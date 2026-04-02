
#include "samplers.h"
#include "../scene/shape.h"
#include "../util/rand.h"

#include <algorithm>
#include <cmath>

constexpr bool IMPORTANCE_SAMPLING = true;

namespace Samplers {

Vec2 Rect::sample(RNG &rng) const {
	//A3T1 - step 2 - supersampling

    // Return a point selected uniformly at random from the rectangle [0,size.x)x[0,size.y)
    // Useful function: rng.unit()

	return Vec2{rng.unit() * size.x, rng.unit() * size.y};
}

float Rect::pdf(Vec2 at) const {
	if (at.x < 0.0f || at.x > size.x || at.y < 0.0f || at.y > size.y) return 0.0f;
	return 1.0f / (size.x * size.y);
}

Vec2 Circle::sample(RNG &rng) const {
	//A3EC - bokeh - circle sampling

    // Return a point selected uniformly at random from a circle defined by its
	// center and radius.
    // Useful function: rng.unit()

    return Vec2{};
}

float Circle::pdf(Vec2 at) const {
	//A3EC - bokeh - circle pdf

	// Return the pdf of sampling the point 'at' for a circle defined by its
	// center and radius.

    return 1.f;
}

Vec3 Point::sample(RNG &rng) const {
	return point;
}

float Point::pdf(Vec3 at) const {
	return at == point ? 1.0f : 0.0f;
}

Vec3 Triangle::sample(RNG &rng) const {
	float u = std::sqrt(rng.unit());
	float v = rng.unit();
	float a = u * (1.0f - v);
	float b = u * v;
	return a * v0 + b * v1 + (1.0f - a - b) * v2;
}

float Triangle::pdf(Vec3 at) const {
	float a = 0.5f * cross(v1 - v0, v2 - v0).norm();
	float u = 0.5f * cross(at - v1, at - v2).norm() / a;
	float v = 0.5f * cross(at - v2, at - v0).norm() / a;
	float w = 1.0f - u - v;
	if (u < 0.0f || v < 0.0f || w < 0.0f) return 0.0f;
	if (u > 1.0f || v > 1.0f || w > 1.0f) return 0.0f;
	return 1.0f / a;
}

Vec3 Hemisphere::Uniform::sample(RNG &rng) const {

	float Xi1 = rng.unit();
	float Xi2 = rng.unit();

	float theta = std::acos(Xi1);
	float phi = 2.0f * PI_F * Xi2;

	float xs = std::sin(theta) * std::cos(phi);
	float ys = std::cos(theta);
	float zs = std::sin(theta) * std::sin(phi);

	return Vec3(xs, ys, zs);
}

float Hemisphere::Uniform::pdf(Vec3 dir) const {
	if (dir.y < 0.0f) return 0.0f;
	return 1.0f / (2.0f * PI_F);
}

Vec3 Hemisphere::Cosine::sample(RNG &rng) const {

	float phi = rng.unit() * 2.0f * PI_F;
	float cos_t = std::sqrt(rng.unit());

	float sin_t = std::sqrt(1 - cos_t * cos_t);
	float x = std::cos(phi) * sin_t;
	float z = std::sin(phi) * sin_t;
	float y = cos_t;

	return Vec3(x, y, z);
}

float Hemisphere::Cosine::pdf(Vec3 dir) const {
	if (dir.y < 0.0f) return 0.0f;
	return dir.y / PI_F;
}

Vec3 Sphere::Uniform::sample(RNG &rng) const {
	//A3T7 - sphere sampler

	// Generate a uniformly random point on the unit sphere.
	// Tip: start with Hemisphere::Uniform

	float u = rng.unit();
	float v = rng.unit();
	float y = 1.0f - 2.0f * u;
	float r = std::sqrt(std::max(0.0f, 1.0f - y * y));
	float phi = 2.0f * PI_F * v;
	return Vec3(r * std::cos(phi), y, r * std::sin(phi));
}

float Sphere::Uniform::pdf(Vec3 dir) const {
	return 1.0f / (4.0f * PI_F);
}

Sphere::Image::Image(const HDR_Image& image) {
	//A3T7 - image sampler init

	// Set up importance sampling data structures for a spherical environment map image.
	// You may make use of the _pdf, _cdf, and total members, or create your own.

	const auto [_w, _h] = image.dimension();
	w = _w;
	h = _h;
	jitter = Rect(Vec2(w > 0 ? 1.0f / float(w) : 1.0f, h > 0 ? 1.0f / float(h) : 1.0f));

	const uint32_t n = w * h;
	_pdf.assign(n, 0.0f);
	float sum = 0.0f;
	for (uint32_t py = 0; py < h; ++py) {
		float vc = (py + 0.5f) / float(h);
		float y = -std::cos(PI_F * vc);
		float st = std::sqrt(std::max(0.0f, 1.0f - y * y));
		for (uint32_t px = 0; px < w; ++px) {
			float wt = image.at(px, py).luma() * st;
			_pdf[py * w + px] = wt;
			sum += wt;
		}
	}
	if (sum <= 1e-20f || !std::isfinite(sum)) {
		float u = 1.0f / float(n);
		std::fill(_pdf.begin(), _pdf.end(), u);
		sum = 1.0f;
	} else {
		for (float& p : _pdf) p /= sum;
	}
	_cdf.resize(n);
	float c = 0.0f;
	for (uint32_t i = 0; i < n; ++i) {
		c += _pdf[i];
		_cdf[i] = c;
	}
}

Vec3 Sphere::Image::sample(RNG &rng) const {
	Sphere::Uniform uni;
	if (!IMPORTANCE_SAMPLING) {
		// Step 1: Uniform sampling
		// Declare a uniform sampler and return its sample
		return uni.sample(rng);
	}
	// Step 2: Importance sampling
	// Use your importance sampling data structure to generate a sample direction.
	// Tip: std::upper_bound
	float u = rng.unit();
	auto it = std::lower_bound(_cdf.begin(), _cdf.end(), u);
	uint32_t idx = static_cast<uint32_t>(it - _cdf.begin());
	if (idx >= w * h) idx = w * h - 1;
	uint32_t px = idx % w;
	uint32_t py = idx / w;
	Vec2 sub = jitter.sample(rng);
	float uf = float(px) / float(w) + sub.x;
	float vf = float(py) / float(h) + sub.y;
	float phi = 2.0f * PI_F * uf;
	float y = -std::cos(PI_F * vf);
	float s = std::sqrt(std::max(0.0f, 1.0f - y * y));
	return Vec3(s * std::cos(phi), y, s * std::sin(phi)).unit();
}

float Sphere::Image::pdf(Vec3 dir) const {
	Sphere::Uniform uni;
	if (!IMPORTANCE_SAMPLING) {
		// Step 1: Uniform sampling
		// Declare a uniform sampler and return its pdf
		return uni.pdf(dir);
	}
	// A3T7 - image sampler importance sampling pdf
	// What is the PDF of this distribution at a particular direction?
	if (!dir.valid()) return 0.0f;
	Vec3 d = dir.unit();
	float sin_theta = std::sqrt(std::max(0.0f, 1.0f - d.y * d.y));
	if (sin_theta < 1e-12f) return 0.0f;
	Vec2 uv = Shapes::Sphere::uv(d);
	float uu = std::clamp(uv.x, 0.0f, 1.0f - 1e-7f);
	float vv = std::clamp(uv.y, 0.0f, 1.0f - 1e-7f);
	uint32_t px = std::min(w - 1, static_cast<uint32_t>(uu * float(w)));
	uint32_t py = std::min(h - 1, static_cast<uint32_t>(vv * float(h)));
	uint32_t idx = py * w + px;
	float p_pix = _pdf[idx];
	float jacobian = (float(w) * float(h)) / (2.0f * PI_F * PI_F * sin_theta);
	return p_pix * jacobian;
}

} // namespace Samplers
