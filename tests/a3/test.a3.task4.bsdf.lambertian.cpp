#include "test.h"
#include "lib/mathlib.h"
#include "scene/material.h"
#include "scene/texture.h"
#include "util/rand.h"

#include <cmath>

Test test_a3_task4_bsdf_lambertian_simple("a3.task4.bsdf.lambertian.simple", []() {
	// This test just checks that the sample function produces a valid sample.

	auto alb_t = std::make_shared<Texture>(Texture{Textures::Constant{Spectrum{1.0f}}});
	auto bsdf = Materials::Lambertian{alb_t};

	Vec3 out;

	RNG rng(1);

	Materials::Scatter s = bsdf.scatter(rng, out, {});
	// Check that the direction is valid
	if (!s.direction.valid() || s.direction.norm() == 0.0f) {
		throw Test::error("BSDF produced invalid sample!");
	}
	float pdf = bsdf.pdf(out, s.direction);
	// Check that the pdf is valid
	if (!std::isfinite(pdf) || pdf < 0.0f) {
		throw Test::error("BSDF produced sample with invalid pdf!");
	}
	// Check the value against the exact attenuation for the first sample from RNG with seed 1
	if (Test::differs(s.attenuation, Spectrum{0.317861f, 0.317861f, 0.317861f})) {
		throw Test::error("BSDF sample attenuation was not equivalent to evaluate!");
	}
});

Test test_a3_task4_bsdf_lambertian_evaluate_normal("a3.task4.bsdf.lambertian.evaluate_normal", []() {
	// White Lambertian, light from straight above (local +Y): cos_i = 1, BRDF term = (albedo/π)*cos = 1/π.
	auto alb_t = std::make_shared<Texture>(Texture{Textures::Constant{Spectrum{1.0f}}});
	Materials::Lambertian bsdf{alb_t};
	Vec3 out{};
	Vec3 in{0.0f, 1.0f, 0.0f};
	Spectrum e = bsdf.evaluate(out, in, {});
	float one_over_pi = 1.0f / PI_F;
	if (Test::differs(e, Spectrum{one_over_pi, one_over_pi, one_over_pi})) {
		throw Test::error("evaluate at normal incidence should equal albedo/π");
	}
});

Test test_a3_task4_bsdf_lambertian_evaluate_grazing("a3.task4.bsdf.lambertian.evaluate_grazing", []() {
	// cos_i = 1/sqrt(2) for direction (1,1,0)/norm — still upper hemisphere.
	auto alb_t = std::make_shared<Texture>(Texture{Textures::Constant{Spectrum{1.0f}}});
	Materials::Lambertian bsdf{alb_t};
	float s = 1.0f / std::sqrt(2.0f);
	Vec3 in{s, s, 0.0f};
	Vec3 out{};
	Spectrum e = bsdf.evaluate(out, in, {});
	float expected = s / PI_F;
	if (Test::differs(e, Spectrum{expected, expected, expected})) {
		throw Test::error("evaluate at 45° should equal (albedo/π)*cos_i");
	}
});

Test test_a3_task4_bsdf_lambertian_evaluate_below_horizon("a3.task4.bsdf.lambertian.evaluate_below_horizon", []() {
	// Incoming direction in the lower hemisphere should contribute nothing (opaque slab).
	auto alb_t = std::make_shared<Texture>(Texture{Textures::Constant{Spectrum{1.0f}}});
	Materials::Lambertian bsdf{alb_t};
	Vec3 out{};
	Vec3 in{0.0f, -1.0f, 0.0f};
	Spectrum e = bsdf.evaluate(out, in, {});
	if (e.r != 0.0f || e.g != 0.0f || e.b != 0.0f) {
		throw Test::error("evaluate should be zero when cos_i <= 0");
	}
});

Test test_a3_task4_bsdf_lambertian_pdf("a3.task4.bsdf.lambertian.pdf", []() {
	// Cosine hemisphere pdf(ω) = cos θ / π = in.y / π for in.y >= 0.
	auto alb_t = std::make_shared<Texture>(Texture{Textures::Constant{Spectrum{1.0f}}});
	Materials::Lambertian bsdf{alb_t};
	Vec3 out{};
	if (Test::differs(bsdf.pdf(out, Vec3{0.0f, 1.0f, 0.0f}), 1.0f / PI_F)) {
		throw Test::error("pdf at normal should be 1/π");
	}
	float cy = 0.25f;
	if (Test::differs(bsdf.pdf(out, Vec3{0.0f, cy, 0.0f}), cy / PI_F)) {
		throw Test::error("pdf should equal in.y/π on the upper hemisphere");
	}
	if (bsdf.pdf(out, Vec3{0.3f, -0.1f, 0.0f}) != 0.0f) {
		throw Test::error("pdf below the horizon should be 0");
	}
});

Test test_a3_task4_bsdf_lambertian_mc_weight("a3.task4.bsdf.lambertian.mc_weight", []() {
	// With cosine sampling: evaluate = (ρ/π)*cos, pdf = cos/π  =>  evaluate/pdf = ρ  (per channel).
	auto alb_t = std::make_shared<Texture>(
	    Texture{Textures::Constant{Spectrum{0.25f, 0.5f, 1.0f}}});
	Materials::Lambertian bsdf{alb_t};
	Spectrum albedo{0.25f, 0.5f, 1.0f};
	Vec3 out{};
	const Vec3 dirs[] = {
	    Vec3{0.0f, 1.0f, 0.0f},
	    Vec3{0.0f, 0.5f, 0.8660254f},
	    Vec3{0.57735f, 0.57735f, 0.57735f},
	    Vec3{-0.910f, 0.35f, 0.22f},
	};
	for (Vec3 in : dirs) {
		in = in.unit();
		if (in.y <= 0.0f) continue;
		Spectrum ev = bsdf.evaluate(out, in, {});
		float p = bsdf.pdf(out, in);
		if (p <= 0.0f || !std::isfinite(p)) {
			throw Test::error("invalid pdf in mc_weight test");
		}
		Spectrum ratio = ev * (1.0f / p);
		if (Test::differs(ratio, albedo)) {
			throw Test::error("evaluate/pdf should equal albedo for cosine-weighted Lambertian");
		}
	}
});

Test test_a3_task4_bsdf_lambertian_scatter_upper_hemisphere("a3.task4.bsdf.lambertian.scatter_upper_hemisphere", []() {
	// Cosine samples must lie strictly in the upper hemisphere (y > 0) for non-degenerate PDF.
	auto alb_t = std::make_shared<Texture>(Texture{Textures::Constant{Spectrum{1.0f}}});
	Materials::Lambertian bsdf{alb_t};
	Vec3 out{};
	for (uint32_t seed = 0; seed < 64; ++seed) {
		RNG rng(seed);
		Materials::Scatter s = bsdf.scatter(rng, out, {});
		if (s.direction.y <= 1.0e-5f) {
			throw Test::error("cosine scatter should have positive y (upper hemisphere)");
		}
		if (Test::differs(s.attenuation, bsdf.evaluate(out, s.direction, {}))) {
			throw Test::error("scatter attenuation must match evaluate(out, scatter.direction)");
		}
	}
});
