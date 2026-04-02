#include "test.h"
#include "lib/mathlib.h"
#include "pathtracer/samplers.h"
#include "util/hdr_image.h"
#include "util/rand.h"

#include <cmath>

Test test_a3_task7_sphere_uniform_pdf("a3.task7.sphere.uniform_pdf", []() {
	Samplers::Sphere::Uniform uni;
	float expected = 1.0f / (4.0f * PI_F);
	Vec3 dirs[] = {
	    Vec3(1.0f, 0.0f, 0.0f).unit(),
	    Vec3(0.0f, 1.0f, 0.0f),
	    Vec3(0.0f, 0.0f, -1.0f),
	    Vec3(0.3f, -0.4f, 0.86f).unit(),
	};
	for (Vec3 d : dirs) {
		if (Test::differs(uni.pdf(d), expected)) {
			throw Test::error("uniform sphere pdf must be 1/(4*pi) everywhere on the sphere");
		}
	}
});

Test test_a3_task7_sphere_uniform_sample_unit("a3.task7.sphere.uniform_unit", []() {
	Samplers::Sphere::Uniform uni;
	RNG rng(0xA37);
	for (int i = 0; i < 600; ++i) {
		Vec3 d = uni.sample(rng);
		if (!d.valid()) {
			throw Test::error("sample direction invalid");
		}
		float n = d.norm();
		if (Test::differs(n, 1.0f)) {
			throw Test::error("uniform sphere samples must be unit length");
		}
	}
});

Test test_a3_task7_sphere_image_sample_valid("a3.task7.sphere.image_sample_valid", []() {
	HDR_Image img(8, 4, Spectrum(0.5f, 0.25f, 0.75f));
	Samplers::Sphere::Image samp(img);
	RNG rng(90210);
	for (int i = 0; i < 200; ++i) {
		Vec3 d = samp.sample(rng);
		if (!d.valid() || d.norm() < 0.99f || d.norm() > 1.01f) {
			throw Test::error("environment importance sample should be a finite unit direction");
		}
		float p = samp.pdf(d);
		if (!std::isfinite(p) || p < 0.0f) {
			throw Test::error("environment pdf should be finite and non-negative");
		}
	}
});

Test test_a3_task7_sphere_image_pdf_matches_sample("a3.task7.sphere.image_pdf_consistency", []() {
	HDR_Image img(6, 3, Spectrum(1.0f, 0.5f, 0.25f));
	Samplers::Sphere::Image samp(img);
	RNG rng(42);
	for (int k = 0; k < 50; ++k) {
		Vec3 d = samp.sample(rng);
		float p1 = samp.pdf(d);
		float p2 = samp.pdf(d.unit());
		if (Test::differs(p1, p2)) {
			throw Test::error("pdf must be invariant to re-normalizing the direction");
		}
	}
});

Test test_a3_task7_sphere_image_black_image_fallback("a3.task7.sphere.image_black_fallback", []() {
	HDR_Image img(3, 3, Spectrum(0.0f));
	Samplers::Sphere::Image samp(img);
	RNG rng(100);
	for (int i = 0; i < 80; ++i) {
		Vec3 d = samp.sample(rng);
		float p = samp.pdf(d);
		if (!std::isfinite(p) || p < 0.0f) {
			throw Test::error("zero-luminance map should still give a valid discrete mixture (uniform over pixels)");
		}
		(void)d;
	}
});
