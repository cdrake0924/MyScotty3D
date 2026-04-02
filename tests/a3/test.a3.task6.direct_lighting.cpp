#include "test.h"
#include "scene/material.h"
#include "scene/texture.h"
#include "util/rand.h"

#include <cmath>

namespace {

float mixture_pdf(float p_bsdf, float p_light) {
	return 0.5f * p_bsdf + 0.5f * p_light;
}

} // namespace

Test test_a3_task6_mixture_pdf_average("a3.task6.mixture_pdf.average", []() {
	float p0 = 0.2f;
	float p1 = 0.8f;
	float pm = mixture_pdf(p0, p1);
	if (Test::differs(pm, 0.5f)) {
		throw Test::error("mixture pdf should be the arithmetic mean of the two strategies");
	}
});

Test test_a3_task6_mixture_pdf_symmetric("a3.task6.mixture_pdf.symmetric", []() {
	float a = 0.37f;
	float b = 0.91f;
	if (Test::differs(mixture_pdf(a, b), mixture_pdf(b, a))) {
		throw Test::error("mixture pdf must be symmetric in its arguments");
	}
});

Test test_a3_task6_mixture_estimator_weight("a3.task6.mixture.estimator_weight", []() {
	Spectrum f(0.3f, 0.6f, 0.1f);
	Spectrum L(2.0f, 2.0f, 2.0f);
	float p_bsdf = 0.1f;
	float p_light = 0.3f;
	float p_mix = mixture_pdf(p_bsdf, p_light);
	Spectrum contrib = f * L * (1.0f / p_mix);
	float expect_scale = 1.0f / p_mix;
	if (Test::differs(contrib.r, f.r * L.r * expect_scale) ||
	    Test::differs(contrib.g, f.g * L.g * expect_scale) ||
	    Test::differs(contrib.b, f.b * L.b * expect_scale)) {
		throw Test::error("single-sample MIS contribution should be (f * L) / p_mix");
	}
});

Test test_a3_task6_mixture_pdf_nonnegative("a3.task6.mixture_pdf.nonnegative", []() {
	if (mixture_pdf(0.0f, 0.0f) < 0.0f) {
		throw Test::error("mixture pdf should be non-negative");
	}
	float pm = mixture_pdf(0.0f, 0.4f);
	if (pm < 0.0f || !std::isfinite(pm)) {
		throw Test::error("mixture pdf with one zero strategy should still be valid");
	}
});

Test test_a3_task6_bsdf_evaluate_matches_task4_lambertian("a3.task6.task4_lambertian_consistency", []() {
	auto alb = std::make_shared<Texture>(Texture{Textures::Constant{Spectrum{0.8f, 0.2f, 0.5f}}});
	Materials::Lambertian bsdf{alb};
	Vec3 out(0.2f, 0.9f, -0.1f);
	out = out.unit();
	RNG rng(42);
	Materials::Scatter s = bsdf.scatter(rng, out, {});
	Spectrum from_scatter = s.attenuation;
	Spectrum from_eval = bsdf.evaluate(out, s.direction, {});
	if (Test::differs(from_scatter, from_eval)) {
		throw Test::error("Lambertian scatter attenuation must match evaluate (task 6 uses evaluate for both strategies)");
	}
});
