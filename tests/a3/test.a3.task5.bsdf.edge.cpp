#include "test.h"
#include "scene/material.h"
#include "scene/texture.h"
#include "util/rand.h"

#include <cmath>
#include <string>

static Vec3 random_unit(RNG &rng) {
	for (;;) {
		float x = rng.unit() * 2.0f - 1.0f;
		float y = rng.unit() * 2.0f - 1.0f;
		float z = rng.unit() * 2.0f - 1.0f;
		Vec3 v(x, y, z);
		float n = v.norm();
		if (n > 1e-4f) return v / n;
	}
}

Test test_a3_task5_reflect_formula("a3.task5.reflect.formula", []() {
	Vec3 a(0.0f, 1.0f, 0.0f);
	Vec3 r = Materials::reflect(a);
	if (Test::differs(r, Vec3(0.0f, 1.0f, 0.0f))) {
		throw Test::error("reflect(+Y) should fix the tangent-parallel components; normal incidence matches.");
	}
	Vec3 b(0.7f, 0.3f, -0.65f);
	b = b.unit();
	r = Materials::reflect(b);
	if (Test::differs(r, Vec3(-b.x, b.y, -b.z))) {
		throw Test::error("reflect must negate x and z, preserve y");
	}
});

Test test_a3_task5_refract_ior_one_inverse("a3.task5.refract.ior_one_inverse", []() {
	const Spectrum tm(0.2f, 0.4f, 0.6f);
	auto tm_t = std::make_shared<Texture>(Texture{Textures::Constant{tm}});
	Materials::Refract bsdf{tm_t, 1.0f};
	for (uint32_t seed = 0; seed < 128; ++seed) {
		RNG rng(seed);
		Vec3 out = random_unit(rng);
		Materials::Scatter s = bsdf.scatter(rng, out, {});
		if (Test::differs(s.direction, out * -1.0f)) {
			throw Test::error("ior=1 refract must give incoming = -outward");
		}
		if (Test::differs(s.attenuation, tm)) {
			throw Test::error("ior=1 refract attenuation should be transmittance * 1^2");
		}
	}
});

Test test_a3_task5_refract_tir("a3.task5.refract.tir", []() {
	const Spectrum tm(0.5f, 0.5f, 0.5f);
	auto tm_t = std::make_shared<Texture>(Texture{Textures::Constant{tm}});
	Materials::Refract bsdf{tm_t, 1.5f};
	Vec3 out = Vec3(0.99f, -0.14142136f, 0.0f);
	out.normalize();
	RNG rng(1);
	Materials::Scatter s = bsdf.scatter(rng, out, {});
	Vec3 expected_reflect = Materials::reflect(out);
	if (Test::differs(s.direction, expected_reflect)) {
		throw Test::error("TIR should mirror-reflect");
	}
	if (Test::differs(s.attenuation, tm)) {
		throw Test::error("TIR path should use transmittance color as stub attenuation");
	}
});

Test test_a3_task5_schlick_edge("a3.task5.schlick.edge", []() {
	float f0 = Materials::schlick(Vec3(0.0f, 1.0f, 0.0f), 1.5f);
	float r0 = (1.0f - 1.5f) / (1.0f + 1.5f);
	r0 = r0 * r0;
	if (Test::differs(f0, r0)) {
		throw Test::error("Schlick at normal incidence should equal R0");
	}
	if (Materials::schlick(Vec3(0.0f, 1.0f, 0.0f), 1.0f) != 0.0f) {
		throw Test::error("Schlick with ior=1 must be 0");
	}
	float grazing = Materials::schlick(Vec3(0.0f, 0.001f, 0.9999995f).unit(), 1.5f);
	if (!(grazing > f0)) {
		throw Test::error("Schlick should increase toward grazing angles");
	}
});

Test test_a3_task5_mirror_stress("a3.task5.mirror.stress", []() {
	auto rfl_t = std::make_shared<Texture>(Texture{Textures::Constant{Spectrum{1.0f}}});
	Materials::Mirror bsdf{rfl_t};
	RNG rng(0xC001D00D);
	for (int i = 0; i < 800; ++i) {
		Vec3 out = random_unit(rng);
		Materials::Scatter s = bsdf.scatter(rng, out, {});
		Vec3 want(-out.x, out.y, -out.z);
		if (Test::differs(s.direction, want)) {
			throw Test::error("mirror scatter direction mismatch at iteration " + std::to_string(i));
		}
		if (Test::differs(s.attenuation, Spectrum{1.0f, 1.0f, 1.0f})) {
			throw Test::error("mirror attenuation should match white reflectance");
		}
	}
});

Test test_a3_task5_refract_stress_snell("a3.task5.refract.stress_snell", []() {
	const Spectrum tm(0.25f, 0.5f, 0.75f);
	auto tm_t = std::make_shared<Texture>(Texture{Textures::Constant{tm}});
	const float ior = 1.33f;
	Materials::Refract bsdf{tm_t, ior};
	RNG rng(0xBEEF);
	for (int i = 0; i < 400; ++i) {
		Vec3 out = random_unit(rng);
		bool internal = false;
		Vec3 wi = Materials::refract(out, ior, internal);
		if (internal) continue;
		float oxy = out.x * out.x + out.z * out.z;
		float wixy = wi.x * wi.x + wi.z * wi.z;
		if (oxy > 1e-8f && wixy > 1e-8f) {
			float sin_o = std::sqrt(oxy);
			float sin_i = std::sqrt(wixy);
			float eta_i = out.y > 0.0f ? ior : 1.0f;
			float eta_t = out.y > 0.0f ? 1.0f : ior;
			float lhs = eta_i * sin_i;
			float rhs = eta_t * sin_o;
			if (std::abs(lhs - rhs) > 0.02f) {
				throw Test::error("Snell invariant failed at iteration " + std::to_string(i));
			}
		}
	}
});

Test test_a3_task5_glass_ior_one_always_refract("a3.task5.glass.ior_one_refract", []() {
	auto tm_t = std::make_shared<Texture>(Texture{Textures::Constant{Spectrum{0.25f, 0.5f, 0.75f}}});
	auto rfl_t = std::make_shared<Texture>(Texture{Textures::Constant{Spectrum{1.0f, 0.0f, 0.0f}}});
	Materials::Glass bsdf{tm_t, rfl_t, 1.0f};
	for (uint32_t seed = 0; seed < 64; ++seed) {
		RNG rng(seed);
		Vec3 out = random_unit(rng);
		Materials::Scatter s = bsdf.scatter(rng, out, {});
		if (Test::differs(s.direction, out * -1.0f)) {
			throw Test::error("ior=1 glass must refract (F=0), never mirror");
		}
	}
});
