
#include "material.h"
#include "../util/rand.h"

#include <cmath>

namespace Materials {

Vec3 reflect(Vec3 dir) {
	//A3T5 Materials - reflect helper

	// Return direction to incoming light that would be
	// reflected out in direction dir from surface
	// with normal (0,1,0)

	return Vec3(-dir.x, dir.y, -dir.z);
}

Vec3 refract(Vec3 out_dir, float index_of_refraction, bool& was_internal) {
	//A3T5 Materials - refract helper

	// Use Snell's Law to refract out_dir through the surface.
	// Return the refracted direction. Set was_internal to true if
	// refraction does not occur due to total internal reflection,
	// and false otherwise.

	// The surface normal is (0,1,0)

	was_internal = false;
	Vec3 wo = out_dir.unit();
	float cos_o = wo.y;
	float sin2_o = std::max(0.0f, 1.0f - cos_o * cos_o);
	float sin_o = std::sqrt(sin2_o);

	if (sin_o < 1e-6f) {
		float cos_i = cos_o > 0.0f ? -1.0f : 1.0f;
		return Vec3(0.0f, cos_i, 0.0f).unit();
	}

	float eta_i = cos_o > 0.0f ? index_of_refraction : 1.0f;
	float eta_t = cos_o > 0.0f ? 1.0f : index_of_refraction;
	float sin2_i = (eta_t * eta_t) / (eta_i * eta_i) * sin2_o;
	if (sin2_i > 1.0f) {
		was_internal = true;
		return Vec3{};
	}
	float cos_i_mag = std::sqrt(std::max(0.0f, 1.0f - sin2_i));
	float sin_i = std::sqrt(sin2_i);
	float cos_i = cos_o > 0.0f ? -cos_i_mag : cos_i_mag;
	Vec3 h(-wo.x / sin_o, 0.0f, -wo.z / sin_o);
	h.normalize();
	return (h * sin_i + Vec3(0.0f, cos_i, 0.0f)).unit();
}

float schlick(Vec3 in_dir, float index_of_refraction) {
	//A3T5 Materials - Schlick's approximation helper

	// Implement Schlick's approximation of the Fresnel reflection factor.

	if (std::abs(index_of_refraction - 1.0f) < 1e-4f) return 0.0f;
	float cos_theta = std::abs(in_dir.y);
	if (cos_theta > 1.0f) cos_theta = 1.0f;
	float r0 = (1.0f - index_of_refraction) / (1.0f + index_of_refraction);
	r0 = r0 * r0;
	return r0 + (1.0f - r0) * std::pow(1.0f - cos_theta, 5.0f);
}

Spectrum Lambertian::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	//A3T4: Materials - Lambertian BSDF evaluation

	//We get a texture from the weak_ptr albedo lock
	//Then use the evaulate function to get the base color at that texture
	//The normal in local space is (0,1,0) so we can just use the y component 
	//base is surface color, divide by pi to diffuse the light, multiply by cos_i to weight the light
	Spectrum base = albedo.lock()->evaluate(uv);
	float cos_i = std::max(0.0f, in.y);
	return base * (cos_i / PI_F);
}

Scatter Lambertian::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	//A3T4: Materials - Lambertian BSDF scattering
	Samplers::Hemisphere::Cosine sampler;
	Scatter ret;
	ret.direction = sampler.sample(rng);
	ret.attenuation = evaluate(out, ret.direction, uv);
	return ret;
}

float Lambertian::pdf(Vec3 out, Vec3 in) const {
	//A3T4: Materials - Lambertian BSDF probability density function
	Samplers::Hemisphere::Cosine sampler;
	return sampler.pdf(in);
}

Spectrum Lambertian::emission(Vec2 uv) const {
	return {};
}

std::weak_ptr<Texture> Lambertian::display() const {
	return albedo;
}

void Lambertian::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(albedo);
}

Spectrum Mirror::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	return {};
}

Scatter Mirror::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	//A3T5: mirror

	// Use reflect to compute the new direction
	// Don't forget that this is a discrete material!
	// Similar to albedo, reflectance represents the ratio of incoming light to reflected light

	Scatter ret;
	ret.direction = reflect(out);
	ret.attenuation = reflectance.lock()->evaluate(uv);
	return ret;
}

float Mirror::pdf(Vec3 out, Vec3 in) const {
	return 0.0f;
}

Spectrum Mirror::emission(Vec2 uv) const {
	return {};
}

std::weak_ptr<Texture> Mirror::display() const {
	return reflectance;
}

void Mirror::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(reflectance);
}

Spectrum Refract::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	return {};
}

Scatter Refract::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	//A3T5 - refract

	// Use refract to determine the new direction - what happens in the total internal reflection case?
	// Be wary of your eta1/eta2 ratio - are you entering or leaving the surface?
	// Don't forget that this is a discrete material!
	// For attenuation, be sure to take a look at the Specular Transimission section of the PBRT textbook for a derivation
	//  You do not need to scale by the Fresnel Coefficient - you'll only need to account for the correct ratio of indices of refraction

	Scatter ret;
	bool internal = false;
	Vec3 wi = refract(out, ior, internal);
	if (internal) {
		ret.direction = reflect(out);
		ret.attenuation = transmittance.lock()->evaluate(uv);
	} else {
		ret.direction = wi;
		float eta_i = out.y > 0.0f ? ior : 1.0f;
		float eta_t = out.y > 0.0f ? 1.0f : ior;
		float s = eta_t / eta_i;
		ret.attenuation = transmittance.lock()->evaluate(uv) * (s * s);
	}
	return ret;
}

float Refract::pdf(Vec3 out, Vec3 in) const {
	return 0.0f;
}

Spectrum Refract::emission(Vec2 uv) const {
	return {};
}

bool Refract::is_emissive() const {
	return false;
}

bool Refract::is_specular() const {
	return true;
}

bool Refract::is_sided() const {
	return true;
}

std::weak_ptr<Texture> Refract::display() const {
	return transmittance;
}

void Refract::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(transmittance);
}

Spectrum Glass::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	return {};
}

Scatter Glass::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	//A3T5 - glass

	// (1) Compute Fresnel coefficient. Tip: Schlick's approximation.
	// (2) Reflect or refract probabilistically based on Fresnel coefficient. Tip: RNG::coin_flip
	// (3) Compute attenuation based on reflectance or transmittance

	// Be wary of your eta1/eta2 ratio - are you entering or leaving the surface?
	// What happens upon total internal reflection?
	// When debugging Glass, it may be useful to compare to a pure-refraction BSDF
	// For attenuation, be sure to take a look at the Specular Transimission section of the PBRT textbook for a derivation
	//  You do not need to scale by the Fresnel Coefficient - you'll only need to account for the correct ratio of indices of refraction

	Scatter ret;
	float F = schlick(out, ior);
	if (rng.coin_flip(F)) {
		ret.direction = reflect(out);
		ret.attenuation = reflectance.lock()->evaluate(uv);
	} else {
		bool internal = false;
		Vec3 wi = refract(out, ior, internal);
		if (internal) {
			ret.direction = reflect(out);
			ret.attenuation = reflectance.lock()->evaluate(uv);
		} else {
			ret.direction = wi;
			float eta_i = out.y > 0.0f ? ior : 1.0f;
			float eta_t = out.y > 0.0f ? 1.0f : ior;
			float s = eta_t / eta_i;
			ret.attenuation = transmittance.lock()->evaluate(uv) * (s * s);
		}
	}
	return ret;
}

float Glass::pdf(Vec3 out, Vec3 in) const {
	return 0.0f;
}

Spectrum Glass::emission(Vec2 uv) const {
	return {};
}

bool Glass::is_emissive() const {
	return false;
}

bool Glass::is_specular() const {
	return true;
}

bool Glass::is_sided() const {
	return true;
}

std::weak_ptr<Texture> Glass::display() const {
	return transmittance;
}

void Glass::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(reflectance);
	f(transmittance);
}

Spectrum Emissive::evaluate(Vec3 out, Vec3 in, Vec2 uv) const {
	return {};
}

Scatter Emissive::scatter(RNG &rng, Vec3 out, Vec2 uv) const {
	Scatter ret;
	ret.direction = {};
	ret.attenuation = {};
	return ret;
}

float Emissive::pdf(Vec3 out, Vec3 in) const {
	return 0.0f;
}

Spectrum Emissive::emission(Vec2 uv) const {
	return emissive.lock()->evaluate(uv);
}

bool Emissive::is_emissive() const {
	return true;
}

bool Emissive::is_specular() const {
	return true;
}

bool Emissive::is_sided() const {
	return false;
}

std::weak_ptr<Texture> Emissive::display() const {
	return emissive;
}

void Emissive::for_each(const std::function<void(std::weak_ptr<Texture>&)>& f) {
	f(emissive);
}

} // namespace Materials

bool operator!=(const Materials::Lambertian& a, const Materials::Lambertian& b) {
	return a.albedo.lock() != b.albedo.lock();
}

bool operator!=(const Materials::Mirror& a, const Materials::Mirror& b) {
	return a.reflectance.lock() != b.reflectance.lock();
}

bool operator!=(const Materials::Refract& a, const Materials::Refract& b) {
	return a.transmittance.lock() != b.transmittance.lock() || a.ior != b.ior;
}

bool operator!=(const Materials::Glass& a, const Materials::Glass& b) {
	return a.reflectance.lock() != b.reflectance.lock() ||
	       a.transmittance.lock() != b.transmittance.lock() || a.ior != b.ior;
}

bool operator!=(const Materials::Emissive& a, const Materials::Emissive& b) {
	return a.emissive.lock() != b.emissive.lock();
}

bool operator!=(const Material& a, const Material& b) {
	if (a.material.index() != b.material.index()) return false;
	return std::visit(
		[&](const auto& material) {
			return material != std::get<std::decay_t<decltype(material)>>(b.material);
		},
		a.material);
}
