
#include "particles.h"

bool Particles::Particle::update(const PT::Aggregate &scene, Vec3 const &gravity, const float radius, const float dt) {

	float dt_remaining = dt;
	const float eps = 1e-6f;
	uint32_t safety = 16;

	while (dt_remaining > eps && safety-- > 0) {
		float speed = velocity.norm();

		if (speed < eps) {
			velocity += gravity * dt_remaining;
			dt_remaining = 0.0f;
			break;
		}

		Ray ray(position, velocity);
		PT::Trace hit = scene.hit(ray);

		if (hit.hit) {
			float cos_angle = std::abs(dot(ray.dir, hit.normal));
			float offset = (cos_angle > eps) ? (radius / cos_angle) : 0.0f;
			float t_contact = hit.distance - offset;

			if (t_contact < 0.0f) {
				velocity = velocity - 2.0f * dot(velocity, hit.normal) * hit.normal;
				continue;
			}

			float t_travel = t_contact / speed;
			if (t_travel >= dt_remaining) {
				position += velocity * dt_remaining;
				velocity += gravity * dt_remaining;
				dt_remaining = 0.0f;
			} else {
				position += ray.dir * t_contact;
				velocity = velocity - 2.0f * dot(velocity, hit.normal) * hit.normal;
				velocity += gravity * t_travel;
				dt_remaining -= t_travel;
			}
		} else {
			position += velocity * dt_remaining;
			velocity += gravity * dt_remaining;
			dt_remaining = 0.0f;
		}
	}

	age -= dt;
	return age > 0.0f;
}

void Particles::advance(const PT::Aggregate& scene, const Mat4& to_world, float dt) {

	if(step_size < EPS_F) return;

	step_accum += dt;

	while(step_accum > step_size) {
		step(scene, to_world);
		step_accum -= step_size;
	}
}

void Particles::step(const PT::Aggregate& scene, const Mat4& to_world) {

	std::vector<Particle> next;
	next.reserve(particles.size());

	for(Particle& p : particles) {
		if(p.update(scene, gravity, radius, step_size)) {
			next.emplace_back(p);
		}
	}

	if(rate > 0.0f) {

		float cos = std::cos(Radians(spread_angle) / 2.0f);

		double begin_t = current_step * double(step_size) * double(rate);
		double end_t = (current_step + 1) * double(step_size) * double(rate);

		uint64_t begin_i = uint64_t(std::max(0.0, std::ceil(begin_t)));
		uint64_t end_i = uint64_t(std::max(0.0, std::ceil(end_t)));

		for (uint64_t i = begin_i; i < end_i; ++i) {

			float y = lerp(cos, 1.0f, rng.unit());
			float t = 2 * PI_F * rng.unit();
			float d = std::sqrt(1.0f - y * y);
			Vec3 dir = initial_velocity * Vec3(d * std::cos(t), y, d * std::sin(t));

			Particle p;
			p.position = to_world * Vec3(0.0f, 0.0f, 0.0f);
			p.velocity = to_world.rotate(dir);
			p.age = lifetime;
			next.push_back(p);
		}
	}

	particles = std::move(next);
	current_step += 1;
}

void Particles::reset() {
	particles.clear();
	step_accum = 0.0f;
	current_step = 0;
	rng.seed(seed);
}

bool operator!=(const Particles& a, const Particles& b) {
	return a.gravity != b.gravity
	|| a.radius != b.radius
	|| a.initial_velocity != b.initial_velocity
	|| a.spread_angle != b.spread_angle
	|| a.lifetime != b.lifetime
	|| a.rate != b.rate
	|| a.step_size != b.step_size
	|| a.seed != b.seed;
}
