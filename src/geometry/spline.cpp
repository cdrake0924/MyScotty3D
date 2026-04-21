
#include "../geometry/spline.h"

template<typename T> T Spline<T>::at(float time) const {

	if (knots.empty()) return T();
	if (knots.size() == 1) return knots.begin()->second;

	auto first_it = knots.begin();
	if (time <= first_it->first) return first_it->second;

	auto last_it = std::prev(knots.end());
	if (time >= last_it->first) return last_it->second;

	auto k2_it = knots.upper_bound(time);
	auto k1_it = std::prev(k2_it);

	float t1 = k1_it->first;
	float t2 = k2_it->first;
	T p1 = k1_it->second;
	T p2 = k2_it->second;

	float t0;
	T p0;
	if (k1_it == first_it) {
		t0 = t1 - (t2 - t1);
		p0 = p1 - (p2 - p1);
	} else {
		auto k0_it = std::prev(k1_it);
		t0 = k0_it->first;
		p0 = k0_it->second;
	}

	float t3;
	T p3;
	auto k3_it = std::next(k2_it);
	if (k3_it == knots.end()) {
		t3 = t2 + (t2 - t1);
		p3 = p2 + (p2 - p1);
	} else {
		t3 = k3_it->first;
		p3 = k3_it->second;
	}

	T m0 = (p2 - p0) * (1.0f / (t2 - t0));
	T m1 = (p3 - p1) * (1.0f / (t3 - t1));

	float dt = t2 - t1;
	float u = (time - t1) / dt;

	return cubic_unit_spline(u, p1, p2, m0 * dt, m1 * dt);
}

template<typename T>
T Spline<T>::cubic_unit_spline(float time, const T& position0, const T& position1,
                               const T& tangent0, const T& tangent1) {

	float t = time;
	float t2 = t * t;
	float t3 = t2 * t;

	float h00 = 2.0f * t3 - 3.0f * t2 + 1.0f;
	float h10 = t3 - 2.0f * t2 + t;
	float h01 = -2.0f * t3 + 3.0f * t2;
	float h11 = t3 - t2;

	return h00 * position0 + h10 * tangent0 + h01 * position1 + h11 * tangent1;
}

template class Spline<float>;
template class Spline<double>;
template class Spline<Vec4>;
template class Spline<Vec3>;
template class Spline<Vec2>;
template class Spline<Mat4>;
template class Spline<Spectrum>;
