#include "shape.h"
#include "../geometry/util.h"

namespace Shapes {

Vec2 Sphere::uv(Vec3 dir) {
	float u = std::atan2(dir.z, dir.x) / (2.0f * PI_F);
	if (u < 0.0f) u += 1.0f;
	float v = std::acos(-1.0f * std::clamp(dir.y, -1.0f, 1.0f)) / PI_F;
	return Vec2{u, v};
}

BBox Sphere::bbox() const {
	BBox box;
	box.enclose(Vec3(-radius));
	box.enclose(Vec3(radius));
	return box;
}

PT::Trace Sphere::hit(Ray ray) const {
	//A3T2 - sphere hit

    // Sphere centered at origin, radius = Sphere::radius
    // Ray: r(t) = o + t*d
    // Substituting into ||x||^2 - r^2 = 0:
    //   ||d||^2 * t^2 + 2*(o.d)*t + ||o||^2 - r^2 = 0
    // Since ray.dir is normalized, a = ||d||^2 = 1, simplifying to:
    //   t^2 + 2*(o.d)*t + ||o||^2 - r^2 = 0

    float b = 2.0f * dot(ray.point, ray.dir);
    float c = dot(ray.point, ray.point) - radius * radius;

    // discriminant = b^2 - 4c  (a=1 so simplified)
    float disc = b * b - 4.0f * c;

    PT::Trace ret;
    ret.origin = ray.point;
    ret.hit = false;
    ret.distance = 0.0f;
    ret.position = Vec3{};
    ret.normal   = Vec3{};
    ret.uv       = Vec2{};

    // Negative discriminant means no real intersection
    if (disc < 0.0f) return ret;

    float sqrt_disc = std::sqrt(disc);
    float t1 = (-b - sqrt_disc) / 2.0f; // closer intersection
    float t2 = (-b + sqrt_disc) / 2.0f; // farther intersection

    // Pick the smallest t within dist_bounds.
    // If t1 is out of bounds but t2 is in bounds, use t2 (ray origin inside sphere).
    float t = -1.0f;
    if (t1 >= ray.dist_bounds.x && t1 <= ray.dist_bounds.y) {
        t = t1;
    } else if (t2 >= ray.dist_bounds.x && t2 <= ray.dist_bounds.y) {
        t = t2;
    } else {
        return ret; // both intersections out of bounds
    }

    ret.hit      = true;
    ret.distance = t;
    ret.position = ray.at(t);
    // Normal points outward from center (origin), so just normalize the position
    ret.normal   = ret.position.unit();
    ret.uv       = uv(ret.normal);

    return ret;
}

Vec3 Sphere::sample(RNG &rng, Vec3 from) const {
	die("Sampling sphere area lights is not implemented yet.");
}

float Sphere::pdf(Ray ray, Mat4 pdf_T, Mat4 pdf_iT) const {
	die("Sampling sphere area lights is not implemented yet.");
}

Indexed_Mesh Sphere::to_mesh() const {
	return Util::closed_sphere_mesh(radius, 2);
}

} // namespace Shapes

bool operator!=(const Shapes::Sphere& a, const Shapes::Sphere& b) {
	return a.radius != b.radius;
}

bool operator!=(const Shape& a, const Shape& b) {
	if (a.shape.index() != b.shape.index()) return false;
	return std::visit(
		[&](const auto& shape) {
			return shape != std::get<std::decay_t<decltype(shape)>>(b.shape);
		},
		a.shape);
}