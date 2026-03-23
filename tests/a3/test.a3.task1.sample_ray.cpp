#include "test.h"
#include "pathtracer/samplers.h"
#include "scene/camera.h"
#include "util/rand.h"


constexpr uint32_t max_depth = 0;

// Create a camera with the following setups
static std::pair<Camera, Mat4> setup_cam(Vec2 wh, Vec3 cent, Vec3 pos, float fov, float ar) {
	Camera c;
	c.aspect_ratio = ar;
	c.vertical_fov = fov;
	c.film.width = static_cast<uint32_t>(wh.x);
	c.film.height = static_cast<uint32_t>(wh.y);
	c.film.samples = 1;
	c.film.max_ray_depth = max_depth;
	c.near_plane = 0.01f;
	Mat4 world_to_camera = Mat4::look_at(pos, cent, Vec3{0, 1, 0});
	return {c, world_to_camera.inverse()};
}

// -----------------------------------------------------------------------
// Original tests
// -----------------------------------------------------------------------

Test test_a3_task1_sample_ray_simple("a3.task1.sample_ray.simple", []() {
	// Create a camera and get the transform matrix from camera to world
	auto [cam, iV] = setup_cam(Vec2(1, 1), Vec3(0, 0, -1), Vec3(), Degrees(2.0f * std::atan(0.5f)), 1.0f);

	// Create a plane from a point and a normal
	Plane p(Vec3(0, 0, -1), Vec3(0, 0, 1));

	// Number of rays to sample
	constexpr uint32_t N = 100000;

	RNG rng;
	for (uint32_t i = 0; i < N; i++) {
		auto [ret, pdf] = cam.sample_ray(rng, 0, 0);
		ret.transform(iV);

		Line l(ret.point, ret.dir);
		Vec3 hitp;
		if (!p.hit(l, hitp)) {
			throw Test::error("Ray did not hit image plane!");
		}

		Vec2 uv = Vec2{hitp.x, hitp.y} + Vec2{0.5f};
		if (uv.x < 0 || uv.x > 1 || uv.y < 0 || uv.y > 1) {
			throw Test::error("Ray hit outside image plane!");
		}
	}
});

Test test_a3_task1_sample_ray_miss("a3.task1.sample_ray.miss", []() {
	// Create a camera and get the transform matrix from camera to world
	auto [cam, iV] = setup_cam(Vec2(1, 1), Vec3(0, 0, 1), Vec3(), Degrees(2.0f * std::atan(0.5f)), 1.0f);

	// Create a plane from a point and a normal
	Plane p(Vec3(0, 0, -1), Vec3(0, 0, 1));

	// Number of rays to sample
	constexpr uint32_t N = 100000;

	RNG rng;
	for (uint32_t i = 0; i < N; i++) {
		auto [ret, pdf] = cam.sample_ray(rng, 0, 0);
		ret.transform(iV);

		Line l(ret.point, ret.dir);
		Vec3 hitp;
		if (p.hit(l, hitp)) {
			throw Test::error("Ray did hit image plane!");
		}
	}
});

// -----------------------------------------------------------------------
// New tests
// -----------------------------------------------------------------------

// TEST: Ray origin is always at (0,0,0) in camera space.
// sample_ray is specified to return a ray starting at the camera origin.
// A common mistake is accidentally offsetting the origin (e.g. using the
// sensor-plane point as the ray origin instead of zero).
Test test_a3_task1_sample_ray_origin("a3.task1.sample_ray.origin", []() {
	auto [cam, iV] = setup_cam(Vec2(4, 4), Vec3(0, 0, -1), Vec3(), 90.0f, 1.0f);

	constexpr uint32_t N = 1000;
	RNG rng;
	for (uint32_t i = 0; i < N; i++) {
		auto [ray, pdf] = cam.sample_ray(rng, 0, 0);
		// Ray point must be exactly the origin in camera space
		if (ray.point != Vec3(0.0f, 0.0f, 0.0f)) {
			throw Test::error("Ray origin is not at (0,0,0)!");
		}
	}
});

// TEST: Ray direction must be a unit vector.
// sample_ray normalizes the direction before returning it. Forgetting to
// normalize (or normalizing incorrectly) breaks all downstream distance
// calculations in the path tracer.
Test test_a3_task1_sample_ray_normalized("a3.task1.sample_ray.normalized", []() {
	auto [cam, iV] = setup_cam(Vec2(8, 8), Vec3(0, 0, -1), Vec3(), 60.0f, 1.0f);

	constexpr uint32_t N = 1000;
	constexpr float eps = 1e-5f;
	RNG rng;
	for (uint32_t i = 0; i < N; i++) {
		for (uint32_t py = 0; py < 8; py++) {
			for (uint32_t px = 0; px < 8; px++) {
				auto [ray, pdf] = cam.sample_ray(rng, px, py);
				float len = ray.dir.norm();
				if (std::abs(len - 1.0f) > eps) {
					throw Test::error("Ray direction is not a unit vector! Length: " +
					                  std::to_string(len));
				}
			}
		}
	}
});

// TEST: Ray always points into the -Z half-space (camera space).
// The camera looks down -Z, so every generated ray direction must have a
// negative Z component. A positive Z means the ray points behind the camera,
// which would produce completely wrong renders.
Test test_a3_task1_sample_ray_forward("a3.task1.sample_ray.forward", []() {
	auto [cam, iV] = setup_cam(Vec2(16, 16), Vec3(0, 0, -1), Vec3(), 90.0f, 1.0f);

	constexpr uint32_t N = 500;
	RNG rng;
	for (uint32_t py = 0; py < 16; py++) {
		for (uint32_t px = 0; px < 16; px++) {
			for (uint32_t i = 0; i < N; i++) {
				auto [ray, pdf] = cam.sample_ray(rng, px, py);
				if (ray.dir.z >= 0.0f) {
					throw Test::error("Ray direction has non-negative Z in camera space! Z=" +
					                  std::to_string(ray.dir.z));
				}
			}
		}
	}
});

// TEST: Rays from different pixels must point in different directions.
// If sample_ray ignores px/py and always returns the same direction,
// the render will look like a single point sample replicated everywhere.
// We check that the center-of-mass directions of opposite corner pixels
// are clearly distinct.
Test test_a3_task1_sample_ray_distinct_pixels("a3.task1.sample_ray.distinct_pixels", []() {
	auto [cam, iV] = setup_cam(Vec2(4, 4), Vec3(0, 0, -1), Vec3(), 90.0f, 1.0f);

	constexpr uint32_t N = 10000;
	constexpr float eps = 0.05f; // pixels should differ by much more than this
	RNG rng;

	// Accumulate mean directions for top-left (0,0) and bottom-right (3,3)
	Vec3 mean_tl, mean_br;
	for (uint32_t i = 0; i < N; i++) {
		auto [r00, p00] = cam.sample_ray(rng, 0, 0);
		auto [r33, p33] = cam.sample_ray(rng, 3, 3);
		mean_tl = mean_tl + r00.dir;
		mean_br = mean_br + r33.dir;
	}
	mean_tl = mean_tl / float(N);
	mean_br = mean_br / float(N);

	float diff_x = std::abs(mean_tl.x - mean_br.x);
	float diff_y = std::abs(mean_tl.y - mean_br.y);
	if (diff_x < eps && diff_y < eps) {
		throw Test::error("Rays from (0,0) and (3,3) point in nearly the same direction! "
		                  "sample_ray may be ignoring px/py.");
	}
});

// TEST: Rays for the same pixel stay within that pixel's footprint on the sensor.
// sample_ray adds a random sub-pixel offset, but that offset must stay inside
// [px, px+1) x [py, py+1). Offsets that spill into neighboring pixels break
// supersampling and produce aliasing artifacts.
// Strategy: for pixel (0,0) of a 1x1 film with FOV = 2*atan(0.5), the sensor
// plane at z=-1 goes from (-0.5,-0.5) to (0.5,0.5). All hit points must stay
// within that square.
Test test_a3_task1_sample_ray_within_pixel("a3.task1.sample_ray.within_pixel", []() {
	// 1x1 film: the entire sensor is one pixel, so any in-bounds ray must hit
	// the unit square centered at the origin on the z=-1 plane.
	auto [cam, iV] = setup_cam(Vec2(1, 1), Vec3(0, 0, -1), Vec3(), Degrees(2.0f * std::atan(0.5f)), 1.0f);

	Plane sensor_plane(Vec3(0, 0, -1), Vec3(0, 0, 1));
	constexpr uint32_t N = 100000;
	constexpr float eps = 1e-4f;
	RNG rng;

	for (uint32_t i = 0; i < N; i++) {
		auto [ray, pdf] = cam.sample_ray(rng, 0, 0);
		// Transform to world (identity here since camera is at origin looking at -z)
		ray.transform(iV);

		Line l(ray.point, ray.dir);
		Vec3 hitp;
		if (!sensor_plane.hit(l, hitp)) {
			throw Test::error("Ray did not hit sensor plane!");
		}
		// hitp should be within [-0.5-eps, 0.5+eps] in x and y
		if (hitp.x < -0.5f - eps || hitp.x > 0.5f + eps ||
		    hitp.y < -0.5f - eps || hitp.y > 0.5f + eps) {
			throw Test::error("Ray hit outside the pixel's sensor footprint! Hit: (" +
			                  std::to_string(hitp.x) + ", " + std::to_string(hitp.y) + ")");
		}
	}
});

// TEST: PDF returned must be positive for all samples.
// The path tracer divides by the PDF when accumulating samples. A zero or
// negative PDF causes division-by-zero or sign errors in the energy estimate.
Test test_a3_task1_sample_ray_pdf_positive("a3.task1.sample_ray.pdf_positive", []() {
	auto [cam, iV] = setup_cam(Vec2(8, 8), Vec3(0, 0, -1), Vec3(), 60.0f, 1.0f);

	constexpr uint32_t N = 1000;
	RNG rng;
	for (uint32_t py = 0; py < 8; py++) {
		for (uint32_t px = 0; px < 8; px++) {
			for (uint32_t i = 0; i < N; i++) {
				auto [ray, pdf] = cam.sample_ray(rng, px, py);
				if (pdf <= 0.0f) {
					throw Test::error("sample_ray returned a non-positive PDF: " +
					                  std::to_string(pdf));
				}
			}
		}
	}
});

// TEST: Aspect ratio is respected.
// A wide-screen camera (ar=2) should produce rays whose X spread is twice
// the Y spread. We measure this by comparing the standard deviation of the
// X and Y hit positions on the sensor plane and checking their ratio.
Test test_a3_task1_sample_ray_aspect_ratio("a3.task1.sample_ray.aspect_ratio", []() {
	// Single-pixel wide camera, aspect ratio 2:1, FOV chosen so
	// sensor half-height = 0.5 and half-width = 1.0 at z=-1.
	float fov = Degrees(2.0f * std::atan(0.5f)); // half_height = 0.5
	auto [cam, iV] = setup_cam(Vec2(1, 1), Vec3(0, 0, -1), Vec3(), fov, 2.0f);

	Plane sensor_plane(Vec3(0, 0, -1), Vec3(0, 0, 1));
	constexpr uint32_t N = 100000;
	RNG rng;

	float sum_x = 0.0f, sum_y = 0.0f;
	float sum_x2 = 0.0f, sum_y2 = 0.0f;

	for (uint32_t i = 0; i < N; i++) {
		auto [ray, pdf] = cam.sample_ray(rng, 0, 0);
		ray.transform(iV);
		Line l(ray.point, ray.dir);
		Vec3 hitp;
		if (!sensor_plane.hit(l, hitp)) {
			throw Test::error("Ray missed sensor plane during aspect ratio test!");
		}
		sum_x  += hitp.x;  sum_y  += hitp.y;
		sum_x2 += hitp.x * hitp.x;
		sum_y2 += hitp.y * hitp.y;
	}

	float mean_x = sum_x / N, mean_y = sum_y / N;
	float var_x  = sum_x2 / N - mean_x * mean_x;
	float var_y  = sum_y2 / N - mean_y * mean_y;
	// std_x / std_y should equal aspect_ratio = 2
	float ratio = std::sqrt(var_x / (var_y + 1e-9f));
	if (std::abs(ratio - 2.0f) > 0.1f) {
		throw Test::error("Aspect ratio not respected! Expected X/Y spread ratio ~2.0, got " +
		                  std::to_string(ratio));
	}
});

// TEST: max_ray_depth is correctly propagated to the returned ray.
// do_trace relies on ray.depth to limit recursion. If sample_ray forgets to
// set it, every ray will have depth 0 and indirect lighting will be skipped.
Test test_a3_task1_sample_ray_depth("a3.task1.sample_ray.depth", []() {
	Camera cam;
	cam.aspect_ratio = 1.0f;
	cam.vertical_fov = 90.0f;
	cam.film.width  = 4;
	cam.film.height = 4;
	cam.film.samples = 1;
	cam.film.max_ray_depth = 7; // non-default depth to catch bugs
	cam.near_plane = 0.01f;

	RNG rng;
	for (uint32_t py = 0; py < 4; py++) {
		for (uint32_t px = 0; px < 4; px++) {
			auto [ray, pdf] = cam.sample_ray(rng, px, py);
			if (ray.depth != 7) {
				throw Test::error("Ray depth not set correctly! Expected 7, got " +
				                  std::to_string(ray.depth));
			}
		}
	}
});