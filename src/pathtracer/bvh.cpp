
#include "bvh.h"
#include "aggregate.h"
#include "instance.h"
#include "tri_mesh.h"

#include <algorithm>
#include <array>
#include <limits>
#include <stack>

namespace PT {

struct BVHBuildData {
	BVHBuildData(size_t start, size_t range, size_t dst) : start(start), range(range), node(dst) {
	}
	size_t start; ///< start index into the primitive array
	size_t range; ///< range of index into the primitive array
	size_t node;  ///< address to update
};

struct SAHBucketData {
	BBox bb;          ///< bbox of all primitives
	size_t num_prims; ///< number of primitives in the bucket
};

template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	//A3T3 - build a bvh

	// Keep these
    nodes.clear();
    primitives = std::move(prims);

    // Construct a BVH from the given vector of primitives and maximum leaf
    // size configuration.

	if (primitives.empty()) {
		root_idx = 0;
		return;
	}
	if (max_leaf_size == 0) max_leaf_size = 1;

	root_idx = new_node();
	std::stack<BVHBuildData> build_stack;
	build_stack.emplace(0, primitives.size(), root_idx);

	constexpr size_t n_buckets = 12;
	constexpr float traversal_cost = 1.0f;

	while (!build_stack.empty()) {
		BVHBuildData build = build_stack.top();
		build_stack.pop();

		BBox node_bbox;
		BBox centroid_bbox;
		for (size_t i = build.start; i < build.start + build.range; i++) {
			BBox prim_bbox = primitives[i].bbox();
			node_bbox.enclose(prim_bbox);
			centroid_bbox.enclose(prim_bbox.center());
		}

		if (build.range <= max_leaf_size || centroid_bbox.empty()) {
			nodes[build.node] = Node{node_bbox, build.start, build.range, build.node, build.node};
			continue;
		}

		float best_cost = std::numeric_limits<float>::infinity();
		size_t best_axis = 0;
		size_t best_bucket_split = 0;
		bool found_split = false;

		for (size_t axis = 0; axis < 3; axis++) {
			uint32_t axis_u = static_cast<uint32_t>(axis);
			float cmin = centroid_bbox.min[axis_u];
			float cmax = centroid_bbox.max[axis_u];
			float extent = cmax - cmin;
			if (extent <= 0.0f) continue;

			std::array<SAHBucketData, n_buckets> buckets;
			for (auto& bucket : buckets) {
				bucket.bb.reset();
				bucket.num_prims = 0;
			}

			for (size_t i = build.start; i < build.start + build.range; i++) {
				BBox prim_bbox = primitives[i].bbox();
				float c = prim_bbox.center()[axis_u];
				size_t b = static_cast<size_t>(n_buckets * ((c - cmin) / extent));
				if (b >= n_buckets) b = n_buckets - 1;
				buckets[b].num_prims++;
				buckets[b].bb.enclose(prim_bbox);
			}

			std::array<BBox, n_buckets - 1> left_bbox;
			std::array<BBox, n_buckets - 1> right_bbox;
			std::array<size_t, n_buckets - 1> left_count{};
			std::array<size_t, n_buckets - 1> right_count{};

			BBox running_left;
			size_t running_left_count = 0;
			for (size_t i = 0; i + 1 < n_buckets; i++) {
				running_left_count += buckets[i].num_prims;
				running_left.enclose(buckets[i].bb);
				left_count[i] = running_left_count;
				left_bbox[i] = running_left;
			}

			BBox running_right;
			size_t running_right_count = 0;
			for (size_t i = n_buckets - 1; i > 0; i--) {
				running_right_count += buckets[i].num_prims;
				running_right.enclose(buckets[i].bb);
				right_count[i - 1] = running_right_count;
				right_bbox[i - 1] = running_right;
			}

			float parent_area = node_bbox.surface_area();
			if (parent_area <= 0.0f) continue;

			for (size_t i = 0; i + 1 < n_buckets; i++) {
				if (left_count[i] == 0 || right_count[i] == 0) continue;
				float sah =
				    traversal_cost +
				    (left_bbox[i].surface_area() * static_cast<float>(left_count[i]) +
				     right_bbox[i].surface_area() * static_cast<float>(right_count[i])) /
				        parent_area;
				if (sah < best_cost) {
					best_cost = sah;
					best_axis = axis;
					best_bucket_split = i;
					found_split = true;
				}
			}
		}

		size_t mid = build.start + build.range / 2;
		if (found_split) {
			uint32_t best_axis_u = static_cast<uint32_t>(best_axis);
			float cmin = centroid_bbox.min[best_axis_u];
			float cmax = centroid_bbox.max[best_axis_u];
			float extent = cmax - cmin;
			if (extent > 0.0f) {
				auto begin = primitives.begin() + static_cast<std::ptrdiff_t>(build.start);
				auto end = begin + static_cast<std::ptrdiff_t>(build.range);
				auto split_at = std::partition(begin, end, [&](const Primitive& prim) {
					size_t bucket = static_cast<size_t>(
					    n_buckets * ((prim.bbox().center()[best_axis_u] - cmin) / extent));
					if (bucket >= n_buckets) bucket = n_buckets - 1;
					return bucket <= best_bucket_split;
				});
				mid = build.start + static_cast<size_t>(split_at - begin);
			}
		}

		// Guard against degenerate partitions where one side gets everything.
		if (mid == build.start || mid == build.start + build.range) {
			Vec3 diag = centroid_bbox.max - centroid_bbox.min;
			size_t axis = 0;
			if (diag.y > diag.x && diag.y >= diag.z) axis = 1;
			else if (diag.z > diag.x && diag.z > diag.y) axis = 2;

			auto begin = primitives.begin() + static_cast<std::ptrdiff_t>(build.start);
			auto nth = begin + static_cast<std::ptrdiff_t>(build.range / 2);
			auto end = begin + static_cast<std::ptrdiff_t>(build.range);
			uint32_t axis_u = static_cast<uint32_t>(axis);
			std::nth_element(begin, nth, end, [&](const Primitive& a, const Primitive& b) {
				return a.bbox().center()[axis_u] < b.bbox().center()[axis_u];
			});
			mid = build.start + build.range / 2;
		}

		size_t l_idx = new_node();
		size_t r_idx = new_node();
		nodes[build.node] = Node{node_bbox, build.start, build.range, l_idx, r_idx};

		size_t l_range = mid - build.start;
		size_t r_range = build.range - l_range;
		build_stack.emplace(mid, r_range, r_idx);
		build_stack.emplace(build.start, l_range, l_idx);
	}

}

template<typename Primitive> Trace BVH<Primitive>::hit(const Ray& ray) const {
	//A3T3 - traverse your BVH

    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.

	if (nodes.empty()) return {};

	Vec2 root_times = ray.dist_bounds;
	if (!nodes[root_idx].bbox.hit(ray, root_times)) return {};

	Trace ret;
	float closest = ray.dist_bounds.y;
	std::stack<size_t> tstack;
	tstack.push(root_idx);

	while (!tstack.empty()) {
		size_t idx = tstack.top();
		tstack.pop();
		const Node& node = nodes[idx];

		Vec2 node_times = ray.dist_bounds;
		node_times.y = std::min(node_times.y, closest);
		if (!node.bbox.hit(ray, node_times)) continue;

		if (node.is_leaf()) {
			for (size_t i = node.start; i < node.start + node.size; i++) {
				Ray clipped = ray;
				clipped.dist_bounds.y = std::min(clipped.dist_bounds.y, closest);
				Trace hit = primitives[i].hit(clipped);
				if (hit.hit && hit.distance < closest) closest = hit.distance;
				ret = Trace::min(ret, hit);
			}
			continue;
		}

		Vec2 l_times = ray.dist_bounds;
		Vec2 r_times = ray.dist_bounds;
		l_times.y = std::min(l_times.y, closest);
		r_times.y = std::min(r_times.y, closest);
		bool l_hit = nodes[node.l].bbox.hit(ray, l_times);
		bool r_hit = nodes[node.r].bbox.hit(ray, r_times);

		if (l_hit && r_hit) {
			// Stack is LIFO, so push farther first to visit nearer child next.
			if (l_times.x <= r_times.x) {
				tstack.push(node.r);
				tstack.push(node.l);
			} else {
				tstack.push(node.l);
				tstack.push(node.r);
			}
		} else if (l_hit) {
			tstack.push(node.l);
		} else if (r_hit) {
			tstack.push(node.r);
		}
	}

	return ret;
}

template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	build(std::move(prims), max_leaf_size);
}

template<typename Primitive> std::vector<Primitive> BVH<Primitive>::destructure() {
	nodes.clear();
	return std::move(primitives);
}

template<typename Primitive>
template<typename P>
typename std::enable_if<std::is_copy_assignable_v<P>, BVH<P>>::type BVH<Primitive>::copy() const {
	BVH<Primitive> ret;
	ret.nodes = nodes;
	ret.primitives = primitives;
	ret.root_idx = root_idx;
	return ret;
}

template<typename Primitive> Vec3 BVH<Primitive>::sample(RNG &rng, Vec3 from) const {
	if (primitives.empty()) return {};
	int32_t n = rng.integer(0, static_cast<int32_t>(primitives.size()));
	return primitives[n].sample(rng, from);
}

template<typename Primitive>
float BVH<Primitive>::pdf(Ray ray, const Mat4& T, const Mat4& iT) const {
	if (primitives.empty()) return 0.0f;
	float ret = 0.0f;
	for (auto& prim : primitives) ret += prim.pdf(ray, T, iT);
	return ret / primitives.size();
}

template<typename Primitive> void BVH<Primitive>::clear() {
	nodes.clear();
	primitives.clear();
}

template<typename Primitive> bool BVH<Primitive>::Node::is_leaf() const {
	// A node is a leaf if l == r, since all interior nodes must have distinct children
	return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
	Node n;
	n.bbox = box;
	n.start = start;
	n.size = size;
	n.l = l;
	n.r = r;
	nodes.push_back(n);
	return nodes.size() - 1;
}
 
template<typename Primitive> BBox BVH<Primitive>::bbox() const {
	if(nodes.empty()) return BBox{Vec3{0.0f}, Vec3{0.0f}};
	return nodes[root_idx].bbox;
}

template<typename Primitive> size_t BVH<Primitive>::n_primitives() const {
	return primitives.size();
}

template<typename Primitive>
uint32_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, uint32_t level,
                                   const Mat4& trans) const {

	std::stack<std::pair<size_t, uint32_t>> tstack;
	tstack.push({root_idx, 0u});
	uint32_t max_level = 0u;

	if (nodes.empty()) return max_level;

	while (!tstack.empty()) {

		auto [idx, lvl] = tstack.top();
		max_level = std::max(max_level, lvl);
		const Node& node = nodes[idx];
		tstack.pop();

		Spectrum color = lvl == level ? Spectrum(1.0f, 0.0f, 0.0f) : Spectrum(1.0f);
		GL::Lines& add = lvl == level ? active : lines;

		BBox box = node.bbox;
		box.transform(trans);
		Vec3 min = box.min, max = box.max;

		auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

		edge(min, Vec3{max.x, min.y, min.z});
		edge(min, Vec3{min.x, max.y, min.z});
		edge(min, Vec3{min.x, min.y, max.z});
		edge(max, Vec3{min.x, max.y, max.z});
		edge(max, Vec3{max.x, min.y, max.z});
		edge(max, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

		if (!node.is_leaf()) {
			tstack.push({node.l, lvl + 1});
			tstack.push({node.r, lvl + 1});
		} else {
			for (size_t i = node.start; i < node.start + node.size; i++) {
				uint32_t c = primitives[i].visualize(lines, active, level - lvl, trans);
				max_level = std::max(c + lvl, max_level);
			}
		}
	}
	return max_level;
}

template class BVH<Triangle>;
template class BVH<Instance>;
template class BVH<Aggregate>;
template BVH<Triangle> BVH<Triangle>::copy<Triangle>() const;

} // namespace PT
