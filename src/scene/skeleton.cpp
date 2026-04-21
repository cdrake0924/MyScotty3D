#include <unordered_set>
#include "skeleton.h"
#include "test.h"
#include <iostream>

void Skeleton::Bone::compute_rotation_axes(Vec3 *x_, Vec3 *y_, Vec3 *z_) const {
	assert(x_ && y_ && z_);
	auto &x = *x_;
	auto &y = *y_;
	auto &z = *z_;

	y = extent.unit();
	if (!y.valid()) {
		y = Vec3{0.0f, 1.0f, 0.0f};
	}

	x = Vec3{1.0f, 0.0f, 0.0f};
	x = (x - dot(x,y) * y).unit();
	if (!x.valid()) {
		x = Vec3{0.0f, 0.0f, 1.0f};
		x = (x - dot(x,y) * y).unit();
	}

	z = cross(x,y);

	float cr = std::cos(roll / 180.0f * PI_F);
	float sr = std::sin(roll / 180.0f * PI_F);
	std::tie(x, z) = std::make_pair(cr * x + sr * -z, cr * z + sr * x);
}

std::vector< Mat4 > Skeleton::bind_pose() const {
	std::vector< Mat4 > bind;
	bind.reserve(bones.size());

	for (size_t i = 0; i < bones.size(); ++i) {
		Bone const &bone = bones[i];
		Mat4 local_to_parent;
		if (bone.parent == -1U) {
			local_to_parent = Mat4::translate(base);
			bind.emplace_back(local_to_parent);
		} else {
			local_to_parent = Mat4::translate(bones[bone.parent].extent);
			bind.emplace_back(bind[bone.parent] * local_to_parent);
		}
	}

	assert(bind.size() == bones.size());
	return bind;
}

std::vector< Mat4 > Skeleton::current_pose() const {
	std::vector< Mat4 > pose;
	pose.reserve(bones.size());

	for (size_t i = 0; i < bones.size(); ++i) {
		Bone const &bone = bones[i];

		Vec3 x_axis, y_axis, z_axis;
		bone.compute_rotation_axes(&x_axis, &y_axis, &z_axis);

		Mat4 R = Mat4::angle_axis(bone.pose.z, z_axis)
		       * Mat4::angle_axis(bone.pose.y, y_axis)
		       * Mat4::angle_axis(bone.pose.x, x_axis);

		if (bone.parent == -1U) {
			pose.emplace_back(Mat4::translate(base + base_offset) * R);
		} else {
			Mat4 local_to_parent = Mat4::translate(bones[bone.parent].extent) * R;
			pose.emplace_back(pose[bone.parent] * local_to_parent);
		}
	}

	return pose;
}

std::vector< Vec3 > Skeleton::gradient_in_current_pose() const {
	std::vector< Vec3 > gradient(bones.size(), Vec3{0.0f, 0.0f, 0.0f});

	std::vector< Mat4 > pose = current_pose();

	for (auto const &handle : handles) {
		if (!handle.enabled) continue;
		if (handle.bone >= bones.size()) continue;

		BoneIndex tip_bone = handle.bone;
		Vec3 tip = pose[tip_bone] * bones[tip_bone].extent;
		Vec3 diff = tip - handle.target;

		BoneIndex b = tip_bone;
		while (b != -1U) {
			Bone const &bone = bones[b];

			Vec3 x_axis_local, y_axis_local, z_axis_local;
			bone.compute_rotation_axes(&x_axis_local, &y_axis_local, &z_axis_local);

			Mat4 base_xform;
			if (bone.parent == -1U) {
				base_xform = Mat4::translate(base + base_offset);
			} else {
				base_xform = pose[bone.parent] * Mat4::translate(bones[bone.parent].extent);
			}

			Vec3 r = base_xform * Vec3{0.0f, 0.0f, 0.0f};

			Mat4 Rz = Mat4::angle_axis(bone.pose.z, z_axis_local);
			Mat4 Ry = Mat4::angle_axis(bone.pose.y, y_axis_local);

			Vec3 axis_z = base_xform.rotate(z_axis_local);
			Vec3 axis_y = (base_xform * Rz).rotate(y_axis_local);
			Vec3 axis_x = (base_xform * Rz * Ry).rotate(x_axis_local);

			Vec3 d = tip - r;
			Vec3 dp_drx = cross(axis_x, d);
			Vec3 dp_dry = cross(axis_y, d);
			Vec3 dp_drz = cross(axis_z, d);

			gradient[b].x += dot(diff, dp_drx);
			gradient[b].y += dot(diff, dp_dry);
			gradient[b].z += dot(diff, dp_drz);

			b = bone.parent;
		}
	}

	assert(gradient.size() == bones.size());
	return gradient;
}

bool Skeleton::solve_ik(uint32_t steps) {
	const float tau = 1.0f;
	const float eps_sq = 1e-8f;

	for (uint32_t s = 0; s < steps; ++s) {
		std::vector< Vec3 > grads = gradient_in_current_pose();

		float mag_sq = 0.0f;
		for (auto const &g : grads) mag_sq += g.norm_squared();
		if (mag_sq < eps_sq) return true;

		for (size_t i = 0; i < bones.size(); ++i) {
			bones[i].pose -= tau * grads[i];
		}
	}

	return false;
}

Vec3 Skeleton::closest_point_on_line_segment(Vec3 const &a, Vec3 const &b, Vec3 const &p) {
	Vec3 ab = b - a;
	float len_sq = dot(ab, ab);
	if (len_sq < EPS_F * EPS_F) return a;
	float t = dot(p - a, ab) / len_sq;
	if (t < 0.0f) t = 0.0f;
	if (t > 1.0f) t = 1.0f;
	return a + t * ab;
}

void Skeleton::assign_bone_weights(Halfedge_Mesh *mesh_) const {
	assert(mesh_);
	auto &mesh = *mesh_;

	std::vector< Mat4 > bind = bind_pose();
	std::vector< std::pair< Vec3, Vec3 > > segments;
	segments.reserve(bones.size());
	for (size_t b = 0; b < bones.size(); ++b) {
		Vec3 start = bind[b] * Vec3{0.0f, 0.0f, 0.0f};
		Vec3 end = bind[b] * bones[b].extent;
		segments.emplace_back(start, end);
	}

	for (auto &v : mesh.vertices) {
		v.bone_weights.clear();

		std::vector< std::pair< uint32_t, float > > weights;
		float total = 0.0f;

		for (size_t b = 0; b < bones.size(); ++b) {
			float r = bones[b].radius;
			if (r <= 0.0f) continue;
			Vec3 closest = closest_point_on_line_segment(segments[b].first, segments[b].second, v.position);
			float d = (v.position - closest).norm();
			float w = std::max(0.0f, r - d) / r;
			if (w > 0.0f) {
				weights.emplace_back(uint32_t(b), w);
				total += w;
			}
		}

		if (total > 0.0f) {
			for (auto const &p : weights) {
				Halfedge_Mesh::Vertex::Bone_Weight bw;
				bw.bone = p.first;
				bw.weight = p.second / total;
				v.bone_weights.emplace_back(bw);
			}
		}
	}
}

Indexed_Mesh Skeleton::skin(Halfedge_Mesh const &mesh, std::vector< Mat4 > const &bind, std::vector< Mat4 > const &current) {
	assert(bind.size() == current.size());

	std::vector< Mat4 > bind_to_current;
	bind_to_current.reserve(bind.size());
	std::vector< Mat4 > bind_to_current_normal;
	bind_to_current_normal.reserve(bind.size());
	for (size_t b = 0; b < bind.size(); ++b) {
		Mat4 T = current[b] * bind[b].inverse();
		bind_to_current.emplace_back(T);
		bind_to_current_normal.emplace_back(Mat4::transpose(T.inverse()));
	}

	std::unordered_map< Halfedge_Mesh::VertexCRef, Vec3 > skinned_positions;
	std::unordered_map< Halfedge_Mesh::HalfedgeCRef, Vec3 > skinned_normals;
	skinned_positions.reserve(mesh.vertices.size());
	skinned_normals.reserve(mesh.halfedges.size());

	for (auto vi = mesh.vertices.begin(); vi != mesh.vertices.end(); ++vi) {
		Vec3 skinned_pos;
		Mat4 normal_xform = Mat4::I;
		bool has_weights = !vi->bone_weights.empty();

		if (has_weights) {
			Mat4 weighted = Mat4::Zero;
			Mat4 weighted_normal = Mat4::Zero;
			for (auto const &bw : vi->bone_weights) {
				if (bw.bone >= bind_to_current.size()) continue;
				weighted += bind_to_current[bw.bone] * bw.weight;
				weighted_normal += bind_to_current_normal[bw.bone] * bw.weight;
			}
			skinned_pos = weighted * vi->position;
			normal_xform = weighted_normal;
		} else {
			skinned_pos = vi->position;
		}

		skinned_positions.emplace(vi, skinned_pos);

		auto h = vi->halfedge;
		do {
			Vec3 n;
			if (has_weights) {
				n = normal_xform.rotate(h->corner_normal);
				if (n.norm_squared() > 0.0f) n = n.unit();
				else n = h->corner_normal;
			} else {
				n = h->corner_normal;
			}
			skinned_normals.emplace(h, n);
			h = h->twin->next;
		} while (h != vi->halfedge);
	}

	std::vector< Indexed_Mesh::Vert > verts;
	std::vector< Indexed_Mesh::Index > idxs;

	for (auto f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
		if (f->boundary) continue;

		uint32_t corners_begin = static_cast<uint32_t>(verts.size());
		auto h = f->halfedge;
		do {
			Indexed_Mesh::Vert vert;
			vert.pos = skinned_positions.at(h->vertex);
			vert.norm = skinned_normals.at(h);
			vert.uv = h->corner_uv;
			vert.id = f->id;
			verts.emplace_back(vert);
			h = h->next;
		} while (h != f->halfedge);
		uint32_t corners_end = static_cast<uint32_t>(verts.size());

		for (size_t i = corners_begin + 1; i + 1 < corners_end; ++i) {
			idxs.emplace_back(corners_begin);
			idxs.emplace_back(static_cast<uint32_t>(i));
			idxs.emplace_back(static_cast<uint32_t>(i + 1));
		}
	}

	return Indexed_Mesh(std::move(verts), std::move(idxs));
}

void Skeleton::for_bones(const std::function<void(Bone&)>& f) {
	for (auto& bone : bones) {
		f(bone);
	}
}


void Skeleton::erase_bone(BoneIndex bone) {
	assert(bone < bones.size());
	for (uint32_t b = 0; b < bones.size(); ++b) {
		if (bones[b].parent == -1U) continue;
		if (bones[b].parent == bone) {
			assert(b > bone);
			bones[b].extent += bones[bone].extent;
			bones[b].parent = bones[bone].parent;
		} else if (bones[b].parent > bone) {
			assert(b > bones[b].parent);
			bones[b].parent -= 1;
		}
	}
	bones.erase(bones.begin() + bone);
	for (uint32_t h = 0; h < handles.size(); ) {
		if (handles[h].bone == bone) {
			erase_handle(h);
		} else if (handles[h].bone > bone) {
			handles[h].bone -= 1;
			++h;
		} else {
			++h;
		}
	}
}

void Skeleton::erase_handle(HandleIndex handle) {
	assert(handle < handles.size());

	handles.erase(handles.begin() + handle);
}


Skeleton::BoneIndex Skeleton::add_bone(BoneIndex parent, Vec3 extent) {
	assert(parent == -1U || parent < bones.size());
	Bone bone;
	bone.extent = extent;
	bone.parent = parent;

	std::unordered_set< uint32_t > used;
	for (auto const &b : bones) {
		used.emplace(b.channel_id);
	}
	while (used.count(next_bone_channel_id)) ++next_bone_channel_id;
	bone.channel_id = next_bone_channel_id++;

	BoneIndex index = BoneIndex(bones.size());
	bones.emplace_back(bone);

	return index;
}

Skeleton::HandleIndex Skeleton::add_handle(BoneIndex bone, Vec3 target) {
	assert(bone < bones.size());
	Handle handle;
	handle.bone = bone;
	handle.target = target;

	std::unordered_set< uint32_t > used;
	for (auto const &h : handles) {
		used.emplace(h.channel_id);
	}
	while (used.count(next_handle_channel_id)) ++next_handle_channel_id;
	handle.channel_id = next_handle_channel_id++;

	HandleIndex index = HandleIndex(handles.size());
	handles.emplace_back(handle);

	return index;
}


Skeleton Skeleton::copy() {
	return *this;
}

void Skeleton::make_valid() {
	for (uint32_t b = 0; b < bones.size(); ++b) {
		if (!(bones[b].parent == -1U || bones[b].parent < b)) {
			warn("bones[%u].parent is %u, which is not < %u; setting to -1.", b, bones[b].parent, b);
			bones[b].parent = -1U;
		}
	}
	if (bones.empty() && !handles.empty()) {
		warn("Have %u handles but no bones. Deleting handles.", uint32_t(handles.size()));
		handles.clear();
	}
	for (uint32_t h = 0; h < handles.size(); ++h) {
		if (handles[h].bone >= HandleIndex(bones.size())) {
			warn("handles[%u].bone is %u, which is not < bones.size(); setting to 0.", h, handles[h].bone);
			handles[h].bone = 0;
		}
	}
}

Indexed_Mesh Skinned_Mesh::bind_mesh() const {
	return Indexed_Mesh::from_halfedge_mesh(mesh, Indexed_Mesh::SplitEdges);
}

Indexed_Mesh Skinned_Mesh::posed_mesh() const {
	return Skeleton::skin(mesh, skeleton.bind_pose(), skeleton.current_pose());
}

Skinned_Mesh Skinned_Mesh::copy() {
	return Skinned_Mesh{mesh.copy(), skeleton.copy()};
}
