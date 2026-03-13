#include "test.h"
#include "geometry/halfedge.h"

#include <map>
#include <set>

static void expect_triangulate(Halfedge_Mesh& mesh) {
	size_t numEdges = mesh.edges.size();
	size_t numFaces = mesh.faces.size();

	// count the number of triangle edges/faces to be generated after triangulation
	size_t c = 0;
	for (Halfedge_Mesh::FaceRef f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		if (!f->boundary && f->degree() > 3) c += f->degree() - 3;
	}

	std::map<uint32_t, Vec3> verts;
	for (Halfedge_Mesh::VertexRef v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
		verts.insert({v->id, v->position});
	}
	std::set<uint32_t> edge_ids;
	for (Halfedge_Mesh::EdgeRef e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
		edge_ids.insert(e->id);
	}
	std::set<uint32_t> face_ids;
	for (Halfedge_Mesh::FaceRef f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		face_ids.insert(f->id);
	}

	//TODO: halfedges + twin relationships should all be preserved
	//TODO: shape of boundary faces should be preserved

	mesh.triangulate();

	if (auto msg = mesh.validate()) {
		throw Test::error("Invalid mesh: " + msg.value().second);
	}

	// check that all faces are degree 3
	for (Halfedge_Mesh::FaceRef f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		if (!f->boundary && f->degree() != 3) {
			throw Test::error("Triangulation created a non-triangular face!");
		}
	}

	// check for expected number of elements
	if (numEdges + c != mesh.edges.size()) {
		throw Test::error("Triangulation did not create the expected number of edges!");
	}
	if (numFaces + c != mesh.faces.size()) {
		throw Test::error("Triangulation did not create the expected number of faces!");
	}

	// check that original elements are preserved
	std::map<uint32_t, Vec3> new_verts;
	for (Halfedge_Mesh::VertexRef v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
		new_verts.insert({v->id, v->position});
	}
	std::set<uint32_t> new_edge_ids;
	for (Halfedge_Mesh::EdgeRef e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
		new_edge_ids.insert(e->id);
	}
	std::set<uint32_t> new_face_ids;
	for (Halfedge_Mesh::FaceRef f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		new_face_ids.insert(f->id);
	}

	if (verts != new_verts) {
		throw Test::error("Triangulation should preserve original vertices!");
	}
	std::set<uint32_t> diff;
	std::set_difference(edge_ids.begin(), edge_ids.end(), new_edge_ids.begin(), new_edge_ids.end(),
	                    std::inserter(diff, diff.end()));
	if (diff.size() != 0) {
		throw Test::error("Triangulation should preserve original edges!");
	}
	std::set_difference(face_ids.begin(), face_ids.end(), new_face_ids.begin(), new_face_ids.end(),
	                    std::inserter(diff, diff.end()));
	if (diff.size() != 0) {
		throw Test::error("Triangulation should preserve original faces!");
	}
}


/*
BASIC CASE

Triangulates a square
*/
Test test_a2_g1_triangulate_basic_square("a2.g1.triangulate.basic.square", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3{-0.5f, 0.0f,-0.5f}, Vec3{-0.5f, 0.0f, 0.5f},
		Vec3{ 0.5f, 0.0f,-0.5f}, Vec3{ 0.5f, 0.0f, 0.5f}
	},{
		{1, 3, 2, 0}
	});

	// Many different implementations of triangulating, so just checks that all the faces are triangles and some other misc things
	expect_triangulate(mesh);
});

/*
BASIC CASE

Triangulates a cube with square faces
*/
Test test_a2_g1_triangulate_basic_quad_cube("a2.g1.triangulate.basic.quad_cube", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3{-1.0f, 1.0f, 1.0f}, 	Vec3{-1.0f, 1.0f, -1.0f},
		Vec3{-1.0f, -1.0f, -1.0f}, 	Vec3{-1.0f, -1.0f, 1.0f},
		Vec3{1.0f, -1.0f, -1.0f}, 	Vec3{1.0f, -1.0f, 1.0f},
		Vec3{1.0f, 1.0f, -1.0f}, 	Vec3{1.0f, 1.0f, 1.0f}
	}, {
		{3, 0, 1, 2}, 
		{5, 3, 2, 4}, 
		{7, 5, 4, 6}, 
		{0, 7, 6, 1},
		{0, 3, 5, 7}, 
		{6, 4, 2, 1}
	});

	// Many different implementations of triangulating, so just checks that all the faces are triangles and some other misc things
	expect_triangulate(mesh);
});


/*
EDGE CASE: mesh that is already fully triangulated.
Triangulating a triangle mesh should be a no-op
(same element counts, same IDs, same positions).
*/
Test test_a2_g1_triangulate_already_triangle("a2.g1.triangulate.already_triangle", []() {
	// A single triangle
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3{0.0f, 0.0f, 0.0f},
		Vec3{1.0f, 0.0f, 0.0f},
		Vec3{0.0f, 1.0f, 0.0f}
	},{
		{0, 1, 2}
	});
	expect_triangulate(mesh);
});

/*
EDGE CASE: pentagon (5-sided face).
A pentagon needs exactly 2 diagonal cuts to become 3 triangles.
Expected: +2 edges, +2 faces.
*/
Test test_a2_g1_triangulate_pentagon("a2.g1.triangulate.pentagon", []() {
	// Regular-ish pentagon in the XZ plane
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3{ 0.0f,  0.0f,  1.0f},
		Vec3{ 0.951f, 0.0f,  0.309f},
		Vec3{ 0.588f, 0.0f, -0.809f},
		Vec3{-0.588f, 0.0f, -0.809f},
		Vec3{-0.951f, 0.0f,  0.309f}
	},{
		{0, 1, 2, 3, 4}
	});
	expect_triangulate(mesh);
});

/*
EDGE CASE: hexagon (6-sided face).
A hexagon needs exactly 3 diagonal cuts to become 4 triangles.
Expected: +3 edges, +3 faces.
*/
Test test_a2_g1_triangulate_hexagon("a2.g1.triangulate.hexagon", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3{ 1.0f, 0.0f,  0.0f},
		Vec3{ 0.5f, 0.0f,  0.866f},
		Vec3{-0.5f, 0.0f,  0.866f},
		Vec3{-1.0f, 0.0f,  0.0f},
		Vec3{-0.5f, 0.0f, -0.866f},
		Vec3{ 0.5f, 0.0f, -0.866f}
	},{
		{0, 1, 2, 3, 4, 5}
	});
	expect_triangulate(mesh);
});

/*
MIXED CASE: mesh with both triangles and quads.
Only the quads should be split; triangles should be untouched.
This is a simple shape: one triangle face and one quad face sharing an edge.
*/
Test test_a2_g1_triangulate_mixed_tri_quad("a2.g1.triangulate.mixed_tri_quad", []() {
	// A "house" shape: a square base with a triangular roof peak sharing the top edge
	//   v4
	//  /  \
	// v2 - v3
	// |    |
	// v0 - v1
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3{-1.0f, 0.0f, 0.0f}, // v0
		Vec3{ 1.0f, 0.0f, 0.0f}, // v1
		Vec3{-1.0f, 1.0f, 0.0f}, // v2
		Vec3{ 1.0f, 1.0f, 0.0f}, // v3
		Vec3{ 0.0f, 2.0f, 0.0f}  // v4
	},{
		{0, 1, 3, 2}, // quad (the square base)
		{2, 3, 4}     // triangle (the roof)
	});
	expect_triangulate(mesh);
});

/*
MIXED CASE: mesh with multiple quads sharing vertices.
Two quads sharing a single vertex (like an "X" seen from above).
Both should be split independently.
*/
Test test_a2_g1_triangulate_two_quads("a2.g1.triangulate.two_quads", []() {
	// Two separate quad faces that share one vertex (v4 in center)
	// but are not connected by an edge (they share a vertex only).
	// Use a flat mesh with two separate quad patches:
	//  v0-v1   v5-v6
	//  |  |    |  |
	//  v2-v3   v7-v8
	// (Two disconnected quads -- still a valid mesh with a single boundary loop)
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3{-2.0f, 0.0f,  1.0f}, // v0
		Vec3{-1.0f, 0.0f,  1.0f}, // v1
		Vec3{-2.0f, 0.0f,  0.0f}, // v2
		Vec3{-1.0f, 0.0f,  0.0f}, // v3
		Vec3{ 1.0f, 0.0f,  1.0f}, // v4
		Vec3{ 2.0f, 0.0f,  1.0f}, // v5
		Vec3{ 1.0f, 0.0f,  0.0f}, // v6
		Vec3{ 2.0f, 0.0f,  0.0f}  // v7
	},{
		{0, 1, 3, 2},
		{4, 5, 7, 6}
	});
	expect_triangulate(mesh);
});

/*
BOUNDARY CASE: a quad with a boundary face.
The quad should be triangulated; the boundary face should be untouched.
Boundary faces must never be triangulated.
*/
Test test_a2_g1_triangulate_boundary_untouched("a2.g1.triangulate.boundary_untouched", []() {
	// A single quad (open mesh, so it has a boundary loop around the outside)
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3{0.0f, 0.0f, 0.0f},
		Vec3{1.0f, 0.0f, 0.0f},
		Vec3{1.0f, 0.0f, 1.0f},
		Vec3{0.0f, 0.0f, 1.0f}
	},{
		{0, 1, 2, 3}
	});

	// Record how many boundary faces exist before
	size_t boundary_faces_before = 0;
	for (auto f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		if (f->boundary) boundary_faces_before++;
	}

	expect_triangulate(mesh);

	// Boundary face count must not change
	size_t boundary_faces_after = 0;
	for (auto f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		if (f->boundary) boundary_faces_after++;
	}
	if (boundary_faces_before != boundary_faces_after) {
		throw Test::error("Triangulation should not change the number of boundary faces!");
	}
});
