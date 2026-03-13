#include "test.h"
#include "geometry/halfedge.h"

#include <set>

static void expect_cc(Halfedge_Mesh& mesh, Halfedge_Mesh const &after) {
	size_t numVerts = mesh.vertices.size();
	size_t numEdges = mesh.edges.size();
	size_t numFaces = mesh.faces.size();

	// count the number of new edges to be generated after catmull-clark subdiv
	// count the number of quad faces to be generated after catmull-clark subdiv
	size_t ce = 0;
	size_t cf = 0;
	for (Halfedge_Mesh::FaceRef f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		if (!f->boundary) {
			ce += f->degree();
			cf += f->degree() - 1;
		}
	}

	mesh.catmark_subdivide();

	if (auto msg = mesh.validate()) {
		throw Test::error("Invalid mesh: " + msg.value().second);
	}

	// check that all faces are degree 4
	for (Halfedge_Mesh::FaceRef f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
		if (!f->boundary && f->degree() != 4) {
			throw Test::error("Catmull-clark subdivision created a non-quad face!");
		}
	}

	// check for expected number of elements
	if (numVerts + numEdges + numFaces - mesh.n_boundaries() != mesh.vertices.size()) {
		throw Test::error("Catmull-clark subdivision did not create the expected number of vertices!");
	}
	if (numEdges * 2 + ce != mesh.edges.size()) {
		throw Test::error("Catmull-clark subdivision did not create the expected number of edges!");
	}
	if (numFaces + cf != mesh.faces.size()) {
		throw Test::error("Catmull-clark subdivision did not create the expected number of faces!");
	}

	// check mesh shape:
	if (auto diff = Test::differs(mesh, after)) {
		throw Test::error("Result does not match expected: " + diff.value());
	}
}

/*
EDGE CASE

Catmull-Clark subdivides a square
*/
Test test_a2_g3_catmull_clark_edge_square("a2.g3.catmull_clark.edge.square", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({ 
		Vec3{-1.0f, 1.0f, 0.0f}, Vec3{ 1.0f, 1.0f, 0.0f},
		Vec3{-1.0f,-1.0f, 0.0f}, Vec3{ 1.0f,-1.0f, 0.0f} 
	},{ 
		{2, 3, 1, 0} 
	});

	Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({ 
		Vec3{-0.75f, 0.75f, 0.0f}, 	Vec3{ 0.0f, 1.0f, 0.0f }, 
		Vec3{ 0.75f, 0.75f, 0.0f}, 	Vec3{-1.0f, 0.0f, 0.0f}, 
		Vec3{ 0.0f, 0.0f, 0.0f }, 	Vec3{ 1.0f, 0.0f, 0.0f},
		Vec3{-0.75f, -0.75f, 0.0f}, Vec3{ 0.0f,-1.0f, 0.0f }, 
		Vec3{ 0.75f, -0.75f, 0.0f} 
	},{ 
		{3, 4, 1, 0}, 
		{4, 5, 2, 1},  
		{6, 7, 4, 3}, 
		{7, 8, 5, 4} 
	});

	expect_cc(mesh, after);
});

/*
BASIC CASE

Catmull-Clark subdivides a cube with square faces
*/
Test test_a2_g3_catmull_clark_basic_quad_cube("a2.g3.catmull_clark.basic.quad_cube", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({ 
		Vec3{-1.0f, 1.0f, 1.0f}, 	Vec3{-1.0f, 1.0f, -1.0f},
		Vec3{-1.0f, -1.0f, -1.0f}, 	Vec3{-1.0f, -1.0f, 1.0f},
		Vec3{1.0f, -1.0f, -1.0f}, 	Vec3{1.0f, -1.0f, 1.0f},
		Vec3{1.0f, 1.0f, -1.0f}, 	Vec3{1.0f, 1.0f, 1.0f}
	},{ 
		{3, 0, 1, 2}, 
		{5, 3, 2, 4}, 
		{7, 5, 4, 6}, 
		{0, 7, 6, 1}, 
		{0, 3, 5, 7}, 
		{6, 4, 2, 1}
	});

	Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({ 
		Vec3{-1.0f, 0.0f, 0.0f}, 					Vec3{-0.75f, 0.0f, 0.75f},
		Vec3{-0.555556f, 0.555556f, 0.555556f}, 	Vec3{-0.75f, 0.75f, 0.0f},
		Vec3{-0.555556f, 0.555556f, -0.555556f}, 	Vec3{-0.75f, 0.0f, -0.75f},
		Vec3{-0.555556f, -0.555556f, -0.555556f}, 	Vec3{-0.75f, -0.75f, 0.0f},
		Vec3{-0.555556f, -0.555556f, 0.555556f}, 	Vec3{0.0f, -1.0f, 0.0f},
		Vec3{0.0f, -0.75f, 0.75f}, 					Vec3{0.0f, -0.75f, -0.75f},
		Vec3{0.555556f, -0.555556f, -0.555556f}, 	Vec3{0.75f, -0.75f, 0.0f},
		Vec3{0.555556f, -0.555556f, 0.555556f}, 	Vec3{1.0f, 0.0f, 0.0f},
		Vec3{0.75f, 0.0f, 0.75f}, 					Vec3{0.75f, 0.0f, -0.75f},
		Vec3{0.555556f, 0.555556f, -0.555556f}, 	Vec3{0.75f, 0.75f, 0.0f},
		Vec3{0.555556f, 0.555556f, 0.555556f}, 		Vec3{0.0f, 1.0f, 0.0f},
		Vec3{0.0f, 0.75f, 0.75f}, 					Vec3{0.0f, 0.75f, -0.75f},
		Vec3{0.0f, 0.0f, 1.0f}, 					Vec3{0.0f, 0.0f, -1.0f}
	},{ 
		{3, 0, 1, 2},     {5, 0, 3, 4},     {7, 0, 5, 6},     {1, 0, 7, 8},     {7, 9, 10, 8},    {11, 9, 7, 6},
		{13, 9, 11, 12},  {10, 9, 13, 14},  {13, 15, 16, 14}, {17, 15, 13, 12}, {19, 15, 17, 18}, {16, 15, 19, 20},
		{19, 21, 22, 20}, {23, 21, 19, 18}, {3, 21, 23, 4},   {22, 21, 3, 2},   {10, 24, 1, 8},   {16, 24, 10, 14},
		{22, 24, 16, 20}, {1, 24, 22, 2},   {11, 25, 17, 12}, {5, 25, 11, 6},   {23, 25, 5, 4},   {17, 25, 23, 18}
	});

	expect_cc(mesh, after);
});

/*
EDGE CASE

Catmull-Clark subdivides a single triangle.
Tests boundary vertex and boundary edge rules with a 3-sided face.

Input: triangle with vertices A(1,0,0), B(-1,0,0), C(0,1,0).
All 3 edges are boundary. All 3 vertices are boundary.

Boundary edge rule:   E = midpoint(v0, v1)
Boundary vertex rule: V = (1/8)*bn1 + (6/8)*v + (1/8)*bn2

Face point (centroid): (0, 1/3, 0)
Edge points: E_AB=(0,0,0), E_BC=(-0.5,0.5,0), E_CA=(0.5,0.5,0)
Vertex points: VA=(0.625,0.125,0), VB=(-0.625,0.125,0), VC=(0,0.75,0)

Result: 3 quads, one per original vertex.
*/
Test test_a2_g3_catmull_clark_edge_triangle("a2.g3.catmull_clark.edge.triangle", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3{ 1.0f, 0.0f, 0.0f},   // 0: A
		Vec3{-1.0f, 0.0f, 0.0f},   // 1: B
		Vec3{ 0.0f, 1.0f, 0.0f}    // 2: C
	},{
		{0, 1, 2}
	});

	// New vertices (7 total after subdividing 1 triangle):
	// [0] VA   = (0.625,    0.125, 0)
	// [1] VB   = (-0.625,   0.125, 0)
	// [2] VC   = (0,        0.75,  0)
	// [3] E_AB = (0,        0,     0)  <- boundary edge midpoint A-B
	// [4] E_BC = (-0.5,     0.5,   0)  <- boundary edge midpoint B-C
	// [5] E_CA = (0.5,      0.5,   0)  <- boundary edge midpoint C-A
	// [6] F    = (0,        1/3,   0)  <- face centroid
	Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
		Vec3{ 0.625f,     0.125f,    0.0f},   // 0: VA
		Vec3{-0.625f,     0.125f,    0.0f},   // 1: VB
		Vec3{ 0.0f,       0.75f,     0.0f},   // 2: VC
		Vec3{ 0.0f,       0.0f,      0.0f},   // 3: E_AB
		Vec3{-0.5f,       0.5f,      0.0f},   // 4: E_BC
		Vec3{ 0.5f,       0.5f,      0.0f},   // 5: E_CA
		Vec3{ 0.0f,       0.333333f, 0.0f}    // 6: F
	},{
		{0, 5, 6, 3},   // quad around A: VA, E_CA, F, E_AB
		{1, 3, 6, 4},   // quad around B: VB, E_AB, F, E_BC
		{2, 4, 6, 5}    // quad around C: VC, E_BC, F, E_CA
	});

	expect_cc(mesh, after);
});

/*
BASIC CASE

Catmull-Clark subdivides two triangles sharing an interior edge (a diamond).
Tests that the interior edge rule uses both adjacent face points, while
boundary edges use simple midpoints, and all 4 boundary vertices use the crease rule.

Vertices: A(0,1,0), B(-1,0,0), C(0,-1,0), D(1,0,0)
Faces: {A,B,C} and {A,C,D}  -- shared interior edge A-C

Interior edge A-C: (A+C+F_ABC+F_ACD)/4 = (0,0,0)
Boundary vertices A,B,C,D all use the crease rule with their 2 boundary neighbors each.
*/
Test test_a2_g3_catmull_clark_basic_two_triangles("a2.g3.catmull_clark.basic.two_triangles", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3{ 0.0f,  1.0f, 0.0f},   // 0: A
		Vec3{-1.0f,  0.0f, 0.0f},   // 1: B
		Vec3{ 0.0f, -1.0f, 0.0f},   // 2: C
		Vec3{ 1.0f,  0.0f, 0.0f}    // 3: D
	},{
		{0, 1, 2},   // face ABC
		{0, 2, 3}    // face ACD
	});

	// Face points:
	//   F_ABC = (-1/3, 0, 0)   [9]
	//   F_ACD = ( 1/3, 0, 0)   [10]
	//
	// Edge points:
	//   [4] E_AC (interior, faces ABC+ACD): (A+C+F_ABC+F_ACD)/4 = (0,0,0)
	//   [5] E_AB (boundary): (-0.5,  0.5, 0)
	//   [6] E_BC (boundary): (-0.5, -0.5, 0)
	//   [7] E_AD (boundary): ( 0.5,  0.5, 0)
	//   [8] E_CD (boundary): ( 0.5, -0.5, 0)
	//
	// Vertex points (all boundary, crease rule):
	//   [0] VA: boundary nbrs B,D -> (0,    0.75, 0)
	//   [1] VB: boundary nbrs A,C -> (-0.75, 0,   0)
	//   [2] VC: boundary nbrs B,D -> (0,   -0.75, 0)
	//   [3] VD: boundary nbrs A,C -> (0.75,  0,   0)
	Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
		Vec3{ 0.0f,      0.75f,  0.0f},   // 0: VA
		Vec3{-0.75f,     0.0f,   0.0f},   // 1: VB
		Vec3{ 0.0f,     -0.75f,  0.0f},   // 2: VC
		Vec3{ 0.75f,     0.0f,   0.0f},   // 3: VD
		Vec3{ 0.0f,      0.0f,   0.0f},   // 4: E_AC (interior)
		Vec3{-0.5f,      0.5f,   0.0f},   // 5: E_AB
		Vec3{-0.5f,     -0.5f,   0.0f},   // 6: E_BC
		Vec3{ 0.5f,      0.5f,   0.0f},   // 7: E_AD
		Vec3{ 0.5f,     -0.5f,   0.0f},   // 8: E_CD
		Vec3{-0.333333f, 0.0f,   0.0f},   // 9:  F_ABC
		Vec3{ 0.333333f, 0.0f,   0.0f}    // 10: F_ACD
	},{
		// 3 quads from face ABC:
		{0, 5,  9, 4},   // around A: VA, E_AB, F_ABC, E_AC
		{1, 6,  9, 5},   // around B: VB, E_BC, F_ABC, E_AB
		{2, 4,  9, 6},   // around C: VC, E_AC, F_ABC, E_BC
		// 3 quads from face ACD:
		{0, 4, 10, 7},   // around A: VA, E_AC, F_ACD, E_AD
		{2, 8, 10, 4},   // around C: VC, E_CD, F_ACD, E_AC
		{3, 7, 10, 8}    // around D: VD, E_AD, F_ACD, E_CD
	});

	expect_cc(mesh, after);
});

/*
BASIC CASE

Catmull-Clark subdivides a 2x1 rectangle (non-square single quad).
Tests that non-uniform boundary edges use simple midpoints and that boundary
vertices with unequal neighbor distances are correctly weighted by the crease rule.

Vertices: A(0,0,0), B(2,0,0), C(2,1,0), D(0,1,0)   face {A,B,C,D}
All 4 edges are boundary (single non-boundary face).

Face point F = (1, 0.5, 0)
E_AB=(1,0,0), E_BC=(2,0.5,0), E_CD=(1,1,0), E_DA=(0,0.5,0)
VA=(0.25,0.125,0), VB=(1.75,0.125,0), VC=(1.75,0.875,0), VD=(0.25,0.875,0)
*/
Test test_a2_g3_catmull_clark_basic_rectangle("a2.g3.catmull_clark.basic.rectangle", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3{0.0f, 0.0f, 0.0f},   // 0: A
		Vec3{2.0f, 0.0f, 0.0f},   // 1: B
		Vec3{2.0f, 1.0f, 0.0f},   // 2: C
		Vec3{0.0f, 1.0f, 0.0f}    // 3: D
	},{
		{0, 1, 2, 3}
	});

	// [0] VA   = (0.25,  0.125, 0)
	// [1] VB   = (1.75,  0.125, 0)
	// [2] VC   = (1.75,  0.875, 0)
	// [3] VD   = (0.25,  0.875, 0)
	// [4] E_AB = (1.0,   0.0,   0)
	// [5] E_BC = (2.0,   0.5,   0)
	// [6] E_CD = (1.0,   1.0,   0)
	// [7] E_DA = (0.0,   0.5,   0)
	// [8] F    = (1.0,   0.5,   0)
	Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
		Vec3{0.25f,  0.125f, 0.0f},   // 0: VA
		Vec3{1.75f,  0.125f, 0.0f},   // 1: VB
		Vec3{1.75f,  0.875f, 0.0f},   // 2: VC
		Vec3{0.25f,  0.875f, 0.0f},   // 3: VD
		Vec3{1.0f,   0.0f,   0.0f},   // 4: E_AB
		Vec3{2.0f,   0.5f,   0.0f},   // 5: E_BC
		Vec3{1.0f,   1.0f,   0.0f},   // 6: E_CD
		Vec3{0.0f,   0.5f,   0.0f},   // 7: E_DA
		Vec3{1.0f,   0.5f,   0.0f}    // 8: F
	},{
		{0, 7, 8, 4},   // quad around A: VA, E_DA, F, E_AB
		{1, 4, 8, 5},   // quad around B: VB, E_AB, F, E_BC
		{2, 5, 8, 6},   // quad around C: VC, E_BC, F, E_CD
		{3, 6, 8, 7}    // quad around D: VD, E_CD, F, E_DA
	});

	expect_cc(mesh, after);
});