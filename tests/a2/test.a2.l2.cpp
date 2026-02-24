#include "test.h"
#include "geometry/halfedge.h"

static void expect_split(Halfedge_Mesh &mesh, Halfedge_Mesh::EdgeRef edge, Halfedge_Mesh const &after) {
	if (auto ret = mesh.split_edge(edge)) {
		if (auto msg = mesh.validate()) {
			throw Test::error("Invalid mesh: " + msg.value().second);
		}
		// check mesh shape:
		if (auto difference = Test::differs(mesh, after, Test::CheckAllBits)) {
			throw Test::error("Resulting mesh did not match expected: " + *difference);
		}
	} else {
		throw Test::error("split_edge rejected operation!");
	}
}

/*
BASIC CASE:

Initial mesh:
0--1\
|  | \
|  |  2
|  | /
3--4/

Split Edge on Edge: 1-4

After mesh:
0--1\
|\ | \
| \2--3
|  | /
4--5/
*/
Test test_a2_l2_split_edge_basic_simple("a2.l2.split_edge.basic.simple", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f),
		                                            Vec3(2.2f, 0.0f, 0.0f),
		Vec3(-1.3f,-0.7f, 0.0f), Vec3(1.4f, -1.0f, 0.0f)
	}, {
		{0, 3, 4, 1}, 
		{1, 4, 2}
	});
	Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->next->edge;

	Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
		Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f),
		                         Vec3(1.25f, 0.0f, 0.0f),  Vec3(2.2f, 0.0f, 0.0f),
		Vec3(-1.3f,-0.7f, 0.0f), Vec3(1.4f, -1.0f, 0.0f)
	}, {
		{0, 4, 5, 2}, 
		{0, 2, 1}, 
		{1, 2, 3}, 
		{2, 5, 3}
	});

	expect_split(mesh, edge, after);
});

/*
EDGE CASE: 

Initial mesh:
0--1\
|  | \
|  |  2
|  | /
3--4/

Split Edge on Edge: 0-1

After mesh:
0--1--2\
|  /  | \
| /   |  3
|/    | /
4-----5/
*/
Test test_a2_l2_split_edge_edge_boundary("a2.l2.split_edge.edge.boundary", []() {
	Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
		Vec3(-1.0f, 1.1f, 0.0f), Vec3(1.1f, 1.0f, 0.0f),
		                                            Vec3(2.2f, 0.0f, 0.0f),
		Vec3(-1.3f,-0.7f, 0.0f), Vec3(1.4f, -1.0f, 0.0f)
	}, {
		{0, 3, 4, 1}, 
		{1, 4, 2}
	});
	Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->next->next->edge;

	Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
		Vec3(-1.0f, 1.1f, 0.0f),  Vec3(0.05f, 1.05f, 0.0f), Vec3(1.1f, 1.0f, 0.0f),
		                                            						Vec3(2.2f, 0.0f, 0.0f),
		Vec3(-1.3f,-0.7f, 0.0f), 							Vec3(1.4f, -1.0f, 0.0f)
	}, {
		{0, 4, 1}, 
		{1, 4, 5, 2}, 
		{2, 5, 3}
	});

	expect_split(mesh, edge, after);
});

/*

CASE 1: TWO TRIANGLES (Interior Edge)

Initial mesh (Square made of 2 triangles):

0-------3

| \ |

| \ |

| \ |

1-------2

Split Edge on Edge: 0-2 (Diagonal)

After mesh (Square made of 4 triangles):

0-------4

| \ / |

| \2/ |

| / \ |

1-------3

*/

Test test_a2_l2_split_edge_two_triangles("a2.l2.split_edge.two_triangles", []() {

Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({

Vec3(-1.0f, 1.0f, 0.0f), Vec3(-1.0f, -1.0f, 0.0f),

Vec3(1.0f, -1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f)

}, {

{0, 1, 2},

{0, 2, 3}

});

Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->next->edge;

Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({

Vec3(-1.0f, 1.0f, 0.0f), Vec3(-1.0f, -1.0f, 0.0f),

Vec3(0.0f, 0.0f, 0.0f), // added vertex 2

Vec3(1.0f, -1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f)

}, {

{0, 1, 2},

{2, 1, 3},

{0, 2, 4},

{2, 3, 4}

});

expect_split(mesh, edge, after);

});

/*

CASE 2: SINGLE TRIANGLE (Boundary Edge)

Initial mesh:

0

|\

| \

| \

1---2

Split Edge on Edge: 1-2 (Bottom boundary)

After mesh:

0

|\ \

| \ \

| \ \

1---2---3

*/

Test test_a2_l2_split_edge_triangle_boundary("a2.l2.split_edge.triangle_boundary", []() {

Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({

Vec3(0.0f, 1.0f, 0.0f), Vec3(-1.0f, -1.0f, 0.0f), Vec3(1.0f, -1.0f, 0.0f)

}, {

{0, 1, 2}

});

Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->edge;

Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({

Vec3(0.0f, 1.0f, 0.0f), Vec3(-1.0f, -1.0f, 0.0f),

Vec3(0.0f, -1.0f, 0.0f), // added vertex 2

Vec3(1.0f, -1.0f, 0.0f)

}, {

{0, 1, 2},

{0, 2, 3}

});

expect_split(mesh, edge, after);

});

/*

CASE 3: TWO QUADS (Interior Edge)

Initial mesh:

0---1---2

| | |

| | |

| | |

3---4---5

Split Edge on Edge: 1-4

After mesh:

0---1---3

| \ | |

| 2 |

| | \ |

4---5---6

*/

Test test_a2_l2_split_edge_two_quads("a2.l2.split_edge.two_quads", []() {

Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({

Vec3(-1.0f, 1.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(1.0f, 1.0f, 0.0f),

Vec3(-1.0f, -1.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f), Vec3(1.0f, -1.0f, 0.0f)

}, {

{0, 3, 4, 1},

{1, 4, 5, 2}

});

Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->next->edge;

Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({

Vec3(-1.0f, 1.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f),

Vec3(0.0f, 0.0f, 0.0f), // added vertex 2

Vec3(1.0f, 1.0f, 0.0f),

Vec3(-1.0f, -1.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f), Vec3(1.0f, -1.0f, 0.0f)

}, {

{0, 4, 5, 2},

{0, 2, 1},

{2, 5, 6},

{1, 2, 6, 3}

});

expect_split(mesh, edge, after);

});

/*
CASE: HIGH DEGREE VERTEX (Star Center)

Initial mesh (4-triangle fan around center):

(-1,1)-------(1,1)
  |  \       / |
  |   \     /  |
  |    \   /   |
  |     (0,0)  |
  |    /   \   |
  |   /     \  |
  |  /       \ |
(-1,-1)-----(1,-1)

Faces: {0,1,4}, {0,4,3}, {0,3,2}, {0,2,1}
where 0=(0,0), 1=(-1,1), 2=(1,1), 3=(1,-1), 4=(-1,-1)

Split Edge on Edge: 0-2  (center to top-right)

After mesh: new vertex 5=(0.5,0.5) inserted on edge 0-2
Diagonals added:
  - In face {0,1,2}: v=5 connects to vertex ccw from 2 = vertex 1 (-1,1)
  - In face {0,2,3}: v=5 connects to vertex ccw from 2 = vertex 3 (1,-1)

(-1,1)-------5-------(1,1)
  |  \      /|\      /|
  |   \    / | \    / |
  |    \  /  |  \  /  |
  |    (0,0) |  (1,1) |
  |    /  \  |  /     |
  |   /    \ | /      |
  |  /      \|/       |
(-1,-1)-----(1,-1)

Resulting faces:
  {0, 1, 5}   -- center, top-left, new midpoint
  {5, 1, 2}   -- new midpoint, top-left, top-right   (diagonal 5->1)
  {0, 5, 3}   -- center, new midpoint, bottom-right
  {5, 2, 3}   -- new midpoint, top-right, bottom-right (diagonal 5->3)
  {0, 3, 4}   -- unchanged
  {0, 4, 1}   -- unchanged
*/
Test test_a2_l2_split_edge_star_center("a2.l2.split_edge.star_center", []() {
    Halfedge_Mesh mesh = Halfedge_Mesh::from_indexed_faces({
        Vec3( 0.0f,  0.0f, 0.0f), // 0 = center
        Vec3(-1.0f,  1.0f, 0.0f), // 1 = top-left
        Vec3( 1.0f,  1.0f, 0.0f), // 2 = top-right
        Vec3( 1.0f, -1.0f, 0.0f), // 3 = bottom-right
        Vec3(-1.0f, -1.0f, 0.0f)  // 4 = bottom-left
    }, {
        {0, 1, 2}, // top face
        {0, 2, 3}, // right face  <-- edge 0-2 is shared between these two
        {0, 3, 4}, // bottom face
        {0, 4, 1}  // left face
    });

    // Edge 0-2: halfedges.begin() is in face {0,1,2}
    // begin() = 0->1, ->next = 1->2, ->next->next = 2->0, ->edge = edge 0-2
    Halfedge_Mesh::EdgeRef edge = mesh.halfedges.begin()->next->next->edge;

    Halfedge_Mesh after = Halfedge_Mesh::from_indexed_faces({
        Vec3( 0.0f,  0.0f, 0.0f), // 0 = center
        Vec3(-1.0f,  1.0f, 0.0f), // 1 = top-left
        Vec3( 1.0f,  1.0f, 0.0f), // 2 = top-right
        Vec3( 1.0f, -1.0f, 0.0f), // 3 = bottom-right
        Vec3(-1.0f, -1.0f, 0.0f), // 4 = bottom-left
        Vec3( 0.5f,  0.5f, 0.0f)  // 5 = new midpoint of edge 0-2
    }, {
        {0, 1, 5},  // center, top-left, midpoint
        {5, 1, 2},  // midpoint, top-left, top-right  (new face from splitting {0,1,2})
        {0, 5, 3},  // center, midpoint, bottom-right
        {5, 2, 3},  // midpoint, top-right, bottom-right (new face from splitting {0,2,3})
        {0, 3, 4},  // unchanged
        {0, 4, 1}   // unchanged
    });

    expect_split(mesh, edge, after);
});
