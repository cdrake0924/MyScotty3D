
#include "halfedge.h"

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <iostream>

/******************************************************************
*********************** Local Operations **************************
******************************************************************/

/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it cannot perform an operation (i.e., because
    the resulting mesh does not have a valid representation).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementation, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/


/*
 * add_face: add a standalone face to the mesh
 *  sides: number of sides
 *  radius: distance from vertices to origin
 *
 * We provide this method as an example of how to make new halfedge mesh geometry.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::add_face(uint32_t sides, float radius) {
	//faces with fewer than three sides are invalid, so abort the operation:
	if (sides < 3) return std::nullopt;


	std::vector< VertexRef > face_vertices;
	//In order to make the first edge point in the +x direction, first vertex should
	// be at -90.0f - 0.5f * 360.0f / float(sides) degrees, so:
	float const start_angle = (-0.25f - 0.5f / float(sides)) * 2.0f * PI_F;
	for (uint32_t s = 0; s < sides; ++s) {
		float angle = float(s) / float(sides) * 2.0f * PI_F + start_angle;
		VertexRef v = emplace_vertex();
		v->position = radius * Vec3(std::cos(angle), std::sin(angle), 0.0f);
		face_vertices.emplace_back(v);
	}

	assert(face_vertices.size() == sides);

	//assemble the rest of the mesh parts:
	FaceRef face = emplace_face(false); //the face to return
	FaceRef boundary = emplace_face(true); //the boundary loop around the face

	std::vector< HalfedgeRef > face_halfedges; //will use later to set ->next pointers

	for (uint32_t s = 0; s < sides; ++s) {
		//will create elements for edge from a->b:
		VertexRef a = face_vertices[s];
		VertexRef b = face_vertices[(s+1)%sides];

		//h is the edge on face:
		HalfedgeRef h = emplace_halfedge();
		//t is the twin, lies on boundary:
		HalfedgeRef t = emplace_halfedge();
		//e is the edge corresponding to h,t:
		EdgeRef e = emplace_edge(false); //false: non-sharp

		//set element data to something reasonable:
		//(most ops will do this with interpolate_data(), but no data to interpolate here)
		h->corner_uv = a->position.xy() / (2.0f * radius) + 0.5f;
		h->corner_normal = Vec3(0.0f, 0.0f, 1.0f);
		t->corner_uv = b->position.xy() / (2.0f * radius) + 0.5f;
		t->corner_normal = Vec3(0.0f, 0.0f,-1.0f);

		//thing -> halfedge pointers:
		e->halfedge = h;
		a->halfedge = h;
		if (s == 0) face->halfedge = h;
		if (s + 1 == sides) boundary->halfedge = t;

		//halfedge -> thing pointers (except 'next' -- will set that later)
		h->twin = t;
		h->vertex = a;
		h->edge = e;
		h->face = face;

		t->twin = h;
		t->vertex = b;
		t->edge = e;
		t->face = boundary;

		face_halfedges.emplace_back(h);
	}

	assert(face_halfedges.size() == sides);

	for (uint32_t s = 0; s < sides; ++s) {
		face_halfedges[s]->next = face_halfedges[(s+1)%sides];
		face_halfedges[(s+1)%sides]->twin->next = face_halfedges[s]->twin;
	}

	return face;
}


/*
 * bisect_edge: split an edge without splitting the adjacent faces
 *  e: edge to split
 *
 * returns: added vertex
 *
 * We provide this as an example for how to implement local operations.
 * (and as a useful subroutine!)
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::bisect_edge(EdgeRef e) {
	// Phase 0: draw a picture
	//
	// before:
	//    ----h--->
	// v1 ----e--- v2
	//   <----t---
	//
	// after:
	//    --h->    --h2->
	// v1 --e-- vm --e2-- v2
	//    <-t2-    <--t--
	//

	// Phase 1: collect existing elements
	HalfedgeRef h = e->halfedge;
	HalfedgeRef t = h->twin;
	VertexRef v1 = h->vertex;
	VertexRef v2 = t->vertex;

	// Phase 2: Allocate new elements, set data
	VertexRef vm = emplace_vertex();
	vm->position = (v1->position + v2->position) / 2.0f;
	interpolate_data({v1, v2}, vm); //set bone_weights

	EdgeRef e2 = emplace_edge();
	e2->sharp = e->sharp; //copy sharpness flag

	HalfedgeRef h2 = emplace_halfedge();
	interpolate_data({h, h->next}, h2); //set corner_uv, corner_normal

	HalfedgeRef t2 = emplace_halfedge();
	interpolate_data({t, t->next}, t2); //set corner_uv, corner_normal

	// The following elements aren't necessary for the bisect_edge, but they are here to demonstrate phase 4
    FaceRef f_not_used = emplace_face();
    HalfedgeRef h_not_used = emplace_halfedge();

	// Phase 3: Reassign connectivity (careful about ordering so you don't overwrite values you may need later!)

	vm->halfedge = h2;

	e2->halfedge = h2;

	assert(e->halfedge == h); //unchanged

	//n.b. h remains on the same face so even if h->face->halfedge == h, no fixup needed (t, similarly)

	h2->twin = t;
	h2->next = h->next;
	h2->vertex = vm;
	h2->edge = e2;
	h2->face = h->face;

	t2->twin = h;
	t2->next = t->next;
	t2->vertex = vm;
	t2->edge = e;
	t2->face = t->face;
	
	h->twin = t2;
	h->next = h2;
	assert(h->vertex == v1); // unchanged
	assert(h->edge == e); // unchanged
	//h->face unchanged

	t->twin = h2;
	t->next = t2;
	assert(t->vertex == v2); // unchanged
	t->edge = e2;
	//t->face unchanged


	// Phase 4: Delete unused elements
    erase_face(f_not_used);
    erase_halfedge(h_not_used);

	// Phase 5: Return the correct iterator
	return vm;
}


/*
 * split_edge: split an edge and adjacent (non-boundary) faces
 *  e: edge to split
 *
 * returns: added vertex. vertex->halfedge should lie along e
 *
 * Note that when splitting the adjacent faces, the new edge
 * should connect to the vertex ccw from the ccw-most end of e
 * within the face.
 *
 * Do not split adjacent boundary faces.
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(EdgeRef e) {
	// A2L2 (REQUIRED): split_edge
    HalfedgeRef h = e->halfedge;
    HalfedgeRef t = h->twin;
    FaceRef f0 = h->face;
    FaceRef f1 = t->face;

    if (f0->boundary && f1->boundary) return std::nullopt;

    HalfedgeRef h_next = h->next;
    HalfedgeRef h_prev = h_next;
    while (h_prev->next != h) h_prev = h_prev->next;

    HalfedgeRef t_next = t->next;
    HalfedgeRef t_prev = t_next;
    while (t_prev->next != t) t_prev = t_prev->next;

    HalfedgeRef h_next_next = h_next->next;
    HalfedgeRef t_next_next = t_next->next;
    VertexRef C0 = h_next_next->vertex;
    VertexRef C1 = t_next_next->vertex;

    auto v_opt = bisect_edge(e);
    if (!v_opt) return std::nullopt;
    VertexRef v = *v_opt;

    HalfedgeRef h_vB = v->halfedge;
    HalfedgeRef h_Av = h;          
    HalfedgeRef t_vA = h_Av->twin;
    HalfedgeRef t_Bv = h_vB->twin;

    HalfedgeRef h_c0_prev = h_next;
    while (h_c0_prev->next != h_next_next) h_c0_prev = h_c0_prev->next;

    EdgeRef   e0   = emplace_edge();
    HalfedgeRef h_vc0 = emplace_halfedge();
    HalfedgeRef h_c0v = emplace_halfedge();
    FaceRef   f2   = emplace_face();

    h_vc0->twin = h_c0v;  
	h_c0v->twin = h_vc0;
    h_vc0->edge = e0;     
	h_c0v->edge = e0;
    e0->halfedge = h_vc0;
    h_vc0->vertex = v;
    h_c0v->vertex = C0;

    h_Av->next = h_vc0;
    h_vc0->next = h_next_next;
    h_vc0->face = f0;
    f0->halfedge = h_Av;

    h_c0_prev->next = h_c0v;
    h_c0v->next = h_vB;
    h_c0v->face = f2;
    h_vB->face = f2;
    {
        HalfedgeRef cur = h_next;
        while (cur != h_c0v) {
            cur->face = f2;
            cur = cur->next;
        }
    }
    f2->halfedge = h_vB;

    if (!f1->boundary) {
        HalfedgeRef t_c1_prev = t_next;
        while (t_c1_prev->next != t_next_next) t_c1_prev = t_c1_prev->next;

        EdgeRef     e1    = emplace_edge();
        HalfedgeRef h_vc1 = emplace_halfedge();
        HalfedgeRef h_c1v = emplace_halfedge();
        FaceRef     f3    = emplace_face();

        h_vc1->twin = h_c1v;  h_c1v->twin = h_vc1;
        h_vc1->edge = e1;     h_c1v->edge = e1;
        e1->halfedge = h_vc1;
        h_vc1->vertex = v;
        h_c1v->vertex = C1;

        t_Bv->next = h_vc1;
        h_vc1->next = t_next_next;
        h_vc1->face = f1;
        t_Bv->face = f1;
        {
            HalfedgeRef cur = t_next_next;
            while (cur != t_Bv) {
                cur->face = f1;
                cur = cur->next;
            }
        }
        f1->halfedge = t_Bv;

        t_c1_prev->next = h_c1v;
        h_c1v->next = t_vA;
        h_c1v->face = f3;
        t_vA->face = f3;
        {
            HalfedgeRef cur = t_next;
            while (cur != h_c1v) {
                cur->face = f3;
                cur = cur->next;
            }
        }
        f3->halfedge = t_vA;
    }

    return v;
}



/*
 * inset_vertex: divide a face into triangles by placing a vertex at f->center()
 *  f: the face to add the vertex to
 *
 * returns:
 *  std::nullopt if insetting a vertex would make mesh invalid
 *  the inset vertex otherwise
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::inset_vertex(FaceRef f) {
	// A2Lx4 (OPTIONAL): inset vertex
	
	(void)f;
    return std::nullopt;
}


/* [BEVEL NOTE] Note on the beveling process:

	Each of the bevel_vertex, bevel_edge, and extrude_face functions do not represent
	a full bevel/extrude operation. Instead, they should update the _connectivity_ of
	the mesh, _not_ the positions of newly created vertices. In fact, you should set
	the positions of new vertices to be exactly the same as wherever they "started from."

	When you click on a mesh element while in bevel mode, one of those three functions
	is called. But, because you may then adjust the distance/offset of the newly
	beveled face, we need another method of updating the positions of the new vertices.

	This is where bevel_positions and extrude_positions come in: these functions are
	called repeatedly as you move your mouse, the position of which determines the
	amount / shrink parameters. These functions are also passed an array of the original
	vertex positions, stored just after the bevel/extrude call, in order starting at
	face->halfedge->vertex, and the original element normal, computed just *before* the
	bevel/extrude call.

	Finally, note that the amount, extrude, and/or shrink parameters are not relative
	values -- you should compute a particular new position from them, not a delta to
	apply.
*/

/*
 * bevel_vertex: creates a face in place of a vertex
 *  v: the vertex to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(VertexRef v) {
	//A2Lx5 (OPTIONAL): Bevel Vertex
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in bevel_vertex_helper (A2Lx5h)

	(void)v;
    return std::nullopt;
}

/*
 * bevel_edge: creates a face in place of an edge
 *  e: the edge to bevel
 *
 * returns: reference to the new face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(EdgeRef e) {
	//A2Lx6 (OPTIONAL): Bevel Edge
	// Reminder: This function does not update the vertex positions.
	// remember to also fill in bevel_edge_helper (A2Lx6h)

	(void)e;
    return std::nullopt;
}

/*
 * extrude_face: creates a face inset into a face
 *  f: the face to inset
 *
 * returns: reference to the inner face
 *
 * see also [BEVEL NOTE] above.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::extrude_face(FaceRef f) {
	//A2L4: Extrude Face
	// Reminder: This function does not update the vertex positions.
	// Remember to also fill in extrude_helper (A2L4h)

	if(f->boundary) return std::nullopt;

	Vec3 base_normal = f->normal();
	if(!std::isfinite(base_normal.x) ||
	   !std::isfinite(base_normal.y) ||
	   !std::isfinite(base_normal.z) ||
	   base_normal.norm() == 0.0f)
	   base_normal = Vec3(0.0f, 0.0f, 0.0f);

	std::vector<HalfedgeRef> old_h;
	std::vector<VertexRef> old_v;

	{
		HalfedgeRef start = f->halfedge;
		if(start == halfedges.end()) return std::nullopt;

		HalfedgeRef curr = start;
		do {
			old_h.emplace_back(curr);
			old_v.emplace_back(curr->vertex);
			curr = curr->next;
			if(old_h.size() > 100000) return std::nullopt;
		} while(curr != start);
	}

	uint32_t n = (uint32_t)old_h.size();
	if(n < 3) std::nullopt;

	std::vector<VertexRef> new_v(n);

	for (uint32_t i = 0; i < n; i++){
		new_v[i] = emplace_vertex();
		new_v[i]->position = old_v[i]->position;
	}

	std::vector<HalfedgeRef> h(n);
	std::vector<EdgeRef> inner_e(n);

	for (uint32_t i = 0; i < n; i++){
		h[i] = emplace_halfedge();
		inner_e[i] = emplace_edge(false);

		h[i]->vertex = new_v[i];
		h[i]->face = f;

		h[i]->edge = inner_e[i];
		inner_e[i]->halfedge = h[i];

		h[i]->corner_normal = base_normal;
	}

	for (uint32_t i = 0; i < n; i++){
		uint32_t i1 = (i + 1) % n;
		h[i]->next = h[i1];
	}

	f->halfedge = h[0];

	for (uint32_t i = 0; i < n; i++) new_v[i]->halfedge = h[i];

	std::vector<FaceRef> quad_f(n);
	std::vector<HalfedgeRef> q1(n), q2(n), q3(n);

	for (uint32_t i = 0; i < n; i++){
		uint32_t i1 = (i + 1) % n;

		quad_f[i] = emplace_face(false);

		HalfedgeRef q0 = old_h[i];

		q1[i] = emplace_halfedge(); 
		q2[i] = emplace_halfedge();
		q3[i] = emplace_halfedge();

		q0->face = quad_f[i];
		q1[i]->face = quad_f[i];
		q2[i]->face = quad_f[i];
		q3[i]->face = quad_f[i];

		quad_f[i]->halfedge = q0;

		q1[i]->vertex = old_v[i1];
		q2[i]->vertex = new_v[i1];
		q3[i]->vertex = new_v[i];

		q0->next = q1[i];
		q1[i]->next = q2[i];
		q2[i]->next = q3[i];
		q3[i]->next = q0;

		q1[i]->corner_normal = base_normal;
		q2[i]->corner_normal = base_normal;
		q3[i]->corner_normal = base_normal;

		q0->corner_normal = base_normal;
	}

	for (uint32_t i = 0; i < n; i++){
		h[i]->twin = q2[i];
		q2[i]->twin = h[i];
		q2[i]->edge = inner_e[i];
	}


	std::vector<EdgeRef> vert_e(n);
	for (uint32_t i = 0; i < n; i++){
		uint32_t ip = (i + n - 1) % n;

		HalfedgeRef up = q1[ip];  
		HalfedgeRef down = q3[i]; 

		vert_e[i] = emplace_edge(false);

		up->twin = down;
		down->twin = up;

		up->edge = vert_e[i];
		down->edge = vert_e[i];
		vert_e[i]->halfedge = up;
	}

	for (uint32_t i = 0; i < n; i++){
		if (old_v[i]->halfedge == halfedges.end() || old_v[i]->halfedge->vertex != old_v[i])
			old_v[i]->halfedge = old_h[i];
	}

	return f;
}

/*
 * flip_edge: rotate non-boundary edge ccw inside its containing faces
 *  e: edge to flip
 *
 * if e is a boundary edge, does nothing and returns std::nullopt
 * if flipping e would create an invalid mesh, does nothing and returns std::nullopt
 *
 * otherwise returns the edge, post-rotation
 *
 * does not create or destroy mesh elements.
 * A2L1
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(EdgeRef e) {
    HalfedgeRef a = e->halfedge;
    HalfedgeRef b = a->twin;

    if (a->face->boundary || b->face->boundary) return std::nullopt;

    HalfedgeRef an = a->next;
    HalfedgeRef ap = an;
    while (ap->next != a) ap = ap->next;

    HalfedgeRef bn = b->next;
    HalfedgeRef bp = bn;
    while (bp->next != b) bp = bp->next;

    HalfedgeRef h = ap->twin;
    do {
        if (h->vertex == bp->vertex) return std::nullopt;
        h = h->next->twin;
    } while (h != ap->twin);

    FaceRef f0 = a->face;
    FaceRef f1 = b->face;

    a->vertex = bn->next->vertex;
    b->vertex = an->next->vertex;

    ap->next = bn;
    bp->next = an;
    a->next  = an->next;
    b->next  = bn->next;
    bn->next = a;
    an->next = b;

    a->face = f0;
	bn->face = f0;
	ap->face = f0;

	b->face = f1;
	an->face = f1;
	bp->face = f1;

    f0->halfedge = a;
    f1->halfedge = b;

    for (HalfedgeRef t : {a, b, an, bn, ap, bp}) {
		t->vertex->halfedge = t;
	}

    return e;
}


/*
 * make_boundary: add non-boundary face to boundary
 *  face: the face to make part of the boundary
 *
 * if face ends up adjacent to other boundary faces, merge them into face
 *
 * if resulting mesh would be invalid, does nothing and returns std::nullopt
 * otherwise returns face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::make_boundary(FaceRef face) {
	//A2Lx7: (OPTIONAL) make_boundary

	return std::nullopt; //TODO: actually write this code!
}

/*
 * dissolve_vertex: merge non-boundary faces adjacent to vertex, removing vertex
 *  v: vertex to merge around
 *
 * if merging would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_vertex(VertexRef v) {
	// A2Lx1 (OPTIONAL): Dissolve Vertex

    return std::nullopt;
}

/*
 * dissolve_edge: merge the two faces on either side of an edge
 *  e: the edge to dissolve
 *
 * merging a boundary and non-boundary face produces a boundary face.
 *
 * if the result of the merge would be an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::dissolve_edge(EdgeRef e) {
	// A2Lx2 (OPTIONAL): dissolve_edge

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data
	
    return std::nullopt;
}

/* collapse_edge: collapse edge to a vertex at its middle
 *  e: the edge to collapse
 *
 * if collapsing the edge would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(EdgeRef e) {
	//A2L3: Collapse Edge

	auto prev_halfedge = [](HalfedgeRef x) -> HalfedgeRef{
		HalfedgeRef prev = x;
		while(prev->next != x) prev = prev->next;
		return prev;
	};

	HalfedgeRef h_ab = e->halfedge;
	HalfedgeRef t_ba = h_ab->twin;
	VertexRef a = h_ab->vertex;
	VertexRef b = t_ba->vertex;

	if(e->on_boundary()){
		HalfedgeRef h_boundary = (h_ab->face->boundary ? h_ab : t_ba);
		HalfedgeRef h_real = h_boundary->twin;

		FaceRef f_boundary = h_boundary->face;
		FaceRef f_real = h_real->face;

		if(f_real->boundary || f_real->degree() != 3) return std::nullopt;

		HalfedgeRef h_bc = h_real->next;
		HalfedgeRef h_ca = h_bc->next;
		VertexRef c = h_ca->vertex;

		a->position = 0.5f * (a->position + b->position);

		{
			HalfedgeRef start = b->halfedge;
			if(start == halfedges.end()) return std::nullopt;
			if(start == h_boundary) start = h_boundary->twin->next;

			HalfedgeRef curr = start;
			for(uint32_t i = 0; i < 100000; i++){
				if(curr == halfedges.end()) break;
				if(curr != h_boundary && curr != h_bc) curr->vertex = a;

				HalfedgeRef n = curr->twin;
				if(n == halfedges.end()) break;
				n = n->next;
				if(n == halfedges.end()) break;

				curr = n;
				if(curr == start) break;
			}
		}
		EdgeRef e_ac = h_ca->edge;
		EdgeRef e_bc = h_bc->edge;
		HalfedgeRef t_ac = h_ca->twin;
		HalfedgeRef t_cb = h_bc->twin;

		t_ac->twin = t_cb;
		t_cb->twin = t_ac;
		t_cb->edge = e_ac;
		e_ac->halfedge = t_ac;

		auto h_boundary_prev = prev_halfedge(h_boundary);
		h_boundary_prev->next = h_boundary->next;
		f_boundary->halfedge = h_boundary_prev->next;

		a->halfedge = t_ac;
		c->halfedge = t_cb;

		erase_face(f_real);
		erase_halfedge(h_real);
		erase_halfedge(h_bc);
		erase_halfedge(h_ca);
		erase_halfedge(h_boundary);
		erase_edge(e);
		erase_edge(e_bc);
		erase_vertex(b);

		return a;
	}

	FaceRef f_h = h_ab->face;
	FaceRef f_t = t_ba->face;
	if(f_h->boundary || f_t->boundary) return std::nullopt;
	if(f_h->degree() <= 3 || f_t->degree() <= 3) return std::nullopt;

	auto h_da = prev_halfedge(h_ab);
	auto t_eb = prev_halfedge(t_ba);

	a->position = 0.5f * (a->position + b->position);

	{
		HalfedgeRef start = b->halfedge;
		if(start == halfedges.end()) return std::nullopt;
		if(start == t_ba) start = t_ba->twin->next;

		HalfedgeRef curr = start;
		for(uint32_t i = 0; i < 100000; i++){
			if(curr == halfedges.end()) break;
			if(curr != t_ba) curr->vertex = a;

			HalfedgeRef n = curr->twin;
			if(n == halfedges.end()) break;
			n = n->next;
			if(n == halfedges.end()) break;

			curr = n;
			if(curr == start) break;
		}
	}

	HalfedgeRef h_bc = h_ab->next;
	HalfedgeRef t_af = t_ba->next;

	h_da->next = h_bc;
	t_eb->next = t_af;

	f_h->halfedge = h_bc;
	f_t->halfedge = t_af;

	a->halfedge = h_bc;
	if(a->halfedge == halfedges.end() || a->halfedge->vertex != a)
		a->halfedge = t_af;
	
	erase_halfedge(h_ab);
	erase_halfedge(t_ba);
	erase_edge(e);
	erase_vertex(b);

	return a;
}

/*
 * collapse_face: collapse a face to a single vertex at its center
 *  f: the face to collapse
 *
 * if collapsing the face would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns the newly collapsed vertex
 */
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(FaceRef f) {
	//A2Lx3 (OPTIONAL): Collapse Face

	//Reminder: use interpolate_data() to merge corner_uv / corner_normal data on halfedges
	// (also works for bone_weights data on vertices!)

    return std::nullopt;
}

/*
 * weld_edges: glue two boundary edges together to make one non-boundary edge
 *  e, e2: the edges to weld
 *
 * if welding the edges would result in an invalid mesh, does nothing and returns std::nullopt
 * otherwise returns e, updated to represent the newly-welded edge
 */
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::weld_edges(EdgeRef e, EdgeRef e2) {
	//A2Lx8: Weld Edges

	//Reminder: use interpolate_data() to merge bone_weights data on vertices!

    return std::nullopt;
}



/*
 * bevel_positions: compute new positions for the vertices of a beveled vertex/edge
 *  face: the face that was created by the bevel operation
 *  start_positions: the starting positions of the vertices
 *     start_positions[i] is the starting position of face->halfedge(->next)^i
 *  direction: direction to bevel in (unit vector)
 *  distance: how far to bevel
 *
 * push each vertex from its starting position along its outgoing edge until it has
 *  moved distance `distance` in direction `direction`. If it runs out of edge to
 *  move along, you may choose to extrapolate, clamp the distance, or do something
 *  else reasonable.
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after bevel_vertex or bevel_edge.
 * (So you can assume the local topology is set up however your bevel_* functions do it.)
 *
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::bevel_positions(FaceRef face, std::vector<Vec3> const &start_positions, Vec3 direction, float distance) {
	//A2Lx5h / A2Lx6h (OPTIONAL): Bevel Positions Helper
	
	// The basic strategy here is to loop over the list of outgoing halfedges,
	// and use the preceding and next vertex position from the original mesh
	// (in the start_positions array) to compute an new vertex position.
	
}

/*
 * extrude_positions: compute new positions for the vertices of an extruded face
 *  face: the face that was created by the extrude operation
 *  move: how much to translate the face
 *  shrink: amount to linearly interpolate vertices in the face toward the face's centroid
 *    shrink of zero leaves the face where it is
 *    positive shrink makes the face smaller (at shrink of 1, face is a point)
 *    negative shrink makes the face larger
 *
 * only changes vertex positions (no connectivity changes!)
 *
 * This is called repeatedly as the user interacts, just after extrude_face.
 * (So you can assume the local topology is set up however your extrude_face function does it.)
 *
 * Using extrude face in the GUI will assume a shrink of 0 to only extrude the selected face
 * Using bevel face in the GUI will allow you to shrink and increase the size of the selected face
 * 
 * see also [BEVEL NOTE] above.
 */
void Halfedge_Mesh::extrude_positions(FaceRef face, Vec3 move, float shrink) {
	//A2L4h: Extrude Positions Helper

	//General strategy:
	// use mesh navigation to get starting positions from the surrounding faces,
	// compute the centroid from these positions + use to shrink,
	// offset by move
	
	if(face->boundary) return;

	HalfedgeRef start = face->halfedge;
	if(start == halfedges.end()) return;

	std::vector<HalfedgeRef> h;
	std::vector<VertexRef> v;
	std::vector<Vec3> pos;

	{
		HalfedgeRef curr = start;
		do {
			h.emplace_back(curr);
			v.emplace_back(curr->vertex);

			HalfedgeRef htn = curr->twin->next;
			VertexRef old_v = htn->twin->vertex;
			pos.emplace_back(old_v->position);

			curr = curr->next;
			if(h.size() > 100000) return;
		} while(curr != start);
	}

	size_t n = h.size();
	if(n < 3) return;

	Vec3 centroid(0.0f);
	for(size_t i = 0; i < n; i++) centroid += pos[i];
	centroid /= float(n);

	float scale = 1.0f - shrink;
	for(size_t i = 0; i < n; i++){
		Vec3 p = centroid + (pos[i] - centroid) * scale + move;
		v[i]->position = p;
	}
}

