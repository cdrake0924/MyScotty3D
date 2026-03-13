#include "halfedge.h"

#include <unordered_map>
#include <unordered_set>

#include <iostream>


/*
 * triangulate: split all non-boundary faces into triangles.
 *
 * Works on all valid meshes.
 */
void Halfedge_Mesh::triangulate() {
	//A2G1: triangulation
	// Collect all non-boundary faces upfront (so iteration is safe while we modify)
	std::vector<FaceRef> to_triangulate;
	for (FaceRef f = faces.begin(); f != faces.end(); ++f) {
		if (!f->boundary) {
			to_triangulate.push_back(f);
		}
	}

	for (FaceRef f : to_triangulate) {
		std::cout << "[triangulate] Processing face id=" << f->id << std::endl;
		// Fan triangulation from the first vertex of the face.
		// For a face with halfedges h0 -> h1 -> h2 -> ... -> h_{n-1} -> h0,
		// we repeatedly "cut off" a triangle at the front: [h0->vertex, h1->vertex, h2->vertex]
		// leaving a smaller polygon.

		// Keep slicing until the face is a triangle.
		while (true) {
			// Count edges in face
			HalfedgeRef h0 = f->halfedge;
			HalfedgeRef h1 = h0->next;
			HalfedgeRef h2 = h1->next;
			if (h2->next == h0) break; // already a triangle

			// Create new diagonal edge: from h0->vertex to h2->vertex
			// This diagonal splits f into:
			//   triangle: h0, h1, new_he_back  (new face tri)
			//   polygon:  new_he_fwd, h2, ..., h_{n-1}  (stays in f)

			EdgeRef new_edge = emplace_edge();
			HalfedgeRef new_he_fwd = emplace_halfedge(); // lives in the new triangle, pointing h2->vertex -> h0->vertex
			HalfedgeRef new_he_back = emplace_halfedge(); // lives in remaining polygon, pointing h0->vertex -> h2->vertex
			FaceRef new_tri = emplace_face();

			// The new triangle face uses: h0, h1, new_he_fwd
			// new_he_fwd points from h2->vertex to h0->vertex (closes the triangle)
			// new_he_back points from h0->vertex to h2->vertex (opening of remaining polygon)

			// Assign new edge halfedges
			new_edge->halfedge = new_he_fwd;

			// Set up new_he_fwd (inside the new triangle)
			new_he_fwd->set_tnvef(new_he_back, h0, h2->vertex, new_edge, new_tri);

			// Set up new_he_back (inside the remaining polygon f)
			new_he_back->set_tnvef(new_he_fwd, h2, h0->vertex, new_edge, f);

			// Update h0 and h1 to be in new_tri
			h0->face = new_tri;
			h1->face = new_tri;
			h1->next = new_he_fwd; // close the triangle

			// new_tri starts at h0
			new_tri->halfedge = h0;

			// Remaining polygon: new_he_back -> h2 -> ... -> h_{n-1} -> h0 (now h0 is NOT in f)
			// Find h_{n-1}: last halfedge that points back to h0
			HalfedgeRef last = h2;
			while (last->next != h0) last = last->next;
			last->next = new_he_back;

			// f now starts at new_he_back
			f->halfedge = new_he_back;

			// std::cout << "  new_tri halfedges: v" << h0->vertex->id
			//           << " -> v" << h1->vertex->id
			//           << " -> v" << new_he_fwd->vertex->id
			//           << " -> (back to v" << h0->vertex->id << ")" << std::endl;
			// std::cout << "  remaining f halfedges: v" << new_he_back->vertex->id
			//           << " -> v" << h2->vertex->id << " -> ..." << std::endl;

			// Interpolate corner data for new halfedges
			interpolate_data({HalfedgeCRef(h1), HalfedgeCRef(h2)}, new_he_fwd);
			interpolate_data({HalfedgeCRef(h0)}, new_he_back); 
		}
	}
}

/*
 * linear_subdivide: split faces into quads without moving anything.
 *
 * Works on all valid meshes.
 *
 * (NOTE: uses catmark_subdivide_helper for subdivision)
 */
void Halfedge_Mesh::linear_subdivide() {
	std::unordered_map< VertexCRef, Vec3 > vertex_positions;
	std::unordered_map< EdgeCRef, Vec3 > edge_vertex_positions;
	std::unordered_map< FaceCRef, Vec3 > face_vertex_positions;

	//A2G2: linear subdivision

	// For every vertex, assign its current position to vertex_positions[v]:
	for (VertexCRef v = vertices.begin(); v != vertices.end(); ++v) {
		vertex_positions[v] = v->position;
	}

    // For every edge, assign the midpoint of its adjacent vertices to edge_vertex_positions[e]:
	for (EdgeCRef e = edges.begin(); e != edges.end(); ++e) {
		edge_vertex_positions[e] = e->center();
	}

    // For every *non-boundary* face, assign the centroid (i.e., arithmetic mean) to face_vertex_positions[f]:
	for (FaceCRef f = faces.begin(); f != faces.end(); ++f) {
		if (!f->boundary) {
			face_vertex_positions[f] = f->center();
		}
	}

	//use the helper function to actually perform the subdivision:
	catmark_subdivide_helper(vertex_positions, edge_vertex_positions, face_vertex_positions);
}

/*
 * catmark_subdivide: split faces into quads with positions calculated by
 *   the Catmull-Clark ruleset.
 *
 * Works on all valid meshes.
 *
 * (NOTE: uses catmark_subdivide_helper for subdivision)
 */
void Halfedge_Mesh::catmark_subdivide() {
	std::unordered_map< VertexCRef, Vec3 > vertex_positions;
	std::unordered_map< EdgeCRef, Vec3 > edge_vertex_positions;
	std::unordered_map< FaceCRef, Vec3 > face_vertex_positions;

	//A2G3: Catmull-Clark Subdivision

	// --- Face points ---
	// Each non-boundary face gets a new vertex at its centroid.
	for (FaceCRef f = faces.begin(); f != faces.end(); ++f) {
		if (!f->boundary) {
			face_vertex_positions[f] = f->center();
		}
	}

	// --- Edge points ---
	// For interior edges: average of the two endpoints and the two adjacent face points.
	// For boundary edges: simple midpoint of the two endpoints.
	for (EdgeCRef e = edges.begin(); e != edges.end(); ++e) {
		HalfedgeCRef h  = e->halfedge;
		HalfedgeCRef ht = h->twin;
		Vec3 v0 = h->vertex->position;
		Vec3 v1 = ht->vertex->position;

		if (e->on_boundary()) {
			// Boundary rule: midpoint of endpoints
			edge_vertex_positions[e] = (v0 + v1) / 2.0f;
		} else {
			// Interior rule: (v0 + v1 + face_point_0 + face_point_1) / 4
			Vec3 f0 = face_vertex_positions.at(h->face);
			Vec3 f1 = face_vertex_positions.at(ht->face);
			edge_vertex_positions[e] = (v0 + v1 + f0 + f1) / 4.0f;
		}
	}

	// --- Vertex points ---
	// For boundary vertices: (1/8)*prev_boundary + (6/8)*v + (1/8)*next_boundary
	// For interior vertices: (F + 2E + (n-3)*P) / n
	//   where F = avg of adjacent face points, E = avg of adjacent edge midpoints,
	//   P = original vertex position, n = valence (number of adjacent faces)
	for (VertexCRef v = vertices.begin(); v != vertices.end(); ++v) {
		if (v->on_boundary()) {
			// Find the two boundary edges incident to v and average their other endpoints
			Vec3 boundary_sum = Vec3{0.0f, 0.0f, 0.0f};
			int boundary_count = 0;
			HalfedgeCRef h = v->halfedge;
			do {
				if (h->edge->on_boundary()) {
					// The other endpoint of this edge
					boundary_sum += h->twin->vertex->position;
					boundary_count++;
				}
				h = h->twin->next;
			} while (h != v->halfedge);
			// Crease rule: 1/8 * sum_of_2_boundary_neighbors + 6/8 * P
			vertex_positions[v] = (boundary_sum / 8.0f) + (6.0f / 8.0f) * v->position;
		} else {
			// Interior vertex: Catmull-Clark weighted average
			// n = valence
			int n = 0;
			Vec3 F{0.0f, 0.0f, 0.0f}; // sum of adjacent face points
			Vec3 E{0.0f, 0.0f, 0.0f}; // sum of adjacent edge midpoints
			HalfedgeCRef h = v->halfedge;
			do {
				n++;
				F += face_vertex_positions.at(h->face);
				E += (v->position + h->twin->vertex->position) / 2.0f;
				h = h->twin->next;
			} while (h != v->halfedge);

			float fn = static_cast<float>(n);
			vertex_positions[v] = (F / fn + 2.0f * E / fn + (fn - 3.0f) * v->position) / fn;
		}
	}

	
	//Now, use the provided helper function to actually perform the subdivision:
	catmark_subdivide_helper(vertex_positions, edge_vertex_positions, face_vertex_positions);

}

/*
 * loop_subdivide: sub-divide non-boundary faces with the Loop subdivision rule
 * 
 * If all non-boundary faces are triangles:
 *   subdivides mesh using the Loop subdivision rule
 *   returns true
 * Otherwise:
 *   does not change mesh
 *   returns false
 *
 * Do note that this requires a working implementation of edge split and edge flip
 */
bool Halfedge_Mesh::loop_subdivide() {

	//preamble: check for any non-triangular non-boundary faces:
	for (FaceCRef f = faces.begin(); f != faces.end(); ++f) {
		if (f->boundary) continue; //ignore boundary faces for this check
		if (f->halfedge->next->next->next != f->halfedge) {
			//found a non-triangular face!
			return false;
		}
	}

	//if execution reaches this point, all non-boundary faces are triangular, so proceed to subdivide:

	// A2Go1: Loop subdivision.
	// ---------------------------------------------------------------
	// Step 1: Compute new positions for ALL OLD vertices using the
	//         Loop subdivision weighting rule, store in vertex_new_pos.
	// ---------------------------------------------------------------
	std::unordered_map<VertexRef, Vec3> vertex_new_pos;
 
	for (VertexRef v = vertices.begin(); v != vertices.end(); ++v) {
		uint32_t n = v->degree();
 
		if (v->on_boundary()) {
			// Boundary rule: 1/8*(left + right) + 3/4*v
			Vec3 boundary_sum{0.0f, 0.0f, 0.0f};
			HalfedgeRef h = v->halfedge;
			do {
				if (h->edge->on_boundary()) {
					boundary_sum += h->twin->vertex->position;
				}
				h = h->twin->next;
			} while (h != v->halfedge);
			vertex_new_pos[v] = (1.0f / 8.0f) * boundary_sum + (3.0f / 4.0f) * v->position;
		} else {
			// Interior rule: (1 - n*beta)*v + beta*sum_of_neighbors
			// Warren's formula: beta = 3/16 if n==3, else 3/(8n)
			float beta = (n == 3) ? (3.0f / 16.0f) : (3.0f / (8.0f * float(n)));
 
			Vec3 neighbor_sum{0.0f, 0.0f, 0.0f};
			HalfedgeRef h = v->halfedge;
			do {
				neighbor_sum += h->twin->vertex->position;
				h = h->twin->next;
			} while (h != v->halfedge);
 
			vertex_new_pos[v] = (1.0f - float(n) * beta) * v->position + beta * neighbor_sum;
		}
	}
 
	// ---------------------------------------------------------------
	// Step 2: Compute positions for new edge-midpoint vertices using
	//         the Loop edge rule, store in edge_new_pos.
	// ---------------------------------------------------------------
	std::unordered_map<EdgeRef, Vec3> edge_new_pos;
 
	for (EdgeRef e = edges.begin(); e != edges.end(); ++e) {
		HalfedgeRef h = e->halfedge;
		HalfedgeRef t = h->twin;
 
		if (e->on_boundary()) {
			// Boundary rule: simple midpoint
			edge_new_pos[e] = (h->vertex->position + t->vertex->position) / 2.0f;
		} else {
			// Interior rule: 3/8*(A+B) + 1/8*(C+D)
			// A,B = endpoints; C,D = opposite vertices of the two adjacent triangles
			Vec3 A = h->vertex->position;
			Vec3 B = t->vertex->position;
			Vec3 C = h->next->next->vertex->position;
			Vec3 D = t->next->next->vertex->position;
			edge_new_pos[e] = (3.0f / 8.0f) * (A + B) + (1.0f / 8.0f) * (C + D);
		}
	}
 
	// ---------------------------------------------------------------
	// Step 3: Split every OLD edge. Snapshot the last old edge first
	//         so newly-created edges are not processed.
	//         Place each new midpoint vertex at its precomputed position.
	// ---------------------------------------------------------------
	std::unordered_set<VertexRef> is_new_vertex;
	std::vector<EdgeRef> new_cross_edges;
 
	// Record all original edge ids BEFORE any splits. After bisect_edge,
	// the original EdgeRef is kept as one half (A->v) and a brand-new EdgeRef
	// is created for the other half (v->B). That new half-edge is NOT in
	// old_edge_ids, so we can use this set to distinguish old-edge halves
	// from true cross-edges.
	std::unordered_set<uint32_t> old_edge_ids;
	for (EdgeRef e = edges.begin(); e != edges.end(); ++e) {
		old_edge_ids.insert(e->id);
	}
 
	EdgeRef last_old_edge = std::prev(edges.end());
 
	for (EdgeRef e = edges.begin(); ; ++e) {
		Vec3 new_pos = edge_new_pos.at(e);
 
		// bisect_edge (called inside split_edge) keeps the original EdgeRef e
		// as the A->v half, and creates a brand-new EdgeRef (call it e2) for
		// the v->B half. vm->halfedge is set to h2, which lies on e2.
		// So after split_edge returns, v_new->halfedge->edge == e2 (new, not
		// in old_edge_ids), while e itself is still in old_edge_ids.
		// We must skip BOTH e and e2 — only e0/e1 (the diagonal cross-cuts) are
		// the edges we want to flip.
		auto result = split_edge(e);
		if (result.has_value()) {
			VertexRef v_new = result.value();
			v_new->position = new_pos;
			is_new_vertex.insert(v_new);
 
			// e2 is the other old-edge half: it lives on v_new->halfedge
			// (bisect_edge sets vm->halfedge = h2 which is on e2).
			EdgeRef e2 = v_new->halfedge->edge;
 
			HalfedgeRef h = v_new->halfedge;
			do {
				EdgeRef ei = h->edge;
				// Skip the two old-edge halves: e (still has its original id,
				// so in old_edge_ids) and e2 (new id but is the bisected half,
				// not a cross-edge — identified explicitly above).
				if (!old_edge_ids.count(ei->id) && ei != e2) {
					new_cross_edges.push_back(ei);
				}
				h = h->twin->next;
			} while (h != v_new->halfedge);
		}
 
		if (e == last_old_edge) break;
	}
 
	// ---------------------------------------------------------------
	// Step 4: Flip every new cross-edge that connects one old vertex
	//         and one new vertex. This gives the correct Loop topology.
	// ---------------------------------------------------------------
	for (EdgeRef e : new_cross_edges) {
		VertexRef va = e->halfedge->vertex;
		VertexRef vb = e->halfedge->twin->vertex;
		bool a_new = is_new_vertex.count(va) > 0;
		bool b_new = is_new_vertex.count(vb) > 0;
		if (a_new != b_new) {
			flip_edge(e);
		}
	}
 
	// ---------------------------------------------------------------
	// Step 5: Write the precomputed positions back to the old vertices.
	//         (New vertices were already positioned in Step 3.)
	// ---------------------------------------------------------------
	for (auto& [v, pos] : vertex_new_pos) {
		v->position = pos;
	}
 
	return true;
}

//isotropic_remesh: improves mesh quality through local operations.
// Do note that this requires a working implementation of EdgeSplit, EdgeFlip, and EdgeCollapse
void Halfedge_Mesh::isotropic_remesh(Isotropic_Remesh_Parameters const &params) {

	//A2Go2: Isotropic Remeshing
	// Optional! Only one of {A2Go1, A2Go2, A2Go3} is required!

	// Compute the mean edge length. This will be the "target length".

    // Repeat the four main steps for `outer_iterations` iterations:

    // -> Split edges much longer than the target length.
	//     ("much longer" means > target length * params.longer_factor)

    // -> Collapse edges much shorter than the target length.
	//     ("much shorter" means < target length * params.shorter_factor)

    // -> Flip each edge if it improves vertex degree.

    // -> Finally, apply some tangential smoothing to the vertex positions.
	//     This means move every vertex in the plane of its normal,
	//     toward the centroid of its neighbors, by params.smoothing_step of
	//     the total distance (so, smoothing_step of 1 would move all the way,
	//     smoothing_step of 0 would not move).
	// -> Repeat the tangential smoothing part params.smoothing_iterations times.

	//NOTE: many of the steps in this function will be modifying the element
	//      lists they are looping over. Take care to avoid use-after-free
	//      or infinite-loop problems.

}

struct Edge_Record {
	Edge_Record() {
	}
	Edge_Record(std::unordered_map<uint32_t, Mat4>& VQ, Halfedge_Mesh::EdgeRef e) : edge(e) {
		// Combined quadric for this edge = sum of the two endpoint quadrics
		Halfedge_Mesh::VertexRef v0 = e->halfedge->vertex;
		Halfedge_Mesh::VertexRef v1 = e->halfedge->twin->vertex;
		Mat4 Q = VQ.at(v0->id) + VQ.at(v1->id);
 
		// Try to solve the 3x3 linear system:
		//   [Q[0][0] Q[1][0] Q[2][0]] [x]   [-Q[3][0]]
		//   [Q[0][1] Q[1][1] Q[2][1]] [y] = [-Q[3][1]]
		//   [Q[0][2] Q[1][2] Q[2][2]] [z]   [-Q[3][2]]
		// (Mat4 is column-major: Q[col][row])
		Mat4 A = Q;
		// Zero out the 4th row and column to isolate the 3x3 system,
		// then set the bottom-right to 1 so it's invertible independently.
		A[0][3] = 0.0f; A[1][3] = 0.0f; A[2][3] = 0.0f; A[3][3] = 1.0f;
		A[3][0] = 0.0f; A[3][1] = 0.0f; A[3][2] = 0.0f;
 
		float det = A.det();
		if (std::abs(det) > 1e-8f) {
			// Solvable: optimal position minimizes quadric error
			Mat4 Ainv = A.inverse();
			Vec4 rhs = Vec4{-Q[3][0], -Q[3][1], -Q[3][2], 1.0f};
			Vec4 sol  = Ainv * rhs;
			optimal = Vec3{sol.x, sol.y, sol.z};
		} else {
			// Degenerate: fall back to whichever endpoint (or midpoint) has lower cost
			Vec3 mid = (v0->position + v1->position) / 2.0f;
			auto cost3 = [&](Vec3 p) -> float {
				Vec4 hp{p.x, p.y, p.z, 1.0f};
				return dot(hp, Q * hp);
			};
			float c0  = cost3(v0->position);
			float c1  = cost3(v1->position);
			float cm  = cost3(mid);
			if (c0 <= c1 && c0 <= cm)      optimal = v0->position;
			else if (c1 <= c0 && c1 <= cm) optimal = v1->position;
			else                            optimal = mid;
		}
 
		// Score = quadric error at the optimal position
		Vec4 ho{optimal.x, optimal.y, optimal.z, 1.0f};
		score = dot(ho, Q * ho);
	}
	Halfedge_Mesh::EdgeRef edge;
	Vec3 optimal;
	float score;
};
 
bool operator<(const Edge_Record& r1, const Edge_Record& r2) {
	if (r1.score != r2.score) {
		return (r1.score < r2.score);
	}
	Halfedge_Mesh::EdgeRef e1 = r1.edge;
	Halfedge_Mesh::EdgeRef e2 = r2.edge;
	return &*e1 < &*e2;
}
 
template<class T> struct MutablePriorityQueue {
	void insert(const T& item) {
		queue.insert(item);
	}
	void remove(const T& item) {
		if (queue.find(item) != queue.end()) {
			queue.erase(item);
		}
	}
	const T& top() const {
		return *(queue.begin());
	}
	void pop() {
		queue.erase(queue.begin());
	}
	size_t size() {
		return queue.size();
	}
 
	std::set<T> queue;
};
 
/*
 * simplify: reduce edge count through collapses
 *  ratio: proportion of original faces to retain
 *
 * you may choose to have your implementation work only on triangle meshes,
 *  in which case it may return 'false' if there are non-triangular
 *  non-boundary faces
 *
 * returns false if it ran out of edges to collapse
 * returns true otherwise
 * 
 * Do note that this requires a working implementation of EdgeCollapse
 */
bool Halfedge_Mesh::simplify(float ratio) {
 
	//A2Go3: simplification
 
	std::unordered_map<uint32_t, Mat4> face_quadrics;
	std::unordered_map<uint32_t, Mat4> vertex_quadrics;
	std::unordered_map<uint32_t, Edge_Record> edge_records;
	MutablePriorityQueue<Edge_Record> queue;
 
	std::cerr << "[simplify] START ratio=" << ratio
	          << " faces=" << faces.size()
	          << " edges=" << edges.size()
	          << " verts=" << vertices.size() << std::endl;
 
	// ---------------------------------------------------------------
	// Step 1: Compute a quadric for each non-boundary face.
	// ---------------------------------------------------------------
	for (FaceRef f = faces.begin(); f != faces.end(); ++f) {
		if (f->boundary) continue;
 
		Vec3 n = f->normal().unit();
		float d = -dot(n, f->halfedge->vertex->position);
		Vec4 p{n.x, n.y, n.z, d};
 
		Mat4 K;
		for (int col = 0; col < 4; ++col) {
			for (int row = 0; row < 4; ++row) {
				K[col][row] = p[row] * p[col];
			}
		}
		face_quadrics[f->id] = K;
	}
	std::cerr << "[simplify] Step 1 done: " << face_quadrics.size() << " face quadrics" << std::endl;
 
	// ---------------------------------------------------------------
	// Step 2: Compute each vertex's quadric = sum of incident face quadrics.
	// ---------------------------------------------------------------
	for (VertexRef v = vertices.begin(); v != vertices.end(); ++v) {
		Mat4 Q = Mat4::Zero;
		HalfedgeRef h = v->halfedge;
		do {
			if (!h->face->boundary) {
				Q = Q + face_quadrics.at(h->face->id);
			}
			h = h->twin->next;
		} while (h != v->halfedge);
		vertex_quadrics[v->id] = Q;
	}
	std::cerr << "[simplify] Step 2 done: " << vertex_quadrics.size() << " vertex quadrics" << std::endl;
 
	// ---------------------------------------------------------------
	// Step 3: Build the priority queue — one Edge_Record per edge.
	// ---------------------------------------------------------------
	for (EdgeRef e = edges.begin(); e != edges.end(); ++e) {
		Edge_Record rec(vertex_quadrics, e);
		edge_records[e->id] = rec;
		queue.insert(rec);
	}
	std::cerr << "[simplify] Step 3 done: " << queue.size() << " edges in queue" << std::endl;
 
	// ---------------------------------------------------------------
	// Step 4: Collapse loop.
	// ---------------------------------------------------------------
	size_t current_faces = 0;
	for (FaceRef f = faces.begin(); f != faces.end(); ++f) {
		if (!f->boundary) current_faces++;
	}
	size_t target_faces = std::max(size_t(1),
	    static_cast<size_t>(std::round(ratio * float(current_faces))));
 
	std::cerr << "[simplify] Step 4 start: current_faces=" << current_faces
	          << " target_faces=" << target_faces << std::endl;
 
	size_t iter = 0;
	size_t collapses = 0;
	size_t rejections = 0;
 
	while (current_faces > target_faces) {
		if (queue.size() == 0) return false;

		Edge_Record best = queue.top();
		queue.pop();

		EdgeRef e = best.edge;

		VertexRef v0 = e->halfedge->vertex;
		VertexRef v1 = e->halfedge->twin->vertex;

		uint32_t id0 = v0->id;
		uint32_t id1 = v1->id;

		// skip stale records
		auto it = edge_records.find(e->id);
		if (it == edge_records.end()) continue;
		if (it->second.score != best.score) continue;

		auto result = collapse_edge(e);
		if (!result.has_value()) {
			++rejections;
			continue;
		}

		++collapses;

		VertexRef v_new = result.value();
		v_new->position = best.optimal;

		vertex_quadrics[v_new->id] =
			vertex_quadrics[id0] + vertex_quadrics[id1];

		current_faces -= 2;

		// rebuild incident edges
		HalfedgeRef h = v_new->halfedge;
		do {
			EdgeRef incident = h->edge;
			// Remove the OLD record from the queue before inserting a new one
			auto old_it = edge_records.find(incident->id);
			if (old_it != edge_records.end()) {
				queue.remove(old_it->second);
			}
			Edge_Record rec(vertex_quadrics, incident);
			edge_records[incident->id] = rec;
			queue.insert(rec);
			h = h->twin->next;
		} while (h != v_new->halfedge);
	}

	std::cerr << "[simplify] DONE after " << iter << " iters ("
	          << collapses << " collapses, " << rejections << " rejections). "
	          << "final_faces=" << current_faces << std::endl;
 
	return true;
}

/*
 * catmark_subdivide_helper: add vertex in every edge and non-boundary face, set positions from parameters
 *
 * Works on all valid meshes.
 */
void Halfedge_Mesh::catmark_subdivide_helper(
	std::unordered_map< VertexCRef, Vec3 > const &vertex_positions, //positions for vertices after subdivision
	std::unordered_map< EdgeCRef, Vec3 > const &edge_vertex_positions, //positions for new vertices added in each edge
	std::unordered_map< FaceCRef, Vec3 > const &face_vertex_positions //positions for new vertices added in each face
	) {

	//check that positions were supplied for every vertex:
	for (VertexCRef v = vertices.begin(); v != vertices.end(); ++v) {
		if (!vertex_positions.count(v)) {
			throw std::runtime_error("No vertex position supplied for vertex with id " + std::to_string(v->id) + ".");
		}
	}

	//check that positions were supplied for every edge:
	for (EdgeCRef e = edges.begin(); e != edges.end(); ++e) {
		if (!edge_vertex_positions.count(e)) {
			throw std::runtime_error("No edge vertex position supplied for edge with id " + std::to_string(e->id) + ".");
		}
	}

	//check that positions were supplied for every (non-boundary) face:
	for (FaceCRef f = faces.begin(); f != faces.end(); ++f) {
		if (f->boundary) {
			if (face_vertex_positions.count(f)) {
				throw std::runtime_error("Extraneous vertex position was supplied for boundary face with id " + std::to_string(f->id) + ".");
			}
		} else {
			if (!face_vertex_positions.count(f)) {
				throw std::runtime_error("No vertex position supplied for face with id " + std::to_string(f->id) + ".");
			}
		}
	}


	{ //check that mesh is in a valid state to start with:
		auto error = validate();
		if (error) {
			throw std::runtime_error("catmark_subdivide_helper called on invalid mesh: " + error.value().second);
		}
	}

	if (vertices.empty() || edges.empty() || faces.empty()) {
		//empty mesh must be empty:
		assert(vertices.empty() && edges.empty() && faces.empty());
		return;
	}

	//store the old last vertex, face, and edge to allow iterating over only the old elements later:
	//(this works because the emplace_* functions add to the end of the element lists)
	VertexRef last_old_vertex = std::prev(vertices.end());
	EdgeRef last_old_edge = std::prev(edges.end());
	FaceRef last_old_face = std::prev(faces.end());
	
	//(can't store .end() iterators because emplace_back puts things "before the end")

	//------------------------
	//split every edge:
	//old halfedges stay connected to their vertices
	//old edge stays connected to e->halfedge->vertex

	//before:
	//     -----h---->
	//  v1 -----e----- v2
	//     <----t-----
	//after:
	//     --h->    --h2->
	//  v1 --e-- vm --e2-- v2
	//     <-t2-    <--t--

	for (EdgeRef e = edges.begin(); e != std::next(last_old_edge); ++e) {
		HalfedgeRef h = e->halfedge;
		HalfedgeRef t = h->twin; assert(t->edge == e);
		VertexRef v1 = h->vertex;
		VertexRef v2 = t->vertex;

		//new elements:
		VertexRef vm = emplace_vertex();
		HalfedgeRef h2 = emplace_halfedge();
		HalfedgeRef t2 = emplace_halfedge();
		EdgeRef e2 = emplace_edge(e->sharp);

		//middle vertex:
		vm->halfedge = h2; //could also use t2
		vm->position = edge_vertex_positions.at(e);
		interpolate_data({v1, v2}, vm);

		//second edge:
		e2->halfedge = h2;

		//second halfedge:
		h2->next = h->next;
		h2->twin = t;
		h2->vertex = vm;
		h2->edge = e2;
		h2->face = h->face;
		interpolate_data({h, h->next}, h2);

		//second twin halfedge:
		t2->next = t->next;
		t2->twin = h;
		t2->vertex = vm;
		t2->edge = e;
		t2->face = t->face;
		interpolate_data({t, t->next}, t2);

		//fix up pointers for existing halfedges:
		h->next = h2;
		h->twin = t2;

		t->next = t2;
		t->twin = h2;
		t->edge = e2;
	}


	//---------------------------
	//split (non-boundary) faces:

	//before:
	//
	//  v0 <-h7- v7 <-h6- v6
	//  |                 ^
	//  h0                h5
	//  v                 |
	//  v1       f        v5
	//  |                 ^
	//  h1                h4
	//  v                 |
	//  v2 -h2-> v3 -h3-> v4
	//
	//after:
	//  v0 <-h7- v7 <-h6- v6
	//  |        |        ^
	//  h0   f   e3   f3  h5
	//  v  --c-> |        |
	//  v1 --e0- vm --e2- v5 
	//  |  <-t-- |        ^
	//  h1   f1  e1  f2   h4
	//  v        |        |
	//  v2 -h2-> v3 -h3-> v4
	//
	// (each new eN has new halfedges as you'd expect,
	//  with eN->halfedge being directed toward the central vertex.)

	for (FaceRef f = faces.begin(); f != std::next(last_old_face); ++f) {
		if (f->boundary) continue; //skip boundary faces

		//get face halfedges:
		std::vector< HalfedgeRef > face_halfedges;
		{
			HalfedgeRef h = f->halfedge;
			do {
				face_halfedges.emplace_back(h);
				h = h->next;
			} while (h != f->halfedge);
			assert(face_halfedges.size() % 2 == 0); //should always be pairs of halfedges along subdivided edges!
		}

		//get face vertices and corners to interpolate data from:
		// (skip the odd vertices/halfedges -- they were just added)
		std::vector< HalfedgeCRef > from_corners;
		std::vector< VertexCRef > from_vertices;
		for (uint32_t i = 0; i < face_halfedges.size(); i += 2) {
			from_corners.emplace_back(face_halfedges[i]);
			from_vertices.emplace_back(face_halfedges[i]->vertex);
		}

		//add central vertex:
		VertexRef vm = emplace_vertex();
		vm->position = face_vertex_positions.at(f);
		interpolate_data(from_vertices, vm);

		//add halfedges and edges around the central vertex:
		std::vector< EdgeRef > inner_edges;
		for (uint32_t i = 0; i + 1 < face_halfedges.size(); i += 2) {
			EdgeRef e = emplace_edge(false);
			HalfedgeRef c = emplace_halfedge();
			HalfedgeRef t = emplace_halfedge();

			e->halfedge = c;

			//halfedge coming from the side:
			c->twin = t;
			//c->next will be set later
			c->vertex = face_halfedges[i+1]->vertex;
			c->edge = e;
			//c->face will be set later
			interpolate_data({face_halfedges[i+1]}, c); //just copy the data

			//halfedge coming from the center:
			t->twin = c;
			//t->next will be set later
			t->vertex = vm;
			t->edge = e;
			//t->face will be set later
			interpolate_data(from_corners, t);

			if (i == 0) vm->halfedge = t;

			//save edge for later connection:
			inner_edges.emplace_back(e);
		}

		//hook up pointers for all the quads:
		for (uint32_t i = 0; i + 1 < face_halfedges.size(); i += 2) {
			HalfedgeRef h0 = face_halfedges[i];
			HalfedgeRef h1 = inner_edges.at(i/2)->halfedge;
			HalfedgeRef h2 = inner_edges.at((i/2 == 0 ? inner_edges.size()-1 : i/2-1))->halfedge->twin;
			HalfedgeRef h3 = face_halfedges[(i == 0 ? face_halfedges.size()-1 : i-1)];

			//connect halfedges around the face:
			h0->next = h1;
			h1->next = h2;
			h2->next = h3;
			assert(h3->next == h0); //already connected and part of the face

			//connect halfedges to the face:
			if (i == 0) {
				//first face re-uses f:
				assert(f->halfedge == h0);
				assert(h0->face == f);
				h1->face = f;
				h2->face = f;
				assert(h3->face == f);
			} else {
				//other faces made fresh:
				FaceRef n = emplace_face(false);
				n->halfedge = h0;
				h0->face = n;
				h1->face = n;
				h2->face = n;
				h3->face = n;
			}
		}

	}

	//--------------------------
	//update positions for vertices
	for (VertexRef v = vertices.begin(); v != std::next(last_old_vertex); ++v) {
		v->position = vertex_positions.at(v);
	}

	{ //PARANOIA: sanity check:
		auto ret = validate();
		if (ret) {
			warn("After subdivide, validate says:\n  %s", ret.value().second.c_str());
		}
		assert(!ret && "subdivide helper should never break topology");
	}
}

/*
 * flip_orientation: flip direction of all halfedges
 *
 * works on all valid meshes.
 */
void Halfedge_Mesh::flip_orientation() {

	//store new h->vertex and v->halfedge pointers:
	std::unordered_map<Halfedge const *, VertexRef> he_to_v;
	std::unordered_map<Vertex const *, HalfedgeRef> v_to_he;
	for (auto &he : halfedges) {
		he_to_v[&he] = he.twin->vertex;
	}
	for (auto &v : vertices) {
		v_to_he[&v] = v.halfedge->twin;
	}

	//reverse all face loops:
	for (auto &face : faces) {
		//read off halfedges around face:
		std::vector<HalfedgeRef> hs;
		std::vector<Vec2> uvs;
		std::vector<Vec3> normals;

		HalfedgeRef h = face.halfedge;
		do {
			hs.emplace_back(h);
			uvs.emplace_back(h->corner_uv);
			normals.emplace_back(h->corner_normal);
			h = h->next;
		} while (h != face.halfedge);

		//reverse face ordering:
		for (uint32_t i = 0; i < hs.size(); ++i) {
			hs[(i+1)%hs.size()]->next = hs[i];
			hs[i]->corner_uv = uvs[(i+1)%hs.size()];
			hs[i]->corner_normal = normals[(i+1)%hs.size()];
		}
	}

	//update h->vertex and v->halfedge pointers:
	for (auto &he : halfedges) {
		he.vertex = he_to_v.at(&he);
	}
	for (auto &v : vertices) {
		v.halfedge = v_to_he.at(&v);
	}
}

/*
 * set_corner_normals: compute face-corner normals based on `sharp` flag and smoothing threshold
 *
 * works on all valid meshes.
 */
void Halfedge_Mesh::set_corner_normals(float threshold) {
	//first, figure out which edges to consider sharp for this operation:
	std::unordered_set< Edge const * > sharp_edges;
	sharp_edges.reserve(edges.size());

	//all edges between boundary and non-boundary get marked sharp regardless of mode:
	for (auto const &edge : edges) {
		if (edge.halfedge->face->boundary != edge.halfedge->twin->face->boundary) {
			sharp_edges.emplace(&edge);
		}
	}

	if (threshold >= 180.0f) {
		//"smooth mode" -- all other edges are considered smooth
	} else {
		//"flat mode" / "auto mode" -- any edges which are marked sharp or have face angle <= threshold get marked sharp:
		float cos_threshold = std::cos( Radians( std::clamp(threshold, 0.0f, 180.0f) ) );
		if (threshold <= 0.0f) cos_threshold = 2.0f; //make sure everything is sharp
		for (auto const &edge : edges) {
			//get adjacent halfedges:
			HalfedgeRef h1 = edge.halfedge;
			HalfedgeRef h2 = h1->twin;
			if (h1->face->boundary || h2->face->boundary) {
				//don't care about edges boundary-boundary, and inside-boundary already marked.
				//thus: nothing to do here
			} else if (edge.sharp) {
				//flagged as sharp, so mark it sharp:
				sharp_edges.emplace(&edge);
			} else {
				//inside-inside edge, non-marked, check angle:
				Vec3 n1 = h1->face->normal();
				Vec3 n2 = h2->face->normal();
				float cos = dot(n1,n2);
				if (cos <= cos_threshold) {
					//treat as sharp:
					sharp_edges.emplace(&edge);
				}
			}
		}
	}

	//clear current corner normals:
	for (auto h = halfedges.begin(); h != halfedges.end(); ++h) {
		h->corner_normal = Vec3{0.0f, 0.0f, 0.0f};
	}

	//now circulate all vertices to set normals:
	for (auto const &v : vertices) {
		//get halfedge leaving this vertex:
		HalfedgeRef begin = v.halfedge;
		assert(&*begin->vertex == &v);

		//circulate begin until it is at a sharp edge (thus, the next corner starts a smoothing group):
		do {
			if (sharp_edges.count(&*begin->edge)) break;
			begin = begin->twin->next;
		} while (begin != v.halfedge); //could be all one big happy smoothing group

		//store all corners around the vertex:
		struct Corner {
			HalfedgeRef in; //halfedge pointing to v
			HalfedgeRef out; //halfedge pointing away from v
			Vec3 weighted_normal; //face normal at corner, weighted... somehow (see below)
		};

		std::vector< std::vector< Corner > > groups;
		HalfedgeRef h = begin;
		do {
			//start a new smoothing group on sharp edges (or at the very first edge):
			if (h == begin || sharp_edges.count(&*h->edge)) {
				groups.emplace_back();
			}
			//add corner after h to current smoothing group:
			Corner corner;
			corner.in = h->twin;
			corner.out = h->twin->next;
			assert(corner.in->face == corner.out->face); //PARANOIA
			{ //compute some sort of weighted normal:
				assert(&*corner.in->vertex != &v);
				assert(&*corner.in->twin->vertex == &v);
				Vec3 from = corner.in->vertex->position - v.position;

				assert(&*corner.out->vertex == &v);
				assert(&*corner.out->twin->vertex != &v);
				Vec3 to = corner.out->twin->vertex->position - v.position;
				/*
				//basic area weighting (weird with non-flat faces and reflex vertices):
				corner.weighted_normal = cross(to - v.position, from - v.position);
				*/
				/*//sort sort of angle weighting thing -- this never works as well as one would hope:
				//...still needs work for reflex angles also
				float angle = std::atan2(cross(from,to).norm(), dot(from, to));
				corner.weighted_normal = angle * corner.in->face->normal();
				*/
				//some other sort of slightly fancy area weighting:
				corner.weighted_normal = cross(to - v.position, from - v.position).norm() * corner.in->face->normal();
			}
			groups.back().emplace_back(corner);

			//advance h:
			h = h->twin->next;
		} while (h != begin);

		//compute weighted normals per-corner:
		for (auto const &group : groups) {
			assert(!group.empty());
			if (group[0].in->face->boundary) {
				//boundary group.
				//PARANOIA:
				for (auto const &corner : group) {
					assert(corner.in->face->boundary);
				}
				//no need for normals on boundary corners
				continue;
			}
			//compute weighted normal:
			Vec3 sum = Vec3{0.0f, 0.0f, 0.0f};
			for (auto const &corner : group) {
				sum += corner.weighted_normal;
			}
			//normalize:
			sum = sum.unit();
			//assign to all corners in group:
			for (auto const &corner : group) {
				assert(&*corner.out->vertex == &v);
				corner.out->corner_normal = sum;
			}
		}
	}

	//normals computed!
}

/*
 * set_corner_uvs_per_face: set uv coordinates to map texture per-face
 */
void Halfedge_Mesh::set_corner_uvs_per_face() {
	//clear existing UVs:
	for (auto &halfedge : halfedges) {
		halfedge.corner_uv = Vec2(0.0f, 0.0f);
	}
	
	//set UVs per-face:
	for (auto const &face : faces) {
		if (face.boundary) continue;

		//come up with a plane perpendicular-ish to the face:
		Vec3 n = face.normal();
		Vec3 p1;
		if (std::abs(n.x) < std::abs(n.y) && std::abs(n.x) < std::abs(n.z)) {
			p1 = Vec3(1.0f, 0.0f, 0.0f);
		} else if (std::abs(n.y) < std::abs(n.z)) {
			p1 = Vec3(0.0f, 1.0f, 0.0f);
		} else {
			p1 = Vec3(0.0f, 0.0f, 1.0f);
		}
		p1 = (p1 - dot(p1, n) * n).unit();
		Vec3 p2 = cross(n, p1);

		//find bounds of face on plane:
		Vec2 min = Vec2(std::numeric_limits< float >::infinity(), std::numeric_limits< float >::infinity());
		Vec2 max = Vec2(-std::numeric_limits< float >::infinity(), -std::numeric_limits< float >::infinity());
		HalfedgeRef v = face.halfedge;
		do {
			Vec2 pt = Vec2(dot(p1, v->vertex->position), dot(p2, v->vertex->position));
			min = hmin(min, pt);
			max = hmax(max, pt);
			v = v->next;
		} while (v != face.halfedge);

		//set corner uvs based on position within bounds:
		do {
			Vec2 pt = Vec2(dot(p1, v->vertex->position), dot(p2, v->vertex->position));
			v->corner_uv = Vec2(
				(pt.x - min.x) / (max.x - min.x),
				(pt.y - min.y) / (max.y - min.y)
			);
			v = v->next;
		} while (v != face.halfedge);
	}
}

/*
 * set_corner_uvs_project: set uv coordinates to map texture by projection to a plane
 */
void Halfedge_Mesh::set_corner_uvs_project(Vec3 origin, Vec3 u_axis, Vec3 v_axis) {

	u_axis /= u_axis.norm_squared();
	v_axis /= v_axis.norm_squared();

	for (auto &halfedge : halfedges) {
		if (halfedge.face->boundary) {
			halfedge.corner_uv = Vec2(0.0f, 0.0f);
		} else {
			halfedge.corner_uv = Vec2(
				dot(halfedge.vertex->position - origin, u_axis),
				dot(halfedge.vertex->position - origin, v_axis)
			);
		}
	}
}
