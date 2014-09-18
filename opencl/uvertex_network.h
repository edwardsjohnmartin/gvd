#ifndef __UVERTEX_NETWORK_H__
#define __UVERTEX_NETWORK_H__

//----------------------------------------------------------------------
// Unmanaged vertex network struct
//----------------------------------------------------------------------

#include "./vec.h"
#include "./vertex.h"
#include "./geometry.h"

NAMESPACE_OCT_BEGIN

//----------------------------------------
// struct UVertexNetwork
//----------------------------------------
typedef struct {
  __GLOBAL__ int* header;
  // int num_vertices;
  __GLOBAL__ Vertex* vertices;

  // int num_cpoints;
  __GLOBAL__ GeomPoint* cpoints;
} UVertexNetwork;

//------------------------------------------------------------
// Initialization
//------------------------------------------------------------

void Initialize(UVertexNetwork* vn);

inline UVertexNetwork make_uvertex_network(
    __GLOBAL__ int* header,
    // int num_vertices, __GLOBAL__ Vertex* vertices,
    // int num_cpoints, __GLOBAL__ GeomPoint* cpoints) {
    __GLOBAL__ Vertex* vertices,
    __GLOBAL__ GeomPoint* cpoints) {
  // UVertexNetwork uvn = { header, num_vertices, vertices, num_cpoints, cpoints };
  UVertexNetwork uvn = { header, vertices, cpoints };
  return uvn;
}

//------------------------------------------------------------
// Header functions
//------------------------------------------------------------

inline int NumVerticesFromHeader(__GLOBAL__ int* header) {
  return header[0];
}

inline int NumCPointsFromHeader(__GLOBAL__ int* header) {
  return header[1];
}

inline int NumVertices(UVertexNetwork vn) {
  return NumVerticesFromHeader(vn.header);
  // return vn.header[0];
  // return vn.num_vertices;
}

inline int NumCPoints(UVertexNetwork vn) {
  return NumCPointsFromHeader(vn.header);
  // return vn.header[1];
  // return vn.num_cpoints;
}

inline int IncNumVertices(UVertexNetwork* vn) {
  // return vn->header[0]++;
#ifdef OPEN_CL
  return atomic_inc(&vn->header[0]);
#else
  return vn->header[0]++;
#endif
  // return vn->num_vertices++;
}

inline int IncNumCPoints(UVertexNetwork* vn) {
  // return vn->header[1]++;
#ifdef OPEN_CL
  return atomic_inc(&vn->header[1]);
#else
  return vn->header[1]++;
#endif
  // return vn->num_cpoints++;
}

//------------------------------------------------------------
// Modifiers
//------------------------------------------------------------

// Constructs 2^D new cells by subdividing the cell base_vi.
//
// Returns the new subcells with lvi indexing (see bit.h)
// slvi2vi is filled with full subdivided cell indices if non-null.
// void Subdivide(
//     const int base_vi, int* subcells_, int* slvi2vi, UVertexNetwork* vn);
void Subdivide(
    const int base_vi, int* slvi2vi, UVertexNetwork* vn);

void Subdivide_A(const int base_vi, UVertexNetwork vn);
// void Subdivide_B(const int base_vi, UVertexNetwork vn);
// void Subdivide_C(const int base_vi, UVertexNetwork vn);

bool SetClosestPoint(const int vi, const int candidate, UVertexNetwork* vn);

// Returns the new index
inline int CreateNewClosestPoint(
    const intn p, const int label, UVertexNetwork* vn) {
  // const int p_idx = vn->_num_cpoints++;
  const int p_idx = IncNumCPoints(vn);
  vn->cpoints[p_idx].p = p;
  vn->cpoints[p_idx].l = label;
  return p_idx;
}

//------------------------------------------------------------
// Accessor functions
//------------------------------------------------------------

inline int Neighbor(const int vi, const Direction d, UVertexNetwork vn) {
  return v_neighbor_index(d, vn.vertices[vi]);
}

int FindNeighbor(const int a_vi, const Direction d,
                 const level_t level, UVertexNetwork vn);

inline int Label(const int vi, UVertexNetwork vn) {
  const int p_idx = v_closest_point(vn.vertices[vi]);
  if (p_idx == -1) return -1;
  return vn.cpoints[p_idx].l;
}

// Return a copy.  Returning a reference is dangerous since
// the vector could be reallocated upon resize.
inline intn ClosestPoint(const int vi, UVertexNetwork vn) {
  const int p_idx = v_closest_point(vn.vertices[vi]);
  assert(p_idx != -1);
  return vn.cpoints[p_idx].p;
}

inline int ClosestPointIndex(const int vi, UVertexNetwork vn) {
  return v_closest_point(vn.vertices[vi]);
}

// Returns the distance from the vertex to its closest geometry point.
// Returns -1 if there is no closest geometry point.
inline int Dist(const int vi, UVertexNetwork vn) {
  const int p_idx = v_closest_point(vn.vertices[vi]);
  if (p_idx == -1) return -1;
  const intn a = vn.cpoints[p_idx].p;
  const intn b = vn.vertices[vi].position;
  return (int)length(convert_floatn(a-b));
}

inline uchar IsBase(const int vi, UVertexNetwork vn) {
  return v_is_base(vn.vertices[vi]);
}

// These only apply if the vertex is the base of a leaf cell
inline level_t CellLevel(const int vi, UVertexNetwork vn) {
  assert(v_is_base(vn.vertices[vi]));
  return v_cell_level(vn.vertices[vi]);
}

inline index_t CellWidth(const int vi, UVertexNetwork vn) {
  assert(v_is_base(vn.vertices[vi]));
  return Level2CellWidth(CellLevel(vi, vn));
}

// Returns true if cells vi0 and vi1 are incident to each other
bool AreAdjacent(const int vi0, const int vi1, UVertexNetwork vn);

// Returns true if vi is incident to the cell described by base_vi.
// Check is done by looking at the positions.
bool IsIncident(const int vi, const int base_vi, UVertexNetwork vn);

// This function is linear!
int NumCells(UVertexNetwork vn);

// nbrs must be of size at least kNumIncidentVertices
// Returns the number of neighbors
int GetNeighborCells(
    const int base_vi, __GLOBAL__ int* nbrs, UVertexNetwork vn);

// Finds the index of the vertex at the given slvi index.  This is
// inefficient.
int Slvi2Vi(const int base_vi, const int slvi, const int level,
            UVertexNetwork vn);

void BuildSlvi2Vi(const int base_vi, int* slvi2vi, const int level,
                  UVertexNetwork vn);

NAMESPACE_OCT_END

#endif
