/*******************************************************
 ** Generalized Voronoi Diagram Project               **
 ** Copyright (c) 2015 John Martin Edwards            **
 ** Scientific Computing and Imaging Institute        **
 ** 72 S Central Campus Drive, Room 3750              **
 ** Salt Lake City, UT 84112                          **
 **                                                   **
 ** For information about this project contact        **
 ** John Edwards at                                   **
 **    edwardsjohnmartin@gmail.com                    **
 ** or visit                                          **
 **    sci.utah.edu/~jedwards/research/gvd/index.html **
 *******************************************************/

#include <stdio.h>
#include <limits.h>

#include "./cl_octree.h"
#include "./distance3.h"
#include "./cl_ambiguous.h"

NAMESPACE_OCT_BEGIN

//--------------------------------------------------
//--------------------------------------------------
// Declarations
//--------------------------------------------------
//--------------------------------------------------

bool IsSubdividable(
    const int vi,
    const int filter,
    const UVertexNetwork vn,
    const __GLOBAL__ uchar* to_subdivide);

//--------------------------------------------------
//--------------------------------------------------
// Definitions
//--------------------------------------------------
//--------------------------------------------------

//------------------------------------------------------------------------------
// FindToSubdivide code
//------------------------------------------------------------------------------

//------------------------------------------------------------
// NeighborVisitor struct and callback
// These are used in FindToSubdivide2
//------------------------------------------------------------

inline void CheckNeighborCompatibility(
    int base_vi, int nbr_vi,
    UVertexNetwork uvn, Vi2Geometries vi2geometries,
    __GLOBAL__ uchar* to_subdivide) {
  if (to_subdivide[base_vi]) {
    // Already subdividing
    return;
  }

  if (!gcell_is_empty(base_vi, vi2geometries) &&
      !gcell_is_empty(nbr_vi, vi2geometries)) {
    Geometries base_geoms = get_geometries(base_vi, vi2geometries);
    Geometries nbr_geoms = get_geometries(nbr_vi, vi2geometries);
    assert(!cell_has_multiple_labels(base_geoms));
    const int base_label = g_label(get_geometry(0, base_geoms));
    const int nbr_label = g_label(get_geometry(0, nbr_geoms));
    if (base_label != nbr_label ||
        cell_has_multiple_labels(nbr_geoms)) {
      // Make buffer between differently-labeled cells
      to_subdivide[base_vi] = true;
    }
  }

  if (!to_subdivide[base_vi]) {
    assert(CellLevel(nbr_vi, uvn) - CellLevel(base_vi, uvn) <= kGradation);
    if (CellLevel(nbr_vi, uvn) - CellLevel(base_vi, uvn) == kGradation) {
      if (to_subdivide[nbr_vi]) {
        // Gradation
        to_subdivide[base_vi] = true;
      }
    }
  }
}

//------------------------------------------------------------
// FindToSubdivide1
// Checks if a cell intersects with multiple geometries
//------------------------------------------------------------
void FindToSubdivide1(
    int vi,
    __GLOBAL__ int* nv_header, __GLOBAL__ Vertex* vertices,
    __GLOBAL__ int* vi2geometries_array,
    int max_level,
    __GLOBAL__ uchar* to_subdivide) {
  to_subdivide[vi] = false;

  UVertexNetwork uvn = make_uvertex_network(nv_header, vertices, 0);

  if (!IsBase(vi, uvn)) {
    return;
  }
  assert(CellLevel(vi, uvn) <= max_level);
  if (CellLevel(vi, uvn) == max_level) {
    return;
  }

  Vi2Geometries vi2geometries = make_vi2geometries(vi2geometries_array);
  if (gcell_is_empty(vi, vi2geometries)) {
    return;
  }
  Geometries geometries =
      get_geometries(vi, vi2geometries);

  const uchar subdivide = 
      v_is_base(uvn.vertices[vi]) &&
      !cell_is_empty(geometries) &&
      cell_has_multiple_labels(geometries);
  to_subdivide[vi] = subdivide;
}

//------------------------------------------------------------
// FindToSubdivide2
// Checks if a cell's neighbors either have different labels
// (in order to have a buffer region) or are subdivided more
// than kGradation than this cell (in order to have bounded
// gradation).
//------------------------------------------------------------
void FindToSubdivide2(
    int vi,
    __GLOBAL__ int* nv_header, __GLOBAL__ Vertex* vertices,
    __GLOBAL__ int* vi2nbrs,
    __GLOBAL__ int* vi2geometries_array,
    int max_level,
    __GLOBAL__ uchar* to_subdivide,
    __GLOBAL__ uchar* changed) {
  // FindSubdivide1 decided already that this cell must be subdivided.
  // We may need to set divide of neighbors.
  if (to_subdivide[vi])
    return;

  UVertexNetwork uvn = make_uvertex_network(
      nv_header, vertices, 0);

  if (!IsBase(vi, uvn)) {
    return;
  }
  assert(CellLevel(vi, uvn) <= max_level);
  if (CellLevel(vi, uvn) == max_level) {
    return;
  }

  Vi2Geometries vi2geometries = make_vi2geometries(vi2geometries_array);
#ifdef OPEN_CL
  __GLOBAL__ int* nbrs = vi2nbrs + vi * kNumIncidentCells;
  int num_nbrs = 0;
  while (nbrs[num_nbrs] > -1 && num_nbrs < kNumIncidentCells) ++num_nbrs;
#else
  int nbrs[kNumIncidentCells];
  const int num_nbrs = GetNeighborCells(vi, nbrs, uvn);
#endif
  assert(num_nbrs <= kNumIncidentCells);
  for (int i = 0; i < num_nbrs; ++i) {
    CheckNeighborCompatibility(vi, nbrs[i], uvn, vi2geometries, to_subdivide);
  }
  if (to_subdivide[vi]) {
    to_subdivide[vi] = true;
    *changed = true;
  }
}

//------------------------------------------------------------------------------
// IsSubdividable
// Checks only to_subdivide and filter.  Does not check octree level.
//------------------------------------------------------------------------------
bool IsSubdividable(
    const int vi,
    const int filter,
    const UVertexNetwork vn,
    const __GLOBAL__ uchar* to_subdivide) {

  if (!to_subdivide[vi])
    return false;

  // Use filter to see if we can subdivide
  const intn position = vn.vertices[vi].position;
  const intn pos = position / CellWidth(vi, vn);
  for (int dim = 0, f = filter; dim < DIM; ++dim, f /= 2) {
    if (intn_comp(dim, pos) % 2 != f % 2) {
      return false;
    }
  }

  return true;
}

//------------------------------------------------------------------------------
// SubdivideCell
// filter describes which cell should be subdivided if adjacent cells
// of identical width are both to be subdivided.
//------------------------------------------------------------------------------
void SubdivideCell(
    const int vi,
    const int filter,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    const __GLOBAL__ uchar* to_subdivide) {
  UVertexNetwork vn = make_uvertex_network(
      vn_header, vertex_array, 0);

  if (!IsSubdividable(vi, filter, vn, to_subdivide))
    return;

  // Subdivide
  // Subdivide(vi, 0, &vn);
  Subdivide_A(vi, vn);
  // Subdivide_B(vi, vn);
  // Subdivide_C(vi, vn);
}

void SubdivideCell_A(
    const int vi,
    const int filter,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    const __GLOBAL__ uchar* to_subdivide) {
  UVertexNetwork vn = make_uvertex_network(
      vn_header, vertex_array, 0);
  if (!IsSubdividable(vi, filter, vn, to_subdivide))
    return;
  Subdivide_A(vi, vn);
}

// void SubdivideCell_B(
//     const int vi,
//     const int filter,
//     __GLOBAL__ int* vn_header,
//     __GLOBAL__ Vertex* vertex_array,
//     const __GLOBAL__ uchar* to_subdivide) {
//   UVertexNetwork vn = make_uvertex_network(
//       vn_header, vertex_array, 0);
//   if (!IsSubdividable(vi, filter, vn, to_subdivide))
//     return;
//   Subdivide_B(vi, vn);
// }

// void SubdivideCell_C(
//     const int vi,
//     const int filter,
//     __GLOBAL__ int* vn_header,
//     __GLOBAL__ Vertex* vertex_array,
//     const __GLOBAL__ uchar* to_subdivide) {
//   UVertexNetwork vn = make_uvertex_network(
//       vn_header, vertex_array, 0);
//   if (!IsSubdividable(vi, filter, vn, to_subdivide))
//     return;
//   Subdivide_C(vi, vn);
// }

void ClipGeometriesAfterSubdivide(
    const int vi,
    const int filter,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    const __GLOBAL__ uchar* to_subdivide,
    __GLOBAL__ int* vi2geometries_array,
    int num_gvertices, __GLOBAL__ intn* gvertices,
    int num_goffsets, __GLOBAL__ int* goffsets) {

  UVertexNetwork vn = make_uvertex_network(
      vn_header, vertex_array, 0);

  if (!IsSubdividable(vi, filter, vn, to_subdivide))
    return;

  const index_t cur_width = Level2CellWidth(CellLevel(vi, vn));

  Vi2Geometries vi2geometries = make_vi2geometries(vi2geometries_array);
  GeomVertices geom_vertices = make_geom_vertices(
      num_gvertices, gvertices, num_goffsets, goffsets);

  if(!gcell_is_empty(vi, vi2geometries)) {
    Geometries geometries = get_geometries(vi, vi2geometries);
    // Clip geometries
    for (int lvi = 1; lvi < (1<<DIM); ++lvi) {
      const int sub_vi = vn.vertices[vi].corners[lvi];
      if (sub_vi == -1) {
        printf("Error: SubdivideCell vi = %d, lvi = %d, sub_vi = -1\n",
               vi, lvi);
        return;
      }
      allocate(geometries, sub_vi, vi2geometries);
      const intn center = vn.vertices[sub_vi].position + (cur_width>>1);
      ClipGeometries(geometries, get_geometries(sub_vi, vi2geometries),
                     center, cur_width, geom_vertices);
    }
    const intn center = vn.vertices[vi].position + (cur_width>>1);
    ClipGeometries(geometries, get_geometries(vi, vi2geometries),
                   center, cur_width, geom_vertices);
  }
}

void ComputeNonEmptyVertexDistances(
    const int base_vi,
    const int lvi,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    __GLOBAL__ GeomPoint* cpoints_array,
    __GLOBAL__ int* vi2geometries_array,
    int num_gvertices, __GLOBAL__ intn* gvertices,
    int num_goffsets, __GLOBAL__ int* goffsets) {

#ifdef OCT2D
  throw std::logic_error("GPU-based geometry clipping not supported in 2D. "
                         "Run with --cpu option.");
#else
  UVertexNetwork vn = make_uvertex_network(
      vn_header, vertex_array, cpoints_array);
  Vi2Geometries vi2g = make_vi2geometries(vi2geometries_array);
  GeomVertices geom_vertices = make_geom_vertices(
      num_gvertices, gvertices, num_goffsets, goffsets);

  if (!IsBase(base_vi, vn) || gcell_is_empty(base_vi, vi2g))
    return;

  // Only compute distance for vertex at local vertex index lvi
  // in relation to base_vi.
  const int vi = vn.vertices[base_vi].corners[lvi];
  const intn v_point = vn.vertices[vi].position;
  Geometries geometries = get_geometries(base_vi, vi2g);
  const int n = g_n(geometries);
  Distance3i min_dist = make_dist3(INT_MAX, make_int3(0, 0, 0));
  Geometry min_geometry;
  for (int i = 0; i < n; ++i) {
    Geometry geometry = get_geometry(i, geometries);
    __GLOBAL__ intn* gverts =
        get_geom_vertices(g_label(geometry), geom_vertices);
    Distance3i dist = distance_geometry(v_point, geometry, gverts);
    if (dist.d < min_dist.d) {
      min_dist = dist;
      min_geometry = geometry;
    }
  }

  assert(min_dist.d < INT_MAX);
  const int label = g_label(min_geometry);
  const intn point = min_dist.p;
  const int dist = min_dist.d;

  // If the vertex network already has a closest point for vi, then use
  // that.  This is because no other vi will be using that closest point.
  // Otherwise, create and new closest point.
  int cp_idx = ClosestPointIndex(vi, vn);
  if (cp_idx == -1) {
    cp_idx = CreateNewClosestPoint(point, label, &vn);
    SetClosestPoint(vi, cp_idx, &vn);
  } else if (Dist(vi, vn) > dist) {
    SetClosestPoint(vi, cp_idx, &vn);
  }

#endif
}

void PullDistances(
    const int vi,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    __GLOBAL__ GeomPoint* cpoints_array,
    __GLOBAL__ uchar* changed) {
  UVertexNetwork vn = make_uvertex_network(
      vn_header, vertex_array, cpoints_array);

  for (int axis = 0; axis < DIM; ++axis) {
    for (int pos = 0; pos < 2; ++pos) {
      const Direction dir = DirectionFromAxis(axis, pos);
      const int n_vi = Neighbor(vi, dir, vn);
      if (n_vi > -1) {
        const int n_cp_idx = ClosestPointIndex(n_vi, vn);
        if (n_cp_idx > -1) {
          if (SetClosestPoint(vi, n_cp_idx, &vn))
            *changed = true;
        }
      }
    }
  }
}

void FindAmbiguous(
    int vi,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    __GLOBAL__ GeomPoint* cpoints_array,
    int max_level,
    __GLOBAL__ uchar* ambiguous,
    __GLOBAL__ int* count) {
  ambiguous[vi] = 0;
  if (vertex_array[vi].level < max_level &&
      vertex_array[vi].base &&
      IsAmbiguous3(vi, vn_header, vertex_array, cpoints_array)) {
    ambiguous[vi] = 1;
#ifdef OPEN_CL
    atomic_inc(count);
#else
    (*count)++;
#endif
  }
}

NAMESPACE_OCT_END
