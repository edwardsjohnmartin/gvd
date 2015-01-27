#include <cstring>

#include "./geometry.h"
#include "./tribox.h"

NAMESPACE_OCT_BEGIN

// Suppose m cells are candidates for subdivision.  Then each of those
// cells must have 2^D * geometries_size(cell[i]) space to write
// new geometries to.
//
// candidates are the cells that are candidates for subdivision.
int compute_sparse_vi2geometries_size(
    const Vi2Geometries dense_geoms,
    __GLOBAL__ uchar* candidates) {
  const int n = g_N(dense_geoms);
  int size = vi2g_offset + n;
  for (int vi = 0; vi < n; ++vi) {
    int cell_size;
    if (gcell_is_empty(vi, dense_geoms)) {
      cell_size = 1;
    } else {
      cell_size = geometries_size(get_geometries(vi, dense_geoms));
    }
    if (candidates[vi]) {
      size += cell_size * kNumNewVertices;
    } else {
      size += cell_size;
    }
  }
  return size;
}

// Suppose m cells are candidates for subdivision.  Then the number of
// new cells will be M = m * (2^D - 1).  Each of these new cells is allocated
// space at the end of the array, with enough space to have the same
// geometries as its parent cell.
//
// candidates is of size dense_geoms.n (the number of vertices) and
// candidates[vi] is true if vi is a candidate for subdivision.
Vi2Geometries make_sparse_vi2geometries(
    const Vi2Geometries dense_geoms,
    __GLOBAL__ uchar* candidates,
    __GLOBAL__ int* array,
    int array_size) {
  const int dense_N = g_N(dense_geoms);

  // Find m
  int m = 0;
  for (int vi = 0; vi < dense_N; ++vi) {
    if (candidates[vi]) {
      ++m;
    }
  }
  // M is the number of new vertices to be added
  const int M = m * kNumNewVertices;
  const int sparse_N = dense_N + M;
  
  // Vi2Geometries sparse_geoms =
  //     { array, dense_geoms.free_offset+M, array, array_size };
  Vi2Geometries sparse_geoms =
      { array+vi2g_offset, array };
  array[0] = array_size;
  array[1] = g_free_offset(dense_geoms)+M;
  // Iterate through the old geometries incrementing their offsets by M.
  int gsize = 0;
  for (int vi = 0; vi < dense_N; ++vi) {
    const int dense_offset = dense_geoms.offsets[vi];
    if (dense_offset == -1) {
      sparse_geoms.offsets[vi] = -1;
    } else {
      sparse_geoms.offsets[vi] = dense_offset + M;
      gsize += geometries_size(get_geometries(vi, dense_geoms));
    }
  }
  // assert(gsize == dense_geoms.free_offset-dense_N);
  assert(gsize == g_free_offset(dense_geoms)-(dense_N+vi2g_offset));
  // Copy the existing geometries
#ifdef OPEN_CL
  const int num = g_free_offset(dense_geoms)-(vi2g_offset+dense_N);
  for (int i = 0; i < num; ++i) {
    sparse_geoms.array[vi2g_offset+sparse_N+i] =
        dense_geoms.array[vi2g_offset+dense_N+i];
  }
#else
  memcpy(sparse_geoms.array+vi2g_offset+sparse_N,
         dense_geoms.array+vi2g_offset+dense_N,
         (g_free_offset(dense_geoms)-(vi2g_offset+dense_N)) * sizeof(int));
#endif

  // All remaining vertices are empty
  for (int vi = dense_N; vi < sparse_N; ++vi) {
    sparse_geoms.offsets[vi] = -1;
  }
  // sparse_geoms.free_offset = sparse_N + (dense_geoms.free_offset-dense_N);
  // assert(sparse_geoms.free_offset <= array_size);
  assert(g_free_offset(sparse_geoms) <= array_size);
  return sparse_geoms;
}

// Given a sparse array of octree geometries, defragment (condense)
// so that all free space is at the end.
// N is the current number of vertices.  Truncate the number of
// vertices in geoms to N.
Vi2Geometries condense_vi2geometries(
    const Vi2Geometries sparse, const int dense_N, __GLOBAL__ int* array,
    int array_size) {
  // Vi2Geometries dense = { array, array_size, array, array_size };
  Vi2Geometries dense = { array+vi2g_offset, array };
  array[0] = array_size;
  array[1] = array_size;
  assert(dense_N <= g_N(sparse));
  int dense_offset = vi2g_offset + dense_N;
  for (int vi = 0; vi < dense_N; ++vi) {
    if (gcell_is_empty(vi, sparse) ||
        g_n(get_geometries(vi, sparse)) == 0) {
      dense.offsets[vi] = -1;
    } else {
      dense.offsets[vi] = dense_offset;
      Geometries sparse_geoms = get_geometries(vi, sparse);
      Geometries dense_geoms =
          { dense.array+dense.offsets[vi], dense.array+dense.offsets[vi] };

      const int n = g_n(sparse_geoms);
      int geom_offset = n;
      g_set_n(n, dense_geoms);
      for (int gi = 0; gi < n; ++gi) {
        dense_geoms.offsets[gi] = geom_offset;
        Geometry sparse_geom = get_geometry(gi, sparse_geoms);
        Geometry dense_geom = get_geometry(gi, dense_geoms);
        const int size = geometry_size(sparse_geom);
#ifdef OPEN_CL
        for (int k = 0; k < size; ++k) {
          dense_geom.array[k] = sparse_geom.array[k];
        }
#else
        memcpy(dense_geom.array,
               sparse_geom.array,
               size * sizeof(int));
#endif
        geom_offset += size;
      }
      dense_offset += geometries_size(dense_geoms);
    }
  }
  // dense.free_offset = dense_offset;
  g_set_free_offset(dense_offset, dense);
  // assert(dense.free_offset <= array_size);
  assert(g_free_offset(dense) <= array_size);
  return dense;
}

// Allocates space to copy, at most, the number of geometries in source,
// to the cell at vi.
void allocate(const Geometries source, const int vi,
              Vi2Geometries vi2geometries) {
  assert(vi < g_N(vi2geometries));
  assert(vi2geometries.offsets[vi] == -1);

  const int size = geometries_size(source);
  // vi2geometries.offsets[vi] = g_free_offset(vi2geometries);
  // g_set_free_offset(g_free_offset(vi2geometries) + size, vi2geometries);
  vi2geometries.offsets[vi] = g_add_free_offset(size, vi2geometries);

  __GLOBAL__ int* target_array =
      vi2geometries.array + vi2geometries.offsets[vi];
#ifdef OPEN_CL
  for (int i = 0; i < size; ++i) {
    target_array[i] = source.array[i];
  }
#else
  memcpy(target_array, source.array, size * sizeof(int));
#endif
}

//------------------------------------------------------------
// ClipGeometries
//------------------------------------------------------------
void ClipGeometries(
    const Geometries source, Geometries target, const intn center,
    const index_t cell_width, const GeomVertices geom_vertices) {
  const bool old = false;

#ifdef OCT2D
  throw std::logic_error("GPU-based geometry clipping not supported in 2D. "
                         "Run with --cpu option.");
#else
  if (old) {
    // original
    const int n = g_n(source);
    for (int i = 0; i < n; ++i) {
      Geometry source_g = get_geometry(i, source);
      Geometry target_g = get_geometry(i, target);
      const int label = g_label(source_g);
      const __GLOBAL__ intn* vertices = get_geom_vertices(label, geom_vertices);
      const int m = g_m(source_g);
      g_set_m(0, target_g);
      // Loop through each face
      for (int j = 0; j < m; ++j) {
        Face t = source_g.faces[j];
        if (TriBoxOverlap(center, cell_width >> 1, vertices, t)) {
          // Add face to geometry
          const int target_m = g_m(target_g);
          target_g.faces[target_m] = t;
          g_set_m(target_m+1, target_g);
        }
      }
    }
  } else {
    // new
    const int n = g_n(source);
    // Count the number of non-empty geometries
    int target_n = 0;
    for (int i = 0; i < n; ++i) {
      Geometry source_g = get_geometry(i, source);
      const int label = g_label(source_g);
      const __GLOBAL__ intn* vertices = get_geom_vertices(label, geom_vertices);
      const int m = g_m(source_g);
      // Loop through each face
      bool empty = true;
      for (int j = 0; empty && j < m; ++j) {
        Face t = source_g.faces[j];
        empty = !TriBoxOverlap(center, cell_width >> 1, vertices, t);
      }
      if (!empty)
        ++target_n;
    }

    // Source
    //   |           |g0                   |g1                   |g2
    //    o0 o1 o2 o3 l0 m0 g00.xyz g01.xyz l1 m1 g10.xyz g11.xyz l2...
    // Target1
    //   |           |g1                   
    //    o1 o1 o2 o3 l1 m1 g10.xyz g11.xyz
    // Target2
    //   |  |g1                   
    //    o1 l1 m1 g10.xyz g11.xyz

    // Copy non-empty geometries
    assert(target_n <= n);
    int offset1 = n;
    int target_i = 0;
    for (int i = 0; target_i < target_n && i < n; ++i) {
      Geometry source_g = get_geometry(i, source);
      const int label = g_label(source_g);
      const __GLOBAL__ intn* vertices = get_geom_vertices(label, geom_vertices);
      const int m = g_m(source_g);
      target.array[target_i] = offset1;
      Geometry target_g = get_geometry(target_i, target);
      g_set_m(0, target_g);
      g_set_label(label, target_g);
      assert(g_label(target_g) == label);
      // Loop through each face
      for (int j = 0; j < m; ++j) {
        Face t = source_g.faces[j];
        if (TriBoxOverlap(center, cell_width >> 1, vertices, t)) {
          // Add face to geometry
          const int target_m = g_m(target_g);
          target_g.faces[target_m] = t;
          g_set_m(target_m+1, target_g);
        }
      }
      if (g_m(target_g) > 0) {
        offset1 += geometry_size(target_g);
        ++target_i;
      }
    }
    // We wrote everything to offset1.  Change offsets and block copy.
    const int size = offset1 - n;
    const int diff_n = n - target_n;
    for (int gi = 0; gi < target_n; ++gi) {
      target.array[gi] = target.array[gi] - diff_n;
    }
    // todo: this is inefficient
    for (int k = 0; k < size; ++k) {
      target.array[target_n+k] = target.array[n+k];
    }
    // Can't use memcpy because the memory overlaps.
    // memcpy(target.array+target_n,
    //        target.array+n,
    //        size);

    if (target_n == 0)
      target.array[0] = 0;

    assert(geometries_size(target) <= geometries_size(source));
  }

#endif
}

NAMESPACE_OCT_END
