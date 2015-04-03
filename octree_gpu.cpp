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

#include <sys/time.h>

#include "./octree_gpu.h"
#include "./geometry_cpp.h"
#include "./timer.h"
#include "./wall_timer.h"
#include "./octree.h"
#include "./search.h"
#include "./vertices_gpu_state.h"
#include "./ambiguous.h"
#include "./gpu.h"

#include "./opencl/geometry.h"
#include "./opencl/cl_octree.h"

using namespace std;

NAMESPACE_OCT_BEGIN

//------------------------------------------------------------------------------
// SubdivideCellGpu
//------------------------------------------------------------------------------
void SubdivideCellGpu_cpu(
    const int vi,
    UVertexNetwork& uvn,
    const uchar* to_subdivide,
    int* vi2geometries_array,
    int num_gvertices, __GLOBAL__ intn* gvertices,
    int num_goffsets, __GLOBAL__ int* goffsets) {

  if (!to_subdivide[vi])
    return;

  Vi2Geometries vi2geometries = make_vi2geometries(vi2geometries_array);
  GeomVertices geom_vertices = make_geom_vertices(
      num_gvertices, gvertices, num_goffsets, goffsets);

  const level_t cur_level = CellLevel(vi, uvn);
  const index_t cur_width = Level2CellWidth(cur_level);

  // Subdivide
  // int subcells[1<<DIM];
  // int* slvi2vi = 0;
  // Subdivide(vi, subcells, slvi2vi, &uvn);
  Subdivide(vi, 0, &uvn);
  int* subcells = uvn.vertices[vi].corners;

  if(!gcell_is_empty(vi, vi2geometries)) {
    Geometries geometries = get_geometries(vi, vi2geometries);
    // Clip geometries
    for (int lvi = 1; lvi < (1<<DIM); ++lvi) {
      const int sub_vi = subcells[lvi];
      allocate(geometries, sub_vi, vi2geometries);
      const intn center = uvn.vertices[sub_vi].position + (cur_width>>2);
      ClipGeometries(geometries, get_geometries(sub_vi, vi2geometries),
                     center, cur_width>>1, geom_vertices);
    }
    const intn center = uvn.vertices[vi].position + (cur_width>>2);
    ClipGeometries(geometries, get_geometries(vi, vi2geometries),
                   center, cur_width>>1, geom_vertices);
  }

  // *num_vertices = NumVertices(uvn);
}

// //------------------------------------------------------------------------------
// // SubdivideCellsGpu
// // Subdivide a level of the octree.
// // base_vis are the base vertices to subdivide.
// //------------------------------------------------------------------------------
// void SubdivideCellsGpu(
//     Vi2Geometries vi2geometries,
//     const uchar* to_subdivide,
//     UVertexNetwork& vertices,
//     const GeomVertices geom_vertices,
//     const OctreeOptions& o) {

//   // vertices.num_vertices will change during subdivisions,
//   // so cache a copy for iteration.
//   const int n = NumVertices(vertices);
//   for (int filter = 0; filter < (1<<DIM); ++filter) {
//     for (int vi = 0; vi < n; ++vi) {
//       SubdivideCell(
//           vi,
//           filter,
//           vertices.header, vertices.vertices,
//           to_subdivide);//,
//           // vi2geometries.array,
//           // geom_vertices.num_vertices, geom_vertices.vertices,
//           // geom_vertices.num_offsets, geom_vertices.offsets);
//       ClipGeometriesAfterSubdivide(
//           vi,
//           filter,
//           vertices.header, vertices.vertices,
//           to_subdivide,
//           vi2geometries.array,
//           geom_vertices.num_vertices, geom_vertices.vertices,
//           geom_vertices.num_offsets, geom_vertices.offsets);
//     }
//   }
// }

// //------------------------------------------------------------------------------
// // SubdivideCellsGpu
// // Subdivide a level of the octree.
// // base_vis are the base vertices to subdivide.
// //------------------------------------------------------------------------------
// void SubdivideCellsGpu_cpu(
//     Vi2Geometries vi2geometries,
//     const uchar* to_subdivide,
//     UVertexNetwork& vertices,
//     const GeomVertices geom_vertices,
//     const OctreeOptions& o) {

//   // vertices.num_vertices will change during subdivisions,
//   // so cache a copy for iteration.
//   const int n = NumVertices(vertices);
//   for (int vi = 0; vi < n; ++vi) {
//     SubdivideCellGpu_cpu(
//         vi,
//         vertices,
//         to_subdivide,
//         vi2geometries.array,
//         geom_vertices.num_vertices, geom_vertices.vertices,
//         geom_vertices.num_offsets, geom_vertices.offsets);
//   }
// }

//------------------------------------------------------------
// GetToSubdivideGpu
// Makes gpu calls.  Must retrieve array from gpu memory
// after calling this function.
//------------------------------------------------------------
void GetToSubdivideGpu(
    int num_vertices,
    const int max_level,
    Gpu& gpu) {

  int header[2] = { num_vertices, 0 };
  cl_event h_event = gpu.EnqueueWriteVNHeader(header);
  gpu.EnqueueFindToSubdivide1(num_vertices, max_level, GpuEvents(h_event));
  gpu.Finish();

  uchar changed = true;
  int changed_count = 0;
  while (changed) {
    changed_count++;
    changed = false;

    cl_event write_event =
        gpu.EnqueueWriteChanged(&changed);
    int header[2] = { num_vertices, 0 };
    cl_event he = gpu.EnqueueWriteVNHeader(header);
    cl_event gpu_event =
        gpu.EnqueueFindToSubdivide2(num_vertices, max_level,
                                    GpuEvents(write_event, he));
    gpu.EnqueueReadChanged(&changed, GpuEvents(gpu_event));
    gpu.Finish();
  }
  // cout << "Gradation required " << changed_count << " iterations" << endl;
}

void ReadFromGpu(
    MVertexNetwork& mvertices,
    shared_array<int>& vi2geometries_array,
    Gpu& gpu) {

  mvertices = make_mvertex_network();

  // Read vertices and cpoints from gpu
  gpu.EnqueueReadVNHeader(mvertices.header.get());
  gpu.Finish();
  // Reserve space
  // const int num_vertices = NumVerticesFromHeader(mvertices.header.get());
  // mvn_reserve(num_vertices, &mvertices);
  mvn_reserve(gpu.VerticesSize(), &mvertices);
  // const int num_cpoints = NumCPoints(mvertices);
  // mvn_cp_reserve(num_cpoints, &mvertices);
  mvn_cp_reserve(gpu.CPointsSize(), &mvertices);
  vi2geometries_array.reset(new int[gpu.Vi2GeometriesSize()]);
  // Read from gpu
  // gpu.EnqueueReadVertices(num_vertices, mvertices.vertices.get());
  gpu.EnqueueReadVertices(gpu.VerticesSize(), mvertices.vertices.get());
  // if (num_cpoints > 0)
  if (gpu.CPointsSize() > 0)
    // gpu.EnqueueReadCPoints(num_cpoints, mvertices.cpoints.get());
    gpu.EnqueueReadCPoints(gpu.CPointsSize(), mvertices.cpoints.get());
  gpu.EnqueueReadVi2Geometries(vi2geometries_array.get());
  gpu.Finish();
}

void WriteToGpu(
    MVertexNetwork& mvertices,
    shared_array<int>& vi2geometries_array,
    Gpu& gpu) {

  Vi2Geometries vi2geometries = make_vi2geometries(vi2geometries_array.get());

  gpu.CreateVertices(mvertices.vertex_array_capacity);
  if (mvertices.cpoint_array_capacity > 0)
    gpu.CreateCPoints(mvertices.cpoint_array_capacity);
  gpu.CreateVi2Geometries(g_array_size(vi2geometries));

  gpu.EnqueueWriteVNHeader(mvertices.header.get());
  gpu.EnqueueWriteVertices(NumVertices(mvertices), mvertices.vertices.get());
  if (mvertices.cpoint_array_capacity > 0)
    gpu.EnqueueWriteCPoints(NumCPoints(mvertices), mvertices.cpoints.get());
  gpu.EnqueueWriteVi2Geometries(vi2geometries.array);
  gpu.Finish();
}

// Initialize gpu memory
void Stage0(
    Vi2Geometries& vi2geometries,
    shared_array<int>& vi2geometries_array,
    GeomVertices& geom_vertices,
    MVertexNetwork& mvertices,
    Gpu& gpu,
    const OctreeOptions& o) {

  // Initialize gpu memory
  gpu.CreateVertices(NumVertices(mvertices));
  gpu.CreateGeomVertices(geom_vertices);
  gpu.CreateVi2Geometries(g_array_size(vi2geometries));
  gpu.EnqueueWriteVNHeader(mvertices.header.get());
  gpu.EnqueueWriteVertices(NumVertices(mvertices), mvertices.vertices.get());
  gpu.EnqueueWriteGeomVertices(geom_vertices);
  gpu.EnqueueWriteVi2Geometries(vi2geometries.array);
  gpu.Finish();
}

// Ensure the vertex array is large enough to accommodate subdivided vertices
void Stage1a(
    const int num_vertices,
    Gpu& gpu,
    const OctreeOptions& o,
    uchar& added) {

  // Get to_subdivide array
  gpu.CreateToSubdivide(num_vertices);
  GetToSubdivideGpu(num_vertices, o.max_level, gpu);

  // Count the number of vertices to subdivide
  cl_event count_event = gpu.EnqueueCountToSubdivide(num_vertices);
  int scount;
  gpu.EnqueueReadCount(&scount, GpuEvents(count_event));
  gpu.Finish();
  added = (scount > 0);

  // Ensure the vertex array is large enough to accommodate subdivided
  // vertices
  const int threshold = num_vertices + scount * kSubAdded;
  const int new_n =
      num_vertices + scount * kSubAdded * o.verts_alloc_factor;
  gpu.EnsureVertices(threshold, new_n);
}

// prepare sparse array
void Stage1b(
    const int num_vertices,
    Gpu& gpu,
    const OctreeOptions& o,
    uchar& added) {
  int sparse_size;
  cl_event csvs_e = gpu.EnqueueComputeSparseVi2GeometriesSize();
  gpu.EnqueueReadSize(&sparse_size, GpuEvents(csvs_e));
  gpu.Finish();
  gpu.MakeSparseVi2Geometries(sparse_size);
  cout << "sparse_size_gpu = " << sparse_size << endl;
}

// Subdivide
void Stage1ci(
    const int num_vertices,
    Gpu& gpu,
    const OctreeOptions& o) {
  GpuEvents e = gpu.EnqueueSubdivideCell(num_vertices);
  // gpu.EnqueueClipGeometries(num_vertices);
  gpu.Finish();
}

// Clip geometries
void Stage1cii(
    const int num_vertices,
    Gpu& gpu,
    const OctreeOptions& o) {
  gpu.EnqueueClipGeometries(num_vertices);
  gpu.Finish();
}

// condense the sparse geometry array
void Stage1d(
    const int num_vertices,
    Gpu& gpu,
    const OctreeOptions& o) {
  gpu.CondenseVi2Geometries();
}

// Initialize the heap
void Stage2(
    const int num_vertices,
    Gpu& gpu,
    const OctreeOptions& o) {
  // Initialize the heap
  gpu.EnsureCPoints();
  gpu.EnqueueComputeNonEmptyVertexDistances(num_vertices);
  gpu.Finish();
}

// Wavefront expansion
void Stage3(
    const int num_vertices,
    Gpu& gpu,
    const OctreeOptions& o) {
  uchar changed = true;
  while (changed) {
    changed = false;
    cl_event write_event = gpu.EnqueueWriteChanged(&changed);
    cl_event gpu_event =
        gpu.EnqueuePullDistances(num_vertices, GpuEvents(write_event));
    gpu.EnqueueReadChanged(&changed, GpuEvents(gpu_event));
    gpu.Finish();
  }
}

// Subdivide ambiguous cells
int Stage4(
    int& num_vertices,
    Gpu& gpu,
    const OctreeOptions& o) {
  for (bool changed = true; changed;) {
    // Get ambiguous cells
    gpu.CreateToSubdivide(num_vertices);
    int tmp = 0;
    cl_event wse = gpu.EnqueueWriteSize(&tmp);
    gpu.Finish();
    cl_event fae =
        gpu.EnqueueFindAmbiguous(num_vertices, o.ambiguous_max_level,
                                 GpuEvents(wse));
    gpu.Finish();
    int scount;
    gpu.EnqueueReadSize(&scount, GpuEvents(fae));
    gpu.Finish();
    cl_event abcde = gpu.EnqueueCountToSubdivide(num_vertices);
    gpu.EnqueueReadCount(&scount, GpuEvents(abcde));
    gpu.Finish();
    changed = (scount > 0);

    if (scount > 0) {
      // Allocate space for new vertices
      const int threshold = num_vertices + scount * kSubAdded;
      const int new_n =
          num_vertices + scount * kSubAdded * o.verts_alloc_factor;
      gpu.EnsureVertices(threshold, new_n);

      // Subdivide
      GpuEvents se = gpu.EnqueueSubdivideCell(num_vertices);
      GpuEvents se2 = gpu.EnqueueClipGeometries(num_vertices, se);
      int header[2];
      gpu.EnqueueReadVNHeader(header, se2);
      gpu.Finish();
      num_vertices = NumVerticesFromHeader(header);

      // Expand wavefront to new cells
      uchar wchanged = true;
      while (wchanged) {
        wchanged = false;
        cl_event write_event = gpu.EnqueueWriteChanged(&wchanged);
        cl_event gpu_event =
            gpu.EnqueuePullDistances(num_vertices, GpuEvents(write_event));
        gpu.EnqueueReadChanged(&wchanged, GpuEvents(gpu_event));
        gpu.Finish();
      }
    }
    if (o.BoolValue("AMBIGUOUS_OUTPUT", false)) {
      cout << "  Ambiguous (gpu2) num: " << scount << endl;
    }
  }
  return num_vertices;
}

// Ensure the vertex array is large enough to accommodate subdivided vertices
void Stage1a_cpu(
    const int num_vertices,
    Gpu& gpu,
    const OctreeOptions& o,
    uchar& added) {

  MVertexNetwork mvertices;
  shared_array<int> vi2geometries_array;
  ReadFromGpu(mvertices, vi2geometries_array, gpu);
  Vi2Geometries vi2geometries = make_vi2geometries(vi2geometries_array.get());
  assert(num_vertices == NumVertices(mvertices));
  assert(0 == mvertices.cpoint_array_capacity);

  // shared_array<Vertex> vertices(new Vertex[num_vertices]);
  // int header[2] = { num_vertices, 0 };

  // prepare to_subdivide array
  int scount = 0;
  // shared_array<uchar> to_subdivide(new uchar[NumVertices(mvertices)]);
  shared_array<uchar> to_subdivide(new uchar[num_vertices]);
  // memset(to_subdivide.get(), NumVertices(mvertices) * sizeof(uchar), 0);
  memset(to_subdivide.get(), num_vertices * sizeof(uchar), 0);

  g_to_subdivide = to_subdivide.get();
  // for (int vi = 0; vi < NumVertices(mvertices); ++vi) {
  for (int vi = 0; vi < num_vertices; ++vi) {
    FindToSubdivide1(
        vi,
        mvertices.header.get(), mvertices.vertices.get(),
        // header, vertices.get(),
        vi2geometries.array,
        o.max_level,
        g_to_subdivide);
  }
  uchar changed2 = true;
  int changed_count = 0;
  while (changed2) {
    changed_count++;
    changed2 = false;
    // for (int vi = 0; vi < NumVertices(mvertices); ++vi) {
    for (int vi = 0; vi < num_vertices; ++vi) {
      FindToSubdivide2(
          vi,
          mvertices.header.get(), mvertices.vertices.get(),
          // header, vertices.get(),
          0,
          vi2geometries.array,
          o.max_level,
          g_to_subdivide,
          &changed2);
    }
  }
  cout << "Gradation required " << changed_count << " iterations" << endl;

  // scount = CountToSubdivide(
  //     NumVertices(mvertices), g_to_subdivide);
  scount = CountToSubdivide(
      num_vertices, g_to_subdivide);
  added = (scount > 0);

  // const int threshold = NumVertices(mvertices) + scount * kSubAdded;
  const int threshold = num_vertices + scount * kSubAdded;
  // const int new_n =
  //     NumVertices(mvertices) + scount * kSubAdded * o.verts_alloc_factor;
  const int new_n =
      num_vertices + scount * kSubAdded * o.verts_alloc_factor;
  mvn_ensure(threshold, new_n, &mvertices);
  // gpu.EnsureVertices(threshold, new_n);

  gpu.CreateToSubdivide(num_vertices);
  gpu.EnqueueWriteToSubdivide(num_vertices, to_subdivide.get());
  gpu.Finish();
  WriteToGpu(mvertices, vi2geometries_array, gpu);
}

// prepare sparse array
void Stage1b_cpu(
    const int num_vertices,
    Gpu& gpu,
    const OctreeOptions& o,
    uchar& added) {
  // Read
  MVertexNetwork mvertices;
  shared_array<int> vi2geometries_array;
  ReadFromGpu(mvertices, vi2geometries_array, gpu);
  Vi2Geometries vi2geometries = make_vi2geometries(vi2geometries_array.get());
  assert(num_vertices == NumVertices(mvertices));
  shared_array<uchar> to_subdivide(new uchar[num_vertices]);
  gpu.EnqueueReadToSubdivide(num_vertices, to_subdivide.get());
  gpu.Finish();

  // Execute
  const int sparse_size = compute_sparse_vi2geometries_size(
      vi2geometries, to_subdivide.get());
  int* sparse_array = new int[sparse_size];
  vi2geometries = make_sparse_vi2geometries(
      vi2geometries, to_subdivide.get(), sparse_array, sparse_size);
  vi2geometries_array.reset(sparse_array);

  // Write
  WriteToGpu(mvertices, vi2geometries_array, gpu);
}

// Subdivide
void Stage1ci_cpu(
    const int num_vertices,
    GeomVertices& geom_vertices,
    Gpu& gpu,
    const OctreeOptions& o) {
  // Read
  MVertexNetwork mvertices;
  shared_array<int> vi2geometries_array;
  ReadFromGpu(mvertices, vi2geometries_array, gpu);
  Vi2Geometries vi2geometries = make_vi2geometries(vi2geometries_array.get());
  assert(num_vertices == NumVertices(mvertices));
  shared_array<uchar> to_subdivide(new uchar[num_vertices]);
  gpu.EnqueueReadToSubdivide(num_vertices, to_subdivide.get());
  gpu.Finish();

  // Execute
  UVertexNetwork vertices = make_vertex_network(mvertices);
  // vertices.num_vertices will change during subdivisions,
  // so cache a copy for iteration.
  const int n = NumVertices(vertices);
  for (int vi = 0; vi < n; ++vi) {
    for (int filter = 0; filter < (1<<DIM); ++filter) {
      // SubdivideCell(
      //     vi, filter, vertices.header, vertices.vertices, to_subdivide.get());
      SubdivideCell_A(
          vi, filter, vertices.header, vertices.vertices, to_subdivide.get());
      // SubdivideCell_B(
      //     vi, filter, vertices.header, vertices.vertices, to_subdivide.get());
      // SubdivideCell_C(
      //     vi, filter, vertices.header, vertices.vertices, to_subdivide.get());
    }
  }

  // Write
  WriteToGpu(mvertices, vi2geometries_array, gpu);
}

// Clip geometries
void Stage1cii_cpu(
    const int num_vertices,
    GeomVertices& geom_vertices,
    Gpu& gpu,
    const OctreeOptions& o) {
  // Read
  MVertexNetwork mvertices;
  shared_array<int> vi2geometries_array;
  ReadFromGpu(mvertices, vi2geometries_array, gpu);
  Vi2Geometries vi2geometries = make_vi2geometries(vi2geometries_array.get());
  // assert(num_vertices == NumVertices(mvertices));
  shared_array<uchar> to_subdivide(new uchar[num_vertices]);
  gpu.EnqueueReadToSubdivide(num_vertices, to_subdivide.get());
  gpu.Finish();

  // Execute
  UVertexNetwork vertices = make_vertex_network(mvertices);
  // vertices.num_vertices will change during subdivisions,
  // so cache a copy for iteration.
  const int n = num_vertices;//NumVertices(vertices);
  for (int vi = 0; vi < n; ++vi) {
    for (int filter = 0; filter < (1<<DIM); ++filter) {
      ClipGeometriesAfterSubdivide(
          vi, filter, vertices.header, vertices.vertices,
          to_subdivide.get(),
          vi2geometries.array,
          geom_vertices.num_vertices, geom_vertices.vertices,
          geom_vertices.num_offsets, geom_vertices.offsets);
    }
  }

  // Write
  WriteToGpu(mvertices, vi2geometries_array, gpu);
}

// condense the sparse geometry array
void Stage1d_cpu(
    const int num_vertices,
    Gpu& gpu,
    const OctreeOptions& o) {
  // Read
  MVertexNetwork mvertices;
  shared_array<int> vi2geometries_array;
  ReadFromGpu(mvertices, vi2geometries_array, gpu);
  Vi2Geometries vi2geometries = make_vi2geometries(vi2geometries_array.get());
  assert(num_vertices == NumVertices(mvertices));
  gpu.Finish();

  // Execute
  UVertexNetwork vertices = make_vertex_network(mvertices);
  const int cur_size = g_array_size(vi2geometries);
  int* dense_array = new int[cur_size];
  vi2geometries = condense_vi2geometries(
      vi2geometries, NumVertices(vertices), dense_array, cur_size);
  vi2geometries_array.reset(dense_array);

  update_mvertex_network(vertices, mvertices);

  // Write
  WriteToGpu(mvertices, vi2geometries_array, gpu);
}

// Initialize the heap
void Stage2_cpu(
    const int num_vertices,
    GeomVertices& geom_vertices,
    Gpu& gpu,
    const OctreeOptions& o) {
  // Read
  MVertexNetwork mvertices;
  shared_array<int> vi2geometries_array;
  ReadFromGpu(mvertices, vi2geometries_array, gpu);
  Vi2Geometries vi2geometries = make_vi2geometries(vi2geometries_array.get());
  assert(num_vertices == NumVertices(mvertices));
  gpu.Finish();

  // Execute
  UVertexNetwork vertices = make_vertex_network(mvertices);
  mvn_cp_reserve(NumVertices(mvertices), &mvertices);
  UVertexNetwork vn = make_vertex_network(mvertices);
  // Compute distances on corners of non-empty cells
  for (int lvi = 0; lvi < (1<<DIM); ++lvi) {
    for (int vi = 0; vi < NumVertices(vn); ++vi) {
      ComputeNonEmptyVertexDistances(
          vi, lvi, vn.header, vn.vertices, vn.cpoints,
          vi2geometries.array,
          geom_vertices.num_vertices, geom_vertices.vertices,
          geom_vertices.num_offsets, geom_vertices.offsets);
    }
  }

  vn = make_vertex_network(mvertices);

  // Write
  WriteToGpu(mvertices, vi2geometries_array, gpu);
}

// Wavefront expansion
void Stage3_cpu(
    const int num_vertices,
    Gpu& gpu,
    const OctreeOptions& o) {
  // Read
  MVertexNetwork mvertices;
  shared_array<int> vi2geometries_array;
  ReadFromGpu(mvertices, vi2geometries_array, gpu);
  Vi2Geometries vi2geometries = make_vi2geometries(vi2geometries_array.get());
  assert(num_vertices == NumVertices(mvertices));
  gpu.Finish();

  // Execute
  UVertexNetwork vn = make_vertex_network(mvertices);
  uchar changed = true;
  while (changed) {
    changed = false;
    for (int vi = 0; vi < NumVertices(vn); ++vi) {
      PullDistances(vi, vn.header, vn.vertices, vn.cpoints, &changed);
    }
  }

  // Write
  WriteToGpu(mvertices, vi2geometries_array, gpu);
}

// Subdivide ambiguous cells
int Stage4_cpu(
    int& num_vertices,
    Gpu& gpu,
    const OctreeOptions& o) {
  // Read
  MVertexNetwork mvertices;
  shared_array<int> vi2geometries_array;
  ReadFromGpu(mvertices, vi2geometries_array, gpu);
  Vi2Geometries vi2geometries = make_vi2geometries(vi2geometries_array.get());
  assert(num_vertices == NumVertices(mvertices));
  gpu.Finish();

  // Execute
  UVertexNetwork vn = make_vertex_network(mvertices);
  assert(num_vertices == NumVertices(vn));
  for (bool changed = true; changed;) {
    // Get ambiguous cells
    // shared_array<uchar> ambiguous(new uchar[NumVertices(vn)]);
    shared_array<uchar> ambiguous(new uchar[num_vertices]);
    int scount = 0;
    // for (int vi = 0; vi < NumVertices(vn); ++vi) {
    for (int vi = 0; vi < num_vertices; ++vi) {
      FindAmbiguous(
          vi, vn.header, vn.vertices, vn.cpoints,
          o.ambiguous_max_level, ambiguous.get(), &scount);
    }
    changed = (scount > 0);
    // Allocate memory
    // const int threshold = NumVertices(mvertices) + scount * kSubAdded;
    const int threshold = num_vertices + scount * kSubAdded;
    // const int new_n =
    //     NumVertices(mvertices) + scount * kSubAdded * o.verts_alloc_factor;
    const int new_n =
        num_vertices + scount * kSubAdded * o.verts_alloc_factor;
    mvn_ensure(threshold, new_n, &mvertices);
    vn = make_vertex_network(mvertices);
    // Do subdivide
    // const int n = NumVertices(mvertices);
    const int n = num_vertices;
    for (int filter = 0; filter < (1<<DIM); ++filter) {
      for (int vi = 0; vi < n; ++vi) {
        SubdivideCell(
            vi,
            filter,
            vn.header, vn.vertices,
            ambiguous.get());//, 0, 0, 0, 0, 0);
      }
    }
    num_vertices = NumVertices(vn);
    uchar wchanged = true;
    while (wchanged) {
      wchanged = false;
      // for (int vi = 0; vi < NumVertices(vn); ++vi) {
      for (int vi = 0; vi < num_vertices; ++vi) {
        PullDistances(vi, vn.header, vn.vertices, vn.cpoints, &wchanged);
      }
    }
    if (o.BoolValue("AMBIGUOUS_OUTPUT", false)) {
      cout << "  Ambiguous (gpu2) num: " << scount << endl;
    }
  }

  // Write
  WriteToGpu(mvertices, vi2geometries_array, gpu);
  return NumVertices(mvertices);
}

//------------------------------------------------------------------------------
// BuildOctreeGpu1
// Builds an octree considering only the geometries and does not perform
// a wavefront expansion.  End result is an octree such that there is a buffer
// of empty cells between cells containing differently-labeled geometries.
//------------------------------------------------------------------------------
void BuildOctreeGpu1(
    Vi2Geometries& vi2geometries,
    shared_array<int>& vi2geometries_array,
    GeomVertices& geom_vertices,
    MVertexNetwork& mvertices,
    Gpu& gpu,
    const OctreeOptions& o) {

  Timer t("BuildOctreeGpu1");
  WallTimer wt("BuildOctreeGpu1 (wall)");

  Stage0(vi2geometries, vi2geometries_array, geom_vertices, mvertices, gpu, o);

  int num_vertices = NumVertices(mvertices);

  // subdivide level by level
  uchar added = true;
  int iter = 0;
  while (added) {
    added = false;
    cout << "Octree iteration " << iter++ << endl;

    // Ensure vertex array is large enough to accommodate subdivided vertices
    if (o.BoolValue("STAGE_1A_GPU", true)) {
      Stage1a(num_vertices, gpu, o, added);
    } else {
      Stage1a_cpu(num_vertices, gpu, o, added);
    }

    // prepare sparse array
    if (o.BoolValue("STAGE_1B_GPU", true)) {
      Stage1b(num_vertices, gpu, o, added);
    } else {
      Stage1b_cpu(num_vertices, gpu, o, added);
    }

    // Subdivide
    if (o.BoolValue("STAGE_1Ci_GPU", true)) {
      Stage1ci(num_vertices, gpu, o);
    } else {
      Stage1ci_cpu(num_vertices, geom_vertices, gpu, o);
    }

    // DO NOT update num_vertices until after geometries are clipped

    // Clip geometries
    if (o.BoolValue("STAGE_1Cii_GPU", true)) {
      Stage1cii(num_vertices, gpu, o);
    } else {
      Stage1cii_cpu(num_vertices, geom_vertices, gpu, o);
    }

    // Read num_vertices from gpu
    int header[2];
    gpu.EnqueueReadVNHeader(header);
    gpu.Finish();
    num_vertices = NumVerticesFromHeader(header);

    // condense the sparse geometry array
    if (o.BoolValue("STAGE_1D_GPU", true)) {
      Stage1d(num_vertices, gpu, o);
    } else {
      Stage1d_cpu(num_vertices, gpu, o);
    }
  }

  // Initialize the heap
  if (o.BoolValue("STAGE_2_GPU", true)) {
    Stage2(num_vertices, gpu, o);
  } else {
    Stage2_cpu(num_vertices, geom_vertices, gpu, o);
  }

  gpu.EnqueueReadVNHeader(mvertices.header.get());
  gpu.Finish();
  mvn_reserve(num_vertices, &mvertices);

  t.restart("Wavefront expansion");

  // Wavefront expansion
  if (o.BoolValue("STAGE_3_GPU", true)) {
    Stage3(num_vertices, gpu, o);
  } else {
    Stage3_cpu(num_vertices, gpu, o);
  }

  // Subdivide ambiguous cells
  t.restart("Ambiguous cell subdivision");
  if (o.BoolValue("STAGE_4_GPU", true)) {
    num_vertices = Stage4(num_vertices, gpu, o);
  } else {
    num_vertices = Stage4_cpu(num_vertices, gpu, o);
  }

  // Read vertices and cpoints from gpu
  gpu.EnqueueReadVNHeader(mvertices.header.get());
  gpu.Finish();
  // Reserve space
  num_vertices = NumVerticesFromHeader(mvertices.header.get());
  mvn_reserve(num_vertices, &mvertices);
  const int num_cpoints = NumCPoints(mvertices);
  mvn_cp_reserve(num_cpoints, &mvertices);
  // Read from gpu
  gpu.EnqueueReadVertices(num_vertices, mvertices.vertices.get());
  gpu.EnqueueReadCPoints(num_cpoints, mvertices.cpoints.get());
  gpu.Finish();
}

//------------------------------------------------------------------------------
// BuildOctreeGpu1_cpu
// Builds an octree considering only the geometries and does not perform
// a wavefront expansion.  End result is an octree such that there is a buffer
// of empty cells between cells containing differently-labeled geometries.
//------------------------------------------------------------------------------
void BuildOctreeGpu1_cpu(
    Vi2Geometries& vi2geometries,
    shared_array<int>& vi2geometries_array,
    GeomVertices& geom_vertices,
    MVertexNetwork& mvertices,
    const OctreeOptions& o) {

  Timer t("BuildOctreeGpu1");
  WallTimer wt("BuildOctreeGpu1 (wall)");

  // cout << "Initial: " << vi2geometries << endl;

  // subdivide level by level
  uchar added = true;
  int iter = 0;
  while (added) {
    added = false;
    cout << "Octree iteration " << iter++ << endl;

    // prepare to_subdivide array
    int scount = 0;
    shared_array<uchar> to_subdivide(new uchar[NumVertices(mvertices)]);
    memset(to_subdivide.get(), NumVertices(mvertices) * sizeof(uchar), 0);

    g_to_subdivide = to_subdivide.get();
    // for (int vi = 0; vi < mvertices.num_vertices; ++vi) {
    //   InitToSubdivide(vi, g_to_subdivide);
    // }
    for (int vi = 0; vi < NumVertices(mvertices); ++vi) {
      FindToSubdivide1(
          vi,
          mvertices.header.get(), mvertices.vertices.get(),
          vi2geometries.array,
          o.max_level,
          g_to_subdivide);
    }
    uchar changed2 = true;
    int changed_count = 0;
    cout << "num_vertices = " << NumVertices(mvertices) << endl;
    while (changed2) {
      changed_count++;
      changed2 = false;
      for (int vi = 0; vi < NumVertices(mvertices); ++vi) {
        FindToSubdivide2(
            vi,
            // NumVertices(mvertices), mvertices.vertices.get(),
            mvertices.header.get(), mvertices.vertices.get(),
            // mvertices.num_cpoints, mvertices.cpoints.get(),
            0,
            vi2geometries.array,
            o.max_level,
            g_to_subdivide,
            &changed2);
      }
    }
    cout << "Gradation required " << changed_count << " iterations" << endl;

    scount = CountToSubdivide(
        NumVertices(mvertices), g_to_subdivide);
    added = (scount > 0);

    const int threshold = NumVertices(mvertices) + scount * kSubAdded;
    const int new_n =
        NumVertices(mvertices) + scount * kSubAdded * o.verts_alloc_factor;
    // mvn_ensure(scount, &mvertices);
    mvn_ensure(threshold, new_n, &mvertices);
    UVertexNetwork vertices = make_vertex_network(mvertices);

    // prepare sparse array
    const int sparse_size = compute_sparse_vi2geometries_size(
        vi2geometries, to_subdivide.get());
    int* sparse_array = new int[sparse_size];
    vi2geometries = make_sparse_vi2geometries(
        vi2geometries, to_subdivide.get(), sparse_array, sparse_size);
    vi2geometries_array.reset(sparse_array);

    // cout << "After sparsify: " << vi2geometries << endl;

    // Subdivide
    const int n = NumVertices(vertices);
    for (int vi = 0; vi < n; ++vi) {
      for (int filter = 0; filter < (1<<DIM); ++filter) {
        SubdivideCell(
            vi, filter, vertices.header, vertices.vertices,
            to_subdivide.get());
      }
    }
    for (int vi = 0; vi < n; ++vi) {
      for (int filter = 0; filter < (1<<DIM); ++filter) {
        ClipGeometriesAfterSubdivide(
            vi, filter, vertices.header, vertices.vertices,
            to_subdivide.get(),
            vi2geometries.array,
            geom_vertices.num_vertices, geom_vertices.vertices,
            geom_vertices.num_offsets, geom_vertices.offsets);
      }
    }
    // SubdivideCellsGpu_cpu(
    //     vi2geometries, to_subdivide.get(), vertices, geom_vertices, o);

    // cout << "After subdivision: " << vi2geometries << endl;

    // condense the sparse geometry array
    // const int cur_size = vi2geometries.array_size;
    const int cur_size = g_array_size(vi2geometries);
    int* dense_array = new int[cur_size];
    vi2geometries = condense_vi2geometries(
        vi2geometries, NumVertices(vertices), dense_array, cur_size);
    vi2geometries_array.reset(dense_array);

    // cout << "After condense: " << vi2geometries << endl;

    update_mvertex_network(vertices, mvertices);
    cout << "num_vertices2 = " << NumVertices(mvertices) << endl;
  }

  // Initialize the heap
  mvn_cp_reserve(NumVertices(mvertices), &mvertices);
  UVertexNetwork vn = make_vertex_network(mvertices);
  // Compute distances on corners of non-empty cells
  for (int lvi = 0; lvi < (1<<DIM); ++lvi) {
    for (int vi = 0; vi < NumVertices(vn); ++vi) {
      ComputeNonEmptyVertexDistances(
          vi, lvi, vn.header, vn.vertices, vn.cpoints,
          vi2geometries.array,
          geom_vertices.num_vertices, geom_vertices.vertices,
          geom_vertices.num_offsets, geom_vertices.offsets);
    }
  }

  vn = make_vertex_network(mvertices);
  t.restart("Wavefront expansion");

  // Wavefront expansion
  {
    uchar changed = true;
    while (changed) {
      changed = false;
      for (int vi = 0; vi < NumVertices(vn); ++vi) {
        PullDistances(vi, vn.header, vn.vertices, vn.cpoints, &changed);
      }
    }
  }

  vn = make_vertex_network(mvertices);

  // Subdivide ambiguous cells
  if (o.ambiguous_max_level > 0) {
    t.restart("Ambiguous cell subdivision");
    int num_vertices = NumVertices(vn);
    for (bool changed = true; changed;) {
      // Get ambiguous cells
      // shared_array<uchar> ambiguous(new uchar[NumVertices(vn)]);
      shared_array<uchar> ambiguous(new uchar[num_vertices]);
      int scount = 0;
      // for (int vi = 0; vi < NumVertices(vn); ++vi) {
      for (int vi = 0; vi < num_vertices; ++vi) {
        FindAmbiguous(
            vi, vn.header, vn.vertices, vn.cpoints,
            o.ambiguous_max_level, ambiguous.get(), &scount);
      }
      changed = (scount > 0);
      // Allocate memory
      // const int threshold = NumVertices(mvertices) + scount * kSubAdded;
      const int threshold = num_vertices + scount * kSubAdded;
      // const int new_n =
      //     NumVertices(mvertices) + scount * kSubAdded * o.verts_alloc_factor;
      const int new_n =
          num_vertices + scount * kSubAdded * o.verts_alloc_factor;
      mvn_ensure(threshold, new_n, &mvertices);
      vn = make_vertex_network(mvertices);
      // Do subdivide
      // const int n = NumVertices(mvertices);
      const int n = num_vertices;
      for (int filter = 0; filter < (1<<DIM); ++filter) {
        for (int vi = 0; vi < n; ++vi) {
          SubdivideCell(
              vi,
              filter,
              vn.header, vn.vertices,
              ambiguous.get());//, 0, 0, 0, 0, 0);
        }
      }
      num_vertices = NumVertices(vn);
      uchar wchanged = true;
      while (wchanged) {
        wchanged = false;
        // for (int vi = 0; vi < NumVertices(vn); ++vi) {
        for (int vi = 0; vi < num_vertices; ++vi) {
          PullDistances(vi, vn.header, vn.vertices, vn.cpoints, &wchanged);
        }
      }
      if (o.BoolValue("AMBIGUOUS_OUTPUT", false)) {
        cout << "  Ambiguous (gpu2) num: " << scount << endl;
      }
    }
  }
}

//------------------------------------------------------------------------------
// BuildOctreeGpu
//------------------------------------------------------------------------------
void BuildOctreeGpu(
    const vector<LabeledGeometry>& geometries,
    GeomVertices& geom_vertices,
    MVertexNetwork& mvertices,
    const OctreeOptions& o) {

  Gpu gpu;

  Vi2Geometries vi2geometries;
  shared_array<int> vi2geometries_array =
      Convert(NumVertices(mvertices), geometries, vi2geometries);
  
  Timer t("Building octree", "*BuildOctree*");
  t.set_output(o.timings || o.report_statistics);

  // Build octree according to geometries -- no wavefront expansion yet
  if (!o.BoolValue("CPU_GPU_IMPL", true)) {
    BuildOctreeGpu1(vi2geometries, vi2geometries_array,
                    geom_vertices, mvertices, gpu, o);
  } else {
    BuildOctreeGpu1_cpu(vi2geometries, vi2geometries_array,
                    geom_vertices, mvertices, o);
  }

  t.stop();
  if (o.report_statistics) {
    UVertexNetwork vn = make_vertex_network(mvertices);
    cout << "Number of octree cells: " << NumCells(vn) << endl;
    cout << "Number of octree vertices: " << NumVertices(vn) << endl;
  }

  // update_mvertex_network(vn, mvertices);
}

NAMESPACE_OCT_END
