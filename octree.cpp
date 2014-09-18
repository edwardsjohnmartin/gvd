#include <algorithm>
#include <iomanip>

#include "./triangle_cpp.h"
#include "./octree.h"
#include "./octree_gpu.h"
#include "./ambiguous.h"
#include "./search.h"
#include "./vector2.h"
#include "./vertices_gpu_state.h"
#include "./geometry_cpp.h"
#include "./opencl/uvertex_network.h"
#include "./opencl/cl_octree.h"

namespace oct {

// See comments to Index() in bit.h.
const int kFaceIndices[4][4][12] = {
  { // 0D
    {  0, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00 }, // the point
    { 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00 },
    { 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00 },
    { 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00 } },
  { // 1D
    {  0,  2, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00 }, // endpoints
    {  1, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00 }, // edge
    { 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00 },
    { 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00 } },
  { // 2D
    {  0,  2,  6,  8, 00, 00, 00, 00, 00, 00, 00, 00 }, // corners
    {  1,  3,  5,  7, 00, 00, 00, 00, 00, 00, 00, 00 }, // edges
    {  4, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00 }, // center
    { 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00 } },
  { // 3D
    {  0,  2,  6,  8, 18, 20, 24, 26, 00, 00, 00, 00 }, // corners
    {  1,  3,  5,  7,  9, 11, 15, 17, 19, 21, 23, 25 }, // edges
    {  4, 10, 12, 14, 16, 22, 00, 00, 00, 00, 00, 00 }, // faces
    { 13, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00 } }// center
};

// To access:
//   kNumFaceIndices[D][dim]
// where dim is the subspace dimension
const int kNumFaceIndices[4][4] = {
  {  1,  0,  0,  0 }, // 0D
  {  2,  1,  0,  0 }, // 1D
  {  4,  4,  1,  0 }, // 2D
  {  8, 12,  6,  1 }  // 3D
};

const int kNumSubdividedArray[4] = { 0, 1, 9, 27 };

//------------------------------------------------------------------------------
// Geometry/array conversion routines
//------------------------------------------------------------------------------

// ** Format for N geometries.  Call the array A.
// offset     description
// ------     -----------
// 0          index offset to geometry 0.  Equal to N.  This value can double
//            as the number of geometries.
// 1..N-1     index offsets to geometries 1..N-1
// N          geometry 0 record
// A[1..N-1]  geometry 1..N-1 record
// 
// ** Geometry record
// offset     description
// ------     -----------
// 0          label
// 1          # triangles = T
// 2..4       triangle 0 - (int, int, int) triple
// 5..2+3*T   triangles 1 - T-1

// The size of an array containing geometries belonging to a
// single octree cell.
template <typename LabeledGeometry>
size_t GeometryArraySize(std::vector<LabeledGeometry>& geoms) {
  const int D = DIM;//LabeledGeometry::D;

  int total_primitive_count = 0;
  for (int i = 0; i < geoms.size(); ++i) {
    total_primitive_count += geoms[i].size();
  }
  // 1 int for each geometry offset
  // 1 int for each geometry label
  // 1 int for each geometry primitive count
  // D*size ints for each geometry
  return std::max((size_t)1, 3*geoms.size() + D*total_primitive_count);
}

template <typename LabeledGeometry>
size_t GeometryArraySize(std::vector<std::vector<LabeledGeometry> >& geoms) {
  int geom_size_sum = 0;
  for (int i = 0; i < geoms.size(); ++i) {
    geom_size_sum += GeometryArraySize(geoms[i]);
  }
  return geom_size_sum + geoms.size();
}

//----------------------------------------
// Geometry -> array
//----------------------------------------

// 0       label
// 1       # primitives (n)
// 2..D*n      primitives
//
// total size: 2 + D*n
// Returns the number of integers written to the array
template <typename LabeledGeometry>
int ConvertGeometryToArray(const LabeledGeometry& geom, int* geom_array){
  const int D = DIM;
  geom_array[0] = geom.GetLabel();
  geom_array[1] = geom.size();
  std::copy(geom.GetPrimitives().begin(), geom.GetPrimitives().end(),
            (Face*)(geom_array+2));
  return 2 + geom.size() * D;
}

// Convert the geometries belonging to an octree cell to an array.
// Returns the number of integers written to the array.
template <typename LabeledGeometry>
int ConvertGeometriesToArray(vector<LabeledGeometry>& geoms, int* geom_array){
  // const int D = LabeledGeometry::D;
  // const int D = DIM;
  int offset = geoms.size();
  for (int i = 0; i < geoms.size(); ++i) {
    // offset to the geometry record
    geom_array[i] = offset;
    offset += ConvertGeometryToArray(geoms[i], geom_array+offset);
  }
  if (offset == 0) {
    // no geometries
    // todo - change back?
    // geom_array[0] = -1;
    geom_array[0] = 0;
    offset = 1;
  }
  return offset;
}

template <typename LabeledGeometry>
shared_array<int> ConvertGeometriesToArray(vector<LabeledGeometry>& geoms){
  const size_t n = GeometryArraySize(geoms);
  shared_array<int> geom_array(new int[n]);
  ConvertGeometriesToArray(geoms, geom_array.get());
  return geom_array;
}

// 2-dimensional array of geometries
template <typename LabeledGeometry>
shared_array<int> ConvertGeometriesToArray(
    std::vector<std::vector<LabeledGeometry> >& geoms){
  const int total_size = GeometryArraySize(geoms);
  shared_array<int> geom_array(new int[total_size]);
  int offset = geoms.size();
  for (int i = 0; i < geoms.size(); ++i) {
    geom_array[i] = offset;
    offset += ConvertGeometriesToArray(geoms[i], geom_array.get()+offset);
  }
  return geom_array;
}

//----------------------------------------
// Array -> geometry
//----------------------------------------

template <typename LabeledGeometry>
vector<LabeledGeometry> ConvertArrayToGeometries(const int* geom_array) {
  const int D = DIM;

  std::vector<LabeledGeometry> geoms;
  int num_geoms = geom_array[0];
  for (int i = 0; i < num_geoms; i++) {
    int offset = geom_array[i];
    int label = geom_array[offset];
    int num_tris = geom_array[offset+1];
    std::vector<Face> tris;
    std::copy((Face*)(geom_array+offset+2),
              (Face*)(geom_array+offset+2+num_tris*D),
              std::back_inserter(tris));
    LabeledGeometry lgeom(tris, label);
    geoms.push_back(lgeom);
  }
  return geoms;
}

// 2-dimensional array of geometries
template <typename LabeledGeometry>
std::vector<std::vector<LabeledGeometry> > ConvertArrayToGeometries2(
    const int* geom_array) {
  std::vector<std::vector<LabeledGeometry> > geoms;
  const int num_arrays = geom_array[0];
  for (int i = 0; i < num_arrays; i++) {
    int offset = geom_array[i];
    geoms.push_back(ConvertArrayToGeometries<LabeledGeometry>(
        geom_array+offset));
  }
  return geoms;
}

// //------------------------------------------------------------
// // InitWaveGpu
// // Initialize the wavefront using the gpu
// //------------------------------------------------------------
// void InitWaveGpu_old(
//     vector<int>& nonempty_vertices,
//     vector<intn >& nonempty_points,
//     shared_array<level_t> all_levels,
//     vector<vector<LabeledGeometry> >& base2geometries,
//     GeomVertices& geom_vertices,
//     UVertexNetwork& vertices,
//     multiset<HeapVertex>& heap,
//     const OctreeOptions& o) {
// #ifdef __OPEN_CL_SUPPORT__
//   static const int kNumCorners = (1 << DIM);

//   intn* all_geom_vertices = geom_vertices.vertices;
//   int* all_geom_vertex_offsets = geom_vertices.offsets;
//   const int num_geom_vertices = geom_vertices.num_vertices;
//   const int num_geom_offsets = geom_vertices.num_offsets;

//   // Number of cells to compute distances for
//   const int N = nonempty_vertices.size();
//   shared_array<intn> base_points(new intn[N]);
//   shared_array<int> all_corner_vertices(new int[N*kNumCorners]);
//   for (int i = 0; i < N; ++i) {
//     const int vi = nonempty_vertices[i];
//     const intn p = nonempty_points[i];
//     base_points[i] = p;
//     memcpy(all_corner_vertices.get()+i*kNumCorners,
//            vertices.vertices[vi].corners, (1<<DIM) * sizeof(int));
//   }
//   shared_array<intn> all_points_gpu(new intn[N*kNumCorners]);
//   shared_array<int> all_labels_gpu(new int[N*kNumCorners]);

//   shared_array<int> geom_array =
//       ConvertGeometriesToArray(base2geometries);
//   const int geom_array_size = GeometryArraySize(base2geometries);
//   const level_t max_level = o.max_level;

//   // Setup OpenCL
//   // input/output buffers
//   int error;
//   cl_mem gpu_base_points =
//       clCreateBuffer(context_cl, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
//                      sizeof(intn)*N, base_points.get(), &error);
//   cl_mem gpu_all_cell_levels =
//       clCreateBuffer(context_cl, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
//                      sizeof(level_t)*N, all_levels.get(), &error);
//   cl_mem gpu_all_geom_array =
//       clCreateBuffer(context_cl, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
//                      sizeof(int)*geom_array_size, geom_array.get(), &error);
//   cl_mem gpu_all_corner_vertices =
//       clCreateBuffer(context_cl, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
//                      sizeof(int)*N*kNumCorners, all_corner_vertices.get(), &error);
//   cl_mem gpu_points =
//       clCreateBuffer(context_cl, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR,
//                      sizeof(intn)*N*kNumCorners, all_points_gpu.get(), &error);
//   cl_mem gpu_labels =
//       clCreateBuffer(context_cl, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR,
//                      sizeof(int)*N*kNumCorners, all_labels_gpu.get(), &error);
//   cl_mem gpu_all_geom_vertices =
//       clCreateBuffer(context_cl, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
//                      sizeof(intn)*num_geom_vertices, all_geom_vertices,
//                      &error);
//   cl_mem gpu_all_geom_vertex_offsets =
//       clCreateBuffer(context_cl, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
//                      sizeof(int)*num_geom_offsets,
//                      all_geom_vertex_offsets,
//                      &error);
//   CheckError (error);

//   // Setup the kernel_init_wave arguments
//   clSetKernelArg(
//       kernel_init_wave, 0, sizeof(cl_mem), &gpu_base_points);
//   clSetKernelArg(
//       kernel_init_wave, 1, sizeof(cl_mem), &gpu_all_cell_levels);
//   clSetKernelArg(
//       kernel_init_wave, 2, sizeof(cl_mem), &gpu_all_geom_array);
//   clSetKernelArg(
//       kernel_init_wave, 3, sizeof(cl_mem), &gpu_all_corner_vertices);
//   clSetKernelArg(
//       kernel_init_wave, 4, sizeof(cl_mem), &gpu_points);
//   clSetKernelArg(
//       kernel_init_wave, 5, sizeof(cl_mem), &gpu_labels);
//   clSetKernelArg(
//       kernel_init_wave, 6, sizeof(level_t), &max_level);
//   clSetKernelArg(
//       kernel_init_wave, 7, sizeof(cl_mem), &gpu_all_geom_vertices);
//   clSetKernelArg(
//       kernel_init_wave, 8, sizeof(cl_mem), &gpu_all_geom_vertex_offsets);

//   cl_command_queue queue_cl = CreateCommandQueue(context_cl, &error);
//   CheckError (error);

//   error = clEnqueueWriteBuffer(
//       queue_cl, gpu_base_points, CL_TRUE, 0,
//       sizeof(intn)*N, base_points.get(), 0, nullptr, nullptr);
//   CheckError(error, "clEnqueueWriteBuffer");
//   error = clEnqueueWriteBuffer(
//       queue_cl, gpu_all_cell_levels, CL_TRUE, 0,
//       sizeof(level_t)*N, all_levels.get(), 0, nullptr, nullptr);
//   CheckError(error, "clEnqueueWriteBuffer");
//   error = clEnqueueWriteBuffer(
//       queue_cl, gpu_all_geom_array, CL_TRUE, 0,
//       sizeof(int)*geom_array_size, geom_array.get(), 0, nullptr, nullptr);
//   CheckError(error, "clEnqueueWriteBuffer");
//   error = clEnqueueWriteBuffer(
//       queue_cl, gpu_all_corner_vertices, CL_TRUE, 0,
//       sizeof(int)*N*kNumCorners, all_corner_vertices.get(), 0, nullptr, nullptr);
//   CheckError(error, "clEnqueueWriteBuffer");
//   error = clEnqueueWriteBuffer(
//       queue_cl, gpu_all_geom_vertices, CL_TRUE, 0,
//       sizeof(intn)*num_geom_vertices, all_geom_vertices, 0, nullptr,
//       nullptr);
//   CheckError(error, "clEnqueueWriteBuffer");
//   error = clEnqueueWriteBuffer(
//       queue_cl, gpu_all_geom_vertex_offsets, CL_TRUE, 0,
//       sizeof(int)*num_geom_offsets, all_geom_vertex_offsets, 0, nullptr,
//       nullptr);
//   CheckError(error, "clEnqueueWriteBuffer");

//   // Run the processing
//   size_t globalWorkSize[] = { N };
//   cl_event gpu_event;
//   error = clEnqueueNDRangeKernel(
//       queue_cl, kernel_init_wave, 1, nullptr, globalWorkSize,
//       nullptr, 0, nullptr, &gpu_event);
//   CheckError(error, "clEnqueueNDRangeKernel");
    
//   // Read results
//   error = clEnqueueReadBuffer(
//       queue_cl, gpu_points, CL_TRUE, 0,
//       sizeof(intn)*N*kNumCorners, all_points_gpu.get(), 1, &gpu_event, nullptr);
//   CheckError(error, "clEnqueueReadBuffer");
//   error = clEnqueueReadBuffer(
//       queue_cl, gpu_labels, CL_TRUE, 0,
//       sizeof(int)*N*kNumCorners, all_labels_gpu.get(), 1, &gpu_event, nullptr);
//   CheckError(error, "clEnqueueReadBuffer");

//   clFinish(queue_cl);

//   clReleaseMemObject(gpu_base_points);
//   clReleaseMemObject(gpu_all_cell_levels);
//   clReleaseMemObject(gpu_all_geom_array);
//   clReleaseMemObject(gpu_all_corner_vertices);
//   clReleaseMemObject(gpu_points);
//   clReleaseMemObject(gpu_labels);
//   clReleaseMemObject(gpu_all_geom_vertices);
//   clReleaseMemObject(gpu_all_geom_vertex_offsets);

//   clReleaseCommandQueue(queue_cl);
//   // End OpenCL

//   if (o.BoolValue("TEST_INIT_WAVE_GPU", false)) {
//     cout << "Testing gpu wave initialization" << endl;
//     shared_array<intn> all_points_cpu(new intn[N*kNumCorners]);
//     shared_array<int> all_labels_cpu(new int[N*kNumCorners]);
//     for (int i = 0; i < N; ++i) {
//       ComputeCornerDistancesParallel3<DIM>(
//           i,
//           base_points.get(), all_levels.get(), geom_array.get(),
//           all_corner_vertices.get(), all_points_cpu.get(),
//           all_labels_cpu.get(), o.max_level,
//           // all_geom_vertices.get(), all_geom_vertex_offsets.get());
//           all_geom_vertices, all_geom_vertex_offsets);
//     }

//     int done_count = 0;
//     for (int i = 0; done_count < 3*kNumCorners && i < N; ++i) {
//       for (int j = 0; j < kNumCorners; ++j) {
//         const int k = i*kNumCorners+j;
//         bool bad = false;
//         for (int l = 0; l < DIM; ++l) {
//           // Don't worry about round off
//           if (abs(all_points_cpu[k].s[l]-all_points_gpu[k].s[l]) > 2) {
//             bad = true;
//           }
//         }
//         if (bad) {
//           cout << i << "," << j << " "
//                << all_points_cpu[k] << "\t" << all_points_gpu[k] << endl;
//           done_count++;
//         }
//       }
//     }
//     done_count = 0;
//     for (int i = 0; done_count < 3*kNumCorners && i < N; ++i) {
//       for (int j = 0; j < kNumCorners; ++j) {
//         const int k = i*kNumCorners+j;
//         if (all_labels_cpu[k] != all_labels_gpu[k]) {
//           cout << i << "/" << N << "," << j << " "
//                << all_labels_cpu[k] << "\t" << all_labels_gpu[k] << endl;
//           done_count++;
//         }
//       }
//     }
//   }

//   intn* all_points = all_points_gpu.get();
//   int* all_labels = all_labels_gpu.get();

//   for (int i = 0; i < N; ++i) {
//     intn* points = all_points+i*kNumCorners;
//     int* labels = all_labels+i*kNumCorners;
//     const intn& base_point = base_points[i];
//     const level_t level = all_levels[i];
//     const index_t width = Level2CellWidth(level);
//     int* corner_vertices = all_corner_vertices.get()+i*kNumCorners;
//     for (int lvi = 0; lvi < kNumCorners; ++lvi) {
//       const int vi = corner_vertices[lvi];
//       const intn v_point = Position(lvi, base_point, width);
//       const intn& point = points[lvi];
//       const int dist = length(point - v_point);
//       const int label = labels[lvi];
//       const int curve_pointi = CreateNewClosestPoint(point, label, &vertices);
//       const bool closest_changed =
//           SetClosestPoint(vi, curve_pointi, &vertices);
//       if (closest_changed) {
//         heap.insert(HeapVertex(dist, vi));
//       }
//     }
//   }
// #endif //__OPEN_CL_SUPPORT__
// }

template <int D>
vector<vector<int> > ComputeBase2Incident(const VertexNetwork& vertices) {
  // Given a vertex, find the base vertices of the up to 2^D cells that
  // it is incident to.
  FindBasesVisitor<D> fbv(&vertices);
  VisitVerticesBFS<D>(vertices, fbv);

  // For each base vertex, store all vertices incident to its cell.
  return fbv.GetBase2Incident();
}

//------------------------------------------------------------------------------
// GetNeighborsToSubdivide
//
// Returns a vector of all neighbors that need to be subdivided.
//------------------------------------------------------------------------------
template <int D>
vector<int> GetNeighborsToSubdivide(
    const int vi,
    const intn& base_point,
    const char level,
    VertexNetwork& vertices) {
  typedef Direction Dir;
  vector<int> nbrs;

  for (int axis = 0; axis < D; ++axis) {
    for (int pos = 0; pos < 2; ++pos) {
      const Dir d = DirectionFromAxis(axis, pos);
      const int n_vi = vertices.FindNeighbor(vi, d, level);
      if (n_vi > -1 &&
          vertices.IsBase(n_vi) &&
          vertices.CellLevel(n_vi) == level) {
        nbrs.push_back(n_vi);
      }
    }
  }
  return nbrs;
}

// Subdivide a level of the octree.
// base_vis are the base vertices to subdivide.
template <typename LabeledGeometry, typename AddedIter>
void SubdivideCells(
    const vector<int>& candidates,
    vector<vector<LabeledGeometry> >& base2geometries,
    VertexNetwork& vertices,
    const GeomVertices& geom_vertices,
    vector<set<int> >& adjacent_cells,
    const OctreeOptions& o,
    AddedIter added) {

  for (int i = 0; i < candidates.size(); ++i) {
    const int vi = candidates[i];
    const vector<LabeledGeometry>& geometries = base2geometries[vi];

    const bool subdivide = 
        !geometries.empty() &&
        (HasMultipleLabels(geometries) || o.full_subdivide);
    
    if (subdivide) {
      SubdivideCell(
          vi, base2geometries, geom_vertices, 0, adjacent_cells,
          added, vertices);
    }
  }
}

//------------------------------------------------------------------------------
// Subdivide cells until at least one empty buffer cell exists between
// objects.
//------------------------------------------------------------------------------
// template <int D, typename LabeledGeometry>
void MakeBuffer(
    VertexNetwork& vertices,
    const GeomVertices& geom_vertices,
    std::vector<std::vector<LabeledGeometry> >& base2geometries,
    vector<set<int> >& adjacent_cells,
    const OctreeOptions& o) {
  typedef Direction Dir;
  // static const int kNumSubdivided = kNumSubdividedArray[DIM];

  vector<int> candidates;
  for (int vi = 0; vi < vertices.size(); ++vi) {
    candidates.push_back(vi);
  }

  bool subdivided = true;
  while (subdivided) {
    subdivided = false;

    vector<int> added;
    // Perform subdivisions
    const int N = vertices.size();
    for (int vi = 0; vi < N; ++vi) {
      if (vertices.IsBase(vi) &&
          vertices.CellLevel(vi) < o.max_level &&
          !adjacent_cells[vi].empty()) {
        subdivided = true;
        list<int> queue;
        int slvi2vi[kNumSubdivided];
        SubdivideCell(
            vi, base2geometries, geom_vertices,
            slvi2vi, adjacent_cells, back_inserter(queue), vertices);
      }
    }
  }
}

void MakeBufferGpu(
    MVertexNetwork& mvertices,
    const GeomVertices& geom_vertices,
    std::vector<std::vector<LabeledGeometry> >& base2geometries,
    vector<set<int> >& adjacent_cells,
    const OctreeOptions& o) {
  typedef Direction Dir;
  // static const int kNumSubdivided = kNumSubdividedArray[DIM];

  bool subdivided = true;
  while (subdivided) {
    subdivided = false;

    // Count subdivisions
    const int N = NumVertices(mvertices);
    int scount = 0;
    for (int vi = 0; vi < N; ++vi) {
      if (v_is_base(mvertices.vertices[vi]) &&
          v_cell_level(mvertices.vertices[vi]) < o.max_level &&
          !adjacent_cells[vi].empty()) {
        ++scount;
      }
    }
    // mvn_ensure(scount, &mvertices);
    const int threshold = NumVertices(mvertices) + scount * kSubAdded;
    const int new_n =
        NumVertices(mvertices) + scount * kSubAdded * o.verts_alloc_factor;
    mvn_ensure(threshold, new_n, &mvertices);
    UVertexNetwork vertices = make_vertex_network(mvertices);
  
    vector<int> added;
    // Perform subdivisions
    for (int vi = 0; vi < N; ++vi) {
      if (IsBase(vi, vertices) &&
          CellLevel(vi, vertices) < o.max_level &&
          !adjacent_cells[vi].empty()) {
        subdivided = true;
        list<int> queue;
        int slvi2vi[kNumSubdivided];
        SubdivideCellGpu(
            vi, base2geometries, geom_vertices,
            slvi2vi, adjacent_cells, back_inserter(queue), vertices);
      }
    }

    update_mvertex_network(vertices, mvertices);
  }

}

//------------------------------------------------------------------------------
// SetDistances
//
// Performs a wavefront expansion out from the cells in the heap, setting
// distances on vertices as it goes.
//------------------------------------------------------------------------------
void SetDistances(
    VertexNetwork& vertices, std::multiset<HeapVertex>& heap,
                  const OctreeOptions& o) {
  using namespace std;

  const bool output = o.BoolValue("SET_DISTANCES_OUTPUT", false);
  Timer t("SetDistances");
  if (!output) {
    t.kill();
  }

  int num_visited = 0;
  int heap_size = heap.size();
  if (output)
    cout << "SetDistances";

  // Set the distances on the vertices belonging to leaf cells that
  // contain a portion of the surface.
  // The heap represents the wavefront expansion.
  std::vector<bool> processed(vertices.size(), false);
  while (!heap.empty()) {
    HeapVertex v = *heap.begin();
    heap.erase(heap.begin());
    --heap_size;

    const int vi = v.vi;
    // const intn v_point = v.v_point;
    const intn v_point = vertices.Position(vi);
    if (!processed[vi]) {
      processed[vi] = true;
      // output
      ++num_visited;
      if (output) {
        if (num_visited % 1000 == 0)
          cout << "." << flush;
        if (num_visited % 20000 == 0)
          cout << (num_visited/1000) << "K" << flush;
      }
      // end output

      // Get all directions that have a valid neighbor
      // vector<Direction<D> > nbr_dirs;
      vector<Direction> nbr_dirs;
      for (int i = 0; i < DIM; ++i) {
        for (int pos = 0; pos < 2; ++pos) {
          // const Direction<D> d = Direction<D>::FromAxis(i, pos);
          const Direction d = DirectionFromAxis(i, pos);
          const int n_vi = vertices.Neighbor(vi, d);
          if (n_vi != -1) {
            nbr_dirs.push_back(d);
          }
        }
      }
      const int label = vertices.Label(vi);
      // Iterate through the neighbors
      for (int i = 0; i < nbr_dirs.size(); ++i) {
        const HeapVertex hv =
            PushDistance<DIM>(vertices, vi, nbr_dirs[i], v_point, label, o);
        if (hv.dist > -1) {
          heap.insert(hv);
          ++heap_size;
          processed[hv.vi] = false;
        }
      }
    }
  }
  if (output) {
    cout << endl;
    cout << "SetDistances num_visited = " << num_visited << endl;
  }
}

//------------------------------------------------------------------------------
// SetDistances
//
// Performs a wavefront expansion out from the cells in the heap, setting
// distances on vertices as it goes.
//------------------------------------------------------------------------------
void SetDistances(
    UVertexNetwork& vertices, std::multiset<HeapVertex>& heap,
                  const OctreeOptions& o) {
  using namespace std;

  const bool output = o.BoolValue("SET_DISTANCES_OUTPUT", false);
  Timer t("SetDistances");
  if (!output) {
    t.kill();
  }

  int num_visited = 0;
  int heap_size = heap.size();
  if (output)
    cout << "SetDistances";

  // Set the distances on the vertices belonging to leaf cells that
  // contain a portion of the surface.
  // The heap represents the wavefront expansion.
  std::vector<bool> processed(NumVertices(vertices), false);//size(), false);
  while (!heap.empty()) {
    HeapVertex v = *heap.begin();
    heap.erase(heap.begin());
    --heap_size;

    const int vi = v.vi;
    // const intn v_point = v.v_point;
    // const intn v_point = vertices.Position(vi);
    const intn v_point = vertices.vertices[vi].position;
    if (!processed[vi]) {
      processed[vi] = true;
      // output
      ++num_visited;
      if (output) {
        if (num_visited % 1000 == 0)
          cout << "." << flush;
        if (num_visited % 20000 == 0)
          cout << (num_visited/1000) << "K" << flush;
      }
      // end output

      // Get all directions that have a valid neighbor
      // vector<Direction<D> > nbr_dirs;
      vector<Direction> nbr_dirs;
      for (int i = 0; i < DIM; ++i) {
        for (int pos = 0; pos < 2; ++pos) {
          // const Direction<D> d = Direction<D>::FromAxis(i, pos);
          const Direction d = DirectionFromAxis(i, pos);
          // const int n_vi = vertices.Neighbor(vi, d);
          const int n_vi = Neighbor(vi, d, vertices);
          if (n_vi != -1) {
            nbr_dirs.push_back(d);
          }
        }
      }
      // const int label = vertices.Label(vi);
      const int label = Label(vi, vertices);
      // Iterate through the neighbors
      for (int i = 0; i < nbr_dirs.size(); ++i) {
        const HeapVertex hv =
            PushDistance(vertices, vi, nbr_dirs[i], v_point, label, o);
        if (hv.dist > -1) {
          heap.insert(hv);
          ++heap_size;
          processed[hv.vi] = false;
        }
      }
    }
  }
  if (output) {
    cout << endl;
    cout << "SetDistances num_visited = " << num_visited << endl;
  }
}

template <int D>
void DisplayLevelHistogram(const VertexNetwork& vertices,
                           const OctreeOptions& o) {
  const int N = o.max_level;
  vector<int> count(N+1, 0);
  for (int vi = 0; vi < vertices.size(); ++vi) {
    if (vertices.IsBase(vi))
      count[vertices.CellLevel(vi)]++;
  }
  double max_count = 0;
  for (int i = 0; i < N+1; ++i) {
    max_count = max(max_count, (double)count[i]);
  }
  cout << "Level histogram:" << endl;
  for (int i = 0; i < N+1; ++i) {
    cout << "  " << i << ": " << setw(7) << count[i] << " ";
    const int c = 40 * (count[i] / max_count);
    for (int j = 0; j < c; ++j) {
      cout << "*";
    }
    cout << endl;
  }
  
}

//------------------------------------------------------------------------------
// BuildOctreeCpu
//------------------------------------------------------------------------------
// template <typename LabeledGeometry>
void BuildOctreeCpu(
    // const vector<vector<floatn> >& all_vertices,
    const vector<LabeledGeometry>& geometries,
    GeomVertices& geom_vertices,
    // const vector<vector<Face> >& all_faces,
    const BoundingBox<floatn>& bb,
    VertexNetwork& vertices,
    const OctreeOptions& o) {
  static const int D = DIM;

  // const int* corner_vertices = vertices.GetCorners(0);

  std::multiset<HeapVertex> heap;
  // Accessed as base2geometries[base_vi].  This gives all geometries
  // that intersect with the cell.
  std::vector<std::vector<LabeledGeometry> > base2geometries;
  base2geometries.push_back(geometries);
  std::vector<LabeledGeometry> leaf_geometries;

  Timer t("Building octree", "*BuildOctree*");
  t.set_output(o.timings || o.report_statistics);

  // Build the octree
  // if (o.BoolValue("TEST", false)) {
    
  // } else {
    // Original implementation - no gpu
    // Subdivide level by level
    level_t level = 0;
    list<int> added;
    added.push_back(0);
    vector<set<int> > adjacent_cells(1);
    while (!added.empty() && level < o.max_level) {
      const vector<int> candidates(added.begin(), added.end());
      added.clear();
      SubdivideCells(
          candidates, base2geometries, vertices, geom_vertices,
          adjacent_cells, o, back_inserter(added));
      ++level;
    }
    // cout << "num std levels = " << (int) level << endl;

    // Subdivide until we have a buffer of empty cells between objects
    if (!o.full_subdivide && o.make_buffer) {
      t.restart("Making buffer");
      // cout << "  Number of octree cells: " << vertices.NumCells() << endl;
      // DisplayLevelHistogram(vertices, o);
      MakeBuffer(
          vertices, geom_vertices, base2geometries, adjacent_cells, o);
    }
  // }

  t.restart("Initialize wavefront - compute corner distances");
  // cout << "  Number of octree cells: " << vertices.NumCells() << endl;
  // DisplayLevelHistogram(vertices, o);

  // Initialize wavefront with non-empty vertices
  std::vector<int> nonempty_vertices;
  std::vector<intn> nonempty_points;
  NonEmptyVertexCollector<D, LabeledGeometry> nonempty(
      vertices, nonempty_vertices, nonempty_points, base2geometries);
  VisitVertices<D>(vertices, nonempty);

  const int N = nonempty_vertices.size();
  shared_array<level_t> all_levels(new level_t[N]);
  for (int i = 0; i < N; ++i) {
    const int vi = nonempty_vertices[i];
    all_levels[i] = vertices.CellLevel(vi);
  } 
  // if (o.gpu) {
  //   InitWaveGpu<D>(
  //       nonempty_vertices, nonempty_points, all_levels, base2geometries,
  //       geom_vertices, vertices, heap, o);
  // } else {
    for (int i = 0; i < N; ++i) {
      const int vi = nonempty_vertices[i];
      const intn p = nonempty_points[i];
      const level_t level = all_levels[i];
      const int* corner_vertices = vertices.GetCorners(vi);
      ComputeCornerDistances<D>(
          p, level, base2geometries[vi], vertices,
          corner_vertices, geom_vertices, heap, o);
    }
  // }

  t.restart("Wavefront expansion");

  // Wavefront expansion
  SetDistances(vertices, heap, o);

  // Subdivide ambiguous cells
  if (o.ambiguous_max_level > 0) {
    t.restart("Ambiguous cell subdivision");
    // SubdivideAmbiguous(vertices, gpu_state, o);
    // if (o.gpu) {
    //   VerticesGpuState<D> gpu_state(kernel_update_vertices, vertices);
    //   SubdivideAmbiguousGpu(vertices, gpu_state, o);
    // } else {
      SubdivideAmbiguousCpu(vertices, o);
    // }
  }

  t.stop();
  if (o.report_statistics) {
    cout << "Number of octree cells: " << vertices.NumCells() << endl;
    // DisplayLevelHistogram(vertices, o);
    cout << "Number of octree vertices: " << vertices.size() << endl;
  }
}

std::vector<std::vector<LabeledGeometry> > test_base2geometries;

void GpuTest(
    const vector<LabeledGeometry>& geometries,
    GeomVertices& geom_vertices,
    const BoundingBox<floatn>& bb,
    VertexNetwork& vertices,
    const OctreeOptions& o) {
  std::vector<std::vector<LabeledGeometry> > base2geometries;
  base2geometries.push_back(geometries);

  level_t level = 0;
  list<int> added;
  added.push_back(0);
  vector<set<int> > adjacent_cells(1);
  while (!added.empty() && level < o.max_level) {
    const vector<int> candidates(added.begin(), added.end());
    added.clear();
    SubdivideCells(
        candidates, base2geometries, vertices, geom_vertices,
        adjacent_cells, o, back_inserter(added));
    ++level;
  }
  test_base2geometries = base2geometries;
}

//------------------------------------------------------------------------------
// BuildOctree
//
// Given polylines or triangles, build an octree.
//------------------------------------------------------------------------------
// ManagedVertexNetwork BuildOctree(
VertexNetwork BuildOctree(
// void BuildOctree(
    const vector<vector<floatn> >& all_vertices,
    const vector<vector<Face> >& all_faces,
    const BoundingBox<floatn>& bb,
    // VertexNetwork& vertices,
    // ManagedVertexNetwork& mvertices,
    const OctreeOptions& o) {
  static const int D = DIM;
  // const int kNumCorners = (1 << D);

  // ManagedVertexNetwork mvertices = make_mvertex_network();
  // VertexNetwork vertices = make_vertex_network(mvertices);
  // VertexNetwork vertices = make_vertex_network();

  // if (vertices.size() > (1<<D)) {
  //   cerr << vertices.size() << endl;
  //   throw logic_error("vertices passed to BuildOctree must be empty");
  // }

  // if (all_faces.size() != all_vertices.size()) return make_vertex_network();
  // if (all_faces.empty()) return make_vertex_network();
  // if (all_faces[0].empty()) return make_vertex_network();
  if (all_faces.size() != all_vertices.size()) return VertexNetwork();
  if (all_faces.empty()) return VertexNetwork();
  if (all_faces[0].empty()) return VertexNetwork();

  Timer t("Convert vertices", "*BuildOctree*");
  t.set_output(o.timings || o.report_statistics);

  // Find bounding box of vertices
  BoundingBox<floatn> vert_bb;
  for (int j = 0; j < all_vertices.size(); ++j) {
    const std::vector<floatn>& fvertices = all_vertices[j];
    for (int i = 0; i < fvertices.size(); ++i) {
      vert_bb(fvertices[i]);
    }
  }

  if (!bb.IsSquare()) {
    // cerr << "Bounding box must be square: " << bb.size()
    //      << " (but it's probably close enough)" << endl;
    // throw logic_error("Bounding box must be square");
  }

  const float dwidth = bb.size().s[0];

  // if (dwidth == 0) return make_vertex_network();
  if (dwidth == 0) return VertexNetwork();

  // Count total number of geometry vertices
  int num_geom_vertices = 0;
  for (int j = 0; j < all_vertices.size(); ++j) {
    num_geom_vertices += all_vertices[j].size();
  }

  // Convert vertices to integer coordinates
  std::vector<LabeledGeometry> geometries;
  shared_array<intn> all_geom_vertices(new intn[num_geom_vertices]);
  shared_array<int> all_geom_vertex_offsets(new int[all_vertices.size()]);
  GeomVertices geom_vertices = {
      num_geom_vertices, all_geom_vertices.get(),
      all_vertices.size(), all_geom_vertex_offsets.get() };
  int geom_vert_idx = 0;
  for (int j = 0; j < all_vertices.size(); ++j) {
    all_geom_vertex_offsets[j] = geom_vert_idx;
    const std::vector<floatn>& fvertices = all_vertices[j];
    const int n = fvertices.size();
    for (int i = 0; i < n; ++i) {
      intn p = make_intn(0);
      for (int k = 0; k < D; ++k) {
        const double d =
            (kWidth-1) * ((fvertices[i].s[k] - bb.min().s[k]) / dwidth);
        int v = static_cast<int>(d+0.5);
        if (v < 0) {
          cerr << "Coordinate in dimension " << k << " is less than zero.  d = "
              << d << " v = " << v << endl;
          cerr << "  fvertices[i][k] = " << fvertices[i].s[k]
               << " bb.min()[k] = " << bb.min().s[k] << endl;
          cerr << "  dwidth = " << dwidth << " kwidth = " << kWidth << endl;
          v = 0;
        }
        p.s[k] = v;
      }
      all_geom_vertices[geom_vert_idx++] = p;
    }
    if (n > 0) {
#ifdef OCT3D
      LabeledGeometry lg(all_faces[j], j);
#else
      // LabeledGeometry lg(geom_vertices.Vertices(j), all_faces[j], j);
      // todo not sure if .vertices should be dependent on j
      LabeledGeometry lg(geom_vertices.vertices, all_faces[j], j);
#endif
      geometries.push_back(lg);
    }
  }

  VertexNetwork vertices;
  if (vertices.size() > (1<<D)) {
    cerr << vertices.size() << endl;
    throw logic_error("vertices passed to BuildOctree must be empty");
  }

  if (o.gpu) {
    GpuTest(geometries, geom_vertices, bb, vertices, o);

    MVertexNetwork mvertices = make_mvertex_network();
    BuildOctreeGpu(geometries, geom_vertices, mvertices, o);
    vertices = VertexNetwork(
        NumVertices(mvertices), mvertices.vertices.get(),
        NumCPoints(mvertices), mvertices.cpoints.get());
  } else {
    BuildOctreeCpu(geometries, geom_vertices, bb, vertices, o);
  }

  // return make_mvertex_network(vertices);
  return vertices;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Template instantiations
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifdef OCT2D
template size_t GeometryArraySize(std::vector<LabeledGeometry2>& geoms);
template size_t GeometryArraySize(
    std::vector<std::vector<LabeledGeometry2> >& geoms);
template std::vector<LabeledGeometry2> ConvertArrayToGeometries(
    const int* geom_array);
template int ConvertGeometriesToArray(
    std::vector<LabeledGeometry2>& geoms, int* geom_array);
template shared_array<int> ConvertGeometriesToArray(
    std::vector<LabeledGeometry2>& geoms);
template std::vector<std::vector<LabeledGeometry2> > ConvertArrayToGeometries2(
    const int* geom_array);
template shared_array<int> ConvertGeometriesToArray(
    std::vector<std::vector<LabeledGeometry2> >& geoms);
template vector<vector<int> > ComputeBase2Incident<2>(
    const VertexNetwork& vertices);
template vector<int> GetNeighborsToSubdivide<2>(
    const int vi,
    const int2& base_point,
    const char level,
    VertexNetwork& vertices);
#endif

#ifdef OCT3D
template size_t GeometryArraySize(std::vector<LabeledGeometry3>& geoms);
template size_t GeometryArraySize(
    std::vector<std::vector<LabeledGeometry3> >& geoms);
template std::vector<LabeledGeometry3> ConvertArrayToGeometries(
    const int* geom_array);
template int ConvertGeometriesToArray(
    std::vector<LabeledGeometry3>& geoms, int* geom_array);
template shared_array<int> ConvertGeometriesToArray(
    std::vector<LabeledGeometry3>& geoms);
template std::vector<std::vector<LabeledGeometry3> > ConvertArrayToGeometries2(
    const int* geom_array);
template shared_array<int> ConvertGeometriesToArray(
    std::vector<std::vector<LabeledGeometry3> >& geoms);
template vector<vector<int> > ComputeBase2Incident<3>(
    const VertexNetwork& vertices);
template vector<int> GetNeighborsToSubdivide<3>(
    const int vi,
    const int3& base_point,
    const char level,
    VertexNetwork& vertices);
#endif

}
