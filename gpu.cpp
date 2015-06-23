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

#ifdef __OPEN_CL_SUPPORT__

#include "./gpu.h"

using namespace std;

NAMESPACE_OCT_BEGIN

Gpu::Gpu()
    : mpoints(0),
      to_subdivide(0), vertices(0), vi2nbrs(0), cpoints(0), vi2geometries(0),
      gvertices(0), goffsets(0),
      to_subdivide_n(0), vertices_n(0), cpoints_n(0), vi2geometries_n(0),
      gvertices_n(0), goffsets_n(0) {
  int error;
  queue = CreateCommandQueue(context_cl, &error);
  CheckError(error);

  // Create vn_header buffer
  vn_header = clCreateBuffer(
      context_cl, CL_MEM_READ_WRITE, sizeof(int) * 2, nullptr, &error);
  CheckError(error, "CreateVNHeader");

  // Create changed buffer
  changed = clCreateBuffer(
      context_cl, CL_MEM_WRITE_ONLY, sizeof(uchar), nullptr, &error);
  CheckError(error, "CreateChanged");

  // Create count buffer
  count = clCreateBuffer(
      context_cl, CL_MEM_WRITE_ONLY, sizeof(int), nullptr, &error);
  CheckError(error, "CreateCount");

  // Create size buffer
  size = clCreateBuffer(
      context_cl, CL_MEM_WRITE_ONLY, sizeof(int), nullptr, &error);
  CheckError(error, "CreateSize");
}

Gpu::~Gpu() {
  if (to_subdivide)
    clReleaseMemObject(to_subdivide);
  if (vn_header)
    clReleaseMemObject(vn_header);
  if (vertices)
    clReleaseMemObject(vertices);
  if (vi2nbrs)
    clReleaseMemObject(vi2nbrs);
  if (cpoints)
    clReleaseMemObject(cpoints);
  if (vi2geometries)
    clReleaseMemObject(vi2geometries);
  if (changed)
    clReleaseMemObject(changed);
  if (count)
    clReleaseMemObject(count);
  if (size)
    clReleaseMemObject(size);
  if (gvertices)
    clReleaseMemObject(gvertices);
  if (goffsets)
    clReleaseMemObject(goffsets);

  clReleaseCommandQueue(queue);
}

//------------------------------------------------------------
// Create buffer functions
//------------------------------------------------------------

void Gpu::CreateToSubdivide(int n) {
  if (to_subdivide)
    clReleaseMemObject(to_subdivide);

  // to_subdivide_size = sizeof(uchar) * n;
  to_subdivide_n = n;

  int error;
  to_subdivide = clCreateBuffer(
      // context_cl, CL_MEM_READ_WRITE, to_subdivide_size, nullptr, &error);
      context_cl, CL_MEM_READ_WRITE, sizeof(uchar) * to_subdivide_n,
      nullptr, &error);
  CheckError(error, "CreateToSubdivide");
}

void Gpu::CreateMPoints(const int n) {
  if (mpoints)
    clReleaseMemObject(mpoints);

  mpoints_n = n;

  int error;
  mpoints = clCreateBuffer(
      context_cl, CL_MEM_READ_ONLY, sizeof(int)*mpoints_n,
      nullptr, &error);
  CheckError(error, "CreateMPoints");
}

void Gpu::CreateVertices(int n) {
  if (vertices)
    clReleaseMemObject(vertices);
  if (vi2nbrs)
    clReleaseMemObject(vi2nbrs);

  vertices_n = n;

  int error;
  vertices = clCreateBuffer(
      context_cl, CL_MEM_READ_WRITE, sizeof(Vertex)*vertices_n,
      nullptr, &error);
  CheckError(error, "CreateVertices");

  vi2nbrs = clCreateBuffer(
      context_cl, CL_MEM_READ_WRITE, sizeof(int)*vertices_n*kNumIncidentCells,
      nullptr, &error);
  CheckError(error, "CreateVi2Nbrs");
}

void Gpu::CreateCPoints(int n) {
  if (cpoints)
    clReleaseMemObject(cpoints);

  // cpoints_size = sizeof(GeomPoint) * n;
  cpoints_n = n;

  int error;
  cpoints = clCreateBuffer(
      // context_cl, CL_MEM_READ_WRITE, cpoints_size, nullptr, &error);
      context_cl, CL_MEM_READ_WRITE, sizeof(GeomPoint) * cpoints_n,
      nullptr, &error);
  CheckError(error, "CreateCPoints");
}

void Gpu::CreateVi2Geometries(int n) {
  if (vi2geometries)
    clReleaseMemObject(vi2geometries);

  // vi2geometries_size = sizeof(int) * n;
  vi2geometries_n = n;

  int error;
  vi2geometries = clCreateBuffer(
      // context_cl, CL_MEM_READ_WRITE, vi2geometries_size, nullptr, &error);
      context_cl, CL_MEM_READ_WRITE, sizeof(int) * vi2geometries_n,
      nullptr, &error);
  CheckError(error, "CreateVi2Geometries");
}

void Gpu::EnsureVertices(int threshold, int new_n, GpuEvents wait_events) {
  if (!vertices) {
    CreateVertices(new_n);
    return;
  }

  const int old_n = vertices_n;
  // if (old_n < new_n) {
  if (old_n < threshold) {
    const int new_vertices_n = new_n;

    int error;
    cl_mem old_vertices = vertices;
    cl_mem new_vertices = clCreateBuffer(
        context_cl, CL_MEM_READ_WRITE, sizeof(Vertex) * new_vertices_n,
        nullptr, &error);
    CheckError(error, "CreateVertices");

    cl_mem old_vi2nbrs = vi2nbrs;
    cl_mem new_vi2nbrs = clCreateBuffer(
        context_cl, CL_MEM_READ_WRITE,
        sizeof(int) * new_vertices_n * kNumIncidentCells,
        nullptr, &error);
    CheckError(error, "CreateVi2Nbrs");

    // Call kernel to copy data
    cl_kernel k = kernel_copy_vertex_array;
    clSetKernelArg(k, 0, sizeof(cl_mem), &new_vertices);
    clSetKernelArg(k, 1, sizeof(cl_mem), &old_vertices);
    clSetKernelArg(k, 2, sizeof(int), &old_n);

    size_t globalWorkSize = old_n;
    cl_event gpu_event;
    error = clEnqueueNDRangeKernel(
        queue, k, 1, nullptr, &globalWorkSize,
        nullptr, wait_events.size(), wait_events.get(), &gpu_event);
    CheckError(error, "clEnqueueNDRangeKernel");
    Finish();

    clReleaseMemObject(old_vertices);
    vertices = new_vertices;
    vertices_n = new_vertices_n;

    clReleaseMemObject(old_vi2nbrs);
    vi2nbrs = new_vi2nbrs;
  }
}

// Reserves one cpoint space for each vertex.  If there already exists an
// array, the elements are NOT copied.
void Gpu::EnsureCPoints(GpuEvents wait_events) {
  const int new_n = vertices_n;
  if (!cpoints) {
    CreateCPoints(new_n);
    return;
  }

  const int old_n = cpoints_n;
  if (old_n < new_n) {
    const int new_cpoints_n = new_n;

    int error;
    cl_mem old_cpoints = cpoints;
    cl_mem new_cpoints = clCreateBuffer(
        context_cl, CL_MEM_READ_WRITE, sizeof(GeomPoint) * new_cpoints_n,
        nullptr, &error);
    CheckError(error, "CreateCPoints");

    // // Call kernel to copy data
    // cl_kernel k = kernel_copy_vertex_array;
    // clSetKernelArg(k, 0, sizeof(cl_mem), &new_vertices);
    // clSetKernelArg(k, 1, sizeof(cl_mem), &old_vertices);
    // clSetKernelArg(k, 2, sizeof(int), &old_n);

    // size_t globalWorkSize = old_n;
    // cl_event gpu_event;
    // error = clEnqueueNDRangeKernel(
    //     queue, k, 1, nullptr, &globalWorkSize,
    //     nullptr, wait_events.size(), wait_events.get(), &gpu_event);
    // CheckError(error, "clEnqueueNDRangeKernel");
    // Finish();

    clReleaseMemObject(old_cpoints);
    cpoints = new_cpoints;
    cpoints_n = new_cpoints_n;
  }
}

void Gpu::CreateGeomVertices(GeomVertices geom_vertices) {
  if (gvertices)
    clReleaseMemObject(gvertices);
  if (goffsets)
    clReleaseMemObject(goffsets);

  gvertices_n = geom_vertices.num_vertices;
  goffsets_n = geom_vertices.num_offsets;

  int error;
  gvertices = clCreateBuffer(
      context_cl, CL_MEM_READ_WRITE, sizeof(intn)*gvertices_n,
      nullptr, &error);
  CheckError(error, "CreateGeomVertices");
  goffsets = clCreateBuffer(
      context_cl, CL_MEM_READ_WRITE, sizeof(int)*goffsets_n,
      nullptr, &error);
  CheckError(error, "CreateGeomVertices");
}

//------------------------------------------------------------
// Write functions
//------------------------------------------------------------

cl_event Gpu::EnqueueWriteToSubdivide(
    int n, uchar* array, GpuEvents wait_events) {
  cl_event event;
  int error = clEnqueueWriteBuffer(
      queue, to_subdivide, CL_FALSE, 0,
      sizeof(uchar) * n, array, wait_events.size(), wait_events.get(), &event);
  CheckError(error, "WriteToSubdivide");
  return event;
}

cl_event Gpu::EnqueueWriteVNHeader(
    int* header, GpuEvents wait_events) {
  cl_event event;
  int error = clEnqueueWriteBuffer(
      queue, vn_header, CL_FALSE, 0,
      sizeof(int) * 2, header, wait_events.size(), wait_events.get(), &event);
  CheckError(error, "WriteVNHeader");
  return event;
}

cl_event Gpu::EnqueueWriteVertices(
    int n, Vertex* array, GpuEvents wait_events) {
  cl_event event;
  int error = clEnqueueWriteBuffer(
      queue, vertices, CL_FALSE, 0,
      sizeof(Vertex) * n, array, wait_events.size(), wait_events.get(), &event);
  CheckError(error, "WriteVertices");
  return event;
}

cl_event Gpu::EnqueueWriteCPoints(
    int n, GeomPoint* array, GpuEvents wait_events) {
  cl_event event;
  int error = clEnqueueWriteBuffer(
      queue, cpoints, CL_FALSE, 0,
      sizeof(GeomPoint) * n, array,
      wait_events.size(), wait_events.get(), &event);
  CheckError(error, "WriteCPoints");
  return event;
}

cl_event Gpu::EnqueueWriteVi2Geometries(int* array, GpuEvents wait_events) {
  cl_event event;
  int error = clEnqueueWriteBuffer(
      queue, vi2geometries, CL_FALSE, 0,
      sizeof(int) * vi2geometries_n, array,
      wait_events.size(), wait_events.get(), &event);
  CheckError(error, "WriteVi2Geometries");
  return event;
}

cl_event Gpu::EnqueueWriteChanged(uchar* value, GpuEvents wait_events) {
  cl_event event;
  int error = clEnqueueWriteBuffer(
      queue, changed, CL_FALSE, 0,
      sizeof(uchar), value, wait_events.size(), wait_events.get(), &event);
  CheckError(error, "WriteChanged");
  return event;
}

cl_event Gpu::EnqueueWriteSize(int* value, GpuEvents wait_events) {
  cl_event event;
  int error = clEnqueueWriteBuffer(
      queue, size, CL_FALSE, 0,
      sizeof(int), value, wait_events.size(), wait_events.get(), &event);
  CheckError(error, "WriteSize");
  return event;
}

GpuEvents Gpu::EnqueueWriteGeomVertices(
    GeomVertices geom_vertices, GpuEvents wait_events) {
  cl_event e0, e1;
  int error;

  error = clEnqueueWriteBuffer(
      queue, gvertices, CL_FALSE, 0,
      sizeof(intn)*gvertices_n, geom_vertices.vertices,
      wait_events.size(), wait_events.get(), &e0);
  CheckError(error, "WriteGeomVertices");

  error = clEnqueueWriteBuffer(
      queue, goffsets, CL_FALSE, 0,
      sizeof(int)*goffsets_n, geom_vertices.offsets,
      wait_events.size(), wait_events.get(), &e1);
  CheckError(error, "WriteGeomVertices");

  return GpuEvents(e0, e1);
}


//------------------------------------------------------------
// Read functions
//------------------------------------------------------------

cl_event Gpu::EnqueueReadToSubdivide(
    int n, uchar* array, GpuEvents wait_events) {
  cl_event event;
  int error = clEnqueueReadBuffer(
      queue, to_subdivide, CL_FALSE, 0,
      sizeof(uchar) * n, array, wait_events.size(), wait_events.get(), &event);
  CheckError(error, "ReadToSubdivide");
  return event;
}

cl_event Gpu::EnqueueReadVNHeader(
    int* header, GpuEvents wait_events) {
  cl_event event;
  int error = clEnqueueReadBuffer(
      queue, vn_header, CL_FALSE, 0,
      sizeof(int) * 2, header, wait_events.size(), wait_events.get(), &event);
  CheckError(error, "ReadVNHeader");
  return event;
}

cl_event Gpu::EnqueueReadVertices(
    int n__, Vertex* array, GpuEvents wait_events) {
  assert(n__ <= vertices_n);
  cl_event event;
  int error = clEnqueueReadBuffer(
      queue, vertices, CL_FALSE, 0,
      // sizeof(Vertex) * vertices_n, array,
      sizeof(Vertex) * n__, array,
      wait_events.size(), wait_events.get(), &event);
  CheckError(error, "ReadVertices");
  return event;
}

cl_event Gpu::EnqueueReadCPoints(
    int n, GeomPoint* array, GpuEvents wait_events) {
  assert(n <= cpoints_n);
  cl_event event;
  int error = clEnqueueReadBuffer(
      queue, cpoints, CL_FALSE, 0,
      sizeof(GeomPoint) * n, array,
      wait_events.size(), wait_events.get(), &event);
  CheckError(error, "ReadCPoints");
  return event;
}

cl_event Gpu::EnqueueReadVi2Geometries(int* array, GpuEvents wait_events) {
  cl_event event;
  int error = clEnqueueReadBuffer(
      queue, vi2geometries, CL_FALSE, 0,
      sizeof(int) * vi2geometries_n, array,
      wait_events.size(), wait_events.get(), &event);
  CheckError(error, "ReadVi2Geometries");
  return event;
}

cl_event Gpu::EnqueueReadChanged(
    uchar* value, GpuEvents wait_events) {
  cl_event event;
  int error = clEnqueueReadBuffer(
      queue, changed, CL_FALSE, 0,
      sizeof(uchar), value, wait_events.size(), wait_events.get(), &event);
  CheckError(error, "ReadChanged");
  return event;
}

cl_event Gpu::EnqueueReadCount(int* value, GpuEvents wait_events) {
  cl_event event;
  int error = clEnqueueReadBuffer(
      queue, count, CL_FALSE, 0,
      sizeof(int), value, wait_events.size(), wait_events.get(), &event);
  CheckError(error, "ReadCount");
  return event;
}

cl_event Gpu::EnqueueReadSize(int* value, GpuEvents wait_events) {
  cl_event event;
  int error = clEnqueueReadBuffer(
      queue, size, CL_FALSE, 0,
      sizeof(int), value, wait_events.size(), wait_events.get(), &event);
  CheckError(error, "ReadSize");
  return event;
}

//------------------------------------------------------------
// Kernel functions
//------------------------------------------------------------
// cl_event Gpu::EnqueueFindNbrs(
//     int num_vertices, GpuEvents wait_events) {

//   cl_event vi2nbrs_event;
//   {
//     cl_kernel k = kernel_find_nbrs;
//     int ai = 0;
//     clSetKernelArg(k, ai++, sizeof(cl_mem), &vn_header);
//     clSetKernelArg(k, ai++, sizeof(cl_mem), &vertices);
//     clSetKernelArg(k, ai++, sizeof(cl_mem), &vi2nbrs);

//     int error;
//     size_t globalWorkSize[] = { num_vertices };
//     cl_event vi2nbrs_event;
//     error = clEnqueueNDRangeKernel(
//         queue, k, 1, nullptr, globalWorkSize,
//         nullptr, wait_events.size(), wait_events.get(), &vi2nbrs_event);
//     CheckError(error, "clEnqueueNDRangeKernel");
//   }
// }

cl_event Gpu::EnqueueFindToSubdivide1(
    int num_vertices, int max_level, GpuEvents wait_events) {
  cl_kernel k = kernel_find_to_subdivide1;
  // clSetKernelArg(k, 0, sizeof(int), &num_vertices);
  clSetKernelArg(k, 0, sizeof(cl_mem), &vn_header);
  clSetKernelArg(k, 1, sizeof(cl_mem), &vertices);
  clSetKernelArg(k, 2, sizeof(cl_mem), &vi2geometries);
  clSetKernelArg(k, 3, sizeof(int), &max_level);
  clSetKernelArg(k, 4, sizeof(cl_mem), &to_subdivide);

  int error;
  size_t globalWorkSize[] = { static_cast<size_t>(num_vertices) };
  cl_event gpu_event;
  error = clEnqueueNDRangeKernel(
      queue, k, 1, nullptr, globalWorkSize,
      nullptr, wait_events.size(), wait_events.get(), &gpu_event);
  CheckError(error, "clEnqueueNDRangeKernel");
  return gpu_event;
}

cl_event Gpu::EnqueueFindToSubdivide2(
    int num_vertices, int max_level, GpuEvents wait_events) {
  // Compute vi2nbrs
  cl_event vi2nbrs_event;
  {
    cl_kernel k = kernel_find_nbrs;
    int ai = 0;
    clSetKernelArg(k, ai++, sizeof(cl_mem), &vn_header);
    clSetKernelArg(k, ai++, sizeof(cl_mem), &vertices);
    clSetKernelArg(k, ai++, sizeof(cl_mem), &vi2nbrs);

    int error;
    size_t globalWorkSize[] = { static_cast<size_t>(num_vertices) };
    error = clEnqueueNDRangeKernel(
        queue, k, 1, nullptr, globalWorkSize,
        nullptr, wait_events.size(), wait_events.get(), &vi2nbrs_event);
    CheckError(error, "clEnqueueNDRangeKernel");
  }

  cl_kernel k = kernel_find_to_subdivide2;
  int ai = 0;
  clSetKernelArg(k, ai++, sizeof(cl_mem), &vn_header);
  clSetKernelArg(k, ai++, sizeof(cl_mem), &vertices);
  clSetKernelArg(k, ai++, sizeof(cl_mem), &vi2nbrs);
  clSetKernelArg(k, ai++, sizeof(cl_mem), &vi2geometries);
  clSetKernelArg(k, ai++, sizeof(int), &max_level);
  clSetKernelArg(k, ai++, sizeof(cl_mem), &to_subdivide);
  clSetKernelArg(k, ai++, sizeof(cl_mem), &changed);

  int error;
  size_t globalWorkSize[] = { static_cast<size_t>(num_vertices) };
  cl_event gpu_event;
  GpuEvents wait_events2(wait_events, vi2nbrs_event);
  error = clEnqueueNDRangeKernel(
      queue, k, 1, nullptr, globalWorkSize,
      nullptr, wait_events2.size(), wait_events2.get(), &gpu_event);
  CheckError(error, "clEnqueueNDRangeKernel");
  return gpu_event;
}

cl_event Gpu::EnqueueCountToSubdivide(int num_vertices, GpuEvents wait_events) {
  cl_kernel k = kernel_count_to_subdivide;
  clSetKernelArg(k, 0, sizeof(int), &num_vertices);
  clSetKernelArg(k, 1, sizeof(cl_mem), &to_subdivide);
  clSetKernelArg(k, 2, sizeof(cl_mem), &count);

  int error;
  size_t globalWorkSize = 1;
  cl_event gpu_event;
  error = clEnqueueNDRangeKernel(
      queue, kernel_count_to_subdivide, 1, nullptr, &globalWorkSize,
      nullptr, wait_events.size(), wait_events.get(), &gpu_event);
  CheckError(error, "clEnqueueNDRangeKernel");

  return gpu_event;
}

cl_event Gpu::EnqueueComputeSparseVi2GeometriesSize(
    GpuEvents wait_events) {
  cl_kernel k = kernel_compute_sparse_vi2geometries_size;
  clSetKernelArg(k, 0, sizeof(cl_mem), &vi2geometries);
  clSetKernelArg(k, 1, sizeof(cl_mem), &to_subdivide);
  clSetKernelArg(k, 2, sizeof(cl_mem), &size);

  int error;
  size_t globalWorkSize = 1;
  cl_event gpu_event;
  error = clEnqueueNDRangeKernel(
      queue, k, 1, nullptr, &globalWorkSize,
      nullptr, wait_events.size(), wait_events.get(), &gpu_event);
  CheckError(error, "clEnqueueNDRangeKernel");

  return gpu_event;
}

void Gpu::MakeSparseVi2Geometries(
    const int new_n, GpuEvents wait_events) {

  const int old_n = vi2geometries_n;

  int error;
  cl_mem old_vi2geometries = vi2geometries;
  cl_mem new_vi2geometries = clCreateBuffer(
      context_cl, CL_MEM_READ_WRITE, sizeof(int) * new_n,
      nullptr, &error);
  CheckError(error, "MakeSparseVi2Geometries");

  // Call kernel to copy data
  cl_kernel k = kernel_make_sparse_vi2geometries;
  clSetKernelArg(k, 0, sizeof(cl_mem), &old_vi2geometries);
  clSetKernelArg(k, 1, sizeof(cl_mem), &to_subdivide);
  clSetKernelArg(k, 2, sizeof(cl_mem), &new_vi2geometries);
  clSetKernelArg(k, 3, sizeof(int), &new_n);

  size_t globalWorkSize = 1;
  cl_event gpu_event;
  error = clEnqueueNDRangeKernel(
      queue, k, 1, nullptr, &globalWorkSize,
      nullptr, wait_events.size(), wait_events.get(), &gpu_event);
  CheckError(error, "clEnqueueNDRangeKernel");
  Finish();

  clReleaseMemObject(old_vi2geometries);
  vi2geometries = new_vi2geometries;
  vi2geometries_n = new_n;
}

void Gpu::CondenseVi2Geometries(GpuEvents wait_events) {
  const int old_n = vi2geometries_n;

  int error;
  size_t globalWorkSize;
  cl_event gpu_event;

  // Get new_n
  int new_n;
  cl_kernel kc = kernel_compute_dense_vi2geometries_size;
  clSetKernelArg(kc, 0, sizeof(cl_mem), &vi2geometries);
  clSetKernelArg(kc, 1, sizeof(cl_mem), &size);
  globalWorkSize = 1;
  error = clEnqueueNDRangeKernel(
      queue, kc, 1, nullptr, &globalWorkSize,
      nullptr, wait_events.size(), wait_events.get(), &gpu_event);
  CheckError(error, "clEnqueueNDRangeKernel");
  EnqueueReadSize(&new_n, GpuEvents(gpu_event));
  Finish();

  // Allocate space on gpu for new array
  // int error;
  cl_mem old_vi2geometries = vi2geometries;
  cl_mem new_vi2geometries = clCreateBuffer(
      context_cl, CL_MEM_READ_WRITE, sizeof(int) * new_n,
      nullptr, &error);
  CheckError(error, "CondenseVi2Geometries");

  // Call kernel to copy data
  cl_kernel k = kernel_condense_vi2geometries;
  clSetKernelArg(k, 0, sizeof(cl_mem), &old_vi2geometries);
  clSetKernelArg(k, 1, sizeof(cl_mem), &vn_header);
  clSetKernelArg(k, 2, sizeof(cl_mem), &new_vi2geometries);
  clSetKernelArg(k, 3, sizeof(int), &new_n);

  globalWorkSize = 1;
  error = clEnqueueNDRangeKernel(
      queue, k, 1, nullptr, &globalWorkSize,
      nullptr, wait_events.size(), wait_events.get(), &gpu_event);
  CheckError(error, "clEnqueueNDRangeKernel");
  Finish();

  clReleaseMemObject(old_vi2geometries);
  vi2geometries = new_vi2geometries;
  vi2geometries_n = new_n;
}

// void Gpu::CondenseVi2Geometries(
//     int num_vertices, GpuEvents wait_events) {
//   const int n = vi2geometries_n;

//   int error;
//   cl_mem old_vi2geometries = vi2geometries;
//   cl_mem new_vi2geometries = clCreateBuffer(
//       context_cl, CL_MEM_READ_WRITE, sizeof(int) * n,
//       nullptr, &error);
//   CheckError(error, "CondenseVi2Geometries");

//   // Call kernel to copy data
//   cl_kernel k = kernel_condense_vi2geometries;
//   clSetKernelArg(k, 0, sizeof(cl_mem), &old_vi2geometries);
//   clSetKernelArg(k, 1, sizeof(int), &num_vertices);
//   clSetKernelArg(k, 2, sizeof(cl_mem), &new_vi2geometries);
//   clSetKernelArg(k, 3, sizeof(int), &n);

//   size_t globalWorkSize = 1;
//   cl_event gpu_event;
//   error = clEnqueueNDRangeKernel(
//       queue, k, 1, nullptr, &globalWorkSize,
//       nullptr, wait_events.size(), wait_events.get(), &gpu_event);
//   CheckError(error, "clEnqueueNDRangeKernel");
//   Finish();

//   clReleaseMemObject(old_vi2geometries);
//   vi2geometries = new_vi2geometries;
// }

// GpuEvents Gpu::EnqueueSubdivideCell(
//     int num_vertices, GpuEvents wait_events) {
//   cl_kernel k = kernel_subdivide_cell;
//   clSetKernelArg(k, 1, sizeof(cl_mem), &vn_header);
//   clSetKernelArg(k, 2, sizeof(cl_mem), &vertices);
//   clSetKernelArg(k, 3, sizeof(cl_mem), &to_subdivide);

//   int error;
//   size_t globalWorkSize[] = { num_vertices };
//   GpuEvents w_events = wait_events;
//   GpuEvents events;
//   for (int filter = 0; filter < (1<<DIM); ++filter) {
//     clSetKernelArg(k, 0, sizeof(int), &filter);
//     cl_event gpu_event;
//     error = clEnqueueNDRangeKernel(
//         queue, k, 1, nullptr, globalWorkSize,
//         nullptr, w_events.size(), w_events.get(), &gpu_event);
//     w_events = GpuEvents(w_events, gpu_event);
//     events = GpuEvents(events, gpu_event);
//     CheckError(error, "clEnqueueNDRangeKernel");
//   }
//   return events;
// }

GpuEvents Gpu::EnqueueSubdivideCell(
    int num_vertices, GpuEvents wait_events) {
  cl_kernel k[] = { kernel_subdivide_cell_A,
                    kernel_subdivide_cell_B,
                    kernel_subdivide_cell_C };
  for (int i = 0; i < 3; ++i) {
    clSetKernelArg(k[i], 1, sizeof(cl_mem), &vn_header);
    clSetKernelArg(k[i], 2, sizeof(cl_mem), &vertices);
    clSetKernelArg(k[i], 3, sizeof(cl_mem), &to_subdivide);
    clSetKernelArg(k[i], 4, sizeof(cl_mem), &size);
  }

  int error;
  size_t globalWorkSize[] = { static_cast<size_t>(num_vertices) };
  GpuEvents w_events = wait_events;
  GpuEvents events;
  for (int filter = 0; filter < (1<<DIM); ++filter) {
    for (int i = 0; i < 1; ++i) {
      // int size = 0;
      // Finish();
      // EnqueueWriteSize(&size);
      // Finish();

      clSetKernelArg(k[i], 0, sizeof(int), &filter);
      cl_event gpu_event;
      error = clEnqueueNDRangeKernel(
          queue, k[i], 1, nullptr, globalWorkSize,
          nullptr, w_events.size(), w_events.get(), &gpu_event);
      w_events = GpuEvents(w_events, gpu_event);
      events = GpuEvents(events, gpu_event);
      CheckError(error, "clEnqueueNDRangeKernel");

      // if (filter == 0) {
      //   Finish();
      //   EnqueueReadSize(&size);
      //   Finish();
      //   cout << "i = " << i << " size = " << size << endl;
      // }
    }
    // exit(0);
  }
  return events;
}

GpuEvents Gpu::EnqueueClipGeometries(
    int num_vertices, GpuEvents wait_events) {
  cl_kernel k = kernel_clip_geometries;
  clSetKernelArg(k, 1, sizeof(cl_mem), &vn_header);
  clSetKernelArg(k, 2, sizeof(cl_mem), &vertices);
  clSetKernelArg(k, 3, sizeof(cl_mem), &to_subdivide);
  clSetKernelArg(k, 4, sizeof(cl_mem), &vi2geometries);
  clSetKernelArg(k, 5, sizeof(int), &gvertices_n);
  clSetKernelArg(k, 6, sizeof(cl_mem), &gvertices);
  clSetKernelArg(k, 7, sizeof(int), &goffsets_n);
  clSetKernelArg(k, 8, sizeof(cl_mem), &goffsets);

  int error;
  size_t globalWorkSize[] = { static_cast<size_t>(num_vertices) };
  GpuEvents w_events = wait_events;
  GpuEvents events;
  for (int filter = 0; filter < (1<<DIM); ++filter) {
    clSetKernelArg(k, 0, sizeof(int), &filter);
    cl_event gpu_event;
    error = clEnqueueNDRangeKernel(
        queue, k, 1, nullptr, globalWorkSize,
        nullptr, w_events.size(), w_events.get(), &gpu_event);
    w_events = GpuEvents(w_events, gpu_event);
    events = GpuEvents(events, gpu_event);
    CheckError(error, "clEnqueueNDRangeKernel");
  }
  return events;
}

GpuEvents Gpu::EnqueueComputeNonEmptyVertexDistances(
    int num_vertices, GpuEvents wait_events) {
  cl_kernel k = kernel_compute_non_empty_vertex_distances;
  clSetKernelArg(k, 1, sizeof(cl_mem), &vn_header);
  clSetKernelArg(k, 2, sizeof(cl_mem), &vertices);
  clSetKernelArg(k, 3, sizeof(cl_mem), &cpoints);
  clSetKernelArg(k, 4, sizeof(cl_mem), &vi2geometries);
  clSetKernelArg(k, 5, sizeof(int), &gvertices_n);
  clSetKernelArg(k, 6, sizeof(cl_mem), &gvertices);
  clSetKernelArg(k, 7, sizeof(int), &goffsets_n);
  clSetKernelArg(k, 8, sizeof(cl_mem), &goffsets);

  int error;
  size_t globalWorkSize[] = { static_cast<size_t>(num_vertices) };
  GpuEvents w_events = wait_events;
  GpuEvents events;
  for (int filter = 0; filter < (1<<DIM); ++filter) {
    clSetKernelArg(k, 0, sizeof(int), &filter);
    cl_event gpu_event;
    error = clEnqueueNDRangeKernel(
        queue, k, 1, nullptr, globalWorkSize,
        nullptr, w_events.size(), w_events.get(), &gpu_event);
    w_events = GpuEvents(w_events, gpu_event);
    events = GpuEvents(events, gpu_event);
    CheckError(error, "clEnqueueNDRangeKernel");
  }
  return events;
}

cl_event Gpu::EnqueuePullDistances(
    int num_vertices, GpuEvents wait_events) {
  cl_kernel k = kernel_pull_distances;
  clSetKernelArg(k, 0, sizeof(cl_mem), &vn_header);
  clSetKernelArg(k, 1, sizeof(cl_mem), &vertices);
  clSetKernelArg(k, 2, sizeof(cl_mem), &cpoints);
  clSetKernelArg(k, 3, sizeof(cl_mem), &changed);

  int error;
  size_t globalWorkSize[] = { static_cast<size_t>(num_vertices) };
  cl_event gpu_event;
  error = clEnqueueNDRangeKernel(
      queue, k, 1, nullptr, globalWorkSize,
      nullptr, wait_events.size(), wait_events.get(), &gpu_event);
  CheckError(error, "clEnqueueNDRangeKernel");
  return gpu_event;
}

// Be sure to CreateSubdivide first
cl_event Gpu::EnqueueFindAmbiguous(
    int num_vertices, int max_level, GpuEvents wait_events) {
  cl_kernel k = kernel_find_ambiguous;
  clSetKernelArg(k, 0, sizeof(cl_mem), &vn_header);
  clSetKernelArg(k, 1, sizeof(cl_mem), &vertices);
  clSetKernelArg(k, 2, sizeof(cl_mem), &cpoints);
  clSetKernelArg(k, 3, sizeof(int), &max_level);
  clSetKernelArg(k, 4, sizeof(cl_mem), &to_subdivide);
  clSetKernelArg(k, 5, sizeof(cl_mem), &size);

  int error;
  size_t globalWorkSize[] = { static_cast<size_t>(num_vertices) };
  cl_event gpu_event;
  error = clEnqueueNDRangeKernel(
      queue, k, 1, nullptr, globalWorkSize,
      nullptr, wait_events.size(), wait_events.get(), &gpu_event);
  CheckError(error, "clEnqueueNDRangeKernel");
  return gpu_event;
}

NAMESPACE_OCT_END

#endif // __OPEN_CL_SUPPORT__
