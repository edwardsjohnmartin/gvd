#ifdef __OPEN_CL_SUPPORT__

#ifdef __APPLE__
#include "OpenCL/opencl.h"
#else
#include "CL/cl.h"
#endif

#include "./options.h"

#ifndef nullptr
#define nullptr 0
#endif

extern cl_context context_cl;
extern cl_program program_cl;
extern cl_kernel kernel_init_wave;
extern cl_kernel kernel_ambiguous_orig;
extern cl_kernel kernel_update_vertices;
extern cl_kernel kernel_count_ambiguous;

extern cl_kernel kernel_find_nbrs;
extern cl_kernel kernel_find_to_subdivide1;
extern cl_kernel kernel_find_to_subdivide2;
extern cl_kernel kernel_count_to_subdivide;
extern cl_kernel kernel_copy_vertex_array;
extern cl_kernel kernel_compute_sparse_vi2geometries_size;
extern cl_kernel kernel_make_sparse_vi2geometries;
extern cl_kernel kernel_compute_dense_vi2geometries_size;
extern cl_kernel kernel_condense_vi2geometries;
extern cl_kernel kernel_subdivide_cell;
extern cl_kernel kernel_subdivide_cell_A;
extern cl_kernel kernel_subdivide_cell_B;
extern cl_kernel kernel_subdivide_cell_C;
extern cl_kernel kernel_clip_geometries;
extern cl_kernel kernel_compute_non_empty_vertex_distances;
extern cl_kernel kernel_pull_distances;
extern cl_kernel kernel_find_ambiguous;

std::string LoadKernel (const char* name);
cl_program CreateProgram(const std::string& source, cl_context context);
cl_int BuildProgram(cl_program program, bool print_log);
cl_command_queue CreateCommandQueue(cl_context context, int* error);

size_t GetKernelWorkgroupSize(cl_kernel kernel);

void CheckError (cl_int error);
void CheckError(cl_int error, const std::string& call);

void OpenCLInit(const int D, oct::OctreeOptions o, const bool safe_log);
void OpenCLCleanup();

#endif // __OPEN_CL_SUPPORT__
