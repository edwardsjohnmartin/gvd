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

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>

#include "./opencl.h"
#include "./shared_array.h"

using namespace std;

cl_context context_cl;
cl_uint deviceIdCount;
vector<cl_device_id> deviceIds;

// cl_context context;
cl_program program_cl;
cl_kernel kernel_init_wave;
cl_kernel kernel_ambiguous_orig;
cl_kernel kernel_update_vertices;
cl_kernel kernel_count_ambiguous;

cl_kernel kernel_find_nbrs;
cl_kernel kernel_find_to_subdivide1;
cl_kernel kernel_find_to_subdivide2;
cl_kernel kernel_count_to_subdivide;
cl_kernel kernel_copy_vertex_array;
cl_kernel kernel_compute_sparse_vi2geometries_size;
cl_kernel kernel_make_sparse_vi2geometries;
cl_kernel kernel_compute_dense_vi2geometries_size;
cl_kernel kernel_condense_vi2geometries;
cl_kernel kernel_subdivide_cell;
cl_kernel kernel_subdivide_cell_A;
cl_kernel kernel_subdivide_cell_B;
cl_kernel kernel_subdivide_cell_C;
cl_kernel kernel_clip_geometries;
cl_kernel kernel_compute_non_empty_vertex_distances;
cl_kernel kernel_pull_distances;
cl_kernel kernel_find_ambiguous;

cl_context GetCLContext() { return context_cl; }

std::string GetPlatformName (cl_platform_id id)
{
  size_t size = 0;
  clGetPlatformInfo (id, CL_PLATFORM_NAME, 0, nullptr, &size);

  std::string result;
  result.resize (size);
  clGetPlatformInfo (id, CL_PLATFORM_NAME, size,
                     const_cast<char*> (result.data ()), nullptr);

  return result;
}

std::string GetDeviceName (cl_device_id id)
{
  size_t size = 0;
  clGetDeviceInfo (id, CL_DEVICE_NAME, 0, nullptr, &size);

  std::string result;
  result.resize (size);
  clGetDeviceInfo (id, CL_DEVICE_NAME, size,
                   const_cast<char*> (result.data ()), nullptr);

  return result;
}

std::string GetDeviceVersion(cl_device_id id)
{
  size_t size = 0;
  clGetDeviceInfo (id, CL_DEVICE_OPENCL_C_VERSION, 0, nullptr, &size);

  std::string result;
  result.resize (size);
  clGetDeviceInfo (id, CL_DEVICE_OPENCL_C_VERSION, size,
                   const_cast<char*> (result.data ()), nullptr);

  return result;
}

cl_uint GetDeviceMaxComputeUnits(cl_device_id id)
{
  // size_t size = 0;
  // clGetDeviceInfo (id, CL_DEVICE_MAX_COMPUTE_UNITS, 0, nullptr, &size);

  cl_uint result;
  // result.resize (size);
  clGetDeviceInfo (id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint),
                   &result, nullptr);

  return result;
}

vector<size_t> GetDeviceMaxWorkItemSizes(cl_device_id id) {
  cl_uint D;
  clGetDeviceInfo (id, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint),
                   &D, nullptr);

  size_t result[32];
  clGetDeviceInfo (id, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t)*D,
                   result, nullptr);

  vector<size_t> ret;
  for (int i = 0; i < D; ++i) {
    ret.push_back(result[i]);
  }

  return ret;
}

cl_uint GetDeviceMaxWorkGroupSize(cl_device_id id)
{
  cl_uint result;
  clGetDeviceInfo (id, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(cl_uint),
                   &result, nullptr);

  return result;
}

bool GetDeviceLittleEndian(cl_device_id id)
{
  cl_bool result;
  clGetDeviceInfo(id, CL_DEVICE_ENDIAN_LITTLE, sizeof(cl_uint),
                   &result, nullptr);

  return (result == CL_TRUE);
}

// size_t GetKernelWorkgroupSize(cl_kernel kernel) {
// cl_int clGetKernelWorkGroupInfo (cl_kernel kernel, cl_device_id device,
// cl_kernel_work_group_info param_name, size_t param_value_size,
// void *param_value,
// size_t *param_value_size_ret)

string ErrorString(cl_int err) {
  switch (err) {
    case CL_SUCCESS:
      return "Success!";
    case CL_DEVICE_NOT_FOUND:
      return "Device not found.";
    case CL_DEVICE_NOT_AVAILABLE:
      return "Device not available";
    case CL_COMPILER_NOT_AVAILABLE:
      return "Compiler not available";
    case CL_MEM_OBJECT_ALLOCATION_FAILURE:
      return "Memory object allocation failure";
    case CL_OUT_OF_RESOURCES:
      return "Out of resources";
    case CL_OUT_OF_HOST_MEMORY:
      return "Out of host memory";
    case CL_PROFILING_INFO_NOT_AVAILABLE:
      return "Profiling information not available";
    case CL_MEM_COPY_OVERLAP:
      return "Memory copy overlap";
    case CL_IMAGE_FORMAT_MISMATCH:
      return "Image format mismatch";
    case CL_IMAGE_FORMAT_NOT_SUPPORTED:
      return "Image format not supported";
    case CL_BUILD_PROGRAM_FAILURE:
      return "Program build failure";
    case CL_MAP_FAILURE:
      return "Map failure";
    case CL_INVALID_VALUE:
      return "Invalid value";
    case CL_INVALID_DEVICE_TYPE:
      return "Invalid device type";
    case CL_INVALID_PLATFORM:
      return "Invalid platform";
    case CL_INVALID_DEVICE:
      return "Invalid device";
    case CL_INVALID_CONTEXT:
      return "Invalid context";
    case CL_INVALID_QUEUE_PROPERTIES:
      return "Invalid queue properties";
    case CL_INVALID_COMMAND_QUEUE:
      return "Invalid command queue";
    case CL_INVALID_HOST_PTR:
      return "Invalid host pointer";
    case CL_INVALID_MEM_OBJECT:
      return "Invalid memory object";
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
      return "Invalid image format descriptor";
    case CL_INVALID_IMAGE_SIZE:
      return "Invalid image size";
    case CL_INVALID_SAMPLER:
      return "Invalid sampler";
    case CL_INVALID_BINARY:
      return "Invalid binary";
    case CL_INVALID_BUILD_OPTIONS:
      return "Invalid build options";
    case CL_INVALID_PROGRAM:
      return "Invalid program";
    case CL_INVALID_PROGRAM_EXECUTABLE:
      return "Invalid program executable";
    case CL_INVALID_KERNEL_NAME:
      return "Invalid kernel name";
    case CL_INVALID_KERNEL_DEFINITION:
      return "Invalid kernel definition";
    case CL_INVALID_KERNEL:
      return "Invalid kernel";
    case CL_INVALID_ARG_INDEX:
      return "Invalid argument index";
    case CL_INVALID_ARG_VALUE:
      return "Invalid argument value";
    case CL_INVALID_ARG_SIZE:
      return "Invalid argument size";
    case CL_INVALID_KERNEL_ARGS:
      return "Invalid kernel arguments";
    case CL_INVALID_WORK_DIMENSION:
      return "Invalid work dimension";
    case CL_INVALID_WORK_GROUP_SIZE:
      return "Invalid work group size";
    case CL_INVALID_WORK_ITEM_SIZE:
      return "Invalid work item size";
    case CL_INVALID_GLOBAL_OFFSET:
      return "Invalid global offset";
    case CL_INVALID_EVENT_WAIT_LIST:
      return "Invalid event wait list";
    case CL_INVALID_EVENT:
      return "Invalid event";
    case CL_INVALID_OPERATION:
      return "Invalid operation";
    case CL_INVALID_GL_OBJECT:
      return "Invalid OpenGL object";
    case CL_INVALID_BUFFER_SIZE:
      return "Invalid buffer size";
    case CL_INVALID_MIP_LEVEL:
      return "Invalid mip-map level";
    default:
      return "Unknown";
  }
}

void CheckError (cl_int error)
{
  if (error != CL_SUCCESS) {
    cerr << "OpenCL call failed with error " << error << ": "
         << ErrorString(error) << endl;
    exit(1);
  }
}

void CL_CALLBACK callback(const char *errinfo,
                          const void *private_info,
                          size_t cb, void *user_data) {
  cerr << "cl_callback errinfo = " << errinfo << endl;
}

void CheckError(cl_int error, const string& call)
{
  if (error != CL_SUCCESS) {
    cerr << "OpenCL call " << call << " failed with error " << error << ": "
         << ErrorString(error) << endl;
    exit(1);
  }
}

std::string LoadKernel (const char* name)
{
  ifstream in (name);
  if (!in) {
    cerr << "Failed to open " << name << endl;
    exit(1);
  }
  string result(
      (istreambuf_iterator<char>(in)),
      istreambuf_iterator<char>());
  return result;
}

// cl_program CreateProgram(const std::string& source, cl_context context) {
cl_program CreateProgram(cl_context context) {

  vector<string> files;
  files.push_back("./opencl/octree.cl");
  files.push_back("./opencl/opencl_cpps.cl");

  const int count = files.size();
  oct::shared_array<string> ssources(new string[count]);
  oct::shared_array<const char*> sources(new const char*[count]);
  oct::shared_array<size_t> lengths(new size_t[count]);

  for (int i = 0; i < count; ++i) {
    const string k = LoadKernel(files[i].c_str());
    ssources[i] = k;
    sources[i] = ssources[i].data();
    lengths[i] = k.size();
  }

  cl_int error = 0;
  cl_program program = clCreateProgramWithSource(
      context, count, sources.get(), lengths.get(), &error);
  CheckError(error, "clCreateProgramWithSource");
  return program;
}

cl_int BuildProgram(cl_program program, bool print_log) {
#ifdef OCT2D
  const int error = clBuildProgram (
      program, deviceIdCount, deviceIds.data(),
      "-D FILTER_SIZE=1 -D OPEN_CL -D DIM=2", nullptr, nullptr);
#else
  const int error = clBuildProgram (
      program, deviceIdCount, deviceIds.data(),
      "-D FILTER_SIZE=1 -D OPEN_CL -D DIM=3", nullptr, nullptr);
#endif
  if (error != CL_SUCCESS && print_log) {
    // First call to know the proper size
    size_t log_size;
    int err = clGetProgramBuildInfo(
        // program, devices_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
        program, deviceIds[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
    char* build_log = (char*)malloc(log_size+1);

    // Second call to get the log
    err = clGetProgramBuildInfo(
        program, deviceIds[0], CL_PROGRAM_BUILD_LOG, log_size, build_log, NULL);
    build_log[log_size] = '\0';
    printf("--- Build log ---\n ");
    fprintf(stderr, "%s\n", build_log);

    free(build_log);
  }

  return error;
}

cl_command_queue CreateCommandQueue(cl_context context, int* error) {
  return clCreateCommandQueue(context, deviceIds[0], 0, error);
}

bool is_big_endian() {
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};
    return bint.c[0] == 1; 
}

cl_kernel CreateKernel(const string& name) {
  int err;
  cl_kernel k =
      clCreateKernel(program_cl, name.c_str(), &err);
  CheckError(err, name);
  return k;
}

// bool opencl_initialized = false;
void OpenCLInit(const int D, oct::OctreeOptions o, const bool safe_log) {
  // if (opencl_initialized) return;
  // opencl_initialized = true;

  // Discover our platform
  cl_uint platformIdCount = 0;
  clGetPlatformIDs (0, nullptr, &platformIdCount);

  if (platformIdCount == 0) {
    std::cerr << "No OpenCL platform found" << std::endl;
    // return 1;
    throw logic_error("No OpenCL platform found");
  } else {
    cout << "Found " << platformIdCount << " platform(s)" << std::endl;
  }

  std::vector<cl_platform_id> platformIds (platformIdCount);
  clGetPlatformIDs (platformIdCount, platformIds.data(), nullptr);

  for (cl_uint i = 0; i < platformIdCount; ++i) {
    cout << "\t (" << (i+1) << ") : "
         << GetPlatformName (platformIds [i]) << std::endl;
  }

  cl_device_type device_type = CL_DEVICE_TYPE_GPU;
  // if (o.Value("CL_DEVICE_TYPE", "GPU") == "CPU") {
  if (o.Value("CL_DEVICE_TYPE", "CPU") == "CPU") {
    device_type = CL_DEVICE_TYPE_CPU;
  }

  // Get the available devices
  // cl_uint deviceIdCount = 0;
  deviceIdCount = 0;
  // clGetDeviceIDs (platformIds [0], CL_DEVICE_TYPE_ALL, 0, nullptr,
  //                 &deviceIdCount);
  // clGetDeviceIDs (platformIds [0], CL_DEVICE_TYPE_GPU, 0, nullptr,
  //                 &deviceIdCount);
  clGetDeviceIDs (platformIds [0], device_type, 0, nullptr,
                  &deviceIdCount);

  if (deviceIdCount == 0) {
    std::cerr << "No OpenCL devices found" << std::endl;
    // return 1;
    throw logic_error("No OpenCL devices found");
  } else {
    cout << "Found " << deviceIdCount << " device(s)" << std::endl;
  }

  // vector<cl_device_id> deviceIds (deviceIdCount);
  deviceIds.resize(deviceIdCount);
  // clGetDeviceIDs (platformIds [0], CL_DEVICE_TYPE_ALL, deviceIdCount,
  //                 deviceIds.data (), nullptr);
  // clGetDeviceIDs (platformIds [0], CL_DEVICE_TYPE_GPU, deviceIdCount,
  //                 deviceIds.data (), nullptr);
  clGetDeviceIDs (platformIds [0], device_type, deviceIdCount,
                  deviceIds.data (), nullptr);

  for (cl_uint i = 0; i < deviceIdCount; ++i) {
    cout << "\t (" << (i+1) << ") : " 
         << GetDeviceName (deviceIds [i]) << endl;
    cout << "\t\t" << GetDeviceVersion(deviceIds[i]) << endl;
    cout << "\t\tMaxComputeUnits = " << GetDeviceMaxComputeUnits(deviceIds[i]) 
         << endl;
    vector<size_t> wis = GetDeviceMaxWorkItemSizes(deviceIds[i]);
    cout << "\t\tMaxWorkItemSizes = ";
    for (int j = 0; j < wis.size(); ++j) {
      cout << wis[j] << " ";
    }
    cout << endl;
    cout << "\t\tMaxWorkGroupSize = "
         << GetDeviceMaxWorkGroupSize(deviceIds[i]) 
         << endl;
    cout << "\t\tLittleEndian = "
         << GetDeviceLittleEndian(deviceIds[i]) 
         << " (cpu = " << !is_big_endian() << ")"
         << endl;
  }

  if (GetDeviceLittleEndian(deviceIds[0]) != !is_big_endian()) {
    cerr << "Endianness of cpu and gpu differ.  This has not been tested."
         << endl;
    throw logic_error(
        "Endianness of cpu and gpu differ.  This has not been tested.");
  }

  // Set up a context
  const cl_context_properties contextProperties [] = {
    CL_CONTEXT_PLATFORM,
    reinterpret_cast<cl_context_properties> (platformIds [0]),
    0, 0
  };

  cl_int error = CL_SUCCESS;
  context_cl = clCreateContext(
      contextProperties, deviceIdCount,
      deviceIds.data (), callback, nullptr, &error);
  CheckError (error);
  cout << "Context created" << endl;

  // Load program and set up kernels
  time_t cl_time = time(0);
  program_cl = CreateProgram(context_cl);
  cl_int err = BuildProgram(program_cl, safe_log);
  CheckError(err, "clBuildProgram");

  // Create kernels
  if (D == 2) {
    kernel_init_wave = CreateKernel("ComputeLeafVertexDistances2");
  } else {
    kernel_init_wave = CreateKernel("ComputeLeafVertexDistances3");
  }
  kernel_ambiguous_orig = CreateKernel("k_GetAmbiguous3_orig");
  kernel_update_vertices = CreateKernel("UpdateVertices3");
  kernel_count_ambiguous = CreateKernel("CountAmbiguous");

  kernel_find_nbrs = CreateKernel("k_FindNbrs");
  kernel_find_to_subdivide1 = CreateKernel("k_FindToSubdivide1");
  kernel_find_to_subdivide2 = CreateKernel("k_FindToSubdivide2");
  kernel_count_to_subdivide = CreateKernel("k_CountToSubdivide");
  kernel_copy_vertex_array = CreateKernel("k_CopyVertexArray");
  kernel_compute_sparse_vi2geometries_size =
      CreateKernel("k_ComputeSparseVi2GeometriesSize");
  kernel_make_sparse_vi2geometries = CreateKernel("k_MakeSparseVi2Geometries");
  kernel_compute_dense_vi2geometries_size =
      CreateKernel("k_ComputeDenseVi2GeometriesSize");
  kernel_condense_vi2geometries = CreateKernel("k_CondenseVi2Geometries");
  kernel_subdivide_cell = CreateKernel("k_SubdivideCell");
  kernel_subdivide_cell_A = CreateKernel("k_SubdivideCell_A");
  kernel_subdivide_cell_B = CreateKernel("k_SubdivideCell_B");
  kernel_subdivide_cell_C = CreateKernel("k_SubdivideCell_C");
  kernel_clip_geometries = CreateKernel("k_ClipGeometries");
  kernel_compute_non_empty_vertex_distances =
      CreateKernel("k_ComputeNonEmptyVertexDistances");
  kernel_pull_distances = CreateKernel("k_PullDistances");
  kernel_find_ambiguous = CreateKernel("k_FindAmbiguous");

  // output
  cout << "OpenCL kernel initialization time: "
       << difftime(time(0), cl_time) << endl;
}

void OpenCLCleanup() {
  clReleaseKernel(kernel_init_wave);
  clReleaseKernel(kernel_ambiguous_orig);
  clReleaseKernel(kernel_count_ambiguous);
  clReleaseKernel(kernel_update_vertices);

  clReleaseKernel(kernel_find_nbrs);
  clReleaseKernel(kernel_find_to_subdivide1);
  clReleaseKernel(kernel_find_to_subdivide2);
  clReleaseKernel(kernel_count_to_subdivide);
  clReleaseKernel(kernel_copy_vertex_array);
  clReleaseKernel(kernel_compute_sparse_vi2geometries_size);
  clReleaseKernel(kernel_make_sparse_vi2geometries);
  clReleaseKernel(kernel_compute_dense_vi2geometries_size);
  clReleaseKernel(kernel_condense_vi2geometries);
  clReleaseKernel(kernel_subdivide_cell);
  clReleaseKernel(kernel_subdivide_cell_A);
  clReleaseKernel(kernel_subdivide_cell_B);
  clReleaseKernel(kernel_subdivide_cell_C);
  clReleaseKernel(kernel_clip_geometries);
  clReleaseKernel(kernel_compute_non_empty_vertex_distances);
  clReleaseKernel(kernel_pull_distances);
  clReleaseKernel(kernel_find_ambiguous);

  clReleaseProgram(program_cl);
  clReleaseContext(context_cl);
}

#endif //__OPEN_CL_SUPPORT__
