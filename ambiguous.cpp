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

#include "./ambiguous.h"
#include "./search.h"
#include "./vertices_gpu_state.h"
#include <set>

namespace oct {

void EnsureSize(vector<vector<int> >& array, const int idx) {
  if (array.size() <= idx)
    array.resize(idx+1);
}

// template <int D>
// void VerticesGpuState<D>::InitializeVertices() {
//   if (!_vertices) return;

//   // typedef Vertex<D> Vertex_t;
//   typedef Vertex Vertex_t;

//   const int N = _vertices->size();
//   _v_array_size = N*2;
//   int error;

//   // shared_array<Vertex_t> vertex_array(new Vertex_t[_v_array_size]);
//   shared_array<Vertex_t> vertex_array(new Vertex_t[N]);
//   _vertices->CopyVertices(vertex_array.get());

//   // cl_vertex_array =
//   //     clCreateBuffer(context_cl, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
//   //                    sizeof(Vertex_t)*_v_array_size, vertex_array.get(),
//   //                    &error);
//   cl_vertex_array =
//       clCreateBuffer(context_cl, CL_MEM_READ_WRITE,
//                      sizeof(Vertex_t)*_v_array_size, 0,
//                      &error);
//   CheckError(error, "clCreateBuffer");

//   cl_command_queue queue_cl = CreateCommandQueue(context_cl, &error);
//   CheckError(error, "CreateCommandQueue");

//   error = clEnqueueWriteBuffer(
//       queue_cl, cl_vertex_array, CL_TRUE, 0,
//       sizeof(Vertex_t)*N, vertex_array.get(), 0, nullptr, nullptr);
//   CheckError(error, "clEnqueueWriteBuffer");

//   error = clFinish(queue_cl);
//   CheckError(error, "clFinish");
// }

// template <int D>
// void VerticesGpuState<D>::InitializeClosestPoints() {
//   if (!_vertices) return;

//   // typedef Vertex<D> Vertex_t;
//   typedef Vertex Vertex_t;

//   const int M = _vertices->NumClosestPoints();
//   int error;

//   shared_array<int> point_label_array(new int[M]);
//   // _vertices.CopyPointLabels(point_label_array.get());
//   for (int i = 0; i < M; ++i) {
//     point_label_array[i] = _vertices->LabelAtIdx(i);
//   }

//   // cl_point_label_array =
//   //     clCreateBuffer(context_cl, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
//   //                    sizeof(int)*M, point_label_array.get(), &error);
//   cl_point_label_array =
//       clCreateBuffer(context_cl, CL_MEM_READ_WRITE,
//                      sizeof(int)*M, 0, &error);
//   CheckError(error, "clCreateBuffer");

//   cl_command_queue queue_cl = CreateCommandQueue(context_cl, &error);
//   CheckError(error, "CreateCommandQueue");

//   error = clEnqueueWriteBuffer(
//       queue_cl, cl_point_label_array, CL_TRUE, 0,
//       sizeof(int)*M, point_label_array.get(), 0, nullptr, nullptr);
//   CheckError(error, "clEnqueueWriteBuffer");

//   error = clFinish(queue_cl);
//   CheckError(error, "clFinish");
// }

// Returns an array that gives ambiguities of each vertex.
//  0  unambiguous
//  1  ambiguous
// -1  error in computing ambiguity
// template <int D>
shared_array<int> GetAmbiguousGpu(
    const UVertexNetwork& vertices, int& ambiguous_count,
    VerticesGpuState<DIM>& gpu_state,
    const OctreeOptions& o) {
  // const int N = vertices.size();
  const int N = NumVertices(vertices);
  // const int M = vertices.NumClosestPoints();
  const int M = NumCPoints(vertices);

  shared_array<int> ambiguities(new int[N]);
  // shared_array<Vertex> vertex_array(new Vertex[N]);
  // vertices.CopyVertices(vertex_array.get());
  Vertex* vertex_array = vertices.vertices;
  shared_array<int> point_label_array(new int[M]);
  for (int i = 0; i < M; ++i) {
    // point_label_array[i] = vertices.LabelAtIdx(i);
    point_label_array[i] = vertices.cpoints[i].l;
  }

#ifdef __OPEN_CL_SUPPORT__

  // const cl_mem* cl_vertex_array = gpu_state.VertexArray();
  // const cl_mem* cl_point_label_array = gpu_state.PointLabelArray();

  // shared_array<Vertex_t> vertex_array(new Vertex_t[_v_array_size]);

  // typedef Vertex Vertex_t;
  int error;

  cl_mem cl_vertex_array =
      clCreateBuffer(context_cl, CL_MEM_READ_WRITE,
                     // sizeof(Vertex)*_v_array_size, 0,
                     sizeof(Vertex)*N, 0,
                     &error);
  cl_mem cl_point_label_array =
      clCreateBuffer(context_cl, CL_MEM_READ_WRITE,
                     sizeof(int)*M, 0, &error);
  cl_mem cl_ambiguities =
      clCreateBuffer(context_cl, CL_MEM_WRITE_ONLY,
                     sizeof(int)*N, nullptr, &error);

  cl_mem cl_ambiguous_count =
      clCreateBuffer(context_cl, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
                     sizeof(int), &ambiguous_count, &error);

  // Setup the kernel arguments
  clSetKernelArg(
      kernel_ambiguous_orig, 0, sizeof(cl_mem), &cl_vertex_array);
  clSetKernelArg(
      kernel_ambiguous_orig, 1, sizeof(cl_mem), &cl_point_label_array);
  clSetKernelArg(
      kernel_ambiguous_orig, 2, sizeof(cl_mem), &cl_ambiguities);
  clSetKernelArg(
      kernel_count_ambiguous, 0, sizeof(cl_mem), &cl_ambiguities);
  clSetKernelArg(
      kernel_count_ambiguous, 1, sizeof(int), &N);
  clSetKernelArg(
      kernel_count_ambiguous, 2, sizeof(cl_mem), &cl_ambiguous_count);
  
  cl_command_queue queue_cl = CreateCommandQueue(context_cl, &error);

  error = clEnqueueWriteBuffer(
      queue_cl, cl_vertex_array, CL_TRUE, 0,
      // sizeof(Vertex)*N, vertex_array.get(), 0, nullptr, nullptr);
      sizeof(Vertex)*N, vertex_array, 0, nullptr, nullptr);
  CheckError(error, "clEnqueueWriteBuffer");
  error = clEnqueueWriteBuffer(
      queue_cl, cl_point_label_array, CL_TRUE, 0,
      sizeof(int)*M, point_label_array.get(), 0, nullptr, nullptr);
  CheckError(error, "clEnqueueWriteBuffer");

  // Run the processing
  size_t globalWorkSize[] = { static_cast<size_t>(N) };
  size_t globalWorkSize1[] = { static_cast<size_t>(1) };
  cl_event gpu_event;
  error = clEnqueueNDRangeKernel(
      queue_cl, kernel_ambiguous_orig, 1, nullptr, globalWorkSize,
      nullptr, 0, nullptr, &gpu_event);
  CheckError(error, "clEnqueueNDRangeKernel");
  error = clEnqueueNDRangeKernel(
      queue_cl, kernel_count_ambiguous, 1, nullptr, globalWorkSize1,
      nullptr, 0, nullptr, &gpu_event);
  CheckError(error, "clEnqueueNDRangeKernel");

  // Execute
  clFinish(queue_cl);

  if (ambiguous_count > 0) {
    // Read result array
    clReleaseCommandQueue(queue_cl);
    queue_cl = CreateCommandQueue(context_cl, &error);
    // error = clEnqueueReadBuffer(
    //     queue_cl, cl_ambiguities, CL_TRUE, 0,
    //     sizeof(int)*N, ambiguities.get(), 1, &gpu_event, nullptr);
    error = clEnqueueReadBuffer(
        queue_cl, cl_ambiguities, CL_TRUE, 0,
        sizeof(int)*ambiguous_count, ambiguities.get(), 1, &gpu_event, nullptr);
    CheckError(error, "clEnqueueReadBuffer");
    clFinish(queue_cl);

    // Cleanup
    clReleaseMemObject(cl_ambiguities);
    clReleaseMemObject(cl_ambiguous_count);
    clReleaseCommandQueue(queue_cl);
  }

  // cout << "Checking for errors in ambiguity computation" << endl;
  // for (int vi = 0; vi < N; ++vi) {
  //   if (ambiguities[vi] == -1) {
  //     cerr << "ERROR computing ambiguity of cell " << vi
  //          << ".  Setting to ambiguous." << endl;
  //     ambiguities[vi] = 1;
  //   } else if (ambiguities[vi] < 0 || ambiguities[vi] > 1) {
  //     cerr << "ERROR: unexpected ambiguity of cell " << vi
  //          << ".  Setting to ambiguous." << endl;
  //     ambiguities[vi] = 1;
  //   }
  // }
  // cout << "Checking for errors in ambiguity computation" << endl;

  // int count = 0;
  // for (int vi = 0; vi < N; ++vi) {
  //   if (ambiguities[vi]) {
  //     ++count;
  //   }
  // }
  // cout << "count = " << count << " ambiguous_count = " << ambiguous_count
  //      << endl;

  // test for correctness
  // if (o.BoolValue("TEST_AMBIGUOUS_GPU", false)) {
  //   cout << "Testing ambiguous gpu implementation" << endl << flush;
  //   set<int> ambiguous_set(ambiguities.get(),
  //                          ambiguities.get()+ambiguous_count);
  //   vector<vector<int> > base2incident = ComputeBase2Incident<DIM>(vertices);
  //   for (int vi = 0; vi < N; ++vi) {
  //     if (vertices.IsBase(vi)) {
  //       const vector<int>& incident = base2incident[vi];
  //       const bool cpu = IsAmbiguous<DIM>(vi, incident, vertices, o);
  //       // const bool gpu = ambiguities[vi];
  //       const bool gpu = (ambiguous_set.find(vi) != ambiguous_set.end());
  //       if (cpu != gpu) {
  //         cout << "Difference in ambiguity.  vi = " << vi
  //              << " level = " << (int)vertices.CellLevel(vi)
  //              << " cpu = " << cpu << " gpu = " << gpu << endl;
  //       // } else if (gpu) {
  //       //   cout << "Correctly found gpu ambiguous.  vi = " << vi
  //       //        << " level = " << (int)vertices.CellLevel(vi) << endl;
  //       }
  //     }
  //   }
  // }
#endif //__OPEN_CL_SUPPORT__

  return ambiguities;
}

// check is initially filled with -1.
template <int D>
void SubdivideCell(
    const int vi, const intn& center,
    VertexNetwork* vertices,
    const OctreeOptions& o, IndexAndPoint<D>* check,
    vector<vector<int> >& base2incident,
    multiset<HeapVertex>& heap) {
  typedef Direction Direction_t;
  typedef Direction Dir;

  const level_t level = vertices->CellLevel(vi);
  const index_t width = Level2CellWidth(level);
  const index_t width_2 = (width >> 1);
  // static const int kNumSubdivided = kNumSubdividedArray[D];

  // Subdivide the cell
  int slvi2vi[kNumSubdivided];
  vertices->Subdivide(vi, slvi2vi);

  int idx = 0;
  set<int> incident_update(
      base2incident[vi].begin(), base2incident[vi].end());

  base2incident[vi].clear();
  for (int i = 0; i < kNumSubdivided; ++i) {
    const int new_vi = slvi2vi[i];
    // const Direction_t d = DirectionFromSlvi(i);
    // const intn& new_v_point = vertices->Position(new_vi);
    EnsureSize(base2incident, new_vi);
    incident_update.insert(new_vi);
  }

  const int* bases = vertices->GetCorners(vi);
  for (int lvi = 0; lvi < (1<<D); ++lvi) {
    const int new_vi = bases[lvi];
    const intn& new_v_point = vertices->Position(new_vi);
    if (level < o.ambiguous_max_level-1) {
      check[idx++] = IndexAndPoint<D>(new_vi, new_v_point);
    }
  }

  // Update base2incident
  for (set<int>::const_iterator it = incident_update.begin();
       it != incident_update.end(); ++it) {
    const int new_vi = *it;
    const intn new_p = vertices->Position(new_vi);
    vector<int> incident_bases;
    incident_bases.push_back(0);
    for (int dim = 0; dim < D; ++dim) {
      if (new_p.s[dim] == center.s[dim]) {
        const int n = incident_bases.size();
        for (int j = 0; j < n; ++j) {
          incident_bases.push_back(incident_bases[j] | (1<<dim));
        }
      } else if (new_p.s[dim] > center.s[dim]) {
        const int n = incident_bases.size();
        for (int j = 0; j < n; ++j) {
          incident_bases[j] = incident_bases[j] | (1<<dim);
        }
      }
    }
    for (int j = 0; j < incident_bases.size(); ++j) {
      base2incident[bases[incident_bases[j]]].push_back(new_vi);
    }
  }

  // Add appropriate vertices to the wavefront heap.
  // todo: this will become unnecessary if the code below works.
  for (int i = 0; i < kNumSubdivided; ++i) {
    const int new_vi = slvi2vi[i];
    const Direction_t d = DirectionFromSlvi(i);
    const intn new_v_point =
        // Position(d.Pos(), d.Neg(), center, width_2);
        Position(d.pos, d.neg, center, width_2);
    if (vertices->ClosestPointIndex(new_vi) == -1)
      LoadWavefront<D>(new_vi, new_v_point, *vertices, heap);
    else {
      const index_t dist =
          length(vertices->ClosestPoint(new_vi)-new_v_point);
      heap.insert(HeapVertex(dist, new_vi));
    }
  }
}

//------------------------------------------------------------------------------
// Parallel
//------------------------------------------------------------------------------

// A 1-edge array, or directed graph.  Vertices point only to
// their neighbors with larger indices.
//
// Suppose vis = { 1, 5, 7, 8, 11 }
//
//  8__________11
//  |          |
//  |          |
//  |          |7
//  |          |
//  |__________|
// 1            5
//
// Builds array such as
// 1:  { 5, 8 }
// 7:  { 11 }
// ...
//
// The actual representation will be more like
// 0:  { 1, 3 }
// 1:  { 0, -1 }
// ...
// (using indices into the incident array, not the vertex network).
//
// The array needs only be 2^(D-1) for each vertex.
// edge_arrays must already be allocated.
template <int D>
void BuildEdgeArray1(int* edge_arrays,
    const vector<int>& incident, const VertexNetwork& vertices) {
  // i, j, k \in { 0,1,...,incident.size } is the index into incident.
  // vi, vj, vk is the index into vertices

  const int m = incident.size();
  const int n = 1<<(D-1);
  fill(edge_arrays, edge_arrays+m*n, -1);
  for (int i = 0; i < m; ++i) {
    int* edge_array = edge_arrays + i*n;
    const int vi = incident[i];
    int count = 0;
    for (int j = i+1; j < m; ++j) {
      const int vj = incident[j];
      if (vertices.IsNeighbor(vi, vj)) {
        edge_array[count] = j;
        ++count;
      }
    }
  }
}

template <int D>
void PrintColors(
    const vector<int>& incident, const vector<int>& color,
    const VertexNetwork& vertices) {
  // const int n = 1<<(D-1);
  for (int i = 0; i < incident.size(); ++i) {
    const int vi = incident[i];
    cout << "  " << vi << ": " << color[i] << endl;
  }
}

template <int D>
void PrintEdgeArrays(
    const vector<int>& incident, const int* edge_arrays,
    const VertexNetwork& vertices) {
  // const int m = incident.size();
  const int n = 1<<(D-1);
  for (int i = 0; i < incident.size(); ++i) {
    const int vi = incident[i];
    cout << "  " << vi << ": ";
    for (int a = 0; a < n; ++a) {
      const int j = edge_arrays[i*n+a];
      if (j > -1) {
        const int vj = incident[j];
        if (vertices.Label(vi) == vertices.Label(vj))
          cout << "*";
        cout << vj << " ";
      }
    }
    cout << endl;
  }
}

// This is in-place and will modify arr.
int NumUnique(vector<int>& arr) {
  // In-place n^2 check
  const int n = arr.size();
  int count = 1;
  for (int i = 1; i < n; ++i) {
    const int value = arr[i];
    bool found = false;
    for (int j = 0; !found && j < count; ++j) {
      found = (arr[j] == value);
    }
    if (!found) {
      arr[count++] = value;
    }
  }
  return count;
}

template <int D>
vector<int> ColorComponents(const int base_vi,
    const vector<int>& incident, int* edge_arrays,
    const VertexNetwork& vertices, const OctreeOptions& o) {

  const int m = incident.size();
  const int n = 1<<(D-1);
  vector<bool> parent(m, true);

  // Color is different from the label.  Each vertex v starts uninitialized.
  // When v is visited in the DAG search it is assigned the color of the
  // origin vertex.  If v is already initialized, then it's color, together
  // with all other vertices that share that color, are changed to the origin
  // vertex's color.
  vector<int> colors(m, -1);

  for (int i = 0; i < m; ++i) {
    if (colors[i] == -1) {
      colors[i] = i;
    }
    const int ci = colors[i];
    const int vi = incident[i];
    int* edge_array = edge_arrays + i*n;
    // int target = -1;
    for (int a = 0; a < n; ++a) {
      const int j = edge_array[a];
      if (j > -1) {
        if (i >= j) throw logic_error("edge_array not ordered");
        const int vj = incident[j];
        if (vertices.Label(vi) == vertices.Label(vj)) {
          const int cj = colors[j];
          if (cj == -1) {
            // j's color is uninitialized
            colors[j] = ci;
          } else if (cj != ci) {
            // j has a different color than i.  Change all vertices that
            // share j's color to i's color.
            for (int k = 0; k < m; ++k) {
              if (colors[k] == cj) {
                colors[k] = ci;
              }
            }
          }
        }
      }
    }
  }
  return colors;
}

template <int D>
int GetLabelIndex(const int li, const int num_labels, vector<int>& label2idx) {
  int i = 0;
  for (; label2idx[i] > -1; ++i) {
    if (label2idx[i] == li) return i;
  }
  label2idx[i] = li;
  return i;
}

template <int D>
bool ReducedIsClique(
    const int num_labels, const vector<int>& incident,
    const VertexNetwork& vertices) {

  const int l = num_labels;
  vector<int> label2idx(l, -1);
  vector<bool> matrix(l*l, false);
  // Set diagonal to true.
  for (int i = 0; i < l; ++i) {
    matrix[i*l + i] = true;
  }

  const int m = incident.size();
  vector<unsigned char> parent(m, true);
  for (int i = 0; i < m && parent[i]; ++i) {
    const int vi = incident[i];
    // int count = 0;
    for (int j = i+1; j < m; ++j) {
      const int vj = incident[j];
      if (vertices.IsNeighbor(vi, vj) &&
          vertices.Label(vi) != vertices.Label(vj)) {
        const int li = vertices.Label(vi);
        const int lj = vertices.Label(vj);
        const int li_idx = GetLabelIndex<D>(li, num_labels, label2idx);
        const int lj_idx = GetLabelIndex<D>(lj, num_labels, label2idx);
        matrix[li_idx*l + lj_idx] = true;
        matrix[lj_idx*l + li_idx] = true;
      }
    }
  }

  for (int i = 0; i < l*l; ++i) {
    if (!matrix[i]) return false;
  }
  return true;
}

template <int D>
int NumLabels(
    const vector<int>& incident, const VertexNetwork& vertices) {
  // i, j, k \in { 0,1,...,incident.size } is the index into incident.
  // vi, vj, vk is the index into vertices

  const int m = incident.size();
  vector<int> labels(m, -1);
  int num_labels = 0;
  for (int i = 0; i < m; ++i) {
    const int vi = incident[i];
    const int label = vertices.Label(vi);
    if (label == -1) {
      cerr << "Invalid label in NumLabels.  vi = " << vi << endl;
      throw logic_error("Invalid label in NumLabels");
    }
    int added = false;
    for (int j = 0; !added; ++j) {
      if (j == m) {
        throw logic_error("Too many labels");
      }
      if (labels[j] > -1) {
        if (label == labels[j])
          added = true;
      } else {
        if (j != num_labels) throw logic_error("Bad counting in labels");
        labels[j] = label;
        ++num_labels;
        added = true;
      }
    }
  }
  return num_labels;
}

// incident is an array of vertex indices that all are incident to
// the cell with base vertex base_vi.
template <int D>
bool IsAmbiguous(
    const int base_vi, const vector<int>& incident,
    const VertexNetwork& vertices, const OctreeOptions& o) {

  // n is the maximum number of neighbors a vertex can have
  static const int n = 1<<(D-1);

  if (incident.size() < (1<<D)) {
    // If base_vi is not a base vertex then incident will be empty.
    // Otherwise there must be at least 2^D incident vertices.
    cerr << "error incident base_vi = " << base_vi << endl;
    cerr << "  incident.size() = " << incident.size() << endl;
    throw logic_error("error incident");
  }
  const int m = incident.size();
  shared_array<int> edge_arrays(new int[m * n]);
  BuildEdgeArray1<D>(edge_arrays.get(), incident, vertices);

  const int num_labels = NumLabels<D>(incident, vertices);
  vector<int> colors =
      ColorComponents<D>(base_vi, incident, edge_arrays.get(), vertices, o);
  const int num_components = NumUnique(colors);
  if (num_components > num_labels) {
    // If a label is duplicated but is in a different connected component
    // then the reduced graph cannot be a clique -- there is no edge
    // between vertices with identical labels.
    return true;
  }

  // Build a matrix of label-label edges.  Must be a clique.
  if (!ReducedIsClique<D>(num_labels, incident, vertices)) {
    return true;
  }

  return false;
}

void SubdivideAmbiguousCpu(VertexNetwork& vertices, const OctreeOptions& o) {
  int iter = 0;
  // Add all vertices to the to_check queue.
  vector<int> to_check(vertices.size());
  for (int base_vi = 0; base_vi < vertices.size(); ++base_vi) {
    to_check[base_vi] = base_vi;
  }
  int num_subdivided = 0;
  while (!to_check.empty()) {
    iter++;
    // cout << "Ambiguous iteration " << iter << endl;

    // todo: remove expensive recomputation.  Except that my base2incident
    // update appears to be incorrect.  Recomputing base2incident every time
    // fixed it.
    vector<vector<int> > base2incident = ComputeBase2Incident<DIM>(vertices);

    // Iterate through cells and compute ambiguity of each cell
    vector<int> ambiguous;
    const int N = vertices.size();
    for (int base_vi = 0; base_vi < N; ++base_vi) {
      if (vertices.IsBase(base_vi)
          && vertices.CellLevel(base_vi) < o.ambiguous_max_level) {
        const vector<int>& incident = base2incident[base_vi];
        if (IsAmbiguous<DIM>(base_vi, incident, vertices, o)) {
          ambiguous.push_back(base_vi);
        }
      }
    }
    to_check.clear();

    num_subdivided += ambiguous.size();
    multiset<HeapVertex> heap;
    for (int i = 0; i < ambiguous.size(); ++i) {
      const int vi = ambiguous[i];

      const level_t level = vertices.CellLevel(vi);
      const index_t width = Level2CellWidth(level);
      const index_t width_2 = (width >> 1);
      intn center = vertices.Position(vi) + width_2;

      IndexAndPoint<DIM> check[27];
      fill(check, check+27, IndexAndPoint<DIM>());
      SubdivideCell<DIM>(
          vi, center, &vertices, o, check, base2incident, heap);
      for (int j = 0; j < 27 && check[j].vi != -1; ++j) {
        to_check.push_back(check[j].vi);
      }
    }
    SetDistances(vertices, heap, o);
  }
  if (o.BoolValue("AMBIGUOUS_OUTPUT", false)) {
    cout << "  Ambiguous (cpu) num subdivided: " << num_subdivided << endl;
    cout << "  Ambiguous (cpu) num iterations: " << iter << endl;
  }
}

// template <int D>
void SubdivideAmbiguousGpu(
    MVertexNetwork& mvertices, VerticesGpuState<DIM>& gpu_state,
    const OctreeOptions& o) {
  // static const int kNumSubdivided = kNumSubdividedArray[DIM];

  bool subdivide = true;
  int iter = 0;
  int num_subdivided = 0;
  while (subdivide) {
    subdivide = false;
    ++iter;
    vector<int> ambiguous;
    int ambiguous_count;
    int scount = 0;
    UVertexNetwork vertices = make_vertex_network(mvertices);
    shared_array<int> ambiguities = GetAmbiguousGpu(
        vertices, ambiguous_count, gpu_state, o);
    for (int i = 0; i < ambiguous_count; ++i) {
      const int vi = ambiguities[i];
      // if (vertices.CellLevel(vi) < o.ambiguous_max_level) {
      if (CellLevel(vi, vertices) < o.ambiguous_max_level) {
        ambiguous.push_back(ambiguities[i]);
        subdivide = true;
        ++num_subdivided;
        ++scount;
        assert(vi < NumVertices(vertices));
      }
    }
    if (o.BoolValue("AMBIGUOUS_OUTPUT", false)) {
      cout << "    Ambiguous (gpu) count: " << scount
           << " (iteration " << iter << ")"
           << endl;
    }

    // mvertices = make_mvertex_network(vertices);
    update_mvertex_network(vertices, mvertices);
    // mvn_ensure(scount, &mvertices);
    const int threshold = NumVertices(mvertices) + scount * kSubAdded;
    const int new_n =
        NumVertices(mvertices) + scount * kSubAdded * o.verts_alloc_factor;
    mvn_ensure(threshold, new_n, &mvertices);
    vertices = make_vertex_network(mvertices);

    if (subdivide) {
      multiset<HeapVertex> heap;
      int slvi2vi[kNumSubdivided];
      for (int i = 0; i < ambiguous.size(); ++i) {
        const int vi = ambiguous[i];
        // Do the subdivide
        // shared_array<int> subcells(new int[1<<DIM]);
        // Subdivide(vi, subcells.get(), slvi2vi, &vertices);
        Subdivide(vi, slvi2vi, &vertices);
        shared_array<int> subcells(new int[1<<DIM]);
        memcpy(subcells.get(), vertices.vertices[vi].corners,
               sizeof(int)*(1<<DIM));
        // Add incident vertices to heap
        for (int j = 0; j < kNumSubdivided; ++j) {
          const int subvi = slvi2vi[j];
          // if (vertices.ClosestPointIndex(subvi) != -1) {
          if (ClosestPointIndex(subvi, vertices) != -1) {
            // heap.insert(HeapVertex(vertices.Dist(subvi), subvi));
            heap.insert(HeapVertex(Dist(subvi, vertices), subvi));
          }
        }
      }
      SetDistances(vertices, heap, o);
    }
    // mvertices = make_mvertex_network(vertices);
    update_mvertex_network(vertices, mvertices);
  }  
  if (o.BoolValue("AMBIGUOUS_OUTPUT", false)) {
    cout << "  Ambiguous (gpu) num subdivided: " << num_subdivided << endl;
    cout << "  Ambiguous (gpu) num iterations: " << iter << endl;
  }

  // mvertices = make_mvertex_network(vertices);
}

// template <int D>
// void SubdivideAmbiguous(
//     VertexNetwork& vertices, VerticesGpuState<D>& gpu_state,
//     const OctreeOptions& o) {
//   if (o.gpu && o.BoolValue("TEST", false)) {
//     SubdivideAmbiguousGpu(vertices, gpu_state, o);
//   } else {
//     SubdivideAmbiguousCpu(vertices, o);
//   }
// }

//------------------------------------------------------------------------------
// FindBasesVisitor
//------------------------------------------------------------------------------

template <int N>
bool FindBasesVisitor<N>::operator()(
    const int vi, const intn& p) {
  typedef oct::Direction Direction_t;

  // Important variables
  //   axis_p := axis in subspace (in the plane if D==2)
  //   axis_n := axis in native space

  if (_visited[vi]) return true;

  // Make sure all neighbors in negative cardinal directions are complete.
  // Also find a direction dir that has a neighbor.
  int n_dir = -1;
  for (int axis_p = 1; axis_p < M; axis_p=(axis_p<<1)) {
    const int axis_n = _p2n[axis_p];
    const Direction_t d = DirectionFromPosNeg(0, axis_n);
    const int n_vi = _vertices->Neighbor(vi, d);
    if (n_vi != -1) {
      if (!_visited[n_vi]) {
        intn n_p = p;
        // n_p[axis_p] -= _vertices->NeighborDist(vi, d);
        n_p.s[AXIS_IDX[axis_n]] -= _vertices->NeighborDist(vi, d);
        (*this)(n_vi, n_p);
      }
      n_dir = axis_n;
    }
  }
  if (n_dir == -1 && vi != 0 && D == N) {
    cerr << "No neighbor in a negative direction.  vi = " << vi << endl;
    throw logic_error("No neighbor in a negative direction?");
  }

  // 0 direction will be set only if vi is the base of a cell.
  if (IsBase(vi)) {
    _bases[vi*M] = vi;
  }

  // Process cardinal directions
  for (int axis_p = 1; axis_p < M; axis_p=(axis_p<<1)) {
    if (!Valid2Direction(vi, p, axis_p))
      continue;
    const int axis_n = _p2n[axis_p];
    const int axis_n_idx = AXIS_IDX[_p2n[axis_p]];
    const Direction_t d = DirectionFromPosNeg(0, axis_n);
    const int n_vi = _vertices->Neighbor(vi, d);
    if (n_vi != -1) {
      const int dist = _vertices->NeighborDist(vi, d);
      if (IsBase(n_vi)) {
        _bases[vi*M+axis_p] = n_vi;
        intn d = make_intn(0);
        d.s[axis_n_idx] = dist;
        _distances[vi*M+axis_p] = d;
      } else if (!IsBoundary(n_vi, axis_n)) {
        _bases[vi*M+axis_p] = _bases[n_vi*M+axis_p];
        intn d = _distances[n_vi*M+axis_p];
        d.s[axis_n_idx] += dist;
        _distances[vi*M+axis_p] = d;
      }
    }    
  }

  // 2D diagonals
  // When D == 2, Iterate 1 time: 01/10
  // When D == 3, Iterate 3 times: 001/010, 001/100, 010/100
  for (int i_p = 1; i_p < M; i_p = i_p<<1) {
    for (int j_p = i_p<<1; j_p < M; j_p = j_p<<1) {
      // Find one that has a neighbor in that direction
      int iaxis_p = i_p;
      int jaxis_p = j_p;
      int iaxis_n = _p2n[i_p];
      int jaxis_n = _p2n[j_p];
      const int diag_p = iaxis_p | jaxis_p;
      const int diag_n = iaxis_n | jaxis_n;
      if (diag_n != _p2n[diag_p]) throw logic_error("diag space bad");
      int n_vi = _vertices->Neighbor(vi, DirectionFromPosNeg(0, iaxis_n));
      int dist_to_n = _vertices->NeighborDist(vi, DirectionFromPosNeg(0, iaxis_n));
      if (n_vi == -1) {
        swap(iaxis_n, jaxis_n);
        swap(iaxis_p, jaxis_p);
        n_vi = _vertices->Neighbor(vi, DirectionFromPosNeg(0, iaxis_n));
        dist_to_n = _vertices->NeighborDist(vi, DirectionFromPosNeg(0, iaxis_n));
      }
      if (n_vi != -1) {
        // n_vi is in the iaxis_n direction
        int base = _bases[n_vi*M+jaxis_p];
        intn dist_n_to_base = _distances[n_vi*M+jaxis_p];
        if (base == -1 && !IsBoundary2(n_vi, iaxis_n) &&
            (D == N || !_vertices->HasNeighbor(n_vi, DirectionFromPosNeg(0, jaxis_n)))) {
          base = _bases[n_vi*M+diag_p];
          dist_n_to_base = _distances[n_vi*M+diag_p];
        }
        if (base != -1) {
          _bases[vi*M+diag_p] = base;
          intn d = dist_n_to_base;
          d.s[AXIS_IDX[iaxis_n]] += dist_to_n;
          _distances[vi*M+diag_p] = d;
        }
      } else {
        if (D == 2 && vi != 0 && D == N) {
          cerr << "In 2D at least one neighbor must exist.  vi = " << vi
               << " i = " << iaxis_n << " j = " << jaxis_n << endl;
          throw logic_error("In 2D at least one neighbor must exist");
        }
      }
    }
  }

  if (D == 3 && vi != 0) {
    // 3D direction
    const int diag_n = 7;
    const int n_vi = _vertices->Neighbor(vi, DirectionFromPosNeg(0, n_dir));
    if (n_vi == -1) throw logic_error("n_vi unexpectedly -1");
    // int base = Get(n_vi, diag_n - n_dir);
    int base = _bases[n_vi*M + (diag_n - n_dir)];
    if (base == -1) {
      // base = Get(n_vi, diag_n);
      base = _bases[n_vi*M + diag_n];
    }
    if (base != -1) {
      _bases[vi*M+diag_n] = base;
    }
  }

  _visited[vi] = true;
  return true;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Template instantiations
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifdef OCT2D
template bool IsAmbiguous<2>(
    const int base_vi, const vector<int>& incident,
    const VertexNetwork& vertices, const OctreeOptions& o);
// template void SubdivideAmbiguous(
//     VertexNetwork& vertices,
//     VerticesGpuState<2>& gpu_state,
//     const OctreeOptions& o);
template bool FindBasesVisitor<2>::operator()(
    const int vi, const int2& p);
#endif

#ifdef OCT3D
template bool IsAmbiguous<3>(
    const int base_vi, const vector<int>& incident,
    const VertexNetwork& vertices, const OctreeOptions& o);
template bool FindBasesVisitor<3>::operator()(
    const int vi, const int3& p);
#endif

}
