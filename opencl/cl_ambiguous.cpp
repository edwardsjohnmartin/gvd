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

#include "./vertex.h"
#include "./vec.h"
#include "./uvertex_network.h"

NAMESPACE_OCT_BEGIN

// index and location
typedef struct {
  int idx;
  int3 p;
  int label;
} LVertex;

// The maximum supported number of vertices on a cell
#define MAX_SIZE 128
// The maximum supported number of labels on vertices on a cell
#define MAX_LABELS 16

//--------------------------------------------------
//--------------------------------------------------
// Declarations
//--------------------------------------------------
//--------------------------------------------------
bool Inside(const int3 p, const int width);
int FindLVertex(const LVertex lv, LVertex* array, const int size);
int FindLVertexWithLabel(
    const int cur_label, LVertex* array, const int size);
bool AddUniqueLVertex(
    const LVertex lv, LVertex* array, int* size);
bool RemoveLVertex(
    const LVertex lv, LVertex* array, int* size);
LVertex PopLVertex(LVertex* array, int* size, const int label);
int FindInt(const int val, int* array, const int size);
bool AddInt(const int val, int* array, int* size);
bool AddMEdge(
    const int label0, const int label1,
    bool* matrix, int* labels, int* num_labels);
int3 NeighborPosition(
    __GLOBAL__ const Vertex* v, const int3 v_p, const Direction d);
int IsAmbiguous3_orig(
    const int base_vi,
    __GLOBAL__ const Vertex* vertices,
    __GLOBAL__ const int* point_labels);

//--------------------------------------------------
//--------------------------------------------------
// Definitions
//--------------------------------------------------
//--------------------------------------------------

bool Inside(const int3 p, const int width) {
  return (p.x >= 0 && p.y >= 0 && p.z >= 0 &&
          p.x <= width && p.y <= width && p.z <= width);
}

int FindLVertex(const LVertex lv, LVertex* array, const int size) {
  for (int i = 0; i < size; ++i) {
    if (array[i].idx == lv.idx)
      return i;
  }
  return -1;
}

int FindLVertexWithLabel(
    const int cur_label, LVertex* array, const int size) {
  for (int i = 0; i < size; ++i) {
    if (array[i].label == cur_label)
      return i;
  }
  return -1;
}

// Only adds if the vertex doesn't already exist in the array.
// Given lv doesn't exist in the array, it is added if lv.label == cur_label
// or if no vertex with label lv.label is already in the array.
bool AddUniqueLVertex(
    const LVertex lv, LVertex* array, int* size) {
  if (FindLVertex(lv, array, *size) > -1) return true;
  if (*size == MAX_SIZE) return false;
  array[(*size)++] = lv;
  return true;
}

// Returns true if the vertex existed and was removed.  Returns false
// if the vertex was not in the array.
bool RemoveLVertex(
    const LVertex lv, LVertex* array, int* size) {
  const int i = FindLVertex(lv, array, *size);
  if (i == -1)
    return false;

  LVertex temp = array[i];
  array[i] = array[(*size)-1];
  array[(*size)-1] = temp;
  --(*size);
  return true;
}

// Takes an element out of the array and returns it.  If there are
// any elements with label label, then an element with label
// is guaranteed to be returned.  If no such element exists, an arbitrary
// element is returned.
LVertex PopLVertex(LVertex* array, int* size, const int label) {
  if (array[(*size)-1].label != label) {
    // Find the first element with label "label"
    for (int i = 0; i < *size; ++i) {
      if (array[i].label == label) {
        // Swap element at idx with last element
        LVertex temp = array[i];
        array[i] = array[(*size)-1];
        array[(*size)-1] = temp;
        break;
      }
    }
  }
  // Return last element
  return array[--(*size)];
}

int FindInt(const int val, int* array, const int size) {
  for (int i = 0; i < size; ++i) {
    if (array[i] == val)
      return i;
  }
  return -1;
}

bool AddInt(const int val, int* array, int* size) {
  if (*size == MAX_SIZE) return false;
  array[(*size)++] = val;
  return true;
}

// Adds an edge between labels.
bool AddMEdge(
    const int label0, const int label1,
    bool* matrix, int* labels, int* num_labels) {
  if (label0 != label1) {
    int label0idx=-1, label1idx=-1;
    for (int i = 0; i < *num_labels; ++i) {
      if (labels[i] == label0)
        label0idx = i;
      if (labels[i] == label1)
        label1idx = i;
    }
    if (label0idx == -1) {
      if (*num_labels == MAX_LABELS)
        return false;
      label0idx = (*num_labels)++;
      labels[label0idx] = label0;
    }
    if (label1idx == -1) {
      if (*num_labels == MAX_LABELS)
        return false;
      label1idx = (*num_labels)++;
      labels[label1idx] = label1;
    }
    matrix[label0idx*MAX_LABELS+label1idx] = true;
    matrix[label1idx*MAX_LABELS+label0idx] = true;
  }

  return true;
}

int3 NeighborPosition(
    __GLOBAL__ const Vertex* v, const int3 v_p, const Direction d) {
  const char level = v_neighbor_level(d, *v);
  const int width = Level2CellWidth(level);
  int3 p = v_p;
  if (d.pos & 1)
    // p = p + (int3)(width, 0, 0);
    p = p + make_int3(width, 0, 0);
  if (d.pos & 2)
    p = p + make_int3(0, width, 0);
  if (d.pos & 4)
    p = p + make_int3(0, 0, width);
  if (d.neg) {
    if (d.neg & 1)
      p = p + make_int3(-width, 0, 0);
    if (d.neg & 2)
      p = p + make_int3(0, -width, 0);
    if (d.neg & 4)
      p = p + make_int3(0, 0, -width);
  }
  return p;
}

// Returns 1 for true, 0 for false, -1 for failure
int IsAmbiguous3_orig(
    const int base_vi,
    __GLOBAL__ const Vertex* vertices,
    __GLOBAL__ const int* point_labels) {

  // Initialization
  const unsigned char D = 3;
  __GLOBAL__ const Vertex* base_v = &vertices[base_vi];
  if (!v_is_base(*base_v)) return 0;
  const int width = Level2CellWidth(v_cell_level(*base_v));

  // Vertices that need to be visited
  LVertex queue[MAX_SIZE];
  int queue_size = 0;

  // Vertices that need to be visited but aren't added to the queue because
  // a vertex with the same label already exists in the queue.
  LVertex backlog[MAX_SIZE];
  int backlog_size = 0;

  // All vertices that have been visited already
  int visited[MAX_SIZE];
  int visited_size = 0;

  // All known labels
  bool edge_matrix[MAX_LABELS * MAX_LABELS];
  for (int i = 0; i < MAX_LABELS * MAX_LABELS; ++i)
    edge_matrix[i] = false;
  int labels[MAX_LABELS];
  int labels_size = 0;

  // Labels that have been completed
  int completed_labels[MAX_SIZE];
  int completed_labels_size = 0;

  /* const bool debug = true; */
  const bool debug = false;

  // Visit all vertices incident to the cell
  LVertex base_lv =
      { base_vi, make_int3(0, 0, 0), point_labels[v_closest_point(*base_v)] };
  int cur_label = base_lv.label;
  AddUniqueLVertex(base_lv, queue, &queue_size);
  int ambiguous = 0;
  while (ambiguous == 0 && queue_size > 0) {
    // debug
    if (debug) {
      printf("queue: ");
      for (int k = 0; k < queue_size; ++k) {
        printf("%d ", queue[k].idx);
      }
      printf("\n");
      printf("completed labels: ");
      for (int k = 0; k < completed_labels_size; ++k) {
        printf("%d ", completed_labels[k]);
      }
      printf("\n");
    }

    // Pop from queue
    const LVertex lv = PopLVertex(queue, &queue_size, cur_label);
    RemoveLVertex(lv, backlog, &backlog_size);
    const int vi = lv.idx;
    if (lv.label != cur_label) {
      // Add cur_label to completed completed_labels array
      AddInt(cur_label, completed_labels, &completed_labels_size);
      cur_label = lv.label;
    }
    if (debug) {
      printf("vertex = %d, label = %d\n", vi, cur_label);
    }

    __GLOBAL__ const Vertex* v = &vertices[vi];

    // Update visited array
    if (!AddInt(vi, visited, &visited_size)) {
      if (debug) {
        printf("failed to add visited %d: visited_size = %d\n",
               vi, visited_size);
      }
      ambiguous = -1;
    }

    // Visit neighbors
    for (int axis = 0; ambiguous == 0 && axis < D; ++axis) {
      for (int pos = 0; pos < 2; ++pos) {
        const Direction dir = DirectionFromAxis(axis, pos == 1);
        const int n_vi = v_neighbor_index(dir, *v);
        if (n_vi > -1) {
          // A neighbor exists
          const int3 n_p = NeighborPosition(v, lv.p, dir);
          __GLOBAL__ const Vertex* n_v = &vertices[n_vi];
          const LVertex n_lv =
              { n_vi, n_p, point_labels[v_closest_point(*n_v)] };
          const int n_label = n_lv.label;

          if (Inside(n_p, width)) {
            if (!AddMEdge(lv.label, n_label, edge_matrix,
                         labels, &labels_size)) {
              printf("Failed to add edge");
              ambiguous = -1;
            }
          }

          if (Inside(n_p, width) &&
              FindInt(n_vi, visited, visited_size) == -1) {
            // The neighbor is incident to the cell and hasn't already
            // been visited
            if (n_label != cur_label &&
                FindLVertexWithLabel(n_label, queue, queue_size) > -1) {
              AddUniqueLVertex(n_lv, backlog, &backlog_size);
            } else {
              // Either the neighbor has the current label or no vertex with
              // the neighbor's label has already been added to the queue
              if (FindInt(
                      n_label, completed_labels, completed_labels_size) > -1) {
                // We have a neighbor that has a label that has already been
                // completed, so label n_label is disconnected and the cell
                // is ambiguous.
                ambiguous = 1;
              }

              if (!AddUniqueLVertex(n_lv, queue, &queue_size)) {
                printf("failed to add unique vertex %d: queue_size = %d\n",
                       n_vi, queue_size);
                ambiguous = -1;
              }
            }
          }
        }
      }
    }
  }

  // The edge matrix must be all 1s in order for the collapsed
  // cell to be a clique.
  for (int i = 0; ambiguous == 0 && i < labels_size; ++i) {
    for (int j = 0; ambiguous == 0 && j < labels_size; ++j) {
      if (i != j && !edge_matrix[i*MAX_LABELS+j])
        ambiguous = 1;
    }
  }

  if (ambiguous == 0 && backlog_size > 0) {
    // This is an additional disconnected case
    ambiguous = 1;
  }

  return ambiguous;
}

// Returns 1 for true, 0 for false, -1 for failure
int IsAmbiguous3(
    const int base_vi,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    __GLOBAL__ GeomPoint* cpoints_array) {
  UVertexNetwork vn =
      make_uvertex_network(vn_header, vertex_array, cpoints_array);

  if (!IsBase(base_vi, vn))
    return 0;

  // Initialization
  const unsigned char D = 3;
  __GLOBAL__ const Vertex* base_v = &vn.vertices[base_vi];
  if (!v_is_base(*base_v)) return 0;
  const int width = Level2CellWidth(v_cell_level(*base_v));

  // Vertices that need to be visited
  LVertex queue[MAX_SIZE];
  int queue_size = 0;

  // Vertices that need to be visited but aren't added to the queue because
  // a vertex with the same label already exists in the queue.
  LVertex backlog[MAX_SIZE];
  int backlog_size = 0;

  // All vertices that have been visited already
  int visited[MAX_SIZE];
  int visited_size = 0;

  // All known labels
  bool edge_matrix[MAX_LABELS * MAX_LABELS];
  for (int i = 0; i < MAX_LABELS * MAX_LABELS; ++i)
    edge_matrix[i] = false;
  int labels[MAX_LABELS];
  int labels_size = 0;

  // Labels that have been completed
  int completed_labels[MAX_SIZE];
  int completed_labels_size = 0;

  // const bool debug = true;
  const bool debug = false;
  // const bool debug = (base_vi == 0);

  // Visit all vertices incident to the cell
  LVertex base_lv =
      { base_vi, make_int3(0, 0, 0), Label(base_vi, vn) };
  int cur_label = base_lv.label;
  AddUniqueLVertex(base_lv, queue, &queue_size);
  int ambiguous = 0;
  while (ambiguous == 0 && queue_size > 0) {
    // debug
    if (debug) {
      printf("queue: ");
      for (int k = 0; k < queue_size; ++k) {
        printf("%d ", queue[k].idx);
      }
      printf("\n");
      printf("completed labels: ");
      for (int k = 0; k < completed_labels_size; ++k) {
        printf("%d ", completed_labels[k]);
      }
      printf("\n");
    }

    // Pop from queue
    const LVertex lv = PopLVertex(queue, &queue_size, cur_label);
    RemoveLVertex(lv, backlog, &backlog_size);
    const int vi = lv.idx;
    if (lv.label != cur_label) {
      // Add cur_label to completed completed_labels array
      AddInt(cur_label, completed_labels, &completed_labels_size);
      cur_label = lv.label;
    }
    if (debug) {
      printf("    vertex = %d, label = %d\n", vi, cur_label);
    }

    __GLOBAL__ const Vertex* v = &vn.vertices[vi];

    // Update visited array
    if (!AddInt(vi, visited, &visited_size)) {
      if (debug) {
        printf("failed to add visited %d: visited_size = %d\n",
               vi, visited_size);
      }
      ambiguous = -1;
    }

    const int voi = 854;

    // Visit neighbors
    for (int axis = 0; ambiguous == 0 && axis < D; ++axis) {
      for (int pos = 0; pos < 2; ++pos) {
        const Direction dir = DirectionFromAxis(axis, pos == 1);
        const int n_vi = v_neighbor_index(dir, *v);
        if (n_vi > -1) {
          // if (debug && n_vi == voi) printf("2 %d\n", voi);
          // A neighbor exists
          const int3 n_p = NeighborPosition(v, lv.p, dir);
          // __GLOBAL__ const Vertex* n_v = &vn.vertices[n_vi];
          // const LVertex n_lv =
          //     { n_vi, n_p, point_labels[v_closest_point(*n_v)] };
          const LVertex n_lv =
              { n_vi, n_p, Label(n_vi, vn) };
          const int n_label = n_lv.label;
          if (debug && n_vi == voi) printf("n_vi=%d label=%d\n",
                                           n_vi, (int)n_label);

          if (Inside(n_p, width)) {
            // if (debug && n_vi == voi) printf("3 %d\n", voi);
            if (!AddMEdge(lv.label, n_label, edge_matrix,
                         labels, &labels_size)) {
              printf("Failed to add edge");
              ambiguous = -1;
            }
          }

          if (Inside(n_p, width) &&
              FindInt(n_vi, visited, visited_size) == -1) {
            // if (debug && n_vi == voi) printf("4 %d\n", voi);
            // The neighbor is incident to the cell and hasn't already
            // been visited
            if (n_label != cur_label &&
                FindLVertexWithLabel(n_label, queue, queue_size) > -1) {
              // if (debug && n_vi == voi)
              //   printf("5 %d %d %d\n", n_vi, cur_label, n_label);
              // if (FindLVertex(n_lv, backlog, backlog_size) > -1) {
              //   if (debug && n_vi == voi) printf("9 %d\n", voi);
              // }
              if (!AddUniqueLVertex(n_lv, backlog, &backlog_size)) {
                // if (debug && n_vi == voi) printf("8 %d\n", voi);
              }
            } else {
              // if (debug && n_vi == voi) printf("6 %d\n", voi);
              // Either the neighbor has the current label or no vertex with
              // the neighbor's label has already been added to the queue
              if (FindInt(
                      n_label, completed_labels, completed_labels_size) > -1) {
                // if (debug && n_vi == voi) printf("7 %d\n", voi);
                // We have a neighbor that has a label that has already been
                // completed, so label n_label is disconnected and the cell
                // is ambiguous.
                ambiguous = 1;
              }

              if (!AddUniqueLVertex(n_lv, queue, &queue_size)) {
                printf("failed to add unique vertex %d: queue_size = %d\n",
                       n_vi, queue_size);
                ambiguous = -1;
              }
            }
          }
        }
      }
    }
  }

  // The edge matrix must be all 1s in order for the collapsed
  // cell to be a clique.
  for (int i = 0; ambiguous == 0 && i < labels_size; ++i) {
    for (int j = 0; ambiguous == 0 && j < labels_size; ++j) {
      if (i != j && !edge_matrix[i*MAX_LABELS+j])
        ambiguous = 1;
    }
  }

  if (ambiguous == 0 && backlog_size > 0) {
    // This is an additional disconnected case
    ambiguous = 1;
  }

  return ambiguous;
}

NAMESPACE_OCT_END
