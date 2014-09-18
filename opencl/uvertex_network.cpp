#include "./uvertex_network.h"
#include "./bit.h"

NAMESPACE_OCT_BEGIN

//------------------------------------------------------------
// Utility declarations
//------------------------------------------------------------

int FindInIntArray(int size, int* array, int value);
int FirstBaseNeighbor(const int vi, Direction d, UVertexNetwork vn);

//------------------------------------------------------------
// Private modifier declarations
//------------------------------------------------------------

// Returns the index of the vertex
int CreateVertex(const intn position, UVertexNetwork* vn);

void AddEdge(const int a_vi, const int b_vi, const Direction d,
             const level_t level, UVertexNetwork* vn);

int Break(const int a_vi, const int b_vi, const Direction d,
          const level_t level, UVertexNetwork* vn);

void SubdivideImpl1(
    int slvi2vi[], const int offset, const int axes[], const int level,
    UVertexNetwork* vn);

void SubdivideImpl2(
    int slvi2vi[], const int offset, const int axes[], const int level,
    UVertexNetwork* vn);

void SubdivideEdges2(
    int slvi2vi[], const int offset, const int axes[], const int level,
    UVertexNetwork* vn);

void SubdivideFaces2(
    int slvi2vi[], const int offset, const int axes[], const int level,
    UVertexNetwork* vn);

void SubdivideImpl3(
    int slvi2vi[], const int offset, const int axes[], const int level,
    UVertexNetwork* vn);

void SubdivideEdges3(const int vi,
    int slvi2vi[], const int offset, const int axes[], const int level,
    UVertexNetwork* vn);

void SubdivideFaces3(
    int slvi2vi[], const int offset, const int axes[], const int level,
    UVertexNetwork* vn);

void SubdivideVolumes3(
    int slvi2vi[], const int offset, const int axes[], const int level,
    UVertexNetwork* vn);

void SubdivideImplAddEdgesFromCenter(
    int slvi2vi[], const int offset, const int axes[], const int level,
    const int dim, UVertexNetwork* vn);

//------------------------------------------------------------
// Definitions
//------------------------------------------------------------

// Indexed by negative directions from vertex towards cell center,
// e.g. 01 = positive y, negative x
//
//     01 \_         _/ 00
//          \_     _/
//            \_ _/
//             _*_
//           _/   \_
//         _/       \_
//     11 /           \ 10
//
//
//      2-----------------3
//      |                 |
//      |                 |
//      |                 | 0: 0 x x x
//      |                 | 1: x 0 x x
//      |                 | 2: x x 0 x
//      |                 | 3: x x x 0
//      |                 |
//      0-----------------1
//
//      2--------5--------3 0: 0 x x x
//      |        |        | 1: x 4 x x
//      |        |        | 2: x x 6 x
//      |        |        | 3: x x x 8
//      6--------8--------7 4: 4 0 x x
//      |        |        | 5: x x 8 6
//      |        |        | 6: 6 x 0 x
//      |        |        | 7: x 8 x 4
//      0--------4--------1 8: 8 6 4 0
//

void Initialize(UVertexNetwork* vn) {
  const int base_vi = CreateVertex(make_intn(0), vn);
  assert(base_vi == 0);
  __GLOBAL__ int* lvi2vi = vn->vertices[base_vi].corners;
  lvi2vi[0] = 0;
  for (int i = 1; i < (1<<DIM); ++i) {
    intn p = make_intn(0);
    for (int j = 0; j < DIM; ++j) {
      if ((i & (1<<j))) {
        // p.s[j] = Level2CellWidth(0);
        set_intn_comp(j, &p, Level2CellWidth(0));
      }
    }
    lvi2vi[i] = CreateVertex(p, vn);
  }

  // Add edges
  for (int lvi = 0; lvi < (1<<DIM); ++lvi) {
    const int vi = lvi2vi[lvi];
    for (int axis = 0; axis < DIM; ++axis) {
      const int n_lvi = NbrVertexPositive(lvi, axis);
      if (n_lvi != -1) {
        const int n_vi = lvi2vi[n_lvi];
        const Direction d = DirectionFromAxis(axis, true);
        AddEdge(vi, n_vi, d, 0, vn);
      }
    }
  }

  v_set_is_base(true, &vn->vertices[0]);
  v_set_cell_level(0, &vn->vertices[0]);
}

bool SetClosestPoint(const int vi, const int candidate, UVertexNetwork* vn) {
  bool changed = false;
  const int cur = v_closest_point(vn->vertices[vi]);
  if (cur == -1) {
    // No closest point
    v_set_closest_point(candidate, &vn->vertices[vi]);
    changed = true;
  } else {
    const intn vp = vn->vertices[vi].position;
    const float cur_dist = length(convert_floatn(vn->cpoints[cur].p-vp));
    const float new_dist = length(convert_floatn(vn->cpoints[candidate].p-vp));
    if (new_dist < cur_dist) {
      v_set_closest_point(candidate, &vn->vertices[vi]);
      changed = true;
    }
  }
  return changed;
}

int FindNeighbor(const int a_vi, const Direction d,
                 const level_t level, UVertexNetwork vn) {
  const int goal_dist = Level2CellWidth(level);
  int cur_dist = 0;
  int cur_vi = a_vi;
  while (cur_dist < goal_dist) {
    if (v_neighbor_level(d, vn.vertices[cur_vi]) == -1) return -1;
    cur_dist += Level2CellWidth(v_neighbor_level(d, vn.vertices[cur_vi]));
    cur_vi = v_neighbor_index(d, vn.vertices[cur_vi]);
  }
  if (cur_dist > goal_dist) return -1;
  return cur_vi;
}

int CreateVertex(const intn position, UVertexNetwork* vn) {
  // const int vi = vn->_num_vertices++;
  const int vi = IncNumVertices(vn);
  vn->vertices[vi] = make_vertex(position);
  return vi;
}

void Subdivide(
    const int base_vi, int* slvi2vi_, UVertexNetwork* vn) {
  assert(base_vi < NumVertices(*vn));

  const int level = CellLevel(base_vi, *vn);
  const __GLOBAL__ int* lvi2vi = vn->vertices[base_vi].corners;

  // Set corners of slvi2vi
  int slvi2vi[kNumSubdivided];
  for (int i = 0; i < kNumSubdivided; ++i) {
    slvi2vi[i] = -1; // initialize to -1
  }
  // BuildSlvi2Vi(base_vi, slvi2vi, CellLevel(base_vi, *vn), *vn);
  for (int lvi = 0; lvi < 1<<DIM; ++lvi) {
    // assert(slvi2vi[lvi2slvi[lvi]] == lvi2vi[lvi]);
    slvi2vi[lvi2slvi[lvi]] = lvi2vi[lvi];
  }
    
  int axes[DIM];
  for (int i = 0; i < DIM; ++i) {
    axes[i] = i;
  }

  if (DIM == 3) {
    // SubdivideImpl3(slvi2vi, 0, axes, level+1, vn);
    SubdivideEdges3(slvi2vi[0], slvi2vi, 0, axes, level+1, vn);
    SubdivideFaces3(slvi2vi, 0, axes, level+1, vn);
    SubdivideVolumes3(slvi2vi, 0, axes, level+1, vn);
  }
  else {
    // SubdivideImpl2(slvi2vi, 0, axes, level+1, vn);
    SubdivideEdges2(slvi2vi, 0, axes, level+1, vn);
    SubdivideFaces2(slvi2vi, 0, axes, level+1, vn);
  }

  // Set base and level of subcells
  int subcells[1<<DIM];
  for (int lvi = 0; lvi < (1<<DIM); ++lvi) {
    // Set sub_lvi2vi
    const Direction offset = DirectionFromPosNeg(0, lvi^((1<<DIM)-1));
    const int sub_vi = slvi2vi[Slvi(&offset)];
    subcells[lvi] = sub_vi;

    for (int i = 0; i < (1<<DIM); ++i) {
      Direction added = Add(&offset, i, 0);
      vn->vertices[sub_vi].corners[i] = slvi2vi[Slvi(&added)];
    }

    v_set_is_base(true, &vn->vertices[sub_vi]);
    v_set_cell_level(level+1, &vn->vertices[sub_vi]);
  }

  // Update slvi2vi if desired (if non-null)
  if (slvi2vi_) {
    for (int i = 0; i < kNumSubdivided; ++i) {
      slvi2vi_[i] = slvi2vi[i];
    }
  }

  for (int i = 0; i < kNumSubdivided; ++i) {
    assert(slvi2vi[i] == Slvi2Vi(base_vi, i, CellLevel(base_vi, *vn)-1, *vn));
  }
}

// void AddCenter() {
//     c_vi = CreateVertex(c_p, vn);
//     v_set_neighbor(c_vi, d, level, &vn->vertices[a_vi]);
//     v_set_neighbor(b_vi, d, level, &vn->vertices[c_vi]);
//     v_set_neighbor(a_vi, Reversed(&d), level, &vn->vertices[c_vi]);
//     v_set_neighbor(c_vi, Reversed(&d), level, &vn->vertices[b_vi]);
// }

//         ---------v3---------
//        |                    |
//        |                    |
//        |                    |
//        |                    |
//       v0         *          v1
//        |                    |
//    ^   |                    |
//    |   |                    |
// d1 |    ---------v2---------
//    |
//    --------->
//       d0

inline intn FindPosition(const intn base_p, const int axis, const int width) {
  intn p = base_p;
  if (axis == 0)
    p.x += width;
  else if (axis == 1)
    p.y += width;
#ifdef OCT3D
  else if (axis == 2)
    p.z += width;
#endif
  return p;
}

// level is the current level of the cell
// cp is the location of the center
// Returns the center vertex
int AddCenter2(
    const int v0, const int v1, const int axis0,
    const int v2, const int v3, const int axis1,
    const int level, UVertexNetwork vn) {
  // w is half the distance from v0 to v1
  const int w = Level2CellWidth(level+1);
  const Direction d0 = DirectionFromPosAxis(axis0);
  const Direction d1 = DirectionFromPosAxis(axis1);
  int c = -1;
  if (Neighbor(v0, d0, vn) == -1) {
    const intn cp = FindPosition(vn.vertices[v0].position, axis0, w);
    c = CreateVertex(cp, &vn);
    AddEdge(v0, c, d0, level+1, &vn);
    AddEdge(c, v1, d0, level+1, &vn);
    AddEdge(v2, c, d1, level+1, &vn);
    AddEdge(c, v3, d1, level+1, &vn);
  } else {
    c = FindNeighbor(v0, d0, level+1, vn);
  }
  return c;
}

// Center vertices are labeled by their normal axis
int AddCenter3(
    const int vx0, const int vx1, 
    const int vy0, const int vy1, 
    const int vz0, const int vz1, 
    const int level, UVertexNetwork vn) {
  // w is half the distance from v0 to v1
  const int w = Level2CellWidth(level+1);
  const Direction dx = DirectionFromPosAxis(0);
  const Direction dy = DirectionFromPosAxis(1);
  const Direction dz = DirectionFromPosAxis(2);
  const intn cp = FindPosition(vn.vertices[vx0].position, 0, w);
  const int c = CreateVertex(cp, &vn);
  AddEdge(vx0, c, dx, level+1, &vn);
  AddEdge(c, vx1, dx, level+1, &vn);
  AddEdge(vy0, c, dy, level+1, &vn);
  AddEdge(c, vy1, dy, level+1, &vn);
  AddEdge(vz0, c, dz, level+1, &vn);
  AddEdge(c, vz1, dz, level+1, &vn);
  return c;
}

void SetCell(__GLOBAL__ Vertex* v, int level,
             int c0, int c1, int c2, int c3, int c4, int c5, int c6, int c7) {
  v_set_is_base(true, v);
  v_set_cell_level(level, v);
  v->corners[0] = c0;
  v->corners[1] = c1;
  v->corners[2] = c2;
  v->corners[3] = c3;
#ifdef OCT3D
  v->corners[4] = c4;
  v->corners[5] = c5;
  v->corners[6] = c6;
  v->corners[7] = c7;
#endif
}

void Subdivide_A(const int base_vi, UVertexNetwork vn) {
  assert(base_vi < NumVertices(vn));

  const int level = CellLevel(base_vi, vn);

  const Direction x = DirectionFromPosAxis(0);
  const Direction y = DirectionFromPosAxis(1);
  const Direction z = DirectionFromPosAxis(2);
  __GLOBAL__ const int* corners = vn.vertices[base_vi].corners;

  const int b0 = corners[0];
  const int b1 = corners[1];
  const int b2 = corners[2];
  const int b3 = corners[3];
  const int b4 = corners[4];
  const int b5 = corners[5];
  const int b6 = corners[6];
  const int b7 = corners[7];

  const int b01 = Break(b0, b1, x, level+1, &vn);
  const int b23 = Break(b2, b3, x, level+1, &vn);
  const int b45 = Break(b4, b5, x, level+1, &vn);
  const int b67 = Break(b6, b7, x, level+1, &vn);

  const int b02 = Break(b0, b2, y, level+1, &vn);
  const int b13 = Break(b1, b3, y, level+1, &vn);
  const int b46 = Break(b4, b6, y, level+1, &vn);
  const int b57 = Break(b5, b7, y, level+1, &vn);

  const int b04 = Break(b0, b4, z, level+1, &vn);
  const int b15 = Break(b1, b5, z, level+1, &vn);
  const int b26 = Break(b2, b6, z, level+1, &vn);
  const int b37 = Break(b3, b7, z, level+1, &vn);

  const int w2 = Level2CellWidth(level+1);

  const int xaxis = 0;
  const int yaxis = 1;
  const int zaxis = 2;

  const int b0145 = AddCenter2(b04, b15, xaxis, b01, b45, zaxis, level, vn);
  const int b2367 = AddCenter2(b26, b37, xaxis, b23, b67, zaxis, level, vn);
  const int b0246 = AddCenter2(b04, b26, yaxis, b02, b46, zaxis, level, vn);
  const int b1357 = AddCenter2(b15, b37, yaxis, b13, b57, zaxis, level, vn);
  const int b0123 = AddCenter2(b01, b23, yaxis, b02, b13, xaxis, level, vn);
  const int b4567 = AddCenter2(b45, b67, yaxis, b46, b57, xaxis, level, vn);

  const int center =
      AddCenter3(b0246, b1357, b0145, b2367, b0123, b4567, level, vn);

  SetCell(&vn.vertices[b0], level+1,
          b0, b01, b02, b0123, b04, b0145, b0246, center);
  SetCell(&vn.vertices[b01], level+1,
          b01, b1, b0123, b13, b0145, b15, center, b1357);
  SetCell(&vn.vertices[b04], level+1,
          b04, b0145, b0246, center, b4, b45, b46, b4567);
  SetCell(&vn.vertices[b0145], level+1,
          b0145, b15, center, b1357, b45, b5, b4567, b57);
  SetCell(&vn.vertices[b02], level+1,
          b02, b0123, b2, b23, b0246, center, b26, b2367);
  SetCell(&vn.vertices[b0123], level+1,
          b0123, b13, b23, b3, center, b1357, b2367, b37);
  SetCell(&vn.vertices[b0246], level+1,
          b0246, center, b26, b2367, b46, b4567, b6, b67);
  SetCell(&vn.vertices[center], level+1,
          center, b1357, b2367, b37, b4567, b57, b67, b7);
}

void Subdivide_B(const int base_vi, UVertexNetwork vn) {
  assert(base_vi < NumVertices(vn));

  const int level = CellLevel(base_vi, vn);

  const __GLOBAL__ int* lvi2vi = vn.vertices[base_vi].corners;

  // Set corners of slvi2vi
  int slvi2vi[kNumSubdivided];
  // for (int i = 0; i < kNumSubdivided; ++i) {
  //   slvi2vi[i] = -1; // initialize to -1
  // }
  BuildSlvi2Vi(base_vi, slvi2vi, CellLevel(base_vi, vn), vn);
  for (int lvi = 0; lvi < 1<<DIM; ++lvi) {
    assert(slvi2vi[lvi2slvi[lvi]] == lvi2vi[lvi]);
    // slvi2vi[lvi2slvi[lvi]] = lvi2vi[lvi];
  }
    
  int axes[DIM];
  for (int i = 0; i < DIM; ++i) {
    axes[i] = i;
  }

  if (DIM == 3) {
    // SubdivideEdges3(slvi2vi, 0, axes, level+1, vn);
    SubdivideFaces3(slvi2vi, 0, axes, level+1, &vn);
    // SubdivideVolumes3(slvi2vi, 0, axes, level+1, vn);
  }
  // else {
  //   // SubdivideEdges2(slvi2vi, 0, axes, level+1, vn);
  //   SubdivideFaces2(slvi2vi, 0, axes, level+1, &vn);
  // }
  BuildSlvi2Vi(base_vi, slvi2vi, CellLevel(base_vi, vn), vn);
}

void Subdivide_C(const int base_vi, UVertexNetwork vn) {
  assert(base_vi < NumVertices(vn));

  const int level = CellLevel(base_vi, vn);
  const __GLOBAL__ int* lvi2vi = vn.vertices[base_vi].corners;

  // Set corners of slvi2vi
  int slvi2vi[kNumSubdivided];
  // for (int i = 0; i < kNumSubdivided; ++i) {
  //   slvi2vi[i] = -1; // initialize to -1
  // }
  BuildSlvi2Vi(base_vi, slvi2vi, CellLevel(base_vi, vn), vn);
  for (int lvi = 0; lvi < 1<<DIM; ++lvi) {
    assert(slvi2vi[lvi2slvi[lvi]] == lvi2vi[lvi]);
    // slvi2vi[lvi2slvi[lvi]] = lvi2vi[lvi];
  }
  // if (base_vi == 0) {
  //   // printf("aa value = %d\n", vn.vertices[0].corners[6]);
  //   printf("aa value = %d\n", slvi2vi[12]);
  // }
    
  int axes[DIM];
  for (int i = 0; i < DIM; ++i) {
    axes[i] = i;
  }

  if (DIM == 3) {
    // SubdivideEdges3(slvi2vi, 0, axes, level+1, vn);
    // SubdivideFaces3(slvi2vi, 0, axes, level+1, vn);
    SubdivideVolumes3(slvi2vi, 0, axes, level+1, &vn);
  }
  else {
    assert(0);
  }

  // Set base and level of subcells
  int subcells[1<<DIM];
  for (int lvi = 0; lvi < (1<<DIM); ++lvi) {
    // Set sub_lvi2vi
    const Direction offset = DirectionFromPosNeg(0, lvi^((1<<DIM)-1));
    const int sub_vi = slvi2vi[Slvi(&offset)];
    subcells[lvi] = sub_vi;

    for (int i = 0; i < (1<<DIM); ++i) {
      Direction added = Add(&offset, i, 0);
      const int c_vi = slvi2vi[Slvi(&added)];
      vn.vertices[sub_vi].corners[i] = c_vi;
      // if (c_vi == -1) {
      //   printf("-1 corner: %d\n", base_vi);
      //   return;
      // }
      // if (sub_vi == 0 && i == 6) {
      //   printf("corner: %d, %d, %d\n", sub_vi, i, c_vi);
      // }
    }

    v_set_is_base(true, &vn.vertices[sub_vi]);
    v_set_cell_level(level+1, &vn.vertices[sub_vi]);
  }
  // if (base_vi == 0)
  //   printf("value = %d\n", vn.vertices[0].corners[6]);
}

void SubdivideImpl1(
    int slvi2vi[], const int offset, const int axes[], const int level,
    UVertexNetwork* vn) {
  const int dim = 1;

  //     vi0           vi1           vi2
  //     *-------------*-------------*
  const int axis = axes[0];
  const int stride = kSlviStrides[axis];
  const int vi0 = slvi2vi[offset];
  const int vi2 = slvi2vi[offset+2*stride];
  // const Direction d = DirectionFromAxis(axis, true);
  const Direction d = DirectionFromPosAxis(axis);
  const int vi1 = Break(vi0, vi2, d, level, vn);
  slvi2vi[offset+1*stride] = vi1;
}

//------------------------------------------------------------
// Yes, the following functions could be implemented
// recursively, but this code needs to also compile in OpenCL,
// which doesn't support recursion.
//------------------------------------------------------------

void SubdivideImpl2(
    int slvi2vi[], const int offset, const int axes[], const int level,
    UVertexNetwork* vn) {
  const int dim = 2;

  SubdivideEdges2(slvi2vi, offset, axes, level, vn);
  // // Subdivide dim-1 faces
  // for (int normali = 0; normali < dim; ++normali) {
  //   const int normal = axes[normali];
  //   const int stride = kSlviStrides[normal];
  //   int subaxes[DIM];
  //   for (int i = 0; i < dim-1; ++i) {
  //     subaxes[i] = axes[(normali+i+1)%dim];
  //   }
  //   // From base vertex
  //   SubdivideImpl1(slvi2vi,  offset, subaxes, level, vn);
  //   // From neighbor vertex
  //   SubdivideImpl1(slvi2vi, offset+stride*2, subaxes, level, vn);
  // }

  SubdivideImplAddEdgesFromCenter(slvi2vi, offset, axes, level, dim, vn);
}

void SubdivideEdges2(
    int slvi2vi[], const int offset, const int axes[], const int level,
    UVertexNetwork* vn) {
  const int dim = 2;

  // Subdivide dim-1 faces
  for (int normali = 0; normali < dim; ++normali) {
    const int normal = axes[normali];
    const int stride = kSlviStrides[normal];
    int subaxes[DIM];
    for (int i = 0; i < dim-1; ++i) {
      subaxes[i] = axes[(normali+i+1)%dim];
    }
    // From base vertex
    SubdivideImpl1(slvi2vi,  offset, subaxes, level, vn);
    // From neighbor vertex
    SubdivideImpl1(slvi2vi, offset+stride*2, subaxes, level, vn);
  }
}

void SubdivideFaces2(
    int slvi2vi[], const int offset, const int axes[], const int level,
    UVertexNetwork* vn) {
  const int dim = 2;

  SubdivideImplAddEdgesFromCenter(slvi2vi, offset, axes, level, dim, vn);
}

void SubdivideImpl3(
    int slvi2vi[], const int offset, const int axes[], const int level,
    UVertexNetwork* vn) {
  const int dim = 3;

  SubdivideEdges3(slvi2vi[0], slvi2vi, offset, axes, level, vn);
  SubdivideFaces3(slvi2vi, offset, axes, level, vn);
  // // Subdivide dim-1 faces
  // for (int normali = 0; normali < dim; ++normali) {
  //   const int normal = axes[normali];
  //   const int stride = kSlviStrides[normal];
  //   int subaxes[DIM];
  //   for (int i = 0; i < dim-1; ++i) {
  //     subaxes[i] = axes[(normali+i+1)%dim];
  //   }
  //   // From base vertex
  //   SubdivideImpl2(slvi2vi,  offset, subaxes, level, vn);
  //   // From neighbor vertex
  //   SubdivideImpl2(slvi2vi, offset+stride*2, subaxes, level, vn);
  // }

  SubdivideImplAddEdgesFromCenter(slvi2vi, offset, axes, level, dim, vn);
}

void SubdivideEdges3(const int vi,
    int slvi2vi[], const int offset, const int axes[], const int level,
    UVertexNetwork* vn) {

  const Direction x = DirectionFromPosAxis(0);
  const Direction y = DirectionFromPosAxis(1);
  const Direction z = DirectionFromPosAxis(2);
  // __GLOBAL__ const int* corners = vn->vertices[vi].corners;
  // Break(corners[0], corners[1], x, level, vn);
  // Break(corners[2], corners[3], x, level, vn);
  // Break(corners[4], corners[5], x, level, vn);
  // Break(corners[6], corners[7], x, level, vn);

  // Break(corners[0], corners[2], y, level, vn);
  // Break(corners[1], corners[3], y, level, vn);
  // Break(corners[4], corners[6], y, level, vn);
  // Break(corners[5], corners[7], y, level, vn);

  // Break(corners[0], corners[4], z, level, vn);
  // Break(corners[1], corners[5], z, level, vn);
  // Break(corners[2], corners[6], z, level, vn);
  // Break(corners[3], corners[7], z, level, vn);

  slvi2vi[ 1] = Break(slvi2vi[ 0], slvi2vi[ 2], x, level, vn);
  slvi2vi[ 3] = Break(slvi2vi[ 0], slvi2vi[ 6], y, level, vn);
  slvi2vi[ 9] = Break(slvi2vi[ 0], slvi2vi[18], z, level, vn);
  slvi2vi[ 5] = Break(slvi2vi[ 2], slvi2vi[ 8], y, level, vn);
  slvi2vi[11] = Break(slvi2vi[ 2], slvi2vi[20], z, level, vn);
  slvi2vi[ 7] = Break(slvi2vi[ 6], slvi2vi[ 8], x, level, vn);
  slvi2vi[15] = Break(slvi2vi[ 6], slvi2vi[24], z, level, vn);
  slvi2vi[19] = Break(slvi2vi[18], slvi2vi[20], x, level, vn);
  slvi2vi[21] = Break(slvi2vi[18], slvi2vi[24], y, level, vn);
  slvi2vi[25] = Break(slvi2vi[24], slvi2vi[26], x, level, vn);
  slvi2vi[23] = Break(slvi2vi[20], slvi2vi[26], y, level, vn);
  slvi2vi[17] = Break(slvi2vi[ 8], slvi2vi[26], z, level, vn);

  // if (slvi2vi[0] == 0) {
  //   printf("xxx = %d\n", slvi2vi[3]);
  // }
  // const int dim = 3;

  // // Subdivide dim-1 faces
  // for (int normali = 0; normali < dim; ++normali) {
  //   const int normal = axes[normali];
  //   const int stride = kSlviStrides[normal];
  //   int subaxes[DIM];
  //   for (int i = 0; i < dim-1; ++i) {
  //     subaxes[i] = axes[(normali+i+1)%dim];
  //   }
  //   // From base vertex
  //   SubdivideEdges2(slvi2vi,  offset, subaxes, level, vn);
  //   // From neighbor vertex
  //   SubdivideEdges2(slvi2vi, offset+stride*2, subaxes, level, vn);
  // }
}

void SubdivideFaces3(
    int slvi2vi[], const int offset, const int axes[], const int level,
    UVertexNetwork* vn) {
  const int dim = 3;

  // Subdivide dim-1 faces
  for (int normali = 0; normali < dim; ++normali) {
    const int normal = axes[normali];
    const int stride = kSlviStrides[normal];
    int subaxes[DIM];
    for (int i = 0; i < dim-1; ++i) {
      subaxes[i] = axes[(normali+i+1)%dim];
    }
    // From base vertex
    SubdivideFaces2(slvi2vi,  offset, subaxes, level, vn);
    // From neighbor vertex
    SubdivideFaces2(slvi2vi, offset+stride*2, subaxes, level, vn);
  }
}

void SubdivideVolumes3(
    int slvi2vi[], const int offset, const int axes[], const int level,
    UVertexNetwork* vn) {
  const int dim = 3;
  SubdivideImplAddEdgesFromCenter(slvi2vi, offset, axes, level, dim, vn);
}

void SubdivideImplAddEdgesFromCenter(
    int slvi2vi[], const int offset, const int axes[], const int level,
    const int dim, UVertexNetwork* vn) {
  //----------------------
  // Add edges from center
  //----------------------

  int center_slvi = offset;
  int max_slvi = offset;
  for (int j = 0; j < dim; ++j) {
    center_slvi += kSlviStrides[axes[j]];
    max_slvi += 2*kSlviStrides[axes[j]];
  }
  slvi2vi[center_slvi] = -1;

  // In the 2D case it's possible that the cell has already been subdivided.
  // Check for this by attempting to access the center vertex.
  int vi = slvi2vi[offset];
  for (int axisi = 0; vi != -1 && axisi < dim; ++axisi) {
    const int axis = axes[axisi];
    const Direction d = DirectionFromAxis(axis, 1);
    vi = FindNeighbor(vi, d, level, *vn);
  }
  if (vi != -1) {
    slvi2vi[center_slvi] = vi;
  }

  if (slvi2vi[center_slvi] == -1) {
    // Center hasn't already been added.  Add it and associated edges.
    const intn center_p =
        vn->vertices[slvi2vi[offset]].position/2 +
        vn->vertices[slvi2vi[max_slvi]].position/2;
    const int center_vi = CreateVertex(center_p, vn);
    slvi2vi[center_slvi] = center_vi;
    for (int axisi = 0; axisi < dim; ++axisi) {
      const int axis = axes[axisi];
      const int stride = kSlviStrides[axis];
      for (int pos = 0; pos < 2; ++pos) {
        const Direction d = DirectionFromAxis(axis, pos);
        const int n_vi = slvi2vi[center_slvi + stride*(pos?1:-1)];
        AddEdge(center_vi, n_vi, d, level, vn);
      }
    }
  }
}

void AddEdge(
    const int a_vi, const int b_vi, const Direction d,
    const level_t level, UVertexNetwork* vn) {
  assert(v_neighbor_index(d, vn->vertices[a_vi]) == -1);
  assert(v_neighbor_index(Reversed(&d), vn->vertices[b_vi]) == -1);
  v_set_neighbor(b_vi, d, level, &vn->vertices[a_vi]);
  v_set_neighbor(a_vi, Reversed(&d), level, &vn->vertices[b_vi]);
}

// Cases:
//              d ----->
//            _______
//           |       | <- level
//
//    1)     a---------------b
//           a-------c-------b  creates and returns c
//
//    2)     a---*---c---*---b  returns c
//                              * = 0 or more points
//
int Break(const int a_vi, const int b_vi, const Direction d,
          const level_t level, UVertexNetwork* vn) {
  int c_vi;
  // printf("C: %d\n", Neighbor(0, DirectionFromPosAxis(0), *vn));
  // printf("D: a=%d b=%d n=%d\n", a_vi, b_vi, Neighbor(a_vi, d, *vn));
  if (Neighbor(a_vi, d, *vn) == b_vi) {
    // case 1
    // printf("case 1: %d %d\n", a_vi, b_vi);
    const intn c_p =
        vn->vertices[a_vi].position/2 + vn->vertices[b_vi].position/2;
    c_vi = CreateVertex(c_p, vn);
    v_set_neighbor(c_vi, d, level, &vn->vertices[a_vi]);
    v_set_neighbor(b_vi, d, level, &vn->vertices[c_vi]);
    v_set_neighbor(a_vi, Reversed(&d), level, &vn->vertices[c_vi]);
    v_set_neighbor(c_vi, Reversed(&d), level, &vn->vertices[b_vi]);
  } else {
    // case 2 - follow a in direction d until
    c_vi = FindNeighbor(a_vi, d, level, *vn);
    // printf("case 2: %d %d %d %d\n", a_vi, b_vi, Neighbor(a_vi, d, *vn), c_vi);
  }
  return c_vi;
}

// // Cases:
// //              d ----->
// //            _______
// //           |       | <- level
// //
// //    1)     a               b
// //           a-------c-------b  creates and returns c
// //
// //    2)     a---*---c---*---b  returns c
// //                              * = 0 or more points
// //
// int Connect(const int a_vi, const int b_vi, const Direction d,
//           const level_t level, UVertexNetwork* vn) {
//   int c_vi;
//   // printf("C: %d\n", Neighbor(0, DirectionFromPosAxis(0), *vn));
//   // printf("D: a=%d b=%d n=%d\n", a_vi, b_vi, Neighbor(a_vi, d, *vn));
//   if (Neighbor(a_vi, d, *vn) == b_vi) {
//     // case 1
//     // printf("case 1: %d %d\n", a_vi, b_vi);
//     const intn c_p =
//         vn->vertices[a_vi].position/2 + vn->vertices[b_vi].position/2;
//     c_vi = CreateVertex(c_p, vn);
//     v_set_neighbor(c_vi, d, level, &vn->vertices[a_vi]);
//     v_set_neighbor(b_vi, d, level, &vn->vertices[c_vi]);
//     v_set_neighbor(a_vi, Reversed(&d), level, &vn->vertices[c_vi]);
//     v_set_neighbor(c_vi, Reversed(&d), level, &vn->vertices[b_vi]);
//   } else {
//     // case 2 - follow a in direction d until
//     c_vi = FindNeighbor(a_vi, d, level, *vn);
//     // printf("case 2: %d %d %d %d\n", a_vi, b_vi, Neighbor(a_vi, d, *vn), c_vi);
//   }
//   return c_vi;
// }

bool AreAdjacent(const int vi0, const int vi1, UVertexNetwork vn) {
  assert(v_is_base(vn.vertices[vi0]));
  assert(v_is_base(vn.vertices[vi1]));
  const index_t w0 = CellWidth(vi0, vn);
  const index_t w1 = CellWidth(vi1, vn);
  for (int i = 0; i < DIM; ++i) {
    const int min0 = intn_comp(i, vn.vertices[vi0].position);
    const int min1 = intn_comp(i, vn.vertices[vi1].position);
    if (min0 < min1)
      if (min0 + w0 < min1)
        return false;
    if (min1 < min0)
      if (min1 + w1 < min0)
        return false;
  }
  return true;
}

bool IsIncident(const int vi, const int base_vi, UVertexNetwork vn) {
  assert(v_is_base(vn.vertices[base_vi]));
  const index_t width = CellWidth(base_vi, vn);
  const intn p = vn.vertices[vi].position;
  const intn base_p = vn.vertices[base_vi].position;
  for (int i = 0; i < DIM; ++i) {
    if (intn_comp(i, p) < intn_comp(i, base_p)) return false;
    if (intn_comp(i, p) > intn_comp(i, base_p)+width) return false;
  }
  return true;
}

int FindInIntArray(int size, int* array, int value) {
  for (int i = 0; i < size; ++i) {
    if (array[i] == value)
      return i;
  }
  return -1;
}

int FirstBaseNeighbor(const int vi, Direction d, UVertexNetwork vn) {
  int nbr_vi = vi;
  do {
    nbr_vi = Neighbor(nbr_vi, d, vn);
  } while (nbr_vi > -1 && !IsBase(nbr_vi, vn));
  return nbr_vi;
}

// nbrs must be of size at least kNumIncidentVertices
// Returns the number of neighbors
int GetNeighborCells(
    const int base_vi, __GLOBAL__ int* nbrs, UVertexNetwork vn) {

  int num_nbrs = 0;

  int visited[kNumIncidentVertices];
  visited[0] = base_vi;
  int num_visited = 1;

  int stack[kNumIncidentVertices];
  stack[0] = base_vi;
  int stack_idx = 0;
  while (stack_idx > -1) {
    const int vi = stack[stack_idx--];
    for (int axis = 0; axis < DIM; ++axis) {
      for (int pos = 0; pos < 2; ++pos) {
        const int n_vi = FirstBaseNeighbor(
            vi, DirectionFromAxis(axis, pos), vn);
        if (n_vi > -1 && FindInIntArray(num_visited, visited, n_vi) == -1) {
          const bool incident = IsIncident(n_vi, base_vi, vn);
          if (IsBase(n_vi, vn) && (incident || !pos)) {
            assert(abs(CellLevel(n_vi, vn)-CellLevel(base_vi, vn))
                   <= kGradation);
            assert(num_nbrs < kNumIncidentCells);
            nbrs[num_nbrs++] = n_vi;
          }
          if (incident) {
            visited[num_visited++] = n_vi;
            assert(num_visited < kNumIncidentVertices || stack_idx < 0);
            stack[++stack_idx] = n_vi;
          }
        }
      }
    }
  }
  return num_nbrs;
}

int NumCells(UVertexNetwork vn) {
  int count = 0;
  for (int vi = 0; vi < NumVertices(vn); ++vi) {
    if (v_is_base(vn.vertices[vi])) ++count;
  }
  return count;
}

// Returns -1 if doesn't exist
// level is the level of the larger cell.  Subcells are at level+1.
int Slvi2Vi(const int base_vi, int slvi, const int level, UVertexNetwork vn) {
  int vi = base_vi;
  // const int l = CellLevel(base_vi, vn);
  int jumps[DIM];
  for (int axis = DIM-1; vi > -1 && axis >= 0; --axis) {
    // const int jumps = slvi / kSlviStrides[axis];
    jumps[axis] = slvi / kSlviStrides[axis];
    assert(jumps[axis] <= 2);
    // if (jumps > 0) {
    //   vi = FindNeighbor(vi, DirectionFromAxis(axis, true), l-(jumps-1), vn);
    // }
    slvi = slvi % kSlviStrides[axis];
  }
  // 2 jumps first, then 1 jumps.  That way we're not cutting through the
  // center of a partially-constructed cell.
  for (int axis = DIM-1; vi > -1 && axis >= 0; --axis) {
    const int j = jumps[axis];
    if (j == 2) {
      vi = FindNeighbor(vi, DirectionFromAxis(axis, true), level, vn);
    }
  }
  for (int axis = DIM-1; vi > -1 && axis >= 0; --axis) {
    const int j = jumps[axis];
    if (j == 1) {
      vi = FindNeighbor(vi, DirectionFromAxis(axis, true), level+1, vn);
    }
  }
  return vi;
}

void BuildSlvi2Vi(const int base_vi, int* slvi2vi, const int level,
                  UVertexNetwork vn) {
  for (int i = 0; i < kNumSubdivided; ++i) {
    slvi2vi[i] = Slvi2Vi(base_vi, i, level, vn);
  }
}

NAMESPACE_OCT_END

