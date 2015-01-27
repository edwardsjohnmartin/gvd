#include "./vertex_network.h"
#include "./opencl/bit.h"
#include "./octree.h"

namespace oct {

// static const int kSlviStrides[3] = { 1, 3, 9 };

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

VertexNetwork::VertexNetwork() {
  Initialize();
}

VertexNetwork::VertexNetwork(
    int num_vertices, Vertex* vertices, int num_cpoints, GeomPoint* cpoints) {
  for (int i = 0; i < num_vertices; ++i) {
    _vertices.push_back(vertices[i]);
  }
  for (int i = 0; i < num_cpoints; ++i) {
    _closest_points.push_back(cpoints[i]);
  }
}

VertexNetwork::~VertexNetwork() {
  // _changes.reset(0);
}

void VertexNetwork::Clear() {
  // _changes.reset(0);
  _vertices.clear();

  // Must always be in an initialized state
  Initialize();
}

void VertexNetwork::Initialize() {
  const int base_vi = CreateVertex(make_intn(0));
  assert(base_vi == 0);
  // Can't cache because _vertices may be reallocated
  // int* lvi2vi = _vertices[base_vi].base2corners;
  _vertices[base_vi].corners[0] = 0;
  for (int i = 1; i < (1<<DIM); ++i) {
    intn p = make_intn(0);
    for (int j = 0; j < DIM; ++j) {
      if ((i & (1<<j))) {
        p.s[j] = Level2CellWidth(0);
      }
    }
    // TODO!!!
    // Bug: the commented-out code should work but for some reason
    // _vertices[base_vi].corners[i] remains -1 for the first two corners.
    const int c_vi = CreateVertex(p);
    _vertices[base_vi].corners[i] = c_vi;
    // _vertices[base_vi].corners[i] = CreateVertex(p);
  }

  // Add edges
  for (int lvi = 0; lvi < (1<<DIM); ++lvi) {
    const int vi = _vertices[base_vi].corners[lvi];
    for (int axis = 0; axis < DIM; ++axis) {
      const int n_lvi = NbrVertexPositive(lvi, axis);
      if (n_lvi != -1) {
        const int n_vi = _vertices[base_vi].corners[n_lvi];
        const Direction d = DirectionFromAxis(axis, true);
        AddEdge(vi, n_vi, d, 0);
      }
    }
  }

  SetIsBase(0, true);
  SetCellLevel(0, 0);
}

int VertexNetwork::CreateVertex(const intn& position) {
  // if (_changes)
  //   _changes->VertexAdded();

  _vertices.push_back(make_vertex(position));
  return _vertices.size() - 1;
}

shared_array<int> VertexNetwork::Subdivide(
    const int base_vi, int* slvi2vi_) {
  using namespace std;

  const int level = CellLevel(base_vi);

  const int* lvi2vi = GetCorners(base_vi);

  // Set corners of slvi2vi
  int slvi2vi[kNumSubdivided];
  fill(slvi2vi, slvi2vi+kNumSubdivided, -1);
  for (int lvi = 0; lvi < 1<<DIM; ++lvi) {
    slvi2vi[lvi2slvi[lvi]] = lvi2vi[lvi];
  }
    
  int axes[DIM];
  for (int i = 0; i < DIM; ++i) {
    axes[i] = i;
  }
  Subdivide(slvi2vi, 0, axes, level+1, DIM);

  // Set base and level of subcells
  shared_array<int> subcells(new int[1<<DIM]);
  for (int lvi = 0; lvi < (1<<DIM); ++lvi) {
    // Set sub_lvi2vi
    const Direction offset = DirectionFromPosNeg(0, lvi^((1<<DIM)-1));
    const int sub_vi = slvi2vi[Slvi(&offset)];
    subcells[lvi] = sub_vi;

    for (int i = 0; i < (1<<DIM); ++i) {
      Direction added = Add(&offset, i, 0);
      _vertices[sub_vi].corners[i] = slvi2vi[Slvi(&added)];
    }

    SetIsBase(sub_vi, true);
    SetCellLevel(sub_vi, level+1);
  }

  // if (_changes) {
  //   for (int i = 0; i < kNumSubdivided; ++i) {
  //     _changes->VertexChanged(slvi2vi[i]);
  //   }
  // }

  // Update slvi2vi if desired (if non-null)
  if (slvi2vi_) {
    memcpy(slvi2vi_, slvi2vi, kNumSubdivided*sizeof(int));
  }

  return subcells;
}

void VertexNetwork::Subdivide(
    int slvi2vi[], const int offset, const int axes[], const int level,
    const int dim) {

  if (dim == 0)
    throw logic_error("dim == 0 not expected by Subdivide");

  if (dim == 1) {
    //     vi0           vi1           vi2
    //     *-------------*-------------*
    const int axis = axes[0];
    const int stride = kSlviStrides[axis];
    const int vi0 = slvi2vi[offset];
    const int vi2 = slvi2vi[offset+2*stride];
    const Direction d = DirectionFromAxis(axis, true);
    // const bool needs_break = (Neighbor(vi0, d) == vi2);
    assert(vi0 < size());
    assert(vi2 < size());
    const int vi1 = Break(vi0, vi2, d, level);
    slvi2vi[offset+1*stride] = vi1;
  } else {
    // Subdivide dim-1 faces
    for (int normali = 0; normali < dim; ++normali) {
      const int normal = axes[normali];
      const int stride = kSlviStrides[normal];
      shared_array<int> subaxes(new int[dim-1]);
      for (int i = 0; i < dim-1; ++i) {
        subaxes[i] = axes[(normali+i+1)%dim];
      }
      // From base vertex
      Subdivide(slvi2vi,  offset, subaxes.get(), level, dim-1);
      // From neighbor vertex
      Subdivide(slvi2vi, offset+stride*2, subaxes.get(), level, dim-1);
    }

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
      vi = FindNeighbor(vi, d, level);
    }
    if (vi != -1) {
      slvi2vi[center_slvi] = vi;
    }

    if (slvi2vi[center_slvi] == -1) {
      // Center hasn't already been added.  Add it and associated edges.
      assert(slvi2vi[offset] < size());
      assert(slvi2vi[max_slvi] < size());
      const intn center_p =
          // _positions[slvi2vi[offset]]/2 + _positions[slvi2vi[max_slvi]]/2;
          _vertices[slvi2vi[offset]].position/2 +
          _vertices[slvi2vi[max_slvi]].position/2;
      const int center_vi = CreateVertex(center_p);
      slvi2vi[center_slvi] = center_vi;
      for (int axisi = 0; axisi < dim; ++axisi) {
        const int axis = axes[axisi];
        const int stride = kSlviStrides[axis];
        for (int pos = 0; pos < 2; ++pos) {
          const Direction d = DirectionFromAxis(axis, pos);
          const int n_vi = slvi2vi[center_slvi + stride*(pos?1:-1)];
          AddEdge(center_vi, n_vi, d, level);
        }
      }
    }
  }
}

void VertexNetwork::AddEdge(const int a_vi, const int b_vi, const Direction& d,
             const level_t level) {
  assert(v_neighbor_index(d, _vertices[a_vi]) == -1);
  assert(v_neighbor_index(Reversed(&d), _vertices[b_vi]) == -1);
  v_set_neighbor(b_vi, d, level, &_vertices[a_vi]);
  v_set_neighbor(a_vi, Reversed(&d), level, &_vertices[b_vi]);
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
int VertexNetwork::Break(const int a_vi, const int b_vi, const Direction& d,
          const level_t level) {
  int c_vi;
  if (Neighbor(a_vi, d) == b_vi) {
    // case 1
    // const intn c_p = _positions[a_vi]/2 + _positions[b_vi]/2;
    const intn c_p = _vertices[a_vi].position/2 + _vertices[b_vi].position/2;
    c_vi = CreateVertex(c_p);
    v_set_neighbor(c_vi, d, level, &_vertices[a_vi]);
    v_set_neighbor(b_vi, d, level, &_vertices[c_vi]);
    v_set_neighbor(a_vi, Reversed(&d), level, &_vertices[c_vi]);
    v_set_neighbor(c_vi, Reversed(&d), level, &_vertices[b_vi]);
  } else {
    // case 2 - follow a in direction d until
    c_vi = FindNeighbor(a_vi, d, level);
  }
  return c_vi;
}

}

