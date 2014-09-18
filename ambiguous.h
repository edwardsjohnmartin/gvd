#ifndef __AMBIGUOUS_H__
#define __AMBIGUOUS_H__

#include "./octree.h"
#include "./search.h"
#include "./vertices_gpu_state.h"
#include "./opencl/uvertex_network.h"
#include "./mvertex_network.h"

namespace oct {

// Given an axis quickly find its index.
// name | axis | idx
// -----|------|-----
//  x   |  1   |  0
//  y   |  2   |  1
//  z   |  4   |  2
static __CONST__ int AXIS_IDX[] = { -1, 0, 1, -1, 2 };

// Given a normal, find the first orthogonal axis.
// name | axis | oname1 | oaxis1
// -----|------|--------|-------
//  x   |  1   |   y    |   2
//  y   |  2   |   z    |   4
//  z   |  4   |   x    |   1
static __CONST__ int OAXIS1[] = { -1, 2, 4, -1, 1 };

// Given a normal, find the second orthogonal axis.
// name | axis | oname2 | oaxis2
// -----|------|--------|-------
//  x   |  1   |   z    |   4
//  y   |  2   |   x    |   1
//  z   |  4   |   y    |   2
static __CONST__ int OAXIS2[] = { -1, 4, 1, -1, 2 };

// template <int D>
shared_array<int> GetAmbiguousGpu(
    const UVertexNetwork& vertices, VerticesGpuState<DIM>& gpu_state,
    const OctreeOptions& o);

template <int D>
bool IsAmbiguous(
    const int base_vi, const std::vector<int>& incident,
    const VertexNetwork& vertices, const OctreeOptions& o);

// template <int D>
void SubdivideAmbiguousGpu(
    MVertexNetwork& vertices,
    VerticesGpuState<DIM>& gpu_state,
    const OctreeOptions& o);

void SubdivideAmbiguousCpu(
    VertexNetwork& vertices,
    const OctreeOptions& o);

template <int D>
void SubdivideAmbiguous(
    VertexNetwork& vertices,
    VerticesGpuState<D>& gpu_state,
    const OctreeOptions& o);

//------------------------------------------------------------
//------------------------------------------------------------
// Utility
//------------------------------------------------------------
//------------------------------------------------------------

// Indices into the directions
//
//      ------------------
//     |        |         |
//     |   01   |   00    |
//     |        |         |
//      --------*---------
//     |        |         |
//     |   11   |   10    |
//     |        |         |
//      ------------------
//
template <int N>
class FindBasesVisitor : public VertexVisitor<N> {
 public:
  // N is the native dimension of the space.  If Normal != 0 then
  // we're working in a subspace.  In practice, this occurs when
  // N=3 and D=2.

  typedef Direction Direction_t;

  FindBasesVisitor(const VertexNetwork* vertices, int normal = 0,
                   const FindBasesVisitor<N>* fbv3 = 0)
      : Normal(normal),
        D((Normal==0)?N:N-1),
        M(1<<D),
        MASK(M-1),
        _bases(new int[M*vertices->size()]),
        _distances(new intn[M*vertices->size()]),
        _vertices(vertices), _cardinal_only(true),
        _visited(new bool[vertices->size()]),
        _p2n(new int[1<<D]),
        _n2p(new int[1<<N]),
        _fbv3(fbv3) {
    std::fill(_bases.get(), _bases.get()+M*vertices->size(), -1);
    std::fill(_distances.get(), _distances.get()+M*vertices->size(), make_intn());
    std::fill(_visited.get(), _visited.get()+vertices->size(), false);
    // Set the 0 base of vertex 0
    _bases[0] = 0;

    // For Normal = 001:
    //   _p2n[00] = 000
    //   _p2n[01] = 010
    //   _p2n[10] = 100
    //   _p2n[11] = 110
    // For Normal = 010:
    //   _p2n[00] = 000
    //   _p2n[01] = 001
    //   _p2n[10] = 100
    //   _p2n[11] = 101
    // For Normal = 100:
    //   _p2n[00] = 000
    //   _p2n[01] = 001
    //   _p2n[10] = 010
    //   _p2n[11] = 011

    // Initialize _n2p to -1
    for (int i = 0; i < (1<<N); ++i) {
      _n2p[i] = -1;
    }

    if (Normal == 1) {
      for (int i = 0; i < (1<<D); ++i) {
        _p2n[i] = (i<<1);
        _n2p[_p2n[i]] = i;
      }
    } else if (Normal == 2) {
      for (int i = 0; i < (1<<D); ++i) {
        _p2n[i] = ((i&6)<<1) + (i&1);
        _n2p[_p2n[i]] = i;
      }
    } else {
      for (int i = 0; i < M; ++i) {
        _p2n[i] = i;
        _n2p[_p2n[i]] = i;
      }
    }
  }

  // bool operator()(const int vi, const Vec<index_t, N>& p);
  bool operator()(const int vi, const intn& p);

  inline bool IsBoundary(const int vi, const int axis) const {
    if (D == 2) return false;

    // To be a boundary, vi must have a base in a direction orthogonal
    // to axis, and the cell width of that base must be greater than the
    // distance from vi to the base.
    // for (int oaxis = 1; oaxis < M; oaxis=(oaxis<<1)) {
    //   if (oaxis != axis) {
        // const int base_vi = _bases[vi*M+oaxis];
        // const int dist = _distances[vi*M+oaxis];
        // if (base_vi != -1 && dist < _vertices->LeafCellWidth(base_vi)) {
        //   return true;
        // }
      // }
    // }

    // new method
    const int oaxis = axis;
    const int base_vi = _bases[vi*M+oaxis];
    // const int dist = _distances[vi*M+oaxis].L1norm();
    const int dist = L1norm(_distances[vi*M+oaxis]);
    if (base_vi != -1 && dist == _vertices->CellWidth(base_vi)) {
      return true;
    }

    return false;
  }

  inline bool IsBoundary2(const int vi, const int axis) const {
    if (D == 2) return false;

    const int axis_idx = AXIS_IDX[axis];
    const int axis1_idx = (axis_idx+1)%3;
    const int axis2_idx = (axis_idx+2)%3;
    const int axis1 = (1<<axis1_idx);
    const int axis2 = (1<<axis2_idx);
    const int diag = axis1 | axis2;
    const int dirs[] = { axis1, axis2, diag };
    const int factors[] = { 1, 1, 2 };

    for (int i = 0; i < 3; ++i) {
      const int dir = dirs[i];
      const int base_vi = _bases[vi*M+dir];
      // const int dist = _distances[vi*M+dir].L1norm();
      const int dist = L1norm(_distances[vi*M+dir]);
      if (base_vi != -1) {
        if (!_vertices->IsBase(base_vi)) {
          throw logic_error("base_vi is not a base");
        }
        const int cell_width = _vertices->CellWidth(base_vi);
        if (dist < cell_width * factors[i]) {
          return true;
        }
      }
    }
    return false;
  }

  inline bool IsBase(const int vi) const {
    if (D == 3) {
      return _vertices->IsBase(vi);
    }
    return _vertices->IsCorner(vi, Normal);
  }

  void SetCardinalOnly(bool b) {
    _cardinal_only = b;
  }
  
  vector<vector<int> > GetBase2Incident() const {
    // Now, for each base vertex, store all vertices incident to its cell.
    vector<vector<int> > base2incident(_vertices->size());
    for (int vi = 0; vi < _vertices->size(); ++vi) {
      for (int j = 0; j < M; ++j) {
        // const int base_vi = Get(vi, j);
        const int base_vi = _bases[vi*M+j];
        if (base_vi > -1) {
          base2incident[base_vi].push_back(vi);
        }
      }
    }
    return base2incident;
  }

  // dir_N is the direction in native coordinate space
  inline intn GetDist(const int vi, const int dir_N) const {
    if (_n2p[dir_N] == -1) return make_intn(-1);
    return _distances[vi*M + _n2p[dir_N]];
  }

  // dir_N is the direction in native coordinate space
  inline int GetBase(const int vi, const int dir_N) const {
    if (_n2p[dir_N] == -1) return -1;
    return _bases[vi*M + _n2p[dir_N]];
  }

  inline bool HasBase(const int vi, const int dir_N) const {
    return GetBase(vi, dir_N) > -1;
  }

  // Finds the first base available among:
  //   n
  //   n | axis1
  //   n | axis2
  //   n | axis1 | axis2
  inline int GetBase(
      const int vi, const int n, const int axis1, const int axis2) const {
    int base_vi = GetBase(vi, n);
    if (base_vi == -1 && !IsBoundary2(vi, axis1)) {
      base_vi = GetBase(vi, n | axis1);
    }
    if (base_vi == -1 && !IsBoundary2(vi, axis2)) {
      base_vi = GetBase(vi, n | axis2);
    }
    if (base_vi == -1) {
      base_vi = GetBase(vi, n | axis1 | axis2);
    }
    return base_vi;
  }

  inline int GetNormal() const {
    return Normal;
  }

  inline bool Valid2Direction(
      // const int vi, const Vec<index_t, N>& p,
      const int vi, const intn& p,
      const int iaxis_p) const {
    if (Normal == 0) return true;

    const int jaxis_p = MASK - iaxis_p;
    const int iaxis_n = _p2n[iaxis_p];
    const int jaxis_n = _p2n[jaxis_p];

    const int base_axis = _fbv3->GetBase(vi, iaxis_n);
    if (base_axis > -1) {
      return true;
    }
    const int base_diag = _fbv3->GetBase(vi, iaxis_n | jaxis_n);
    if (base_diag > -1) {
      const intn dist = _fbv3->GetDist(vi, iaxis_n | jaxis_n);
      const int cell_width = _vertices->CellWidth(base_diag);
      if (dist.s[AXIS_IDX[jaxis_n]] < cell_width) {
        return true;
      }
    }
    
    if (p.s[AXIS_IDX[Normal]] == kWidth) {
      return true;
    }

    return false;
  }

 private:
  // D == 2:
  //   x    = 001
  //   y    = 010
  //   M    = 100
  //   MASK = 011
  // D == 3:
  //   x    = 0001
  //   y    = 0010
  //   z    = 0100
  //   M    = 1000
  //   MASK = 0111
  // M (Max) is z*2 so as to account for diagonals (e.g. xz=0101)
  // D is the computing dimension
  // N is the native dimension
  // D = N-1 if we're computing bases on a plane in 3D
  int Normal;
  int D;
  int M;
  int MASK;

  shared_array<int> _bases;
  shared_array<intn> _distances;
  const VertexNetwork* _vertices;
  bool _cardinal_only;

  // For use in vertex visitor mode
  shared_array<bool> _visited;

  // n = native
  // p = projected
  shared_array<int> _p2n;
  shared_array<int> _n2p;

  const FindBasesVisitor<N>* _fbv3;
};

} // end namespace

#endif
