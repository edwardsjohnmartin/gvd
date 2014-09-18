#include "./tile3.h"

#include <vector>

#include "./opencl/vec.h"

#include "./octree.h"
#include "./medial.h"
#include "./search.h"

bool TileSurfaceVisitorGpu2::operator()(const int vi, const int3& p) {
  const int l = _vertices->Label(vi);
  // for (int axis = 1; axis < M; axis = (axis<<1)) {
  for (int axis_idx = 0; axis_idx < 3; ++axis_idx) {
    const int axis = (1<<axis_idx);
    const Direction_t dir = oct::DirectionFromPosNeg(0, axis);
    const int n_vi = _vertices->Neighbor(vi, dir);
    if (n_vi != -1) {
      const int n_l = _vertices->Label(n_vi);
      if (l != n_l) {
        // TODO: compute the actual intersection point
        int3 q = p;
        const double d = _vertices->NeighborDist(vi, dir);
        q.s[axis_idx] -= d;
        int3 intersection = make_int3(0);
        int dist;
        oct::tie(intersection, dist) = oct::ComputeIntersection<3>(
            vi, n_vi, p, q, *_vertices, _o);
        if (dist < 0) throw logic_error("computeintersection bad");

        // Add the intersection point in each adjacent plane for the
        // 2D mesh
        oct::shared_array<ViNormal> shared(new ViNormal[4]);
        // Get four cells for which (vi, n_vi) is on the boundary.
        // If axis is x, then there are two shared cells in each of the
        // two containing planes: xy and xz.
        const int axis1_idx = (axis_idx+1)%3;
        const int axis2_idx = (axis_idx+2)%3;
        const int axis1 = (1<<axis1_idx);
        const int axis2 = (1<<axis2_idx);
        // normal is axis2
        GetSharedCells(
            vi, n_vi, axis, axis1, &_fbv_2[axis2_idx], shared.get());
        // normal is axis1
        GetSharedCells(
            vi, n_vi, axis, axis2, &_fbv_2[axis1_idx], shared.get()+2);
        // Map an incident cell/plane to the intersection
        for (int i = 0; i < 4; ++i) {
          const int p_base_vi = shared[i].vi;
          const int p_normal = shared[i].normal;
          if (p_normal == axis) throw logic_error("normal == axis");
          if (p_base_vi != -1) {
            const bool next_axis = (i/2 == 0);
            const bool pos = (i%2 == 0);
            int l0 = l;
            int l1 = n_l;
            // Correct orientation
            if (next_axis == pos)
              swap(l0, l1);

            _p_cell2ints[p_normal][p_base_vi].push_back(
                LabeledIntersection(intersection, l0, l1, dist));
          }
        }
      }
    }
  }
  return true;
}

//                 ---------
//                |         |
//                |         | <-- shared[0]
//                |         |
// oaxis     n_vi |---------| vi
//    ^           |         |
//    |           |         | <-- shared[1]
//    |           |         |
//    |           |---------|
//    |   
//    |   
//  --|------------------>
//    |                 axis
//
// Given two vertices sharing an edge, find the two 2D cells sharing
// the edge.  The first cell is in the positive direction and the second
// is in the negative direction.
//
// n_vi is in the negative axis direction from vi.
// shared must have space for two adjacent cells.
// axis is the primary axis: n_vi is a neighbor of vi along axis.
// oaxis is orthogonal to axis and normal = axis x oaxis.
void TileSurfaceVisitorGpu2::GetSharedCells(
    const int vi, const int n_vi, const int axis, const int oaxis,
    const oct::FindBasesVisitor<3>* fbv,
    ViNormal* shared) {
  const int normal = fbv->GetNormal();
  const int base_pos = fbv->GetBase(vi, axis);
  const int base_neg = fbv->GetBase(vi, axis | oaxis);
  shared[0] = ViNormal(base_pos, normal);
  shared[1] = ViNormal(base_neg, normal);
}

bool TileSurfaceVisitorGpu3::operator()(const int vi, const int3& p) {

  //              |                   |
  //          ____|____*__________    |
  //         |    |               |   |
  //         |    |               *   |
  //         |    +_______________|___|
  //         *   /    ^  ^ normal |  /
  //         |  /axis1| /         | /
  //         | /      |/___> axis2|/
  //     |   vi______________|____|
  //     |   /               |   /
  //     |  /                |  /
  //     | /                 | /
  //     |/__________________|/
  //     -
  //
  // * = intersections
  // + = base_vi_pos
  // - = base_vi_neg

  for (int normal_idx = 0; normal_idx < 3; ++normal_idx) {
    const int axis1_idx = (normal_idx+1)%3;
    const int axis2_idx = (normal_idx+2)%3;
    const int normal = (1<<normal_idx);
    const int axis1 = (1<<axis1_idx);
    const int axis2 = (1<<axis2_idx);

    // intersections is all intersections on edges lying on the 2D cell
    // with base vertex vi and the given normal.
    const vector<LabeledIntersection>& intersections =
        _p_cell2ints[normal][vi];
    const int n = intersections.size();
    if (n > 0) {
      // Find the center of the 2D cell's intersections.  A 2D cell is
      // defined by a base vertex and a normal.

      double3 sum = make_double3(0);
      double dist_sum = 0;
      for (int i = 0; i < n; ++i) {
        sum += convert_double3(intersections[i].p);
        dist_sum += intersections[i].dist;
      }
      const int3 center = convert_int3(sum / n);
      const int center_dist = (int)(dist_sum / n);
      
      // Find out which 3D cells vi is incident to.  There will
      // be at most two: one on the negative normal side and one on
      // the positive normal side.

      // On the negative side, labels are ordered in clockwise direction.
      // So the first label should be 1,2,3 and opposite for the second.
      // On the positive side, labels are ordered in anticlockwise direction.
      // So the first label should be 3,2,1 and opposite for the second.

      const int base_vi_neg = _fbv->GetBase(vi, normal, axis1, axis2);
      const int base_vi_pos = _fbv->GetBase(vi, 0, axis1, axis2);

      if (n > 2) {
        // Add intersections to tri_verts and count the number of
        // intersections for each label
        vector<int> label2count(_num_meshes, 0);
        for (int i = 0; i < n; ++i) {
          const int3 p = intersections[i].p;
          const int dist = intersections[i].dist;
          for (int j = 0; j < 2; ++j) {
            const int label = intersections[i].label[j];
            if (label >= 0 && label < _num_meshes) {
              _tri_verts[label].push_back(p);
              _vert_dist[label].push_back(dist);
              if (dist < 0) throw logic_error("dist bad");
              ++label2count[label];
            } else {
              // This is a consequence of "issue 1", noted in SetDistances
              // (found in octree.cpp)
              cerr << "Bad label" << endl;
            }
          }
        }
        // Add center to any tri_verts arrays that have an intersection
        for (int label = 0; label < _num_meshes; ++label) {
          if (label2count[label] > 0) {
            _tri_verts[label].push_back(center);
            _vert_dist[label].push_back(center_dist);
            if (center_dist < 0) throw logic_error("center_dist bad");
          }
        }

        // Add the intersection to cell2edges
        vector<int> label2added(_num_meshes, 0);
        for (int i = 0; i < n; ++i) {
          const int l0 = intersections[i].label[0];
          const int l1 = intersections[i].label[1];
          if (l0 >= 0 && l0 < _num_meshes &&
              l1 >= 0 && l1 < _num_meshes) {
            _gvd_graph->AddEdge(l0, l1);
            const int idx0 = _tri_verts[l0].size() - label2count[l0] - 1;
            const int idx1 = _tri_verts[l1].size() - label2count[l1] - 1;
            const int off0 = label2added[l0]++;
            const int off1 = label2added[l1]++;
            // intersection index and center index, respectively
            const int idxi0 = idx0+off0;
            const int idxi1 = idx1+off1;
            const int idxc0 = _tri_verts[l0].size() - 1;
            const int idxc1 = _tri_verts[l1].size() - 1;

            if (base_vi_neg > -1) {
              // Edge is (2d_center, 1d_center)
              // The first label will be (3d_center, 2d_center, 1d_center)
              // The second label will be the opposite
              _cell2edges[l0][base_vi_neg].push_back(
                  LabeledEdge(make_edge(idxi0, idxc0), l1));
              _cell2edges[l1][base_vi_neg].push_back(
                  LabeledEdge(make_edge(idxc1, idxi1), l0));
            }
            if (base_vi_pos > -1) {
              // Edge is (2d_center, 1d_center)
              // The second label will be (3d_center, 2d_center, 1d_center)
              // The first label will be the opposite
              _cell2edges[l0][base_vi_pos].push_back(
                  LabeledEdge(make_edge(idxc0, idxi0), l1));
              _cell2edges[l1][base_vi_pos].push_back(
                  LabeledEdge(make_edge(idxi1, idxc1), l0));
            }
          }
        }
      } else if (n == 2) {
        const int l0 = intersections[0].label[0];
        const int l1 = intersections[1].label[0];
        if (l0 >= 0 && l0 < _num_meshes &&
            l1 >= 0 && l1 < _num_meshes) {
          _gvd_graph->AddEdge(l0, l1);

          const int3 p0 = intersections[0].p;
          const int3 p1 = intersections[1].p;
          const int dist0 = intersections[0].dist;
          const int dist1 = intersections[1].dist;
          const int idx0 = _tri_verts[l0].size();
          const int idx1 = _tri_verts[l1].size();
          _tri_verts[l0].push_back(p0);
          _tri_verts[l0].push_back(p1);
          _tri_verts[l1].push_back(p0);
          _tri_verts[l1].push_back(p1);
          if (dist0 < 0) throw logic_error("dist0 bad");
          if (dist1 < 0) throw logic_error("dist1 bad");
          _vert_dist[l0].push_back(dist0);
          _vert_dist[l0].push_back(dist1);
          _vert_dist[l1].push_back(dist0);
          _vert_dist[l1].push_back(dist1);

          if (base_vi_neg > -1) {
            _cell2edges[l0][base_vi_neg].push_back(
                LabeledEdge(make_edge(idx0, idx0+1), l1));
            _cell2edges[l1][base_vi_neg].push_back(
                LabeledEdge(make_edge(idx1+1, idx1), l0));
          }
          if (base_vi_pos > -1) {
            _cell2edges[l0][base_vi_pos].push_back(
                LabeledEdge(make_edge(idx0+1, idx0), l1));
            _cell2edges[l1][base_vi_pos].push_back(
                LabeledEdge(make_edge(idx1, idx1+1), l0));
          }
        } else {
          // This is a consequence of "issue 1", noted in SetDistances
          // (found in octree.cpp)
          cerr << "Bad label" << endl;
        }
      }
    }
  }
  return true;
}
