#include "./OctreeUtils.h"
#include "./karras.h"

#include <iostream>
#include <fstream>
#include <set>

#include "./bb.h"
#include "./opencl/geom.h"

using std::vector;
using std::cout;
using std::endl;
using std::logic_error;

namespace Karras {

OctCell FindLeaf(
    const intn& p, const vector<OctNode>& octree, const Resln& resln) {
  const Morton z = xyz2z(p, resln);

  // Set up mask
  int mask = 0;
  for (int i = 0; i < DIM; ++i) {
    mask |= (0x1 << i);
  }

  intn origin = make_uni_intn(0);
  int width = resln.width;
  OctNode const * node = &octree[0];
  int idx = 0;
  for (int i = resln.mbits-DIM; i >= 0; i-=DIM) {
    const int octant = (z >> i).getBlock(0) & mask;
    width /= 2;

    if (octant % 2 == 1)
      origin += make_intn(width, 0);
    if (octant / 2 == 1)
      origin += make_intn(0, width);

    if (node->is_leaf(octant)) {
      // Empty leaf node
      return OctCell(origin, width, idx, node, octant, 0,
                     (*node)[octant]);
    }
    idx = (*node)[octant];
    node = &octree[idx];
  }

  ofstream out("find.err");
  out << p << endl;
  out << resln << endl;
  out << octree << endl;
  out.close();

  cerr << "Didn't find leaf node" << endl;
  cerr << "p = " << p << endl;
  cerr << "resln = " << resln << endl;
  // cerr << "octree = " << octree.size() << " " << octree.back() << endl;
  cerr << "octree = " << octree << endl;

  throw logic_error("Didn't find leaf node");
}

OctCell FindNeighbor(
    const OctCell& cell, const int intersection,
    const std::vector<OctNode>& octree, const Resln& resln) {
  throw logic_error("Not implemented");
}

std::vector<CellIntersection> FindIntersections(
    const floatn& a, const floatn& b, const OctCell& cell,
    const std::vector<OctNode>& octree, const Resln& resln) {
  // static const float EPSILON = 1e-6;
  floatn af = a;
  floatn bf = b;
  // v is normalized (3 lines down)
  floatn v = bf - af;
  float len = length(v);
  v = v / len;

  const intn& origin = cell.get_origin();
  const int width = cell.get_width();
  BoundingBox<floatn> bb;
  bb(convert_floatn(origin));
  bb(convert_floatn(origin) + make_uni_floatn(width));

  vector<float> t_values;
  vector<floatn> p_values;
  {
    // Bottom edge
    const int y = origin.y;
    const float t = (y - af.y) / v.y;
    t_values.push_back(t);
    const floatn p = make_floatn(af.x + v.x * t, y);
    p_values.push_back(p);
  } {
    // Right edge
    const int x = origin.x + width;
    const float t = (x - af.x) / v.x;
    t_values.push_back(t);
    const floatn p = make_floatn(x, af.y + v.y * t);
    p_values.push_back(p);
  } {
    // Top edge
    const int y = origin.y + width;
    const float t = (y - af.y) / v.y;
    t_values.push_back(t);
    const floatn p = make_floatn(af.x + v.x * t, y);
    p_values.push_back(p);
  } {
    // Left edge
    const int x = origin.x;
    const float t = (x - af.x) / v.x;
    t_values.push_back(t);
    const floatn p = make_floatn(x, af.y + v.y * t);
    p_values.push_back(p);
  }

  vector<CellIntersection> ret;
  std::set<floatn> point_set;
  for (int i = 0; i < t_values.size(); ++i) {
    const float t = t_values[i];
    const floatn p = p_values[i];
    // const float x = p.s[0];
    // const float y = p.s[1];
    if (t >= -EPSILON && t < len+EPSILON) {
      if (bb.in_closed(p, EPSILON)) {
        ret.push_back(CellIntersection(t, p));
        point_set.insert(p);
      }
    }
  }
  
  if (ret.size() == 2) {
    if (dist(ret[0].p, ret[1].p) < EPSILON) {
      ret.resize(1);
    }
  }

  if (ret.size() > 2) {
    // If there are 3 or more intersections then choose the two that
    // most closely match the direction v of the intersecting line
    CellIntersection best[2];
    float best_dot = 0;
    for (int i = 0; i < ret.size(); ++i) {
      for (int j = i+1; j < ret.size(); ++j) {
        const floatn c = ret[i].p;
        const floatn d = ret[j].p;
        const floatn u = normalize(d-c);
        const float cur_dot = fabs(dot(u, v));
        if (cur_dot > best_dot) {
          best_dot = cur_dot;
          best[0] = ret[i];
          best[1] = ret[j];
        }
      }
    }
    ret.clear();
    ret.push_back(best[0]);
    ret.push_back(best[1]);
  }

  return ret;
}

} // namespace
