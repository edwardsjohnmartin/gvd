#include "./OctreeUtils.h"
#include "./karras.h"

#include <iostream>

#include "./bb.h"

using std::vector;
using std::cout;
using std::endl;
using std::logic_error;

namespace Karras {

OctCell FindLeaf(
    const intn& p, const vector<OctNode>& octree, const Resln& resln) {
  const int z = xyz2z(p, resln);

  // Set up mask
  int mask = 0;
  for (int i = 0; i < DIM; ++i) {
    mask |= (0x1 << i);
  }

  intn origin = make_uni_intn(0);
  int width = resln.width;
  OctNode const * node = &octree[0];
  for (int i = resln.mbits-DIM; i > 0; i-=DIM) {
    const int octant = (z >> i) & mask;
    width /= 2;

    if (octant % 2 == 1)
      origin += make_intn(width, 0);
    if (octant / 2 == 1)
      origin += make_intn(0, width);

    if (node->is_leaf(octant)) {
      // Empty leaf node
      return OctCell(origin, width, node, octant, 0, (*node)[octant]);
    }
    const int idx = (*node)[octant];
    node = &octree[idx];
  }
  throw logic_error("Didn't find leaf node");
  // return node;
}

OctCell FindNeighbor(
    const OctCell& cell, const int intersection,
    const std::vector<OctNode>& octree, const Resln& resln) {
  throw logic_error("Not implemented");
}

std::vector<intn> FindIntersections(
    const intn& a, const intn& b, const OctCell& cell,
    const std::vector<OctNode>& octree, const Resln& resln) {
  floatn af = convert_floatn(a);
  floatn bf = convert_floatn(b);
  floatn vf = bf - af;
  float len = length(vf);
  vf = vf / len;

  vector<intn> ret;

  BoundingBox<intn> bb;
  bb(cell.get_origin());
  bb(cell.get_origin() + make_uni_intn(cell.get_width()));

  vector<float> t_values;
  vector<intn> p_values;
  {
    // Bottom edge
    const int y = cell.get_origin().y;
    const float t = (y - af.y) / vf.y;
    t_values.push_back(t);
    const intn p = make_intn(af.x + vf.x * t, y);
    p_values.push_back(p);
  } {
    // Right edge
    const int x = cell.get_origin().x + cell.get_width();
    const float t = (x - af.x) / vf.x;
    t_values.push_back(t);
    const intn p = make_intn(x, af.y + vf.y * t);
    p_values.push_back(p);
  } {
    // Top edge
    const int y = cell.get_origin().y + cell.get_width();
    const float t = (y - af.y) / vf.y;
    t_values.push_back(t);
    const intn p = make_intn(af.x + vf.x * t, y);
    p_values.push_back(p);
  } {
    // Left edge
    const int x = cell.get_origin().x;
    const float t = (x - af.x) / vf.x;
    t_values.push_back(t);
    const intn p = make_intn(x, af.y + vf.y * t);
    p_values.push_back(p);
  }

  for (int i = 0; i < t_values.size(); ++i) {
    const float t = t_values[i];
    const intn p = p_values[i];
    if (t >= 0 && t < len) {
      // const intn p = convert_intn(af + vf * t);
      if (bb.in_closed(p)) {
        ret.push_back(p);
      } else if (i == 1) {
        cout << bb << " " << p << endl;
      }
    }
  }

  // p = af + t * vf
  // (p.x - af.x) / vf.x = t

  return ret;
}

} // namespace
