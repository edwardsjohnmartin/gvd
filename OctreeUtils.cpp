#include "./OctreeUtils.h"
#include "./karras.h"

#include <iostream>
#include <fstream>
#include <set>

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
  int idx = 0;
  for (int i = resln.mbits-DIM; i >= 0; i-=DIM) {
    const int octant = (z >> i) & mask;
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
    const intn& a, const intn& b, const OctCell& cell,
    const std::vector<OctNode>& octree, const Resln& resln) {
  floatn af = convert_floatn(a);
  floatn bf = convert_floatn(b);
  floatn vf = bf - af;
  float len = length(vf);
  vf = vf / len;

  vector<CellIntersection> ret;

  BoundingBox<floatn> bb;
  bb(convert_floatn(cell.get_origin()));
  bb(convert_floatn(cell.get_origin()) + make_uni_floatn(cell.get_width()));

  vector<float> t_values;
  vector<floatn> p_values;
  {
    // Bottom edge
    const int y = cell.get_origin().y;
    const float t = (y - af.y) / vf.y;
    t_values.push_back(t);
    const floatn p = make_floatn(af.x + vf.x * t, y);
    p_values.push_back(p);
  } {
    // Right edge
    const int x = cell.get_origin().x + cell.get_width();
    const float t = (x - af.x) / vf.x;
    t_values.push_back(t);
    const floatn p = make_floatn(x, af.y + vf.y * t);
    p_values.push_back(p);
  } {
    // Top edge
    const int y = cell.get_origin().y + cell.get_width();
    const float t = (y - af.y) / vf.y;
    t_values.push_back(t);
    const floatn p = make_floatn(af.x + vf.x * t, y);
    p_values.push_back(p);
  } {
    // Left edge
    const int x = cell.get_origin().x;
    const float t = (x - af.x) / vf.x;
    t_values.push_back(t);
    const floatn p = make_floatn(x, af.y + vf.y * t);
    p_values.push_back(p);
  }

  std::set<intn> point_set;
  for (int i = 0; i < t_values.size(); ++i) {
    const float t = t_values[i];
    const floatn p = p_values[i];
    if (t >= 0 && t < len) {
      if (bb.in_closed(p)) {
        const intn qp = convert_intn(p);
        if (point_set.find(qp) == point_set.end()) {
          ret.push_back(CellIntersection(t, qp));
          point_set.insert(qp);
        }
      }
    }
  }

  // p = af + t * vf
  // (p.x - af.x) / vf.x = t

  return ret;
}

} // namespace
