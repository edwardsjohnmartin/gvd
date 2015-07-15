#ifndef __OCTREE_UTILS_H__
#define __OCTREE_UTILS_H__

#include <vector>

#include "./opencl/defs.h"
#include "./opencl/vec.h"
#include "./OctNode.h"
#include "./Resln.h"

namespace Karras {

OctCell FindLeaf(
    const intn& p, const std::vector<OctNode>& octree, const Resln& resln);

OctCell FindNeighbor(
    const OctCell& cell, const int intersection,
    const std::vector<OctNode>& octree, const Resln& resln);

struct CellIntersection {
  CellIntersection(const float t_, const intn p_)
      : t(t_), p(p_) {}
  float t;
  intn p;
};

// Find intersections of the line segment ab with an octree cell.
std::vector<CellIntersection> FindIntersections(
    const intn& a, const intn& b, const OctCell& cell,
    const std::vector<OctNode>& octree, const Resln& resln);

} // namespace

#endif
