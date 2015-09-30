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
  CellIntersection() {}
  CellIntersection(const float t_, const floatn p_)
      : t(t_), p(p_) {}
  float t;
  floatn p;
};

// Find intersections of the line segment ab with an octree cell.
std::vector<CellIntersection> FindIntersections(
    // const intn& a, const intn& b, const OctCell& cell,
    const floatn& a, const floatn& b, const OctCell& cell,
    const std::vector<OctNode>& octree, const Resln& resln);

} // namespace

#endif
