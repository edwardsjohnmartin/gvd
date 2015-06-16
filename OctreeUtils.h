#ifndef __OCTREE_UTILS_H__
#define __OCTREE_UTILS_H__

#include <vector>

#include "./opencl/defs.h"
#include "./opencl/vec.h"
#include "./OctNode.h"
#include "./Resln.h"

namespace Karras {

OctNode const* FindNode(
    const intn& p, const std::vector<OctNode>& octree, const Resln& resln);

} // namespace

#endif
