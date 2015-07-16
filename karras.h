#ifndef __KARRAS_H__
#define __KARRAS_H__

#include <vector>
#include <stdexcept>

#include "./opencl/defs.h"
#include "./opencl/vec.h"
#include "./Resln.h"
#include "./OctNode.h"

namespace Karras {

Morton xyz2z(intn p, const Resln& r);
intn z2xyz(const Morton z, const Resln& r);

std::vector<intn> Quantize(const std::vector<floatn>& points, const Resln& r);

std::vector<OctNode> BuildOctree(
    const std::vector<intn>& opoints, const Resln& r, const bool verbose=false);

// Debug output
void OutputOctree(const std::vector<OctNode>& octree);

} // namespace

// inline std::ostream& operator<<(std::ostream& out, const Karras::Resln& resln) {
//   out << "width=" << resln.width << ", volume=" << resln.volume
//       << ", bits=" << resln.bits
//       << ", mbits=" << resln.mbits;
//   return out;
// }

#endif
