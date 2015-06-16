#include "./OctreeUtils.h"
#include "./karras.h"

#include <iostream>

using std::vector;
using std::cout;
using std::endl;

namespace Karras {

OctNode const* FindNode(
    const intn& p, const vector<OctNode>& octree, const Resln& resln) {
  const int z = xyz2z(p, resln);
  int mask = 0;
  for (int i = 0; i < DIM; ++i) {
    mask |= (0x1 << i);
  }
  cout << "mask = " << mask << endl;
  OctNode const * node = &octree[0];
  for (int i = resln.mbits-DIM; i > 0; i-=DIM) {
    const int octant = (z >> i) & mask;
    const int idx = (*node)[octant];
    if (idx == -1) {
      // Empty leaf node
      return node;
    }
    node = &octree[idx];
  }
  return node;
}

} // namespace
