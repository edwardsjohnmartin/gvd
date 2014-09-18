#ifndef __OCTREE_GPU_H__
#define __OCTREE_GPU_H__

#include <vector>

#include "./opencl/defs.h"

#include "./mvertex_network.h"
#include "./vectorn.h"
#include "./options.h"

NAMESPACE_OCT_BEGIN

void BuildOctreeGpu(
    const std::vector<LabeledGeometry>& geometries,
    GeomVertices& geom_vertices,
    MVertexNetwork& mvertices,
    const OctreeOptions& o);

NAMESPACE_OCT_END

#endif
