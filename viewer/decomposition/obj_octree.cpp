#include "obj_octree.h"

using namespace oct;


ObjectOctree::ObjectOctree(OctreeOptions& options)
{
    //this->options = options;
}


void ObjectOctree::build(
    const std::vector<LabeledGeometry>& vertices,
    GeomVertices& geom_vertices,
    const BoundingBox<floatn>& bounding_box)
{
}
