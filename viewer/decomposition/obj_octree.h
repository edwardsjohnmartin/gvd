#ifndef OBJ_OCTREE_H
#define OBJ_OCTREE_H


#include <vector>

#include "../../vertex_network.h"
#include "../../options.h"
#include "../../bb.h"
#include "../../octree.h"
#include "../../geometry_cpp.h"
#include "../../vectorn.h"


class ObjectOctree : public oct::VertexNetwork
{

  private:
    oct::OctreeOptions options;

  public:
    ObjectOctree(oct::OctreeOptions& options);

    void build(
        const std::vector<oct::LabeledGeometry>& vertices,
        oct::GeomVertices& geom_vertices,
        const BoundingBox<floatn>& bounding_box);

};



#endif
