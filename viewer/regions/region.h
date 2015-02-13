/*
 */

#ifndef REGION_H
#define REGION_H

#include <vector>
#include "../../opencl/vec_cpp.h"


class Region {

  private:
    float2 average_normal;
    int unique_id;
    std::vector<int> vertices; // TODO - use different data struct
    int left_neighbor, right_neighbor;

  public:
    Region(float2 normal, int id);

    void addVertex(int vert_id);
    void addVerts(const std::vector<int> &verts);

    void setNeighbors(int left_neighbor, int right_neighbor);

    Region merge(Region &other, int new_id);
    double angleTo(Region &other);

    float2 getNormal();
    int getNumVerts();
    int getId();
    int getLeftNeighbor();
    int getRightNeighbor();
};


#endif
