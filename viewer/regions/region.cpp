#include "region.h"


/* Constructor:
 * Create a new region with the normal (which should be the average normal
 * spanning all vertices of this region), and the given unique id.
 * Vertices should be inserted separetely using the addVertex or addVerts
 * functions.
 */
Region::Region(float2 normal, int id)
{
    average_normal = normal;
    unique_id = id;
}


/* Add a vertex to this region, making it a part of the region.
 */
void Region::addVertex(int vert_id)
{
    vertices.push_back(vert_id);
}


/* Add a list of vertices to this region, making all of them a part of this
 * region.
 */
void Region::addVerts(const std::vector<int> &verts)
{
    vertices.insert(vertices.end(), verts.begin(), verts.end());
}


/* Combines this region with the other region, assigning it the average normal
 * between the two regions (weighted appropriately by number of vertices),
 * and adds the vertices from both regions into the new one.
 * The given new_id is assigned as the new region's unique id, and should be
 * different from the IDs of any other regions.
 * Returns the newly created region.
 */
Region Region::merge(Region &other, int new_id)
{
    float2 n1 = getNormal() * getNumVerts();
    float2 n2 = other.getNormal() * other.getNumVerts();
    
    Region mergedRegion(((n1 + n2) / 2.0), new_id);
    mergedRegion.addVerts(vertices);
    mergedRegion.addVerts(other.vertices);

    return mergedRegion;
}


/* Getter: normal
 * Returns the average normal of this region.
 */
float2 Region::getNormal()
{
    return average_normal;
}


/* Getter: number of vertices
 * Returns the number of vertices encompassed by this region.
 */
int Region::getNumVerts()
{
    return vertices.size();
}


/* Getter: unique id
 * Returns the unique id of this region.
 */
int Region::getId()
{
    return unique_id;
}
