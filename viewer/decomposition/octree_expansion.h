/**
 *
 */
#ifndef OCTREE_EXPANSION_H_
#define OCTREE_EXPANSION_H_


template <int DIM>
class OctreeExpansion
{

  public:
    
    OctreeExpansion()
    {
        assert(DIM == 2 || DIM == 3);
    }


    /**
     * Returns true if the provided N-dimensional vertex belongs to the provided
     * N-dimensional octree cell.
     */
    bool vertexBelongsToCell(;)
    {
        ;
    }
    
    
    /**
     * Given a GVD octree, expands the tree into D dimensiona
     */
    void expandCells(VertexNetwork octree_verts,
                     std::vector<MyVec<float, DIM> > object_verts)
    {
        ;
    }

}


#endif
