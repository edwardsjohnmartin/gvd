/*******************************************************
 ** Generalized Voronoi Diagram Project               **
 ** Copyright (c) 2015 John Martin Edwards            **
 ** Scientific Computing and Imaging Institute        **
 ** 72 S Central Campus Drive, Room 3750              **
 ** Salt Lake City, UT 84112                          **
 **                                                   **
 ** For information about this project contact        **
 ** John Edwards at                                   **
 **    edwardsjohnmartin@gmail.com                    **
 ** or visit                                          **
 **    sci.utah.edu/~jedwards/research/gvd/index.html **
 *******************************************************/


// //------------------------------------------------------------
// // ConnCompEdges
// //
// // Finds all pairs of touching labels.
// //------------------------------------------------------------
// template <int D>
// std::set<Edge> ConnCompEdges(const VertexNetwork& vertices,
//                              const int start_vi,
//                              const intn& start_p,
//                              // const BoundingBox<index_t, D>& bb,
//                              const BoundingBox<intn>& bb,
//                              int& num_labels,
//                              const std::vector<int>& axes);

// //------------------------------------------------------------
// // ConnCompEdges
// //------------------------------------------------------------
// template <int D>
// std::set<Edge> ConnCompEdges(const VertexNetwork& vertices,
//                              const int start_vi,
//                              const intn& start_p,
//                              // const BoundingBox<index_t, D>& bb,
//                              const BoundingBox<intn>& bb,
//                              int& num_labels);

// //------------------------------------------------------------
// // ConnComp
// //
// // Performs a breadth-first search over the vertices finding
// // the largest connected component starting from start_vi.
// // VertexVisitor is of the form
// //    bool Visit(const int vi, const int3& p);
// // EdgeVisitor is of the form
// //    bool Visit(const int vi, const int n_vi,
// //               const int3& p, const int3& q, const Direction<D>& d);
// //------------------------------------------------------------
// template <int D>
// std::vector<int> ConnComp(const VertexNetwork& vertices,
//                           const int start_vi,
//                           const intn& start_p,
//                           // const BoundingBox<index_t, D>& bb,
//                           const BoundingBox<intn>& bb,
//                           const std::vector<int>& axes);

// //------------------------------------------------------------
// // ConnComp
// //
// // Performs a breadth-first search over the vertices finding
// // the largest connected component starting from start_vi.
// // VertexVisitor is of the form
// //    bool Visit(const int vi, const int3& p);
// // EdgeVisitor is of the form
// //    bool Visit(const int vi, const int n_vi,
// //               const int3& p, const int3& q, const Direction<D>& d);
// //------------------------------------------------------------
// template <int D>
// std::vector<int> ConnComp(const VertexNetwork& vertices,
//                           const int start_vi,
//                           const intn& start_p,
//                           // const BoundingBox<index_t, D>& bb);
//                           const BoundingBox<intn>& bb);
