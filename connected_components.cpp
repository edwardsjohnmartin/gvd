
// //------------------------------------------------------------------------------
// // Connected components functions
// //------------------------------------------------------------------------------

// //------------------------------------------------------------
// // ConnCompEdges
// //
// // Finds all pairs of touching labels.
// //------------------------------------------------------------
// template <int D>
// std::set<Edge> ConnCompEdges(const VertexNetwork& vertices,
//                              const int start_vi,
//                              const intn& start_p,
//                              const BoundingBox<intn>& bb,
//                              int& num_labels,
//                              const std::vector<int>& axes) {
//   typedef VisitedSet VisitedStruct;

//   assert(bb.in_closed(start_p));

//   std::set<int> all_labels;
//   std::set<Edge> pairs;

//   if (vertices.empty()) return pairs;

//   const int num_axes = axes.size();
//   VisitedStruct visited(vertices.size());
//   std::list<std::pair<int, intn> > new_search;
//   new_search.push_back(make_pair(start_vi, start_p));
//   while (!new_search.empty()) {
//     // Outer search.  Each one starts a connected component
//     const int s_vi = new_search.front().first;
//     const intn s_p = new_search.front().second;
//     new_search.pop_front();
//     if (!visited[s_vi]) {
//       // components.push_back(label);
//       std::list<std::pair<int, intn> > queue;
//       queue.push_back(make_pair(s_vi, s_p));
//       while (!queue.empty()) {
//         const int vi = queue.front().first;
//         const int label = vertices.Label(vi);
//         const intn p = queue.front().second;
//         queue.pop_front();
//         if (!visited[vi]) {
//           visited.SetVisited(vi);
//           for (int i = 0; i < num_axes; ++i) {
//             const int axis = axes[i];
//             for (int pos = 0; pos < 2; ++pos) {
//               // const oct::Direction<D> d = Direction::FromAxis(axis, pos);
//               const oct::Direction d = DirectionFromAxis(axis, pos);
//               const int n_vi = vertices.Neighbor(vi, d);
//               if (n_vi != -1) {
//                 intn q = p;
//                 if (pos)
//                   q.s[axis] += vertices.NeighborDist(vi, d);
//                 else
//                   q.s[axis] -= vertices.NeighborDist(vi, d);
//                 if (bb.in_closed(q)) {
//                   const int n_label = vertices.Label(n_vi);
//                   if (n_label != label) {
//                     pairs.insert(Edge(n_label, label));
//                   }
//                   all_labels.insert(n_label);
//                   queue.push_back(make_pair(n_vi, q));
//                   // else
//                   //   new_search.push_back(make_pair(n_vi, q));
//                 }
//               }
//             }
//           }
//         }
//       }
//     }
//   }
//   num_labels = all_labels.size();
//   return pairs;
// }

// //------------------------------------------------------------
// // ConnCompEdges
// //------------------------------------------------------------
// template <int D>
// std::set<Edge> ConnCompEdges(const VertexNetwork& vertices,
//                              const int start_vi,
//                              const intn& start_p,
//                              const BoundingBox<intn>& bb,
//                              int& num_labels) {
//   typedef VisitedSet VisitedStruct;

//   assert(bb.in_closed(start_p));

//   std::vector<int> axes;
//   for (int i = 0; i < D; ++i) axes.push_back(i);

//   return ConnCompEdges<D>(vertices, start_vi, start_p, bb, num_labels, axes);
// }

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
//                           const BoundingBox<intn>& bb,
//                           const std::vector<int>& axes) {
//   typedef VisitedSet VisitedStruct;

//   assert(bb.in_closed(start_p));

//   std::vector<int> components;
//   if (vertices.empty()) return components;

//   const int num_axes = axes.size();
//   VisitedStruct visited(vertices.size());
//   std::list<std::pair<int, intn> > new_search;
//   new_search.push_back(make_pair(start_vi, start_p));
//   // Outer search.  Each one starts a connected component.
//   while (!new_search.empty()) {
//     const int s_vi = new_search.front().first;
//     const intn s_p = new_search.front().second;
//     new_search.pop_front();
//     if (!visited[s_vi]) {
//       const int label = vertices.Label(s_vi);
//       components.push_back(label);
//       std::list<std::pair<int, intn> > queue;
//       queue.push_back(make_pair(s_vi, s_p));
//       // Inner search.  Sticks with a label.
//       while (!queue.empty()) {
//         const int vi = queue.front().first;
//         const intn p = queue.front().second;
//         queue.pop_front();
//         if (!visited[vi]) {
//           visited.SetVisited(vi);
//           // Search out over each axis
//           for (int i = 0; i < num_axes; ++i) {
//             const int axis = axes[i];
//             // Each direction
//             for (int pos = 0; pos < 2; ++pos) {
//               // const oct::Direction<D> d = Direction::FromAxis(axis, pos);
//               const oct::Direction d = DirectionFromAxis(axis, pos);
//               const int n_vi = vertices.Neighbor(vi, d);
//               if (n_vi != -1) {
//                 intn q = p;
//                 if (pos)
//                   q.s[axis] += vertices.NeighborDist(vi, d);
//                 else
//                   q.s[axis] -= vertices.NeighborDist(vi, d);
//                 if (bb.in_closed(q)) {
//                   if (vertices.Label(n_vi) == label)
//                     queue.push_back(make_pair(n_vi, q));
//                   else
//                     new_search.push_back(make_pair(n_vi, q));
//                 }
//               }
//             }
//           }
//         }
//       }
//     }
//   }
//   return components;
// }

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
//                           const BoundingBox<intn>& bb) {
//   typedef VisitedSet VisitedStruct;

//   assert(bb.in_closed(start_p));

//   std::vector<int> axes;
//   for (int i = 0; i < D; ++i) axes.push_back(i);

//   return ConnComp<D>(vertices, start_vi, start_p, bb, axes);
// }


#ifdef OCT2D
// template std::set<Edge> ConnCompEdges<2>(
//     const VertexNetwork& vertices,
//     const int start_vi,
//     const int2& start_p,
//     const BoundingBox<int2>& bb,
//     int& num_labels,
//     const std::vector<int>& axes);
// template std::set<Edge> ConnCompEdges<2>(
//     const VertexNetwork& vertices,
//     const int start_vi,
//     const int2& start_p,
//     const BoundingBox<int2>& bb,
//     int& num_labels);
// template std::vector<int> ConnComp<2>(
//     const VertexNetwork& vertices,
//     const int start_vi,
//     const int2& start_p,
//     const BoundingBox<int2>& bb,
//     const std::vector<int>& axes);
// template std::vector<int> ConnComp<2>(
//     const VertexNetwork& vertices,
//     const int start_vi,
//     const int2& start_p,
//     const BoundingBox<int2>& bb);
#endif

#ifdef OCT3D
// template std::set<Edge> ConnCompEdges<3>(
//     const VertexNetwork& vertices,
//     const int start_vi,
//     const int3& start_p,
//     const BoundingBox<int3>& bb,
//     int& num_labels,
//     const std::vector<int>& axes);
// template std::set<Edge> ConnCompEdges<3>(
//     const VertexNetwork& vertices,
//     const int start_vi,
//     const int3& start_p,
//     const BoundingBox<int3>& bb,
//     int& num_labels);
// template std::vector<int> ConnComp<3>(
//     const VertexNetwork& vertices,
//     const int start_vi,
//     const int3& start_p,
//     const BoundingBox<int3>& bb,
//     const std::vector<int>& axes);
// template std::vector<int> ConnComp<3>(
//     const VertexNetwork& vertices,
//     const int start_vi,
//     const int3& start_p,
//     const BoundingBox<int3>& bb);
#endif
