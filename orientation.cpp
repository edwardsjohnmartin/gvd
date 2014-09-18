#include <vector>
#include <queue>
#include <set>
#include <iomanip>

#include "./opencl/triangle.h"
#include "./opencl/edge.h"
#include "./opencl/vec.h"

#include "./triangle_cpp.h"
#include "./edge_cpp.h"

using namespace std;

vector<Triangle> adj_triangles(const Edge& e,
                                const vector<vector<Triangle> >& v2t)
{
  std::vector<Triangle> adj;
  const std::vector<Triangle>& triangles = v2t[e.s[0]];
  // foreach(const Triangle& t, v2t[e[0]]) {
  typedef vector<Triangle>::const_iterator Iter;
  for (Iter it = triangles.begin(); it != triangles.end(); ++it) {
    const Triangle& t = *it;
    if (triangle_has_edge(e.s[0], e.s[1], t))
      adj.push_back(t);
  }
  return adj;
}

Triangle adj_triangle(const Triangle& t, int vi,
                      const vector<vector<Triangle> >& v2t)
{
  vector<Triangle> tris = adj_triangles(triangle_opposite(vi, t), v2t);
  if (tris[0] != t) return tris[0];
  return tris[1];
}

vector<Triangle> adj_triangles(const Triangle& t,
                                const vector<vector<Triangle> >& v2t)
{
  set<Triangle> tris;
  // foreach(const Edge& e, t.edges()) {
  for (int i = 0; i < 3; ++i) {
    Edge e = make_edge(t.s[i], t.s[(i+1)%3]);
    vector<Triangle> ts = adj_triangles(e, v2t);
    tris.insert(ts.begin(), ts.end());
  }
  tris.erase(t);
  return vector<Triangle>(tris.begin(), tris.end());
}

// Testing:
// 0        1
//
//
// 2        3
//
//
// 4        5
// vector<Triangle> triangles;
// triangles.push_back(Triangle(0, 2, 1));
// triangles.push_back(Triangle(2, 1, 3));
// triangles.push_back(Triangle(2, 3, 4));
// triangles.push_back(Triangle(2, 3, 5));
// for (int i = 0; i < triangles.size(); ++i) {
//   cout << triangles[i] << endl;
// }
// triangles = orient(triangles);
// cout << "fixed" << endl;
// for (int i = 0; i < triangles.size(); ++i) {
//   cout << triangles[i] << endl;
// }
// return 0;

// Surface is assumed to be orientable.
// triangles will likely change order.
vector<Triangle> orient_new(const vector<Triangle>& triangles)
{
  if (triangles.empty()) return vector<Triangle>();

  const vector<vector<Triangle> > v2t =
      build_v2t(triangles.begin(), triangles.end());
  //const int n = v2t.size();

  vector<Triangle> o_triangles;
  set<Triangle> visited;
  // oriented edges
  set<Edge> o_edges;
  vector<Edge> tempe;
  get_edges(triangles[0], back_inserter(tempe));
  // o_edges.insert(tempe.begin(), tempe.end());
  // breadth-first search
  queue<Triangle> q;
  q.push(triangles[0]);
  while (!q.empty()) {
    Triangle t = q.front();
    q.pop();
    if (visited.find(t) == visited.end()) {
      visited.insert(t);

      // find a visited edge and reverse if necessary
      Triangle new_t = t;
      vector<Edge> edges;
      get_edges(t, back_inserter(edges));
      bool done = false;
      for (int i = 0; i < 3 && !done; ++i) {
        const Edge& e = edges[i];
        if (o_edges.find(e) != o_edges.end()) {
          const Edge& oe = *o_edges.find(e);
          if (oe.s[0] == e.s[0]) {
            // need to reverse this triangle
            new_t = make_triangle(t.s[2], t.s[1], t.s[0]);
          }
          done = true;
        }
      }

      // add edges to o_edges and triangle to o_triangles
      // edges = new_t.edges();
      get_edges(new_t, back_inserter(edges));
      o_edges.insert(edges.begin(), edges.end());
      o_triangles.push_back(new_t);

      // add neighbors to the queue
      vector<Triangle> nbrs = adj_triangles(t, v2t);
      for (int i = 0; i < nbrs.size(); ++i) {
        const Triangle& nbr = nbrs[i];
        if (visited.find(nbr) == visited.end()) {
          q.push(nbr);
        }
      }
      // foreach (const Triangle& nbr, nbrs) {
      //   if (visited.find(nbr) == visited.end()) {
      //     q.push(nbr);
      //   }
      // }
    }
  }
  return o_triangles;
}
// Surface is assumed to be orientable.
// triangles will likely change order.
vector<Triangle> orient(const vector<Triangle>& triangles) {
  if (triangles.empty()) return vector<Triangle>();

  const vector<vector<Triangle> > v2t =
      build_v2t(triangles.begin(), triangles.end());
  const map<Edge, vector<Triangle> > e2t =
      build_e2t(triangles.begin(), triangles.end());
  //const int n = v2t.size();

  // Remove triangles as we complete them
  set<Triangle> t_set(triangles.begin(), triangles.end());
  // oriented triangles
  vector<Triangle> o_triangles;
  while (!t_set.empty()) {
    // cout << "starting manifold section" << endl;
    set<Triangle> visited;
    // oriented edges
    set<Edge> o_edges;
    // breadth-first search
    queue<Triangle> q;
    // q.push(triangles[0]);
    q.push(*t_set.begin());
    while (!q.empty()) {
      Triangle t = q.front();
      q.pop();
      if (visited.find(t) == visited.end()) {
        // cout << "visiting " << t << endl;
        visited.insert(t);
        t_set.erase(t);

        // find a visited edge and reverse if necessary
        Triangle new_t = t;
        vector<Edge> edges;
        get_edges(t, back_inserter(edges));
        bool done = false;
        for (int i = 0; i < 3 && !done; ++i) {
          const Edge& e = edges[i];
          if (o_edges.find(e) != o_edges.end()) {
            const Edge& oe = *o_edges.find(e);
            if (oe.s[0] == e.s[0]) {
              // need to reverse this triangle
              new_t = make_triangle(t.s[2], t.s[1], t.s[0]);
            }
            done = true;
          }
        }

        // add edges to o_edges and triangle to o_triangles
        get_edges(new_t, inserter(o_edges, o_edges.begin()));
        o_triangles.push_back(new_t);

        // add neighbors to the queue
        for (int i = 0; i < 3; ++i) {
          Edge e = triangle_opposite(t.s[i], t);
          const vector<Triangle>& nbrs = e2t.find(e)->second;
          if (nbrs.size() == 2) {
            // only traverse manifold edges
            const Triangle& nbr = (nbrs[0]==t)?nbrs[1]:nbrs[0];
            if (visited.find(nbr) == visited.end()) {
              q.push(nbr);
            }
          }
        }

        //   vector<Triangle> nbrs = adj_triangles(t, v2t);
        //   typedef vector<Triangle>::const_iterator Iter;
        //   for (Iter it = nbrs.begin(); it != nbrs.end(); ++it) {
        //     const Triangle& nbr = *it;
        //     if (visited.find(nbr) == visited.end()) {
        //       q.push(nbr);
        //     }
        //   }
      }
    }
  }
  return o_triangles;
}

// Surface is assumed to be orientable.
// triangles will likely change order.
void orient_lean(vector<Triangle>& triangles) {
  const vector<vector<Triangle> > v2t =
      build_v2t(triangles.begin(), triangles.end());
  const map<Edge, vector<Triangle> > e2t =
      build_e2t(triangles.begin(), triangles.end());
  //const int n = v2t.size();

  for (int i = 0; i < triangles.size(); ++i) {
    Triangle& t = triangles[i];
    vector<Edge> edges;
    get_edges(t, back_inserter(edges));
    int num_diffs = 0;
    for (int j = 0; j < edges.size(); ++j) {
      // For each edge
      const Edge& e = edges[j];
      const vector<Triangle>& nbrs = e2t.find(e)->second;
      for (int k = 0; k < nbrs.size(); ++k) {
        const Triangle& nbr = nbrs[k];
        if (nbr != t) {
          if (triangle_next(e.s[0], nbr) == e.s[1]) {
            ++num_diffs;
          }
        }
      }
    }
    if (num_diffs > 1) {
      t = triangle_inverted(t);
    }
  }
}

bool adjacent(const Triangle& T, const Edge& e,
              const vector<vector<Triangle> >& v2t, Triangle& T_)
{
  typedef vector<Triangle>::const_iterator Iter;
  for (Iter it = v2t[e.s[0]].begin();
       it != v2t[e.s[0]].end();
       ++it) {
    const Triangle& t = *it;
  // foreach(const Triangle& t, v2t[e[0]]) {
    if (triangle_has_edge(e.s[0], e.s[1], t) && t != T) {
      T_ = t;
      return true;
    }
  }
  return false;
}

bool same_orientation(const Triangle& t1, const Triangle& t2) {
  int t2_0 = triangle_g2l(t1.s[0], t2);
  if (t1.s[0] == t2.s[t2_0] && t1.s[1] == t2.s[(t2_0+1)%3] && t1.s[1] == t2.s[(t2_0+1)%3])
    return true;
  else
    return false;
}

void check_adj_orientation(const Triangle& cur_tri,
                           set<Triangle>& tri_set,
                           vector<vector<Triangle> >& v2t) {
  for (int li = 0; li < 3; ++li) {
    int startV = cur_tri.s[li]; //global index
    int endV = cur_tri.s[(li+1)%3];  //global index
    Edge e = make_edge(startV, endV);
    Triangle adj_t;
    if (adjacent(cur_tri, e, v2t, adj_t)) {
      if (tri_set.find(adj_t) == tri_set.end()){
        int startV_a = triangle_g2l(startV, adj_t); //local index of adj
        int endV_a = triangle_g2l(endV, adj_t);  //local index of adj
        
        if ( ((endV_a+1)%3) == startV_a) {
          //good      
          tri_set.insert(adj_t);
          check_adj_orientation(adj_t, tri_set, v2t);
        }
        else {
          //reverse
          Triangle tri_rev = make_triangle(adj_t.s[2], adj_t.s[1], adj_t.s[0]);
          //update v2t
          for (int v = 0; v<3; ++v){
            typedef vector<Triangle>::iterator Iter;
            for (Iter it = v2t[adj_t.s[v]].begin();
                 it != v2t[adj_t.s[v]].end();
                 ++it) {
            // BOOST_FOREACH(Triangle& t, v2t[adj_t[v]]) {
              Triangle& t = *it;
              if (t == adj_t){
                t = tri_rev;
              }
            }
          }
          tri_set.insert(tri_rev);
          check_adj_orientation(tri_rev, tri_set, v2t);
        }
      }
      else {
        if (!same_orientation(adj_t, *tri_set.find(adj_t))){
          throw logic_error("unknown error");
        }
      }
    }
    else {
      //No Adjacent Triangle on edge
    }
  }
}

Triangle replace(const Triangle& t, int vi, int vi_new) {
  int v[3];
  for (int i = 0; i < 3; ++i) {
    if (t.s[i] == vi)
      v[i] = vi_new;
    else
      v[i] = t.s[i];
  }
  return make_triangle(v);
}

// Splits a surface into multiple orientable surfaces
vector<Triangle> split(const vector<Triangle>& triangles,
                       vector<float3>& vertices) {
  if (triangles.empty()) return vector<Triangle>();

  const vector<vector<Triangle> > v2t =
      build_v2t(triangles.begin(), triangles.end());
  map<Edge, vector<Triangle> > e2t =
      build_e2t(triangles.begin(), triangles.end());
  //const int n = v2t.size();

  // Remove triangles as we complete them
  set<Triangle> t_set(triangles.begin(), triangles.end());
  // oriented triangles
  vector<Triangle> all_o_triangles;
  while (!t_set.empty()) {
    // cout << "starting manifold section" << endl;
    map<int, int> replaced;
    vector<Triangle> o_triangles;
    set<Triangle> visited;
    // oriented edges
    set<Edge> o_edges;
    // breadth-first search
    queue<Triangle> q;
    // q.push(triangles[0]);
    q.push(*t_set.begin());
    while (!q.empty()) {
      Triangle t = q.front();
      q.pop();
      if (visited.find(t) == visited.end()) {
        // cout << "visiting " << t << endl;
        visited.insert(t);
        t_set.erase(t);

        // find a visited edge and reverse if necessary
        Triangle new_t = t;
        vector<Edge> edges;
        get_edges(t, back_inserter(edges));
        bool done = false;
        for (int i = 0; i < 3 && !done; ++i) {
          const Edge& e = edges[i];
          if (o_edges.find(e) != o_edges.end()) {
            const Edge& oe = *o_edges.find(e);
            if (oe.s[0] == e.s[0]) {
              // need to reverse this triangle
              new_t = make_triangle(t.s[2], t.s[1], t.s[0]);
            }
            done = true;
          }
        }

        // add neighbors to the queue and break non-manifold edges
        for (int i = 0; i < 3; ++i) {
          Edge e = triangle_opposite(t.s[i], t);
          vector<Triangle>& nbrs = e2t.find(e)->second;
          if (nbrs.size() == 2) {
            // only traverse manifold edges
            const Triangle& nbr = (nbrs[0]==t)?nbrs[1]:nbrs[0];
            if (visited.find(nbr) == visited.end()) {
              q.push(nbr);
            }
          } else if (nbrs.size() > 2) {
            // non-manifold edge
            // find this triangle and remove it from the list
            for (int j = 0; j < nbrs.size()-1; ++j) {
              if (nbrs[j] == t)
                swap(nbrs[j], nbrs[nbrs.size()-1]);
            }
            if (nbrs[nbrs.size()-1] != t) {
              cout << "t not found in neighbor list" << endl;
              cout << "  t = " << endl;
              for (int j = 0; j < 3; ++j) {
                cout << "    " << setprecision(15) << vertices[t.s[j]] << endl;
              }
              cout << endl;
              throw logic_error("t not found in neighbor list");
            }
            nbrs.pop_back();

            for (int i = 0; i < 2; ++i) {
              const int oldi = e.s[i];
              int newi;
              if (replaced.find(oldi) != replaced.end()) {
                newi = replaced.find(oldi)->second;
              } else {
                newi = vertices.size();
                replaced[oldi] = newi;
                vertices.push_back(vertices[oldi]);
              }
              for (int j = 0; j < 3; ++j) {
                new_t = replace(new_t, oldi, newi);
              }
            }
          }
        }

        // add edges to o_edges and triangle to o_triangles
        get_edges(new_t, inserter(o_edges, o_edges.begin()));
        o_triangles.push_back(new_t);
      }
    }
    for (int i = 0; i < o_triangles.size(); ++i) {
      Triangle t = o_triangles[i];
      for (int j = 0; j < 3; ++j) {
        if (replaced.find(t.s[j]) != replaced.end()) {
          t = replace(t, t.s[j], replaced.find(t.s[j])->second);
        }
      }
      all_o_triangles.push_back(t);
    }
  }
  return all_o_triangles;
}

