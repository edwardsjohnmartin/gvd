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

#ifndef __OCT_GVD_H__
#define __OCT_GVD_H__

#include <vector>
#include <utility>
#include <set>
#include <algorithm>
#include <cfloat>

#include "./opencl/edge.h"
#include "./opencl/vertex.h"
#include "./search.h"

namespace oct {

//----------------------------------------
// class LabeledVertex
//----------------------------------------
template <int D>
class LabeledVertex {
 public:
  LabeledVertex() {}
  LabeledVertex(const intn& point, const int vi, const int label)
      : _point(point), _vi(vi), _label(label) {}

  const intn& Point() const { return _point; }
  int Index() const { return _vi; }
  const int Label() const { return _label; }

 private:
  intn _point;
  int _vi;
  int _label;
  // 16 bytes in 3D
};

template <int D>
struct LabeledSegment {
  LabeledSegment() {}
  LabeledSegment(const LabeledVertex<D>& p1, const LabeledVertex<D>& p2) {
    endpoints[0] = p1;
    endpoints[1] = p2;
  }
  const LabeledVertex<D>& operator[](int i) const { return endpoints[i]; }
  LabeledVertex<D> endpoints[2];
};

// Takes intersection points around a 2D face and connects them
// in some meaningful way.  If midpoint is true, connects to the
// centroid of the face.  Otherwise, recursively connects the closest
// pair.
//    |-----------*------|
//    |           |      |
//    |          /       |
//    |         |        |
//    *--------/         |
//    |        |         |
//    |       / \        |
//    |      /   \       |
//    |-----*-----*------|
//
// Each line is also given two labels that it separates.
template <int D, typename Out_iter, typename Label_iter>
// void Connect(const std::vector<Vec<Constants::index_t, D> >& points,
void Connect(const std::vector<intn >& points,
             // const std::vector<Edge>& segments,
             const std::vector<LabeledSegment<D> >& segments,
             bool midpoint, Out_iter lines, Label_iter line_segments) {

  const int n = points.size();
  if (n == 2) {
    *lines++ = std::make_pair(points[0], points[1]);
    // We return only the segment that corresponds to the first point.
    *line_segments++ = segments[0];
    return;
  }

  if (midpoint) {
    // Connect at barycenter of all points
    const double nd = static_cast<double>(n);
    intn c = make_intn(0);
    for (int k = 0; k < n; ++k) {
      c += points[k]/nd;
    }
    for (int k = 0; k < n; ++k) {
      *lines++ = std::make_pair(points[k], c);
      *line_segments++ = segments[k];
    }
  } else {
    // find the closest two adjacent points
    int dist = INT_MAX;
    int mini = -1;
    for (int i = 0; i < n; ++i) {
      // const int d = (points[i]-points[(i+1)%n]).norm();
      const int d = length(points[i]-points[(i+1)%n]);
      if (d < dist) {
        dist = d;
        mini = i;
      }
    }
    // find the point furthest from mini
    dist = -1;
    int maxi = -1;
    for (int i = 0; i < n; ++i) {
      // const int d = (points[mini]-points[i]).norm();
      const int d = length(points[mini]-points[i]);
      if (d > dist) {
        dist = d;
        maxi = i;
      }
    }
    // branch at barycenter of mini, mini+1, maxi
    intn c = make_intn(0);
    c += points[mini] / 3;
    c += points[(mini+1)%n] / 3;
    c += points[maxi] / 3;
    *lines++ = std::make_pair(points[mini], c);
    *lines++ = std::make_pair(points[(mini+1)%n], c);
    *line_segments++ = segments[mini];
    *line_segments++ = segments[(mini+1)%n];
    if (n == 3) {
      *lines++ = std::make_pair(points[maxi], c);
      *line_segments++ = segments[maxi];
    }

    if (n > 3) {
      vector<intn> new_points;
      for (int i = 0; i < n; ++i) {
        if (i != mini && i != (mini+1)%n) {
          new_points.push_back(points[i]);
        }
        if (i == mini) {
          new_points.push_back(c);
        }
      }
      // Connect<D>(new_points, midpoint, lines);
      Connect<D>(new_points, segments, midpoint, lines, line_segments);
    }
  }
}

template <typename A, typename B>
struct tie_t {
  tie_t(A& a_, B& b_) : a(a_), b(b_) {}
  tie_t& operator=(const std::pair<A, B>& p) {
    a = p.first;
    b = p.second;
    return *this;
  }
  A& a;
  B& b;
};

template <typename A, typename B>
tie_t<A, B> tie(A& a, B& b) { return tie_t<A, B>(a, b); }


//------------------------------------------------------------
// ComputeQ
//
//   \\                   ||
//    |*                  ||
//   // \__               ||
//  //     \         _____*|
// //       a-------b     ||
//          |       |
//          |       |
//          |_______|
//
// * = pa, pb (closest points on curves to a and b)
// This function uses floating point computations rather
// than integer.
//
// Returns the computed point q and a distance measure from
// pa and pb.
//------------------------------------------------------------
template <typename Point>
std::pair<Point, float> ComputeQ(
    const Point& a, const Point& b,
    const Point& pa, const Point& pb, const bool error_msg) {

  using namespace std;
  static const int D = DIM;

  double alpha = 0;
  for (int i = 0; i < D; ++i) {
    alpha += (pb.s[i]*pb.s[i]-pa.s[i]*pa.s[i]);
  }

  bool bad = false;
  Point q = a;
  for (int i = 0; i < D; ++i) {
    if (a.s[i] != b.s[i]) {
      if (pa.s[i] == pb.s[i]) {
        // q.s[i] = pa.s[i];
        if (fabs(pa.s[i] - a.s[i]) < fabs(pa.s[i] - b.s[i]))
          q.s[i] = a.s[i];
        else
          q.s[i] = b.s[i];
      } else {
        // const double tf = (-(2*(a*(pa-pb)) + alpha)) / (2*(pa[i]-pb[i]));
        const double tf = (-(2*(dot(a, (pa-pb))) + alpha)) /
            (2*(pa.s[i]-pb.s[i]));
        //const index_t t = static_cast<index_t>(tf);
        q.s[i] += tf;
        // if (q[i] < std::min(a[i], b[i]) || q[i] > std::max(a[i], b[i])) {
        //   bad = true;
        // }
        // Clamp to endpoints if necessary.  This would be because the GVD
        // passes through a or b and there is numerical error in computing q.
        if (q.s[i] < std::min(a.s[i], b.s[i])) {
          if ((int)(std::min(a.s[i], b.s[i])-q.s[i]) > 1 && error_msg) {
            // More serious if the error is greater than 1.
            bad = true;
            cerr << "ERROR: bad computeq" << endl;
            cerr << "  i = " << i << endl;
            cerr << "  q = " << convert_intn(q) << endl;
            cerr << "  a = " << convert_intn(a) << " b = " << convert_intn(b) << endl;
            cerr << "  pa = " << convert_intn(pa) << " pb = " << convert_intn(pb) << endl;
          }
          q.s[i] = std::min(a.s[i], b.s[i]);
        } else if (q.s[i] > std::max(a.s[i], b.s[i])) {
          if ((int)(q.s[i]-std::max(a.s[i], b.s[i])) > 1 && error_msg) {
            // More serious if the error is greater than 1.
            bad = true;
            cerr << "ERROR: bad computeq" << endl;
            cerr << "  i = " << i << endl;
            cerr << "  q = " << convert_intn(q) << endl;
            cerr << "  a = " << convert_intn(a) << " b = " << convert_intn(b) << endl;
            cerr << "  pa = " << convert_intn(pa) << " pb = " << convert_intn(pb) << endl;
          }
          q.s[i] = std::max(a.s[i], b.s[i]);
        }
      }
    }
  }

  // double dist = (pa-q).norm() + (pb-q).norm();
  double dist = length(pa-q) + length(pb-q);
  if (bad) {
    dist = FLT_MAX;
  }
  return make_pair(q, dist);
}

template <int D>
std::pair<intn, double> ComputeIntersection(
    const int avi, const int bvi,
    const intn& a_, const intn& b_,
    const VertexNetwork& vertices, const OctreeOptions& o) {

  const int coi = o.cell_of_interest;

  const int alabel = vertices.Label(avi);
  const int blabel = vertices.Label(bvi);
  if (alabel == blabel) throw logic_error("Same labels");

  const doublen ac = convert_doublen(vertices.ClosestPoint(avi));
  const doublen bc = convert_doublen(vertices.ClosestPoint(bvi));

  const doublen a = convert_doublen(a_);
  const doublen b = convert_doublen(b_);
  doublen m;
  float d;
  tie(m, d) = oct::ComputeQ(a, b, ac, bc, avi==coi||bvi==coi);
  bool simple = o.simple_q;
  if (d < 0) {
    cerr << "Bad distance.  d = " << d << endl;
    cerr << "avi = " << avi << " bvi = " << bvi << endl;
    cerr << "a = " << a << " b = " << b << endl;
    cerr << "ac = " << ac << " bc = " << bc
         << endl;
    simple = true;
  }

  intn q = make_intn(0);
  if (simple) {
    q = convert_intn((a+b)/2);
    // d = ((a-ac).norm()+(b-bc).norm())/2;
    d = (length(a-ac)+length(b-bc))/2;
  } else {
    q = convert_intn(m);
  }

  return std::make_pair(q, (int)d);
}

//------------------------------------------------------------
// IntersectionCollector
//
// Edge visitor that computes intersections and stores them.
//------------------------------------------------------------
template <int D, typename Out_iter>
struct IntersectionCollector {
  IntersectionCollector(const VertexNetwork& verts, Out_iter points_iter,
                        const OctreeOptions& o_)
      : vertices(verts), points(points_iter), o(o_) {}

  bool operator()(int avi, int bvi,
                  intn a_, intn b_,
                  // oct::Direction<D> d) {
                  oct::Direction d) {
    using namespace std;
    // Ensure that q computations are done identically for every
    // edge since each one is computed twice.
    if (bvi < avi) {
      swap(avi, bvi);
      swap(a_, b_);
      // d = d.Opposite();  // this didn't compile for some reason
      // d = d.Reversed();
      d = Reversed(&d);
      // d = oct::Direction<D>(d.Neg(), d.Pos());
    }
    const int alabel = vertices.Label(avi);
    const int blabel = vertices.Label(bvi);
    if (alabel != blabel) {
      intn q = make_intn(0);
      int dist;
      tie(q, dist) = ComputeIntersection<D>(avi, bvi, a_, b_, vertices, o);

      // const doublen a = convert_doublen(a_);
      // const doublen b = convert_doublen(b_);
      LabeledVertex<D> a_labeled(a_, avi, vertices.Label(avi));
      LabeledVertex<D> b_labeled(b_, bvi, vertices.Label(bvi));
      *points++ = std::make_pair(q, LabeledSegment<D>(a_labeled, b_labeled));
    }
    return true;
  }
  const VertexNetwork& vertices;
  Out_iter points;
  OctreeOptions o;
};

template <int D, typename Out_iter>
void GetIntersectionsAroundFace(
    const int vi,
    const intn& base_point,
    const level_t level,
    const int primary_axis,
    const int secondary_axis,
    const VertexNetwork& vertices,
    Out_iter points,
    const OctreeOptions& o) {
  // const oct::Direction<D> prime_d =
  const oct::Direction prime_d =
      // oct::Direction<D>::FromAxis(primary_axis, true);
      oct::DirectionFromAxis(primary_axis, true);
  // const oct::Direction<D> sec_d =
  const oct::Direction sec_d =
      // oct::Direction<D>::FromAxis(secondary_axis, true);
      oct::DirectionFromAxis(secondary_axis, true);
  int v_id_p = vertices.FindNeighbor(vi, prime_d, level);
  int v_id_s = vertices.FindNeighbor(vi, sec_d, level);
  int v_id_ps = vertices.FindNeighbor(v_id_p, sec_d, level);
  int v_id_sp = vertices.FindNeighbor(v_id_s, prime_d, level);
  if (v_id_p == -1 || v_id_s == -1 || v_id_ps == -1 ||
      v_id_sp == -1 || v_id_ps != v_id_sp) {
    throw logic_error("Unexpected problem");
  }
  IntersectionCollector<D, Out_iter> v(vertices, points, o);
  WalkAroundFace<D>(vi, base_point, level,
                    primary_axis, secondary_axis,
                    vertices, v);
}
}

#endif
