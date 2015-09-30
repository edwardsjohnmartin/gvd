#include "./geom.h"

using std::vector;
using std::logic_error;
using std::swap;
using std::cout;
using std::endl;

floatn centroid(const float_seg& s) {
  return (s.a() + s.b()) / 2;
}

floatn closest(const floatn& p, const float_seg& s) {
  if (s.is_degenerate())
    return s.a();
  static const float EPSILON = 1e-6;
  floatn v = s.vector();
  float len = length(v);
  v = v / len;
  floatn u = p - s.a();
  float t = dot(u, v);
  floatn q = s.a() + v * t;
  if (t < 0 || dist(q, s.a()) < EPSILON)
    q = s.a();
  else if (t > len || dist(q, s.b()) < EPSILON)
    q = s.b();
  return q;
}

bool x_given_y(const float_seg& s, const float& y, float* x) {
  // p + tv
  // p.y + t*v.y = y
  // t = (y - p.y) / v.y
  // x = p.x + t*v.x
  const floatn p = s.p();
  const floatn v = s.vector();
  const float t = (y - p.y) / v.y;
  *x = p.x + t * v.x;
  return t <= 1;
}

bool y_given_x(const float_seg& s, const float& x, float* y) {
  // p + tv
  // p.x + t*v.x = x
  // t = (x - p.x) / v.x
  // y = p.y + t*v.y
  const floatn p = s.p();
  const floatn v = s.vector();
  const float t = (x - p.x) / v.x;
  *y = p.y + t * v.y;
  return t <= 1;
}

// Given a box and a line segment s such that s.a() is the center of the box:
//      _    _____   ____
//     |    |   __|_/
//   d |    |  /  |
//     |_   |_____|
//
//
// compute
//             _____   ____
//            |     |_/
//            |     |
//            |_____|
bool clip_to_box(const float_seg& s, const float& d, float_seg* new_s) {
  // c is the center of the box
  const floatn c = s.a();
  const float d2 = d/2;
  const floatn s_vec = s.vector();
  bool valid;
  float x, y;
  if (fabs(s_vec.x) > fabs(s_vec.y)) {
    //           ____
    //      ____/
    // ____/
    x = c.x + d2;
    if (s_vec.x < 0)
      x = c.x - d2;
    valid = y_given_x(s, x, &y);
  } else {
    //      _/
    //    _/
    //  _/
    // /
    y = c.y + d2;
    if (s_vec.y < 0)
      y = c.y - d2;
    valid = x_given_y(s, y, &x);
  }
  if (valid)
    *new_s = float_seg(make_floatn(x, y), s.b());
  return valid;
}

// Returns true if two line segments intersect at any more than one point.
//
//       p1
//         *_
//           \_
//             \_
//  p3 *---------\_-----------* p4
//                 \_
//                   \_
//                     \_* p2
//
bool multi_intersection(
    const floatn& p1, const floatn& p2, const floatn& p3, const floatn& p4) {
  // Store the values for fast access and easy
  // equations-to-code conversion
  const float x1 = p1.x, x2 = p2.x, x3 = p3.x, x4 = p4.x;
  const float y1 = p1.y, y2 = p2.y, y3 = p3.y, y4 = p4.y;
 
  const float d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
  if (d == 0) {
    // Lines are parallel (determinant is zero)
    if (x2-x1 != 0) {
      // Lines are not vertical
      const double m = (y2-y1)/(x2-x1);
      if (y1-y3 == m * (x1-x3)) {
        // Lines are collinear
        const float len1 = fabs(x2-x1);
        const float len2 = fabs(x4-x3);
        const float len =
            fmax(x1, fmax(x2, fmax(x3, x4))) - fmin(x1, fmin(x2, fmin(x3, x4)));
        if (len < len1 + len2) {
          return true;
        }
      }
    } else {
      // Lines are vertical
      if (x1 == x3) {
        // Lines are collinear
        const float len1 = fabs(y2-y1);
        const float len2 = fabs(y4-y3);
        const float len =
            fmax(y1, fmax(y2, fmax(y3, y4))) - fmin(y1, fmin(y2, fmin(y3, y4)));
        if (len < len1 + len2) {
          return true;
        }
      }
    }
  }
  return false;
}

bool multi_intersection(const float_seg& A, const float_seg& B) {
  return multi_intersection(A.a(), A.b(), B.a(), B.b());
}

// line-line intersection
//
//       p1
//         *_
//           \_
//             \_
//  p3 *---------\_-----------* p4
//                 \_
//                   \_
//                     \_* p2
//
// segs_intersect is true if the line segments intersect. False otherwise.
floatn line_intersection(
    const floatn& p1, const floatn& p2, const floatn& p3, const floatn& p4,
    bool* lines_intersect, bool* segs_intersect) {

  // TODO: remove dependence on double precision. Can we do everything with
  // integer arithmetic?

  // Store the values for fast access and easy
  // equations-to-code conversion
  const double x1 = p1.x, x2 = p2.x, x3 = p3.x, x4 = p4.x;
  const double y1 = p1.y, y2 = p2.y, y3 = p3.y, y4 = p4.y;
 
  const double d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
  // If d is zero, there is no intersection
  if (d == 0) {
    *lines_intersect = false;
    *segs_intersect = false;
    return make_floatn(0);
  }
  *lines_intersect = true;
 
  // Get the x and y
  const double pre = (x1*y2 - y1*x2), post = (x3*y4 - y3*x4);
  const double x = ( pre * (x3 - x4) - (x1 - x2) * post ) / d;
  const double y = ( pre * (y3 - y4) - (y1 - y2) * post ) / d;
 
  // cout << x1 << " " << y2 << " " << y1 << " " << x2 << endl;
  // cout << "p1 = " << p1 << " p2 = " << p2 << " p3 = " << p3 << " p4 = " << p4
  //      << endl;
  // cout << "d = " << d << " pre = " << pre << " post = " << post << endl;

  // Check if the x and y coordinates are within both lines
  *segs_intersect = true;
  if ( x < fmin(x1, x2) || x > fmax(x1, x2) ||
       x < fmin(x3, x4) || x > fmax(x3, x4) ) {
    *segs_intersect = false;
  }
  if ( y < fmin(y1, y2) || y > fmax(y1, y2) ||
       y < fmin(y3, y4) || y > fmax(y3, y4) ) {
    *segs_intersect = false;
  }
 
  return make_floatn(x, y);
}

// floatn line_intersection(
//     const floatn& p1, const floatn& p2, const floatn& p3, const floatn& p4,
//     bool* lines_intersect, bool* segs_intersect) {
//   const float x1 = p1.x;
//   const float x2 = p2.x;
//   const float y1 = p1.y;
//   const float y2 = p2.y;
//   float A2 = p4.y-p3.y;//other.Y2() - other.Y1();
//   float B2 = p4.x-p3.x;//other.X2() - other.X1();
//   float C2 = A2*p3.x + B2*p3.y;//A2*other.X1() + B2*other.Y1();

//   float A1 = y2 - y1;
//   float B1 = x2 - x1;
//   float C1 = A1 * x1 + B1 * y1;

//   float det = A1*B2 - A2*B1;
//   if (det == 0) {
//     *lines_intersect = false;
//     *segs_intersect = false;
//     return make_floatn(0);
//   }
//   *lines_intersect = true;

//   cout << "d = " << det << endl;

//   const float x = (B2 * C1 - B1 * C2) / det;
//   const float y = (A1 * C2 - A2 * C1) / det;
//   return make_floatn(x, y);
// }

floatn line_intersection(const float_seg& a, const float_seg& b,
                         bool* lines_intersect, bool* segs_intersect) {
  return line_intersection(a.a(), a.b(), b.a(), b.b(),
                           lines_intersect, segs_intersect);
}

// Given two segments
//
//             a  ___*
//            ___/
//           /
//          *            
//  *------------------------*
//             b
//
// find the two closest points
//
//             a  ___*
//            ___/ 
//           / 
//          x            
//  *-------x-----------------*
//           b
//
// Also works for intersecting line segments.
void closest(const float_seg& a, const float_seg& b,
            floatn* ca, floatn* cb,
             bool* ca_end, bool* cb_end) {
  // Check for intersection
  bool lines_intersect, segs_intersect;
  const floatn l_intersection = line_intersection(
      a, b, &lines_intersect, &segs_intersect);
  if (segs_intersect) {
    *ca = *cb = l_intersection;
    *ca_end = *cb_end = false;
  } else {
    // No intersection
    floatn c = closest(a.a(), b);
    *ca = a.a();
    *cb = c;
    *ca_end = true;
    *cb_end = (c == b.a() || c == b.b());
    float dist = length(a.a()-c);

    floatn c_ = closest(a.b(), b);
    float dist_ = length(a.b()-c_);
    if (dist_ < dist) {
      c = c_;
      *ca = a.b();
      *cb = c;
      *ca_end = true;
      *cb_end = (c == b.a() || c == b.b());
      dist = dist_;
    }

    c_ = closest(b.a(), a);
    dist_ = length(b.a()-c_);
    if (dist_ < dist) {
      c = c_;
      *ca = c;
      *cb = b.a();
      *ca_end = (c == a.a() || c == a.b());
      *cb_end = true;
      dist = dist_;
    }

    c_ = closest(b.b(), a);
    dist_ = length(b.b()-c_);
    if (dist_ < dist) {
      c = c_;
      *ca = c;
      *cb = b.b();
      *ca_end = (c == a.a() || c == a.b());
      *cb_end = true;
      dist = dist_;
    }

    // Now consider this case (and the corresponding case where endpoints
    // have identical y values):
    //
    //                 ____/
    //            ____/
    //       ____/
    //      /
    //     *
    //     *____
    //          \____
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        const floatn& a0 = a[i];
        const floatn& a1 = a[(i+1)%2];
        const floatn& b0 = b[j];
        const floatn& b1 = b[(j+1)%2];
        if (a0.x == b0.x) {
          const int dir_ay = sgn(a1.y - a0.y);
          const int dir_by = sgn(b1.y - b0.y);
          const int dir_aby = sgn(b0.y - a0.y);
          if (dir_ay != dir_by && dir_ay != dir_aby) {
            *ca = a[i];
            *cb = b[j];
            *ca_end = *cb_end = true;
          }
        } else if (a0.y == b0.y) {
          const int dir_ax = sgn(a1.x - a0.x);
          const int dir_bx = sgn(b1.x - b0.x);
          const int dir_abx = sgn(b0.x - a0.x);
          if (dir_ax != dir_bx && dir_ax != dir_abx) {
            *ca = a[i];
            *cb = b[j];
            *ca_end = *cb_end = true;
          }
        }
      }
    }
  }
}

// Given two segments
//
//             a  ___*
//            ___/
//           /
//          *            
//  *------------------------*
//             b
//
// compute three new segments
//
//             aa ___*
//            ___/ 
//           / --> 
//          *            
//  *-------*-----------------*
//   <-- bb0      bb1 -->
//
// and choose the bb that has a positive dot product with aa
//
//             aa ___*
//            ___/ 
//           / --> 
//          *            
//          *-----------------*
//                bb -->
//
void v_segs(const float_seg& a, const float_seg& b,
            float_seg* aa_, float_seg* bb_) {
  floatn c = closest(a.a(), b);
  float dist = length(a.a()-c);
  float_seg aa = a;
  float_seg bb0(c, b.a());
  float_seg bb1(c, b.b());

  floatn c_ = closest(a.b(), b);
  float dist_ = length(a.b()-c_);
  if (dist_ < dist) {
    c = c_;
    dist = dist_;
    aa = float_seg(a.b(), a.a());
    bb0 = float_seg(c, b.a());
    bb1 = float_seg(c, b.b());
  }

  c_ = closest(b.a(), a);
  dist_ = length(b.a()-c_);
  if (dist_ < dist) {
    c = c_;
    dist = dist_;
    aa = float_seg(b.a(), b.b());
    bb0 = float_seg(c, a.a());
    bb1 = float_seg(c, a.b());
  }

  c_ = closest(b.b(), a);
  dist_ = length(b.b()-c_);
  if (dist_ < dist) {
    c = c_;
    dist = dist_;
    aa = float_seg(b.b(), b.a());
    bb0 = float_seg(c, a.a());
    bb1 = float_seg(c, a.b());
  }

  float_seg bb = bb0;
  if (dot(aa.vector(), bb1.vector()) > 0)
    bb = bb1;

  *aa_ = aa;
  *bb_ = bb;
}

// Given two segments
//
//             b  ___*
//            ___/
//           /
//          *            
//  *___
//      \______
//             \______
//          a          \___*
//
// orient and swap as necessary so that the point closest to the other
//
//             a  ___*
//            ___/ 
//           / --> 
//          *            
//  *___      c
//      \___*__
//             \______
//          b-->      \___*
//
// Point c is the closest point on b to segment a.
//
//             a  ___*
//            ___/ |
//           / --> |
//          *      | d'    
//  *___  d |      |
//      \___*______|____________
//             \___|__
//          b-->      \___*
//
//
vector<floatn> v_sample(const float_seg& a, const float_seg& b) {
  vector<floatn> samples;
  return samples;
}

// Pass A and B by value! They'll change in the function.
void FitBoxesNoIntersection(float_seg A, float_seg B, const float& min_d,
                            std::vector<floatn>* samples,
                            vector<floatn>* origins, vector<float>* lengths) {
  static const int OFFSET_FACTOR = 4;

  const float_seg A_orig(A);
  const float_seg B_orig(B);
  bool done = false;
  while (!done) {
    // Closest points between A and B
    floatn ca, cb;
    // Closest points between B and C
    floatn c_cb;
    // Whether A's (resp. B's) closest point is an endpoint of the segment
    bool ca_end, cb_end;
    // Whether B's (resp. C's) closest point is an endpoint of the segment
    bool c_cb_end;

    closest(A, B, &ca, &cb, &ca_end, &cb_end);

    if (ca == cb) {
      throw logic_error("FixBoxes doesn't support intersections");
    }

    if (!ca_end && !cb_end) throw logic_error("no end but not intersecting");
    // Make A the segment with the end intersection
    if (cb_end) {
      swap(ca, cb);
      swap(ca_end, cb_end);
      swap(A, B);
    }
    if (!ca_end) {
      throw logic_error("Closest points between A and B must "
                        "include one endpoint");
    }
    // Reverse A if necessary so that A.a() == ca
    if (A.a() != ca) {
      A.reverse();
    }
    // If cb is an endpoint of B, reverse B if necessary so that B.a() == cb
    if (cb_end && B.a() != cb) {
      B.reverse();
    }

    float d;
    floatn o;
    // floatn new_samples[2];
    floatn sample_offset = make_floatn(0, 0);
    if (dist(ca, cb) < min_d) {
      o = (ca+cb)/2 - make_uni_floatn(min_d/2);
      d = min_d;
      c_cb = cb;
      c_cb_end = cb_end;
      // new_samples[0] = ca;
      // new_samples[1] = cb;
    } else if (cb_end && (A.a().x == B.a().x || A.a().y == B.a().y)) {
      // Special case: closest points are axis aligned
      d = dist(A.a(), B.a());
      o = make_floatn(fmin(A.a().x, B.a().x), fmin(A.a().y, B.a().y));
      if (A.a().x == B.a().x) {
        sample_offset = make_floatn(d/OFFSET_FACTOR, 0);
        if (A.vector().x < 0) {
          o = o - make_floatn(d, 0);
          sample_offset = sample_offset * -1;
        }
      } else if (A.a().y == B.a().y) {
        sample_offset = make_floatn(0, d/16);
        if (A.vector().y < 0) {
          o = o - make_floatn(0, d);
          sample_offset = sample_offset * -1;
        }
      }
      c_cb = cb;
      c_cb_end = true;
    } else {
      // Given the vector v between closest points ca and cb, compute the
      // best fit for a square between the closest points. The square doesn't
      // necessarily intersect the closest points. Intersecting points are
      // ca and c_cb (which intersects with line segment B).
      const floatn v = cb - ca;
      d = length(v);
      const float cx_dir = (v.x >= 0) ? 2*d : -2*d;
      const float cy_dir = (v.y >= 0) ? 2*d : -2*d;
      const float_seg C(ca, ca + make_floatn(cx_dir, cy_dir));
      floatn c_cc;
      bool c_cc_end;
      closest(B, C, &c_cb, &c_cc, &c_cb_end, &c_cc_end);
      float d_ = dist(ca, c_cb) / (float)sqrt(2);
      if (c_cb != c_cc) {
        // There's no intersection between B and C
        d_ = fmax(fabs(c_cb.x-ca.x), fabs(c_cb.y-ca.y));
      }
      if (d_ < d) {
        d = d_;
      }
      const float box_x_dir = (v.x >= 0) ? 0 : -d;
      const float box_y_dir = (v.y >= 0) ? 0 : -d;
      o = ca + make_floatn(box_x_dir, box_y_dir);
    }
    // Add samples between ca and cb to avoid numerical problems when
    // ca or cb are on cell boundary
    samples->push_back(ca + (cb-ca)/OFFSET_FACTOR + sample_offset);
    samples->push_back(cb + (ca-cb)/OFFSET_FACTOR + sample_offset);
    origins->push_back(o);
    lengths->push_back(d);

    // If c_cb is not an endpoint, split B into two line segments
    // and discard the one with negative dot product with A
    if (c_cb_end) {
      if (B.a() != c_cb) {
        B.reverse();
      }
    } else {
      const float_seg B0(c_cb, B.a());
      const float_seg B1(c_cb, B.b());
      if (dot(A.unit(), B0.unit()) > 0) {
        B = B0;
      } else {
        B = B1;
      }
    }

    if (dot(A.unit(), B.unit()) <= 0) {
      done = true;
    } else {
      const floatn bisect = unit(A.unit()+B.unit());
      bool horizontal;
      if (A.a().x == B.a().x) {
        horizontal = true;
      } else if (A.a().y == B.a().y) {
        horizontal = false;
      } else {
        horizontal = (fabs(bisect.x) > fabs(bisect.y));
      }

      float A_x, A_y, B_x, B_y;
      bool A_intersects, B_intersects;
      if (horizontal) {
        if (bisect.x < 0) {
          A_x = B_x = o.x;
        } else {
          A_x = B_x = o.x + d;
        }
        A_intersects = y_given_x(A, A_x, &A_y);
        B_intersects = y_given_x(B, B_x, &B_y);
      } else {
        if (bisect.y < 0) {
          A_y = B_y = o.y;
        } else {
          A_y = B_y = o.y + d;
        }
        A_intersects = x_given_y(A, A_y, &A_x);
        B_intersects = x_given_y(B, B_y, &B_x);
      }
      done = (!A_intersects || !B_intersects);
      if (A.a() == make_floatn(A_x, A_y) && B.a() == make_floatn(B_x, B_y)) {
        done = true;
        cout << "Infinite loop" << endl;
        cout << "A_orig = " << A_orig << " B_orig = " << B_orig << endl;
      }
      A = float_seg(make_floatn(A_x, A_y), A.b());
      B = float_seg(make_floatn(B_x, B_y), B.b());
    }
  }
}

// Pass A and B by value! They'll change in the function.
void FitBoxes(float_seg A, float_seg B, const float& min_d,
              std::vector<floatn>* samples,
              vector<floatn>* origins, vector<float>* lengths) {
  // First check to see if there's an intersection. If so, split the
  // segments into multiple segments.
  floatn ca, cb;
  bool ca_end, cb_end;
  closest(A, B, &ca, &cb, &ca_end, &cb_end);
  //if (dist(ca, cb) < min_d / 2) {
  if (ca == cb) {
    // Intersection or very close approach
    samples->push_back(ca);
    origins->push_back(ca - make_floatn(min_d/2, min_d/2));
    lengths->push_back(min_d);
    // Reorder A and B so they are moving left to right
    if (A.a().x > A.b().x)
      A.reverse();
    if (B.a().x > B.b().x)
      B.reverse();
    float_seg A0_orig(ca, A.a());
    float_seg A1_orig(ca, A.b());
    float_seg B0_orig(cb, B.a());
    float_seg B1_orig(cb, B.b());
    float_seg A0, A1, B0, B1;
    bool A0_valid = clip_to_box(A0_orig, min_d, &A0);
    bool A1_valid = clip_to_box(A1_orig, min_d, &A1);
    bool B0_valid = clip_to_box(B0_orig, min_d, &B0);
    bool B1_valid = clip_to_box(B1_orig, min_d, &B1);
    const float dotp = dot(A0.unit(), B0.unit());
    if (dotp < 0) {
      swap(A0, A1);
      swap(A0_valid, A1_valid);
    }
    if (A0_valid && B0_valid)
      FitBoxesNoIntersection(A0, B0, min_d, samples, origins, lengths);
    if (A1_valid && B1_valid)
      FitBoxesNoIntersection(A1, B1, min_d, samples, origins, lengths);
  } else {
    FitBoxesNoIntersection(A, B, min_d, samples, origins, lengths);
  }
}
