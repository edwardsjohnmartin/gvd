#ifndef __GEOM_H__
#define __GEOM_H__

#include <vector>

#include "./segment.h"

inline int sgn(const float& f) {
  return (f < 0) ? -1 : ((f > 0) ? 1 : 0);
}

floatn centroid(const float_seg& s);

inline float dist(const floatn& a, const floatn& b) {
  return length(a-b);
}

// Finds the point on segment s that is closest to p.
floatn closest(const floatn& p, const float_seg& s);

// Finds the point on segment s that intersects with the horizontal line at y
// Returns true if such a point exists.
bool x_given_y(const float_seg& s, const float& y, float* x);

// Finds the point on segment s that intersects with the vertical line at x.
// Returns true if such a point exists.
bool y_given_x(const float_seg& s, const float& x, float* y);

// s - segment to be clipped
// d - edge length of box
// returns false if s is entirely contained within the box, true otherwise
// Precondition: s.a() is the center of the box
bool clip_to_box(const float_seg& s, const float& d, float_seg* new_s);

void closest(const float_seg& a, const float_seg& b,
             float2* ca, float2* cb,
             bool* ca_end, bool* cb_end);

void v_segs(const float_seg& a, const float_seg& b,
            float_seg* aa_, float_seg* bb_);

std::vector<floatn> v_sample(const float_seg& a, const float_seg& b);

// Pass A and B by value! They'll change in the function.
void FitBoxes(float_seg A, float_seg B, const float& min_d,
              std::vector<floatn>* origins, std::vector<float>* lengths);

bool multi_intersection(const float_seg& A, const float_seg& B);

floatn line_intersection(
    const floatn& p1, const floatn& p2, const floatn& p3, const floatn& p4,
    bool* lines_intersect, bool* segs_intersect);

#endif
