#ifndef __VEC_CPP_H__
#define __VEC_CPP_H__

#include <cmath>
#include <iostream>

#ifdef __OPEN_CL_SUPPORT__
//------------------------------------------------------------
// Use OpenCL's types
//------------------------------------------------------------
#define CL_VEC
#ifdef __APPLE__
#include "OpenCL/opencl.h"
#else
#include "CL/cl.h"
#endif
typedef cl_uchar2 bool2;
typedef cl_int2 int2;
typedef cl_float2 float2;
typedef cl_double2 double2;
typedef cl_uchar3 bool3;
typedef cl_int3 int3;
typedef cl_float3 float3;
typedef cl_double3 double3;
typedef cl_uchar4 bool4;
typedef cl_int4 int4;
typedef cl_float4 float4;
typedef cl_double4 double4;
#else
//------------------------------------------------------------
// Define our own type
//------------------------------------------------------------
template <class NumType, int NumDims>
struct MyVec {
  NumType s[NumDims];
};
// union MyVec {
//   struct { NumType s[NumDims]; };
//   struct { NumType x, y, z, w; };
// };

// #define __extension__
// typedef union {
//     int s[4];
// #if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
//    __extension__ struct{ cl_int  x, y, z, w; };
//    __extension__ struct{ cl_int  s0, s1, s2, s3; };
//    __extension__ struct{ cl_int2 lo, hi; };
// #endif
// } MyVec;

typedef MyVec<unsigned char, 2> bool2;
typedef MyVec<int, 2> int2;
typedef MyVec<float, 2> float2;
typedef MyVec<double, 2> double2;

typedef MyVec<int, 3> int3;
typedef MyVec<float, 3> float3;
typedef MyVec<double, 3> double3;
typedef MyVec<unsigned char, 3> bool3;

typedef MyVec<unsigned char, 4> bool4;
typedef MyVec<int, 4> int4;
typedef MyVec<float, 4> float4;
typedef MyVec<double, 4> double4;
#endif

//------------------------------------------------------------
// Define typen
//------------------------------------------------------------
#ifdef OCT2D
typedef bool2 booln;
typedef int2 intn;
typedef float2 floatn;
typedef double2 doublen;
#else
typedef bool3 booln;
typedef int3 intn;
typedef float3 floatn;
typedef double3 doublen;
#endif

//------------------------------------------------------------
// PointAndLabel struct
//------------------------------------------------------------
struct PointAndLabel {
  PointAndLabel(const intn& p_, const int l_) : p(p_), l(l_) {}
  intn p;
  int l;
};

//------------------------------------------------------------
// Functions
//------------------------------------------------------------

// Construction
#define MAKE_VECN(typen,type)                                         \
  inline typen make_##typen(type a=0, type b=0, type c=0, type d=0) {   \
    type arr[] = { a, b, c, d };                                        \
    return *(typen*)(arr);                                              \
    }
MAKE_VECN(bool2,bool)
MAKE_VECN(int2,int)
MAKE_VECN(float2,float)
MAKE_VECN(double2,double)
MAKE_VECN(bool3,bool)
MAKE_VECN(int3,int)
MAKE_VECN(float3,float)
MAKE_VECN(double3,double)
MAKE_VECN(bool4,bool)
MAKE_VECN(int4,int)
MAKE_VECN(float4,float)
MAKE_VECN(double4,double)

// #ifdef OCT2D
MAKE_VECN(booln,bool)
MAKE_VECN(intn,int)
MAKE_VECN(floatn,float)
MAKE_VECN(doublen,double)
// #else
// MAKE_VECN(booln,bool)
// MAKE_VECN(intn,int)
// MAKE_VECN(floatn,float)
// MAKE_VECN(doublen,double)
// #endif

#define MAKE_VECN_(typen,typen_,type,n)                         \
  inline typen make_##typen(const typen_& v_, type a=0) {       \
    type arr[] = { a, a, a, a };                                \
    for (int i = 0; i < n-1; ++i) {                             \
      arr[i] = v_.s[i];                                         \
    }                                                           \
    return *(typen*)(arr);                                      \
  }
MAKE_VECN_(bool3,bool2,bool,3)
MAKE_VECN_(int3,int2,int,3)
MAKE_VECN_(float3,float2,float,3)
MAKE_VECN_(double3,double2,double,3)
MAKE_VECN_(bool4,bool3,bool,4)
MAKE_VECN_(int4,int3,int,4)
MAKE_VECN_(float4,float3,float,4)
MAKE_VECN_(double4,double3,double,4)
// //------------------------------------------------------------
// // cast
// #define VEC_2_POINTER(typen,type)                           \
//   inline operator const type*(const typen& v) {      \
//     return &v.s[0];                                     \
//   }
// VEC_2_POINTER(int2,int)
// VEC_2_POINTER(float2,float)
// VEC_2_POINTER(double2,double)
// VEC_2_POINTER(int3,int)
// VEC_2_POINTER(float3,float)
// VEC_2_POINTER(double3,double)

//------------------------------------------------------------
// output operators
#define OUT_VEC(typen,n)                                        \
  inline std::ostream& operator<<(std::ostream& out, const typen& v) { \
    for (int i = 0; i < n; i++) {                               \
      out << v.s[i];                                            \
      if (i < n-1) out << " ";                                  \
    }                                                           \
    return out;                                                 \
  }
OUT_VEC(bool2,2)
OUT_VEC(int2,2)
OUT_VEC(float2,2)
OUT_VEC(double2,2)
#ifndef CL_VEC
OUT_VEC(bool3,3)
OUT_VEC(int3,3)
OUT_VEC(float3,3)
OUT_VEC(double3,3)
#endif
OUT_VEC(bool4,4)
OUT_VEC(int4,4)
OUT_VEC(float4,4)
OUT_VEC(double4,4)
#define IN_VEC(typen,n)                                         \
  inline std::istream& operator>>(std::istream& in, typen& v) {        \
    for (int i = 0; i < n; i++) {                               \
      in >> v.s[i];                                             \
    }                                                           \
    return in;                                                  \
  }
IN_VEC(bool2,2)
IN_VEC(int2,2)
IN_VEC(float2,2)
IN_VEC(double2,2)
#ifndef CL_VEC
IN_VEC(bool3,3)
IN_VEC(int3,3)
IN_VEC(float3,3)
IN_VEC(double3,3)
#endif
IN_VEC(bool4,4)
IN_VEC(int4,4)
IN_VEC(float4,4)
IN_VEC(double4,4)

//------------------------------------------------------------
// min/max
#define VEC_VEC_MIN(typen,n)                                    \
  inline typen vec_min(const typen& a, const typen& b) {                \
    typen result = make_##typen(); \
    for (int i = 0; i < n; i++) {                               \
      result.s[i] = (a.s[i] < b.s[i] ? a.s[i] : b.s[i]);        \
    }                                                           \
    return result;                                              \
  }
VEC_VEC_MIN(int2,2)
VEC_VEC_MIN(float2,2)
VEC_VEC_MIN(double2,2)
VEC_VEC_MIN(int3,3)
VEC_VEC_MIN(float3,3)
VEC_VEC_MIN(double3,3)
#define VEC_VEC_MAX(typen,n)                                    \
  inline typen vec_max(const typen& a, const typen& b) {            \
    typen result = make_##typen();                                    \
    for (int i = 0; i < n; i++) {                               \
      result.s[i] = (a.s[i] > b.s[i] ? a.s[i] : b.s[i]);        \
    }                                                           \
    return result;                                              \
  }
VEC_VEC_MAX(int2,2)
VEC_VEC_MAX(float2,2)
VEC_VEC_MAX(double2,2)
VEC_VEC_MAX(int3,3)
VEC_VEC_MAX(float3,3)
VEC_VEC_MAX(double3,3)

//------------------------------------------------------------
// Binary operators
// All binary operators work component-wise on vectors
//------------------------------------------------------------

#define VEC_OP_SCALAR(typen,type,n,op)                          \
  inline typen operator op (const typen v, const type a) {      \
  typen result = make_##typen();                                \
    for (int i = 0; i < n; i++) result.s[i] = v.s[i] op a;      \
    return result;                                              \
  }
#define VEC_OP_VEC(typen,n,op)                                  \
  inline typen operator op (const typen u, const typen v) {     \
    typen result = make_##typen();                                               \
    for (int i = 0; i < n; i++) result.s[i] = u.s[i] op v.s[i]; \
    return result;                                              \
  }
#define VEC_OP_ASSIGN_SCALAR(typen,type,n,op)           \
  inline typen& operator op (typen& v, const type a) {  \
    for (int i = 0; i < n; i++) v.s[i] op a;            \
    return v;                                           \
  }
#define VEC_OP_ASSIGN_VEC(typen,n,op)                   \
  inline typen& operator op (typen& u, const typen v) { \
    for (int i = 0; i < n; i++) u.s[i] op v.s[i];       \
    return u;                                           \
  }
// v * a
VEC_OP_SCALAR(int2, int, 2, *)
VEC_OP_SCALAR(float2, float, 2, *)
VEC_OP_SCALAR(double2, double, 2, *)
VEC_OP_SCALAR(int3, int, 3, *)
VEC_OP_SCALAR(float3, float, 3, *)
VEC_OP_SCALAR(double3, double, 3, *)
// v *= a
VEC_OP_ASSIGN_SCALAR(int2, int, 2, *=)
VEC_OP_ASSIGN_SCALAR(float2, float, 2, *=)
VEC_OP_ASSIGN_SCALAR(double2, double, 2, *=)
VEC_OP_ASSIGN_SCALAR(int3, int, 3, *=)
VEC_OP_ASSIGN_SCALAR(float3, float, 3, *=)
VEC_OP_ASSIGN_SCALAR(double3, double, 3, *=)
// u * v
VEC_OP_VEC(int2, 2, *)
VEC_OP_VEC(float2, 2, *)
VEC_OP_VEC(double2, 2, *)
VEC_OP_VEC(int3, 3, *)
VEC_OP_VEC(float3, 3, *)
VEC_OP_VEC(double3, 3, *)
// v / a
VEC_OP_SCALAR(int2, int, 2, /)
VEC_OP_SCALAR(float2, float, 2, /)
VEC_OP_SCALAR(double2, double, 2, /)
VEC_OP_SCALAR(int3, int, 3, /)
VEC_OP_SCALAR(float3, float, 3, /)
VEC_OP_SCALAR(double3, double, 3, /)
// v /= a
VEC_OP_ASSIGN_SCALAR(int2, int, 2, /=)
VEC_OP_ASSIGN_SCALAR(float2, float, 2, /=)
VEC_OP_ASSIGN_SCALAR(double2, double, 2, /=)
VEC_OP_ASSIGN_SCALAR(int3, int, 3, /=)
VEC_OP_ASSIGN_SCALAR(float3, float, 3, /=)
VEC_OP_ASSIGN_SCALAR(double3, double, 3, /=)
// v + a
VEC_OP_SCALAR(double2, double, 2, +)
VEC_OP_SCALAR(float2, float, 2, +)
VEC_OP_SCALAR(int2, int, 2, +)
VEC_OP_SCALAR(double3, double, 3, +)
VEC_OP_SCALAR(float3, float, 3, +)
VEC_OP_SCALAR(int3, int, 3, +)
// v += a
VEC_OP_ASSIGN_SCALAR(int2, int, 2, +=)
VEC_OP_ASSIGN_SCALAR(float2, float, 2, +=)
VEC_OP_ASSIGN_SCALAR(double2, double, 2, +=)
VEC_OP_ASSIGN_SCALAR(int3, int, 3, +=)
VEC_OP_ASSIGN_SCALAR(float3, float, 3, +=)
VEC_OP_ASSIGN_SCALAR(double3, double, 3, +=)
// u + v
VEC_OP_VEC(double2, 2, +)
VEC_OP_VEC(float2, 2, +)
VEC_OP_VEC(int2, 2, +)
VEC_OP_VEC(double3, 3, +)
VEC_OP_VEC(float3, 3, +)
VEC_OP_VEC(int3, 3, +)
// u += v
VEC_OP_ASSIGN_VEC(double2, 2, +=)
VEC_OP_ASSIGN_VEC(float2, 2, +=)
VEC_OP_ASSIGN_VEC(int2, 2, +=)
VEC_OP_ASSIGN_VEC(double3, 3, +=)
VEC_OP_ASSIGN_VEC(float3, 3, +=)
VEC_OP_ASSIGN_VEC(int3, 3, +=)
// v - a
VEC_OP_SCALAR(double2, double, 2, -)
VEC_OP_SCALAR(float2, float, 2, -)
VEC_OP_SCALAR(int2, int, 2, -)
VEC_OP_SCALAR(double3, double, 3, -)
VEC_OP_SCALAR(float3, float, 3, -)
VEC_OP_SCALAR(int3, int, 3, -)
// u - v
VEC_OP_VEC(double2, 2, -)
VEC_OP_VEC(float2, 2, -)
VEC_OP_VEC(int2, 2, -)
VEC_OP_VEC(double3, 3, -)
VEC_OP_VEC(float3, 3, -)
VEC_OP_VEC(int3, 3, -)
// - (negate)
#define VEC_NEGATE(typen,n)                             \
  inline typen operator-(const typen v) {               \
    typen result = make_##typen();                                       \
    for (int i = 0; i < n; i++) result.s[i] = -v.s[i];  \
    return result;                                      \
  }
VEC_NEGATE(double2, 2)
VEC_NEGATE(float2, 2)
VEC_NEGATE(int2, 2)
VEC_NEGATE(double3, 3)
VEC_NEGATE(float3, 3)
VEC_NEGATE(int3, 3)

//------------------------------------------------------------
// equality
#define VEC_EQUALS(typen,n)                                     \
  inline bool operator==(const typen& a, const typen& b) {      \
    bool eq = true;                                             \
    for (int i = 0; i < n && (eq = (a.s[i] == b.s[i])); ++i);  \
    return eq;                                                  \
  }
VEC_EQUALS(int3,3)
VEC_EQUALS(float3,3)
VEC_EQUALS(double3,3)
VEC_EQUALS(int2,2)
VEC_EQUALS(float2,2)
VEC_EQUALS(double2,2)
#define VEC_NOT_EQUALS(typen,n)                                 \
  inline bool operator!=(const typen& a, const typen& b) {      \
    return !(a == b);                                           \
  }
VEC_NOT_EQUALS(int3,3)
VEC_NOT_EQUALS(float3,3)
VEC_NOT_EQUALS(double3,3)
VEC_NOT_EQUALS(int2,2)
VEC_NOT_EQUALS(float2,2)
VEC_NOT_EQUALS(double2,2)
// <
#define VEC_LESS_SCALAR(typen,type,n)                           \
  inline bool operator<(const typen& v, const type& a) {        \
    bool lt = true;                                             \
    for (int i = 0; i < n && (lt = (v.s[i] < a)); i++);         \
    return lt;                                                  \
  }
VEC_LESS_SCALAR(int3,int,3)
VEC_LESS_SCALAR(float3,float,3)
VEC_LESS_SCALAR(double3,double,3)
VEC_LESS_SCALAR(int2,int,2)
VEC_LESS_SCALAR(float2,float,2)
VEC_LESS_SCALAR(double2,double,2)
#define VEC_LESS_VEC(typen,n)                                   \
  inline bool operator<(const typen& u, const typen& v) {       \
    for (int i = 0; i < n; ++i) {                               \
      if (u.s[i] < v.s[i]) return true;                         \
      if (u.s[i] > v.s[i]) return false;                        \
    }                                                           \
    return false;                                               \
  }
VEC_LESS_VEC(int3,3)
VEC_LESS_VEC(float3,3)
VEC_LESS_VEC(double3,3)
VEC_LESS_VEC(int2,2)
VEC_LESS_VEC(float2,2)
VEC_LESS_VEC(double2,2)

//------------------------------------------------------------
// Type conversions
#define CONVERT_VEC3(from,to,type)                              \
  inline to convert_##to(const from v) {                        \
    return make_##to((type)v.s[0], (type)v.s[1], (type)v.s[2]); \
  }
#define CONVERT_VEC2(from,to,type)                      \
  inline to convert_##to(const from v) {                \
    return make_##to((type)v.s[0], (type)v.s[1]);       \
  }

#ifdef OCT2D
#define CONVERT_VECN(from,to,type)                      \
  inline to convert_##to(const from v) {                \
    return make_##to((type)v.s[0], (type)v.s[1]);       \
  }
#else
#define CONVERT_VECN(from,to,type)                              \
  inline to convert_##to(const from v) {                        \
    return make_##to((type)v.s[0], (type)v.s[1], (type)v.s[2]); \
  }
#endif

CONVERT_VEC3(float3,int3,int)
CONVERT_VEC3(double3,int3,int)
CONVERT_VEC3(int3,float3,float)
CONVERT_VEC3(double3,float3,float)
CONVERT_VEC3(int3,double3,double)
CONVERT_VEC3(float3,double3,double)

CONVERT_VEC2(float2,int2,int)
CONVERT_VEC2(double2,int2,int)
CONVERT_VEC2(int2,float2,float)
CONVERT_VEC2(double2,float2,float)
CONVERT_VEC2(int2,double2,double)
CONVERT_VEC2(float2,double2,double)

CONVERT_VECN(floatn,intn,int)
CONVERT_VECN(doublen,intn,int)
CONVERT_VECN(intn,floatn,float)
CONVERT_VECN(doublen,floatn,float)
CONVERT_VECN(intn,doublen,double)
CONVERT_VECN(floatn,doublen,double)

//------------------------------------------------------------
// L1norm
inline int L1norm(const int2 a) {
  int result = 0;
  for (int i = 0; i < 2; i++) result += a.s[i];
  return result;
}
inline int L1norm(const int3 a) {
  int result = 0;
  for (int i = 0; i < 3; i++) result += a.s[i];
  return result;
}
//------------------------------------------------------------
// dot
inline double dot(const double3 a, const double3 b) {
  double result = 0;
  for (int i = 0; i < 3; i++) result += a.s[i]*b.s[i];
  return result;
}
inline float dot(const float3 a, const float3 b) {
  float result = 0;
  for (int i = 0; i < 3; i++) result += a.s[i]*b.s[i];
  return result;
}
inline double dot(const double2 a, const double2 b) {
  double result = 0;
  for (int i = 0; i < 2; i++) result += a.s[i]*b.s[i];
  return result;
}
inline float dot(const float2 a, const float2 b) {
  float result = 0;
  for (int i = 0; i < 2; i++) result += a.s[i]*b.s[i];
  return result;
}
//------------------------------------------------------------
// length2
inline double length2(const double3 v) {
  return dot(v, v);
}
inline float length2(const float3 v) {
  return dot(v, v);
}
inline int length2(const int3 v) {
  throw std::logic_error("can't get norm2 of integer type");
}
inline double length2(const double2 v) {
  return dot(v, v);
}
inline float length2(const float2 v) {
  return dot(v, v);
}
inline int length2(const int2 v) {
  throw std::logic_error("can't get norm2 of integer type");
}
//------------------------------------------------------------
// length
#define fast_length length
inline double length(const double3 v) {
  return sqrt(length2(v));
}
inline float length(const float3 v) {
  return sqrt(length2(v));
}
inline int length(const int3 v) {
  float3 vf = make_float3();
  for (int i = 0; i < 3; ++i) {
    vf.s[i] = v.s[i];
  }
  float n = length(vf);
  return (int)(n+0.5);
}
inline double length(const double2 v) {
  return sqrt(length2(v));
}
inline float length(const float2 v) {
  return sqrt(length2(v));
}
inline int length(const int2 v) {
  float2 vf = make_float2();
  for (int i = 0; i < 2; ++i) {
    vf.s[i] = v.s[i];
  }
  float n = length(vf);
  return (int)(n+0.5);
}
//------------------------------------------------------------
// normalize
#define fast_normalize normalize
inline double3 normalize(const double3 v) {
  return v / length(v);
}
inline float3 normalize(const float3 v) {
  return v / length(v);
}
inline int3 normalize(const int3 v) {
  return v / length(v);
}
inline double2 normalize(const double2 v) {
  return v / length(v);
}
inline float2 normalize(const float2 v) {
  return v / length(v);
}
inline int2 normalize(const int2 v) {
  return v / length(v);
}
//------------------------------------------------------------
// cross
#define VEC_CROSS_VEC(type3)                            \
  inline type3 cross(const type3 a, const type3 b) {    \
    type3 result = make_##type3();                      \
    result.s[0] = a.s[1]*b.s[2]-a.s[2]*b.s[1];          \
    result.s[1] = a.s[2]*b.s[0]-a.s[0]*b.s[2];          \
    result.s[2] = a.s[0]*b.s[1]-a.s[1]*b.s[0];          \
    return result;                                      \
  }
VEC_CROSS_VEC(int3)
VEC_CROSS_VEC(float3)
VEC_CROSS_VEC(double3)

#endif
