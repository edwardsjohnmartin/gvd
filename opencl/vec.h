#ifndef _VEC_H_
#define _VEC_H_

#ifdef OPEN_CL
#include "./vec_cl.h"
#else
#include "./vec_cpp.h"
#endif

inline int intn_comp(int i, intn p) {
  if (i == 0)
    return p.x;
  if (i == 1)
    return p.y;
#ifdef OCT3D
  if (i == 2)
    return p.z;
#endif
  return -1;
}

inline void set_intn_comp(int i, intn* p, int value) {
  if (i == 0)
    p->x = value;
  else if (i == 1)
    p->y = value;
#ifdef OCT3D
  else if (i == 2)
    p->z = value;
#endif
}

#endif
