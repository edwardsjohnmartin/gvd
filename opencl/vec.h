#ifndef _VEC_H_
#define _VEC_H_

#ifdef OPEN_CL
#include "./vec_cl.h"
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

#else
#include "./vec_cpp.h"
inline int intn_comp(int i, intn p) {
  if (i == 0)
    return p.s[0];
  if (i == 1)
    return p.s[1];
#ifdef OCT3D
  if (i == 2)
    return p.s[2];
#endif
  return -1;
}
inline void set_intn_comp(int i, intn* p, int value) {
  if (i == 0)
    p->s[0] = value;
  else if (i == 1)
    p->s[1] = value;
#ifdef OCT3D
  else if (i == 2)
    p->s[2] = value;
#endif
}

#endif


#endif
