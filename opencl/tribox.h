#ifndef __TRI_BOX_H__
#define __TRI_BOX_H__

#include "./defs.h"
#include "./triangle.h"
#include "./vec.h"

// int triBoxOverlap(const float boxcenter[3],
//                   const float boxhalfsize[3],
//                   const float triverts[3][3]);
//                   // const float** triverts);

bool TriBoxOverlap(const int3 boxcenter,
                   const int boxhalfsize,
                   const __GLOBAL__ int3* vertices,
                   const Triangle tri);

#endif
