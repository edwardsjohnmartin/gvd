#ifndef __CL_AMBIGUOUS_H__
#define __CL_AMBIGUOUS_H__

#include "./defs.h"
#include "./vertex.h"

NAMESPACE_OCT_BEGIN

int IsAmbiguous3_orig(
    const int base_vi,
    global const Vertex* vertices,
    global const int* point_labels);

int IsAmbiguous3(
    const int base_vi,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    __GLOBAL__ GeomPoint* cpoints_array);

NAMESPACE_OCT_END

#endif
