#ifndef __GEOMETRY_CPP_H__
#define __GEOMETRY_CPP_H__

//------------------------------------------------------------
// Test code to be used only in c++
//------------------------------------------------------------

#include <vector>
#include "./opencl/geometry.h"
#include "./vectorn.h"

NAMESPACE_OCT_BEGIN

LabeledGeometry3 Convert(Geometry geometry);
std::vector<LabeledGeometry3> Convert(Geometries geometries);
std::vector<std::vector<LabeledGeometry3> > Convert(
    Vi2Geometries vi2geometries);
// shared_array<int> Convert(
//     const std::vector<std::vector<LabeledGeometry3> >& lg,
//     Vi2Geometries& vi2geometries);
shared_array<int> Convert(
    const int num_vertices,
    const std::vector<LabeledGeometry3>& geometries,
    Vi2Geometries& vi2geometries);
// void Insert(const std::vector<LabeledGeometry3>& old_geometries,
//                    const int vi, Vi2Geometries& vi2geometries);
std::ostream& operator<<(
    std::ostream& out, Vi2Geometries vi2geometries);
void Compare(
    const std::vector<std::vector<LabeledGeometry3> >& base2geometries,
    Vi2Geometries vi2geometries);

NAMESPACE_OCT_END

#endif
