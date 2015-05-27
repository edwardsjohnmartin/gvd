#ifndef __KARRAS_H__
#define __KARRAS_H__

#include <vector>

#include "./opencl/vec.h"

int xyz2z(intn p);
intn z2xyz(const int z);

void ktest(const std::vector<floatn>& points);
void ktest(const std::vector<intn>& points);


#endif
