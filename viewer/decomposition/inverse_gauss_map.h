/**
 */

#ifndef INVERSE_GAUSS_MAP_H
#define INVERSE_GAUSS_MAP_H


#include <assert.h>
#include "../../opencl/vec_cpp.h"


class InverseGaussMap {

  private:
    int dimension;
    int resolution;

    int getBin2D(float x, float y);
    int getBin3D(float x, float y, float z);
    int getBinFromAngle(float angle, int k);


  public:
    InverseGaussMap(int dimension, int resolution);

    void setResolution(int resolution);

    /**
     * Returns the bin for the correct dimension.
     */
    template<int D> int getBin(MyVec<float, D> normal) {
        assert(D == dimension);
        if(dimension == 2)
            return getBin2D(normal.x, normal.y);
        else
            return getBin3D(normal.x, normal.y, normal.z);
    }

};


#endif
