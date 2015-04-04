/**
 * The InverseGaussMap provides binning of arbitrary size (this value can be
 * changed) that groups angles defined by normal vectors in 2D or 3D. This
 * object provides constant-time lookups and resolution updates to the Gauss
 * map, and only needs the 2D or 3D normal vector to return a mapping to the
 * circle or sphere.
 *
 * The spherical Gauss map is based on Lambert mapping.
 */

#ifndef INVERSE_GAUSS_MAP_H
#define INVERSE_GAUSS_MAP_H


#include <assert.h>
#include <math.h>

#include "../../opencl/vec_cpp.h"
#include <iostream>


/** class InverseGaussMap
 * Defines the inverse Gauss map object with constant-time bin functions and
 * updatable resolution. The dimension of this object must be specified (2
 * or 3 dimensions only). If D is invalid, the object will throw an assertion
 * error during construction.
 */
template<int DIM> class InverseGaussMap
{

  private:

    // The resolution
    int resolution;
    int resolution_sqrt;

    /**
     * Returns the bin number for a 2D circular Gauss map.
     */
    int getBin2D(float x, float y)
    {
        assert(DIM == 2);
        float theta = acos(x);
        if(y < 0)
            theta += M_PI;
        return getBinFromAngle(theta, resolution);
    }

    float getAngle(float x, float y)
    {
        if(x == 0 && y == 0)
            return 0;
        float norm = sqrt(x*x + y*y);
        x = x / norm;
        y = y / norm;
        float theta = acos(x);
        if(y < 0)
            theta += M_PI;
        return theta;
    }

    /**
     * Returns the bin number for a 3D spherical Gauss map.
     */
    int getBin3D(float x, float y, float z)
    {
        assert(DIM == 3);
        float theta = getAngle(x, y);
        float phi = getAngle(y, z);
        theta = acos(z);
        phi = atan2(y, x);
        int x_bin = getBinFromAngle(theta, resolution_sqrt);
        int y_bin = getBinFromAngle(phi, resolution_sqrt);
        std::cout << x_bin << " / " << y_bin << std::endl;
        return (y_bin * resolution_sqrt) + x_bin;
    }

    /**
     * Returns the integer binning of the angle in between 0 and k (exclusive),
     * provided the angle is between 0 and 2 Pi. If angle is larger than 2 Pi
     * or smaller than 0, the resulting bin will be converted to the correct
     * value. Assumes k is correctly larger than 0.
     */
    int getBinFromAngle(float theta, int k)
    {
        int bin = (int)floor((theta / (2*M_PI)) * k);
        bin = bin % k;
        if(bin < 0)
            bin += 4;
        return bin;
    }


  public:

    /**
     * Initialize the Gauss map in 2D or 3D (circle or sphere, respectively),
     * and set up the initial resolution. Dimension must be 2 or 3, and
     * resolution must be greater than 0.
     */
    InverseGaussMap(int resolution)
    {
        assert(DIM == 2 || DIM == 3);
        setResolution(resolution);
    }

    /**
     * Set the resolution of the Gauss map. Resolution must be greater than 0.
     * In 3D, resolution must be a perfect square, and will be set to the next
     * largest perfect square if such a value is not provided.
     */
    void setResolution(int resolution)
    {
        assert(resolution > 0);
        if(DIM == 3)
        {
            resolution_sqrt = (int)ceil(pow(resolution, 0.5));
            resolution = resolution_sqrt * resolution_sqrt;
        }
        this->resolution = resolution;
    }

    /**
     * Returns the bin for the correct dimension by calling the appropriate
     * 2D or 3D bin function.
     */
    int getBin(MyVec<float, DIM> normal) {
        if(DIM == 2)
            return getBin2D(normal.x, normal.y);
        else
            return getBin3D(normal.x, normal.y, normal.z);
    }

};


#endif
