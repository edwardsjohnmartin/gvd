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
#include <iostream>

#include "../../opencl/vec_cpp.h"


/** class InverseGaussMap
 * Defines the inverse Gauss map object with constant-time bin functions and
 * updatable resolution. The dimension of this object must be specified (2
 * or 3 dimensions only). If D is invalid, the object will throw an assertion
 * error during construction.
 */
template<int DIMENSION>
class InverseGaussMap
{

  private:

    // The resolution
    int resolution;
    int resolution_sqrt;

    /**
     * Returns the bin number for a 2D circular Gauss map.
     */
    int getBin2D(float x, float y) const
    {
        assert(DIMENSION == 2);
        float theta = atan2(y, x);
        if (theta < 0)
            theta += 2*M_PI;
        return getBinFromAngle(theta, resolution);
    }

    /**
     * Returns the bin number for a 3D spherical Gauss map.
     */
    int getBin3D(float x, float y, float z) const
    {
        assert(DIMENSION == 3);
        // TODO - algorithm for finding angles might be wrong?
        float theta = atan2(y, x);
        float phi = atan2(z, sqrt(x*x + y*y));
        int x_bin = getBinFromAngle(theta, resolution_sqrt);
        int y_bin = getBinFromAngle(phi, resolution_sqrt);
        return (y_bin * resolution_sqrt) + x_bin;
    }

    /**
     * Returns the integer binning of the angle in between 0 and k (exclusive),
     * provided the angle is between 0 and 2 Pi. If angle is larger than 2 Pi
     * or smaller than 0, the resulting bin will be converted to the correct
     * value. Assumes k is correctly larger than 0.
     */
    int getBinFromAngle(float theta, int k) const
    {
        int bin = (int)floor((theta / (2*M_PI)) * k);
        bin = bin % k;
        if(bin < 0)
            bin += k;
        return bin;
    }


  public:

    /**
     * Initialize the Gauss map in 2D or 3D (circle or sphere, respectively),
     * and set up the initial resolution. Dimension must be 2 or 3, and
     * resolution must be greater than 0.
     */
    InverseGaussMap(int resolution = 4)
    {
        assert(DIMENSION == 2 || DIMENSION == 3);
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
        if(DIMENSION == 3)
        {
            resolution_sqrt = (int)ceil(pow(resolution, 0.5));
            resolution = resolution_sqrt * resolution_sqrt;
        }
        this->resolution = resolution;
    }

    /**
     * Returns the resulution (number of bins).
     */
    int getResolution() const
    {
        return resolution;
    }

    /**
     * Returns the bin for the correct dimension by calling the appropriate
     * 2D or 3D bin function.
     */
    int getBin(MyVec<float, DIMENSION> normal) const
    {
        if(DIMENSION == 2)
            return getBin2D(normal.x, normal.y);
        else
            return getBin3D(normal.x, normal.y, normal.z);
    }

};


#endif
