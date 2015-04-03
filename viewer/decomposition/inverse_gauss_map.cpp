#include "inverse_gauss_map.h"

#include <assert.h>
#include <math.h>

#include <iostream>
using namespace std;


/**
 * Initialize the Gauss map in 2D or 3D (circle or sphere, respectively), and
 * set up the initial resolution. Dimension must be 2 or 3, and resolution
 * must be greater than 0.
 */
InverseGaussMap::InverseGaussMap(int dimension, int resolution)
    : dimension(dimension)
{
    assert(dimension == 2 || dimension == 3);
    setResolution(resolution);
}


/**
 * Set the resolution of the Gauss map. Resolution must be greater than 0. In
 * 3D, resolution must be a perfect square, and will be set to the next largest
 * perfect square if such a value is not provided.
 */
void InverseGaussMap::setResolution(int resolution)
{
    assert(resolution > 0);
    if(dimension == 3)
    {
        int sqroot = (int)ceil(pow(resolution, 0.5));
        resolution = sqroot * sqroot;
    }
    this->resolution = resolution;
}


/**
 * Returns the bin number (depending on 2D or 3D binning).
 */
int InverseGaussMap::getBin(int normal)
{
    return 0;
}


int InverseGaussMap::getBin2D(int normal)
{
    return 0;
}


/**
 * Returns the integer binning of the angle in between 0 and k (exclusive),
 * provided the angle is between 0 and 2 Pi. If angle is larger than 2 Pi or
 * smaller than 0, the resulting bin will be converted to the correct value.
 */
int InverseGaussMap::getBinFromAngle(float theta, int k)
{
    int bin = (int)floor((theta / (2*M_PI)) * k);
    bin = bin % k;
    if(bin < 0)
        bin += 4;
    cout << theta << " -> " << bin << endl;
    return bin;
}
