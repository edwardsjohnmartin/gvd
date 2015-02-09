/*
 */

#ifndef REGION_H
#define REGION_H


class Region {

  private:
    float2 min_norm;
    float2 max_norm;

  public:
    Region(float2 normal) {
        min_norm = normal;
        max_norm = normal;
    }

    void merge(Region *other) {
        // TODO
        // if other.max_norm > max_norm, max_norm = other.max_norm
        // if other.min_norm < min_norm, min_norm = other.min_norm
    }

    float2 getNormal() {
        // TODO - left or right? average??
        return min_norm;
    }
};


#endif
