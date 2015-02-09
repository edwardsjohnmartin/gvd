/*
 */

#ifndef REGION_H
#define REGION_H


class Region {

  private:
    float2 average_norm;
    int num_regions;
    int unique_id;

  public:
    Region(float2 normal, int id) {
        average_norm = normal;
        num_regions = 1;
        unique_id = id;
    }

    void merge(Region *other) {
        float2 weighted_average = getNormal() * numRegions();
        weighted_average += other->getNormal() * other->numRegions();
        num_regions += other->numRegions();
        average_norm = weighted_average / num_regions;
    }

    float2 getNormal() {
        return average_norm;
    }

    int numRegions() {
        return num_regions;
    }

    int getId() {
        return unique_id;
    }
};


#endif
