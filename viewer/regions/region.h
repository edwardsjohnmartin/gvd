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
    Region(float2 normal, int id, int num_parts = 1) {
        average_norm = normal;
        num_regions = num_parts;
        unique_id = id;
    }

    Region * merge(Region *other, int new_id) {
        float2 weighted_average = getNormal() * numRegions();
        weighted_average += other->getNormal() * other->numRegions();
        average_norm = weighted_average / num_regions;
        int num_parts = numRegions() + other->numRegions();
        return new Region(average_norm, new_id, num_parts);
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
