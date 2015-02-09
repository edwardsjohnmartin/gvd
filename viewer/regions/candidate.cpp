#include <iostream>
#include <cmath>

#include "candidate.h"


/* Builds the candidate between the two regions identified by their unique
 * IDs and their respective normals.
 */
Candidate::Candidate(std::vector<Region *> &regions, int region1, int region2)
{
    this->regions = regions;
    this->region1 = region1;
    this->region2 = region2;
}


/* Computes the score (angle between the two regions).
 */
double Candidate::getScore() const
{
    float2 normal1 = regions[region1]->getNormal();
    float2 normal2 = regions[region2]->getNormal();
    double dotp = normal1.x * normal2.x + normal1.y * normal2.y;
    double score = acos(dotp);
    //std::cout << score << std::endl;
    return score;
    return 0;
}
