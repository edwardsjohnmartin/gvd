/* The candidate class for region merging is used for the priority queue
 * to score regions based on their angle similarity.
 */

#ifndef CANDIDATE_H
#define CANDIDATE_H


#include "../../opencl/vec_cpp.h"
#include "region.h"


class Candidate {

  public:
    int region1;
    int region2;

  private:
    Region ** regions;

  public:

    // Initialize the regions by IDs and computes their score by their normals.
    Candidate(Region ** regions, int region1, int region2);

    // Returns the current score of this candidate pair.
    double getScore() const;

    // Compares the score of this object with the other object.
    // Returns true if this score is lower than the other score.
    bool operator < (const Candidate &other) const
    {
        return this->getScore() < other.getScore();
    }

};


#endif
