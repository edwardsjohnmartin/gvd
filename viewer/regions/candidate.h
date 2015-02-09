/* The candidate class for region merging is used for the priority queue
 * to score regions based on their angle similarity.
 */

#ifndef CANDIDATE_H
#define CANDIDATE_H


#include "../../opencl/vec_cpp.h"


class Candidate {

  private:

    // Unique IDs of the regions that this candidate pair encompasses.
    int region1;
    int region2;

    // Score (angle) between the two region normals.
    double score;

    // Helper function to compute the angle (and score) between two normals.
    void computeScore(float2 normal1, float2 normal2);


  public:

    // Initialize the regions by IDs and computes their score by their normals.
    Candidate(int region1, int region2, float2 normal1, float2 normal2);

    // Compares the score of this object with the other object.
    // Returns true if this score is lower than the other score.
    bool operator < (const Candidate &other) const
    {
        return score < other.score;
    }

    // Returns the current score of this candidate pair.
    double getScore();

};


#endif
