/* The candidate class for region merging is used for the priority queue
 * to score regions based on their angle similarity.
 */

#ifndef CANDIDATE_H
#define CANDIDATE_H


struct Candidate {

    int region1;
    int region2;
    double score;

    // Constructor: assigns the two regions and their score (angle)..
    Candidate(int region1, int region2, double score)
       : region1(region1), region2(region2), score(score) { }

    // Compares the score of this object with the other object.
    // Returns true if this score is lower than the other score.
    bool operator < (const Candidate &other) const
    {
        return score < other.score;
    }

};


#endif
