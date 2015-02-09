#include <iostream>
#include <cmath>

#include "candidate.h"


/* Builds the candidate between the two regions identified by their unique
 * IDs and their respective normals.
 */
Candidate::Candidate(int region1, int region2, float2 normal1, float2 normal2)
{
    this->region1 = region1;
    this->region2 = region2;
    computeScore(normal1, normal2);
}


/* Evaluates the score based on the angle between the two normals.
 * Assumes the vectors are normalized.
 * A lower score implies the vectors are closer together.
 */
void Candidate::computeScore(float2 normal1, float2 normal2)
{
    double dotp = normal1.x * normal2.x + normal1.y * normal2.y;
    score = acos(dotp);
    std::cout << score << std::endl;
}


/* Returns the score (angle between the two regions).
 */
double Candidate::getScore()
{
    return this->score;
}
