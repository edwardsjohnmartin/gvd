#include "float2.h"


class Candidate {

  private:
    int region1;
    int region2;
    double score;

    void computeScore(float2 normal1, float2 normal2);

  public:
    Candidate(int region1, int region2, float2 normal1, float2 normal2);

    bool operator < (const Candidate &other) const
    {
        return score < other.score;
    }

    double getScore();

};
