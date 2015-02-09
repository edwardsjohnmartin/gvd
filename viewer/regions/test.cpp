#include <iostream>
#include <set>

#include "candidate.h"

using namespace std;


int main()
{
    float2 a = make_float2(0, 1);
    float2 b = make_float2(0.5, 0.866);
    float2 c = make_float2(1, 0);
    set<Candidate> s;
    s.insert(Candidate(1, 2, b, c));
    s.insert(Candidate(0, 1, a, b));
    s.insert(Candidate(0, 2, a, c));
    set<Candidate>::iterator it = s.begin();
    Candidate x = *(it);
    s.erase(it);
    cout << "Top score: " << x.getScore() << endl;
    it = s.begin();
    x = *(it);
    s.erase(it);
    cout << "Next best score: " << x.getScore() << endl;
    return 0;
}
