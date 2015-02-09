#include <iostream>
#include <set>

#include "candidate.h"

using namespace std;


int main()
{
    set<Candidate> s;
    float2 a;
    float2 b;
    float2 c;
    a.x = 0; a.y = 1;
    b.x = 0.5; b.y = 0.866;
    c.x = 1; c.y = 0;
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
