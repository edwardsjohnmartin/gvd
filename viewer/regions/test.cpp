#include <iostream>
#include <set>
#include <vector>

#include "candidate.h"
#include "region_merging.h"

using namespace std;


void test_merge_algorithm()
{
    vector<float2> normals;
    normals.push_back(make_float2(0, 1));
    normals.push_back(make_float2(0.5, 0.866));
    normals.push_back(make_float2(1, 0));
    normals.push_back(make_float2(-1, 0));
    mergeRegions(normals);
}


void test_candidate()
{
    /*float2 a = make_float2(0, 1);
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
    cout << "Next best score: " << x.getScore() << endl;*/
}


int main()
{
    //test_candidate();
    test_merge_algorithm();
    return 0;
}
