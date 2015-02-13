#include <iostream>
#include <set>
#include <vector>

#include "candidate.h"
#include "region_merging.h"

using namespace std;


void test_merge_algorithm()
{
    vector<float2> normals;
    // sqrt(3)/2 = 0.866
    // sqrt(2)/2 = 0.707

    // R1 (0, 1, 2)
    normals.push_back(make_float2(1, 0));
    normals.push_back(make_float2(0.866, 0.5));
    normals.push_back(make_float2(0.707, 0.707));

    // R2 (3)
    normals.push_back(make_float2(-0.866, 0.5));

    // R3 (4, 5, 6)
    normals.push_back(make_float2(-0.5, -0.866));
    normals.push_back(make_float2(0, -1));
    normals.push_back(make_float2(0.5, -0.866));
    vector<Region> regions = mergeRegions(normals);
    for(vector<Region>::iterator it = regions.begin();
        it != regions.end(); it++)
    {
        vector<int> verts = it->getVerts();
        cout << "Region " << it->getId() << ": ";
        for(vector<int>::iterator it2 = verts.begin();
            it2 != verts.end(); it2++)
        {
            cout << (*it2) << ", ";
        }
        cout << endl;
    }
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
