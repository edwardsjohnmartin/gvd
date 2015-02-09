/* Region.
 */

#include <set>
#include <vector>
#include <iostream>

#include "region.h"
#include "candidate.h"


#define THRESHOLD 1.0


void mergeRegions(std::vector<float2> normals)
{
    // build the regions table
    const int num_regions = normals.size();
    Region * regions[num_regions];
    for(int i=0; i<num_regions; i++)
        regions[i] = new Region(normals[i]);

    // set up the priority queue with all region neighbors
    std::set<Candidate> queue;
    for(int i=0; i<num_regions - 1; i++)
        queue.insert(Candidate(regions, i, i+1));
    if(num_regions > 1)
        queue.insert(Candidate(regions, num_regions-1, 0));

    // now run the algorithm
    std::set<Candidate>::iterator it;
    while(!queue.empty())
    {
        it = queue.begin();
        Candidate best = *(it);
        if(best.getScore() < THRESHOLD)
        {
            regions[best.region1]->merge(regions[best.region2]);
            regions[best.region2] = regions[best.region1];
        }
        else
        {
            break;
        }
        queue.erase(it);
    }

    // now collect unique regions from the map...

    /*
    // set up the queue with all region neighbors
    std::set<Candidate> queue;
    for(int i=0; i<verts.size() - 1; i++)
        queue.insert(Candidate(i, i+1, verts[i], verts[i+1]));
    int last_idx = verts.size() - 1;
    if(last_idx > 0)
        queue.insert(Candidate(last_idx, 0, verts[last_idx], verts[0]));
    
    // run the merging algorithm
    std::set<Candidate>::iterator it;
    while(!queue.empty())
    {
        it = queue.begin();
        Candidate best = *(it);
        std::cout << "Score: " << best.getScore() << std::endl;
        // TODO - if region is within angle, merge them
        queue.erase(it);
    }
    */
}
