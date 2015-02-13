/* Region.
 */

#include <algorithm>
#include <set>
#include <unordered_set>
#include <vector>
#include <iostream>

#include "region.h"
#include "candidate.h"


//#define THRESHOLD 1.571 // half pi (quarter circle)
#define THRESHOLD 0.7854 // quarter pi (either circle)


std::vector<Region> mergeRegions(std::vector<float2> normals)
{
    const int num_verts = normals.size();

    // Create a list of regions, and add a new region for each vertex.
    std::vector<Region> regions;
    regions.reserve(num_verts);
    for(int vert_id = 0; vert_id < num_verts; vert_id++)
    {
        regions.push_back(Region(normals[vert_id], vert_id));
        regions[vert_id].addVertex(vert_id);
    }

    // If only one or two vertices, done. Return the list as is.
    // TODO - necessary?
    if(num_verts < 2)
        return regions;
    
    // Build the queue of candidates.
    std::set<Candidate> queue;
    for(int i = 0; i < regions.size()-1; i++)
    {
        queue.insert(Candidate(regions[i].getId(), regions[i+1].getId(),
                               regions[i].angleTo(regions[i+1])));
    }
    queue.insert(Candidate(regions[regions.size()-1].getId(),
                           regions[0].getId(),
                           regions[regions.size()-1].angleTo(regions[0])));

    std::unordered_set<int> merged;
    int next_id = regions.size();

    // Run the merging algorithm loop.
    std::set<Candidate>::iterator it;
    while(!queue.empty())
    {
        // Pop the top of the queue.
        it = queue.begin();
        Candidate best = *(it);
        queue.erase(it);

        // If either region of this candidate was already merged, skip it.
        if(merged.find(best.region1) != merged.end() ||
           merged.find(best.region2) != merged.end())
        {
            continue;
        }

        // If the best candidate's angle (score) exceeds the threshold, stop.
        if(best.score > THRESHOLD)
            break;
    }

    return regions;
}
