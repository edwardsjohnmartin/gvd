/* Region.
 */

#include <algorithm>
#include <set>
#include <unordered_set>
#include <vector>
#include <iostream>

#include "region.h"
#include "candidate.h"


#define THRESHOLD 1.571 // half pi (quarter circle)
//#define THRESHOLD 0.7854 // quarter pi (either circle)


std::vector<Region> mergeRegions(std::vector<float2> normals)
{
    const int num_verts = normals.size();

    // Create a list of regions, and add a new region for each vertex.
    std::vector<Region> regions;
    regions.reserve(num_verts);
    for(int vert_id = 0; vert_id < num_verts; vert_id++)
    {
        Region region = Region(normals[vert_id], vert_id);
        region.addVertex(vert_id);
        int left_neighbor = (vert_id > 0) ? (vert_id - 1) : (num_verts - 1);
        int right_neighbor = (vert_id < num_verts - 1) ? (vert_id + 1) : 0;
        region.setNeighbors(left_neighbor, right_neighbor);
        regions.push_back(region);
    }

    // If only one or two vertices, done. Return the list as is.
    // TODO - necessary?
    if(num_verts < 2)
        return regions;
    
    // Build the queue of candidates. Candidates are built with the assumption
    // that the regions are ordered left to right. This allows for more
    // efficient neighbor search later on.
    std::set<Candidate> queue;
    for(int i = 0; i < regions.size()-1; i++)
    {
        queue.insert(Candidate(regions[i].getId(), regions[i+1].getId(),
                               regions[i].angleTo(regions[i+1])));
    }
    queue.insert(Candidate(regions[regions.size()-1].getId(),
                           regions[0].getId(),
                           regions[regions.size()-1].angleTo(regions[0])));

    // Keep track of merged regions.
    std::unordered_set<int> merged;

    // Run the merging algorithm loop.
    while(!queue.empty())
    {
        // Pop the top of the queue.
        std::set<Candidate>::iterator top = queue.begin();
        Candidate best = *(top);
        queue.erase(top);

        // If either region of this candidate was already merged, skip it.
        if(merged.find(best.region1) != merged.end() ||
           merged.find(best.region2) != merged.end())
        {
            continue;
        }

        // If the best candidate's angle (score) exceeds the threshold, stop.
        if(best.score > THRESHOLD) { std::cout << "done: " << best.score << std::endl;
            break;}

        // If this candidate is valid, merge the two regions and add new
        // candidates for all of the region's neighbors.
        Region new_region = regions[best.region1].merge(regions[best.region2],
                                                        regions.size());
        std::cout << "Merged regions " << best.region1 << " and " << best.region2 << std::endl;

        int left_neighbor = regions[best.region1].getLeftNeighbor();
        int right_neighbor = regions[best.region2].getRightNeighbor();

        queue.insert(Candidate(left_neighbor, new_region.getId(),
                               regions[left_neighbor].angleTo(new_region)));
        queue.insert(Candidate(new_region.getId(), right_neighbor,
                               new_region.angleTo(regions[right_neighbor])));
        std::cout << "New scores: " << regions[left_neighbor].angleTo(new_region)
                  << ", " << new_region.angleTo(regions[right_neighbor]) << std::endl;

        new_region.setNeighbors(left_neighbor, right_neighbor);

        // Add the old regions to the merged pile, and add the new region.
        regions.push_back(new_region);
        merged.insert(best.region1);
        merged.insert(best.region2);
    }

    // Remove all regions that have been merged to other regions.
    std::vector<Region>::iterator region_itr = regions.begin();
    while(region_itr != regions.end())
    {
        if(merged.find(region_itr->getId()) != merged.end())
            region_itr = regions.erase(region_itr);
        else
            region_itr++;
    }

    return regions;
}
