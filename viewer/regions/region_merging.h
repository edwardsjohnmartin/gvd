/* Region.
 */

#include <algorithm>
#include <set>
#include <vector>
#include <iostream>

#include "region.h"
#include "candidate.h"


//#define THRESHOLD 1.571 // half pi (quarter circle)
#define THRESHOLD 0.7854 // quarter pi (either circle)


void mergeRegions(std::vector<float2> normals)
{
    // build the regions table
    const int num_parts = normals.size();
    Region * regions[num_parts];
    for(int i=0; i<num_parts; i++)
        regions[i] = new Region(normals[i], i);

    // set up the priority queue with all region neighbors
    std::set<Candidate> queue;
    for(int i=0; i<num_parts - 1; i++)
        queue.insert(Candidate(regions, i, i+1));
    if(num_parts > 1)
        queue.insert(Candidate(regions, num_parts-1, 0));

    // housekeeping: keep track of indices we can delete later
    std::vector<int> merged_indices;

    // now run the algorithm
    std::set<Candidate>::iterator it;
    while(!queue.empty())
    {
        it = queue.begin();
        Candidate best = *(it);
        // if best.region1 in merged or best.region2 in merged,
        //  continue
        if(best.getScore() < THRESHOLD)
        {
            regions[best.region1]->merge(regions[best.region2]);
            regions[best.region2] = regions[best.region1];
            merged_indices.push_back(best.region2);
            std::cout << "Merged " << best.region2 << " with "
                      << best.region1 << std::endl;
            // merged = new Region(region1, region2, next id);
            // for all neighbors of merged:
            //  queue.insert(Candidate(merged, neighbor);
            // seen.push_back(region1.id);
            // seen.push_back(region2.id);
        }
        else {
            std::cout << "Failed: " << best.getScore() << std::endl;
            break;
        }
        queue.erase(it);
    }

    // remove the regions that were merged to other regions
    std::vector<Region *> remaining_regions;
    remaining_regions.assign(regions, regions + num_parts);
    std::sort(remaining_regions.begin(), remaining_regions.end());
    for(int offset=0; offset<merged_indices.size(); offset++)
    {
        int remove_index = merged_indices[offset] - offset;
        remaining_regions.erase(remaining_regions.begin() + remove_index);
    }

    std::cout << "Num regions remaining: " << remaining_regions.size() << std::endl;
    for(int i=0; i<remaining_regions.size(); i++)
    {
        std::cout << "idx = " << remaining_regions[i]->getId() << std::endl;
    }
}
