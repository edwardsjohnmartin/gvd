/* Region.
 */

#include <set>
#include <vector>
#include <iostream>

#include "candidate.h"


void mergeRegions(std::vector<float2> verts)
{
    std::set<Candidate> queue;
    for(int i=0; i<verts.size() - 1; i++)
        queue.insert(Candidate(i, i+1, verts[i], verts[i+1]));
    int last_vert_idx = verts.size() - 1;
    if(last_vert_idx > 0)
        queue.insert(Candidate(last_vert_idx, 0, verts[last_vert_idx], verts[0]));
    
    std::set<Candidate>::iterator it;
    while(!queue.empty())
    {
        it = queue.begin();
        Candidate c = *(it);
        std::cout << "Score: " << c.getScore() << std::endl;
        // TODO - if region is within angle, merge them
        queue.erase(it);
    }
}
