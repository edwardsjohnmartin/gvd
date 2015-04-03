#include <iostream>

#include "inverse_gauss_map.h"

using namespace std;


int main()
{
    InverseGaussMap<2> g1(4);
    InverseGaussMap<3> g2(4);

    float2 x;
    x.x = 1;
    x.y = 0;
    cout << g1.getBin(x) << endl;

    float3 y;
    y.x = 3;
    y.y = 2;
    y.z = 0;
    cout << g2.getBin(y) << endl;
}
