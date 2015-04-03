#include <iostream>

#include "inverse_gauss_map.h"

using namespace std;


int main()
{
    InverseGaussMap g1(2, 8);
    InverseGaussMap g2(3, 8);

    float2 x;
    x.x = 2;
    x.y = 1;
    cout << g1.getBin<2>(x) << endl;

    float3 y;
    y.x = 3;
    y.y = 2;
    y.z = 0;
    cout << g2.getBin<3>(y) << endl;
}
