#include <iostream>

#include "inverse_gauss_map.h"

using namespace std;


void test() {
    InverseGaussMap<3> g2(16);
    float3 v;
v.x = -0.200825; v.y = 1.034244; v.z = 0.003236;
    cout << g2.getBin(v) << endl;
v.x = -0.388418; v.y = 0.977338; v.z = 0.003236;
    cout << g2.getBin(v) << endl;
v.x = -0.561304; v.y = 0.884928; v.z = 0.003236;
    cout << g2.getBin(v) << endl;
v.x = -0.712841; v.y = 0.760565; v.z = 0.003236;
    cout << g2.getBin(v) << endl;
v.x = -0.837204; v.y = 0.609029; v.z = 0.003236;
    cout << g2.getBin(v) << endl;
v.x = -0.929614; v.y = 0.436142; v.z = 0.003236;
    cout << g2.getBin(v) << endl;
v.x = -0.986519; v.y = 0.248549; v.z = 0.003236;
    cout << g2.getBin(v) << endl;
v.x = -1.005734; v.y = 0.053459; v.z = 0.003236;
    cout << g2.getBin(v) << endl;
v.x = -0.197076; v.y = 1.034244; v.z = -0.034824;
    cout << g2.getBin(v) << endl;
v.x = -0.381065; v.y = 0.977338; v.z = -0.071422;
    cout << g2.getBin(v) << endl;
v.x = -0.550629; v.y = 0.884928; v.z = -0.105151;
    cout << g2.getBin(v) << endl;
v.x = -0.699254; v.y = 0.760565; v.z = -0.134714;
    cout << g2.getBin(v) << endl;
v.x = -0.821227; v.y = 0.609029; v.z = -0.158976;
    cout << g2.getBin(v) << endl;
v.x = -0.911862; v.y = 0.436142; v.z = -0.177004;
    cout << g2.getBin(v) << endl;
v.x = -0.967674; v.y = 0.248549; v.z = -0.188106;
    cout << g2.getBin(v) << endl;
v.x = -0.986519; v.y = 0.053459; v.z = -0.191855;
    cout << g2.getBin(v) << endl;
}
int main()
{
    InverseGaussMap<2> g1(16);
    InverseGaussMap<3> g2(16);

    float2 x;
    x.x = 1;
    x.y = 0;
    //cout << g1.getBin(x) << endl;

    float3 y;
    y.x = -1;
    y.y = 0;
    y.z = 0;
    float norm = sqrt(y.x*y.x+y.y*y.y+y.z*y.z);
    y.x /= norm;
    y.y /= norm;
    y.z /= norm;
    cout << g2.getBin(y) << endl;
    test();
}
