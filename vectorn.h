#ifdef OCT2D
#include "./vector2.h"
namespace oct {
typedef LabeledGeometry2 LabeledGeometry;
}

#else
#include "./vector3.h"
namespace oct {
typedef LabeledGeometry3 LabeledGeometry;
}
#endif

