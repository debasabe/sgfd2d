#ifndef MACROSHDR
#define MACROSHDR

#define _USE_MATH_DEFINES
#include <cmath>
#include "cmatrix.h"

// Physical Domain (1km x 1km)
#define XMIN 0.0
#define XMAX 1000.0
#define ZMIN 0.0
#define ZMAX 1000.0
// Time domain
#define TMAX 1.5
// Source parameters
#define PKFREQ 50.0
#define XSRC 500.0
#define ZSRC 500.0
#define NSX 0.0
#define NSZ 1.0
#define XSW 10.0
#define ZSW 10.0
#define VMAX 2000.0
// Interval between saves
#define SAVEINT 200

#endif
