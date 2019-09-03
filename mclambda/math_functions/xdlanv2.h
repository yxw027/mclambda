#ifndef XDLANV2_H
#define XDLANV2_H

// Include Files
#include <cmath>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "..\mclambda_types.h"

// Function Declarations
extern void xdlanv2(double *a, double *b, double *c, double *d, double *rt1r,
                    double *rt1i, double *rt2r, double *rt2i, double *cs, double
                    *sn);

#endif
