#ifndef ERFC_H
#define ERFC_H

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
extern double b_erfc(double x);
static double rt_powd_snf(double u0, double u1);

#endif
