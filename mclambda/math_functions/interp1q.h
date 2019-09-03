#ifndef INTERP1Q_H
#define INTERP1Q_H

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
extern double interp1q(const double x[31], const double y[31], double xi);

#endif
