#ifndef XZGETRF_H
#define XZGETRF_H

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
extern void xzgetrf(int n, double A[], int ipiv[], int *info);

#endif
