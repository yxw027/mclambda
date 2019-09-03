#ifndef XZLASCL_H
#define XZLASCL_H

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
extern void b_xzlascl(int n, double cfrom, double cto, creal_T A[]);
extern void xzlascl(int n, double cfrom, double cto, creal_T A[]);

#endif
