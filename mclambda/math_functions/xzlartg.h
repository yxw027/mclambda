#ifndef XZLARTG_H
#define XZLARTG_H

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
extern void b_xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn);
extern void xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn,
                    creal_T *r);

#endif
