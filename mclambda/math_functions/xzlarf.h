#ifndef XZLARF_H
#define XZLARF_H

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
extern void b_xzlarf(int n1, int m, int n, int iv0, double tau, double C[], int ic0, double
              work[]);
extern void xzlarf(int n1, int n, int iv0, double tau, double C[], int ic0, double
                   work[]);

#endif
