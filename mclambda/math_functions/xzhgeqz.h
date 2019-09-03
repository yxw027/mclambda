#ifndef XZHGEQZ_H
#define XZHGEQZ_H

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
extern void xzhgeqz(int n, const creal_T A[], int ilo, int ihi, int *info, creal_T alpha1
             [], creal_T beta1[]);

#endif
