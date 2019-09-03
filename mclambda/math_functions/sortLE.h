#ifndef SORTLE_H
#define SORTLE_H

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
extern boolean_T sortLE(const emxArray_real_T *v, const int col_data[], const
  int col_size[2], int irow1, int irow2);

#endif
