// --------------------------------------------------------------------------
//                           SSEARCH HEADER FILE
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
//  Release date  : AUG-2019
//  Authors       : Hernández Olcina, Jorge
//
//  Master in Geomatics Engineering and Geoinformation
//  Universidad Politécnica de Valencia
//  Hochschule Karlsruhe - Technik und Wirtschaft University of Applied Sciences
// --------------------------------------------------------------------------
#ifndef SSEARCH_H
#define SSEARCH_H

// Include Files
#include <cmath>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "..\math_functions\rt_nonfinite.h"
#include "..\math_functions\rtwtypes.h"
#include "..\LAMBDA_types.h"

// --------------------------------------------------------------------------
//                          Function Declarations
// --------------------------------------------------------------------------
extern void b_ssearch(const double ahat_data[], const int ahat_size[1], const
                      double L_data[], const int L_size[2], const double D_data[],
                      double ncands, emxArray_real_T *afixed, emxArray_real_T
                      *sqnorm);
extern void ssearch(const double ahat[12], const double L[144], const double D
                    [12], double ncands, emxArray_real_T *afixed,
                    emxArray_real_T *sqnorm);

#endif
