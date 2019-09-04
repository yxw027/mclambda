// --------------------------------------------------------------------------
//                          MC-LAMBDA HEADER FILE
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
//  Release date  : AUG-2019
//  Authors       : Hernández Olcina, Jorge
//
//  Master in Geomatics Engineering and Geoinformation
//  Universidad Politécnica de Valencia
//  Hochschule Karlsruhe - Technik und Wirtschaft University of Applied Sciences
// --------------------------------------------------------------------------
#ifndef LAMBDA_H
#define LAMBDA_H

// Include Files
#include <cmath>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "math_functions\rt_nonfinite.h"
#include "math_functions\rtwtypes.h"
#include "mclambda_types.h"

// --------------------------------------------------------------------------
//                          Function Declarations
// --------------------------------------------------------------------------
static double rt_remd_snf(double u0, double u1);

extern void mclambda(int n, double ahat[], const double Qahat[], double method, double
            param, const emxArray_char_T *type, double value, emxArray_real_T
            *afixed, emxArray_real_T *sqnorm, double *Ps, double Qzhat[],
            double Z[], double *nfixed, double *mu);

#endif
