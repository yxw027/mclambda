// --------------------------------------------------------------------------
//                        MC-LAMBDA rtwutil HEADER FILE
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
//  Release date  : AUG-2019
//  Authors       : Hernández Olcina, Jorge
//
//  Master in Geomatics Engineering and Geoinformation
//  Universidad Politécnica de Valencia
//  Hochschule Karlsruhe - Technik und Wirtschaft University of Applied Sciences
// --------------------------------------------------------------------------
#ifndef LAMBDA_RTWUTIL_H
#define LAMBDA_RTWUTIL_H

// Include Files
#include <cmath>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "math_functions\rt_nonfinite.h"
#include "math_functions\rtwtypes.h"
#include "LAMBDA_types.h"

// --------------------------------------------------------------------------
//                          Function Declarations
// --------------------------------------------------------------------------
extern double rt_hypotd_snf(double u0, double u1);

extern double rt_roundd_snf(double u);

#endif
