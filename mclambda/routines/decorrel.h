// --------------------------------------------------------------------------
//                           DECORREL HEADER FILE
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
//  Release date  : AUG-2019
//  Authors       : Hernández Olcina, Jorge
//
//  Master in Geomatics Engineering and Geoinformation
//  Universidad Politécnica de Valencia
//  Hochschule Karlsruhe - Technik und Wirtschaft University of Applied Sciences
// --------------------------------------------------------------------------
#ifndef DECORREL_H
#define DECORREL_H

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
extern void decorrel(const double Qahat[144], const double ahat[12], double
                     Qzhat[144], double Z[144], double L[144], double D[12],
                     double zhat[12], double iZt[144]);

#endif
