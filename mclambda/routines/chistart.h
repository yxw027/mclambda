// --------------------------------------------------------------------------
//                           CHISTART HEADER FILE
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
//  Release date  : AUG-2019
//  Authors       : Hernández Olcina, Jorge
//
//  Master in Geomatics Engineering and Geoinformation
//  Universidad Politécnica de Valencia
//  Hochschule Karlsruhe - Technik und Wirtschaft University of Applied Sciences
// --------------------------------------------------------------------------
#ifndef CHISTART_H
#define CHISTART_H

// Include Files
#include <cmath>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "..\math_functions\rt_nonfinite.h"
#include "..\math_functions\rtwtypes.h"
#include "..\mclambda_types.h"

// --------------------------------------------------------------------------
//                          Function Declarations
// --------------------------------------------------------------------------
double chistart(int n, const double D[], const double L[], const double ahat[], double ncands);

#endif

