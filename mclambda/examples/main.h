// --------------------------------------------------------------------------
//                   MAIN MC-LAMBDA SOFTWARE HEADER FILE
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
//  Release date  : AUG-2019
//  Authors       : Hernández Olcina, Jorge
//
//  Master in Geomatics Engineering and Geoinformation
//  Universidad Politécnica de Valencia
//  Hochschule Karlsruhe - Technik und Wirtschaft University of Applied Sciences
// --------------------------------------------------------------------------
#ifndef MAIN_H
#define MAIN_H

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
static void argInit_12x12_real_T(double result[144]);

static void argInit_12x1_real_T(double result[12]);

static char argInit_char_T();

static double argInit_real_T();

static emxArray_char_T *c_argInit_UnboundedxUnbounded_c();

static void main_LAMBDA();

extern int main(int argc, const char * const argv[]);

#endif
