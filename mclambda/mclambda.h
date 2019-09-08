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
#include <stdio.h>
#include <stdarg.h>
#include <time.h>
#include <ctype.h>
#include "math_functions\rt_nonfinite.h"
#include "math_functions\rtwtypes.h"
#include "mclambda_types.h"
#ifdef WIN32
#include <winsock2.h>
#include <windows.h>
#else
#include <pthread.h>
#endif
#ifdef __cplusplus
extern "C" {
#endif

#ifdef WIN_DLL
#define EXPORT __declspec(dllexport) /* for Windows DLL */
#else
#define EXPORT
#endif

// --------------------------------------------------------------------------
//                          Function Declarations
// --------------------------------------------------------------------------
static double rt_remd_snf(double u0, double u1);

EXPORT void mclambda(int n, double ahat[], double Qahat[], double method, double
            param, const emxArray_char_T *type, double value, emxArray_real_T
            *afixed, emxArray_real_T *sqnorm, double *Ps, double Qzhat[],
            double Z[], double *nfixed, double *mu);

#ifdef __cplusplus
}
#endif
#endif
