// --------------------------------------------------------------------------
//                        MC-LAMBDA emxAPI HEADER FILE
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
//  Release date  : AUG-2019
//  Authors       : Hernández Olcina, Jorge
//
//  Master in Geomatics Engineering and Geoinformation
//  Universidad Politécnica de Valencia
//  Hochschule Karlsruhe - Technik und Wirtschaft University of Applied Sciences
// --------------------------------------------------------------------------
#ifndef LAMBDA_EMXAPI_H
#define LAMBDA_EMXAPI_H

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
EXPORT emxArray_char_T *emxCreateND_char_T(int numDimensions, int *size);

EXPORT emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size);

EXPORT emxArray_char_T *emxCreateWrapperND_char_T(char *data, int numDimensions,
  int *size);

EXPORT emxArray_real_T *emxCreateWrapperND_real_T(double *data, int
  numDimensions, int *size);

EXPORT emxArray_char_T *emxCreateWrapper_char_T(char *data, int rows, int cols);

EXPORT emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols);

EXPORT emxArray_char_T *emxCreate_char_T(int rows, int cols);

EXPORT emxArray_real_T *emxCreate_real_T(int rows, int cols);

EXPORT void emxDestroyArray_char_T(emxArray_char_T *emxArray);

EXPORT void emxDestroyArray_real_T(emxArray_real_T *emxArray);

EXPORT void emxInitArray_char_T(emxArray_char_T **pEmxArray, int numDimensions);

EXPORT void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions);

#ifdef __cplusplus
}
#endif
#endif
