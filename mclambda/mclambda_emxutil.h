// --------------------------------------------------------------------------
//                        MC-LAMBDA emxutil HEADER FILE
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
//  Release date  : AUG-2019
//  Authors       : Hernández Olcina, Jorge
//
//  Master in Geomatics Engineering and Geoinformation
//  Universidad Politécnica de Valencia
//  Hochschule Karlsruhe - Technik und Wirtschaft University of Applied Sciences
// --------------------------------------------------------------------------
#ifndef LAMBDA_EMXUTIL_H
#define LAMBDA_EMXUTIL_H

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
extern void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, unsigned
  int elementSize);

extern void emxFreeMatrix_cell_wrap_3(cell_wrap_3 pMatrix[2]);

static void emxFreeStruct_cell_wrap_3(cell_wrap_3 *pStruct);

static void emxInitStruct_cell_wrap_3(cell_wrap_3 *pStruct);

extern void emxFree_char_T(emxArray_char_T **pEmxArray);

extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);

extern void emxFree_real_T(emxArray_real_T **pEmxArray);

extern void emxInitMatrix_cell_wrap_3(cell_wrap_3 pMatrix[2]);

extern void emxInit_char_T(emxArray_char_T **pEmxArray, int numDimensions);

extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);

extern void emxInit_int32_T1(emxArray_int32_T **pEmxArray, int numDimensions);

extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);

extern void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions);

#endif
