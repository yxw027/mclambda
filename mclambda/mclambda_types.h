// --------------------------------------------------------------------------
//                        MC-LAMBDA TYPES HEADER FILE
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
//  Release date  : AUG-2019
//  Authors       : Hernández Olcina, Jorge
//
//  Master in Geomatics Engineering and Geoinformation
//  Universidad Politécnica de Valencia
//  Hochschule Karlsruhe - Technik und Wirtschaft University of Applied Sciences
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
//
//  DESCRIPTION:
//
//  Here you can find the different variable formats defined for the use 
//  of the MC-LAMBDA routine
//
// --------------------------------------------------------------------------
#ifndef LAMBDA_TYPES_H
#define LAMBDA_TYPES_H

// Include Files
#include "math_functions\rtwtypes.h"

// --------------------------------------------------------------------------
//                          struct_emxArray_real_T
// --------------------------------------------------------------------------
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T


typedef struct {
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
} emxArray_real_T;

#endif

typedef struct {
  emxArray_real_T *f1;
} cell_wrap_3;
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
//                          struct_emxArray__common
// --------------------------------------------------------------------------
#ifndef struct_emxArray__common
#define struct_emxArray__common

typedef struct {
  void *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
} emxArray__common;

#endif
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
//                          struct_emxArray_char_T
// --------------------------------------------------------------------------

#ifndef struct_emxArray_char_T
#define struct_emxArray_char_T

typedef struct {
  char *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
} emxArray_char_T;

#endif
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
//                          struct_emxArray_int32_T
// --------------------------------------------------------------------------

#ifndef struct_emxArray_int32_T
#define struct_emxArray_int32_T

typedef struct { 
  int *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
} emxArray_int32_T;

#endif
// --------------------------------------------------------------------------

#endif

