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


struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

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

struct emxArray__common
{
  void *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
//                          struct_emxArray_char_T
// --------------------------------------------------------------------------

#ifndef struct_emxArray_char_T
#define struct_emxArray_char_T

struct emxArray_char_T
{
  char *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
//                          struct_emxArray_int32_T
// --------------------------------------------------------------------------

#ifndef struct_emxArray_int32_T
#define struct_emxArray_int32_T

struct emxArray_int32_T
{
  int *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif
// --------------------------------------------------------------------------

#endif

