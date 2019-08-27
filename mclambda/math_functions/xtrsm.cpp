//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xtrsm.cpp
//
// MATLAB Coder version            : 3.3
// C/C++ source code generated on  : 08-Aug-2019 14:38:13
//

// Include Files
#include "rt_nonfinite.h"
#include "..\LAMBDA.h"
#include "xtrsm.h"

// Function Definitions

//
// Arguments    : const double A[144]
//                double B[144]
// Return Type  : void
//
void xtrsm(const double A[144], double B[144])
{
  int j;
  int jBcol;
  int k;
  int kAcol;
  int i;
  for (j = 0; j < 12; j++) {
    jBcol = 12 * j;
    for (k = 11; k >= 0; k += -1) {
      kAcol = 12 * k;
      if (B[k + jBcol] != 0.0) {
        B[k + jBcol] /= A[k + kAcol];
        for (i = 0; i + 1 <= k; i++) {
          B[i + jBcol] -= B[k + jBcol] * A[i + kAcol];
        }
      }
    }
  }
}

//
// File trailer for xtrsm.cpp
//
// [EOF]
//
