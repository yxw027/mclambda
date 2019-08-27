//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: inv.cpp
//
// MATLAB Coder version            : 3.3
// C/C++ source code generated on  : 08-Aug-2019 14:38:13
//

// Include Files
#include "rt_nonfinite.h"
#include "..\LAMBDA.h"
#include "inv.h"
#include "xtrsm.cpp"
#include "xzgetrf.cpp"

// Function Definitions

//
// Arguments    : const double x[144]
//                double y[144]
// Return Type  : void
//
void inv(const double x[144], double y[144])
{
  int c;
  double b_x[144];
  int ipiv[12];
  int k;
  signed char p[12];
  int j;
  int i;
  for (c = 0; c < 144; c++) {
    y[c] = 0.0;
    b_x[c] = x[c];
  }

  xzgetrf(b_x, ipiv, &c);
  for (c = 0; c < 12; c++) {
    p[c] = (signed char)(1 + c);
  }

  for (k = 0; k < 11; k++) {
    if (ipiv[k] > 1 + k) {
      c = p[ipiv[k] - 1];
      p[ipiv[k] - 1] = p[k];
      p[k] = (signed char)c;
    }
  }

  for (k = 0; k < 12; k++) {
    c = p[k] - 1;
    y[k + 12 * (p[k] - 1)] = 1.0;
    for (j = k; j + 1 < 13; j++) {
      if (y[j + 12 * c] != 0.0) {
        for (i = j + 1; i + 1 < 13; i++) {
          y[i + 12 * c] -= y[j + 12 * c] * b_x[i + 12 * j];
        }
      }
    }
  }

  xtrsm(b_x, y);
}

//
// File trailer for inv.cpp
//
// [EOF]
//
