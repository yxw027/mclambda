//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xzgetrf.cpp
//
// MATLAB Coder version            : 3.3
// C/C++ source code generated on  : 08-Aug-2019 14:38:13
//

// Include Files
#include "rt_nonfinite.h"
#include "..\LAMBDA.h"
#include "xzgetrf.h"

// Function Definitions

//
// Arguments    : double A[144]
//                int ipiv[12]
//                int *info
// Return Type  : void
//
void xzgetrf(double A[144], int ipiv[12], int *info)
{
  int i9;
  int j;
  int c;
  int iy;
  int ix;
  double smax;
  int jy;
  double s;
  int b_j;
  int ijA;
  for (i9 = 0; i9 < 12; i9++) {
    ipiv[i9] = 1 + i9;
  }

  *info = 0;
  for (j = 0; j < 11; j++) {
    c = j * 13;
    iy = 0;
    ix = c;
    smax = std::abs(A[c]);
    for (jy = 2; jy <= 12 - j; jy++) {
      ix++;
      s = std::abs(A[ix]);
      if (s > smax) {
        iy = jy - 1;
        smax = s;
      }
    }

    if (A[c + iy] != 0.0) {
      if (iy != 0) {
        ipiv[j] = (j + iy) + 1;
        ix = j;
        iy += j;
        for (jy = 0; jy < 12; jy++) {
          smax = A[ix];
          A[ix] = A[iy];
          A[iy] = smax;
          ix += 12;
          iy += 12;
        }
      }

      i9 = (c - j) + 12;
      for (iy = c + 1; iy + 1 <= i9; iy++) {
        A[iy] /= A[c];
      }
    } else {
      *info = j + 1;
    }

    iy = c + 13;
    jy = c + 12;
    for (b_j = 1; b_j <= 11 - j; b_j++) {
      smax = A[jy];
      if (A[jy] != 0.0) {
        ix = c + 1;
        i9 = (iy - j) + 11;
        for (ijA = iy; ijA + 1 <= i9; ijA++) {
          A[ijA] += A[ix] * -smax;
          ix++;
        }
      }

      jy += 12;
      iy += 12;
    }
  }

  if ((*info == 0) && (!(A[143] != 0.0))) {
    *info = 12;
  }
}

//
// File trailer for xzgetrf.cpp
//
// [EOF]
//
