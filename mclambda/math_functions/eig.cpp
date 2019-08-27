//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: eig.cpp
//
// MATLAB Coder version            : 3.3
// C/C++ source code generated on  : 08-Aug-2019 14:38:13
//

// Include Files
/*#include "rt_nonfinite.h"
#include "LAMBDA.h"*/
#include "eig.h"
#include "schur.cpp"
#include "xzgeev.cpp"

// Function Definitions

//
// Arguments    : const double A[144]
//                creal_T V[12]
// Return Type  : void
//
void eig(const double A[144], creal_T V[12])
{
  boolean_T p;
  int info;
  int i;
  boolean_T exitg2;
  creal_T beta1[12];
  double b_A[144];
  creal_T T[144];
  int exitg1;
  double V_re;
  double brm;
  double bim;
  double d;
  p = false;
  for (info = 0; info < 144; info++) {
    if (p || rtIsInf(A[info]) || rtIsNaN(A[info])) {
      p = true;
    } else {
      p = false;
    }
  }

  if (p) {
    for (i = 0; i < 12; i++) {
      V[i].re = rtNaN;
      V[i].im = 0.0;
    }
  } else {
    p = true;
    info = 0;
    exitg2 = false;
    while ((!exitg2) && (info < 12)) {
      i = 0;
      do {
        exitg1 = 0;
        if (i <= info) {
          if (!(A[i + 12 * info] == A[info + 12 * i])) {
            p = false;
            exitg1 = 1;
          } else {
            i++;
          }
        } else {
          info++;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }

    if (p) {
      memcpy(&b_A[0], &A[0], 144U * sizeof(double));
      schur(b_A, T);
      for (info = 0; info < 12; info++) {
        V[info] = T[info + 12 * info];
      }
    } else {
      xzgeev(A, &info, V, beta1);
      for (info = 0; info < 12; info++) {
        V_re = V[info].re;
        if (beta1[info].im == 0.0) {
          if (V[info].im == 0.0) {
            V[info].re /= beta1[info].re;
            V[info].im = 0.0;
          } else if (V[info].re == 0.0) {
            V[info].re = 0.0;
            V[info].im /= beta1[info].re;
          } else {
            V[info].re /= beta1[info].re;
            V[info].im /= beta1[info].re;
          }
        } else if (beta1[info].re == 0.0) {
          if (V[info].re == 0.0) {
            V[info].re = V[info].im / beta1[info].im;
            V[info].im = 0.0;
          } else if (V[info].im == 0.0) {
            V[info].re = 0.0;
            V[info].im = -(V_re / beta1[info].im);
          } else {
            V[info].re = V[info].im / beta1[info].im;
            V[info].im = -(V_re / beta1[info].im);
          }
        } else {
          brm = std::abs(beta1[info].re);
          bim = std::abs(beta1[info].im);
          if (brm > bim) {
            bim = beta1[info].im / beta1[info].re;
            d = beta1[info].re + bim * beta1[info].im;
            V[info].re = (V[info].re + bim * V[info].im) / d;
            V[info].im = (V[info].im - bim * V_re) / d;
          } else if (bim == brm) {
            if (beta1[info].re > 0.0) {
              bim = 0.5;
            } else {
              bim = -0.5;
            }

            if (beta1[info].im > 0.0) {
              d = 0.5;
            } else {
              d = -0.5;
            }

            V[info].re = (V[info].re * bim + V[info].im * d) / brm;
            V[info].im = (V[info].im * bim - V_re * d) / brm;
          } else {
            bim = beta1[info].re / beta1[info].im;
            d = beta1[info].im + bim * beta1[info].re;
            V[info].re = (bim * V[info].re + V[info].im) / d;
            V[info].im = (bim * V[info].im - V_re) / d;
          }
        }
      }
    }
  }
}

//
// File trailer for eig.cpp
//
// [EOF]
//
