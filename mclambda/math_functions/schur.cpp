//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: schur.cpp
//
// MATLAB Coder version            : 3.3
// C/C++ source code generated on  : 08-Aug-2019 14:38:13
//

// Include Files
#include "rt_nonfinite.h"
#include "..\LAMBDA.h"
#include "schur.h"
#include "xdlanv2.cpp"
#include "xdhseqr.cpp"
#include "xgehrd.cpp"
#include "..\LAMBDA_rtwutil.h"

// Function Definitions

//
// Arguments    : double A[144]
//                creal_T V[144]
// Return Type  : void
//
void schur(double A[144], creal_T V[144])
{
  boolean_T p;
  int istart;
  double Vr[144];
  int j;
  int i;
  double r;
  double s;
  double t1_re;
  double t1_im;
  double mu1_im;
  double rt1i;
  double mu1_re;
  double rt2i;
  double cs;
  double sn;
  p = false;
  for (istart = 0; istart < 144; istart++) {
    if (p || rtIsInf(A[istart]) || rtIsNaN(A[istart])) {
      p = true;
    } else {
      p = false;
    }
  }

  if (p) {
    for (istart = 0; istart < 144; istart++) {
      V[istart].re = rtNaN;
      V[istart].im = 0.0;
    }

    istart = 3;
    for (j = 0; j < 10; j++) {
      for (i = istart; i < 13; i++) {
        V[(i + 12 * j) - 1].re = 0.0;
        V[(i + 12 * j) - 1].im = 0.0;
      }

      istart++;
    }
  } else {
    xgehrd(A);
    memcpy(&Vr[0], &A[0], 144U * sizeof(double));
    eml_dlahqr(Vr);
    istart = 4;
    for (j = 0; j < 9; j++) {
      for (i = istart; i < 13; i++) {
        Vr[(i + 12 * j) - 1] = 0.0;
      }

      istart++;
    }

    for (istart = 0; istart < 144; istart++) {
      V[istart].re = Vr[istart];
      V[istart].im = 0.0;
    }

    for (istart = 10; istart >= 0; istart += -1) {
      if (Vr[(istart + 12 * istart) + 1] != 0.0) {
        r = Vr[istart + 12 * istart];
        s = Vr[istart + 12 * (istart + 1)];
        t1_re = Vr[(istart + 12 * istart) + 1];
        t1_im = Vr[(istart + 12 * (istart + 1)) + 1];
        xdlanv2(&r, &s, &t1_re, &t1_im, &mu1_im, &rt1i, &mu1_re, &rt2i, &cs, &sn);
        mu1_re = mu1_im - Vr[(istart + 12 * (istart + 1)) + 1];
        r = rt_hypotd_snf(rt_hypotd_snf(mu1_re, rt1i), Vr[(istart + 12 * istart)
                          + 1]);
        if (rt1i == 0.0) {
          mu1_re /= r;
          mu1_im = 0.0;
        } else if (mu1_re == 0.0) {
          mu1_re = 0.0;
          mu1_im = rt1i / r;
        } else {
          mu1_re /= r;
          mu1_im = rt1i / r;
        }

        s = Vr[(istart + 12 * istart) + 1] / r;
        for (j = istart; j + 1 < 13; j++) {
          t1_re = V[istart + 12 * j].re;
          t1_im = V[istart + 12 * j].im;
          r = V[istart + 12 * j].re;
          V[istart + 12 * j].re = (mu1_re * V[istart + 12 * j].re + mu1_im *
            V[istart + 12 * j].im) + s * V[(istart + 12 * j) + 1].re;
          V[istart + 12 * j].im = (mu1_re * V[istart + 12 * j].im - mu1_im * r)
            + s * V[(istart + 12 * j) + 1].im;
          r = mu1_re * V[(istart + 12 * j) + 1].im + mu1_im * V[(istart + 12 * j)
            + 1].re;
          V[(istart + 12 * j) + 1].re = (mu1_re * V[(istart + 12 * j) + 1].re -
            mu1_im * V[(istart + 12 * j) + 1].im) - s * t1_re;
          V[(istart + 12 * j) + 1].im = r - s * t1_im;
        }

        for (i = 0; i + 1 <= istart + 2; i++) {
          t1_re = V[i + 12 * istart].re;
          t1_im = V[i + 12 * istart].im;
          r = mu1_re * V[i + 12 * istart].im + mu1_im * V[i + 12 * istart].re;
          V[i + 12 * istart].re = (mu1_re * V[i + 12 * istart].re - mu1_im * V[i
            + 12 * istart].im) + s * V[i + 12 * (istart + 1)].re;
          V[i + 12 * istart].im = r + s * V[i + 12 * (istart + 1)].im;
          r = V[i + 12 * (istart + 1)].re;
          V[i + 12 * (istart + 1)].re = (mu1_re * V[i + 12 * (istart + 1)].re +
            mu1_im * V[i + 12 * (istart + 1)].im) - s * t1_re;
          V[i + 12 * (istart + 1)].im = (mu1_re * V[i + 12 * (istart + 1)].im -
            mu1_im * r) - s * t1_im;
        }

        V[(istart + 12 * istart) + 1].re = 0.0;
        V[(istart + 12 * istart) + 1].im = 0.0;
      }
    }
  }
}

//
// File trailer for schur.cpp
//
// [EOF]
//
