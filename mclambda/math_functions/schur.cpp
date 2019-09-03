// Include Files
#include "rt_nonfinite.h"
#include "..\mclambda.h"
#include "schur.h"
#include "xdlanv2.cpp"
#include "xdhseqr.cpp"
#include "xgehrd.cpp"
#include "..\mclambda_rtwutil.h"

// --------------------------------------------------------------------------
void schur(int n, double A[], creal_T V[])
{
  boolean_T p;
  int istart;
  double Vr[n*n];
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
  for (istart = 0; istart < n*n; istart++) {
    if (p || rtIsInf(A[istart]) || rtIsNaN(A[istart])) {
      p = true;
    } else {
      p = false;
    }
  }

  if (p) {
    for (istart = 0; istart < n*n; istart++) {
      V[istart].re = rtNaN;
      V[istart].im = 0.0;
    }

    istart = 3;
    for (j = 0; j < 10; j++) {
      for (i = istart; i < n+1; i++) {
        V[(i + n * j) - 1].re = 0.0;
        V[(i + n * j) - 1].im = 0.0;
      }

      istart++;
    }
  } else {
    xgehrd(n, A);
    memcpy(&Vr[0], &A[0], n*n * sizeof(double));
    eml_dlahqr(n, Vr);
    istart = 4;
    for (j = 0; j < 9; j++) {
      for (i = istart; i < n+1; i++) {
        Vr[(i + n * j) - 1] = 0.0;
      }

      istart++;
    }

    for (istart = 0; istart < n*n; istart++) {
      V[istart].re = Vr[istart];
      V[istart].im = 0.0;
    }

    for (istart = 10; istart >= 0; istart += -1) {
      if (Vr[(istart + n * istart) + 1] != 0.0) {
        r = Vr[istart + n * istart];
        s = Vr[istart + n * (istart + 1)];
        t1_re = Vr[(istart + n * istart) + 1];
        t1_im = Vr[(istart + n * (istart + 1)) + 1];
        xdlanv2(&r, &s, &t1_re, &t1_im, &mu1_im, &rt1i, &mu1_re, &rt2i, &cs, &sn);
        mu1_re = mu1_im - Vr[(istart + n * (istart + 1)) + 1];
        r = rt_hypotd_snf(rt_hypotd_snf(mu1_re, rt1i), Vr[(istart + n * istart)
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

        s = Vr[(istart + n * istart) + 1] / r;
        for (j = istart; j + 1 < n+1; j++) {
          t1_re = V[istart + n * j].re;
          t1_im = V[istart + n * j].im;
          r = V[istart + n * j].re;
          V[istart + n * j].re = (mu1_re * V[istart + n * j].re + mu1_im *
            V[istart + n * j].im) + s * V[(istart + n * j) + 1].re;
          V[istart + n * j].im = (mu1_re * V[istart + n * j].im - mu1_im * r)
            + s * V[(istart + n * j) + 1].im;
          r = mu1_re * V[(istart + n * j) + 1].im + mu1_im * V[(istart + n * j)
            + 1].re;
          V[(istart + n * j) + 1].re = (mu1_re * V[(istart + n * j) + 1].re -
            mu1_im * V[(istart + n * j) + 1].im) - s * t1_re;
          V[(istart + n * j) + 1].im = r - s * t1_im;
        }

        for (i = 0; i + 1 <= istart + 2; i++) {
          t1_re = V[i + n * istart].re;
          t1_im = V[i + n * istart].im;
          r = mu1_re * V[i + n * istart].im + mu1_im * V[i + n * istart].re;
          V[i + n * istart].re = (mu1_re * V[i + n * istart].re - mu1_im * V[i
            + n * istart].im) + s * V[i + n * (istart + 1)].re;
          V[i + n * istart].im = r + s * V[i + n * (istart + 1)].im;
          r = V[i + n * (istart + 1)].re;
          V[i + n * (istart + 1)].re = (mu1_re * V[i + n * (istart + 1)].re +
            mu1_im * V[i + n * (istart + 1)].im) - s * t1_re;
          V[i + n * (istart + 1)].im = (mu1_re * V[i + n * (istart + 1)].im -
            mu1_im * r) - s * t1_im;
        }

        V[(istart + n * istart) + 1].re = 0.0;
        V[(istart + n * istart) + 1].im = 0.0;
      }
    }
  }
}
// --------------------------------------------------------------------------
