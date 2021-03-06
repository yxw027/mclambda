// Include Files
#include "rt_nonfinite.h"
#include "eig.h"
#include "schur.cpp"
#include "xzgeev.cpp"

// --------------------------------------------------------------------------
void eig(int n, const double A[], creal_T V[])
{
  boolean_T p;
  int info;
  int i;
  boolean_T exitg2;
  creal_T beta1[9999];
  double b_A[9999];
  creal_T T[9999];
  int exitg1;
  double V_re;
  double brm;
  double bim;
  double d;
  p = false;
  for (info = 0; info < n*n; info++) {
    if (p || rtIsInf(A[info]) || rtIsNaN(A[info])) {
      p = true;
    } else {
      p = false;
    }
  }

  if (p) {
    for (i = 0; i < n; i++) {
      V[i].re = rtNaN;
      V[i].im = 0.0;
    }
  } else {
    p = true;
    info = 0;
    exitg2 = false;
    while ((!exitg2) && (info < n)) {
      i = 0;
      do {
        exitg1 = 0;
        if (i <= info) {
          if (!(A[i + n * info] == A[info + n * i])) {
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
      memcpy(&b_A[0], &A[0], n*n * sizeof(double));
      schur(n, b_A, T);
      for (info = 0; info < n; info++) {
        V[info] = T[info + n * info];
      }
    } else {
      xzgeev(n, A, &info, V, beta1);
      for (info = 0; info < n; info++) {
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
// --------------------------------------------------------------------------
