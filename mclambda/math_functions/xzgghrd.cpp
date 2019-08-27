//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xzgghrd.cpp
//
// MATLAB Coder version            : 3.3
// C/C++ source code generated on  : 08-Aug-2019 14:38:13
//

// Include Files
#include "rt_nonfinite.h"
#include "..\LAMBDA.h"
#include "xzgghrd.h"
#include "xzlartg.h"

// Function Definitions

//
// Arguments    : int ilo
//                int ihi
//                creal_T A[144]
// Return Type  : void
//
void xzgghrd(int ilo, int ihi, creal_T A[144])
{
  int jcol;
  int jrow;
  double c;
  creal_T s;
  int j;
  double stemp_re;
  double stemp_im;
  double A_im;
  double A_re;
  if (!(ihi < ilo + 2)) {
    for (jcol = ilo - 1; jcol + 1 < ihi - 1; jcol++) {
      for (jrow = ihi - 1; jrow + 1 > jcol + 2; jrow--) {
        xzlartg(A[(jrow + 12 * jcol) - 1], A[jrow + 12 * jcol], &c, &s, &A[(jrow
                 + 12 * jcol) - 1]);
        A[jrow + 12 * jcol].re = 0.0;
        A[jrow + 12 * jcol].im = 0.0;
        for (j = jcol + 1; j + 1 < 13; j++) {
          stemp_re = c * A[(jrow + 12 * j) - 1].re + (s.re * A[jrow + 12 * j].re
            - s.im * A[jrow + 12 * j].im);
          stemp_im = c * A[(jrow + 12 * j) - 1].im + (s.re * A[jrow + 12 * j].im
            + s.im * A[jrow + 12 * j].re);
          A_im = A[(jrow + 12 * j) - 1].im;
          A_re = A[(jrow + 12 * j) - 1].re;
          A[jrow + 12 * j].re = c * A[jrow + 12 * j].re - (s.re * A[(jrow + 12 *
            j) - 1].re + s.im * A[(jrow + 12 * j) - 1].im);
          A[jrow + 12 * j].im = c * A[jrow + 12 * j].im - (s.re * A_im - s.im *
            A_re);
          A[(jrow + 12 * j) - 1].re = stemp_re;
          A[(jrow + 12 * j) - 1].im = stemp_im;
        }

        s.re = -s.re;
        s.im = -s.im;
        for (j = 0; j + 1 <= ihi; j++) {
          stemp_re = c * A[j + 12 * jrow].re + (s.re * A[j + 12 * (jrow - 1)].re
            - s.im * A[j + 12 * (jrow - 1)].im);
          stemp_im = c * A[j + 12 * jrow].im + (s.re * A[j + 12 * (jrow - 1)].im
            + s.im * A[j + 12 * (jrow - 1)].re);
          A_im = A[j + 12 * jrow].im;
          A_re = A[j + 12 * jrow].re;
          A[j + 12 * (jrow - 1)].re = c * A[j + 12 * (jrow - 1)].re - (s.re *
            A[j + 12 * jrow].re + s.im * A[j + 12 * jrow].im);
          A[j + 12 * (jrow - 1)].im = c * A[j + 12 * (jrow - 1)].im - (s.re *
            A_im - s.im * A_re);
          A[j + 12 * jrow].re = stemp_re;
          A[j + 12 * jrow].im = stemp_im;
        }
      }
    }
  }
}

//
// File trailer for xzgghrd.cpp
//
// [EOF]
//
