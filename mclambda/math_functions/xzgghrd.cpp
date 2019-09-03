// Include Files
#include "rt_nonfinite.h"
#include "..\mclambda.h"
#include "xzgghrd.h"
#include "xzlartg.h"

// --------------------------------------------------------------------------
void xzgghrd(int n, int ilo, int ihi, creal_T A[])
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
        xzlartg(A[(jrow + n * jcol) - 1], A[jrow + n * jcol], &c, &s, &A[(jrow
                 + n * jcol) - 1]);
        A[jrow + n * jcol].re = 0.0;
        A[jrow + n * jcol].im = 0.0;
        for (j = jcol + 1; j + 1 < n+1; j++) {
          stemp_re = c * A[(jrow + n * j) - 1].re + (s.re * A[jrow + n * j].re
            - s.im * A[jrow + n * j].im);
          stemp_im = c * A[(jrow + n * j) - 1].im + (s.re * A[jrow + n * j].im
            + s.im * A[jrow + n * j].re);
          A_im = A[(jrow + n * j) - 1].im;
          A_re = A[(jrow + n * j) - 1].re;
          A[jrow + n * j].re = c * A[jrow + n * j].re - (s.re * A[(jrow + n *
            j) - 1].re + s.im * A[(jrow + n * j) - 1].im);
          A[jrow + n * j].im = c * A[jrow + n * j].im - (s.re * A_im - s.im *
            A_re);
          A[(jrow + n * j) - 1].re = stemp_re;
          A[(jrow + n * j) - 1].im = stemp_im;
        }

        s.re = -s.re;
        s.im = -s.im;
        for (j = 0; j + 1 <= ihi; j++) {
          stemp_re = c * A[j + n * jrow].re + (s.re * A[j + n * (jrow - 1)].re
            - s.im * A[j + n * (jrow - 1)].im);
          stemp_im = c * A[j + n * jrow].im + (s.re * A[j + n * (jrow - 1)].im
            + s.im * A[j + n * (jrow - 1)].re);
          A_im = A[j + n * jrow].im;
          A_re = A[j + n * jrow].re;
          A[j + n * (jrow - 1)].re = c * A[j + n * (jrow - 1)].re - (s.re *
            A[j + n * jrow].re + s.im * A[j + n * jrow].im);
          A[j + n * (jrow - 1)].im = c * A[j + n * (jrow - 1)].im - (s.re *
            A_im - s.im * A_re);
          A[j + n * jrow].re = stemp_re;
          A[j + n * jrow].im = stemp_im;
        }
      }
    }
  }
}
// --------------------------------------------------------------------------
