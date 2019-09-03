// Include Files
#include "rt_nonfinite.h"
#include "..\mclambda.h"
#include "xgehrd.h"
#include "xzlarf.cpp"
#include "schur.h"
#include "xnrm2.cpp"
#include "..\mclambda_rtwutil.h"

// --------------------------------------------------------------------------
void xgehrd(int n, double a[])
{
  double work[n];
  int i;
  int im1n;
  int in;
  int ia0;
  double alpha1;
  double d2;
  double xnorm;
  double tau[n-1];
  int knt;
  int i4;
  int k;
  memset(&work[0], 0, n * sizeof(double));
  for (i = 0; i < n-1; i++) {
    im1n = i * n + 2;
    in = (i + 1) * n;
    ia0 = i + 3;
    if (!(ia0 < n)) {
      ia0 = n;
    }

    ia0 += i * n;
    alpha1 = a[(i + n * i) + 1];
    d2 = 0.0;
    xnorm = xnrm2(10 - i, a, ia0);
    if (xnorm != 0.0) {
      xnorm = rt_hypotd_snf(a[(i + n * i) + 1], xnorm);
      if (a[(i + n * i) + 1] >= 0.0) {
        xnorm = -xnorm;
      }

      if (std::abs(xnorm) < 1.0020841800044864E-292) {
        knt = 0;
        i4 = (ia0 - i) + 9;
        do {
          knt++;
          for (k = ia0; k <= i4; k++) {
            a[k - 1] *= 9.9792015476736E+291;
          }

          xnorm *= 9.9792015476736E+291;
          alpha1 *= 9.9792015476736E+291;
        } while (!(std::abs(xnorm) >= 1.0020841800044864E-292));

        xnorm = rt_hypotd_snf(alpha1, xnrm2(10 - i, a, ia0));
        if (alpha1 >= 0.0) {
          xnorm = -xnorm;
        }

        d2 = (xnorm - alpha1) / xnorm;
        alpha1 = 1.0 / (alpha1 - xnorm);
        i4 = (ia0 - i) + 9;
        while (ia0 <= i4) {
          a[ia0 - 1] *= alpha1;
          ia0++;
        }

        for (k = 1; k <= knt; k++) {
          xnorm *= 1.0020841800044864E-292;
        }

        alpha1 = xnorm;
      } else {
        d2 = (xnorm - a[(i + n * i) + 1]) / xnorm;
        alpha1 = 1.0 / (a[(i + n * i) + 1] - xnorm);
        i4 = (ia0 - i) + 9;
        while (ia0 <= i4) {
          a[ia0 - 1] *= alpha1;
          ia0++;
        }

        alpha1 = xnorm;
      }
    }

    tau[i] = d2;
    a[(i + n * i) + 1] = 1.0;
    xzlarf(n, n-1 - i, i + im1n, tau[i], a, in + 1, work);
    b_xzlarf(n, n-1 - i, n-1 - i, i + im1n, tau[i], a, (i + in) + 2, work);
    a[(i + n * i) + 1] = alpha1;
  }
}
// --------------------------------------------------------------------------
