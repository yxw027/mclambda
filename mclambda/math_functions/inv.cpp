// Include Files
#include "rt_nonfinite.h"
#include "..\mclambda.h"
#include "inv.h"
#include "xtrsm.cpp"
#include "xzgetrf.cpp"

// --------------------------------------------------------------------------
void inv(int n, const double x[], double y[])
{
  int c;
  double b_x[n*n];
  int ipiv[n];
  int k;
  signed char p[n];
  int j;
  int i;
  for (c = 0; c < n*n; c++) {
    y[c] = 0.0;
    b_x[c] = x[c];
  }

  xzgetrf(n, b_x, ipiv, &c);
  for (c = 0; c < n; c++) {
    p[c] = (signed char)(1 + c);
  }

  for (k = 0; k < n-1; k++) {
    if (ipiv[k] > 1 + k) {
      c = p[ipiv[k] - 1];
      p[ipiv[k] - 1] = p[k];
      p[k] = (signed char)c;
    }
  }

  for (k = 0; k < n; k++) {
    c = p[k] - 1;
    y[k + n * (p[k] - 1)] = 1.0;
    for (j = k; j + 1 < n+1; j++) {
      if (y[j + n * c] != 0.0) {
        for (i = j + 1; i + 1 < n+1; i++) {
          y[i + n * c] -= y[j + n * c] * b_x[i + n * j];
        }
      }
    }
  }

  xtrsm(n, b_x, y);
}
// --------------------------------------------------------------------------
