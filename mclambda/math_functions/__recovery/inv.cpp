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
  double b_x[999];
  double y_1[999];
  int ipiv[999];
  int k;
  signed char p[999];
  int j;
  int i;
  for (c = 0; c < n*n; c++) {
    y_1[c] = 0.0;
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
    y_1[k + n * (p[k] - 1)] = 1.0;
    for (j = k; j + 1 < n+1; j++) {
      if (y_1[j + n * c] != 0.0) {
        for (i = j + 1; i + 1 < n+1; i++) {
          y_1[i + n * c] -= y_1[j + n * c] * b_x[i + n * j];
        }
      }
    }
  }
  xtrsm(n, b_x, y_1);
  for (int j = 0; j < n*n; j++){
    y[j] = y_1[j];
  }
}
// --------------------------------------------------------------------------
