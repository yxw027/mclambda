// Include Files
#include "rt_nonfinite.h"
#include "..\mclambda.h"
#include "xzgetrf.h"

// --------------------------------------------------------------------------
void xzgetrf(int n, double A[], int ipiv[], int *info)
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
  for (i9 = 0; i9 < n; i9++) {
    ipiv[i9] = 1 + i9;
  }

  *info = 0;
  for (j = 0; j < n-1; j++) {
    c = j * (n+1);
    iy = 0;
    ix = c;
    smax = std::abs(A[c]);
    for (jy = 2; jy <= n - j; jy++) {
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
        for (jy = 0; jy < n; jy++) {
          smax = A[ix];
          A[ix] = A[iy];
          A[iy] = smax;
          ix += n;
          iy += n;
        }
      }

      i9 = (c - j) + n;
      for (iy = c + 1; iy + 1 <= i9; iy++) {
        A[iy] /= A[c];
      }
    } else {
      *info = j + 1;
    }

    iy = c + n+1;
    jy = c + n;
    for (b_j = 1; b_j <= n-1 - j; b_j++) {
      smax = A[jy];
      if (A[jy] != 0.0) {
        ix = c + 1;
        i9 = (iy - j) + n-1;
        for (ijA = iy; ijA + 1 <= i9; ijA++) {
          A[ijA] += A[ix] * -smax;
          ix++;
        }
      }

      jy += n;
      iy += n;
    }
  }

  if ((*info == 0) && (!(A[n*n-1] != 0.0))) {
    *info = n;
  }
}
// --------------------------------------------------------------------------
