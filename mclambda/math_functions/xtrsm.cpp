// Include Files
#include "rt_nonfinite.h"
#include "..\mclambda.h"
#include "xtrsm.h"

// --------------------------------------------------------------------------
void xtrsm(int n, const double A[], double B[])
{
  int j;
  int jBcol;
  int k;
  int kAcol;
  int i;
  for (j = 0; j < n; j++) {
    jBcol = n * j;
    for (k = n-1; k >= 0; k += -1) {
      kAcol = n * k;
      if (B[k + jBcol] != 0.0) {
        B[k + jBcol] /= A[k + kAcol];
        for (i = 0; i + 1 <= k; i++) {
          B[i + jBcol] -= B[k + jBcol] * A[i + kAcol];
        }
      }
    }
  }
}
// --------------------------------------------------------------------------
