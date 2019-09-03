// Include Files
#include "rt_nonfinite.h"
#include "..\mclambda.h"
#include "xzlangeM.h"
#include "schur.h"
#include "..\mclambda_rtwutil.h"

// --------------------------------------------------------------------------
double xzlangeM(int n, const creal_T x[])
{
  double y;
  int k;
  boolean_T exitg1;
  double absxk;
  y = 0.0;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < n*n)) {
    absxk = rt_hypotd_snf(x[k].re, x[k].im);
    if (rtIsNaN(absxk)) {
      y = rtNaN;
      exitg1 = true;
    } else {
      if (absxk > y) {
        y = absxk;
      }

      k++;
    }
  }

  return y;
}
// --------------------------------------------------------------------------
