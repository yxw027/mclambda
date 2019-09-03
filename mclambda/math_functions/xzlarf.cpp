// Include Files
#include "rt_nonfinite.h"
#include "..\mclambda.h"
#include "xzlarf.h"

// --------------------------------------------------------------------------
void b_xzlarf(int n1, int m, int n, int iv0, double tau, double C[], int ic0, double
              work[])
{
  int lastv;
  int lastc;
  int i;
  boolean_T exitg2;
  int jy;
  int j;
  int i6;
  int ia;
  int exitg1;
  double c;
  int ix;
  if (tau != 0.0) {
    lastv = m;
    i = iv0 + m;
    while ((lastv > 0) && (C[i - 2] == 0.0)) {
      lastv--;
      i--;
    }

    lastc = n;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      i = ic0 + (lastc - 1) * n1;
      ia = i;
      do {
        exitg1 = 0;
        if (ia <= (i + lastv) - 1) {
          if (C[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    if (lastc != 0) {
      for (i = 1; i <= lastc; i++) {
        work[i - 1] = 0.0;
      }

      i = 0;
      i6 = ic0 + n1 * (lastc - 1);
      for (jy = ic0; jy <= i6; jy += n1) {
        ix = iv0;
        c = 0.0;
        j = (jy + lastv) - 1;
        for (ia = jy; ia <= j; ia++) {
          c += C[ia - 1] * C[ix - 1];
          ix++;
        }

        work[i] += c;
        i++;
      }
    }

    if (!(-tau == 0.0)) {
      i = ic0 - 1;
      jy = 0;
      for (j = 1; j <= lastc; j++) {
        if (work[jy] != 0.0) {
          c = work[jy] * -tau;
          ix = iv0;
          i6 = lastv + i;
          for (ia = i; ia + 1 <= i6; ia++) {
            C[ia] += C[ix - 1] * c;
            ix++;
          }
        }

        jy++;
        i += n1;
      }
    }
  }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void xzlarf(int n1, int n, int iv0, double tau, double C[], int ic0, double work[])
{
  int lastv;
  int lastc;
  int i;
  boolean_T exitg2;
  int jy;
  int ix;
  int j;
  int i5;
  int ia;
  int exitg1;
  double temp;
  if (tau != 0.0) {
    lastv = n;
    i = iv0 + n;
    while ((lastv > 0) && (C[i - 2] == 0.0)) {
      lastv--;
      i--;
    }

    lastc = n1;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      i = (ic0 + lastc) - 1;
      ia = i;
      do {
        exitg1 = 0;
        if (ia <= i + (lastv - 1) * n1) {
          if (C[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia += n1;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    if (lastc != 0) {
      for (i = 1; i <= lastc; i++) {
        work[i - 1] = 0.0;
      }

      ix = iv0;
      i5 = ic0 + n1 * (lastv - 1);
      for (jy = ic0; jy <= i5; jy += n1) {
        i = 0;
        j = (jy + lastc) - 1;
        for (ia = jy; ia <= j; ia++) {
          work[i] += C[ia - 1] * C[ix - 1];
          i++;
        }

        ix++;
      }
    }

    if (!(-tau == 0.0)) {
      i = ic0 - 1;
      jy = iv0 - 1;
      for (j = 1; j <= lastv; j++) {
        if (C[jy] != 0.0) {
          temp = C[jy] * -tau;
          ix = 0;
          i5 = lastc + i;
          for (ia = i; ia + 1 <= i5; ia++) {
            C[ia] += work[ix] * temp;
            ix++;
          }
        }

        jy++;
        i += n1;
      }
    }
  }
}
// --------------------------------------------------------------------------
