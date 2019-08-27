//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xzlarf.cpp
//
// MATLAB Coder version            : 3.3
// C/C++ source code generated on  : 08-Aug-2019 14:38:13
//

// Include Files
#include "rt_nonfinite.h"
#include "..\LAMBDA.h"
#include "xzlarf.h"

// Function Definitions

//
// Arguments    : int m
//                int n
//                int iv0
//                double tau
//                double C[144]
//                int ic0
//                double work[12]
// Return Type  : void
//
void b_xzlarf(int m, int n, int iv0, double tau, double C[144], int ic0, double
              work[12])
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
      i = ic0 + (lastc - 1) * 12;
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
      i6 = ic0 + 12 * (lastc - 1);
      for (jy = ic0; jy <= i6; jy += 12) {
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
        i += 12;
      }
    }
  }
}

//
// Arguments    : int n
//                int iv0
//                double tau
//                double C[144]
//                int ic0
//                double work[12]
// Return Type  : void
//
void xzlarf(int n, int iv0, double tau, double C[144], int ic0, double work[12])
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

    lastc = 12;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      i = (ic0 + lastc) - 1;
      ia = i;
      do {
        exitg1 = 0;
        if (ia <= i + (lastv - 1) * 12) {
          if (C[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia += 12;
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
      i5 = ic0 + 12 * (lastv - 1);
      for (jy = ic0; jy <= i5; jy += 12) {
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
        i += 12;
      }
    }
  }
}

//
// File trailer for xzlarf.cpp
//
// [EOF]
//
