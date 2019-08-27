//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xzgeev.cpp
//
// MATLAB Coder version            : 3.3
// C/C++ source code generated on  : 08-Aug-2019 14:38:13
//

// Include Files
#include "rt_nonfinite.h"
#include "..\LAMBDA.h"
#include "xzgeev.h"
#include "xzlascl.cpp"
#include "xzhgeqz.cpp"
#include "xzgghrd.cpp"
#include "xzggbal.cpp"
#include "isfinite.cpp"
#include "xzlangeM.cpp"

// Function Definitions

//
// Arguments    : const double A[144]
//                int *info
//                creal_T alpha1[12]
//                creal_T beta1[12]
// Return Type  : void
//
void xzgeev(const double A[144], int *info, creal_T alpha1[12], creal_T beta1[12])
{
  int i;
  creal_T At[144];
  double anrm;
  boolean_T ilascl;
  double anrmto;
  int ihi;
  int rscale[12];
  for (i = 0; i < 144; i++) {
    At[i].re = A[i];
    At[i].im = 0.0;
  }

  *info = 0;
  anrm = xzlangeM(At);
  if (!b_isfinite(anrm)) {
    for (i = 0; i < 12; i++) {
      alpha1[i].re = rtNaN;
      alpha1[i].im = 0.0;
      beta1[i].re = rtNaN;
      beta1[i].im = 0.0;
    }
  } else {
    ilascl = false;
    anrmto = anrm;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      anrmto = 6.7178761075670888E-139;
      ilascl = true;
    } else {
      if (anrm > 1.4885657073574029E+138) {
        anrmto = 1.4885657073574029E+138;
        ilascl = true;
      }
    }

    if (ilascl) {
      xzlascl(anrm, anrmto, At);
    }

    xzggbal(At, &i, &ihi, rscale);
    xzgghrd(i, ihi, At);
    xzhgeqz(At, i, ihi, info, alpha1, beta1);
    if ((*info == 0) && ilascl) {
      b_xzlascl(anrmto, anrm, alpha1);
    }
  }
}

//
// File trailer for xzgeev.cpp
//
// [EOF]
//
