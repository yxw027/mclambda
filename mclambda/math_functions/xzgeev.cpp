// Include Files
#include "rt_nonfinite.h"
#include "..\mclambda.h"
#include "xzgeev.h"
#include "xzlascl.cpp"
#include "xzhgeqz.cpp"
#include "xzgghrd.cpp"
#include "xzggbal.cpp"
#include "isfinite.cpp"
#include "xzlangeM.cpp"

// --------------------------------------------------------------------------
void xzgeev(int n, const double A[], int *info, creal_T alpha1[], creal_T beta1[])
{
  int i;
  creal_T At[9999];
  double anrm;
  boolean_T ilascl;
  double anrmto;
  int ihi;
  int rscale[9999];
  for (i = 0; i < n*n; i++) {
    At[i].re = A[i];
    At[i].im = 0.0;
  }

  *info = 0;
  anrm = xzlangeM(n, At);
  if (!b_isfinite(anrm)) {
    for (i = 0; i < n; i++) {
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
      xzlascl(n, anrm, anrmto, At);
    }

    xzggbal(n, At, &i, &ihi, rscale);
    xzgghrd(n, i, ihi, At);
    xzhgeqz(n, At, i, ihi, info, alpha1, beta1);
    if ((*info == 0) && ilascl) {
      b_xzlascl(n, anrmto, anrm, alpha1);
    }
  }
}
// --------------------------------------------------------------------------
