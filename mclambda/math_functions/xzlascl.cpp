// Include Files
#include "rt_nonfinite.h"
#include "..\mclambda.h"
#include "xzlascl.h"

// --------------------------------------------------------------------------
void b_xzlascl(int n, double cfrom, double cto, creal_T A[])
{
  double cfromc;
  double ctoc;
  boolean_T notdone;
  double cfrom1;
  double cto1;
  double mul;
  int i8;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    cfrom1 = cfromc * 2.0041683600089728E-292;
    cto1 = ctoc / 4.9896007738368E+291;
    if ((std::abs(cfrom1) > std::abs(ctoc)) && (ctoc != 0.0)) {
      mul = 2.0041683600089728E-292;
      cfromc = cfrom1;
    } else if (std::abs(cto1) > std::abs(cfromc)) {
      mul = 4.9896007738368E+291;
      ctoc = cto1;
    } else {
      mul = ctoc / cfromc;
      notdone = false;
    }

    for (i8 = 0; i8 < n; i8++) {
      A[i8].re *= mul;
      A[i8].im *= mul;
    }
  }
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void xzlascl(int n, double cfrom, double cto, creal_T A[])
{
  double cfromc;
  double ctoc;
  boolean_T notdone;
  double cfrom1;
  double cto1;
  double mul;
  int i7;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    cfrom1 = cfromc * 2.0041683600089728E-292;
    cto1 = ctoc / 4.9896007738368E+291;
    if ((std::abs(cfrom1) > std::abs(ctoc)) && (ctoc != 0.0)) {
      mul = 2.0041683600089728E-292;
      cfromc = cfrom1;
    } else if (std::abs(cto1) > std::abs(cfromc)) {
      mul = 4.9896007738368E+291;
      ctoc = cto1;
    } else {
      mul = ctoc / cfromc;
      notdone = false;
    }

    for (i7 = 0; i7 < n*n; i7++) {
      A[i7].re *= mul;
      A[i7].im *= mul;
    }
  }
}
// --------------------------------------------------------------------------