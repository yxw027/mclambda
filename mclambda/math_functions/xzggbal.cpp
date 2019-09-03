// Include Files
#include "rt_nonfinite.h"
#include "..\mclambda.h"
#include "xzggbal.h"

// --------------------------------------------------------------------------
void xzggbal(int n, creal_T A[], int *ilo, int *ihi, int rscale[])
{
  int i;
  int exitg2;
  int j;
  boolean_T found;
  int ii;
  boolean_T exitg3;
  int nzcount;
  int jj;
  double atmp_re;
  boolean_T exitg4;
  double atmp_im;
  int exitg1;
  for (i = 0; i < n; i++) {
    rscale[i] = 1;
  }

  *ilo = 1;
  *ihi = n;
  do {
    exitg2 = 0;
    i = 0;
    j = 0;
    found = false;
    ii = *ihi;
    exitg3 = false;
    while ((!exitg3) && (ii > 0)) {
      nzcount = 0;
      i = ii;
      j = *ihi;
      jj = 1;
      exitg4 = false;
      while ((!exitg4) && (jj <= *ihi)) {
        if ((A[(ii + n * (jj - 1)) - 1].re != 0.0) || (A[(ii + n * (jj - 1)) -
             1].im != 0.0) || (ii == jj)) {
          if (nzcount == 0) {
            j = jj;
            nzcount = 1;
            jj++;
          } else {
            nzcount = 2;
            exitg4 = true;
          }
        } else {
          jj++;
        }
      }

      if (nzcount < 2) {
        found = true;
        exitg3 = true;
      } else {
        ii--;
      }
    }

    if (!found) {
      exitg2 = 2;
    } else {
      if (i != *ihi) {
        for (ii = 0; ii < n; ii++) {
          atmp_re = A[(i + n * ii) - 1].re;
          atmp_im = A[(i + n * ii) - 1].im;
          A[(i + n * ii) - 1] = A[(*ihi + n * ii) - 1];
          A[(*ihi + n * ii) - 1].re = atmp_re;
          A[(*ihi + n * ii) - 1].im = atmp_im;
        }
      }

      if (j != *ihi) {
        for (ii = 0; ii + 1 <= *ihi; ii++) {
          atmp_re = A[ii + n * (j - 1)].re;
          atmp_im = A[ii + n * (j - 1)].im;
          A[ii + n * (j - 1)] = A[ii + n * (*ihi - 1)];
          A[ii + n * (*ihi - 1)].re = atmp_re;
          A[ii + n * (*ihi - 1)].im = atmp_im;
        }
      }

      rscale[*ihi - 1] = j;
      (*ihi)--;
      if (*ihi == 1) {
        rscale[0] = 1;
        exitg2 = 1;
      }
    }
  } while (exitg2 == 0);

  if (exitg2 == 1) {
  } else {
    do {
      exitg1 = 0;
      i = 0;
      j = 0;
      found = false;
      jj = *ilo;
      exitg3 = false;
      while ((!exitg3) && (jj <= *ihi)) {
        nzcount = 0;
        i = *ihi;
        j = jj;
        ii = *ilo;
        exitg4 = false;
        while ((!exitg4) && (ii <= *ihi)) {
          if ((A[(ii + n * (jj - 1)) - 1].re != 0.0) || (A[(ii + n * (jj - 1))
               - 1].im != 0.0) || (ii == jj)) {
            if (nzcount == 0) {
              i = ii;
              nzcount = 1;
              ii++;
            } else {
              nzcount = 2;
              exitg4 = true;
            }
          } else {
            ii++;
          }
        }

        if (nzcount < 2) {
          found = true;
          exitg3 = true;
        } else {
          jj++;
        }
      }

      if (!found) {
        exitg1 = 1;
      } else {
        if (i != *ilo) {
          for (ii = *ilo - 1; ii + 1 < n+1; ii++) {
            atmp_re = A[(i + n * ii) - 1].re;
            atmp_im = A[(i + n * ii) - 1].im;
            A[(i + n * ii) - 1] = A[(*ilo + n * ii) - 1];
            A[(*ilo + n * ii) - 1].re = atmp_re;
            A[(*ilo + n * ii) - 1].im = atmp_im;
          }
        }

        if (j != *ilo) {
          for (ii = 0; ii + 1 <= *ihi; ii++) {
            atmp_re = A[ii + n * (j - 1)].re;
            atmp_im = A[ii + n * (j - 1)].im;
            A[ii + n * (j - 1)] = A[ii + n * (*ilo - 1)];
            A[ii + n * (*ilo - 1)].re = atmp_re;
            A[ii + n * (*ilo - 1)].im = atmp_im;
          }
        }

        rscale[*ilo - 1] = j;
        (*ilo)++;
        if (*ilo == *ihi) {
          rscale[*ilo - 1] = *ilo;
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}
// --------------------------------------------------------------------------
