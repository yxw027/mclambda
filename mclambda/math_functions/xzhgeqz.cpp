// Include Files
#include "rt_nonfinite.h"
#include "..\mclambda.h"
#include "xzhgeqz.h"
#include "xzlartg.cpp"
#include "sqrt.cpp"

// --------------------------------------------------------------------------
void xzhgeqz(int n, const creal_T A[], int ilo, int ihi, int *info, creal_T alpha1
             [], creal_T beta1[])
{
  creal_T b_A[n*n];
  int i;
  double eshift_re;
  double eshift_im;
  creal_T ctemp;
  double anorm;
  double scale;
  double reAij;
  double sumsq;
  double b_atol;
  boolean_T firstNonZero;
  int j;
  int jp1;
  double ascale;
  boolean_T failed;
  boolean_T guard1 = false;
  double imAij;
  boolean_T guard2 = false;
  int ifirst;
  int istart;
  int ilast;
  double temp2;
  int ilastm1;
  int ifrstm;
  int ilastm;
  int iiter;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  int jiter;
  int exitg1;
  boolean_T exitg2;
  creal_T b_ascale;
  creal_T shift;
  double ad22_re;
  double ad22_im;
  double t1_im;
  memcpy(&b_A[0], &A[0], n*n * sizeof(creal_T));
  *info = 0;
  for (i = 0; i < n; i++) {
    alpha1[i].re = 0.0;
    alpha1[i].im = 0.0;
    beta1[i].re = 1.0;
    beta1[i].im = 0.0;
  }

  eshift_re = 0.0;
  eshift_im = 0.0;
  ctemp.re = 0.0;
  ctemp.im = 0.0;
  anorm = 0.0;
  if (!(ilo > ihi)) {
    scale = 0.0;
    sumsq = 0.0;
    firstNonZero = true;
    for (j = ilo; j <= ihi; j++) {
      jp1 = j + 1;
      if (ihi < j + 1) {
        jp1 = ihi;
      }

      for (i = ilo; i <= jp1; i++) {
        reAij = A[(i + n * (j - 1)) - 1].re;
        imAij = A[(i + n * (j - 1)) - 1].im;
        if (reAij != 0.0) {
          anorm = std::abs(reAij);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = anorm;
            firstNonZero = false;
          } else if (scale < anorm) {
            temp2 = scale / anorm;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            scale = anorm;
          } else {
            temp2 = anorm / scale;
            sumsq += temp2 * temp2;
          }
        }

        if (imAij != 0.0) {
          anorm = std::abs(imAij);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = anorm;
            firstNonZero = false;
          } else if (scale < anorm) {
            temp2 = scale / anorm;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            scale = anorm;
          } else {
            temp2 = anorm / scale;
            sumsq += temp2 * temp2;
          }
        }
      }
    }

    anorm = scale * std::sqrt(sumsq);
  }

  reAij = 2.2204460492503131E-16 * anorm;
  b_atol = 2.2250738585072014E-308;
  if (reAij > 2.2250738585072014E-308) {
    b_atol = reAij;
  }

  reAij = 2.2250738585072014E-308;
  if (anorm > 2.2250738585072014E-308) {
    reAij = anorm;
  }

  ascale = 1.0 / reAij;
  failed = true;
  for (j = ihi; j + 1 < n+1; j++) {
    alpha1[j] = A[j + n * j];
  }

  guard1 = false;
  guard2 = false;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    ifrstm = ilo;
    ilastm = ihi;
    iiter = 0;
    goto60 = false;
    goto70 = false;
    goto90 = false;
    jiter = 1;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1)) {
        if (ilast + 1 == ilo) {
          goto60 = true;
        } else if (std::abs(b_A[ilast + n * ilastm1].re) + std::abs(b_A[ilast +
                    n * ilastm1].im) <= b_atol) {
          b_A[ilast + n * ilastm1].re = 0.0;
          b_A[ilast + n * ilastm1].im = 0.0;
          goto60 = true;
        } else {
          j = ilastm1;
          exitg2 = false;
          while ((!exitg2) && (j + 1 >= ilo)) {
            if (j + 1 == ilo) {
              firstNonZero = true;
            } else if (std::abs(b_A[j + n * (j - 1)].re) + std::abs(b_A[j + n *
                        (j - 1)].im) <= b_atol) {
              b_A[j + n * (j - 1)].re = 0.0;
              b_A[j + n * (j - 1)].im = 0.0;
              firstNonZero = true;
            } else {
              firstNonZero = false;
            }

            if (firstNonZero) {
              ifirst = j + 1;
              goto70 = true;
              exitg2 = true;
            } else {
              j--;
            }
          }
        }

        if (goto60 || goto70) {
          firstNonZero = true;
        } else {
          firstNonZero = false;
        }

        if (!firstNonZero) {
          for (i = 0; i < n; i++) {
            alpha1[i].re = rtNaN;
            alpha1[i].im = 0.0;
            beta1[i].re = rtNaN;
            beta1[i].im = 0.0;
          }

          *info = 1;
          exitg1 = 1;
        } else if (goto60) {
          goto60 = false;
          alpha1[ilast] = b_A[ilast + n * ilast];
          ilast = ilastm1;
          ilastm1--;
          if (ilast + 1 < ilo) {
            failed = false;
            guard2 = true;
            exitg1 = 1;
          } else {
            iiter = 0;
            eshift_re = 0.0;
            eshift_im = 0.0;
            ilastm = ilast + 1;
            if (ifrstm > ilast + 1) {
              ifrstm = ilo;
            }

            jiter++;
          }
        } else {
          if (goto70) {
            goto70 = false;
            iiter++;
            ifrstm = ifirst;
            if (iiter - iiter / 10 * 10 != 0) {
              anorm = ascale * b_A[ilastm1 + n * ilastm1].re;
              reAij = ascale * b_A[ilastm1 + n * ilastm1].im;
              if (reAij == 0.0) {
                shift.re = anorm / 0.28867513459481292;
                shift.im = 0.0;
              } else if (anorm == 0.0) {
                shift.re = 0.0;
                shift.im = reAij / 0.28867513459481292;
              } else {
                shift.re = anorm / 0.28867513459481292;
                shift.im = reAij / 0.28867513459481292;
              }

              anorm = ascale * b_A[ilast + n * ilast].re;
              reAij = ascale * b_A[ilast + n * ilast].im;
              if (reAij == 0.0) {
                ad22_re = anorm / 0.28867513459481292;
                ad22_im = 0.0;
              } else if (anorm == 0.0) {
                ad22_re = 0.0;
                ad22_im = reAij / 0.28867513459481292;
              } else {
                ad22_re = anorm / 0.28867513459481292;
                ad22_im = reAij / 0.28867513459481292;
              }

              temp2 = 0.5 * (shift.re + ad22_re);
              t1_im = 0.5 * (shift.im + ad22_im);
              anorm = ascale * b_A[ilastm1 + n * ilast].re;
              reAij = ascale * b_A[ilastm1 + n * ilast].im;
              if (reAij == 0.0) {
                sumsq = anorm / 0.28867513459481292;
                imAij = 0.0;
              } else if (anorm == 0.0) {
                sumsq = 0.0;
                imAij = reAij / 0.28867513459481292;
              } else {
                sumsq = anorm / 0.28867513459481292;
                imAij = reAij / 0.28867513459481292;
              }

              anorm = ascale * b_A[ilast + n * ilastm1].re;
              reAij = ascale * b_A[ilast + n * ilastm1].im;
              if (reAij == 0.0) {
                scale = anorm / 0.28867513459481292;
                anorm = 0.0;
              } else if (anorm == 0.0) {
                scale = 0.0;
                anorm = reAij / 0.28867513459481292;
              } else {
                scale = anorm / 0.28867513459481292;
                anorm = reAij / 0.28867513459481292;
              }

              reAij = shift.re * ad22_im + shift.im * ad22_re;
              shift.re = ((temp2 * temp2 - t1_im * t1_im) + (sumsq * scale -
                imAij * anorm)) - (shift.re * ad22_re - shift.im * ad22_im);
              shift.im = ((temp2 * t1_im + t1_im * temp2) + (sumsq * anorm +
                imAij * scale)) - reAij;
              b_sqrt(&shift);
              if ((temp2 - ad22_re) * shift.re + (t1_im - ad22_im) * shift.im <=
                  0.0) {
                shift.re += temp2;
                shift.im += t1_im;
              } else {
                shift.re = temp2 - shift.re;
                shift.im = t1_im - shift.im;
              }
            } else {
              anorm = ascale * b_A[ilast + n * ilastm1].re;
              reAij = ascale * b_A[ilast + n * ilastm1].im;
              if (reAij == 0.0) {
                sumsq = anorm / 0.28867513459481292;
                imAij = 0.0;
              } else if (anorm == 0.0) {
                sumsq = 0.0;
                imAij = reAij / 0.28867513459481292;
              } else {
                sumsq = anorm / 0.28867513459481292;
                imAij = reAij / 0.28867513459481292;
              }

              eshift_re += sumsq;
              eshift_im += imAij;
              shift.re = eshift_re;
              shift.im = eshift_im;
            }

            j = ilastm1;
            jp1 = ilastm1 + 1;
            exitg2 = false;
            while ((!exitg2) && (j + 1 > ifirst)) {
              istart = j + 1;
              ctemp.re = ascale * b_A[j + n * j].re - shift.re *
                0.28867513459481292;
              ctemp.im = ascale * b_A[j + n * j].im - shift.im *
                0.28867513459481292;
              anorm = std::abs(ctemp.re) + std::abs(ctemp.im);
              temp2 = ascale * (std::abs(b_A[jp1 + n * j].re) + std::abs
                                (b_A[jp1 + n * j].im));
              reAij = anorm;
              if (temp2 > anorm) {
                reAij = temp2;
              }

              if ((reAij < 1.0) && (reAij != 0.0)) {
                anorm /= reAij;
                temp2 /= reAij;
              }

              if ((std::abs(b_A[j + n * (j - 1)].re) + std::abs(b_A[j + n * (j
                     - 1)].im)) * temp2 <= anorm * b_atol) {
                goto90 = true;
                exitg2 = true;
              } else {
                jp1 = j;
                j--;
              }
            }

            if (!goto90) {
              istart = ifirst;
              ctemp.re = ascale * b_A[(ifirst + n * (ifirst - 1)) - 1].re -
                shift.re * 0.28867513459481292;
              ctemp.im = ascale * b_A[(ifirst + n * (ifirst - 1)) - 1].im -
                shift.im * 0.28867513459481292;
              goto90 = true;
            }
          }

          if (goto90) {
            goto90 = false;
            b_ascale.re = ascale * b_A[istart + n * (istart - 1)].re;
            b_ascale.im = ascale * b_A[istart + n * (istart - 1)].im;
            b_xzlartg(ctemp, b_ascale, &anorm, &shift);
            j = istart;
            jp1 = istart - 2;
            while (j < ilast + 1) {
              if (j > istart) {
                xzlartg(b_A[(j + n * jp1) - 1], b_A[j + n * jp1], &anorm,
                        &shift, &b_A[(j + n * jp1) - 1]);
                b_A[j + n * jp1].re = 0.0;
                b_A[j + n * jp1].im = 0.0;
              }

              for (jp1 = j - 1; jp1 + 1 <= ilastm; jp1++) {
                ad22_re = anorm * b_A[(j + n * jp1) - 1].re + (shift.re * b_A[j
                  + n * jp1].re - shift.im * b_A[j + n * jp1].im);
                ad22_im = anorm * b_A[(j + n * jp1) - 1].im + (shift.re * b_A[j
                  + n * jp1].im + shift.im * b_A[j + n * jp1].re);
                reAij = b_A[(j + n * jp1) - 1].re;
                b_A[j + n * jp1].re = anorm * b_A[j + n * jp1].re - (shift.re *
                  b_A[(j + n * jp1) - 1].re + shift.im * b_A[(j + n * jp1) - 1]
                  .im);
                b_A[j + n * jp1].im = anorm * b_A[j + n * jp1].im - (shift.re *
                  b_A[(j + n * jp1) - 1].im - shift.im * reAij);
                b_A[(j + n * jp1) - 1].re = ad22_re;
                b_A[(j + n * jp1) - 1].im = ad22_im;
              }

              shift.re = -shift.re;
              shift.im = -shift.im;
              jp1 = j;
              if (ilast + 1 < j + 2) {
                jp1 = ilast - 1;
              }

              for (i = ifrstm - 1; i + 1 <= jp1 + 2; i++) {
                ad22_re = anorm * b_A[i + n * j].re + (shift.re * b_A[i + n *
                  (j - 1)].re - shift.im * b_A[i + n * (j - 1)].im);
                ad22_im = anorm * b_A[i + n * j].im + (shift.re * b_A[i + n *
                  (j - 1)].im + shift.im * b_A[i + n * (j - 1)].re);
                reAij = b_A[i + n * j].re;
                b_A[i + n * (j - 1)].re = anorm * b_A[i + n * (j - 1)].re -
                  (shift.re * b_A[i + n * j].re + shift.im * b_A[i + n * j].im);
                b_A[i + n * (j - 1)].im = anorm * b_A[i + n * (j - 1)].im -
                  (shift.re * b_A[i + n * j].im - shift.im * reAij);
                b_A[i + n * j].re = ad22_re;
                b_A[i + n * j].im = ad22_im;
              }

              jp1 = j - 1;
              j++;
            }
          }

          jiter++;
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }

  if (guard2) {
    if (failed) {
      *info = ilast + 1;
      for (jp1 = 0; jp1 + 1 <= ilast + 1; jp1++) {
        alpha1[jp1].re = rtNaN;
        alpha1[jp1].im = 0.0;
        beta1[jp1].re = rtNaN;
        beta1[jp1].im = 0.0;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1) {
    for (j = 0; j + 1 < ilo; j++) {
      alpha1[j] = b_A[j + n * j];
    }
  }
}
// --------------------------------------------------------------------------
