// Include Files
#include "rt_nonfinite.h"
#include "..\mclambda.h"
#include "xdhseqr.h"
#include "xzlarfg.cpp"

// --------------------------------------------------------------------------
int eml_dlahqr(int n, double h[])
{
  int info;
  int j;
  int i;
  boolean_T exitg1;
  int L;
  boolean_T goto150;
  int its;
  boolean_T exitg2;
  int k;
  boolean_T exitg3;
  double tst;
  double htmp1;
  double aa;
  double ab;
  double ba;
  double rt2r;
  double rt1r;
  double s;
  double sn;
  int nr;
  int hoffset;
  int m;
  double d3;
  int b_k;
  double v[3];
  info = 0;
  for (j = 0; j < 9; j++) {
    h[(j + n * j) + 2] = 0.0;
    h[(j + n * j) + 3] = 0.0;
  }

  h[119] = 0.0;
  i = n-1;
  exitg1 = false;
  while ((!exitg1) && (i + 1 >= 1)) {
    L = 1;
    goto150 = false;
    its = 0;
    exitg2 = false;
    while ((!exitg2) && (its < 31)) {
      k = i;
      exitg3 = false;
      while ((!exitg3) && ((k + 1 > L) && (!(std::abs(h[k + n * (k - 1)]) <=
                1.2025010160053837E-291)))) {
        tst = std::abs(h[(k + n * (k - 1)) - 1]) + std::abs(h[k + n * k]);
        if (tst == 0.0) {
          if (k - 1 >= 1) {
            tst = std::abs(h[(k + n * (k - 2)) - 1]);
          }

          if (k + 2 <= n) {
            tst += std::abs(h[(k + n * k) + 1]);
          }
        }

        if (std::abs(h[k + n * (k - 1)]) <= 2.2204460492503131E-16 * tst) {
          htmp1 = std::abs(h[k + n * (k - 1)]);
          tst = std::abs(h[(k + n * k) - 1]);
          if (htmp1 > tst) {
            ab = htmp1;
            ba = tst;
          } else {
            ab = tst;
            ba = htmp1;
          }

          htmp1 = std::abs(h[k + n * k]);
          tst = std::abs(h[(k + n * (k - 1)) - 1] - h[k + n * k]);
          if (htmp1 > tst) {
            aa = htmp1;
            htmp1 = tst;
          } else {
            aa = tst;
          }

          s = aa + ab;
          tst = 2.2204460492503131E-16 * (htmp1 * (aa / s));
          if ((1.2025010160053837E-291 > tst) || rtIsNaN(tst)) {
            d3 = 1.2025010160053837E-291;
          } else {
            d3 = tst;
          }

          if (ba * (ab / s) <= d3) {
            exitg3 = true;
          } else {
            k--;
          }
        } else {
          k--;
        }
      }

      L = k + 1;
      if (k + 1 > 1) {
        h[k + n * (k - 1)] = 0.0;
      }

      if (k + 1 >= i) {
        goto150 = true;
        exitg2 = true;
      } else {
        if (its == 10) {
          s = std::abs(h[(k + n * k) + 1]) + std::abs(h[(k + n * (k + 1)) + 2]);
          tst = 0.75 * s + h[k + n * k];
          aa = -0.4375 * s;
          htmp1 = s;
          ba = tst;
        } else if (its == 20) {
          s = std::abs(h[i + n * (i - 1)]) + std::abs(h[(i + n * (i - 2)) - 1]);
          tst = 0.75 * s + h[i + n * i];
          aa = -0.4375 * s;
          htmp1 = s;
          ba = tst;
        } else {
          tst = h[(i + n * (i - 1)) - 1];
          htmp1 = h[i + n * (i - 1)];
          aa = h[(i + n * i) - 1];
          ba = h[i + n * i];
        }

        s = ((std::abs(tst) + std::abs(aa)) + std::abs(htmp1)) + std::abs(ba);
        if (s == 0.0) {
          rt1r = 0.0;
          ab = 0.0;
          rt2r = 0.0;
          ba = 0.0;
        } else {
          tst /= s;
          htmp1 /= s;
          aa /= s;
          ba /= s;
          ab = (tst + ba) / 2.0;
          tst = (tst - ab) * (ba - ab) - aa * htmp1;
          htmp1 = std::sqrt(std::abs(tst));
          if (tst >= 0.0) {
            rt1r = ab * s;
            rt2r = rt1r;
            ab = htmp1 * s;
            ba = -ab;
          } else {
            rt1r = ab + htmp1;
            rt2r = ab - htmp1;
            if (std::abs(rt1r - ba) <= std::abs(rt2r - ba)) {
              rt1r *= s;
              rt2r = rt1r;
            } else {
              rt2r *= s;
              rt1r = rt2r;
            }

            ab = 0.0;
            ba = 0.0;
          }
        }

        m = i - 2;
        exitg3 = false;
        while ((!exitg3) && (m + 1 >= k + 1)) {
          s = (std::abs(h[m + n * m] - rt2r) + std::abs(ba)) + std::abs(h[(m +
            n * m) + 1]);
          tst = h[(m + n * m) + 1] / s;
          v[0] = (tst * h[m + n * (m + 1)] + (h[m + n * m] - rt1r) * ((h[m +
                    n * m] - rt2r) / s)) - ab * (ba / s);
          v[1] = tst * (((h[m + n * m] + h[(m + n * (m + 1)) + 1]) - rt1r) -
                        rt2r);
          v[2] = tst * h[(m + n * (m + 1)) + 2];
          s = (std::abs(v[0]) + std::abs(v[1])) + std::abs(v[2]);
          tst = v[0] / s;
          v[0] /= s;
          htmp1 = v[1] / s;
          v[1] /= s;
          aa = v[2] / s;
          v[2] /= s;
          if ((m + 1 == k + 1) || (std::abs(h[m + n * (m - 1)]) * (std::abs
                (htmp1) + std::abs(aa)) <= 2.2204460492503131E-16 * std::abs(tst)
               * ((std::abs(h[(m + n * (m - 1)) - 1]) + std::abs(h[m + n * m]))
                  + std::abs(h[(m + n * (m + 1)) + 1])))) {
            exitg3 = true;
          } else {
            m--;
          }
        }

        for (b_k = m; b_k + 1 <= i; b_k++) {
          nr = (i - b_k) + 1;
          if (3 < nr) {
            nr = 3;
          }

          if (b_k + 1 > m + 1) {
            hoffset = b_k + n * (b_k - 1);
            for (j = 1; j <= nr; j++) {
              v[j - 1] = h[(j + hoffset) - 1];
            }
          }

          tst = v[0];
          rt2r = xzlarfg(nr, &tst, v);
          v[0] = tst;
          if (b_k + 1 > m + 1) {
            h[b_k + n * (b_k - 1)] = tst;
            h[(b_k + n * (b_k - 1)) + 1] = 0.0;
            if (b_k + 1 < i) {
              h[(b_k + n * (b_k - 1)) + 2] = 0.0;
            }
          } else {
            if (m + 1 > k + 1) {
              h[b_k + n * (b_k - 1)] *= 1.0 - rt2r;
            }
          }

          tst = v[1];
          htmp1 = rt2r * v[1];
          if (nr == 3) {
            ab = v[2];
            ba = rt2r * v[2];
            for (j = b_k; j + 1 < n+1; j++) {
              aa = (h[b_k + n * j] + tst * h[(b_k + n * j) + 1]) + ab * h[(b_k
                + n * j) + 2];
              h[b_k + n * j] -= aa * rt2r;
              h[(b_k + n * j) + 1] -= aa * htmp1;
              h[(b_k + n * j) + 2] -= aa * ba;
            }

            if (b_k + 4 < i + 1) {
              nr = b_k;
            } else {
              nr = i - 3;
            }

            for (j = 0; j + 1 <= nr + 4; j++) {
              aa = (h[j + n * b_k] + tst * h[j + n * (b_k + 1)]) + ab * h[j +
                n * (b_k + 2)];
              h[j + n * b_k] -= aa * rt2r;
              h[j + n * (b_k + 1)] -= aa * htmp1;
              h[j + n * (b_k + 2)] -= aa * ba;
            }
          } else {
            if (nr == 2) {
              for (j = b_k; j + 1 < n+1; j++) {
                aa = h[b_k + n * j] + tst * h[(b_k + n * j) + 1];
                h[b_k + n * j] -= aa * rt2r;
                h[(b_k + n * j) + 1] -= aa * htmp1;
              }

              for (j = 0; j + 1 <= i + 1; j++) {
                aa = h[j + n * b_k] + tst * h[j + n * (b_k + 1)];
                h[j + n * b_k] -= aa * rt2r;
                h[j + n * (b_k + 1)] -= aa * htmp1;
              }
            }
          }
        }

        its++;
      }
    }

    if (!goto150) {
      info = i + 1;
      exitg1 = true;
    } else {
      if ((L != i + 1) && (L == i)) {
        tst = h[(i + n * i) - 1];
        htmp1 = h[i + n * (i - 1)];
        aa = h[i + n * i];
        xdlanv2(&h[(i + n * (i - 1)) - 1], &tst, &htmp1, &aa, &ab, &ba, &rt2r,
                &rt1r, &s, &sn);
        h[(i + n * i) - 1] = tst;
        h[i + n * (i - 1)] = htmp1;
        h[i + n * i] = aa;
        if (n > i + 1) {
          nr = (i + (i + 1) * n) - 1;
          hoffset = i + (i + 1) * n;
          for (k = 1; k <= n-1 - i; k++) {
            tst = s * h[nr] + sn * h[hoffset];
            h[hoffset] = s * h[hoffset] - sn * h[nr];
            h[nr] = tst;
            hoffset += n;
            nr += n;
          }
        }

        if (!(i - 1 < 1)) {
          nr = (i - 1) * n;
          hoffset = i * n;
          for (k = 1; k < i; k++) {
            tst = s * h[nr] + sn * h[hoffset];
            h[hoffset] = s * h[hoffset] - sn * h[nr];
            h[nr] = tst;
            hoffset++;
            nr++;
          }
        }
      }

      i = L - 2;
    }
  }

  return info;
}

//
// File trailer for xdhseqr.cpp
//
// [EOF]
//
