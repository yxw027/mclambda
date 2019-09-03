// --------------------------------------------------------------------------
//                             LSEARCH ROUTINE
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
//  Release date  : AUG-2019
//  Authors       : Hernández Olcina, Jorge
//
//  Master in Geomatics Engineering and Geoinformation
//  Universidad Politécnica de Valencia
//  Hochschule Karlsruhe - Technik und Wirtschaft University of Applied Sciences
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
//
//  DESCRIPTION:
//
//  This routine finds the integer vector which is closest to a given
//  float vector, in a least squares sence. This is the search-step in
//  integer ambiguity resolution. It is best to perform this search only
//  on ambiguities which have been decorrelated using MCLAMBDA.
//
//  INPUTS:
//
//     n      : Number of float ambiguities
//     ahat   : Float ambiguities (should be decorrelated for computational
//              efficiency)
//     L,D    : LtDL-decomposition of the variance-covariance matrix of the
//              float ambiguities ahat
//     ncands : Number of requested candidates
//
//  OUTPUTS:
//
//     afixed : Estimated integers
//     sqnorm : Corresponding squared norms
//
// --------------------------------------------------------------------------

// Include Files
#include "..\math_functions\sortLE.cpp"
#include "..\math_functions\inv.cpp"
#include "chistart.cpp"
#include "lsearch.h"
#include "..\mclambda_rtwutil.h"

using namespace std;
// --------------------------------------------------------------------------
//                          Function Definitions
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
//                               lsearch
// --------------------------------------------------------------------------
//
// Arguments    : int n                   -> Number of float ambiguities
//                const double ahat       -> Vector of float ambiguities
//                const double L          -> Matrix D of LDL-decomposition
//                const double D          -> Matrix L of LDL-decomposition
//                double ncands           -> Requested number of candidates
//                emxArray_real_T *afixed -> Output afixed
//                emxArray_real_T *sqnorm -> Output sqnorm
//
// Return       : void
//
// --------------------------------------------------------------------------
void lsearch(int n1, const double ahat[], const double L[], const double D[],
             double ncands, emxArray_real_T *afixed, emxArray_real_T *sqnorm)
{
  // ============================= VARIABLES ================================
  double Chi2;
  double Linv[n1*n1];
  int i;
  double right[n1+1];
  double Dinv[n1];
  double left[n1+1];
  int qEnd;
  double dq[n1];
  int endsearch;
  unsigned int ncan;
  int iold;
  int i2;
  double lef[n1];
  double distl[n1];
  double endd[n1];
  emxArray_real_T *tmp;
  emxArray_real_T *ndx;
  int j;
  double reach;
  emxArray_real_T *varargin_2;
  int pEnd;
  int p;
  boolean_T empty_non_axis_sizes;
  int n;
  double mtmp;
  boolean_T exitg1;
  int col_size[2];
  int k;
  emxArray_int32_T *idx;
  int col_data[n1+1];
  emxArray_int32_T *iwork;
  emxArray_real_T *ycol;

  // ============================ START PROGRAM ===============================
  // Computes the initial size of the search ellipsoid (USING CHISTART ROUTINE)
  Chi2 = chistart(n1, D, L, ahat, ncands+1);
  // Initialization
  inv(n1, L, Linv);
  for (i = 0; i < n1; i++) {
    Dinv[i] = 1.0 / D[i];
    right[i] = 0.0;
  }

  right[n1] = Chi2;
  memset(&left[0], 0, n1+1 * sizeof(double));
  for (qEnd = 0; qEnd < n1-1; qEnd++) {
    dq[qEnd] = Dinv[1 + qEnd] / Dinv[qEnd];
  }

  dq[n1-1] = 1.0 / Dinv[n1-1];
  endsearch = 0;
  ncan = 0U;
  i = n1;
  iold = n1;
  qEnd = afixed->size[0] * afixed->size[1];
  afixed->size[0] = n1;
  afixed->size[1] = (int)ncands;
  emxEnsureCapacity((emxArray__common *)afixed, qEnd, sizeof(double));
  i2 = n1 * (int)ncands;
  for (qEnd = 0; qEnd < i2; qEnd++) {
    afixed->data[qEnd] = 0.0;
  }

  qEnd = sqnorm->size[0] * sqnorm->size[1];
  sqnorm->size[0] = 1;
  sqnorm->size[1] = (int)ncands;
  emxEnsureCapacity((emxArray__common *)sqnorm, qEnd, sizeof(double));
  i2 = (int)ncands;
  for (qEnd = 0; qEnd < i2; qEnd++) {
    sqnorm->data[qEnd] = 0.0;
  }

  // Start the main search-loop
  memset(&lef[0], 0, n1 * sizeof(double));
  memset(&distl[0], 0, n1 * sizeof(double));
  memset(&endd[0], 0, n1 * sizeof(double));
  while (!(endsearch != 0)) {
    i--;
    if (iold + 1 <= i + 1) {
      lef[i] += Linv[(i + n1 * i) + 1];
    } else {
      lef[i] = 0.0;
      for (j = 0; j <= 10 - i; j++) {
        i2 = (i + j) + 1;
        lef[i] += Linv[i2 + n1 * i] * distl[i2];
      }
    }

    iold = i;
    right[i] = (right[i + 1] - left[i + 1]) * dq[i];
    reach = std::sqrt(right[i]);
    distl[i] = std::ceil((ahat[i] - reach) - lef[i]) - ahat[i];
    if (distl[i] > reach - lef[i]) {
      // There is nothing at this level, so backtrack
      pEnd = 0;
      p = 0;
      while ((!(p != 0)) && (i + 1 < n1)) {
        i++;
        if (distl[i] < endd[i]) {
          distl[i]++;
          reach = distl[i] + lef[i];
          left[i] = reach * reach;
          p = 1;
          if (i + 1 == n1) {
            pEnd = 1;
          }
        }
      }

      if ((i + 1 == n1) && (!(pEnd != 0))) {
        endsearch = 1;
      }
    } else {
      // Set the right border
      endd[i] = (reach - lef[i]) - 1.0;
      reach = distl[i] + lef[i];
      left[i] = reach * reach;
    }

    if (i + 1 == 1) {
      // Collect the integer vectors and corresponding
      // squared distances, add to vectors "afixed" and "sqnorm" if:
      //   - Less then "ncands" candidates found so far
      //   - The squared norm is smaller than one of the previous ones
      reach = Chi2 - (right[0] - left[0]) * Dinv[0];
      endd[0]++;
      while (distl[0] <= endd[0]) {
        if (ncan < ncands) {
          ncan++;
          for (qEnd = 0; qEnd < n1; qEnd++) {
            afixed->data[qEnd + afixed->size[0] * ((int)ncan - 1)] = distl[qEnd]
              + ahat[qEnd];
          }

          sqnorm->data[(int)ncan - 1] = reach;
        } else {
          i2 = 1;
          n = sqnorm->size[1];
          mtmp = sqnorm->data[0];
          p = 0;
          if (sqnorm->size[1] > 1) {
            if (rtIsNaN(sqnorm->data[0])) {
              pEnd = 1;
              exitg1 = false;
              while ((!exitg1) && (pEnd + 1 <= n)) {
                i2 = pEnd + 1;
                if (!rtIsNaN(sqnorm->data[pEnd])) {
                  mtmp = sqnorm->data[pEnd];
                  p = pEnd;
                  exitg1 = true;
                } else {
                  pEnd++;
                }
              }
            }

            if (i2 < sqnorm->size[1]) {
              while (i2 + 1 <= n) {
                if (sqnorm->data[i2] > mtmp) {
                  mtmp = sqnorm->data[i2];
                  p = i2;
                }

                i2++;
              }
            }
          }

          if (reach < mtmp) {
            for (qEnd = 0; qEnd < n1; qEnd++) {
              afixed->data[qEnd + afixed->size[0] * p] = distl[qEnd] + ahat[qEnd];
            }

            sqnorm->data[p] = reach;
          }
        }

        reach += (2.0 * (distl[0] + lef[0]) + 1.0) * Dinv[0];
        distl[0]++;
      }

      // And backtrack
      pEnd = 0;
      p = 0;
      while ((!(p != 0)) && (i + 1 < n1)) {
        i++;
        if (distl[i] < endd[i]) {
          distl[i]++;
          reach = distl[i] + lef[i];
          left[i] = reach * reach;
          p = 1;
          if (i + 1 == n1) {
            pEnd = 1;
          }
        }
      }

      if ((i + 1 == n1) && (!(pEnd != 0))) {
        endsearch = 1;
      }
    }
  }

  emxInit_real_T(&tmp, 2);
  emxInit_real_T1(&ndx, 1);

  // Sort the resulting candidates, according to the norm
  qEnd = ndx->size[0];
  ndx->size[0] = sqnorm->size[1];
  emxEnsureCapacity((emxArray__common *)ndx, qEnd, sizeof(double));
  i2 = sqnorm->size[1];
  for (qEnd = 0; qEnd < i2; qEnd++) {
    ndx->data[qEnd] = sqnorm->data[sqnorm->size[0] * qEnd];
  }

  emxInit_real_T(&varargin_2, 2);
  qEnd = varargin_2->size[0] * varargin_2->size[1];
  varargin_2->size[0] = afixed->size[1];
  varargin_2->size[1] = n1;
  emxEnsureCapacity((emxArray__common *)varargin_2, qEnd, sizeof(double));
  for (qEnd = 0; qEnd < n1; qEnd++) {
    i2 = afixed->size[1];
    for (iold = 0; iold < i2; iold++) {
      varargin_2->data[iold + varargin_2->size[0] * qEnd] = afixed->data[qEnd +
        afixed->size[0] * iold];
    }
  }

  if (!(ndx->size[0] == 0)) {
    i2 = ndx->size[0];
  } else if (!(varargin_2->size[0] == 0)) {
    i2 = varargin_2->size[0];
  } else {
    i2 = 0;
  }

  empty_non_axis_sizes = (i2 == 0);
  if (empty_non_axis_sizes || (!(ndx->size[0] == 0))) {
    pEnd = 1;
  } else {
    pEnd = 0;
  }

  if (empty_non_axis_sizes || (!(varargin_2->size[0] == 0))) {
    p = n1;
  } else {
    p = 0;
  }

  qEnd = tmp->size[0] * tmp->size[1];
  tmp->size[0] = i2;
  tmp->size[1] = pEnd + p;
  emxEnsureCapacity((emxArray__common *)tmp, qEnd, sizeof(double));
  for (qEnd = 0; qEnd < pEnd; qEnd++) {
    for (iold = 0; iold < i2; iold++) {
      tmp->data[iold + tmp->size[0] * qEnd] = ndx->data[iold + i2 * qEnd];
    }
  }

  for (qEnd = 0; qEnd < p; qEnd++) {
    for (iold = 0; iold < i2; iold++) {
      tmp->data[iold + tmp->size[0] * (qEnd + pEnd)] = varargin_2->data[iold +
        i2 * qEnd];
    }
  }

  col_size[0] = 1;
  col_size[1] = tmp->size[1];
  k = 1;
  emxFree_real_T(&ndx);
  while (k <= tmp->size[1]) {
    col_data[k - 1] = k;
    k++;
  }

  emxInit_int32_T1(&idx, 1);
  n = tmp->size[0];
  qEnd = idx->size[0];
  idx->size[0] = tmp->size[0];
  emxEnsureCapacity((emxArray__common *)idx, qEnd, sizeof(int));
  i2 = tmp->size[0];
  for (qEnd = 0; qEnd < i2; qEnd++) {
    idx->data[qEnd] = 0;
  }

  if ((tmp->size[0] == 0) || (tmp->size[1] == 0)) {
    for (k = 1; k <= n; k++) {
      idx->data[k - 1] = k;
    }
  } else {
    emxInit_int32_T1(&iwork, 1);
    qEnd = iwork->size[0];
    iwork->size[0] = tmp->size[0];
    emxEnsureCapacity((emxArray__common *)iwork, qEnd, sizeof(int));
    for (k = 1; k <= n - 1; k += 2) {
      if (sortLE(tmp, col_data, col_size, k, k + 1)) {
        idx->data[k - 1] = k;
        idx->data[k] = k + 1;
      } else {
        idx->data[k - 1] = k + 1;
        idx->data[k] = k;
      }
    }

    if ((tmp->size[0] & 1) != 0) {
      idx->data[tmp->size[0] - 1] = tmp->size[0];
    }

    i = 2;
    while (i < n) {
      i2 = i << 1;
      j = 1;
      for (pEnd = 1 + i; pEnd < n + 1; pEnd = qEnd + i) {
        p = j;
        iold = pEnd;
        qEnd = j + i2;
        if (qEnd > n + 1) {
          qEnd = n + 1;
        }

        k = 0;
        endsearch = qEnd - j;
        while (k + 1 <= endsearch) {
          if (sortLE(tmp, col_data, col_size, idx->data[p - 1], idx->data[iold -
                     1])) {
            iwork->data[k] = idx->data[p - 1];
            p++;
            if (p == pEnd) {
              while (iold < qEnd) {
                k++;
                iwork->data[k] = idx->data[iold - 1];
                iold++;
              }
            }
          } else {
            iwork->data[k] = idx->data[iold - 1];
            iold++;
            if (iold == qEnd) {
              while (p < pEnd) {
                k++;
                iwork->data[k] = idx->data[p - 1];
                p++;
              }
            }
          }

          k++;
        }

        for (k = 0; k + 1 <= endsearch; k++) {
          idx->data[(j + k) - 1] = iwork->data[k];
        }

        j = qEnd;
      }

      i = i2;
    }

    emxFree_int32_T(&iwork);
  }

  emxInit_real_T1(&ycol, 1);
  i2 = tmp->size[0];
  n = tmp->size[1];
  ncan = (unsigned int)tmp->size[0];
  qEnd = ycol->size[0];
  ycol->size[0] = (int)ncan;
  emxEnsureCapacity((emxArray__common *)ycol, qEnd, sizeof(double));
  for (j = 0; j + 1 <= n; j++) {
    for (i = 0; i + 1 <= i2; i++) {
      ycol->data[i] = tmp->data[(idx->data[i] + tmp->size[0] * j) - 1];
    }

    for (i = 0; i + 1 <= i2; i++) {
      tmp->data[i + tmp->size[0] * j] = ycol->data[i];
    }
  }

  emxFree_real_T(&ycol);
  emxFree_int32_T(&idx);
  i2 = tmp->size[0];
  qEnd = sqnorm->size[0] * sqnorm->size[1];
  sqnorm->size[0] = 1;
  sqnorm->size[1] = i2;
  emxEnsureCapacity((emxArray__common *)sqnorm, qEnd, sizeof(double));
  for (qEnd = 0; qEnd < i2; qEnd++) {
    sqnorm->data[sqnorm->size[0] * qEnd] = tmp->data[qEnd];
  }

  i2 = tmp->size[0];
  qEnd = varargin_2->size[0] * varargin_2->size[1];
  varargin_2->size[0] = i2;
  varargin_2->size[1] = n1;
  emxEnsureCapacity((emxArray__common *)varargin_2, qEnd, sizeof(double));
  for (qEnd = 0; qEnd < n1; qEnd++) {
    for (iold = 0; iold < i2; iold++) {
      varargin_2->data[iold + varargin_2->size[0] * qEnd] = tmp->data[iold +
        tmp->size[0] * (1 + qEnd)];
    }
  }

  qEnd = tmp->size[0];
  i2 = qEnd * n1;
  k = 0;
  emxFree_real_T(&tmp);
  while (k + 1 <= i2) {
    varargin_2->data[k] = rt_roundd_snf(varargin_2->data[k]);
    k++;
  }

  qEnd = afixed->size[0] * afixed->size[1];
  afixed->size[0] = n1;
  afixed->size[1] = varargin_2->size[0];
  emxEnsureCapacity((emxArray__common *)afixed, qEnd, sizeof(double));
  i2 = varargin_2->size[0];
  for (qEnd = 0; qEnd < i2; qEnd++) {
    for (iold = 0; iold < n1; iold++) {
      afixed->data[iold + afixed->size[0] * qEnd] = varargin_2->data[qEnd +
        varargin_2->size[0] * iold];
    }
  }

  emxFree_real_T(&varargin_2);
}
// --------------------------------------------------------------------------
//                              End of lsearch
// --------------------------------------------------------------------------