// --------------------------------------------------------------------------
//                             SSEARCH ROUTINE
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
// Integer ambiguity vector search by employing the search-and-shrink
// technique.
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
#include "..\mclambda_emxutil.h"
#include "..\mclambda_rtwutil.h"
#include "..\math_functions\sort1.cpp"
#include "..\routines\ssearch.h"

using namespace std;
// --------------------------------------------------------------------------
//                          Function Definitions
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
//                                ssearch
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
void ssearch(int n1, const double ahat[], const double L[], const double D[],
             double ncands, emxArray_real_T *afixed, emxArray_real_T *sqnorm)
{
  // ============================= VARIABLES ================================
  int n;
  int loop_ub;
  double Chi2;
  double dist[n1];
  double acond[n1];
  double zcond[n1];
  double step[n1];
  boolean_T endsearch;
  unsigned int count;
  double left;
  double newdist;
  double imax;
  double S[n1*n1];
  double k;
  emxArray_int32_T *iidx;
  emxArray_real_T *b_afixed;
  int ix;
  double S_data[n1-1];

  // ============================ START PROGRAM ===============================
  // Initialization
  n = afixed->size[0] * afixed->size[1];
  afixed->size[0] = n1;
  afixed->size[1] = (int)ncands;
  emxEnsureCapacity((emxArray__common *)afixed, n, sizeof(double));
  loop_ub = n1 * (int)ncands;
  for (n = 0; n < loop_ub; n++) {
    afixed->data[n] = 0.0;
  }

  n = sqnorm->size[0] * sqnorm->size[1];
  sqnorm->size[0] = 1;
  sqnorm->size[1] = (int)ncands;
  emxEnsureCapacity((emxArray__common *)sqnorm, n, sizeof(double));
  loop_ub = (int)ncands;
  for (n = 0; n < loop_ub; n++) {
    sqnorm->data[n] = 0.0;
  }

  // Initializing the variables for searching
  Chi2 = 1.0E+18;

  // Start search with an infinite chi^2
  memset(&dist[0], 0, n1 * sizeof(double));
  memset(&acond[0], 0, n1 * sizeof(double));
  memset(&zcond[0], 0, n1 * sizeof(double));
  memset(&step[0], 0, n1 * sizeof(double));
  dist[n1-1] = 0.0;

  endsearch = false;
  count = 0U;

  // The number of candidates
  acond[n1-1] = ahat[n1-1];
  zcond[n1-1] = rt_roundd_snf(ahat[n1-1]);
  left = ahat[n1-1] - zcond[n1-1];
  newdist = left;
  if (left < 0.0) {
    newdist = -1.0;
  } else if (left > 0.0) {
    newdist = 1.0;
  } else {
    if (left == 0.0) {
      newdist = 0.0;
    }
  }

  step[n1-1] = newdist;

  if (newdist == 0.0) {
    step[n1-1] = 1.0;
  }

  imax = ncands;

  // Initially, the maximum F(z) is at ncands
  memset(&S[0], 0, n1*n1 * sizeof(double));

  // Used to compute conditional ambiguities
  k = (double)(n1);

  // Start the main search-loop
  while (!endsearch) {
    newdist = dist[(int)k - 1] + left * left / D[(int)k - 1];
    if (newdist < Chi2) {
      if (k != 1.0) {
        // Case 1: move down
        k--;
        dist[(int)k - 1] = newdist;
        newdist = zcond[(int)(k + 1.0) - 1] - acond[(int)(k + 1.0) - 1];
        loop_ub = (int)k;
        for (n = 0; n < loop_ub; n++) {
          S_data[n] = S[((int)(k + 1.0) + n1 * n) - 1] + newdist * L[((int)(k +
            1.0) + n1 * n) - 1];
        }

        loop_ub = (int)k;
        for (n = 0; n < loop_ub; n++) {
          S[((int)k + n1 * n) - 1] = S_data[n];
        }

        acond[(int)k - 1] = ahat[(int)k - 1] + S[((int)k + n1 * ((int)k - 1)) -
          1];
        zcond[(int)k - 1] = rt_roundd_snf(acond[(int)k - 1]);
        left = acond[(int)k - 1] - zcond[(int)k - 1];
        newdist = left;
        if (left < 0.0) {
          newdist = -1.0;
        } else if (left > 0.0) {
          newdist = 1.0;
        } else {
          if (left == 0.0) {
            newdist = 0.0;
          }
        }

        step[(int)k - 1] = newdist;

        if (step[(int)k - 1] == 0.0) {
          step[(int)k - 1] = 1.0;
        }

      } else {
        // Case 2: store the found candidate and try next valid integer
        if (count < ncands - 1.0) {
          // Store the first ncands-1 initial points as candidates
          count++;
          for (n = 0; n < n1; n++) {
            afixed->data[n + afixed->size[0] * ((int)count - 1)] = zcond[n];
          }

          sqnorm->data[(int)count - 1] = newdist;

          // Store F(zcond)
        } else {
          for (n = 0; n < n1; n++) {
            afixed->data[n + afixed->size[0] * ((int)imax - 1)] = zcond[n];
          }

          sqnorm->data[(int)imax - 1] = newdist;
          n = sqnorm->size[1];
          Chi2 = sqnorm->data[0];
          loop_ub = 0;
          if (sqnorm->size[1] > 1) {
            for (ix = 1; ix + 1 <= n; ix++) {
              if (sqnorm->data[ix] > Chi2) {
                Chi2 = sqnorm->data[ix];
                loop_ub = ix;
              }
            }
          }

          imax = loop_ub + 1;
        }

        zcond[0] += step[0];

        // Next valid integer
        left = acond[0] - zcond[0];
        newdist = step[0];
        if (step[0] < 0.0) {
          newdist = -1.0;
        } else if (step[0] > 0.0) {
          newdist = 1.0;
        } else {
          if (step[0] == 0.0) {
            newdist = 0.0;
          }
        }

        step[0] = -step[0] - newdist;
      }
    } else {
      // Case 3: exit or move up
      if (k == (double)(n1)) {
        endsearch = true;
      } else {
        k++;

        // Move up
        zcond[(int)k - 1] += step[(int)k - 1];

        // Next valid integer
        left = acond[(int)k - 1] - zcond[(int)k - 1];
        newdist = step[(int)k - 1];
        if (step[(int)k - 1] < 0.0) {
          newdist = -1.0;
        } else if (step[(int)k - 1] > 0.0) {
          newdist = 1.0;
        } else {
          if (step[(int)k - 1] == 0.0) {
            newdist = 0.0;
          }
        }

        step[(int)k - 1] = -step[(int)k - 1] - newdist;
      }
    }
  }

  emxInit_int32_T(&iidx, 2);
  emxInit_real_T(&b_afixed, 2);
  sort(sqnorm, iidx);
  n = b_afixed->size[0] * b_afixed->size[1];
  b_afixed->size[0] = n1;
  b_afixed->size[1] = iidx->size[1];
  emxEnsureCapacity((emxArray__common *)b_afixed, n, sizeof(double));
  loop_ub = iidx->size[1];
  for (n = 0; n < loop_ub; n++) {
    for (ix = 0; ix < n1; ix++) {
      b_afixed->data[ix + b_afixed->size[0] * n] = afixed->data[ix +
        afixed->size[0] * (iidx->data[iidx->size[0] * n] - 1)];
    }
  }

  emxFree_int32_T(&iidx);
  n = afixed->size[0] * afixed->size[1];
  afixed->size[0] = n1;
  afixed->size[1] = b_afixed->size[1];
  emxEnsureCapacity((emxArray__common *)afixed, n, sizeof(double));
  loop_ub = b_afixed->size[1];
  for (n = 0; n < loop_ub; n++) {
    for (ix = 0; ix < n1; ix++) {
      afixed->data[ix + afixed->size[0] * n] = b_afixed->data[ix +
        b_afixed->size[0] * n];
    }
  }

  emxFree_real_T(&b_afixed);
}
// --------------------------------------------------------------------------
//                            End of ssearch
// --------------------------------------------------------------------------