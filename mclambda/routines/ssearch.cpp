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
#include "..\LAMBDA_emxutil.h"
#include "..\LAMBDA_rtwutil.h"
#include "..\math_functions\sort1.cpp"
#include "..\routines\ssearch.h"

using namespace std;
// --------------------------------------------------------------------------
//                          Function Definitions
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
//                               b_ssearch
// --------------------------------------------------------------------------
//
// Arguments    : const double ahat_data  -> Vector of float ambiguities
//                const double ahat_size  -> Size of float ambiguities vector
//                const double L_data     -> Matrix L of LDL-decomposition
//                const double L_size     -> Size of matrix D
//                const double D_data     -> Matrix D of LDL-decomposition
//                double ncands           -> Requested number of candidates
//                emxArray_real_T *afixed -> Output afixed
//                emxArray_real_T *sqnorm -> Output sqnorm
//
// Return       : void
//
// --------------------------------------------------------------------------
void b_ssearch(const double ahat_data[], const int ahat_size[1], const double
               L_data[], const int L_size[2], const double D_data[], double
               ncands, emxArray_real_T *afixed, emxArray_real_T *sqnorm)
{
  // ============================= VARIABLES ================================
  int n;
  int ix;
  int loop_ub;
  double Chi2;
  double dist_data[12];
  boolean_T endsearch;
  unsigned int count;
  double acond_data[12];
  double zcond_data[12];
  double left;
  double newdist;
  double step_data[12];
  double imax;
  int S_size_idx_0;
  double S_data[144];
  int b_n;
  int k;
  int itmp;
  emxArray_int32_T *iidx;
  emxArray_real_T *b_afixed;
  signed char i3;
  signed char tmp_data[12];
  double b_S_data[11];
  double b_zcond_data[12];

  // ============================ START PROGRAM ===============================
  // Initialization
  n = ahat_size[0];
  ix = afixed->size[0] * afixed->size[1];
  afixed->size[0] = ahat_size[0];
  afixed->size[1] = (int)ncands;
  emxEnsureCapacity((emxArray__common *)afixed, ix, sizeof(double));
  loop_ub = ahat_size[0] * (int)ncands;
  for (ix = 0; ix < loop_ub; ix++) {
    afixed->data[ix] = 0.0;
  }

  ix = sqnorm->size[0] * sqnorm->size[1];
  sqnorm->size[0] = 1;
  sqnorm->size[1] = (int)ncands;
  emxEnsureCapacity((emxArray__common *)sqnorm, ix, sizeof(double));
  loop_ub = (int)ncands;
  for (ix = 0; ix < loop_ub; ix++) {
    sqnorm->data[ix] = 0.0;
  }

  // Initializing the variables for searching
  Chi2 = 1.0E+18;

  // Start search with an infinite chi^2
  loop_ub = ahat_size[0];
  for (ix = 0; ix < loop_ub; ix++) {
    dist_data[ix] = 0.0;
  }

  dist_data[ahat_size[0] - 1] = 0.0;

  endsearch = false;
  count = 0U;

  // The number of candidates
  loop_ub = ahat_size[0];
  for (ix = 0; ix < loop_ub; ix++) {
    acond_data[ix] = 0.0;
  }

  acond_data[ahat_size[0] - 1] = ahat_data[ahat_size[0] - 1];
  loop_ub = ahat_size[0];
  for (ix = 0; ix < loop_ub; ix++) {
    zcond_data[ix] = 0.0;
  }

  zcond_data[ahat_size[0] - 1] = rt_roundd_snf(acond_data[ahat_size[0] - 1]);
  left = acond_data[ahat_size[0] - 1] - zcond_data[ahat_size[0] - 1];
  loop_ub = ahat_size[0];
  for (ix = 0; ix < loop_ub; ix++) {
    step_data[ix] = 0.0;
  }

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

  step_data[ahat_size[0] - 1] = newdist;

  if (step_data[ahat_size[0] - 1] == 0.0) {
    step_data[ahat_size[0] - 1] = 1.0;
  }

  imax = ncands;

  // Initially, the maximum F(z) is at ncands
  S_size_idx_0 = ahat_size[0];
  loop_ub = ahat_size[0] * ahat_size[0];
  for (ix = 0; ix < loop_ub; ix++) {
    S_data[ix] = 0.0;
  }

  loop_ub = ahat_size[0];
  for (ix = 0; ix < loop_ub; ix++) {
    b_n = ahat_size[0];
    for (itmp = 0; itmp < b_n; itmp++) {
      S_data[itmp + S_size_idx_0 * ix] = 0.0;
    }
  }

  // Used to compute conditional ambiguities
  k = ahat_size[0] - 1;

  // Start the main search-loop
  while (!endsearch) {
    newdist = dist_data[k] + left * left / D_data[k];
    if (newdist < Chi2) {
      if (k + 1 != 1) {
        // Case 1: move down
        k--;
        dist_data[k] = newdist;
        newdist = zcond_data[k + 1] - acond_data[k + 1];
        for (ix = 0; ix <= k; ix++) {
          b_S_data[ix] = S_data[(k + S_size_idx_0 * ix) + 1] + newdist * L_data
            [(k + L_size[0] * ix) + 1];
        }

        loop_ub = k + 1;
        for (ix = 0; ix < loop_ub; ix++) {
          S_data[k + S_size_idx_0 * ix] = b_S_data[ix];
        }

        acond_data[k] = ahat_data[k] + S_data[k + S_size_idx_0 * k];
        zcond_data[k] = rt_roundd_snf(acond_data[k]);
        left = acond_data[k] - zcond_data[k];
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

        step_data[k] = newdist;

        if (step_data[k] == 0.0) {
          step_data[k] = 1.0;
        }
      } else {
        // Case 2: store the found candidate and try next valid integer
        if (count < ncands - 1.0) {
          // Store the first ncands-1 initial points as candidates
          count++;
          i3 = (signed char)((signed char)afixed->size[0] - 1);
          loop_ub = i3;
          for (ix = 0; ix <= loop_ub; ix++) {
            tmp_data[ix] = (signed char)ix;
          }

          for (ix = 0; ix < n; ix++) {
            b_zcond_data[ix] = zcond_data[ix];
          }

          loop_ub = i3 + 1;
          for (ix = 0; ix < loop_ub; ix++) {
            afixed->data[tmp_data[ix] + afixed->size[0] * ((int)count - 1)] =
              b_zcond_data[ix];
          }

          sqnorm->data[(int)count - 1] = newdist;

          // Store F(zcond)
        } else {
          i3 = (signed char)((signed char)afixed->size[0] - 1);
          loop_ub = i3;
          for (ix = 0; ix <= loop_ub; ix++) {
            tmp_data[ix] = (signed char)ix;
          }

          for (ix = 0; ix < n; ix++) {
            b_zcond_data[ix] = zcond_data[ix];
          }

          loop_ub = i3 + 1;
          for (ix = 0; ix < loop_ub; ix++) {
            afixed->data[tmp_data[ix] + afixed->size[0] * ((int)imax - 1)] =
              b_zcond_data[ix];
          }

          sqnorm->data[(int)imax - 1] = newdist;
          b_n = sqnorm->size[1];
          Chi2 = sqnorm->data[0];
          itmp = 0;
          if (sqnorm->size[1] > 1) {
            for (ix = 1; ix + 1 <= b_n; ix++) {
              if (sqnorm->data[ix] > Chi2) {
                Chi2 = sqnorm->data[ix];
                itmp = ix;
              }
            }
          }

          imax = itmp + 1;
        }

        zcond_data[0] += step_data[0];

        // Next valid integer
        left = acond_data[0] - zcond_data[0];
        newdist = step_data[0];
        if (step_data[0] < 0.0) {
          newdist = -1.0;
        } else if (step_data[0] > 0.0) {
          newdist = 1.0;
        } else {
          if (step_data[0] == 0.0) {
            newdist = 0.0;
          }
        }

        step_data[0] = -step_data[0] - newdist;
      }
    } else {
      // Case 3: exit or move up
      if (k + 1 == n) {
        endsearch = true;
      } else {
        k++;

        // Move up
        zcond_data[k] += step_data[k];

        // Next valid integer
        left = acond_data[k] - zcond_data[k];
        newdist = step_data[k];
        if (step_data[k] < 0.0) {
          newdist = -1.0;
        } else if (step_data[k] > 0.0) {
          newdist = 1.0;
        } else {
          if (step_data[k] == 0.0) {
            newdist = 0.0;
          }
        }

        step_data[k] = -step_data[k] - newdist;
      }
    }
  }

  emxInit_int32_T(&iidx, 2);
  emxInit_real_T(&b_afixed, 2);
  sort(sqnorm, iidx);
  b_n = afixed->size[0];
  ix = b_afixed->size[0] * b_afixed->size[1];
  b_afixed->size[0] = b_n;
  b_afixed->size[1] = iidx->size[1];
  emxEnsureCapacity((emxArray__common *)b_afixed, ix, sizeof(double));
  loop_ub = iidx->size[1];
  for (ix = 0; ix < loop_ub; ix++) {
    for (itmp = 0; itmp < b_n; itmp++) {
      b_afixed->data[itmp + b_afixed->size[0] * ix] = afixed->data[itmp +
        afixed->size[0] * (iidx->data[iidx->size[0] * ix] - 1)];
    }
  }

  emxFree_int32_T(&iidx);
  ix = afixed->size[0] * afixed->size[1];
  afixed->size[0] = b_afixed->size[0];
  afixed->size[1] = b_afixed->size[1];
  emxEnsureCapacity((emxArray__common *)afixed, ix, sizeof(double));
  loop_ub = b_afixed->size[1];
  for (ix = 0; ix < loop_ub; ix++) {
    b_n = b_afixed->size[0];
    for (itmp = 0; itmp < b_n; itmp++) {
      afixed->data[itmp + afixed->size[0] * ix] = b_afixed->data[itmp +
        b_afixed->size[0] * ix];
    }
  }

  emxFree_real_T(&b_afixed);
}
// --------------------------------------------------------------------------
//                            End of b_ssearch
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
//                                ssearch
// --------------------------------------------------------------------------
//
// Arguments    : const double ahat       -> Vector of float ambiguities
//                const double L          -> Matrix D of LDL-decomposition
//                const double D          -> Matrix L of LDL-decomposition
//                double ncands           -> Requested number of candidates
//                emxArray_real_T *afixed -> Output afixed
//                emxArray_real_T *sqnorm -> Output sqnorm
//
// Return       : void
//
// --------------------------------------------------------------------------
void ssearch(const double ahat[12], const double L[144], const double D[12],
             double ncands, emxArray_real_T *afixed, emxArray_real_T *sqnorm)
{
  // ============================= VARIABLES ================================
  int n;
  int loop_ub;
  double Chi2;
  double dist[12];
  double acond[12];
  double zcond[12];
  double step[12];
  boolean_T endsearch;
  unsigned int count;
  double left;
  double newdist;
  double imax;
  double S[144];
  double k;
  emxArray_int32_T *iidx;
  emxArray_real_T *b_afixed;
  int ix;
  double S_data[11];

  // ============================ START PROGRAM ===============================
  // Initialization
  n = afixed->size[0] * afixed->size[1];
  afixed->size[0] = 12;
  afixed->size[1] = (int)ncands;
  emxEnsureCapacity((emxArray__common *)afixed, n, sizeof(double));
  loop_ub = 12 * (int)ncands;
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
  memset(&dist[0], 0, 12U * sizeof(double));
  memset(&acond[0], 0, 12U * sizeof(double));
  memset(&zcond[0], 0, 12U * sizeof(double));
  memset(&step[0], 0, 12U * sizeof(double));
  dist[11] = 0.0;

  endsearch = false;
  count = 0U;

  // The number of candidates
  acond[11] = ahat[11];
  zcond[11] = rt_roundd_snf(ahat[11]);
  left = ahat[11] - zcond[11];
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

  step[11] = newdist;

  if (newdist == 0.0) {
    step[11] = 1.0;
  }

  imax = ncands;

  // Initially, the maximum F(z) is at ncands
  memset(&S[0], 0, 144U * sizeof(double));

  // Used to compute conditional ambiguities
  k = 12.0;

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
          S_data[n] = S[((int)(k + 1.0) + 12 * n) - 1] + newdist * L[((int)(k +
            1.0) + 12 * n) - 1];
        }

        loop_ub = (int)k;
        for (n = 0; n < loop_ub; n++) {
          S[((int)k + 12 * n) - 1] = S_data[n];
        }

        acond[(int)k - 1] = ahat[(int)k - 1] + S[((int)k + 12 * ((int)k - 1)) -
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
          for (n = 0; n < 12; n++) {
            afixed->data[n + afixed->size[0] * ((int)count - 1)] = zcond[n];
          }

          sqnorm->data[(int)count - 1] = newdist;

          // Store F(zcond)
        } else {
          for (n = 0; n < 12; n++) {
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
      if (k == 12.0) {
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
  b_afixed->size[0] = 12;
  b_afixed->size[1] = iidx->size[1];
  emxEnsureCapacity((emxArray__common *)b_afixed, n, sizeof(double));
  loop_ub = iidx->size[1];
  for (n = 0; n < loop_ub; n++) {
    for (ix = 0; ix < 12; ix++) {
      b_afixed->data[ix + b_afixed->size[0] * n] = afixed->data[ix +
        afixed->size[0] * (iidx->data[iidx->size[0] * n] - 1)];
    }
  }

  emxFree_int32_T(&iidx);
  n = afixed->size[0] * afixed->size[1];
  afixed->size[0] = 12;
  afixed->size[1] = b_afixed->size[1];
  emxEnsureCapacity((emxArray__common *)afixed, n, sizeof(double));
  loop_ub = b_afixed->size[1];
  for (n = 0; n < loop_ub; n++) {
    for (ix = 0; ix < 12; ix++) {
      afixed->data[ix + afixed->size[0] * n] = b_afixed->data[ix +
        b_afixed->size[0] * n];
    }
  }

  emxFree_real_T(&b_afixed);
}
// --------------------------------------------------------------------------
//                            End of ssearch
// --------------------------------------------------------------------------
