// --------------------------------------------------------------------------
//                             DECORREL ROUTINE
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
//  This routine creates a decorrelated Q-matrix, by finding the
//  Z-matrix and performing the corresponding transformation.
//
//  INPUTS:
//
//    Qahat: Variance-covariance matrix of ambiguities (original)
//    ahat : Original ambiguities (optional)
//
//  OUTPUTS:
//
//    Qzhat: Variance-covariance matrix of decorrelated ambiguities
//    Z    : Z-transformation matrix
//    L    : L matrix (from LtDL-decomposition of Qzhat)
//    D    : D matrix (from LtDL-decomposition of Qzhat)
//    zhat : Transformed ambiguities (optional)
//    iZt  : inv(Z')-transformation matrix
//
// --------------------------------------------------------------------------

// Include Files
#include "decorrel.h"

using namespace std;

// --------------------------------------------------------------------------
//                          Function Definitions
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
//                               decorrel
// --------------------------------------------------------------------------
//
// Arguments    : const double Qahat -> Variance-covariance matrix of ambiguities
//                const double ahat  -> Vector of float ambiguities
//                double Qzhat       -> Output Qzhat
//                double Z           -> Output Z
//                double L           -> Output L
//                double D           -> Output D
//                double zhat        -> Output zhat
//                double iZt         -> Output iZt
//
// Return       : void
//
// --------------------------------------------------------------------------
void decorrel(const double Qahat[144], const double ahat[12], double Qzhat[144],
              double Z[144], double L[144], double D[12], double zhat[12],
              double iZt[144])
{
  // ============================= VARIABLES ================================
  creal_T unusedExpr[12];
  int i1;
  int sw;
  int ib;
  int i;
  double b_Qahat[144];
  double delta;
  int loop_ub;
  int br;
  double L_data[11];
  double b_iZt[12];
  int ar;
  double lambda;
  double eta;
  double a[4];
  double b_data[20];
  double C_data[20];
  int ic;
  int ia;
  double b_L_data[24];
  double c_iZt[24];

  // ============================ START PROGRAM ===============================
  // Tests on Inputs ahat and Qahat
  // - Is the Q-matrix symmetric?
  // - Is the Q-matrix positive-definite?
  eig(Qahat, unusedExpr);

  // Initialization
  memset(&iZt[0], 0, 144U * sizeof(double));
  for (i1 = 0; i1 < 12; i1++) {
    iZt[i1 + 12 * i1] = 1.0;
  }
  i1 = 10;
  sw = 1;

  // ========================== LDL Decomposition =============================
  for (ib = 0; ib < 144; ib++) {
    b_Qahat[ib] = Qahat[ib];
    L[ib] = 0.0;
  }

  for (i = 0; i < 12; i++) {
    D[11 - i] = b_Qahat[(12 * (11 - i) - i) + 11];
    delta = std::sqrt(b_Qahat[(12 * (11 - i) - i) + 11]);
    loop_ub = 11 - i;
    for (ib = 0; ib <= loop_ub; ib++) {
      L[(12 * ib - i) + 11] = b_Qahat[(12 * ib - i) + 11] / delta;
    }

    for (br = 0; br <= 10 - i; br++) {
      for (ib = 0; ib <= br; ib++) {
        L_data[ib] = b_Qahat[br + 12 * ib] - L[(12 * ib - i) + 11] * L[(12 * br
          - i) + 11];
      }

      loop_ub = br + 1;
      for (ib = 0; ib < loop_ub; ib++) {
        b_Qahat[br + 12 * ib] = L_data[ib];
      }
    }

    loop_ub = 12 - i;
    for (ib = 0; ib < loop_ub; ib++) {
      b_iZt[ib] = L[(12 * ib - i) + 11] / L[(12 * (11 - i) - i) + 11];
    }

    loop_ub = 12 - i;
    for (ib = 0; ib < loop_ub; ib++) {
      L[(12 * ib - i) + 11] = b_iZt[ib];
    }
  }

  // The actual decorrelation procedure
  while (sw != 0) {
    i = 11;
    // Loop for column from n to 1
    sw = 0;
    while ((!(sw != 0)) && (i + 1 > 1)) {
      i--;
      // The i-th column
      if (i + 1 <= i1 + 1) {
        for (br = 0; br <= 10 - i; br++) {
          ar = (i + br) + 1;
          delta = rt_roundd_snf(L[ar + 12 * i]);
          if (delta != 0.0) {
            // If mu not equal to 0
            loop_ub = -ar;
            for (ib = 0; ib <= loop_ub + 11; ib++) {
              L_data[ib] = L[(ar + ib) + 12 * i] - delta * L[(ar + ib) + 12 * ar];
            }

            loop_ub = 12 - ar;
            for (ib = 0; ib < loop_ub; ib++) {
              L[(ar + ib) + 12 * i] = L_data[ib];
            }

            for (ib = 0; ib < 12; ib++) {
              b_iZt[ib] = iZt[ib + 12 * ar] + delta * iZt[ib + 12 * i];
            }

            memcpy(&iZt[ar * 12], &b_iZt[0], 12U * sizeof(double));

            // iZt is inv(Zt) matrix
          }
        }
      }

      delta = D[i] + L[(i + 12 * i) + 1] * L[(i + 12 * i) + 1] * D[i + 1];
      if (delta < D[i + 1]) {
        lambda = D[i + 1] * L[(i + 12 * i) + 1] / delta;
        eta = D[i] / delta;
        D[i] = eta * D[i + 1];
        D[i + 1] = delta;
        if (1 > i) {
          loop_ub = 0;
        } else {
          loop_ub = i;
        }

        a[0] = -L[(i + 12 * i) + 1];
        a[2] = 1.0;
        a[1] = eta;
        a[3] = lambda;
        for (ib = 0; ib < loop_ub; ib++) {
          for (ar = 0; ar < 2; ar++) {
            b_data[ar + (ib << 1)] = L[(ar + i) + 12 * ib];
          }
        }

        i1 = (signed char)loop_ub;
        for (ib = 0; ib < i1; ib++) {
          for (ar = 0; ar < 2; ar++) {
            C_data[ar + (ib << 1)] = 0.0;
          }
        }

        if (loop_ub != 0) {
          i1 = (loop_ub - 1) << 1;
          for (sw = 0; sw <= i1; sw += 2) {
            for (ic = sw + 1; ic <= sw + 2; ic++) {
              C_data[ic - 1] = 0.0;
            }
          }

          br = 0;
          for (sw = 0; sw <= i1; sw += 2) {
            ar = -1;
            for (ib = br; ib + 1 <= br + 2; ib++) {
              if (b_data[ib] != 0.0) {
                ia = ar;
                for (ic = sw; ic + 1 <= sw + 2; ic++) {
                  ia++;
                  C_data[ic] += b_data[ib] * a[ia];
                }
              }

              ar += 2;
            }

            br += 2;
          }
        }

        loop_ub = (signed char)loop_ub;
        for (ib = 0; ib < loop_ub; ib++) {
          for (ar = 0; ar < 2; ar++) {
            L[(ar + i) + 12 * ib] = C_data[ar + (ib << 1)];
          }
        }

        L[(i + 12 * i) + 1] = lambda;

        // Swap rows i and i+1
        if (i + 3 > 12) {
          ib = -2;
          ar = 1;
          i1 = -2;
        } else {
          ib = i;
          ar = 13;
          i1 = i;
        }

        sw = (ar - ib) - 3;
        loop_ub = ar - ib;
        for (ar = 0; ar < 2; ar++) {
          for (br = 0; br <= loop_ub - 4; br++) {
            b_L_data[br + sw * ar] = L[((ib + br) + 12 * ((i - ar) + 1)) + 2];
          }
        }

        for (ib = 0; ib < 2; ib++) {
          for (ar = 0; ar < sw; ar++) {
            L[((i1 + ar) + 12 * (ib + i)) + 2] = b_L_data[ar + sw * ib];
          }

          memcpy(&c_iZt[ib * 12], &iZt[(ib * -12 + i * 12) + 12], 12U * sizeof
                 (double));
        }

        for (ib = 0; ib < 2; ib++) {
          memcpy(&iZt[ib * 12 + i * 12], &c_iZt[ib * 12], 12U * sizeof(double));
        }

        i1 = i;
        sw = 1;
      }
    }
  }

  // Return the transformed Q-matrix and the transformation-matrix
  // Return the decorrelated ambiguities, if they were supplied
  for (ib = 0; ib < 12; ib++) {
    for (ar = 0; ar < 12; ar++) {
      b_Qahat[ar + 12 * ib] = iZt[ib + 12 * ar];
    }
  }

  inv(b_Qahat, Z);
  for (i1 = 0; i1 < 144; i1++) {
    Z[i1] = rt_roundd_snf(Z[i1]);
  }

  for (ib = 0; ib < 12; ib++) {
    for (ar = 0; ar < 12; ar++) {
      b_Qahat[ib + 12 * ar] = 0.0;
      for (i1 = 0; i1 < 12; i1++) {
        b_Qahat[ib + 12 * ar] += Z[i1 + 12 * ib] * Qahat[i1 + 12 * ar];
      }
    }

    zhat[ib] = 0.0;
    for (ar = 0; ar < 12; ar++) {
      Qzhat[ib + 12 * ar] = 0.0;
      for (i1 = 0; i1 < 12; i1++) {
        Qzhat[ib + 12 * ar] += b_Qahat[ib + 12 * i1] * Z[i1 + 12 * ar];
      }

      zhat[ib] += Z[ar + 12 * ib] * ahat[ar];
    }
  }
}
// --------------------------------------------------------------------------
//                              End of decorrel
// --------------------------------------------------------------------------
