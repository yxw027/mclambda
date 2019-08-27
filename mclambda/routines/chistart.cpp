// --------------------------------------------------------------------------
//                            CHISTART ROUTINE
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
//  This routine computes or approximates the initial size of the search
//  ellipsoid. If the requested number of candidates is not more than the
//  dimension + 1, this is done by computing the squared distances of partially
//  conditionally rounded float vectors to the float vector in the metric of the
//  covariance matrix. Otherwise an approximation is used.
//
//  INPUTS:
//
//     L,D   : LtDL-decomposition of the variance-covariance matrix of
//             the float ambiguities (preferably decorrelated)
//     ahat  : float ambiguites (preferably decorrelated)
//     ncands: Requested number of candidates
//     factor: Multiplication factor for the volume of the resulting
//             search ellipsoid
// 
//  OUTPUTS:
//
//     Chi2  : Size of the search ellipsoid
//
// --------------------------------------------------------------------------

// Include Files
#include "..\math_functions\rt_nonfinite.h"
#include "..\LAMBDA.h"
#include "chistart.h"
#include "..\LAMBDA_emxutil.h"
#include "..\math_functions\sort1.h"
#include "..\math_functions\xzgetrf.h"
#include "..\math_functions\xtrsm.h"
#include "..\LAMBDA_rtwutil.h"

using namespace std;

// --------------------------------------------------------------------------
//                          Function Definitions
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
//                               chistart
// --------------------------------------------------------------------------
//
// Arguments    : const double D    -> Matrix D of LtDL-decomposition
//                const double L    -> Matrix L of LtDL-decomposition
//                const double ahat -> Vector of float ambiguities
//                double ncands     -> Requested number of candidates
//
// Return       : double Chi2       -> Size of the search ellipsoid
//
// --------------------------------------------------------------------------
double chistart(const double D[12], const double L[144], const double ahat[12],
                double ncands)
{
  // ============================= VARIABLES ================================
  double Chi2;
  emxArray_real_T *Chi;
  int jBcol;
  double iQ[144];
  double afloat[12];
  int j;
  double A[144];
  int ipiv[12];
  int jp;
  int kBcol;
  int k;
  double dw;
  int i;
  double afixed[12];
  double x;
  double b_ahat[12];

  // ============================ START PROGRAM ===============================
  // Computation depends on the number of candidates to be computed
  int n = sizeof(ahat);
  if (ncands <= n + 1) {
    emxInit_real_T(&Chi, 2);
    // Computation based on the bootstrapping estimator
    jBcol = Chi->size[0] * Chi->size[1];
    Chi->size[0] = 1;
    Chi->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)Chi, jBcol, sizeof(double));
    Chi->data[0] = 0.0;
    for (jBcol = 0; jBcol < 12; jBcol++) {
      afloat[jBcol] = 1.0 / D[jBcol];
    }

    memset(&iQ[0], 0, 144U * sizeof(double));
    for (j = 0; j < 12; j++) {
      iQ[j + 12 * j] = afloat[j];
    }

    memcpy(&A[0], &L[0], 144U * sizeof(double));
    xzgetrf(A, ipiv, &jp);
    for (kBcol = 0; kBcol < 11; kBcol++) {
      if (ipiv[kBcol] != kBcol + 1) {
        jp = ipiv[kBcol] - 1;
        for (jBcol = 0; jBcol < 12; jBcol++) {
          dw = iQ[kBcol + 12 * jBcol];
          iQ[kBcol + 12 * jBcol] = iQ[jp + 12 * jBcol];
          iQ[jp + 12 * jBcol] = dw;
        }
      }
    }

    for (j = 0; j < 12; j++) {
      jBcol = 12 * j;
      for (k = 0; k < 12; k++) {
        jp = 12 * k;
        if (iQ[k + jBcol] != 0.0) {
          for (i = k + 1; i + 1 < 13; i++) {
            iQ[i + jBcol] -= iQ[k + jBcol] * A[i + jp];
          }
        }
      }
    }

    xtrsm(A, iQ);
    for (jBcol = 0; jBcol < 12; jBcol++) {
      for (kBcol = 0; kBcol < 12; kBcol++) {
        A[kBcol + 12 * jBcol] = L[jBcol + 12 * kBcol];
      }
    }

    xzgetrf(A, ipiv, &jp);
    for (j = 0; j < 12; j++) {
      jBcol = 12 * j;
      jp = 12 * j;
      for (k = 1; k <= j; k++) {
        kBcol = 12 * (k - 1);
        if (A[(k + jp) - 1] != 0.0) {
          for (i = 0; i < 12; i++) {
            iQ[i + jBcol] -= A[(k + jp) - 1] * iQ[i + kBcol];
          }
        }
      }

      dw = 1.0 / A[j + jp];
      for (i = 0; i < 12; i++) {
        iQ[i + jBcol] *= dw;
      }
    }

    for (j = 11; j >= 0; j += -1) {
      jBcol = 12 * j;
      jp = 12 * j - 1;
      for (k = j + 2; k < 13; k++) {
        kBcol = 12 * (k - 1);
        if (A[k + jp] != 0.0) {
          for (i = 0; i < 12; i++) {
            iQ[i + jBcol] -= A[k + jp] * iQ[i + kBcol];
          }
        }
      }
    }

    for (jBcol = 10; jBcol >= 0; jBcol += -1) {
      if (ipiv[jBcol] != jBcol + 1) {
        jp = ipiv[jBcol] - 1;
        for (kBcol = 0; kBcol < 12; kBcol++) {
          dw = iQ[kBcol + 12 * jBcol];
          iQ[kBcol + 12 * jBcol] = iQ[kBcol + 12 * jp];
          iQ[kBcol + 12 * jp] = dw;
        }
      }
    }

    for (k = 0; k < 13; k++) {
      memcpy(&afloat[0], &ahat[0], 12U * sizeof(double));
      memcpy(&afixed[0], &ahat[0], 12U * sizeof(double));
      for (i = 0; i < 12; i++) {
        dw = 0.0;
        jBcol = (int)(((12.0 + -(double)i) + -13.0) / -1.0);
        for (j = 0; j < jBcol; j++) {
          dw += L[(12 * (11 - i) - j) + 11] * (afloat[11 - j] - afixed[11 - j]);
        }

        afloat[11 - i] -= dw;
        if (12 - i != 12 - k) {
          afixed[11 - i] = rt_roundd_snf(afloat[11 - i]);
        } else {
          dw = rt_roundd_snf(afloat[11 - i]);
          x = afloat[11 - i] - dw;
          if (x < 0.0) {
            x = -1.0;
          } else if (x > 0.0) {
            x = 1.0;
          } else {
            if (x == 0.0) {
              x = 0.0;
            }
          }

          afixed[11 - i] = dw + x;
        }
      }

      jp = Chi->size[1];
      jBcol = Chi->size[0] * Chi->size[1];
      Chi->size[1] = jp + 1;
      emxEnsureCapacity((emxArray__common *)Chi, jBcol, sizeof(double));
      for (jBcol = 0; jBcol < 12; jBcol++) {
        afloat[jBcol] = ahat[jBcol] - afixed[jBcol];
      }

      dw = 0.0;
      for (jBcol = 0; jBcol < 12; jBcol++) {
        b_ahat[jBcol] = 0.0;
        for (kBcol = 0; kBcol < 12; kBcol++) {
          b_ahat[jBcol] += afloat[kBcol] * iQ[kBcol + 12 * jBcol];
        }

        dw += b_ahat[jBcol] * (ahat[jBcol] - afixed[jBcol]);
      }

      Chi->data[jp] = dw;
    }

    // Sort the results, and return the appropriate number
    // Add an "eps", to make sure there is no boundary problem
    b_sort(Chi);
    Chi2 = Chi->data[(int)ncands] + 1.0E-6;
    emxFree_real_T(&Chi);
  } else {
    // An approximation for the squared norm is computed
    float Vn = (2/n) * (pow(M_PI, (n/2)))/(tgamma(n/2));
    float prodD = 1;
    for (int i=0; i<sizeof(D); i++){
        prodD *=  D[i];
    }
    Chi2  = ncands * pow((ncands/(sqrt(prodD * Vn))), (2/n));

  }

  return Chi2;
}
// --------------------------------------------------------------------------
//                              End of chistart
// --------------------------------------------------------------------------
