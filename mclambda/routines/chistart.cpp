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
//     n     : Number of float ambiguities
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
#include "..\mclambda.h"
#include "chistart.h"
#include "..\mclambda_emxutil.h"
#include "..\math_functions\sort1.h"
#include "..\math_functions\xzgetrf.h"
#include "..\math_functions\xtrsm.h"
#include "..\mclambda_rtwutil.h"

using namespace std;

// --------------------------------------------------------------------------
//                          Function Definitions
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
//                               chistart
// --------------------------------------------------------------------------
//
// Arguments    : int n             -> Number of float ambiguities
//                const double D    -> Matrix D of LtDL-decomposition
//                const double L    -> Matrix L of LtDL-decomposition
//                const double ahat -> Vector of float ambiguities
//                double ncands     -> Requested number of candidates
//
// Return       : double Chi2       -> Size of the search ellipsoid
//
// --------------------------------------------------------------------------
double chistart(int n, const double D[], const double L[], const double ahat[],
                double ncands)
{
  // ============================= VARIABLES ================================
  double Chi2;
  emxArray_real_T *Chi;
  int jBcol;
  double iQ[n*n];
  double afloat[n];
  int j;
  double A[n*n];
  int ipiv[n];
  int jp;
  int kBcol;
  int k;
  double dw;
  int i;
  double afixed[n];
  double x;
  double b_ahat[n];

  // ============================ START PROGRAM ===============================
  // Computation depends on the number of candidates to be computed
  if (ncands <= n + 1) {
    emxInit_real_T(&Chi, 2);
    // Computation based on the bootstrapping estimator
    jBcol = Chi->size[0] * Chi->size[1];
    Chi->size[0] = 1;
    Chi->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)Chi, jBcol, sizeof(double));
    Chi->data[0] = 0.0;
    for (jBcol = 0; jBcol < n; jBcol++) {
      afloat[jBcol] = 1.0 / D[jBcol];
    }

    memset(&iQ[0], 0, n*n * sizeof(double));
    for (j = 0; j < n; j++) {
      iQ[j + n * j] = afloat[j];
    }

    memcpy(&A[0], &L[0], n*n * sizeof(double));
    xzgetrf(n, A, ipiv, &jp);
    for (kBcol = 0; kBcol < n-1; kBcol++) {
      if (ipiv[kBcol] != kBcol + 1) {
        jp = ipiv[kBcol] - 1;
        for (jBcol = 0; jBcol < n; jBcol++) {
          dw = iQ[kBcol + n * jBcol];
          iQ[kBcol + n * jBcol] = iQ[jp + n * jBcol];
          iQ[jp + n * jBcol] = dw;
        }
      }
    }

    for (j = 0; j < n; j++) {
      jBcol = n * j;
      for (k = 0; k < n; k++) {
        jp = n * k;
        if (iQ[k + jBcol] != 0.0) {
          for (i = k + 1; i + 1 < n+1; i++) {
            iQ[i + jBcol] -= iQ[k + jBcol] * A[i + jp];
          }
        }
      }
    }

    xtrsm(n, A, iQ);
    for (jBcol = 0; jBcol < n; jBcol++) {
      for (kBcol = 0; kBcol < n; kBcol++) {
        A[kBcol + n * jBcol] = L[jBcol + n * kBcol];
      }
    }

    xzgetrf(n, A, ipiv, &jp);
    for (j = 0; j < n; j++) {
      jBcol = n * j;
      jp = n * j;
      for (k = 1; k <= j; k++) {
        kBcol = n * (k - 1);
        if (A[(k + jp) - 1] != 0.0) {
          for (i = 0; i < n; i++) {
            iQ[i + jBcol] -= A[(k + jp) - 1] * iQ[i + kBcol];
          }
        }
      }

      dw = 1.0 / A[j + jp];
      for (i = 0; i < n; i++) {
        iQ[i + jBcol] *= dw;
      }
    }

    for (j = n-1; j >= 0; j += -1) {
      jBcol = n * j;
      jp = n * j - 1;
      for (k = j + 2; k < n+1; k++) {
        kBcol = n * (k - 1);
        if (A[k + jp] != 0.0) {
          for (i = 0; i < n; i++) {
            iQ[i + jBcol] -= A[k + jp] * iQ[i + kBcol];
          }
        }
      }
    }

    for (jBcol = 10; jBcol >= 0; jBcol += -1) {
      if (ipiv[jBcol] != jBcol + 1) {
        jp = ipiv[jBcol] - 1;
        for (kBcol = 0; kBcol < n; kBcol++) {
          dw = iQ[kBcol + n * jBcol];
          iQ[kBcol + n * jBcol] = iQ[kBcol + n * jp];
          iQ[kBcol + n * jp] = dw;
        }
      }
    }

    for (k = 0; k < n+1; k++) {
      memcpy(&afloat[0], &ahat[0], n * sizeof(double));
      memcpy(&afixed[0], &ahat[0], n * sizeof(double));
      for (i = 0; i < n; i++) {
        dw = 0.0;
        jBcol = (int)((((double)(n) + -(double)i) + -(n+1.0)) / -1.0);
        for (j = 0; j < jBcol; j++) {
          dw += L[(n * (n-1 - i) - j) + n-1] * (afloat[n-1 - j] - afixed[n-1 - j]);
        }

        afloat[n-1 - i] -= dw;
        if (n - i != n - k) {
          afixed[n-1 - i] = rt_roundd_snf(afloat[n-1 - i]);
        } else {
          dw = rt_roundd_snf(afloat[n-1 - i]);
          x = afloat[n-1 - i] - dw;
          if (x < 0.0) {
            x = -1.0;
          } else if (x > 0.0) {
            x = 1.0;
          } else {
            if (x == 0.0) {
              x = 0.0;
            }
          }

          afixed[n-1 - i] = dw + x;
        }
      }

      jp = Chi->size[1];
      jBcol = Chi->size[0] * Chi->size[1];
      Chi->size[1] = jp + 1;
      emxEnsureCapacity((emxArray__common *)Chi, jBcol, sizeof(double));
      for (jBcol = 0; jBcol < n; jBcol++) {
        afloat[jBcol] = ahat[jBcol] - afixed[jBcol];
      }

      dw = 0.0;
      for (jBcol = 0; jBcol < n; jBcol++) {
        b_ahat[jBcol] = 0.0;
        for (kBcol = 0; kBcol < n; kBcol++) {
          b_ahat[jBcol] += afloat[kBcol] * iQ[kBcol + n * jBcol];
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
    for (int i=0; i<n; i++){
        prodD *=  D[i];
    }
    Chi2  = ncands * pow((ncands/(sqrt(prodD * Vn))), (2/n));

  }

  return Chi2;
}
// --------------------------------------------------------------------------
//                              End of chistart
// --------------------------------------------------------------------------
