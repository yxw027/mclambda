//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: lusolve.cpp
//
// MATLAB Coder version            : 3.3
// C/C++ source code generated on  : 08-Aug-2019 14:38:13
//

// Include Files
#include "rt_nonfinite.h"
#include "..\LAMBDA.h"
#include "lusolve.h"
#include "colon.cpp"

// Function Definitions

//
// Arguments    : const double A_data[]
//                const int A_size[2]
//                double B_data[]
//                int B_size[2]
// Return Type  : void
//
void lusolve(const double A_data[], const int A_size[2], double B_data[], int
             B_size[2])
{
  int n;
  int jBcol;
  int k;
  int ipiv_data[12];
  int ipiv_size[2];
  double b_A_data[144];
  int u1;
  int j;
  int nb;
  int mmj;
  int c;
  int ix;
  int jAcol;
  double temp;
  int kBcol;
  int i;
  double s;
  n = A_size[1];
  jBcol = A_size[0] * A_size[1];
  for (k = 0; k < jBcol; k++) {
    b_A_data[k] = A_data[k];
  }

  eml_signed_integer_colon(A_size[1], ipiv_data, ipiv_size);
  jBcol = A_size[1] - 1;
  u1 = A_size[1];
  if (jBcol < u1) {
    u1 = jBcol;
  }

  for (j = 0; j + 1 <= u1; j++) {
    mmj = n - j;
    c = j * (n + 1);
    if (mmj < 1) {
      jBcol = -1;
    } else {
      jBcol = 0;
      if (mmj > 1) {
        ix = c;
        temp = std::abs(b_A_data[c]);
        for (k = 2; k <= mmj; k++) {
          ix++;
          s = std::abs(b_A_data[ix]);
          if (s > temp) {
            jBcol = k - 1;
            temp = s;
          }
        }
      }
    }

    if (b_A_data[c + jBcol] != 0.0) {
      if (jBcol != 0) {
        ipiv_data[j] = (j + jBcol) + 1;
        ix = j;
        jBcol += j;
        for (k = 1; k <= n; k++) {
          temp = b_A_data[ix];
          b_A_data[ix] = b_A_data[jBcol];
          b_A_data[jBcol] = temp;
          ix += n;
          jBcol += n;
        }
      }

      k = c + mmj;
      for (i = c + 1; i + 1 <= k; i++) {
        b_A_data[i] /= b_A_data[c];
      }
    }

    jBcol = n - j;
    jAcol = (c + n) + 1;
    kBcol = c + n;
    for (nb = 1; nb < jBcol; nb++) {
      temp = b_A_data[kBcol];
      if (b_A_data[kBcol] != 0.0) {
        ix = c + 1;
        k = mmj + jAcol;
        for (i = jAcol; i + 1 < k; i++) {
          b_A_data[i] += b_A_data[ix] * -temp;
          ix++;
        }
      }

      kBcol += n;
      jAcol += n;
    }
  }

  nb = B_size[0];
  if (B_size[0] != 0) {
    for (j = 0; j + 1 <= n; j++) {
      jBcol = nb * j;
      jAcol = n * j;
      for (k = 0; k + 1 <= j; k++) {
        kBcol = nb * k;
        if (b_A_data[k + jAcol] != 0.0) {
          for (i = 0; i + 1 <= nb; i++) {
            B_data[i + jBcol] -= b_A_data[k + jAcol] * B_data[i + kBcol];
          }
        }
      }

      temp = 1.0 / b_A_data[j + jAcol];
      for (i = 0; i + 1 <= nb; i++) {
        B_data[i + jBcol] *= temp;
      }
    }
  }

  if (B_size[0] != 0) {
    for (j = A_size[1]; j > 0; j--) {
      jBcol = nb * (j - 1);
      jAcol = n * (j - 1);
      for (k = j; k + 1 <= n; k++) {
        kBcol = nb * k;
        if (b_A_data[k + jAcol] != 0.0) {
          for (i = 0; i + 1 <= nb; i++) {
            B_data[i + jBcol] -= b_A_data[k + jAcol] * B_data[i + kBcol];
          }
        }
      }
    }
  }

  for (jBcol = A_size[1] - 2; jBcol + 1 > 0; jBcol--) {
    if (ipiv_data[jBcol] != jBcol + 1) {
      jAcol = ipiv_data[jBcol] - 1;
      for (kBcol = 0; kBcol + 1 <= nb; kBcol++) {
        temp = B_data[kBcol + B_size[0] * jBcol];
        B_data[kBcol + B_size[0] * jBcol] = B_data[kBcol + B_size[0] * jAcol];
        B_data[kBcol + B_size[0] * jAcol] = temp;
      }
    }
  }
}

//
// File trailer for lusolve.cpp
//
// [EOF]
//
