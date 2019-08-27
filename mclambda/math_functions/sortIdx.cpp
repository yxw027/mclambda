//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: sortIdx.cpp
//
// MATLAB Coder version            : 3.3
// C/C++ source code generated on  : 08-Aug-2019 14:38:13
//

// Include Files
#include "rt_nonfinite.h"
#include "..\LAMBDA.h"
#include "sortIdx.h"
#include "..\LAMBDA_emxutil.h"

// Function Declarations
static void merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int np,
                  int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
  int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork);

// Function Definitions

//
// Arguments    : emxArray_int32_T *idx
//                emxArray_real_T *x
//                int offset
//                int np
//                int nq
//                emxArray_int32_T *iwork
//                emxArray_real_T *xwork
// Return Type  : void
//
static void merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int np,
                  int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
  int n;
  int qend;
  int p;
  int iout;
  int exitg1;
  if (nq != 0) {
    n = np + nq;
    for (qend = 0; qend + 1 <= n; qend++) {
      iwork->data[qend] = idx->data[offset + qend];
      xwork->data[qend] = x->data[offset + qend];
    }

    p = 0;
    n = np;
    qend = np + nq;
    iout = offset - 1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork->data[p] <= xwork->data[n]) {
        idx->data[iout] = iwork->data[p];
        x->data[iout] = xwork->data[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx->data[iout] = iwork->data[n];
        x->data[iout] = xwork->data[n];
        if (n + 1 < qend) {
          n++;
        } else {
          n = (iout - p) + 1;
          while (p + 1 <= np) {
            idx->data[n + p] = iwork->data[p];
            x->data[n + p] = xwork->data[p];
            p++;
          }

          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

//
// Arguments    : emxArray_int32_T *idx
//                emxArray_real_T *x
//                int offset
//                int n
//                int preSortLevel
//                emxArray_int32_T *iwork
//                emxArray_real_T *xwork
// Return Type  : void
//
static void merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
  int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
  int nPairs;
  int bLen;
  int tailOffset;
  int nTail;
  nPairs = n >> preSortLevel;
  bLen = 1 << preSortLevel;
  while (nPairs > 1) {
    if ((nPairs & 1) != 0) {
      nPairs--;
      tailOffset = bLen * nPairs;
      nTail = n - tailOffset;
      if (nTail > bLen) {
        merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
      }
    }

    tailOffset = bLen << 1;
    nPairs >>= 1;
    for (nTail = 1; nTail <= nPairs; nTail++) {
      merge(idx, x, offset + (nTail - 1) * tailOffset, bLen, bLen, iwork, xwork);
    }

    bLen = tailOffset;
  }

  if (n > bLen) {
    merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
  }
}

//
// Arguments    : emxArray_real_T *x
//                emxArray_int32_T *idx
// Return Type  : void
//
void sortIdx(emxArray_real_T *x, emxArray_int32_T *idx)
{
  int ib;
  int m;
  int n;
  int b_n;
  int i;
  emxArray_int32_T *iwork;
  double x4[4];
  int idx4[4];
  emxArray_real_T *xwork;
  int nNaNs;
  int k;
  int wOffset;
  signed char perm[4];
  int nNonNaN;
  int i4;
  int nBlocks;
  int nPairs;
  int b_iwork[256];
  double b_xwork[256];
  int exitg1;
  ib = x->size[0];
  m = idx->size[0];
  idx->size[0] = ib;
  emxEnsureCapacity((emxArray__common *)idx, m, sizeof(int));
  for (m = 0; m < ib; m++) {
    idx->data[m] = 0;
  }

  n = x->size[0];
  b_n = x->size[0];
  for (i = 0; i < 4; i++) {
    x4[i] = 0.0;
    idx4[i] = 0;
  }

  emxInit_int32_T1(&iwork, 1);
  m = iwork->size[0];
  iwork->size[0] = ib;
  emxEnsureCapacity((emxArray__common *)iwork, m, sizeof(int));
  i = iwork->size[0];
  m = iwork->size[0];
  iwork->size[0] = i;
  emxEnsureCapacity((emxArray__common *)iwork, m, sizeof(int));
  for (m = 0; m < i; m++) {
    iwork->data[m] = 0;
  }

  emxInit_real_T1(&xwork, 1);
  i = x->size[0];
  m = xwork->size[0];
  xwork->size[0] = i;
  emxEnsureCapacity((emxArray__common *)xwork, m, sizeof(double));
  i = xwork->size[0];
  m = xwork->size[0];
  xwork->size[0] = i;
  emxEnsureCapacity((emxArray__common *)xwork, m, sizeof(double));
  for (m = 0; m < i; m++) {
    xwork->data[m] = 0.0;
  }

  nNaNs = 1;
  ib = 0;
  for (k = 0; k + 1 <= b_n; k++) {
    if (rtIsNaN(x->data[k])) {
      idx->data[b_n - nNaNs] = k + 1;
      xwork->data[b_n - nNaNs] = x->data[k];
      nNaNs++;
    } else {
      ib++;
      idx4[ib - 1] = k + 1;
      x4[ib - 1] = x->data[k];
      if (ib == 4) {
        i = k - nNaNs;
        if (x4[0] <= x4[1]) {
          ib = 1;
          m = 2;
        } else {
          ib = 2;
          m = 1;
        }

        if (x4[2] <= x4[3]) {
          wOffset = 3;
          i4 = 4;
        } else {
          wOffset = 4;
          i4 = 3;
        }

        if (x4[ib - 1] <= x4[wOffset - 1]) {
          if (x4[m - 1] <= x4[wOffset - 1]) {
            perm[0] = (signed char)ib;
            perm[1] = (signed char)m;
            perm[2] = (signed char)wOffset;
            perm[3] = (signed char)i4;
          } else if (x4[m - 1] <= x4[i4 - 1]) {
            perm[0] = (signed char)ib;
            perm[1] = (signed char)wOffset;
            perm[2] = (signed char)m;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)ib;
            perm[1] = (signed char)wOffset;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)m;
          }
        } else if (x4[ib - 1] <= x4[i4 - 1]) {
          if (x4[m - 1] <= x4[i4 - 1]) {
            perm[0] = (signed char)wOffset;
            perm[1] = (signed char)ib;
            perm[2] = (signed char)m;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)wOffset;
            perm[1] = (signed char)ib;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)m;
          }
        } else {
          perm[0] = (signed char)wOffset;
          perm[1] = (signed char)i4;
          perm[2] = (signed char)ib;
          perm[3] = (signed char)m;
        }

        idx->data[i - 2] = idx4[perm[0] - 1];
        idx->data[i - 1] = idx4[perm[1] - 1];
        idx->data[i] = idx4[perm[2] - 1];
        idx->data[i + 1] = idx4[perm[3] - 1];
        x->data[i - 2] = x4[perm[0] - 1];
        x->data[i - 1] = x4[perm[1] - 1];
        x->data[i] = x4[perm[2] - 1];
        x->data[i + 1] = x4[perm[3] - 1];
        ib = 0;
      }
    }
  }

  wOffset = b_n - nNaNs;
  if (ib > 0) {
    for (i = 0; i < 4; i++) {
      perm[i] = 0;
    }

    if (ib == 1) {
      perm[0] = 1;
    } else if (ib == 2) {
      if (x4[0] <= x4[1]) {
        perm[0] = 1;
        perm[1] = 2;
      } else {
        perm[0] = 2;
        perm[1] = 1;
      }
    } else if (x4[0] <= x4[1]) {
      if (x4[1] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 2;
        perm[2] = 3;
      } else if (x4[0] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 3;
        perm[2] = 2;
      } else {
        perm[0] = 3;
        perm[1] = 1;
        perm[2] = 2;
      }
    } else if (x4[0] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 1;
      perm[2] = 3;
    } else if (x4[1] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 3;
      perm[2] = 1;
    } else {
      perm[0] = 3;
      perm[1] = 2;
      perm[2] = 1;
    }

    for (k = 1; k <= ib; k++) {
      idx->data[(wOffset - ib) + k] = idx4[perm[k - 1] - 1];
      x->data[(wOffset - ib) + k] = x4[perm[k - 1] - 1];
    }
  }

  m = (nNaNs - 1) >> 1;
  for (k = 1; k <= m; k++) {
    i = idx->data[wOffset + k];
    idx->data[wOffset + k] = idx->data[b_n - k];
    idx->data[b_n - k] = i;
    x->data[wOffset + k] = xwork->data[b_n - k];
    x->data[b_n - k] = xwork->data[wOffset + k];
  }

  if (((nNaNs - 1) & 1) != 0) {
    x->data[(wOffset + m) + 1] = xwork->data[(wOffset + m) + 1];
  }

  nNonNaN = (n - nNaNs) + 1;
  i = 2;
  if (nNonNaN > 1) {
    if (n >= 256) {
      nBlocks = nNonNaN >> 8;
      if (nBlocks > 0) {
        for (wOffset = 1; wOffset <= nBlocks; wOffset++) {
          i4 = ((wOffset - 1) << 8) - 1;
          for (nNaNs = 0; nNaNs < 6; nNaNs++) {
            b_n = 1 << (nNaNs + 2);
            n = b_n << 1;
            nPairs = 256 >> (nNaNs + 3);
            for (k = 1; k <= nPairs; k++) {
              ib = i4 + (k - 1) * n;
              for (i = 1; i <= n; i++) {
                b_iwork[i - 1] = idx->data[ib + i];
                b_xwork[i - 1] = x->data[ib + i];
              }

              m = 0;
              i = b_n;
              do {
                exitg1 = 0;
                ib++;
                if (b_xwork[m] <= b_xwork[i]) {
                  idx->data[ib] = b_iwork[m];
                  x->data[ib] = b_xwork[m];
                  if (m + 1 < b_n) {
                    m++;
                  } else {
                    exitg1 = 1;
                  }
                } else {
                  idx->data[ib] = b_iwork[i];
                  x->data[ib] = b_xwork[i];
                  if (i + 1 < n) {
                    i++;
                  } else {
                    i = ib - m;
                    while (m + 1 <= b_n) {
                      idx->data[(i + m) + 1] = b_iwork[m];
                      x->data[(i + m) + 1] = b_xwork[m];
                      m++;
                    }

                    exitg1 = 1;
                  }
                }
              } while (exitg1 == 0);
            }
          }
        }

        i = nBlocks << 8;
        ib = nNonNaN - i;
        if (ib > 0) {
          merge_block(idx, x, i, ib, 2, iwork, xwork);
        }

        i = 8;
      }
    }

    merge_block(idx, x, 0, nNonNaN, i, iwork, xwork);
  }

  emxFree_real_T(&xwork);
  emxFree_int32_T(&iwork);
}

//
// File trailer for sortIdx.cpp
//
// [EOF]
//
