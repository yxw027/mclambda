// Include Files
#include "rt_nonfinite.h"
#include "..\mclambda.h"
#include "sort1.h"
#include "..\mclambda_emxutil.h"
#include "sortIdx.cpp"

// --------------------------------------------------------------------------
void b_sort(emxArray_real_T *x)
{
  emxArray_real_T *vwork;
  int i12;
  int x_idx_0;
  int i13;
  emxArray_int32_T *b_vwork;
  emxInit_real_T1(&vwork, 1);
  i12 = x->size[1];
  x_idx_0 = x->size[1];
  i13 = vwork->size[0];
  vwork->size[0] = x_idx_0;
  emxEnsureCapacity((emxArray__common *)vwork, i13, sizeof(double));
  for (x_idx_0 = 0; x_idx_0 + 1 <= i12; x_idx_0++) {
    vwork->data[x_idx_0] = x->data[x_idx_0];
  }

  emxInit_int32_T1(&b_vwork, 1);
  sortIdx(vwork, b_vwork);
  x_idx_0 = 0;
  emxFree_int32_T(&b_vwork);
  while (x_idx_0 + 1 <= i12) {
    x->data[x_idx_0] = vwork->data[x_idx_0];
    x_idx_0++;
  }

  emxFree_real_T(&vwork);
}
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
void sort(emxArray_real_T *x, emxArray_int32_T *idx)
{
  emxArray_real_T *vwork;
  int i10;
  int x_idx_0;
  int i11;
  emxArray_int32_T *iidx;
  emxInit_real_T1(&vwork, 1);
  i10 = x->size[1];
  x_idx_0 = x->size[1];
  i11 = vwork->size[0];
  vwork->size[0] = x_idx_0;
  emxEnsureCapacity((emxArray__common *)vwork, i11, sizeof(double));
  i11 = idx->size[0] * idx->size[1];
  idx->size[0] = 1;
  idx->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)idx, i11, sizeof(int));
  for (x_idx_0 = 0; x_idx_0 + 1 <= i10; x_idx_0++) {
    vwork->data[x_idx_0] = x->data[x_idx_0];
  }

  emxInit_int32_T1(&iidx, 1);
  sortIdx(vwork, iidx);
  for (x_idx_0 = 0; x_idx_0 + 1 <= i10; x_idx_0++) {
    x->data[x_idx_0] = vwork->data[x_idx_0];
    idx->data[x_idx_0] = iidx->data[x_idx_0];
  }

  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
}
// --------------------------------------------------------------------------
