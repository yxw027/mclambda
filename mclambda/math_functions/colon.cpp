//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: colon.cpp
//
// MATLAB Coder version            : 3.3
// C/C++ source code generated on  : 08-Aug-2019 14:38:13
//

// Include Files
#include "rt_nonfinite.h"
#include "..\LAMBDA.h"
#include "colon.h"

// Function Definitions

//
// Arguments    : int b
//                int y_data[]
//                int y_size[2]
// Return Type  : void
//
void eml_signed_integer_colon(int b, int y_data[], int y_size[2])
{
  int yk;
  int k;
  y_size[0] = 1;
  y_size[1] = b;
  y_data[0] = 1;
  yk = 1;
  for (k = 2; k <= b; k++) {
    yk++;
    y_data[k - 1] = yk;
  }
}

//
// File trailer for colon.cpp
//
// [EOF]
//
