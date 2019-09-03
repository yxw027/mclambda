// Include Files
#include "rt_nonfinite.h"
#include "..\mclambda.h"
#include "isfinite.h"

// --------------------------------------------------------------------------
boolean_T b_isfinite(double x)
{
  return (!rtIsInf(x)) && (!rtIsNaN(x));
}
// --------------------------------------------------------------------------
