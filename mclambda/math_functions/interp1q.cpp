// Include Files
#include "rt_nonfinite.h"
#include "interp1q.h"

// --------------------------------------------------------------------------
double interp1q(const double x[31], const double y[31], double xi)
{
  double yi;
  int low_i;
  int low_ip1;
  int high_i;
  double r;
  int mid_i;
  yi = rtNaN;
  if ((xi >= x[0]) && (xi <= x[30])) {
    low_i = 1;
    low_ip1 = 2;
    high_i = 31;
    while (high_i > low_ip1) {
      mid_i = (low_i + high_i) >> 1;
      if (xi >= x[mid_i - 1]) {
        low_i = mid_i;
        low_ip1 = mid_i + 1;
      } else {
        high_i = mid_i;
      }
    }

    r = (xi - x[low_i - 1]) / (x[low_i] - x[low_i - 1]);
    if (y[low_i - 1] == y[low_i]) {
      yi = y[low_i - 1];
    } else {
      yi = (1.0 - r) * y[low_i - 1] + r * y[low_i];
    }
  }

  return yi;
}
// --------------------------------------------------------------------------
