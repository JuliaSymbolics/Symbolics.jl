#include <math.h>
void diffeqf(double* du, const double* RHS1, const double* RHS2, const double RHS3) {
  du[0] = RHS3 + RHS1[0] + RHS2[0];
}
