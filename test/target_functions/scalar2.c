#include <math.h>
void diffeqf(double* du, const double* RHS1, const double* RHS2, const double RHS3) {
  du[0] = 2 * RHS3 * 1 + pow(RHS1[0], 2) + 1 / RHS1[1] + pow(sin(RHS2[0]), 3.5);
}
