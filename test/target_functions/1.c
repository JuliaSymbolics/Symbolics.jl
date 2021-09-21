#include <math.h>
void diffeqf(double* internal_var___du, const double* internal_var___u, const double* internal_var___p, const double t) {
  internal_var___du[0] = internal_var___p[0] * internal_var___u[0] + -1 * internal_var___u[0] * internal_var___u[1];
  internal_var___du[1] = -3 * internal_var___u[1] + internal_var___u[0] * internal_var___u[1];
}
