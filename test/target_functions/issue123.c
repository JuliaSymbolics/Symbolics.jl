#include <math.h>
void diffeqf(double* du, const double* RHS1) {
  du[0] = (M[0] * qd[0] + M[6] * qd[1]) * qd[0] * u;
}
