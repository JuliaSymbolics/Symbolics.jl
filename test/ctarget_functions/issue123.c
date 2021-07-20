#include <math.h>
void diffeqf(double* du, const double* RHS1) {
  du[0] = u * (M[0] * qd[0] + M[0] * qd[2] + M[0] * qd[3] + M[6] * qd[1] + M[24] * qd[4] + M[30] * qd[5]) * qd[0];
}
