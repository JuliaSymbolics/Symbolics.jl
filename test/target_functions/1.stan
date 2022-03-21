vector diffeqf(real t,vector internal_var___u,vector internal_var___p) {
  vector[2] internal_var___du;
  internal_var___du[1] = internal_var___p[1] * internal_var___u[1] + -1 * internal_var___u[1] * internal_var___u[2];
  internal_var___du[2] = -3 * internal_var___u[2] + internal_var___u[1] * internal_var___u[2];
  return internal_var___du;
}
