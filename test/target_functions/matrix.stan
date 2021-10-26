vector diffeqf(real Any[],vector internal_var___u,vector internal_var___p) {
  vector[12] internal_var___du;
  internal_var___du[1] = internal_var___u[1];
  internal_var___du[2] = internal_var___u[2];
  internal_var___du[3] = internal_var___u[3];
  internal_var___du[4] = internal_var___u[4];
  internal_var___du[5] = internal_var___u[5];
  internal_var___du[6] = internal_var___u[6];
  internal_var___du[7] = internal_var___u[7];
  internal_var___du[8] = internal_var___u[8];
  internal_var___du[9] = internal_var___p[1];
  internal_var___du[10] = internal_var___p[2];
  internal_var___du[11] = internal_var___p[3];
  internal_var___du[12] = internal_var___p[4];
  return internal_var___du;
}
