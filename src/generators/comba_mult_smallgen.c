/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */

/* program emits a NxN comba multiplier for 1x1 to 16x16 */
#include <stdio.h>

int main(int argc, char **argv)
{
   int N, x, y, z;

   /* print out preamble */
printf(
"/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */\n"
"/* SPDX-License-Identifier: Unlicense */\n"
"#ifndef TFM_PRE_GEN_MPI_C\n"
"#define TFM_DEFINES\n"
"#include \"fp_mul_comba.c\"\n"
"#endif\n"
"\n"
"#if defined(TFM_SMALL_SET)\n"
"void fp_mul_comba_small(const fp_int *A, const fp_int *B, fp_int *C)\n"
"{\n"
"   fp_digit c0, c1, c2, at[32];\n"
"   switch (MAX(A->used, B->used)) {\n"
);

for (N = 1; N <= 16; N++) {

printf(
"\n"
"   case %d:\n"
"      memcpy(at, A->dp, %d * sizeof(fp_digit));\n"
"      memcpy(at+%d, B->dp, %d * sizeof(fp_digit));\n"
"      COMBA_START;\n"
"\n"
"      COMBA_CLEAR;\n", N, N, N, N);

   /* now do the rows */
   for (x = 0; x < (N+N-1); x++) {
printf(
"      /* %d */\n", x);
if (x > 0) {
printf(
"      COMBA_FORWARD;\n");
}
printf("   ");
      for (y = 0; y < N; y++) {
      for (z = 0; z < N; z++) {
          if ((y+z)==x) {
             printf("   MULADD(at[%d], at[%d]);", y, z+N);
          }
      }
      }
printf(
"\n"
"      COMBA_STORE(C->dp[%d]);\n", x);
   }
printf(
"      COMBA_STORE2(C->dp[%d]);\n"
"      C->used = %d;\n"
"      C->sign = A->sign ^ B->sign;\n"
"      fp_clamp(C);\n"
"      COMBA_FINI;\n"
"      break;\n", N+N-1, N+N);
}
printf("   }\n}\n\n#endif\n");

  return 0;
}
