/* TomsFastMath, a fast ISO C bignum library.
 *
 * This project is meant to fill in where LibTomMath
 * falls short.  That is speed ;-)
 *
 * This project is public domain and free for all purposes.
 *
 * Tom St Denis, tomstdenis@gmail.com
 */

/* program emits a NxN comba multiplier */
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
   int N, x, y, z;
   N = atoi(argv[1]);

   /* print out preamble */
printf(
"#define TFM_DEFINES\n"
"#include \"fp_mul_comba.c\"\n"
"\n"
"#if defined(TFM_MUL%d) && FP_SIZE >= %d\n"
"void fp_mul_comba%d(fp_int *A, fp_int *B, fp_int *C)\n"
"{\n"
"   fp_digit c0, c1, c2, at[%d];\n", N, N+N, N, N+N);
if (N >= 32) {
printf(
"   int out_size;\n"
"\n"
"   out_size = A->used + B->used;\n");
} else {
printf(
"\n");
}
printf(
"   memcpy(at, A->dp, %d * sizeof(fp_digit));\n"
"   memcpy(at+%d, B->dp, %d * sizeof(fp_digit));\n"
"   COMBA_START;\n"
"\n"
"   COMBA_CLEAR;\n", N, N, N);

   /* now do the rows */
   for (x = 0; x < (N+N-1); x++) {
printf(
"   /* %d */\n", x);
if (x > 0) {
printf(
"   COMBA_FORWARD;\n");
}
      for (y = 0; y < N; y++) {
      for (z = 0; z < N; z++) {
          if ((y+z)==x) {
             printf("   MULADD(at[%d], at[%d]); ", y, z+N);
          }
      }
      }
printf(
"\n"
"   COMBA_STORE(C->dp[%d]);\n", x);
if (N >= 32 && N*2 != (x+2) &&(x+2) >= 40 && (x+2)%8 == 0) {
   printf(
"\n"
"   /* early out at %d digits, %d*32==%d, or two %d bit operands */\n"
"   if (out_size <= %d) { COMBA_STORE2(C->dp[%d]); C->used = %d; C->sign = A->sign ^ B->sign; fp_clamp(C); COMBA_FINI; return; }\n"
"\n"
   , x+2, x+2, (x+2)*32, (x+2)*16, x+2, x+1, x+2);
}
   }
printf(
"   COMBA_STORE2(C->dp[%d]);\n"
"   C->used = %d;\n"
"   C->sign = A->sign ^ B->sign;\n"
"   fp_clamp(C);\n"
"   COMBA_FINI;\n"
"}\n#endif\n\n\n"
"/* $Source$ */\n"
"/* $Revision$ */\n"
"/* $Date$ */\n"
, N+N-1, N+N);

  return 0;
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */
