/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
#ifndef TFM_PRE_GEN_MPI_C
#define TFM_DEFINES
#include "fp_sqr_comba.c"
#endif

#if defined(TFM_SQR8) && FP_SIZE >= 16
void fp_sqr_comba8(fp_int *A, fp_int *B)
{
   fp_digit *a, b[16], c0, c1, c2;
   fp_digit sc0, sc1, sc2;
#ifdef TFM_ISO
   fp_word tt;
#endif

   a = A->dp;
   COMBA_START;

   /* clear carries */
   CLEAR_CARRY;

   /* output 0 */
   SQRADD(a[0],a[0]);
   COMBA_STORE(b[0]);

   /* output 1 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[1]);
   COMBA_STORE(b[1]);

   /* output 2 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[2]); SQRADD(a[1], a[1]);
   COMBA_STORE(b[2]);

   /* output 3 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[3]); SQRADD2(a[1], a[2]);
   COMBA_STORE(b[3]);

   /* output 4 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[4]); SQRADD2(a[1], a[3]); SQRADD(a[2], a[2]);
   COMBA_STORE(b[4]);

   /* output 5 */
   CARRY_FORWARD;
   SQRADDSC(a[0], a[5]); SQRADDAC(a[1], a[4]); SQRADDAC(a[2], a[3]); SQRADDDB;
   COMBA_STORE(b[5]);

   /* output 6 */
   CARRY_FORWARD;
   SQRADDSC(a[0], a[6]); SQRADDAC(a[1], a[5]); SQRADDAC(a[2], a[4]); SQRADDDB; SQRADD(a[3], a[3]);
   COMBA_STORE(b[6]);

   /* output 7 */
   CARRY_FORWARD;
   SQRADDSC(a[0], a[7]); SQRADDAC(a[1], a[6]); SQRADDAC(a[2], a[5]); SQRADDAC(a[3], a[4]); SQRADDDB;
   COMBA_STORE(b[7]);

   /* output 8 */
   CARRY_FORWARD;
   SQRADDSC(a[1], a[7]); SQRADDAC(a[2], a[6]); SQRADDAC(a[3], a[5]); SQRADDDB; SQRADD(a[4], a[4]);
   COMBA_STORE(b[8]);

   /* output 9 */
   CARRY_FORWARD;
   SQRADDSC(a[2], a[7]); SQRADDAC(a[3], a[6]); SQRADDAC(a[4], a[5]); SQRADDDB;
   COMBA_STORE(b[9]);

   /* output 10 */
   CARRY_FORWARD;
   SQRADD2(a[3], a[7]); SQRADD2(a[4], a[6]); SQRADD(a[5], a[5]);
   COMBA_STORE(b[10]);

   /* output 11 */
   CARRY_FORWARD;
   SQRADD2(a[4], a[7]); SQRADD2(a[5], a[6]);
   COMBA_STORE(b[11]);

   /* output 12 */
   CARRY_FORWARD;
   SQRADD2(a[5], a[7]); SQRADD(a[6], a[6]);
   COMBA_STORE(b[12]);

   /* output 13 */
   CARRY_FORWARD;
   SQRADD2(a[6], a[7]);
   COMBA_STORE(b[13]);

   /* output 14 */
   CARRY_FORWARD;
   SQRADD(a[7], a[7]);
   COMBA_STORE(b[14]);
   COMBA_STORE2(b[15]);
   COMBA_FINI;

   B->used = 16;
   B->sign = FP_ZPOS;
   memcpy(B->dp, b, 16 * sizeof(fp_digit));
   fp_clamp(B);
}
#endif
