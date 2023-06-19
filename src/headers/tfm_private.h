/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
#ifndef TFM_PRIVATE_H_
#define TFM_PRIVATE_H_

#include "tfm.h"

/*
 * Private symbols
 * ---------------
 *
 * On Unix symbols can be marked as hidden if tomsfastmath is compiled
 * as a shared object. By default, symbols are visible.
 */
#if defined(__GNUC__) && __GNUC__ >= 4 && !defined(_WIN32) && !defined(__CYGWIN__)
#   define FP_PRIVATE __attribute__ ((visibility ("hidden")))
#else
#   define FP_PRIVATE
#endif

/* VARIOUS LOW LEVEL STUFFS */
FP_PRIVATE void s_fp_add(fp_int *a, fp_int *b, fp_int *c);
FP_PRIVATE void s_fp_sub(fp_int *a, fp_int *b, fp_int *c);
FP_PRIVATE void fp_reverse(unsigned char *s, int len);

FP_PRIVATE void fp_mul_comba(fp_int *A, fp_int *B, fp_int *C);

#ifdef TFM_SMALL_SET
FP_PRIVATE void fp_mul_comba_small(fp_int *A, fp_int *B, fp_int *C);
#endif

#ifdef TFM_MUL3
FP_PRIVATE void fp_mul_comba3(fp_int *A, fp_int *B, fp_int *C);
#endif
#ifdef TFM_MUL4
FP_PRIVATE void fp_mul_comba4(fp_int *A, fp_int *B, fp_int *C);
#endif
#ifdef TFM_MUL6
FP_PRIVATE void fp_mul_comba6(fp_int *A, fp_int *B, fp_int *C);
#endif
#ifdef TFM_MUL7
FP_PRIVATE void fp_mul_comba7(fp_int *A, fp_int *B, fp_int *C);
#endif
#ifdef TFM_MUL8
FP_PRIVATE void fp_mul_comba8(fp_int *A, fp_int *B, fp_int *C);
#endif
#ifdef TFM_MUL9
FP_PRIVATE void fp_mul_comba9(fp_int *A, fp_int *B, fp_int *C);
#endif
#ifdef TFM_MUL12
FP_PRIVATE void fp_mul_comba12(fp_int *A, fp_int *B, fp_int *C);
#endif
#ifdef TFM_MUL17
FP_PRIVATE void fp_mul_comba17(fp_int *A, fp_int *B, fp_int *C);
#endif

#ifdef TFM_MUL20
FP_PRIVATE void fp_mul_comba20(fp_int *A, fp_int *B, fp_int *C);
#endif
#ifdef TFM_MUL24
FP_PRIVATE void fp_mul_comba24(fp_int *A, fp_int *B, fp_int *C);
#endif
#ifdef TFM_MUL28
FP_PRIVATE void fp_mul_comba28(fp_int *A, fp_int *B, fp_int *C);
#endif
#ifdef TFM_MUL32
FP_PRIVATE void fp_mul_comba32(fp_int *A, fp_int *B, fp_int *C);
#endif
#ifdef TFM_MUL48
FP_PRIVATE void fp_mul_comba48(fp_int *A, fp_int *B, fp_int *C);
#endif
#ifdef TFM_MUL64
FP_PRIVATE void fp_mul_comba64(fp_int *A, fp_int *B, fp_int *C);
#endif

FP_PRIVATE void fp_sqr_comba(fp_int *A, fp_int *B);

#ifdef TFM_SMALL_SET
FP_PRIVATE void fp_sqr_comba_small(fp_int *A, fp_int *B);
#endif

#ifdef TFM_SQR3
FP_PRIVATE void fp_sqr_comba3(fp_int *A, fp_int *B);
#endif
#ifdef TFM_SQR4
FP_PRIVATE void fp_sqr_comba4(fp_int *A, fp_int *B);
#endif
#ifdef TFM_SQR6
FP_PRIVATE void fp_sqr_comba6(fp_int *A, fp_int *B);
#endif
#ifdef TFM_SQR7
FP_PRIVATE void fp_sqr_comba7(fp_int *A, fp_int *B);
#endif
#ifdef TFM_SQR8
FP_PRIVATE void fp_sqr_comba8(fp_int *A, fp_int *B);
#endif
#ifdef TFM_SQR9
FP_PRIVATE void fp_sqr_comba9(fp_int *A, fp_int *B);
#endif
#ifdef TFM_SQR12
FP_PRIVATE void fp_sqr_comba12(fp_int *A, fp_int *B);
#endif
#ifdef TFM_SQR17
FP_PRIVATE void fp_sqr_comba17(fp_int *A, fp_int *B);
#endif

#ifdef TFM_SQR20
FP_PRIVATE void fp_sqr_comba20(fp_int *A, fp_int *B);
#endif
#ifdef TFM_SQR24
FP_PRIVATE void fp_sqr_comba24(fp_int *A, fp_int *B);
#endif
#ifdef TFM_SQR28
FP_PRIVATE void fp_sqr_comba28(fp_int *A, fp_int *B);
#endif
#ifdef TFM_SQR32
FP_PRIVATE void fp_sqr_comba32(fp_int *A, fp_int *B);
#endif
#ifdef TFM_SQR48
FP_PRIVATE void fp_sqr_comba48(fp_int *A, fp_int *B);
#endif
#ifdef TFM_SQR64
FP_PRIVATE void fp_sqr_comba64(fp_int *A, fp_int *B);
#endif
FP_PRIVATE extern const char *fp_s_rmap;

#endif
