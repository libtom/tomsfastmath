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
#ifdef TFM_NO_STDLIB
FP_PRIVATE void fp_memcpy(void *restrict dst, const void *restrict src, size_t n);
FP_PRIVATE void fp_memset(void *restrict dst, unsigned char c, size_t n);

#define toupper(C) __extension__({ \
   int _c = (int)(C); \
   'a' <= _c && 'z' >= _c ? _c - ((int)'a' - (int)'A') : _c; \
})

#define memcpy(D, S, N) __extension__({ \
   uint8_t *_dst = (uint8_t*)(D); \
   const uint8_t *_src = (uint8_t*)(S); \
   size_t _n = (size_t)(N); \
   if (__builtin_constant_p(_n)) { \
      for (size_t _i = 0; _i < _n; ++_i) _dst[_i] = _src[_i]; \
   } else { \
      fp_memcpy(_dst, _src, _n); \
   } \
   (void *)_dst; \
})

/* The loop increment of 128 bytes was determined experimentally - it results
 * in good inline code optimizations in both GCC and clang. */
#define memset(D, C, N) __extension__({ \
   uint8_t *_dst = (uint8_t*)(D); \
   uint8_t _c = (uint8_t)(C); \
   size_t _n = (size_t)(N); \
   if (__builtin_constant_p(_n)) { \
      size_t _b = _n >> 7, _r = _n & 127; \
      for (size_t _i = 0; _i < _b; ++_i) { \
         size_t _j = 0; \
         for (;_j <  64; ++_j) *_dst++ = 0; \
         for (;_j < 128; ++_j) *_dst++ = 0; \
      } \
      size_t _i = 0; \
      if (_r >  0) for (;_i <  64 && _i < _r; ++_i) *_dst++ = _c; \
      if (_r > 64) for (;_i < 128 && _i < _r; ++_i) *_dst++ = _c; \
   } else { \
      fp_memset(_dst, _c, _n); \
   } \
   (void *)_dst; \
})
#endif

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
