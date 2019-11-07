/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
#include <tfm_private.h>

int fp_cmp(fp_int *a, fp_int *b)
{
   if (a->sign == FP_NEG && b->sign == FP_ZPOS) {
      return FP_LT;
   } else if (a->sign == FP_ZPOS && b->sign == FP_NEG) {
      return FP_GT;
   } else {
      /* compare digits */
      if (a->sign == FP_NEG) {
         /* if negative compare opposite direction */
         return fp_cmp_mag(b, a);
      } else {
         return fp_cmp_mag(a, b);
      }
   }
}
