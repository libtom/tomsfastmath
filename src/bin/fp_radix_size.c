/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
#include <tfm_private.h>

int fp_radix_size(const fp_int *a, int radix, int *size)
{
  fp_int  t;
  fp_digit d;

  *size = 0;

  /* check range of the radix */
  if (radix < 2 || radix > 64) {
    return FP_VAL;
  }

  /* quick out if its zero */
  if (fp_iszero(a) == 1) {
     *size = 2;
     return FP_OKAY;
  }

  fp_init_copy(&t, a);

  /* if it is negative output a - */
  if (t.sign == FP_NEG) {
    (*size)++;
    t.sign = FP_ZPOS;
  }

  while (fp_iszero (&t) == FP_NO) {
    fp_div_d (&t, (fp_digit) radix, &t, &d);
    (*size)++;
  }

  /* append a NULL so the string is properly terminated */
  (*size)++;
  return FP_OKAY;

}
