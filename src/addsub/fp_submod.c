/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
#include <tfm_private.h>

/* d = a - b (mod c) */
int fp_submod(const fp_int *a, const fp_int *b, const fp_int *c, fp_int *d)
{
  fp_int tmp;
  fp_zero(&tmp);
  fp_sub(a, b, &tmp);
  return fp_mod(&tmp, c, d);
}
