/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
#include <tfm_private.h>

/* c = a * a (mod b) */
int fp_sqrmod(fp_int *a, fp_int *b, fp_int *c)
{
  fp_int tmp;
  fp_zero(&tmp);
  fp_sqr(a, &tmp);
  return fp_mod(&tmp, b, c);
}
