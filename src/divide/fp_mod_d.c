/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
#include <tfm_private.h>

/* c = a mod b, 0 <= c < b  */
int fp_mod_d(fp_int *a, fp_digit b, fp_digit *c)
{
   return fp_div_d(a, b, NULL, c);
}
