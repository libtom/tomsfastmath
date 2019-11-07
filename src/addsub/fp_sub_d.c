/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
#include <tfm_private.h>

/* c = a - b */
void fp_sub_d(fp_int *a, fp_digit b, fp_int *c)
{
   fp_int tmp;
   fp_set(&tmp, b);
   fp_sub(a, &tmp, c);
}
