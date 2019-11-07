/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
#include <tfm_private.h>

int fp_isprime(fp_int *a)
{
  return fp_isprime_ex(a, 8);
}
