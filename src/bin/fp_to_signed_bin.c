/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
#include <tfm_private.h>

void fp_to_signed_bin(fp_int *a, unsigned char *b)
{
  fp_to_unsigned_bin (a, b + 1);
  b[0] = (unsigned char) ((a->sign == FP_ZPOS) ? 0 : 1);
}
