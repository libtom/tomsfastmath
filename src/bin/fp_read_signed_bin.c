/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
#include <tfm_private.h>

void fp_read_signed_bin(fp_int *a, const unsigned char *b, int c)
{
  /* read magnitude */
  fp_read_unsigned_bin (a, b + 1, c - 1);

  /* first byte is 0 for positive, non-zero for negative */
  if (b[0] == 0) {
     a->sign = FP_ZPOS;
  } else {
     a->sign = FP_NEG;
  }
}
