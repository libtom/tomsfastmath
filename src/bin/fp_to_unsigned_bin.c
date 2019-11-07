/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
#include <tfm_private.h>

void fp_to_unsigned_bin(fp_int *a, unsigned char *b)
{
  int     x;
  fp_int  t;

  fp_init_copy(&t, a);

  x = 0;
  while (fp_iszero (&t) == FP_NO) {
      b[x++] = (unsigned char) (t.dp[0] & 255);
      fp_div_2d (&t, 8, &t, NULL);
  }
  fp_reverse (b, x);
}
