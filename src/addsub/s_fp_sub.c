/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
#include <tfm_private.h>

/* unsigned subtraction ||a|| >= ||b|| ALWAYS! */
void s_fp_sub(fp_int *a, fp_int *b, fp_int *c)
{
  int      x, oldbused, oldused;
  fp_word  t;

  oldused  = c->used;
  oldbused = b->used;
  c->used  = a->used;
  t       = 0;
  for (x = 0; x < oldbused; x++) {
     t         = ((fp_word)a->dp[x]) - (((fp_word)b->dp[x]) + t);
     c->dp[x]  = (fp_digit)t;
     t         = (t >> DIGIT_BIT)&1;
  }
  for (; x < a->used; x++) {
     t         = ((fp_word)a->dp[x]) - t;
     c->dp[x]  = (fp_digit)t;
     t         = (t >> DIGIT_BIT)&1;
   }
  for (; x < oldused; x++) {
     c->dp[x] = 0;
  }
  fp_clamp(c);
}
