/* TomsFastMath, a fast ISO C bignum library.
 *
 * This project is meant to fill in where LibTomMath
 * falls short.  That is speed ;-)
 *
 * This project is public domain and free for all purposes.
 *
 * Tom St Denis, tomstdenis@gmail.com
 */
#include <tfm.h>

#define fp_on_bitnum(a, bitnum) \
    a->dp[(bitnum) >> DIGIT_SHIFT] |= (fp_digit)1 << ((bitnum) & (DIGIT_BIT-1))

#define fp_off_bitnum(a, bitnum) \
    a->dp[(bitnum) >> DIGIT_SHIFT] &= ~((fp_digit)1 << ((bitnum) & (DIGIT_BIT-1)))

/* This is possibly the mother of all prime generation functions, muahahahahaha! */
int fp_prime_random_ex(fp_int *a, int size, int flags, tfm_prime_callback cb, void *dat)
{
   fp_digit maskAND_msb, maskOR_lsb;
   int res, err, bsize, dsize;

   /* sanity check the input */
   if (size <= 1) {
      return FP_VAL;
   }

   /* TFM_PRIME_SAFE implies TFM_PRIME_BBS */
   if (flags & TFM_PRIME_SAFE) {
      flags |= TFM_PRIME_BBS;
   }

   /* calc the digit size */
   dsize = (size + DIGIT_BIT - 1) >> DIGIT_SHIFT;

   /* calc the maskAND value for the MSbyte */
   maskAND_msb = FP_MASK >> ((DIGIT_BIT - (size & (DIGIT_BIT-1))) & (DIGIT_BIT-1));

   /* get the maskOR_lsb */
   maskOR_lsb         = 1;
   if (flags & TFM_PRIME_BBS) {
      maskOR_lsb     |= 3;
   }

   do {
      /* read the bytes */
      if (cb((unsigned char*)&a->dp[0], dsize*DIGIT_BIT, dat) != dsize*DIGIT_BIT) {
         return FP_VAL;
      }
      a->used = dsize;

      /* make sure the MSbyte has the required number of bits */
      a->dp[dsize-1]    &= maskAND_msb;

      /* modify the LSbyte as requested */
      a->dp[0]          |= maskOR_lsb;

      /* turn on the MSbit to force the requested magnitude */
      fp_on_bitnum(a, size-1);

      /* modify the 2nd MSBit */
      if (flags & TFM_PRIME_2MSB_ON) {
          fp_on_bitnum(a, size-2);
      } else if (flags & TFM_PRIME_2MSB_OFF) {
          fp_off_bitnum(a, size-2);
      }

      /* is it prime? */
      res = fp_isprime(a);
      if (res == FP_NO) continue;

      if (flags & TFM_PRIME_SAFE) {
         /* see if (a-1)/2 is prime */
         fp_sub_d(a, 1, a);
         fp_div_2(a, a);

         /* is it prime? */
         res = fp_isprime(a);
      }
   } while (res == FP_NO);

   if (flags & TFM_PRIME_SAFE) {
      /* restore a to the original value */
      fp_mul_2(a, a);
      fp_add_d(a, 1, a);
   }

   return FP_OKAY;
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */
