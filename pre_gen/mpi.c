/* Start: fp_2expt.c */
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

/* computes a = 2**b */
void fp_2expt(fp_int *a, int b)
{
   int     z;

   /* zero a as per default */
   fp_zero (a);

   if (b < 0) { 
      return;
   }

   z = b / DIGIT_BIT;
   if (z >= FP_SIZE) {
      return; 
   }

  /* set the used count of where the bit will go */
  a->used = z + 1;

  /* put the single bit in its place */
  a->dp[z] = ((fp_digit)1) << (b % DIGIT_BIT);
}


/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_2expt.c */

/* Start: fp_add.c */
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

void fp_add(fp_int *a, fp_int *b, fp_int *c)
{
  int     sa, sb;

  /* get sign of both inputs */
  sa = a->sign;
  sb = b->sign;

  /* handle two cases, not four */
  if (sa == sb) {
    /* both positive or both negative */
    /* add their magnitudes, copy the sign */
    c->sign = sa;
    s_fp_add (a, b, c);
  } else {
    /* one positive, the other negative */
    /* subtract the one with the greater magnitude from */
    /* the one of the lesser magnitude.  The result gets */
    /* the sign of the one with the greater magnitude. */
    if (fp_cmp_mag (a, b) == FP_LT) {
      c->sign = sb;
      s_fp_sub (b, a, c);
    } else {
      c->sign = sa;
      s_fp_sub (a, b, c);
    }
  }
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_add.c */

/* Start: fp_add_d.c */
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

/* c = a + b */
void fp_add_d(fp_int *a, fp_digit b, fp_int *c)
{
   fp_int tmp;
   fp_set(&tmp, b);
   fp_add(a,&tmp,c);
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_add_d.c */

/* Start: fp_addmod.c */
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

/* d = a + b (mod c) */
int fp_addmod(fp_int *a, fp_int *b, fp_int *c, fp_int *d)
{
  fp_int tmp;
  fp_zero(&tmp);
  fp_add(a, b, &tmp);
  return fp_mod(&tmp, c, d);
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_addmod.c */

/* Start: fp_cmp.c */
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

int fp_cmp(fp_int *a, fp_int *b)
{
   if (a->sign == FP_NEG && b->sign == FP_ZPOS) {
      return FP_LT;
   } else if (a->sign == FP_ZPOS && b->sign == FP_NEG) {
      return FP_GT;
   } else {
      /* compare digits */
      if (a->sign == FP_NEG) {
         /* if negative compare opposite direction */
         return fp_cmp_mag(b, a);
      } else {
         return fp_cmp_mag(a, b);
      }
   }
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_cmp.c */

/* Start: fp_cmp_d.c */
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

/* compare against a single digit */
int fp_cmp_d(fp_int *a, fp_digit b)
{
  /* compare based on sign */
  if ((b && a->used == 0) || a->sign == FP_NEG) {
    return FP_LT;
  }

  /* compare based on magnitude */
  if (a->used > 1) {
    return FP_GT;
  }

  /* compare the only digit of a to b */
  if (a->dp[0] > b) {
    return FP_GT;
  } else if (a->dp[0] < b) {
    return FP_LT;
  } else {
    return FP_EQ;
  }

}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_cmp_d.c */

/* Start: fp_cmp_mag.c */
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

int fp_cmp_mag(fp_int *a, fp_int *b)
{
   int x;

   if (a->used > b->used) {
      return FP_GT;
   } else if (a->used < b->used) {
      return FP_LT;
   } else {
      for (x = a->used - 1; x >= 0; x--) {
          if (a->dp[x] > b->dp[x]) {
             return FP_GT;
          } else if (a->dp[x] < b->dp[x]) {
             return FP_LT;
          }
      }
   }
   return FP_EQ;
}


/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_cmp_mag.c */

/* Start: fp_cnt_lsb.c */
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

static const int lnz[16] = {
   4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0
};

/* Counts the number of lsbs which are zero before the first zero bit */
int fp_cnt_lsb(fp_int *a)
{
   int x;
   fp_digit q, qq;

   /* easy out */
   if (fp_iszero(a) == 1) {
      return 0;
   }

   /* scan lower digits until non-zero */
   for (x = 0; x < a->used && a->dp[x] == 0; x++);
   q = a->dp[x];
   x *= DIGIT_BIT;

   /* now scan this digit until a 1 is found */
   if ((q & 1) == 0) {
      do {
         qq  = q & 15;
         x  += lnz[qq];
         q >>= 4;
      } while (qq == 0);
   }
   return x;
}


/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_cnt_lsb.c */

/* Start: fp_count_bits.c */
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

int fp_count_bits (fp_int * a)
{
  int     r;
  fp_digit q;

  /* shortcut */
  if (a->used == 0) {
    return 0;
  }

  /* get number of digits and add that */
  r = (a->used - 1) * DIGIT_BIT;

  /* take the last digit and count the bits in it */
  q = a->dp[a->used - 1];
  while (q > ((fp_digit) 0)) {
    ++r;
    q >>= ((fp_digit) 1);
  }
  return r;
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_count_bits.c */

/* Start: fp_div.c */
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

/* a/b => cb + d == a */
int fp_div(fp_int *a, fp_int *b, fp_int *c, fp_int *d)
{
  fp_int  q, x, y, t1, t2;
  int     n, t, i, norm, neg;

  /* is divisor zero ? */
  if (fp_iszero (b) == 1) {
    return FP_VAL;
  }

  /* if a < b then q=0, r = a */
  if (fp_cmp_mag (a, b) == FP_LT) {
    if (d != NULL) {
      fp_copy (a, d);
    } 
    if (c != NULL) {
      fp_zero (c);
    }
    return FP_OKAY;
  }

  fp_init(&q);
  q.used = a->used + 2;

  fp_init(&t1);
  fp_init(&t2);
  fp_init_copy(&x, a);
  fp_init_copy(&y, b);

  /* fix the sign */
  neg = (a->sign == b->sign) ? FP_ZPOS : FP_NEG;
  x.sign = y.sign = FP_ZPOS;

  /* normalize both x and y, ensure that y >= b/2, [b == 2**DIGIT_BIT] */
  norm = fp_count_bits(&y) % DIGIT_BIT;
  if (norm < (int)(DIGIT_BIT-1)) {
     norm = (DIGIT_BIT-1) - norm;
     fp_mul_2d (&x, norm, &x);
     fp_mul_2d (&y, norm, &y);
  } else {
     norm = 0;
  }

  /* note hac does 0 based, so if used==5 then its 0,1,2,3,4, e.g. use 4 */
  n = x.used - 1;
  t = y.used - 1;

  /* while (x >= y*b**n-t) do { q[n-t] += 1; x -= y*b**{n-t} } */
  fp_lshd (&y, n - t);                                             /* y = y*b**{n-t} */

  while (fp_cmp (&x, &y) != FP_LT) {
    ++(q.dp[n - t]);
    fp_sub (&x, &y, &x);
  }

  /* reset y by shifting it back down */
  fp_rshd (&y, n - t);

  /* step 3. for i from n down to (t + 1) */
  for (i = n; i >= (t + 1); i--) {
    if (i > x.used) {
      continue;
    }

    /* step 3.1 if xi == yt then set q{i-t-1} to b-1, 
     * otherwise set q{i-t-1} to (xi*b + x{i-1})/yt */
    if (x.dp[i] == y.dp[t]) {
      q.dp[i - t - 1] = ((((fp_word)1) << DIGIT_BIT) - 1);
    } else {
      fp_word tmp;
      tmp = ((fp_word) x.dp[i]) << ((fp_word) DIGIT_BIT);
      tmp |= ((fp_word) x.dp[i - 1]);
      tmp /= ((fp_word) y.dp[t]);
      q.dp[i - t - 1] = (fp_digit) (tmp);
    }

    /* while (q{i-t-1} * (yt * b + y{t-1})) > 
             xi * b**2 + xi-1 * b + xi-2 
     
       do q{i-t-1} -= 1; 
    */
    q.dp[i - t - 1] = (q.dp[i - t - 1] + 1);
    do {
      q.dp[i - t - 1] = (q.dp[i - t - 1] - 1);

      /* find left hand */
      fp_zero (&t1);
      t1.dp[0] = (t - 1 < 0) ? 0 : y.dp[t - 1];
      t1.dp[1] = y.dp[t];
      t1.used = 2;
      fp_mul_d (&t1, q.dp[i - t - 1], &t1);

      /* find right hand */
      t2.dp[0] = (i - 2 < 0) ? 0 : x.dp[i - 2];
      t2.dp[1] = (i - 1 < 0) ? 0 : x.dp[i - 1];
      t2.dp[2] = x.dp[i];
      t2.used = 3;
    } while (fp_cmp_mag(&t1, &t2) == FP_GT);

    /* step 3.3 x = x - q{i-t-1} * y * b**{i-t-1} */
    fp_mul_d (&y, q.dp[i - t - 1], &t1);
    fp_lshd  (&t1, i - t - 1);
    fp_sub   (&x, &t1, &x);

    /* if x < 0 then { x = x + y*b**{i-t-1}; q{i-t-1} -= 1; } */
    if (x.sign == FP_NEG) {
      fp_copy (&y, &t1);
      fp_lshd (&t1, i - t - 1);
      fp_add (&x, &t1, &x);
      q.dp[i - t - 1] = q.dp[i - t - 1] - 1;
    }
  }

  /* now q is the quotient and x is the remainder 
   * [which we have to normalize] 
   */
  
  /* get sign before writing to c */
  x.sign = x.used == 0 ? FP_ZPOS : a->sign;

  if (c != NULL) {
    fp_clamp (&q);
    fp_copy (&q, c);
    c->sign = neg;
  }

  if (d != NULL) {
    fp_div_2d (&x, norm, &x, NULL);

/* the following is a kludge, essentially we were seeing the right remainder but 
   with excess digits that should have been zero
 */
    for (i = b->used; i < x.used; i++) {
        x.dp[i] = 0;
    }
    fp_clamp(&x);
    fp_copy (&x, d);
  }

  return FP_OKAY;
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_div.c */

/* Start: fp_div_2.c */
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

/* b = a/2 */
void fp_div_2(fp_int * a, fp_int * b)
{
  int     x, oldused;

  oldused = b->used;
  b->used = a->used;
  {
    register fp_digit r, rr, *tmpa, *tmpb;

    /* source alias */
    tmpa = a->dp + b->used - 1;

    /* dest alias */
    tmpb = b->dp + b->used - 1;

    /* carry */
    r = 0;
    for (x = b->used - 1; x >= 0; x--) {
      /* get the carry for the next iteration */
      rr = *tmpa & 1;

      /* shift the current digit, add in carry and store */
      *tmpb-- = (*tmpa-- >> 1) | (r << (DIGIT_BIT - 1));

      /* forward carry to next iteration */
      r = rr;
    }

    /* zero excess digits */
    tmpb = b->dp + b->used;
    for (x = b->used; x < oldused; x++) {
      *tmpb++ = 0;
    }
  }
  b->sign = a->sign;
  fp_clamp (b);
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_div_2.c */

/* Start: fp_div_2d.c */
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

/* c = a / 2**b */
void fp_div_2d(fp_int *a, int b, fp_int *c, fp_int *d)
{
  fp_digit D, r, rr;
  int      x;
  fp_int   t;

  /* if the shift count is <= 0 then we do no work */
  if (b <= 0) {
    fp_copy (a, c);
    if (d != NULL) {
      fp_zero (d);
    }
    return;
  }

  fp_init(&t);

  /* get the remainder */
  if (d != NULL) {
    fp_mod_2d (a, b, &t);
  }

  /* copy */
  fp_copy(a, c);

  /* shift by as many digits in the bit count */
  if (b >= (int)DIGIT_BIT) {
    fp_rshd (c, b / DIGIT_BIT);
  }

  /* shift any bit count < DIGIT_BIT */
  D = (fp_digit) (b % DIGIT_BIT);
  if (D != 0) {
    register fp_digit *tmpc, mask, shift;

    /* mask */
    mask = (((fp_digit)1) << D) - 1;

    /* shift for lsb */
    shift = DIGIT_BIT - D;

    /* alias */
    tmpc = c->dp + (c->used - 1);

    /* carry */
    r = 0;
    for (x = c->used - 1; x >= 0; x--) {
      /* get the lower  bits of this word in a temp */
      rr = *tmpc & mask;

      /* shift the current word and mix in the carry bits from the previous word */
      *tmpc = (*tmpc >> D) | (r << shift);
      --tmpc;

      /* set the carry to the carry bits of the current word found above */
      r = rr;
    }
  }
  fp_clamp (c);
  if (d != NULL) {
    fp_copy (&t, d);
  }
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_div_2d.c */

/* Start: fp_div_d.c */
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

static int s_is_power_of_two(fp_digit b, int *p)
{
   int x;

   for (x = 1; x < DIGIT_BIT; x++) {
      if (b == (((fp_digit)1)<<x)) {
         *p = x;
         return 1;
      }
   }
   return 0;
}

/* a/b => cb + d == a */
int fp_div_d(fp_int *a, fp_digit b, fp_int *c, fp_digit *d)
{
  fp_int   q;
  fp_word  w;
  fp_digit t;
  int      ix;

  /* cannot divide by zero */
  if (b == 0) {
     return FP_VAL;
  }

  /* quick outs */
  if (b == 1 || fp_iszero(a) == 1) {
     if (d != NULL) {
        *d = 0;
     }
     if (c != NULL) {
        fp_copy(a, c);
     }
     return FP_OKAY;
  }

  /* power of two ? */
  if (s_is_power_of_two(b, &ix) == 1) {
     if (d != NULL) {
        *d = a->dp[0] & ((((fp_digit)1)<<ix) - 1);
     }
     if (c != NULL) {
        fp_div_2d(a, ix, c, NULL);
     }
     return FP_OKAY;
  }

  /* no easy answer [c'est la vie].  Just division */
  fp_init(&q);
  
  q.used = a->used;
  q.sign = a->sign;
  w = 0;
  for (ix = a->used - 1; ix >= 0; ix--) {
     w = (w << ((fp_word)DIGIT_BIT)) | ((fp_word)a->dp[ix]);
     
     if (w >= b) {
        t = (fp_digit)(w / b);
        w -= ((fp_word)t) * ((fp_word)b);
      } else {
        t = 0;
      }
      q.dp[ix] = (fp_digit)t;
  }
  
  if (d != NULL) {
     *d = (fp_digit)w;
  }
  
  if (c != NULL) {
     fp_clamp(&q);
     fp_copy(&q, c);
  }
 
  return FP_OKAY;
}


/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_div_d.c */

/* Start: fp_exptmod.c */
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

#ifdef TFM_TIMING_RESISTANT

/* timing resistant montgomery ladder based exptmod 

   Based on work by Marc Joye, Sung-Ming Yen, "The Montgomery Powering Ladder", Cryptographic Hardware and Embedded Systems, CHES 2002
*/
static int _fp_exptmod(fp_int * G, fp_int * X, fp_int * P, fp_int * Y)
{
  fp_int   R[2];
  fp_digit buf, mp;
  int      err, bitcnt, digidx, y;

  /* now setup montgomery  */
  if ((err = fp_montgomery_setup (P, &mp)) != FP_OKAY) {
     return err;
  }

  fp_init(&R[0]);   
  fp_init(&R[1]);   
   
  /* now we need R mod m */
  fp_montgomery_calc_normalization (&R[0], P);

  /* now set R[0][1] to G * R mod m */
  if (fp_cmp_mag(P, G) != FP_GT) {
     /* G > P so we reduce it first */
     fp_mod(G, P, &R[1]);
  } else {
     fp_copy(G, &R[1]);
  }
  fp_mulmod (&R[1], &R[0], P, &R[1]);

  /* for j = t-1 downto 0 do
        r_!k = R0*R1; r_k = r_k^2
  */
  
  /* set initial mode and bit cnt */
  bitcnt = 1;
  buf    = 0;
  digidx = X->used - 1;

  for (;;) {
    /* grab next digit as required */
    if (--bitcnt == 0) {
      /* if digidx == -1 we are out of digits so break */
      if (digidx == -1) {
        break;
      }
      /* read next digit and reset bitcnt */
      buf    = X->dp[digidx--];
      bitcnt = (int)DIGIT_BIT;
    }

    /* grab the next msb from the exponent */
    y     = (fp_digit)(buf >> (DIGIT_BIT - 1)) & 1;
    buf <<= (fp_digit)1;

    /* do ops */
    fp_mul(&R[0], &R[1], &R[y^1]); fp_montgomery_reduce(&R[y^1], P, mp);
    fp_sqr(&R[y], &R[y]);          fp_montgomery_reduce(&R[y], P, mp);
  }

   fp_montgomery_reduce(&R[0], P, mp);
   fp_copy(&R[0], Y);
   return FP_OKAY;
}   

#else

/* y = g**x (mod b) 
 * Some restrictions... x must be positive and < b
 */
static int _fp_exptmod(fp_int * G, fp_int * X, fp_int * P, fp_int * Y)
{
  fp_int   M[64], res;
  fp_digit buf, mp;
  int      err, bitbuf, bitcpy, bitcnt, mode, digidx, x, y, winsize;

  /* find window size */
  x = fp_count_bits (X);
  if (x <= 21) {
    winsize = 1;
  } else if (x <= 36) {
    winsize = 3;
  } else if (x <= 140) {
    winsize = 4;
  } else if (x <= 450) {
    winsize = 5;
  } else {
    winsize = 6;
  } 

  /* init M array */
  memset(M, 0, sizeof(M)); 

  /* now setup montgomery  */
  if ((err = fp_montgomery_setup (P, &mp)) != FP_OKAY) {
     return err;
  }

  /* setup result */
  fp_init(&res);

  /* create M table
   *
   * The M table contains powers of the input base, e.g. M[x] = G^x mod P
   *
   * The first half of the table is not computed though accept for M[0] and M[1]
   */

   /* now we need R mod m */
   fp_montgomery_calc_normalization (&res, P);

   /* now set M[1] to G * R mod m */
   if (fp_cmp_mag(P, G) != FP_GT) {
      /* G > P so we reduce it first */
      fp_mod(G, P, &M[1]);
   } else {
      fp_copy(G, &M[1]);
   }
   fp_mulmod (&M[1], &res, P, &M[1]);

  /* compute the value at M[1<<(winsize-1)] by squaring M[1] (winsize-1) times */
  fp_copy (&M[1], &M[1 << (winsize - 1)]);
  for (x = 0; x < (winsize - 1); x++) {
    fp_sqr (&M[1 << (winsize - 1)], &M[1 << (winsize - 1)]);
    fp_montgomery_reduce (&M[1 << (winsize - 1)], P, mp);
  }

  /* create upper table */
  for (x = (1 << (winsize - 1)) + 1; x < (1 << winsize); x++) {
    fp_mul(&M[x - 1], &M[1], &M[x]);
    fp_montgomery_reduce(&M[x], P, mp);
  }

  /* set initial mode and bit cnt */
  mode   = 0;
  bitcnt = 1;
  buf    = 0;
  digidx = X->used - 1;
  bitcpy = 0;
  bitbuf = 0;

  for (;;) {
    /* grab next digit as required */
    if (--bitcnt == 0) {
      /* if digidx == -1 we are out of digits so break */
      if (digidx == -1) {
        break;
      }
      /* read next digit and reset bitcnt */
      buf    = X->dp[digidx--];
      bitcnt = (int)DIGIT_BIT;
    }

    /* grab the next msb from the exponent */
    y     = (fp_digit)(buf >> (DIGIT_BIT - 1)) & 1;
    buf <<= (fp_digit)1;

    /* if the bit is zero and mode == 0 then we ignore it
     * These represent the leading zero bits before the first 1 bit
     * in the exponent.  Technically this opt is not required but it
     * does lower the # of trivial squaring/reductions used
     */
    if (mode == 0 && y == 0) {
      continue;
    }

    /* if the bit is zero and mode == 1 then we square */
    if (mode == 1 && y == 0) {
      fp_sqr(&res, &res);
      fp_montgomery_reduce(&res, P, mp);
      continue;
    }

    /* else we add it to the window */
    bitbuf |= (y << (winsize - ++bitcpy));
    mode    = 2;

    if (bitcpy == winsize) {
      /* ok window is filled so square as required and multiply  */
      /* square first */
      for (x = 0; x < winsize; x++) {
        fp_sqr(&res, &res);
        fp_montgomery_reduce(&res, P, mp);
      }

      /* then multiply */
      fp_mul(&res, &M[bitbuf], &res);
      fp_montgomery_reduce(&res, P, mp);

      /* empty window and reset */
      bitcpy = 0;
      bitbuf = 0;
      mode   = 1;
    }
  }

  /* if bits remain then square/multiply */
  if (mode == 2 && bitcpy > 0) {
    /* square then multiply if the bit is set */
    for (x = 0; x < bitcpy; x++) {
      fp_sqr(&res, &res);
      fp_montgomery_reduce(&res, P, mp);

      /* get next bit of the window */
      bitbuf <<= 1;
      if ((bitbuf & (1 << winsize)) != 0) {
        /* then multiply */
        fp_mul(&res, &M[1], &res);
        fp_montgomery_reduce(&res, P, mp);
      }
    }
  }

  /* fixup result if Montgomery reduction is used
   * recall that any value in a Montgomery system is
   * actually multiplied by R mod n.  So we have
   * to reduce one more time to cancel out the factor
   * of R.
   */
  fp_montgomery_reduce(&res, P, mp);

  /* swap res with Y */
  fp_copy (&res, Y);
  return FP_OKAY;
}

#endif


int fp_exptmod(fp_int * G, fp_int * X, fp_int * P, fp_int * Y)
{
   fp_int tmp;
   int    err;
   
#ifdef TFM_CHECK
   /* prevent overflows */
   if (P->used > (FP_SIZE/2)) {
      return FP_VAL;
   }
#endif

   /* is X negative?  */
   if (X->sign == FP_NEG) {
      /* yes, copy G and invmod it */
      fp_copy(G, &tmp);
      if ((err = fp_invmod(&tmp, P, &tmp)) != FP_OKAY) {
         return err;
      }
      X->sign = FP_ZPOS;
      err =  _fp_exptmod(&tmp, X, P, Y);
      if (X != Y) {
         X->sign = FP_NEG;
      }
      return err;
   } else {
      /* Positive exponent so just exptmod */
      return _fp_exptmod(G, X, P, Y);
   }
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_exptmod.c */

/* Start: fp_gcd.c */
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

/* c = (a, b) */
void fp_gcd(fp_int *a, fp_int *b, fp_int *c)
{
   fp_int u, v, r;

   /* either zero than gcd is the largest */
   if (fp_iszero (a) == 1 && fp_iszero (b) == 0) {
     fp_abs (b, c);
     return;
   }
   if (fp_iszero (a) == 0 && fp_iszero (b) == 1) {
     fp_abs (a, c);
     return;
   }

   /* optimized.  At this point if a == 0 then
    * b must equal zero too
    */
   if (fp_iszero (a) == 1) {
     fp_zero(c);
     return;
   }

   /* sort inputs */
   if (fp_cmp_mag(a, b) != FP_LT) {
      fp_init_copy(&u, a);
      fp_init_copy(&v, b);
   } else {
      fp_init_copy(&u, b);
      fp_init_copy(&v, a);
   }
 
   fp_zero(&r);
   while (fp_iszero(&v) == FP_NO) {
      fp_mod(&u, &v, &r);
      fp_copy(&v, &u);
      fp_copy(&r, &v);
   }
   fp_copy(&u, c);
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_gcd.c */

/* Start: fp_ident.c */
/* TomsFastMath, a fast ISO C bignum library.
 * 
 * This project is meant to fill in where LibTomMath
 * falls short.  That is speed ;-)
 *
 * This project is public domain and free for all purposes.
 * 
 * Tom St Denis, tomstdenis@gmail.com
 */
#include "tfm.h"

const char *fp_ident(void)
{
   static char buf[1024];

   memset(buf, 0, sizeof(buf));
   snprintf(buf, sizeof(buf)-1,
"TomsFastMath (%s)\n"
"\n"
"Sizeofs\n"
"\tfp_digit = %u\n"
"\tfp_word  = %u\n"
"\n"
"FP_MAX_SIZE = %u\n"
"\n"
"Defines: \n"
#ifdef __i386__
" __i386__ "
#endif
#ifdef __x86_64__
" __x86_64__ "
#endif
#ifdef TFM_X86
" TFM_X86 "
#endif
#ifdef TFM_X86_64
" TFM_X86_64 "
#endif
#ifdef TFM_SSE2
" TFM_SSE2 "
#endif
#ifdef TFM_ARM
" TFM_ARM "
#endif
#ifdef TFM_PPC32
" TFM_PPC32 "
#endif
#ifdef TFM_AVR32
" TFM_AVR32 "
#endif
#ifdef TFM_ECC192
" TFM_ECC192 "
#endif
#ifdef TFM_ECC224
" TFM_ECC224 "
#endif
#ifdef TFM_ECC384
" TFM_ECC384 "
#endif
#ifdef TFM_ECC521
" TFM_ECC521 "
#endif

#ifdef TFM_NO_ASM
" TFM_NO_ASM "
#endif
#ifdef FP_64BIT
" FP_64BIT "
#endif
#ifdef TFM_HUGE
" TFM_HUGE "
#endif
"\n", __DATE__, sizeof(fp_digit), sizeof(fp_word), FP_MAX_SIZE);

   if (sizeof(fp_digit) == sizeof(fp_word)) {
      strncat(buf, "WARNING: sizeof(fp_digit) == sizeof(fp_word), this build is likely to not work properly.\n", 
              sizeof(buf)-1);
   }
   return buf;
}

#ifdef STANDALONE

int main(void)
{
   printf("%s\n", fp_ident());
   return 0;
}

#endif


/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_ident.c */

/* Start: fp_invmod.c */
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

static int fp_invmod_slow (fp_int * a, fp_int * b, fp_int * c)
{
  fp_int  x, y, u, v, A, B, C, D;
  int     res;

  /* b cannot be negative */
  if (b->sign == FP_NEG || fp_iszero(b) == 1) {
    return FP_VAL;
  }

  /* init temps */
  fp_init(&x);    fp_init(&y);
  fp_init(&u);    fp_init(&v);
  fp_init(&A);    fp_init(&B);
  fp_init(&C);    fp_init(&D);

  /* x = a, y = b */
  if ((res = fp_mod(a, b, &x)) != FP_OKAY) {
      return res;
  }
  fp_copy(b, &y);

  /* 2. [modified] if x,y are both even then return an error! */
  if (fp_iseven (&x) == 1 && fp_iseven (&y) == 1) {
    return FP_VAL;
  }

  /* 3. u=x, v=y, A=1, B=0, C=0,D=1 */
  fp_copy (&x, &u);
  fp_copy (&y, &v);
  fp_set (&A, 1);
  fp_set (&D, 1);

top:
  /* 4.  while u is even do */
  while (fp_iseven (&u) == 1) {
    /* 4.1 u = u/2 */
    fp_div_2 (&u, &u);

    /* 4.2 if A or B is odd then */
    if (fp_isodd (&A) == 1 || fp_isodd (&B) == 1) {
      /* A = (A+y)/2, B = (B-x)/2 */
      fp_add (&A, &y, &A);
      fp_sub (&B, &x, &B);
    }
    /* A = A/2, B = B/2 */
    fp_div_2 (&A, &A);
    fp_div_2 (&B, &B);
  }

  /* 5.  while v is even do */
  while (fp_iseven (&v) == 1) {
    /* 5.1 v = v/2 */
    fp_div_2 (&v, &v);

    /* 5.2 if C or D is odd then */
    if (fp_isodd (&C) == 1 || fp_isodd (&D) == 1) {
      /* C = (C+y)/2, D = (D-x)/2 */
      fp_add (&C, &y, &C);
      fp_sub (&D, &x, &D);
    }
    /* C = C/2, D = D/2 */
    fp_div_2 (&C, &C);
    fp_div_2 (&D, &D);
  }

  /* 6.  if u >= v then */
  if (fp_cmp (&u, &v) != FP_LT) {
    /* u = u - v, A = A - C, B = B - D */
    fp_sub (&u, &v, &u);
    fp_sub (&A, &C, &A);
    fp_sub (&B, &D, &B);
  } else {
    /* v - v - u, C = C - A, D = D - B */
    fp_sub (&v, &u, &v);
    fp_sub (&C, &A, &C);
    fp_sub (&D, &B, &D);
  }

  /* if not zero goto step 4 */
  if (fp_iszero (&u) == 0)
    goto top;

  /* now a = C, b = D, gcd == g*v */

  /* if v != 1 then there is no inverse */
  if (fp_cmp_d (&v, 1) != FP_EQ) {
    return FP_VAL;
  }

  /* if its too low */
  while (fp_cmp_d(&C, 0) == FP_LT) {
      fp_add(&C, b, &C);
  }
  
  /* too big */
  while (fp_cmp_mag(&C, b) != FP_LT) {
      fp_sub(&C, b, &C);
  }
  
  /* C is now the inverse */
  fp_copy(&C, c);
  return FP_OKAY;
}

/* c = 1/a (mod b) for odd b only */
int fp_invmod(fp_int *a, fp_int *b, fp_int *c)
{
  fp_int  x, y, u, v, B, D;
  int     neg;

  /* 2. [modified] b must be odd   */
  if (fp_iseven (b) == FP_YES) {
    return fp_invmod_slow(a,b,c);
  }

  /* init all our temps */
  fp_init(&x);  fp_init(&y);
  fp_init(&u);  fp_init(&v);
  fp_init(&B);  fp_init(&D);

  /* x == modulus, y == value to invert */
  fp_copy(b, &x);

  /* we need y = |a| */
  fp_abs(a, &y);

  /* 3. u=x, v=y, A=1, B=0, C=0,D=1 */
  fp_copy(&x, &u);
  fp_copy(&y, &v);
  fp_set (&D, 1);

top:
  /* 4.  while u is even do */
  while (fp_iseven (&u) == FP_YES) {
    /* 4.1 u = u/2 */
    fp_div_2 (&u, &u);

    /* 4.2 if B is odd then */
    if (fp_isodd (&B) == FP_YES) {
      fp_sub (&B, &x, &B);
    }
    /* B = B/2 */
    fp_div_2 (&B, &B);
  }

  /* 5.  while v is even do */
  while (fp_iseven (&v) == FP_YES) {
    /* 5.1 v = v/2 */
    fp_div_2 (&v, &v);

    /* 5.2 if D is odd then */
    if (fp_isodd (&D) == FP_YES) {
      /* D = (D-x)/2 */
      fp_sub (&D, &x, &D);
    }
    /* D = D/2 */
    fp_div_2 (&D, &D);
  }

  /* 6.  if u >= v then */
  if (fp_cmp (&u, &v) != FP_LT) {
    /* u = u - v, B = B - D */
    fp_sub (&u, &v, &u);
    fp_sub (&B, &D, &B);
  } else {
    /* v - v - u, D = D - B */
    fp_sub (&v, &u, &v);
    fp_sub (&D, &B, &D);
  }

  /* if not zero goto step 4 */
  if (fp_iszero (&u) == FP_NO) {
    goto top;
  }

  /* now a = C, b = D, gcd == g*v */

  /* if v != 1 then there is no inverse */
  if (fp_cmp_d (&v, 1) != FP_EQ) {
    return FP_VAL;
  }

  /* b is now the inverse */
  neg = a->sign;
  while (D.sign == FP_NEG) {
    fp_add (&D, b, &D);
  }
  fp_copy (&D, c);
  c->sign = neg;
  return FP_OKAY;
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_invmod.c */

/* Start: fp_isprime.c */
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

/* a few primes */
static const fp_digit primes[256] = {
  0x0002, 0x0003, 0x0005, 0x0007, 0x000B, 0x000D, 0x0011, 0x0013,
  0x0017, 0x001D, 0x001F, 0x0025, 0x0029, 0x002B, 0x002F, 0x0035,
  0x003B, 0x003D, 0x0043, 0x0047, 0x0049, 0x004F, 0x0053, 0x0059,
  0x0061, 0x0065, 0x0067, 0x006B, 0x006D, 0x0071, 0x007F, 0x0083,
  0x0089, 0x008B, 0x0095, 0x0097, 0x009D, 0x00A3, 0x00A7, 0x00AD,
  0x00B3, 0x00B5, 0x00BF, 0x00C1, 0x00C5, 0x00C7, 0x00D3, 0x00DF,
  0x00E3, 0x00E5, 0x00E9, 0x00EF, 0x00F1, 0x00FB, 0x0101, 0x0107,
  0x010D, 0x010F, 0x0115, 0x0119, 0x011B, 0x0125, 0x0133, 0x0137,

  0x0139, 0x013D, 0x014B, 0x0151, 0x015B, 0x015D, 0x0161, 0x0167,
  0x016F, 0x0175, 0x017B, 0x017F, 0x0185, 0x018D, 0x0191, 0x0199,
  0x01A3, 0x01A5, 0x01AF, 0x01B1, 0x01B7, 0x01BB, 0x01C1, 0x01C9,
  0x01CD, 0x01CF, 0x01D3, 0x01DF, 0x01E7, 0x01EB, 0x01F3, 0x01F7,
  0x01FD, 0x0209, 0x020B, 0x021D, 0x0223, 0x022D, 0x0233, 0x0239,
  0x023B, 0x0241, 0x024B, 0x0251, 0x0257, 0x0259, 0x025F, 0x0265,
  0x0269, 0x026B, 0x0277, 0x0281, 0x0283, 0x0287, 0x028D, 0x0293,
  0x0295, 0x02A1, 0x02A5, 0x02AB, 0x02B3, 0x02BD, 0x02C5, 0x02CF,

  0x02D7, 0x02DD, 0x02E3, 0x02E7, 0x02EF, 0x02F5, 0x02F9, 0x0301,
  0x0305, 0x0313, 0x031D, 0x0329, 0x032B, 0x0335, 0x0337, 0x033B,
  0x033D, 0x0347, 0x0355, 0x0359, 0x035B, 0x035F, 0x036D, 0x0371,
  0x0373, 0x0377, 0x038B, 0x038F, 0x0397, 0x03A1, 0x03A9, 0x03AD,
  0x03B3, 0x03B9, 0x03C7, 0x03CB, 0x03D1, 0x03D7, 0x03DF, 0x03E5,
  0x03F1, 0x03F5, 0x03FB, 0x03FD, 0x0407, 0x0409, 0x040F, 0x0419,
  0x041B, 0x0425, 0x0427, 0x042D, 0x043F, 0x0443, 0x0445, 0x0449,
  0x044F, 0x0455, 0x045D, 0x0463, 0x0469, 0x047F, 0x0481, 0x048B,

  0x0493, 0x049D, 0x04A3, 0x04A9, 0x04B1, 0x04BD, 0x04C1, 0x04C7,
  0x04CD, 0x04CF, 0x04D5, 0x04E1, 0x04EB, 0x04FD, 0x04FF, 0x0503,
  0x0509, 0x050B, 0x0511, 0x0515, 0x0517, 0x051B, 0x0527, 0x0529,
  0x052F, 0x0551, 0x0557, 0x055D, 0x0565, 0x0577, 0x0581, 0x058F,
  0x0593, 0x0595, 0x0599, 0x059F, 0x05A7, 0x05AB, 0x05AD, 0x05B3,
  0x05BF, 0x05C9, 0x05CB, 0x05CF, 0x05D1, 0x05D5, 0x05DB, 0x05E7,
  0x05F3, 0x05FB, 0x0607, 0x060D, 0x0611, 0x0617, 0x061F, 0x0623,
  0x062B, 0x062F, 0x063D, 0x0641, 0x0647, 0x0649, 0x064D, 0x0653
};

int fp_isprime(fp_int *a)
{
   fp_int   b;
   fp_digit d;
   int      r, res;

   /* do trial division */
   for (r = 0; r < 256; r++) {
       fp_mod_d(a, primes[r], &d);
       if (d == 0) {
          return FP_NO;
       }
   }

   /* now do 8 miller rabins */
   for (r = 0; r < 8; r++) {
       fp_set(&b, primes[r]);
       fp_prime_miller_rabin(a, &b, &res);
       if (res == FP_NO) {
          return FP_NO;
       }
   }
   return FP_YES;
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_isprime.c */

/* Start: fp_lcm.c */
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

/* c = [a, b] */
void fp_lcm(fp_int *a, fp_int *b, fp_int *c)
{
   fp_int  t1, t2;

   fp_init(&t1);
   fp_init(&t2);
   fp_gcd(a, b, &t1);
   if (fp_cmp_mag(a, b) == FP_GT) {
      fp_div(a, &t1, &t2, NULL);
      fp_mul(b, &t2, c);
   } else {
      fp_div(b, &t1, &t2, NULL);
      fp_mul(a, &t2, c);
   }   
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_lcm.c */

/* Start: fp_lshd.c */
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

void fp_lshd(fp_int *a, int x)
{
   int y;

   /* move up and truncate as required */
   y = MIN(a->used + x - 1, (int)(FP_SIZE-1));

   /* store new size */
   a->used = y + 1;

   /* move digits */
   for (; y >= x; y--) {
       a->dp[y] = a->dp[y-x];
   }
 
   /* zero lower digits */
   for (; y >= 0; y--) {
       a->dp[y] = 0;
   }

   /* clamp digits */
   fp_clamp(a);
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_lshd.c */

/* Start: fp_mod.c */
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

/* c = a mod b, 0 <= c < b  */
int fp_mod(fp_int *a, fp_int *b, fp_int *c)
{
   fp_int t;
   int    err;

   fp_zero(&t);
   if ((err = fp_div(a, b, NULL, &t)) != FP_OKAY) {
      return err;
   }
   if (t.sign != b->sign) {
      fp_add(&t, b, c);
   } else {
      fp_copy(&t, c);
  }
  return FP_OKAY;
}



/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_mod.c */

/* Start: fp_mod_2d.c */
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

/* c = a mod 2**d */
void fp_mod_2d(fp_int *a, int b, fp_int *c)
{
   int x;

   /* zero if count less than or equal to zero */
   if (b <= 0) {
      fp_zero(c);
      return;
   }

   /* get copy of input */
   fp_copy(a, c);
 
   /* if 2**d is larger than we just return */
   if (b >= (DIGIT_BIT * a->used)) {
      return;
   }

  /* zero digits above the last digit of the modulus */
  for (x = (b / DIGIT_BIT) + ((b % DIGIT_BIT) == 0 ? 0 : 1); x < c->used; x++) {
    c->dp[x] = 0;
  }
  /* clear the digit that is not completely outside/inside the modulus */
  c->dp[b / DIGIT_BIT] &= ~((fp_digit)0) >> (DIGIT_BIT - b);
  fp_clamp (c);
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_mod_2d.c */

/* Start: fp_mod_d.c */
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

/* c = a mod b, 0 <= c < b  */
int fp_mod_d(fp_int *a, fp_digit b, fp_digit *c)
{
   return fp_div_d(a, b, NULL, c);
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_mod_d.c */

/* Start: fp_mont_small.c */
#ifdef TFM_SMALL_MONT_SET
/* computes x/R == x (mod N) via Montgomery Reduction */
void fp_montgomery_reduce_small(fp_int *a, fp_int *m, fp_digit mp)
{
   fp_digit c[FP_SIZE], *_c, *tmpm, mu, cy;
   int      oldused, x, y, pa;

#if defined(USE_MEMSET)
   /* now zero the buff */
   memset(c, 0, sizeof c);
#endif
   pa = m->used;

   /* copy the input */
   oldused = a->used;
   for (x = 0; x < oldused; x++) {
       c[x] = a->dp[x];
   }
#if !defined(USE_MEMSET)
   for (; x < 2*pa+3; x++) {
       c[x] = 0;
   }
#endif
   MONT_START;

   switch (pa) {
      case 1:
            x = 0; cy   = 0;
            LOOP_START;
            _c   = c + 0;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
         break;
      case 2:
            x = 0; cy   = 0;
            LOOP_START;
            _c   = c + 0;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 1; cy   = 0;
            LOOP_START;
            _c   = c + 1;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
         break;
      case 3:
            x = 0; cy   = 0;
            LOOP_START;
            _c   = c + 0;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 1; cy   = 0;
            LOOP_START;
            _c   = c + 1;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 2; cy   = 0;
            LOOP_START;
            _c   = c + 2;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
         break;
      case 4:
            x = 0; cy   = 0;
            LOOP_START;
            _c   = c + 0;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 1; cy   = 0;
            LOOP_START;
            _c   = c + 1;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 2; cy   = 0;
            LOOP_START;
            _c   = c + 2;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 3; cy   = 0;
            LOOP_START;
            _c   = c + 3;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
         break;
      case 5:
            x = 0; cy   = 0;
            LOOP_START;
            _c   = c + 0;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 1; cy   = 0;
            LOOP_START;
            _c   = c + 1;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 2; cy   = 0;
            LOOP_START;
            _c   = c + 2;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 3; cy   = 0;
            LOOP_START;
            _c   = c + 3;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 4; cy   = 0;
            LOOP_START;
            _c   = c + 4;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
         break;
      case 6:
            x = 0; cy   = 0;
            LOOP_START;
            _c   = c + 0;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 1; cy   = 0;
            LOOP_START;
            _c   = c + 1;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 2; cy   = 0;
            LOOP_START;
            _c   = c + 2;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 3; cy   = 0;
            LOOP_START;
            _c   = c + 3;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 4; cy   = 0;
            LOOP_START;
            _c   = c + 4;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 5; cy   = 0;
            LOOP_START;
            _c   = c + 5;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
         break;
      case 7:
            x = 0; cy   = 0;
            LOOP_START;
            _c   = c + 0;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 1; cy   = 0;
            LOOP_START;
            _c   = c + 1;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 2; cy   = 0;
            LOOP_START;
            _c   = c + 2;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 3; cy   = 0;
            LOOP_START;
            _c   = c + 3;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 4; cy   = 0;
            LOOP_START;
            _c   = c + 4;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 5; cy   = 0;
            LOOP_START;
            _c   = c + 5;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 6; cy   = 0;
            LOOP_START;
            _c   = c + 6;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
         break;
      case 8:
            x = 0; cy   = 0;
            LOOP_START;
            _c   = c + 0;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 1; cy   = 0;
            LOOP_START;
            _c   = c + 1;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 2; cy   = 0;
            LOOP_START;
            _c   = c + 2;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 3; cy   = 0;
            LOOP_START;
            _c   = c + 3;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 4; cy   = 0;
            LOOP_START;
            _c   = c + 4;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 5; cy   = 0;
            LOOP_START;
            _c   = c + 5;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 6; cy   = 0;
            LOOP_START;
            _c   = c + 6;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 7; cy   = 0;
            LOOP_START;
            _c   = c + 7;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
         break;
      case 9:
            x = 0; cy   = 0;
            LOOP_START;
            _c   = c + 0;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 1; cy   = 0;
            LOOP_START;
            _c   = c + 1;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 2; cy   = 0;
            LOOP_START;
            _c   = c + 2;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 3; cy   = 0;
            LOOP_START;
            _c   = c + 3;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 4; cy   = 0;
            LOOP_START;
            _c   = c + 4;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 5; cy   = 0;
            LOOP_START;
            _c   = c + 5;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 6; cy   = 0;
            LOOP_START;
            _c   = c + 6;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 7; cy   = 0;
            LOOP_START;
            _c   = c + 7;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 8; cy   = 0;
            LOOP_START;
            _c   = c + 8;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
         break;
      case 10:
            x = 0; cy   = 0;
            LOOP_START;
            _c   = c + 0;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 1; cy   = 0;
            LOOP_START;
            _c   = c + 1;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 2; cy   = 0;
            LOOP_START;
            _c   = c + 2;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 3; cy   = 0;
            LOOP_START;
            _c   = c + 3;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 4; cy   = 0;
            LOOP_START;
            _c   = c + 4;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 5; cy   = 0;
            LOOP_START;
            _c   = c + 5;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 6; cy   = 0;
            LOOP_START;
            _c   = c + 6;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 7; cy   = 0;
            LOOP_START;
            _c   = c + 7;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 8; cy   = 0;
            LOOP_START;
            _c   = c + 8;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 9; cy   = 0;
            LOOP_START;
            _c   = c + 9;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
         break;
      case 11:
            x = 0; cy   = 0;
            LOOP_START;
            _c   = c + 0;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 1; cy   = 0;
            LOOP_START;
            _c   = c + 1;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 2; cy   = 0;
            LOOP_START;
            _c   = c + 2;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 3; cy   = 0;
            LOOP_START;
            _c   = c + 3;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 4; cy   = 0;
            LOOP_START;
            _c   = c + 4;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 5; cy   = 0;
            LOOP_START;
            _c   = c + 5;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 6; cy   = 0;
            LOOP_START;
            _c   = c + 6;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 7; cy   = 0;
            LOOP_START;
            _c   = c + 7;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 8; cy   = 0;
            LOOP_START;
            _c   = c + 8;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 9; cy   = 0;
            LOOP_START;
            _c   = c + 9;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 10; cy   = 0;
            LOOP_START;
            _c   = c + 10;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
         break;
      case 12:
            x = 0; cy   = 0;
            LOOP_START;
            _c   = c + 0;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 1; cy   = 0;
            LOOP_START;
            _c   = c + 1;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 2; cy   = 0;
            LOOP_START;
            _c   = c + 2;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 3; cy   = 0;
            LOOP_START;
            _c   = c + 3;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 4; cy   = 0;
            LOOP_START;
            _c   = c + 4;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 5; cy   = 0;
            LOOP_START;
            _c   = c + 5;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 6; cy   = 0;
            LOOP_START;
            _c   = c + 6;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 7; cy   = 0;
            LOOP_START;
            _c   = c + 7;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 8; cy   = 0;
            LOOP_START;
            _c   = c + 8;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 9; cy   = 0;
            LOOP_START;
            _c   = c + 9;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 10; cy   = 0;
            LOOP_START;
            _c   = c + 10;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 11; cy   = 0;
            LOOP_START;
            _c   = c + 11;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
         break;
      case 13:
            x = 0; cy   = 0;
            LOOP_START;
            _c   = c + 0;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 1; cy   = 0;
            LOOP_START;
            _c   = c + 1;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 2; cy   = 0;
            LOOP_START;
            _c   = c + 2;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 3; cy   = 0;
            LOOP_START;
            _c   = c + 3;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 4; cy   = 0;
            LOOP_START;
            _c   = c + 4;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 5; cy   = 0;
            LOOP_START;
            _c   = c + 5;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 6; cy   = 0;
            LOOP_START;
            _c   = c + 6;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 7; cy   = 0;
            LOOP_START;
            _c   = c + 7;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 8; cy   = 0;
            LOOP_START;
            _c   = c + 8;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 9; cy   = 0;
            LOOP_START;
            _c   = c + 9;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 10; cy   = 0;
            LOOP_START;
            _c   = c + 10;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 11; cy   = 0;
            LOOP_START;
            _c   = c + 11;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 12; cy   = 0;
            LOOP_START;
            _c   = c + 12;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
         break;
      case 14:
            x = 0; cy   = 0;
            LOOP_START;
            _c   = c + 0;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 1; cy   = 0;
            LOOP_START;
            _c   = c + 1;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 2; cy   = 0;
            LOOP_START;
            _c   = c + 2;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 3; cy   = 0;
            LOOP_START;
            _c   = c + 3;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 4; cy   = 0;
            LOOP_START;
            _c   = c + 4;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 5; cy   = 0;
            LOOP_START;
            _c   = c + 5;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 6; cy   = 0;
            LOOP_START;
            _c   = c + 6;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 7; cy   = 0;
            LOOP_START;
            _c   = c + 7;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 8; cy   = 0;
            LOOP_START;
            _c   = c + 8;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 9; cy   = 0;
            LOOP_START;
            _c   = c + 9;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 10; cy   = 0;
            LOOP_START;
            _c   = c + 10;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 11; cy   = 0;
            LOOP_START;
            _c   = c + 11;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 12; cy   = 0;
            LOOP_START;
            _c   = c + 12;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 13; cy   = 0;
            LOOP_START;
            _c   = c + 13;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
         break;
      case 15:
            x = 0; cy   = 0;
            LOOP_START;
            _c   = c + 0;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 1; cy   = 0;
            LOOP_START;
            _c   = c + 1;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 2; cy   = 0;
            LOOP_START;
            _c   = c + 2;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 3; cy   = 0;
            LOOP_START;
            _c   = c + 3;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 4; cy   = 0;
            LOOP_START;
            _c   = c + 4;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 5; cy   = 0;
            LOOP_START;
            _c   = c + 5;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 6; cy   = 0;
            LOOP_START;
            _c   = c + 6;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 7; cy   = 0;
            LOOP_START;
            _c   = c + 7;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 8; cy   = 0;
            LOOP_START;
            _c   = c + 8;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 9; cy   = 0;
            LOOP_START;
            _c   = c + 9;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 10; cy   = 0;
            LOOP_START;
            _c   = c + 10;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 11; cy   = 0;
            LOOP_START;
            _c   = c + 11;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 12; cy   = 0;
            LOOP_START;
            _c   = c + 12;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 13; cy   = 0;
            LOOP_START;
            _c   = c + 13;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 14; cy   = 0;
            LOOP_START;
            _c   = c + 14;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
         break;
      case 16:
            x = 0; cy   = 0;
            LOOP_START;
            _c   = c + 0;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 1; cy   = 0;
            LOOP_START;
            _c   = c + 1;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 2; cy   = 0;
            LOOP_START;
            _c   = c + 2;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 3; cy   = 0;
            LOOP_START;
            _c   = c + 3;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 4; cy   = 0;
            LOOP_START;
            _c   = c + 4;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 5; cy   = 0;
            LOOP_START;
            _c   = c + 5;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 6; cy   = 0;
            LOOP_START;
            _c   = c + 6;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 7; cy   = 0;
            LOOP_START;
            _c   = c + 7;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 8; cy   = 0;
            LOOP_START;
            _c   = c + 8;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 9; cy   = 0;
            LOOP_START;
            _c   = c + 9;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 10; cy   = 0;
            LOOP_START;
            _c   = c + 10;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 11; cy   = 0;
            LOOP_START;
            _c   = c + 11;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 12; cy   = 0;
            LOOP_START;
            _c   = c + 12;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 13; cy   = 0;
            LOOP_START;
            _c   = c + 13;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 14; cy   = 0;
            LOOP_START;
            _c   = c + 14;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
            x = 15; cy   = 0;
            LOOP_START;
            _c   = c + 15;
            tmpm = m->dp;
#ifdef INNERMUL8
            INNERMUL8; _c += 8; tmpm += 8;
            INNERMUL8; _c += 8; tmpm += 8;
#else
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
            INNERMUL; ++_c;
#endif
            LOOP_END;
            while (cy) {
               PROPCARRY;
               ++_c;
            }
         break;
  }
  /* now copy out */
  _c   = c + pa;
  tmpm = a->dp;
  for (x = 0; x < pa+1; x++) {
     *tmpm++ = *_c++;
  }

  for (; x < oldused; x++)   {
     *tmpm++ = 0;
  }

  MONT_FINI;

  a->used = pa+1;
  fp_clamp(a);

  /* if A >= m then A = A - m */
  if (fp_cmp_mag (a, m) != FP_LT) {
    s_fp_sub (a, m, a);
  }
}

#endif

/* End: fp_mont_small.c */

/* Start: fp_montgomery_calc_normalization.c */
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

/* computes a = B**n mod b without division or multiplication useful for
 * normalizing numbers in a Montgomery system.
 */
void fp_montgomery_calc_normalization(fp_int *a, fp_int *b)
{
  int     x, bits;

  /* how many bits of last digit does b use */
  bits = fp_count_bits (b) % DIGIT_BIT;
  if (!bits) bits = DIGIT_BIT;

  /* compute A = B^(n-1) * 2^(bits-1) */
  if (b->used > 1) {
     fp_2expt (a, (b->used - 1) * DIGIT_BIT + bits - 1);
  } else {
     fp_set(a, 1);
     bits = 1;
  }

  /* now compute C = A * B mod b */
  for (x = bits - 1; x < (int)DIGIT_BIT; x++) {
    fp_mul_2 (a, a);
    if (fp_cmp_mag (a, b) != FP_LT) {
      s_fp_sub (a, b, a);
    }
  }
}


/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_montgomery_calc_normalization.c */

/* Start: fp_montgomery_reduce.c */
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

/******************************************************************/
#if defined(TFM_X86) && !defined(TFM_SSE2) 
/* x86-32 code */

#define MONT_START 
#define MONT_FINI
#define LOOP_END
#define LOOP_START \
   mu = c[x] * mp

#define INNERMUL                                          \
asm(                                                      \
   "movl %5,%%eax \n\t"                                   \
   "mull %4       \n\t"                                   \
   "addl %1,%%eax \n\t"                                   \
   "adcl $0,%%edx \n\t"                                   \
   "addl %%eax,%0 \n\t"                                   \
   "adcl $0,%%edx \n\t"                                   \
   "movl %%edx,%1 \n\t"                                   \
:"=g"(_c[LO]), "=r"(cy)                                   \
:"0"(_c[LO]), "1"(cy), "g"(mu), "g"(*tmpm++)              \
: "%eax", "%edx", "%cc")

#define PROPCARRY                           \
asm(                                        \
   "addl   %1,%0    \n\t"                   \
   "setb   %%al     \n\t"                   \
   "movzbl %%al,%1 \n\t"                    \
:"=g"(_c[LO]), "=r"(cy)                     \
:"0"(_c[LO]), "1"(cy)                       \
: "%eax", "%cc")

/******************************************************************/
#elif defined(TFM_X86_64)
/* x86-64 code */

#define MONT_START 
#define MONT_FINI
#define LOOP_END
#define LOOP_START \
   mu = c[x] * mp

#define INNERMUL                                          \
asm(                                                      \
   "movq %5,%%rax \n\t"                                   \
   "mulq %4       \n\t"                                   \
   "addq %1,%%rax \n\t"                                   \
   "adcq $0,%%rdx \n\t"                                   \
   "addq %%rax,%0 \n\t"                                   \
   "adcq $0,%%rdx \n\t"                                   \
   "movq %%rdx,%1 \n\t"                                   \
:"=g"(_c[LO]), "=r"(cy)                                   \
:"0"(_c[LO]), "1"(cy), "r"(mu), "r"(*tmpm++)              \
: "%rax", "%rdx", "%cc")

#define INNERMUL8 \
 asm(                  \
 "movq 0(%5),%%rax    \n\t"  \
 "movq 0(%2),%%r10    \n\t"  \
 "movq 0x8(%5),%%r11  \n\t"  \
 "mulq %4             \n\t"  \
 "addq %%r10,%%rax    \n\t"  \
 "adcq $0,%%rdx       \n\t"  \
 "movq 0x8(%2),%%r10  \n\t"  \
 "addq %3,%%rax       \n\t"  \
 "adcq $0,%%rdx       \n\t"  \
 "movq %%rax,0(%0)    \n\t"  \
 "movq %%rdx,%1       \n\t"  \
 \
 "movq %%r11,%%rax    \n\t"  \
 "movq 0x10(%5),%%r11 \n\t"  \
 "mulq %4             \n\t"  \
 "addq %%r10,%%rax    \n\t"  \
 "adcq $0,%%rdx       \n\t"  \
 "movq 0x10(%2),%%r10 \n\t"  \
 "addq %3,%%rax       \n\t"  \
 "adcq $0,%%rdx       \n\t"  \
 "movq %%rax,0x8(%0)  \n\t"  \
 "movq %%rdx,%1       \n\t"  \
 \
 "movq %%r11,%%rax    \n\t"  \
 "movq 0x18(%5),%%r11 \n\t"  \
 "mulq %4             \n\t"  \
 "addq %%r10,%%rax    \n\t"  \
 "adcq $0,%%rdx       \n\t"  \
 "movq 0x18(%2),%%r10 \n\t"  \
 "addq %3,%%rax       \n\t"  \
 "adcq $0,%%rdx       \n\t"  \
 "movq %%rax,0x10(%0) \n\t"  \
 "movq %%rdx,%1       \n\t"  \
 \
 "movq %%r11,%%rax    \n\t"  \
 "movq 0x20(%5),%%r11 \n\t"  \
 "mulq %4             \n\t"  \
 "addq %%r10,%%rax    \n\t"  \
 "adcq $0,%%rdx       \n\t"  \
 "movq 0x20(%2),%%r10 \n\t"  \
 "addq %3,%%rax       \n\t"  \
 "adcq $0,%%rdx       \n\t"  \
 "movq %%rax,0x18(%0) \n\t"  \
 "movq %%rdx,%1       \n\t"  \
 \
 "movq %%r11,%%rax    \n\t"  \
 "movq 0x28(%5),%%r11 \n\t"  \
 "mulq %4             \n\t"  \
 "addq %%r10,%%rax    \n\t"  \
 "adcq $0,%%rdx       \n\t"  \
 "movq 0x28(%2),%%r10 \n\t"  \
 "addq %3,%%rax       \n\t"  \
 "adcq $0,%%rdx       \n\t"  \
 "movq %%rax,0x20(%0) \n\t"  \
 "movq %%rdx,%1       \n\t"  \
 \
 "movq %%r11,%%rax    \n\t"  \
 "movq 0x30(%5),%%r11 \n\t"  \
 "mulq %4             \n\t"  \
 "addq %%r10,%%rax    \n\t"  \
 "adcq $0,%%rdx       \n\t"  \
 "movq 0x30(%2),%%r10 \n\t"  \
 "addq %3,%%rax       \n\t"  \
 "adcq $0,%%rdx       \n\t"  \
 "movq %%rax,0x28(%0) \n\t"  \
 "movq %%rdx,%1       \n\t"  \
 \
 "movq %%r11,%%rax    \n\t"  \
 "movq 0x38(%5),%%r11 \n\t"  \
 "mulq %4             \n\t"  \
 "addq %%r10,%%rax    \n\t"  \
 "adcq $0,%%rdx       \n\t"  \
 "movq 0x38(%2),%%r10 \n\t"  \
 "addq %3,%%rax       \n\t"  \
 "adcq $0,%%rdx       \n\t"  \
 "movq %%rax,0x30(%0) \n\t"  \
 "movq %%rdx,%1       \n\t"  \
 \
 "movq %%r11,%%rax    \n\t"  \
 "mulq %4             \n\t"  \
 "addq %%r10,%%rax    \n\t"  \
 "adcq $0,%%rdx       \n\t"  \
 "addq %3,%%rax       \n\t"  \
 "adcq $0,%%rdx       \n\t"  \
 "movq %%rax,0x38(%0) \n\t"  \
 "movq %%rdx,%1       \n\t"  \
 \
:"=r"(_c), "=r"(cy)                    \
: "0"(_c),  "1"(cy), "g"(mu), "r"(tmpm)\
: "%rax", "%rdx", "%r10", "%r11", "%cc")


#define PROPCARRY                           \
asm(                                        \
   "addq   %1,%0    \n\t"                   \
   "setb   %%al     \n\t"                   \
   "movzbq %%al,%1 \n\t"                    \
:"=g"(_c[LO]), "=r"(cy)                     \
:"0"(_c[LO]), "1"(cy)                       \
: "%rax", "%cc")

/******************************************************************/
#elif defined(TFM_SSE2)  
/* SSE2 code (assumes 32-bit fp_digits) */
/* XMM register assignments:
 * xmm0  *tmpm++, then Mu * (*tmpm++)
 * xmm1  c[x], then Mu
 * xmm2  mp
 * xmm3  cy
 * xmm4  _c[LO]
 */

#define MONT_START \
   asm("movd %0,%%mm2"::"g"(mp))

#define MONT_FINI \
   asm("emms")

#define LOOP_START          \
asm(                        \
"movd %0,%%mm1        \n\t" \
"pxor %%mm3,%%mm3     \n\t" \
"pmuludq %%mm2,%%mm1  \n\t" \
:: "g"(c[x]))

/* pmuludq on mmx registers does a 32x32->64 multiply. */
#define INNERMUL               \
asm(                           \
   "movd %1,%%mm4        \n\t" \
   "movd %2,%%mm0        \n\t" \
   "paddq %%mm4,%%mm3    \n\t" \
   "pmuludq %%mm1,%%mm0  \n\t" \
   "paddq %%mm0,%%mm3    \n\t" \
   "movd %%mm3,%0        \n\t" \
   "psrlq $32, %%mm3     \n\t" \
:"=g"(_c[LO]) : "0"(_c[LO]), "g"(*tmpm++) );

#define INNERMUL8 \
asm(                           \
   "movd 0(%1),%%mm4     \n\t" \
   "movd 0(%2),%%mm0     \n\t" \
   "paddq %%mm4,%%mm3    \n\t" \
   "pmuludq %%mm1,%%mm0  \n\t" \
   "movd 4(%2),%%mm5     \n\t" \
   "paddq %%mm0,%%mm3    \n\t" \
   "movd 4(%1),%%mm6     \n\t" \
   "movd %%mm3,0(%0)     \n\t" \
   "psrlq $32, %%mm3     \n\t" \
\
   "paddq %%mm6,%%mm3    \n\t" \
   "pmuludq %%mm1,%%mm5  \n\t" \
   "movd 8(%2),%%mm6     \n\t" \
   "paddq %%mm5,%%mm3    \n\t" \
   "movd 8(%1),%%mm7     \n\t" \
   "movd %%mm3,4(%0)     \n\t" \
   "psrlq $32, %%mm3     \n\t" \
\
   "paddq %%mm7,%%mm3    \n\t" \
   "pmuludq %%mm1,%%mm6  \n\t" \
   "movd 12(%2),%%mm7    \n\t" \
   "paddq %%mm6,%%mm3    \n\t" \
   "movd 12(%1),%%mm5     \n\t" \
   "movd %%mm3,8(%0)     \n\t" \
   "psrlq $32, %%mm3     \n\t" \
\
   "paddq %%mm5,%%mm3    \n\t" \
   "pmuludq %%mm1,%%mm7  \n\t" \
   "movd 16(%2),%%mm5    \n\t" \
   "paddq %%mm7,%%mm3    \n\t" \
   "movd 16(%1),%%mm6    \n\t" \
   "movd %%mm3,12(%0)    \n\t" \
   "psrlq $32, %%mm3     \n\t" \
\
   "paddq %%mm6,%%mm3    \n\t" \
   "pmuludq %%mm1,%%mm5  \n\t" \
   "movd 20(%2),%%mm6    \n\t" \
   "paddq %%mm5,%%mm3    \n\t" \
   "movd 20(%1),%%mm7    \n\t" \
   "movd %%mm3,16(%0)    \n\t" \
   "psrlq $32, %%mm3     \n\t" \
\
   "paddq %%mm7,%%mm3    \n\t" \
   "pmuludq %%mm1,%%mm6  \n\t" \
   "movd 24(%2),%%mm7    \n\t" \
   "paddq %%mm6,%%mm3    \n\t" \
   "movd 24(%1),%%mm5     \n\t" \
   "movd %%mm3,20(%0)    \n\t" \
   "psrlq $32, %%mm3     \n\t" \
\
   "paddq %%mm5,%%mm3    \n\t" \
   "pmuludq %%mm1,%%mm7  \n\t" \
   "movd 28(%2),%%mm5    \n\t" \
   "paddq %%mm7,%%mm3    \n\t" \
   "movd 28(%1),%%mm6    \n\t" \
   "movd %%mm3,24(%0)    \n\t" \
   "psrlq $32, %%mm3     \n\t" \
\
   "paddq %%mm6,%%mm3    \n\t" \
   "pmuludq %%mm1,%%mm5  \n\t" \
   "paddq %%mm5,%%mm3    \n\t" \
   "movd %%mm3,28(%0)    \n\t" \
   "psrlq $32, %%mm3     \n\t" \
:"=r"(_c) : "0"(_c), "g"(tmpm) );

#define LOOP_END \
asm( "movd %%mm3,%0  \n" :"=r"(cy))

#define PROPCARRY                           \
asm(                                        \
   "addl   %1,%0    \n\t"                   \
   "setb   %%al     \n\t"                   \
   "movzbl %%al,%1 \n\t"                    \
:"=g"(_c[LO]), "=r"(cy)                     \
:"0"(_c[LO]), "1"(cy)                       \
: "%eax", "%cc")

/******************************************************************/
#elif defined(TFM_ARM)
   /* ARMv4 code */

#define MONT_START 
#define MONT_FINI
#define LOOP_END
#define LOOP_START \
   mu = c[x] * mp

#define INNERMUL                    \
asm(                                \
    " LDR    r0,%1            \n\t" \
    " ADDS   r0,r0,%0         \n\t" \
    " MOVCS  %0,#1            \n\t" \
    " MOVCC  %0,#0            \n\t" \
    " UMLAL  r0,%0,%3,%4      \n\t" \
    " STR    r0,%1            \n\t" \
:"=r"(cy),"=m"(_c[0]):"0"(cy),"r"(mu),"r"(*tmpm++),"1"(_c[0]):"r0","%cc");

#define PROPCARRY                  \
asm(                               \
    " LDR   r0,%1            \n\t" \
    " ADDS  r0,r0,%0         \n\t" \
    " STR   r0,%1            \n\t" \
    " MOVCS %0,#1            \n\t" \
    " MOVCC %0,#0            \n\t" \
:"=r"(cy),"=m"(_c[0]):"0"(cy),"1"(_c[0]):"r0","%cc");

#elif defined(TFM_PPC32)

/* PPC32 */
#define MONT_START 
#define MONT_FINI
#define LOOP_END
#define LOOP_START \
   mu = c[x] * mp

#define INNERMUL                     \
asm(                                 \
   " mullw    16,%3,%4       \n\t"   \
   " mulhwu   17,%3,%4       \n\t"   \
   " addc     16,16,%0       \n\t"   \
   " addze    17,17          \n\t"   \
   " lwz      18,%1          \n\t"   \
   " addc     16,16,18       \n\t"   \
   " addze    %0,17          \n\t"   \
   " stw      16,%1          \n\t"   \
:"=r"(cy),"=m"(_c[0]):"0"(cy),"r"(mu),"r"(tmpm[0]),"1"(_c[0]):"16", "17", "18","%cc"); ++tmpm;

#define PROPCARRY                    \
asm(                                 \
   " lwz      16,%1         \n\t"    \
   " addc     16,16,%0      \n\t"    \
   " stw      16,%1         \n\t"    \
   " xor      %0,%0,%0      \n\t"    \
   " addze    %0,%0         \n\t"    \
:"=r"(cy),"=m"(_c[0]):"0"(cy),"1"(_c[0]):"16","%cc");

#elif defined(TFM_PPC64)

/* PPC64 */
#define MONT_START 
#define MONT_FINI
#define LOOP_END
#define LOOP_START \
   mu = c[x] * mp

#define INNERMUL                     \
asm(                                 \
   " mulld    16,%3,%4       \n\t"   \
   " mulhdu   17,%3,%4       \n\t"   \
   " addc     16,16,%0       \n\t"   \
   " addze    17,17          \n\t"   \
   " ldx      18,0,%1        \n\t"   \
   " addc     16,16,18       \n\t"   \
   " addze    %0,17          \n\t"   \
   " sdx      16,0,%1        \n\t"   \
:"=r"(cy),"=m"(_c[0]):"0"(cy),"r"(mu),"r"(tmpm[0]),"1"(_c[0]):"16", "17", "18","%cc"); ++tmpm;

#define PROPCARRY                    \
asm(                                 \
   " ldx      16,0,%1       \n\t"    \
   " addc     16,16,%0      \n\t"    \
   " sdx      16,0,%1       \n\t"    \
   " xor      %0,%0,%0      \n\t"    \
   " addze    %0,%0         \n\t"    \
:"=r"(cy),"=m"(_c[0]):"0"(cy),"1"(_c[0]):"16","%cc");

/******************************************************************/

#elif defined(TFM_AVR32)

/* AVR32 */
#define MONT_START 
#define MONT_FINI
#define LOOP_END
#define LOOP_START \
   mu = c[x] * mp

#define INNERMUL                    \
asm(                                \
    " ld.w   r2,%1            \n\t" \
    " add    r2,%0            \n\t" \
    " eor    r3,r3            \n\t" \
    " acr    r3               \n\t" \
    " macu.d r2,%3,%4         \n\t" \
    " st.w   %1,r2            \n\t" \
    " mov    %0,r3            \n\t" \
:"=r"(cy),"=r"(_c):"0"(cy),"r"(mu),"r"(*tmpm++),"1"(_c):"r2","r3");

#define PROPCARRY                    \
asm(                                 \
   " ld.w     r2,%1         \n\t"    \
   " add      r2,%0         \n\t"    \
   " st.w     %1,r2         \n\t"    \
   " eor      %0,%0         \n\t"    \
   " acr      %0            \n\t"    \
:"=r"(cy),"=r"(&_c[0]):"0"(cy),"1"(&_c[0]):"r2","%cc");

#else

/* ISO C code */
#define MONT_START 
#define MONT_FINI
#define LOOP_END
#define LOOP_START \
   mu = c[x] * mp

#define INNERMUL                                      \
   do { fp_word t;                                    \
   _c[0] = t  = ((fp_word)_c[0] + (fp_word)cy) +      \
                (((fp_word)mu) * ((fp_word)*tmpm++)); \
   cy = (t >> DIGIT_BIT);                             \
   } while (0)

#define PROPCARRY \
   do { fp_digit t = _c[0] += cy; cy = (t < cy); } while (0)

#endif
/******************************************************************/


#define LO  0

#ifdef TFM_SMALL_MONT_SET
#include "fp_mont_small.c"
#endif

/* computes x/R == x (mod N) via Montgomery Reduction */
void fp_montgomery_reduce(fp_int *a, fp_int *m, fp_digit mp)
{
   fp_digit c[FP_SIZE], *_c, *tmpm, mu;
   int      oldused, x, y, pa;

   /* bail if too large */
   if (m->used > (FP_SIZE/2)) {
      return;
   }

#ifdef TFM_SMALL_MONT_SET
   if (m->used <= 16) {
      fp_montgomery_reduce_small(a, m, mp);
      return;
   }
#endif

#if defined(USE_MEMSET)
   /* now zero the buff */
   memset(c, 0, sizeof c);
#endif
   pa = m->used;

   /* copy the input */
   oldused = a->used;
   for (x = 0; x < oldused; x++) {
       c[x] = a->dp[x];
   }
#if !defined(USE_MEMSET)
   for (; x < 2*pa+1; x++) {
       c[x] = 0;
   }
#endif
   MONT_START;

   for (x = 0; x < pa; x++) {
       fp_digit cy = 0;
       /* get Mu for this round */
       LOOP_START;
       _c   = c + x;
       tmpm = m->dp;
       y = 0;
       #if (defined(TFM_SSE2) || defined(TFM_X86_64))
        for (; y < (pa & ~7); y += 8) {
              INNERMUL8;
              _c   += 8;
              tmpm += 8;
           }
       #endif

       for (; y < pa; y++) {
          INNERMUL;
          ++_c;
       }
       LOOP_END;
       while (cy) {
           PROPCARRY;
           ++_c;
       }
  }         

  /* now copy out */
  _c   = c + pa;
  tmpm = a->dp;
  for (x = 0; x < pa+1; x++) {
     *tmpm++ = *_c++;
  }

  for (; x < oldused; x++)   {
     *tmpm++ = 0;
  }

  MONT_FINI;

  a->used = pa+1;
  fp_clamp(a);
  
  /* if A >= m then A = A - m */
  if (fp_cmp_mag (a, m) != FP_LT) {
    s_fp_sub (a, m, a);
  }
}


/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_montgomery_reduce.c */

/* Start: fp_montgomery_setup.c */
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

/* setups the montgomery reduction */
int fp_montgomery_setup(fp_int *a, fp_digit *rho)
{
  fp_digit x, b;

/* fast inversion mod 2**k
 *
 * Based on the fact that
 *
 * XA = 1 (mod 2**n)  =>  (X(2-XA)) A = 1 (mod 2**2n)
 *                    =>  2*X*A - X*X*A*A = 1
 *                    =>  2*(1) - (1)     = 1
 */
  b = a->dp[0];

  if ((b & 1) == 0) {
    return FP_VAL;
  }

  x = (((b + 2) & 4) << 1) + b; /* here x*a==1 mod 2**4 */
  x *= 2 - b * x;               /* here x*a==1 mod 2**8 */
  x *= 2 - b * x;               /* here x*a==1 mod 2**16 */
  x *= 2 - b * x;               /* here x*a==1 mod 2**32 */
#ifdef FP_64BIT
  x *= 2 - b * x;               /* here x*a==1 mod 2**64 */
#endif

  /* rho = -1/m mod b */
  *rho = (((fp_word) 1 << ((fp_word) DIGIT_BIT)) - ((fp_word)x));

  return FP_OKAY;
}


/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_montgomery_setup.c */

/* Start: fp_mul.c */
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

/* c = a * b */
void fp_mul(fp_int *A, fp_int *B, fp_int *C)
{
    int   y, yy;

    /* call generic if we're out of range */
    if (A->used + B->used > FP_SIZE) {
       fp_mul_comba(A, B, C);
       return ;
    }

     y  = MAX(A->used, B->used);
     yy = MIN(A->used, B->used);
    /* pick a comba (unrolled 4/8/16/32 x or rolled) based on the size
       of the largest input.  We also want to avoid doing excess mults if the 
       inputs are not close to the next power of two.  That is, for example,
       if say y=17 then we would do (32-17)^2 = 225 unneeded multiplications 
    */

#ifdef TFM_MUL3
        if (y <= 3) {
           fp_mul_comba3(A,B,C);
           return;
        }
#endif
#ifdef TFM_MUL4
        if (y == 4) {
           fp_mul_comba4(A,B,C);
           return;
        }
#endif
#ifdef TFM_MUL6
        if (y <= 6) {
           fp_mul_comba6(A,B,C);
           return;
        }
#endif
#ifdef TFM_MUL7
        if (y == 7) {
           fp_mul_comba7(A,B,C);
           return;
        }
#endif
#ifdef TFM_MUL8
        if (y == 8) {
           fp_mul_comba8(A,B,C);
           return;
        }
#endif
#ifdef TFM_MUL9
        if (y == 9) {
           fp_mul_comba9(A,B,C);
           return;
        }
#endif
#ifdef TFM_MUL12
        if (y <= 12) {
           fp_mul_comba12(A,B,C);
           return;
        }
#endif
#ifdef TFM_MUL17
        if (y <= 17) {
           fp_mul_comba17(A,B,C);
           return;
        }
#endif

#ifdef TFM_SMALL_SET
        if (y <= 16) {
           fp_mul_comba_small(A,B,C);
           return;
        }
#endif        
#if defined(TFM_MUL20)
        if (y <= 20) {
           fp_mul_comba20(A,B,C);
           return;
        }
#endif
#if defined(TFM_MUL24)
        if (yy >= 16 && y <= 24) {
           fp_mul_comba24(A,B,C);
           return;
        }
#endif
#if defined(TFM_MUL28)
        if (yy >= 20 && y <= 28) {
           fp_mul_comba28(A,B,C);
           return;
        }
#endif
#if defined(TFM_MUL32)
        if (yy >= 24 && y <= 32) {
           fp_mul_comba32(A,B,C);
           return;
        }
#endif
#if defined(TFM_MUL48)
        if (yy >= 40 && y <= 48) {
           fp_mul_comba48(A,B,C);
           return;
        }
#endif        
#if defined(TFM_MUL64)
        if (yy >= 56 && y <= 64) {
           fp_mul_comba64(A,B,C);
           return;
        }
#endif
        fp_mul_comba(A,B,C);
}


/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_mul.c */

/* Start: fp_mul_2.c */
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

void fp_mul_2(fp_int * a, fp_int * b)
{
  int     x, oldused;
   
  oldused = b->used;
  b->used = a->used;

  {
    register fp_digit r, rr, *tmpa, *tmpb;

    /* alias for source */
    tmpa = a->dp;
    
    /* alias for dest */
    tmpb = b->dp;

    /* carry */
    r = 0;
    for (x = 0; x < a->used; x++) {
    
      /* get what will be the *next* carry bit from the 
       * MSB of the current digit 
       */
      rr = *tmpa >> ((fp_digit)(DIGIT_BIT - 1));
      
      /* now shift up this digit, add in the carry [from the previous] */
      *tmpb++ = ((*tmpa++ << ((fp_digit)1)) | r);
      
      /* copy the carry that would be from the source 
       * digit into the next iteration 
       */
      r = rr;
    }

    /* new leading digit? */
    if (r != 0 && b->used != (FP_SIZE-1)) {
      /* add a MSB which is always 1 at this point */
      *tmpb = 1;
      ++(b->used);
    }

    /* now zero any excess digits on the destination 
     * that we didn't write to 
     */
    tmpb = b->dp + b->used;
    for (x = b->used; x < oldused; x++) {
      *tmpb++ = 0;
    }
  }
  b->sign = a->sign;
}


/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_mul_2.c */

/* Start: fp_mul_2d.c */
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

/* c = a * 2**d */
void fp_mul_2d(fp_int *a, int b, fp_int *c)
{
   fp_digit carry, carrytmp, shift;
   int x;

   /* copy it */
   fp_copy(a, c);

   /* handle whole digits */
   if (b >= DIGIT_BIT) {
      fp_lshd(c, b/DIGIT_BIT);
   }
   b %= DIGIT_BIT;

   /* shift the digits */
   if (b != 0) {
      carry = 0;   
      shift = DIGIT_BIT - b;
      for (x = 0; x < c->used; x++) {
          carrytmp = c->dp[x] >> shift;
          c->dp[x] = (c->dp[x] << b) + carry;
          carry = carrytmp;
      }
      /* store last carry if room */
      if (carry && x < FP_SIZE) {
         c->dp[c->used++] = carry;
      }
   }
   fp_clamp(c);
}


/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_mul_2d.c */

/* Start: fp_mul_comba.c */
/* TomsFastMath, a fast ISO C bignum library.
 * 
 * This project is meant to fill in where LibTomMath
 * falls short.  That is speed ;-)
 *
 * This project is public domain and free for all purposes.
 * 
 * Tom St Denis, tomstdenis@gmail.com
 */

/* About this file...

*/

#include <tfm.h>

#if defined(TFM_PRESCOTT) && defined(TFM_SSE2)
   #undef TFM_SSE2
   #define TFM_X86
#endif

/* these are the combas.  Worship them. */
#if defined(TFM_X86)
/* Generic x86 optimized code */

/* anything you need at the start */
#define COMBA_START

/* clear the chaining variables */
#define COMBA_CLEAR \
   c0 = c1 = c2 = 0;

/* forward the carry to the next digit */
#define COMBA_FORWARD \
   do { c0 = c1; c1 = c2; c2 = 0; } while (0);

/* store the first sum */
#define COMBA_STORE(x) \
   x = c0;

/* store the second sum [carry] */
#define COMBA_STORE2(x) \
   x = c1;

/* anything you need at the end */
#define COMBA_FINI

/* this should multiply i and j  */
#define MULADD(i, j)                                      \
asm(                                                      \
     "movl  %6,%%eax     \n\t"                            \
     "mull  %7           \n\t"                            \
     "addl  %%eax,%0     \n\t"                            \
     "adcl  %%edx,%1     \n\t"                            \
     "adcl  $0,%2        \n\t"                            \
     :"=r"(c0), "=r"(c1), "=r"(c2): "0"(c0), "1"(c1), "2"(c2), "m"(i), "m"(j)  :"%eax","%edx","%cc");

#elif defined(TFM_X86_64)
/* x86-64 optimized */

/* anything you need at the start */
#define COMBA_START

/* clear the chaining variables */
#define COMBA_CLEAR \
   c0 = c1 = c2 = 0;

/* forward the carry to the next digit */
#define COMBA_FORWARD \
   do { c0 = c1; c1 = c2; c2 = 0; } while (0);

/* store the first sum */
#define COMBA_STORE(x) \
   x = c0;

/* store the second sum [carry] */
#define COMBA_STORE2(x) \
   x = c1;

/* anything you need at the end */
#define COMBA_FINI

/* this should multiply i and j  */
#define MULADD(i, j)                                      \
asm  (                                                    \
     "movq  %6,%%rax     \n\t"                            \
     "mulq  %7           \n\t"                            \
     "addq  %%rax,%0     \n\t"                            \
     "adcq  %%rdx,%1     \n\t"                            \
     "adcq  $0,%2        \n\t"                            \
     :"=r"(c0), "=r"(c1), "=r"(c2): "0"(c0), "1"(c1), "2"(c2), "g"(i), "g"(j)  :"%rax","%rdx","%cc");

#elif defined(TFM_SSE2)
/* use SSE2 optimizations */

/* anything you need at the start */
#define COMBA_START

/* clear the chaining variables */
#define COMBA_CLEAR \
   c0 = c1 = c2 = 0;

/* forward the carry to the next digit */
#define COMBA_FORWARD \
   do { c0 = c1; c1 = c2; c2 = 0; } while (0);

/* store the first sum */
#define COMBA_STORE(x) \
   x = c0;

/* store the second sum [carry] */
#define COMBA_STORE2(x) \
   x = c1;

/* anything you need at the end */
#define COMBA_FINI \
   asm("emms");

/* this should multiply i and j  */
#define MULADD(i, j)                                     \
asm(                                                     \
    "movd  %6,%%mm0     \n\t"                            \
    "movd  %7,%%mm1     \n\t"                            \
    "pmuludq %%mm1,%%mm0\n\t"                            \
    "movd  %%mm0,%%eax  \n\t"                            \
    "psrlq $32,%%mm0    \n\t"                            \
    "addl  %%eax,%0     \n\t"                            \
    "movd  %%mm0,%%eax  \n\t"                            \
    "adcl  %%eax,%1     \n\t"                            \
    "adcl  $0,%2        \n\t"                            \
    :"=r"(c0), "=r"(c1), "=r"(c2): "0"(c0), "1"(c1), "2"(c2), "m"(i), "m"(j)  :"%eax","%cc");

#elif defined(TFM_ARM)
/* ARM code */

#define COMBA_START 

#define COMBA_CLEAR \
   c0 = c1 = c2 = 0;

#define COMBA_FORWARD \
   do { c0 = c1; c1 = c2; c2 = 0; } while (0);

#define COMBA_STORE(x) \
   x = c0;

#define COMBA_STORE2(x) \
   x = c1;

#define COMBA_FINI

#define MULADD(i, j)                                          \
asm(                                                          \
"  UMULL  r0,r1,%6,%7           \n\t"                         \
"  ADDS   %0,%0,r0              \n\t"                         \
"  ADCS   %1,%1,r1              \n\t"                         \
"  ADC    %2,%2,#0              \n\t"                         \
:"=r"(c0), "=r"(c1), "=r"(c2) : "0"(c0), "1"(c1), "2"(c2), "r"(i), "r"(j) : "r0", "r1", "%cc");

#elif defined(TFM_PPC32)
/* For 32-bit PPC */

#define COMBA_START

#define COMBA_CLEAR \
   c0 = c1 = c2 = 0;

#define COMBA_FORWARD \
   do { c0 = c1; c1 = c2; c2 = 0; } while (0);

#define COMBA_STORE(x) \
   x = c0;

#define COMBA_STORE2(x) \
   x = c1;

#define COMBA_FINI 
   
/* untested: will mulhwu change the flags?  Docs say no */
#define MULADD(i, j)              \
asm(                              \
   " mullw  16,%6,%7       \n\t" \
   " addc   %0,%0,16       \n\t" \
   " mulhwu 16,%6,%7       \n\t" \
   " adde   %1,%1,16       \n\t" \
   " addze  %2,%2          \n\t" \
:"=r"(c0), "=r"(c1), "=r"(c2):"0"(c0), "1"(c1), "2"(c2), "r"(i), "r"(j):"16");

#elif defined(TFM_PPC64)
/* For 64-bit PPC */

#define COMBA_START

#define COMBA_CLEAR \
   c0 = c1 = c2 = 0;

#define COMBA_FORWARD \
   do { c0 = c1; c1 = c2; c2 = 0; } while (0);

#define COMBA_STORE(x) \
   x = c0;

#define COMBA_STORE2(x) \
   x = c1;

#define COMBA_FINI 
   
/* untested: will mulhwu change the flags?  Docs say no */
#define MULADD(i, j)              \
asm(                              \
   " mulld  16,%6,%7       \n\t" \
   " addc   %0,%0,16       \n\t" \
   " mulhdu 16,%6,%7       \n\t" \
   " adde   %1,%1,16       \n\t" \
   " addze  %2,%2          \n\t" \
:"=r"(c0), "=r"(c1), "=r"(c2):"0"(c0), "1"(c1), "2"(c2), "r"(i), "r"(j):"16");

#elif defined(TFM_AVR32)

/* ISO C code */

#define COMBA_START

#define COMBA_CLEAR \
   c0 = c1 = c2 = 0;

#define COMBA_FORWARD \
   do { c0 = c1; c1 = c2; c2 = 0; } while (0);

#define COMBA_STORE(x) \
   x = c0;

#define COMBA_STORE2(x) \
   x = c1;

#define COMBA_FINI 
   
#define MULADD(i, j)             \
asm(                             \
   " mulu.d r2,%6,%7        \n\t"\
   " add    %0,r2           \n\t"\
   " adc    %1,%1,r3        \n\t"\
   " acr    %2              \n\t"\
:"=r"(c0), "=r"(c1), "=r"(c2):"0"(c0), "1"(c1), "2"(c2), "r"(i), "r"(j):"r2","r3");

#else
/* ISO C code */

#define COMBA_START

#define COMBA_CLEAR \
   c0 = c1 = c2 = 0;

#define COMBA_FORWARD \
   do { c0 = c1; c1 = c2; c2 = 0; } while (0);

#define COMBA_STORE(x) \
   x = c0;

#define COMBA_STORE2(x) \
   x = c1;

#define COMBA_FINI 
   
#define MULADD(i, j)                                                              \
   do { fp_word t;                                                                \
   t = (fp_word)c0 + ((fp_word)i) * ((fp_word)j); c0 = t;                         \
   t = (fp_word)c1 + (t >> DIGIT_BIT);            c1 = t; c2 += t >> DIGIT_BIT;   \
   } while (0);

#endif


/* generic PxQ multiplier */
void fp_mul_comba(fp_int *A, fp_int *B, fp_int *C)
{
   int       ix, iy, iz, tx, ty, pa;
   fp_digit  c0, c1, c2, *tmpx, *tmpy;
   fp_int    tmp, *dst;

   COMBA_START;
   COMBA_CLEAR;
   
   /* get size of output and trim */
   pa = A->used + B->used;
   if (pa >= FP_SIZE) {
      pa = FP_SIZE-1;
   }

   if (A == C || B == C) {
      fp_zero(&tmp);
      dst = &tmp;
   } else {
      fp_zero(C);
      dst = C;
   }

   for (ix = 0; ix < pa; ix++) {
      /* get offsets into the two bignums */
      ty = MIN(ix, B->used-1);
      tx = ix - ty;

      /* setup temp aliases */
      tmpx = A->dp + tx;
      tmpy = B->dp + ty;

      /* this is the number of times the loop will iterrate, essentially its 
         while (tx++ < a->used && ty-- >= 0) { ... }
       */
      iy = MIN(A->used-tx, ty+1);

      /* execute loop */
      COMBA_FORWARD;
      for (iz = 0; iz < iy; ++iz) {
          MULADD(*tmpx++, *tmpy--);
      }

      /* store term */
      COMBA_STORE(dst->dp[ix]);
  }
  COMBA_FINI;

  dst->used = pa;
  dst->sign = A->sign ^ B->sign;
  fp_clamp(dst);
  fp_copy(dst, C);
}

#include "fp_mul_comba_small_set.i"
#include "fp_mul_comba_3.i"
#include "fp_mul_comba_4.i"
#include "fp_mul_comba_6.i"
#include "fp_mul_comba_7.i"
#include "fp_mul_comba_8.i"
#include "fp_mul_comba_9.i"
#include "fp_mul_comba_12.i"
#include "fp_mul_comba_17.i"
#include "fp_mul_comba_20.i"
#include "fp_mul_comba_24.i"
#include "fp_mul_comba_28.i"
#include "fp_mul_comba_32.i"
#include "fp_mul_comba_48.i"
#include "fp_mul_comba_64.i"

/* $Source$ */
/* $Revision$ */
/* $Date$ */


/* End: fp_mul_comba.c */

/* Start: fp_mul_d.c */
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

/* c = a * b */
void fp_mul_d(fp_int *a, fp_digit b, fp_int *c)
{
   fp_word  w;
   int      x, oldused;

   oldused = c->used;
   c->used = a->used;
   c->sign = a->sign;
   w       = 0;
   for (x = 0; x < a->used; x++) {
       w         = ((fp_word)a->dp[x]) * ((fp_word)b) + w;
       c->dp[x]  = (fp_digit)w;
       w         = w >> DIGIT_BIT;
   }
   if (w != 0 && (a->used != FP_SIZE)) {
      c->dp[c->used++] = w;
      ++x;
   }
   for (; x < oldused; x++) {
      c->dp[x] = 0;
   }
   fp_clamp(c);
}


/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_mul_d.c */

/* Start: fp_mulmod.c */
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
/* d = a * b (mod c) */
int fp_mulmod(fp_int *a, fp_int *b, fp_int *c, fp_int *d)
{
  fp_int tmp;
  fp_zero(&tmp);
  fp_mul(a, b, &tmp);
  return fp_mod(&tmp, c, d);
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_mulmod.c */

/* Start: fp_prime_miller_rabin.c */
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

/* Miller-Rabin test of "a" to the base of "b" as described in 
 * HAC pp. 139 Algorithm 4.24
 *
 * Sets result to 0 if definitely composite or 1 if probably prime.
 * Randomly the chance of error is no more than 1/4 and often 
 * very much lower.
 */
void fp_prime_miller_rabin (fp_int * a, fp_int * b, int *result)
{
  fp_int  n1, y, r;
  int     s, j;

  /* default */
  *result = FP_NO;

  /* ensure b > 1 */
  if (fp_cmp_d(b, 1) != FP_GT) {
     return;
  }     

  /* get n1 = a - 1 */
  fp_init_copy(&n1, a);
  fp_sub_d(&n1, 1, &n1);

  /* set 2**s * r = n1 */
  fp_init_copy(&r, &n1);

  /* count the number of least significant bits
   * which are zero
   */
  s = fp_cnt_lsb(&r);

  /* now divide n - 1 by 2**s */
  fp_div_2d (&r, s, &r, NULL);

  /* compute y = b**r mod a */
  fp_init(&y);
  fp_exptmod(b, &r, a, &y);

  /* if y != 1 and y != n1 do */
  if (fp_cmp_d (&y, 1) != FP_EQ && fp_cmp (&y, &n1) != FP_EQ) {
    j = 1;
    /* while j <= s-1 and y != n1 */
    while ((j <= (s - 1)) && fp_cmp (&y, &n1) != FP_EQ) {
      fp_sqrmod (&y, a, &y);

      /* if y == 1 then composite */
      if (fp_cmp_d (&y, 1) == FP_EQ) {
         return;
      }
      ++j;
    }

    /* if y != n1 then composite */
    if (fp_cmp (&y, &n1) != FP_EQ) {
       return;
    }
  }

  /* probably prime now */
  *result = FP_YES;
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_prime_miller_rabin.c */

/* Start: fp_prime_random_ex.c */
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

/* This is possibly the mother of all prime generation functions, muahahahahaha! */
int fp_prime_random_ex(fp_int *a, int t, int size, int flags, tfm_prime_callback cb, void *dat)
{
   unsigned char *tmp, maskAND, maskOR_msb, maskOR_lsb;
   int res, err, bsize, maskOR_msb_offset;

   /* sanity check the input */
   if (size <= 1 || t <= 0) {
      return FP_VAL;
   }

   /* TFM_PRIME_SAFE implies TFM_PRIME_BBS */
   if (flags & TFM_PRIME_SAFE) {
      flags |= TFM_PRIME_BBS;
   }

   /* calc the byte size */
   bsize = (size>>3)+(size&7?1:0);

   /* we need a buffer of bsize bytes */
   tmp = malloc(bsize);
   if (tmp == NULL) {
      return FP_MEM;
   }

   /* calc the maskAND value for the MSbyte*/
   maskAND = 0xFF >> (8 - (size & 7));

   /* calc the maskOR_msb */
   maskOR_msb        = 0;
   maskOR_msb_offset = (size - 2) >> 3;
   if (flags & TFM_PRIME_2MSB_ON) {
      maskOR_msb     |= 1 << ((size - 2) & 7);
   } else if (flags & TFM_PRIME_2MSB_OFF) {
      maskAND        &= ~(1 << ((size - 2) & 7));
   }

   /* get the maskOR_lsb */
   maskOR_lsb         = 1;
   if (flags & TFM_PRIME_BBS) {
      maskOR_lsb     |= 3;
   }

   do {
      /* read the bytes */
      if (cb(tmp, bsize, dat) != bsize) {
         err = FP_VAL;
         goto error;
      }
 
      /* work over the MSbyte */
      tmp[0]    &= maskAND;
      tmp[0]    |= 1 << ((size - 1) & 7);

      /* mix in the maskORs */
      tmp[maskOR_msb_offset]   |= maskOR_msb;
      tmp[bsize-1]             |= maskOR_lsb;

      /* read it in */
      fp_read_unsigned_bin(a, tmp, bsize);

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

   err = FP_OKAY;
error:
   free(tmp);
   return err;
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_prime_random_ex.c */

/* Start: fp_radix_size.c */
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

int fp_radix_size(fp_int *a, int radix, int *size)
{
  int     digs;
  fp_int  t;
  fp_digit d;
   
  *size = 0;

  /* check range of the radix */
  if (radix < 2 || radix > 64) {
    return FP_VAL;
  }

  /* quick out if its zero */
  if (fp_iszero(a) == 1) {
     *size = 2;
     return FP_OKAY;
  }

  fp_init_copy(&t, a);

  /* if it is negative output a - */
  if (t.sign == FP_NEG) {
    (*size)++;
    t.sign = FP_ZPOS;
  }

  digs = 0;
  while (fp_iszero (&t) == FP_NO) {
    fp_div_d (&t, (fp_digit) radix, &t, &d);
    (*size)++;
  }

  /* append a NULL so the string is properly terminated */
  (*size)++;
  return FP_OKAY;

}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_radix_size.c */

/* Start: fp_read_radix.c */
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

int fp_read_radix(fp_int *a, char *str, int radix)
{
  int     y, neg;
  char    ch;

  /* make sure the radix is ok */
  if (radix < 2 || radix > 64) {
    return FP_VAL;
  }

  /* if the leading digit is a
   * minus set the sign to negative.
   */
  if (*str == '-') {
    ++str;
    neg = FP_NEG;
  } else {
    neg = FP_ZPOS;
  }

  /* set the integer to the default of zero */
  fp_zero (a);

  /* process each digit of the string */
  while (*str) {
    /* if the radix < 36 the conversion is case insensitive
     * this allows numbers like 1AB and 1ab to represent the same  value
     * [e.g. in hex]
     */
    ch = (char) ((radix < 36) ? toupper (*str) : *str);
    for (y = 0; y < 64; y++) {
      if (ch == fp_s_rmap[y]) {
         break;
      }
    }

    /* if the char was found in the map
     * and is less than the given radix add it
     * to the number, otherwise exit the loop.
     */
    if (y < radix) {
      fp_mul_d (a, (fp_digit) radix, a);
      fp_add_d (a, (fp_digit) y, a);
    } else {
      break;
    }
    ++str;
  }

  /* set the sign only if a != 0 */
  if (fp_iszero(a) != FP_YES) {
     a->sign = neg;
  }
  return FP_OKAY;
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_read_radix.c */

/* Start: fp_read_signed_bin.c */
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

void fp_read_signed_bin(fp_int *a, unsigned char *b, int c)
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

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_read_signed_bin.c */

/* Start: fp_read_unsigned_bin.c */
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

void fp_read_unsigned_bin(fp_int *a, unsigned char *b, int c)
{
  /* zero the int */
  fp_zero (a);

  /* read the bytes in */
  for (; c > 0; c--) {
    fp_mul_2d (a, 8, a);
    a->dp[0] |= *b++;
    a->used += 1;
  }
  fp_clamp (a);
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_read_unsigned_bin.c */

/* Start: fp_reverse.c */
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

/* reverse an array, used for radix code */
void bn_reverse (unsigned char *s, int len)
{
  int     ix, iy;
  unsigned char t;

  ix = 0;
  iy = len - 1;
  while (ix < iy) {
    t     = s[ix];
    s[ix] = s[iy];
    s[iy] = t;
    ++ix;
    --iy;
  }
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_reverse.c */

/* Start: fp_rshd.c */
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

void fp_rshd(fp_int *a, int x)
{
  int y;

  /* too many digits just zero and return */
  if (x >= a->used) {
     fp_zero(a);
     return;
  }

   /* shift */
   for (y = 0; y < a->used - x; y++) {
      a->dp[y] = a->dp[y+x];
   }

   /* zero rest */
   for (; y < a->used; y++) {
      a->dp[y] = 0;
   }
   
   /* decrement count */
   a->used -= x;
   fp_clamp(a);
}


/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_rshd.c */

/* Start: fp_s_rmap.c */
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

/* chars used in radix conversions */
const char *fp_s_rmap = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+/";

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_s_rmap.c */

/* Start: fp_set.c */
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

void fp_set(fp_int *a, fp_digit b)
{
   fp_zero(a);
   a->dp[0] = b;
   a->used  = a->dp[0] ? 1 : 0;
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_set.c */

/* Start: fp_signed_bin_size.c */
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

int fp_signed_bin_size(fp_int *a)
{
  return 1 + fp_unsigned_bin_size (a);
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_signed_bin_size.c */

/* Start: fp_sqr.c */
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

/* b = a*a  */
void fp_sqr(fp_int *A, fp_int *B)
{
    int     y;

    /* call generic if we're out of range */
    if (A->used + A->used > FP_SIZE) {
       fp_sqr_comba(A, B);
       return ;
    }

    y = A->used;
#if defined(TFM_SQR3)
        if (y <= 3) {
           fp_sqr_comba3(A,B);
           return;
        }
#endif
#if defined(TFM_SQR4)
        if (y == 4) {
           fp_sqr_comba4(A,B);
           return;
        }
#endif
#if defined(TFM_SQR6)
        if (y <= 6) {
           fp_sqr_comba6(A,B);
           return;
        }
#endif
#if defined(TFM_SQR7)
        if (y == 7) {
           fp_sqr_comba7(A,B);
           return;
        }
#endif
#if defined(TFM_SQR8)
        if (y == 8) {
           fp_sqr_comba8(A,B);
           return;
        }
#endif
#if defined(TFM_SQR9)
        if (y == 9) {
           fp_sqr_comba9(A,B);
           return;
        }
#endif
#if defined(TFM_SQR12)
        if (y <= 12) {
           fp_sqr_comba12(A,B);
           return;
        }
#endif
#if defined(TFM_SQR17)
        if (y <= 17) {
           fp_sqr_comba17(A,B);
           return;
        }
#endif
#if defined(TFM_SMALL_SET)
        if (y <= 16) {
           fp_sqr_comba_small(A,B);
           return;
        }
#endif
#if defined(TFM_SQR20)
        if (y <= 20) {
           fp_sqr_comba20(A,B);
           return;
        }
#endif
#if defined(TFM_SQR24)
        if (y <= 24) {
           fp_sqr_comba24(A,B);
           return;
        }
#endif
#if defined(TFM_SQR28)
        if (y <= 28) {
           fp_sqr_comba28(A,B);
           return;
        }
#endif
#if defined(TFM_SQR32)
        if (y <= 32) {
           fp_sqr_comba32(A,B);
           return;
        }
#endif
#if defined(TFM_SQR48)
        if (y <= 48) {
           fp_sqr_comba48(A,B);
           return;
        }
#endif
#if defined(TFM_SQR64)
        if (y <= 64) {
           fp_sqr_comba64(A,B);
           return;
        }
#endif
       fp_sqr_comba(A, B);
}


/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_sqr.c */

/* Start: fp_sqr_comba.c */
/*
 * 
 * This project is meant to fill in where LibTomMath
 * falls short.  That is speed ;-)
 *
 * This project is public domain and free for all purposes.
 * 
 * Tom St Denis, tomstdenis@gmail.com
 */
#include <tfm.h>

#if defined(TFM_PRESCOTT) && defined(TFM_SSE2)
   #undef TFM_SSE2
   #define TFM_X86
#endif


#if defined(TFM_X86)

/* x86-32 optimized */

#define COMBA_START

#define CLEAR_CARRY \
   c0 = c1 = c2 = 0;

#define COMBA_STORE(x) \
   x = c0;

#define COMBA_STORE2(x) \
   x = c1;

#define CARRY_FORWARD \
   do { c0 = c1; c1 = c2; c2 = 0; } while (0);

#define COMBA_FINI

#define SQRADD(i, j)                                      \
asm(                                            \
     "movl  %6,%%eax     \n\t"                            \
     "mull  %%eax        \n\t"                            \
     "addl  %%eax,%0     \n\t"                            \
     "adcl  %%edx,%1     \n\t"                            \
     "adcl  $0,%2        \n\t"                            \
     :"=r"(c0), "=r"(c1), "=r"(c2): "0"(c0), "1"(c1), "2"(c2), "m"(i) :"%eax","%edx","%cc");

#define SQRADD2(i, j)                                     \
asm(                                            \
     "movl  %6,%%eax     \n\t"                            \
     "mull  %7           \n\t"                            \
     "addl  %%eax,%0     \n\t"                            \
     "adcl  %%edx,%1     \n\t"                            \
     "adcl  $0,%2        \n\t"                            \
     "addl  %%eax,%0     \n\t"                            \
     "adcl  %%edx,%1     \n\t"                            \
     "adcl  $0,%2        \n\t"                            \
     :"=r"(c0), "=r"(c1), "=r"(c2): "0"(c0), "1"(c1), "2"(c2), "m"(i), "m"(j)  :"%eax","%edx","%cc");

#define SQRADDSC(i, j)                                    \
asm(                                                     \
     "movl  %6,%%eax     \n\t"                            \
     "mull  %7           \n\t"                            \
     "movl  %%eax,%0     \n\t"                            \
     "movl  %%edx,%1     \n\t"                            \
     "xorl  %2,%2        \n\t"                            \
     :"=r"(sc0), "=r"(sc1), "=r"(sc2): "0"(sc0), "1"(sc1), "2"(sc2), "g"(i), "g"(j) :"%eax","%edx","%cc");

#define SQRADDAC(i, j)                                    \
asm(                                                     \
     "movl  %6,%%eax     \n\t"                            \
     "mull  %7           \n\t"                            \
     "addl  %%eax,%0     \n\t"                            \
     "adcl  %%edx,%1     \n\t"                            \
     "adcl  $0,%2        \n\t"                            \
     :"=r"(sc0), "=r"(sc1), "=r"(sc2): "0"(sc0), "1"(sc1), "2"(sc2), "g"(i), "g"(j) :"%eax","%edx","%cc");

#define SQRADDDB                                          \
asm(                                                     \
     "addl %6,%0         \n\t"                            \
     "adcl %7,%1         \n\t"                            \
     "adcl %8,%2         \n\t"                            \
     "addl %6,%0         \n\t"                            \
     "adcl %7,%1         \n\t"                            \
     "adcl %8,%2         \n\t"                            \
     :"=r"(c0), "=r"(c1), "=r"(c2) : "0"(c0), "1"(c1), "2"(c2), "r"(sc0), "r"(sc1), "r"(sc2) : "%cc");

#elif defined(TFM_X86_64)
/* x86-64 optimized */

#define COMBA_START

#define CLEAR_CARRY \
   c0 = c1 = c2 = 0;

#define COMBA_STORE(x) \
   x = c0;

#define COMBA_STORE2(x) \
   x = c1;

#define CARRY_FORWARD \
   do { c0 = c1; c1 = c2; c2 = 0; } while (0);

#define COMBA_FINI

#define SQRADD(i, j)                                      \
asm(                                                     \
     "movq  %6,%%rax     \n\t"                            \
     "mulq  %%rax        \n\t"                            \
     "addq  %%rax,%0     \n\t"                            \
     "adcq  %%rdx,%1     \n\t"                            \
     "adcq  $0,%2        \n\t"                            \
     :"=r"(c0), "=r"(c1), "=r"(c2): "0"(c0), "1"(c1), "2"(c2), "g"(i) :"%rax","%rdx","%cc");

#define SQRADD2(i, j)                                     \
asm(                                                     \
     "movq  %6,%%rax     \n\t"                            \
     "mulq  %7           \n\t"                            \
     "addq  %%rax,%0     \n\t"                            \
     "adcq  %%rdx,%1     \n\t"                            \
     "adcq  $0,%2        \n\t"                            \
     "addq  %%rax,%0     \n\t"                            \
     "adcq  %%rdx,%1     \n\t"                            \
     "adcq  $0,%2        \n\t"                            \
     :"=r"(c0), "=r"(c1), "=r"(c2): "0"(c0), "1"(c1), "2"(c2), "g"(i), "g"(j)  :"%rax","%rdx","%cc");

#define SQRADDSC(i, j)                                    \
asm(                                                     \
     "movq  %6,%%rax     \n\t"                            \
     "mulq  %7           \n\t"                            \
     "movq  %%rax,%0     \n\t"                            \
     "movq  %%rdx,%1     \n\t"                            \
     "xorq  %2,%2        \n\t"                            \
     :"=r"(sc0), "=r"(sc1), "=r"(sc2): "0"(sc0), "1"(sc1), "2"(sc2), "g"(i), "g"(j) :"%rax","%rdx","%cc");

#define SQRADDAC(i, j)                                                         \
asm(                                                     \
     "movq  %6,%%rax     \n\t"                            \
     "mulq  %7           \n\t"                            \
     "addq  %%rax,%0     \n\t"                            \
     "adcq  %%rdx,%1     \n\t"                            \
     "adcq  $0,%2        \n\t"                            \
     :"=r"(sc0), "=r"(sc1), "=r"(sc2): "0"(sc0), "1"(sc1), "2"(sc2), "g"(i), "g"(j) :"%rax","%rdx","%cc");

#define SQRADDDB                                          \
asm(                                                     \
     "addq %6,%0         \n\t"                            \
     "adcq %7,%1         \n\t"                            \
     "adcq %8,%2         \n\t"                            \
     "addq %6,%0         \n\t"                            \
     "adcq %7,%1         \n\t"                            \
     "adcq %8,%2         \n\t"                            \
     :"=r"(c0), "=r"(c1), "=r"(c2) : "0"(c0), "1"(c1), "2"(c2), "r"(sc0), "r"(sc1), "r"(sc2) : "%cc");

#elif defined(TFM_SSE2)

/* SSE2 Optimized */
#define COMBA_START

#define CLEAR_CARRY \
   c0 = c1 = c2 = 0;

#define COMBA_STORE(x) \
   x = c0;

#define COMBA_STORE2(x) \
   x = c1;

#define CARRY_FORWARD \
   do { c0 = c1; c1 = c2; c2 = 0; } while (0);

#define COMBA_FINI \
   asm("emms");

#define SQRADD(i, j)                                      \
asm(                                            \
     "movd  %6,%%mm0     \n\t"                            \
     "pmuludq %%mm0,%%mm0\n\t"                            \
     "movd  %%mm0,%%eax  \n\t"                            \
     "psrlq $32,%%mm0    \n\t"                            \
     "addl  %%eax,%0     \n\t"                            \
     "movd  %%mm0,%%eax  \n\t"                            \
     "adcl  %%eax,%1     \n\t"                            \
     "adcl  $0,%2        \n\t"                            \
     :"=r"(c0), "=r"(c1), "=r"(c2): "0"(c0), "1"(c1), "2"(c2), "m"(i) :"%eax","%cc");

#define SQRADD2(i, j)                                     \
asm(                                            \
     "movd  %6,%%mm0     \n\t"                            \
     "movd  %7,%%mm1     \n\t"                            \
     "pmuludq %%mm1,%%mm0\n\t"                            \
     "movd  %%mm0,%%eax  \n\t"                            \
     "psrlq $32,%%mm0    \n\t"                            \
     "movd  %%mm0,%%edx  \n\t"                            \
     "addl  %%eax,%0     \n\t"                            \
     "adcl  %%edx,%1     \n\t"                            \
     "adcl  $0,%2        \n\t"                            \
     "addl  %%eax,%0     \n\t"                            \
     "adcl  %%edx,%1     \n\t"                            \
     "adcl  $0,%2        \n\t"                            \
     :"=r"(c0), "=r"(c1), "=r"(c2): "0"(c0), "1"(c1), "2"(c2), "m"(i), "m"(j)  :"%eax","%edx","%cc");

#define SQRADDSC(i, j)                                                         \
asm(                                            \
     "movd  %6,%%mm0     \n\t"                            \
     "movd  %7,%%mm1     \n\t"                            \
     "pmuludq %%mm1,%%mm0\n\t"                            \
     "movd  %%mm0,%0     \n\t"                            \
     "psrlq $32,%%mm0    \n\t"                            \
     "movd  %%mm0,%1     \n\t"                            \
     "xorl  %2,%2        \n\t"                            \
     :"=r"(sc0), "=r"(sc1), "=r"(sc2): "0"(sc0), "1"(sc1), "2"(sc2), "m"(i), "m"(j));

#define SQRADDAC(i, j)                                                         \
asm(                                            \
     "movd  %6,%%mm0     \n\t"                            \
     "movd  %7,%%mm1     \n\t"                            \
     "pmuludq %%mm1,%%mm0\n\t"                            \
     "movd  %%mm0,%%eax  \n\t"                            \
     "psrlq $32,%%mm0    \n\t"                            \
     "movd  %%mm0,%%edx  \n\t"                            \
     "addl  %%eax,%0     \n\t"                            \
     "adcl  %%edx,%1     \n\t"                            \
     "adcl  $0,%2        \n\t"                            \
     :"=r"(sc0), "=r"(sc1), "=r"(sc2): "0"(sc0), "1"(sc1), "2"(sc2), "m"(i), "m"(j)  :"%eax","%edx","%cc");

#define SQRADDDB                                          \
asm(                                                     \
     "addl %6,%0         \n\t"                            \
     "adcl %7,%1         \n\t"                            \
     "adcl %8,%2         \n\t"                            \
     "addl %6,%0         \n\t"                            \
     "adcl %7,%1         \n\t"                            \
     "adcl %8,%2         \n\t"                            \
     :"=r"(c0), "=r"(c1), "=r"(c2) : "0"(c0), "1"(c1), "2"(c2), "r"(sc0), "r"(sc1), "r"(sc2) : "%cc");

#elif defined(TFM_ARM)

/* ARM code */

#define COMBA_START

#define CLEAR_CARRY \
   c0 = c1 = c2 = 0;

#define COMBA_STORE(x) \
   x = c0;

#define COMBA_STORE2(x) \
   x = c1;

#define CARRY_FORWARD \
   do { c0 = c1; c1 = c2; c2 = 0; } while (0);

#define COMBA_FINI

/* multiplies point i and j, updates carry "c1" and digit c2 */
#define SQRADD(i, j)                                             \
asm(                                                             \
"  UMULL  r0,r1,%6,%6              \n\t"                         \
"  ADDS   %0,%0,r0                 \n\t"                         \
"  ADCS   %1,%1,r1                 \n\t"                         \
"  ADC    %2,%2,#0                 \n\t"                         \
:"=r"(c0), "=r"(c1), "=r"(c2) : "0"(c0), "1"(c1), "2"(c2), "r"(i) : "r0", "r1", "%cc");
	
/* for squaring some of the terms are doubled... */
#define SQRADD2(i, j)                                            \
asm(                                                             \
"  UMULL  r0,r1,%6,%7              \n\t"                         \
"  ADDS   %0,%0,r0                 \n\t"                         \
"  ADCS   %1,%1,r1                 \n\t"                         \
"  ADC    %2,%2,#0                 \n\t"                         \
"  ADDS   %0,%0,r0                 \n\t"                         \
"  ADCS   %1,%1,r1                 \n\t"                         \
"  ADC    %2,%2,#0                 \n\t"                         \
:"=r"(c0), "=r"(c1), "=r"(c2) : "0"(c0), "1"(c1), "2"(c2), "r"(i), "r"(j) : "r0", "r1", "%cc");

#define SQRADDSC(i, j)                                           \
asm(                                                             \
"  UMULL  %0,%1,%6,%7              \n\t"                         \
"  SUB    %2,%2,%2                 \n\t"                         \
:"=r"(sc0), "=r"(sc1), "=r"(sc2) : "0"(sc0), "1"(sc1), "2"(sc2), "r"(i), "r"(j) : "%cc");

#define SQRADDAC(i, j)                                           \
asm(                                                             \
"  UMULL  r0,r1,%6,%7              \n\t"                         \
"  ADDS   %0,%0,r0                 \n\t"                         \
"  ADCS   %1,%1,r1                 \n\t"                         \
"  ADC    %2,%2,#0                 \n\t"                         \
:"=r"(sc0), "=r"(sc1), "=r"(sc2) : "0"(sc0), "1"(sc1), "2"(sc2), "r"(i), "r"(j) : "r0", "r1", "%cc");

#define SQRADDDB                                                 \
asm(                                                             \
"  ADDS  %0,%0,%3                     \n\t"                      \
"  ADCS  %1,%1,%4                     \n\t"                      \
"  ADC   %2,%2,%5                     \n\t"                      \
"  ADDS  %0,%0,%3                     \n\t"                      \
"  ADCS  %1,%1,%4                     \n\t"                      \
"  ADC   %2,%2,%5                     \n\t"                      \
:"=r"(c0), "=r"(c1), "=r"(c2) : "r"(sc0), "r"(sc1), "r"(sc2), "0"(c0), "1"(c1), "2"(c2) : "%cc");

#elif defined(TFM_PPC32)

/* PPC32 */

#define COMBA_START

#define CLEAR_CARRY \
   c0 = c1 = c2 = 0;

#define COMBA_STORE(x) \
   x = c0;

#define COMBA_STORE2(x) \
   x = c1;

#define CARRY_FORWARD \
   do { c0 = c1; c1 = c2; c2 = 0; } while (0);

#define COMBA_FINI

/* multiplies point i and j, updates carry "c1" and digit c2 */
#define SQRADD(i, j)             \
asm(                             \
   " mullw  16,%6,%6       \n\t" \
   " addc   %0,%0,16       \n\t" \
   " mulhwu 16,%6,%6       \n\t" \
   " adde   %1,%1,16       \n\t" \
   " addze  %2,%2          \n\t" \
:"=r"(c0), "=r"(c1), "=r"(c2):"0"(c0), "1"(c1), "2"(c2), "r"(i):"16","%cc");

/* for squaring some of the terms are doubled... */
#define SQRADD2(i, j)            \
asm(                             \
   " mullw  16,%6,%7       \n\t" \
   " mulhwu 17,%6,%7       \n\t" \
   " addc   %0,%0,16       \n\t" \
   " adde   %1,%1,17       \n\t" \
   " addze  %2,%2          \n\t" \
   " addc   %0,%0,16       \n\t" \
   " adde   %1,%1,17       \n\t" \
   " addze  %2,%2          \n\t" \
:"=r"(c0), "=r"(c1), "=r"(c2):"0"(c0), "1"(c1), "2"(c2), "r"(i), "r"(j):"16", "17","%cc");

#define SQRADDSC(i, j)            \
asm(                              \
   " mullw  %0,%6,%7        \n\t" \
   " mulhwu %1,%6,%7        \n\t" \
   " xor    %2,%2,%2        \n\t" \
:"=r"(sc0), "=r"(sc1), "=r"(sc2):"0"(sc0), "1"(sc1), "2"(sc2), "r"(i),"r"(j) : "%cc");

#define SQRADDAC(i, j)           \
asm(                             \
   " mullw  16,%6,%7       \n\t" \
   " addc   %0,%0,16       \n\t" \
   " mulhwu 16,%6,%7       \n\t" \
   " adde   %1,%1,16       \n\t" \
   " addze  %2,%2          \n\t" \
:"=r"(sc0), "=r"(sc1), "=r"(sc2):"0"(sc0), "1"(sc1), "2"(sc2), "r"(i), "r"(j):"16", "%cc");

#define SQRADDDB                  \
asm(                              \
   " addc   %0,%0,%3        \n\t" \
   " adde   %1,%1,%4        \n\t" \
   " adde   %2,%2,%5        \n\t" \
   " addc   %0,%0,%3        \n\t" \
   " adde   %1,%1,%4        \n\t" \
   " adde   %2,%2,%5        \n\t" \
:"=r"(c0), "=r"(c1), "=r"(c2) : "r"(sc0), "r"(sc1), "r"(sc2), "0"(c0), "1"(c1), "2"(c2) : "%cc");

#elif defined(TFM_PPC64)
/* PPC64 */

#define COMBA_START

#define CLEAR_CARRY \
   c0 = c1 = c2 = 0;

#define COMBA_STORE(x) \
   x = c0;

#define COMBA_STORE2(x) \
   x = c1;

#define CARRY_FORWARD \
   do { c0 = c1; c1 = c2; c2 = 0; } while (0);

#define COMBA_FINI

/* multiplies point i and j, updates carry "c1" and digit c2 */
#define SQRADD(i, j)             \
asm(                             \
   " mulld  16,%6,%6       \n\t" \
   " addc   %0,%0,16       \n\t" \
   " mulhdu 16,%6,%6       \n\t" \
   " adde   %1,%1,16       \n\t" \
   " addze  %2,%2          \n\t" \
:"=r"(c0), "=r"(c1), "=r"(c2):"0"(c0), "1"(c1), "2"(c2), "r"(i):"16","%cc");

/* for squaring some of the terms are doubled... */
#define SQRADD2(i, j)            \
asm(                             \
   " mulld  16,%6,%7       \n\t" \
   " mulhdu 17,%6,%7       \n\t" \
   " addc   %0,%0,16       \n\t" \
   " adde   %1,%1,17       \n\t" \
   " addze  %2,%2          \n\t" \
   " addc   %0,%0,16       \n\t" \
   " adde   %1,%1,17       \n\t" \
   " addze  %2,%2          \n\t" \
:"=r"(c0), "=r"(c1), "=r"(c2):"0"(c0), "1"(c1), "2"(c2), "r"(i), "r"(j):"16", "17","%cc");

#define SQRADDSC(i, j)            \
asm(                              \
   " mulld  %0,%6,%7        \n\t" \
   " mulhdu %1,%6,%7        \n\t" \
   " xor    %2,%2,%2        \n\t" \
:"=r"(sc0), "=r"(sc1), "=r"(sc2):"0"(sc0), "1"(sc1), "2"(sc2), "r"(i),"r"(j) : "%cc");

#define SQRADDAC(i, j)           \
asm(                             \
   " mulld  16,%6,%7       \n\t" \
   " addc   %0,%0,16       \n\t" \
   " mulhdu 16,%6,%7       \n\t" \
   " adde   %1,%1,16       \n\t" \
   " addze  %2,%2          \n\t" \
:"=r"(sc0), "=r"(sc1), "=r"(sc2):"0"(sc0), "1"(sc1), "2"(sc2), "r"(i), "r"(j):"16", "%cc");

#define SQRADDDB                  \
asm(                              \
   " addc   %0,%0,%3        \n\t" \
   " adde   %1,%1,%4        \n\t" \
   " adde   %2,%2,%5        \n\t" \
   " addc   %0,%0,%3        \n\t" \
   " adde   %1,%1,%4        \n\t" \
   " adde   %2,%2,%5        \n\t" \
:"=r"(c0), "=r"(c1), "=r"(c2) : "r"(sc0), "r"(sc1), "r"(sc2), "0"(c0), "1"(c1), "2"(c2) : "%cc");


#elif defined(TFM_AVR32)

/* AVR32 */

#define COMBA_START

#define CLEAR_CARRY \
   c0 = c1 = c2 = 0;

#define COMBA_STORE(x) \
   x = c0;

#define COMBA_STORE2(x) \
   x = c1;

#define CARRY_FORWARD \
   do { c0 = c1; c1 = c2; c2 = 0; } while (0);

#define COMBA_FINI

/* multiplies point i and j, updates carry "c1" and digit c2 */
#define SQRADD(i, j)             \
asm(                             \
   " mulu.d r2,%6,%6       \n\t" \
   " add    %0,%0,r2       \n\t" \
   " adc    %1,%1,r3       \n\t" \
   " acr    %2             \n\t" \
:"=r"(c0), "=r"(c1), "=r"(c2):"0"(c0), "1"(c1), "2"(c2), "r"(i):"r2","r3");

/* for squaring some of the terms are doubled... */
#define SQRADD2(i, j)            \
asm(                             \
   " mulu.d r2,%6,%7       \n\t" \
   " add    %0,%0,r2       \n\t" \
   " adc    %1,%1,r3       \n\t" \
   " acr    %2,            \n\t" \
   " add    %0,%0,r2       \n\t" \
   " adc    %1,%1,r3       \n\t" \
   " acr    %2,            \n\t" \
:"=r"(c0), "=r"(c1), "=r"(c2):"0"(c0), "1"(c1), "2"(c2), "r"(i), "r"(j):"r2", "r3");

#define SQRADDSC(i, j)            \
asm(                              \
   " mulu.d r2,%6,%7        \n\t" \
   " mov    %0,r2           \n\t" \
   " mov    %1,r3           \n\t" \
   " eor    %2,%2           \n\t" \
:"=r"(sc0), "=r"(sc1), "=r"(sc2):"0"(sc0), "1"(sc1), "2"(sc2), "r"(i),"r"(j) : "r2", "r3");

#define SQRADDAC(i, j)           \
asm(                             \
   " mulu.d r2,%6,%7       \n\t" \
   " add    %0,%0,r2       \n\t" \
   " adc    %1,%1,r3       \n\t" \
   " acr    %2             \n\t" \
:"=r"(sc0), "=r"(sc1), "=r"(sc2):"0"(sc0), "1"(sc1), "2"(sc2), "r"(i), "r"(j):"r2", "r3");

#define SQRADDDB                  \
asm(                              \
   " add    %0,%0,%3        \n\t" \
   " adc    %1,%1,%4        \n\t" \
   " adc    %2,%2,%5        \n\t" \
   " add    %0,%0,%3        \n\t" \
   " adc    %1,%1,%4        \n\t" \
   " adc    %2,%2,%5        \n\t" \
:"=r"(c0), "=r"(c1), "=r"(c2) : "r"(sc0), "r"(sc1), "r"(sc2), "0"(c0), "1"(c1), "2"(c2) : "%cc");


#else

#define TFM_ISO

/* ISO C portable code */

#define COMBA_START

#define CLEAR_CARRY \
   c0 = c1 = c2 = 0;

#define COMBA_STORE(x) \
   x = c0;

#define COMBA_STORE2(x) \
   x = c1;

#define CARRY_FORWARD \
   do { c0 = c1; c1 = c2; c2 = 0; } while (0);

#define COMBA_FINI

/* multiplies point i and j, updates carry "c1" and digit c2 */
#define SQRADD(i, j)                                 \
   do { fp_word t;                                   \
   t = c0 + ((fp_word)i) * ((fp_word)j);  c0 = t;    \
   t = c1 + (t >> DIGIT_BIT);             c1 = t; c2 += t >> DIGIT_BIT; \
   } while (0);
  

/* for squaring some of the terms are doubled... */
#define SQRADD2(i, j)                                                 \
   do { fp_word t;                                                    \
   t  = ((fp_word)i) * ((fp_word)j);                                  \
   tt = (fp_word)c0 + t;                 c0 = tt;                              \
   tt = (fp_word)c1 + (tt >> DIGIT_BIT); c1 = tt; c2 += tt >> DIGIT_BIT;       \
   tt = (fp_word)c0 + t;                 c0 = tt;                              \
   tt = (fp_word)c1 + (tt >> DIGIT_BIT); c1 = tt; c2 += tt >> DIGIT_BIT;       \
   } while (0);

#define SQRADDSC(i, j)                                                         \
   do { fp_word t;                                                             \
      t =  ((fp_word)i) * ((fp_word)j);                                        \
      sc0 = (fp_digit)t; sc1 = (t >> DIGIT_BIT); sc2 = 0;                      \
   } while (0);

#define SQRADDAC(i, j)                                                         \
   do { fp_word t;                                                             \
   t = sc0 + ((fp_word)i) * ((fp_word)j);  sc0 = t;                            \
   t = sc1 + (t >> DIGIT_BIT);             sc1 = t; sc2 += t >> DIGIT_BIT;     \
   } while (0);

#define SQRADDDB                                                               \
   do { fp_word t;                                                             \
   t = ((fp_word)sc0) + ((fp_word)sc0) + c0; c0 = t;                                                 \
   t = ((fp_word)sc1) + ((fp_word)sc1) + c1 + (t >> DIGIT_BIT); c1 = t;                              \
   c2 = c2 + ((fp_word)sc2) + ((fp_word)sc2) + (t >> DIGIT_BIT);                                     \
   } while (0);

#endif

#include "fp_sqr_comba_generic.c"
#include "fp_sqr_comba_small_set.i"
#include "fp_sqr_comba_3.i"
#include "fp_sqr_comba_4.i"
#include "fp_sqr_comba_6.i"
#include "fp_sqr_comba_7.i"
#include "fp_sqr_comba_8.i"
#include "fp_sqr_comba_9.i"
#include "fp_sqr_comba_12.i"
#include "fp_sqr_comba_17.i"
#include "fp_sqr_comba_20.i"
#include "fp_sqr_comba_24.i"
#include "fp_sqr_comba_28.i"
#include "fp_sqr_comba_32.i"
#include "fp_sqr_comba_48.i"
#include "fp_sqr_comba_64.i"

/* End: fp_sqr_comba.c */

/* Start: fp_sqr_comba_generic.c */
/* TomsFastMath, a fast ISO C bignum library.
 * 
 * This project is meant to fill in where LibTomMath
 * falls short.  That is speed ;-)
 *
 * This project is public domain and free for all purposes.
 * 
 * Tom St Denis, tomstdenis@gmail.com
 */

/* generic comba squarer */
void fp_sqr_comba(fp_int *A, fp_int *B)
{
  int       pa, ix, iz;
  fp_digit  c0, c1, c2;
  fp_int    tmp, *dst;
#ifdef TFM_ISO
  fp_word   tt;
#endif    

  /* get size of output and trim */
  pa = A->used + A->used;
  if (pa >= FP_SIZE) {
     pa = FP_SIZE-1;
  }

  /* number of output digits to produce */
  COMBA_START;
  CLEAR_CARRY;

  if (A == B) {
     fp_zero(&tmp);
     dst = &tmp;
  } else {
     fp_zero(B);
     dst = B;
  }

  for (ix = 0; ix < pa; ix++) { 
      int      tx, ty, iy;
      fp_digit *tmpy, *tmpx;

      /* get offsets into the two bignums */
      ty = MIN(A->used-1, ix);
      tx = ix - ty;

      /* setup temp aliases */
      tmpx = A->dp + tx;
      tmpy = A->dp + ty;

      /* this is the number of times the loop will iterrate,
         while (tx++ < a->used && ty-- >= 0) { ... }
       */
      iy = MIN(A->used-tx, ty+1);

      /* now for squaring tx can never equal ty 
       * we halve the distance since they approach 
       * at a rate of 2x and we have to round because 
       * odd cases need to be executed
       */
      iy = MIN(iy, (ty-tx+1)>>1);

      /* forward carries */
      CARRY_FORWARD;

      /* execute loop */
      for (iz = 0; iz < iy; iz++) {
          SQRADD2(*tmpx++, *tmpy--);
      }

      /* even columns have the square term in them */
      if ((ix&1) == 0) {
          SQRADD(A->dp[ix>>1], A->dp[ix>>1]);
      }

      /* store it */
      COMBA_STORE(dst->dp[ix]);
  }

  COMBA_FINI;

  /* setup dest */
  dst->used = pa;
  fp_clamp (dst);
  if (dst != B) {
     fp_copy(dst, B);
  }
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_sqr_comba_generic.c */

/* Start: fp_sqrmod.c */
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

/* c = a * a (mod b) */
int fp_sqrmod(fp_int *a, fp_int *b, fp_int *c)
{
  fp_int tmp;
  fp_zero(&tmp);
  fp_sqr(a, &tmp);
  return fp_mod(&tmp, b, c);
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_sqrmod.c */

/* Start: fp_sub.c */
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

/* c = a - b */
void fp_sub(fp_int *a, fp_int *b, fp_int *c)
{
  int     sa, sb;

  sa = a->sign;
  sb = b->sign;

  if (sa != sb) {
    /* subtract a negative from a positive, OR */
    /* subtract a positive from a negative. */
    /* In either case, ADD their magnitudes, */
    /* and use the sign of the first number. */
    c->sign = sa;
    s_fp_add (a, b, c);
  } else {
    /* subtract a positive from a positive, OR */
    /* subtract a negative from a negative. */
    /* First, take the difference between their */
    /* magnitudes, then... */
    if (fp_cmp_mag (a, b) != FP_LT) {
      /* Copy the sign from the first */
      c->sign = sa;
      /* The first has a larger or equal magnitude */
      s_fp_sub (a, b, c);
    } else {
      /* The result has the *opposite* sign from */
      /* the first number. */
      c->sign = (sa == FP_ZPOS) ? FP_NEG : FP_ZPOS;
      /* The second has a larger magnitude */
      s_fp_sub (b, a, c);
    }
  }
}


/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_sub.c */

/* Start: fp_sub_d.c */
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

/* c = a - b */
void fp_sub_d(fp_int *a, fp_digit b, fp_int *c)
{
   fp_int tmp;
   fp_set(&tmp, b);
   fp_sub(a, &tmp, c);
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_sub_d.c */

/* Start: fp_submod.c */
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

/* d = a - b (mod c) */
int fp_submod(fp_int *a, fp_int *b, fp_int *c, fp_int *d)
{
  fp_int tmp;
  fp_zero(&tmp);
  fp_sub(a, b, &tmp);
  return fp_mod(&tmp, c, d);
}


/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_submod.c */

/* Start: fp_to_signed_bin.c */
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

void fp_to_signed_bin(fp_int *a, unsigned char *b)
{
  fp_to_unsigned_bin (a, b + 1);
  b[0] = (unsigned char) ((a->sign == FP_ZPOS) ? 0 : 1);
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_to_signed_bin.c */

/* Start: fp_to_unsigned_bin.c */
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
  bn_reverse (b, x);
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_to_unsigned_bin.c */

/* Start: fp_toradix.c */
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

int fp_toradix(fp_int *a, char *str, int radix)
{
  int     digs;
  fp_int  t;
  fp_digit d;
  char   *_s = str;

  /* check range of the radix */
  if (radix < 2 || radix > 64) {
    return FP_VAL;
  }

  /* quick out if its zero */
  if (fp_iszero(a) == 1) {
     *str++ = '0';
     *str = '\0';
     return FP_OKAY;
  }

  fp_init_copy(&t, a);

  /* if it is negative output a - */
  if (t.sign == FP_NEG) {
    ++_s;
    *str++ = '-';
    t.sign = FP_ZPOS;
  }

  digs = 0;
  while (fp_iszero (&t) == FP_NO) {
    fp_div_d (&t, (fp_digit) radix, &t, &d);
    *str++ = fp_s_rmap[d];
    ++digs;
  }

  /* reverse the digits of the string.  In this case _s points
   * to the first digit [exluding the sign] of the number]
   */
  bn_reverse ((unsigned char *)_s, digs);

  /* append a NULL so the string is properly terminated */
  *str = '\0';
  return FP_OKAY;
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_toradix.c */

/* Start: fp_unsigned_bin_size.c */
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

int fp_unsigned_bin_size(fp_int *a)
{
  int     size = fp_count_bits (a);
  return (size / 8 + ((size & 7) != 0 ? 1 : 0));
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: fp_unsigned_bin_size.c */

/* Start: s_fp_add.c */
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

/* unsigned addition */
void s_fp_add(fp_int *a, fp_int *b, fp_int *c)
{
  int      x, y, oldused;
  register fp_word  t;

  y       = MAX(a->used, b->used);
  oldused = c->used;
  c->used = y;
 
  t = 0;
  for (x = 0; x < y; x++) {
      t         += ((fp_word)a->dp[x]) + ((fp_word)b->dp[x]);
      c->dp[x]   = (fp_digit)t;
      t        >>= DIGIT_BIT;
  }
  if (t != 0 && x < FP_SIZE) {
     c->dp[c->used++] = (fp_digit)t;
     ++x;
  }

  c->used = x;
  for (; x < oldused; x++) {
     c->dp[x] = 0;
  }
  fp_clamp(c);
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: s_fp_add.c */

/* Start: s_fp_sub.c */
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
     t         = (t >> DIGIT_BIT);
   }
  for (; x < oldused; x++) {
     c->dp[x] = 0;
  }
  fp_clamp(c);
}

/* $Source$ */
/* $Revision$ */
/* $Date$ */

/* End: s_fp_sub.c */


/* EOF */
