/* TomsFastMath, a fast ISO C bignum library.
 * 
 * This project is meant to fill in where LibTomMath
 * falls short.  That is speed ;-)
 *
 * This project is public domain and free for all purposes.
 * 
 * Tom St Denis, tomstdenis@iahu.ca
 */
#include <tfm.h>

/* c = 1/a (mod b) for odd b only */
int fp_invmod(fp_int *a, fp_int *b, fp_int *c)
{
  fp_int  x, y, u, v, B, D;
  int     neg;

  /* 2. [modified] b must be odd   */
  if (fp_iseven (b) == FP_YES) {
    return FP_VAL;
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
