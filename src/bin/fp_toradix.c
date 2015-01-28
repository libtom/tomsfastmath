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

/**
 * a:		pointer to fp_int representing the input number
 * str:		output buffer
 * radix:	number of character to use for encoding of the number
 * maxlen:	maximum number of the buffer that can be used
 *
 * The radix value can be in the range 2 to 64. This function converts number
 * a into a string str. This function writes at most size bytes (including the
 * terminating null byte to str. It behaves like snprintf(3) in regard to this.
 *
 * Return: If invalid parameter are detected a negative value is returned. On
 * success the function returns the number of bytes that would be written if
 * the function had enough space. Thus a return value of maxlen or more means
 * that the function was not able store all characters and the output is
 * incomplete.
 */
int fp_toradix_n(fp_int *a, char *str, int radix, unsigned int maxlen)
{
  int     digs;
  fp_int  t;
  fp_digit d;
  char   *_s = str;
  unsigned int wrote;

  /* check range of the radix */
  if (radix < 2 || radix > 64)
    return -1;

  /* quick check for zero */
  if (fp_iszero(a) == 1) {
    if (maxlen >= 2)
       *str++ = '0';
    if (maxlen >= 1)
      *str = '\0';
    return 1;
  }

  wrote = 0;

  fp_init_copy(&t, a);

  /* if it is negative output a - */
  if (t.sign == FP_NEG) {
    ++_s;
    wrote++;
    if (wrote < maxlen)
      *str++ = '-';
    t.sign = FP_ZPOS;
  }

  digs = 0;
  while (fp_iszero (&t) == FP_NO) {
    fp_div_d (&t, (fp_digit) radix, &t, &d);
    wrote++;
    if (wrote < maxlen)
      *str++ = fp_s_rmap[d];
    ++digs;
  }

  /* reverse the digits of the string.  In this case _s points
   * to the first digit [exluding the sign] of the number]
   */
  if (wrote < maxlen)
    fp_reverse ((unsigned char *)_s, digs);

    /* append a NULL so the string is properly terminated */
  if (maxlen >= 1)
    *str = '\0';
  return wrote;
}

/**
 * a:		pointer to fp_int representing the input number
 * str:		output buffer
 * radix:	number of character to use for encoding of the number
 *
 * The radix value can be in the range 2 to 64. This function converts number
 * a into a string str. Please don't use this function because a too small
 * chosen str buffer would lead to an overflow which can not be detected.
 * Please use fp_toradix_n() instead.
 *
 * Return: FP_VAL on error, FP_OKAY on success.
 */
int fp_toradix(fp_int *a, char *str, int radix)
{
  int ret;

  ret = fp_toradix_n(a, str, radix, UINT_MAX);
  if (ret < 0)
    return FP_VAL;
  return FP_OKAY;
}
