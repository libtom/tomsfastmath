/* TomsFastMath, a fast ISO C bignum library.
 * 
 * This project is meant to fill in where LibTomMath
 * falls short.  That is speed ;-)
 *
 * This project is public domain and free for all purposes.
 * 
 * Tom St Denis, tomstdenis@iahu.ca
 */
#ifndef TFM_H_
#define TFM_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>

#undef MIN
#define MIN(x,y) ((x)<(y)?(x):(y))
#undef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))

/* do we want large code? */
#define TFM_LARGE

/* do we want huge code (implies large)?  The answer is, yes. */
#define TFM_HUGE

/* imply TFM_LARGE as required */
#if defined(TFM_HUGE)
   #if !defined(TFM_LARGE)
      #define TFM_LARGE
   #endif
#endif

/* Max size of any number in bits.  Basically the largest size you will be multiplying
 * should be half [or smaller] of FP_MAX_SIZE-four_digit
 *
 * You can externally define this or it defaults to 4096-bits.
 */
#ifndef FP_MAX_SIZE
   #define FP_MAX_SIZE           (4096+(4*DIGIT_BIT))
#endif

/* will this lib work? */
#if (CHAR_BIT & 7)
   #error CHAR_BIT must be a multiple of eight.
#endif
#if FP_MAX_SIZE % CHAR_BIT
   #error FP_MAX_SIZE must be a multiple of CHAR_BIT
#endif

/* autodetect x86-64 and make sure we are using 64-bit digits with x86-64 asm */
#if defined(__x86_64__)
   #if defined(TFM_X86) || defined(TFM_SSE2) || defined(TFM_ARM) 
       #error x86-64 detected, x86-32/SSE2/ARM optimizations are not valid!
   #endif
   #if !defined(TFM_X86_64) && !defined(TFM_NO_ASM)
      #define TFM_X86_64
   #endif
#endif
#if defined(TFM_X86_64)
    #if !defined(FP_64BIT)
       #define FP_64BIT
    #endif
#endif

/* try to detect x86-32 */
#if defined(__i386__) && !defined(TFM_SSE2)
   #if defined(TFM_X86_64) || defined(TFM_ARM) 
       #error x86-32 detected, x86-64/ARM optimizations are not valid!
   #endif
   #if !defined(TFM_X86) && !defined(TFM_NO_ASM)
      #define TFM_X86
   #endif
#endif

/* make sure we're 32-bit for x86-32/sse/arm */
#if (defined(TFM_X86) || defined(TFM_SSE2) || defined(TFM_ARM)) && defined(FP_64BIT)
   #warning x86-32, SSE2 and ARM optimizations require 32-bit digits (undefining)
   #undef FP_64BIT
#endif

/* multi asms? */
#ifdef TFM_X86
   #define TFM_ASM
#endif
#ifdef TFM_X86_64
   #ifdef TFM_ASM
      #error TFM_ASM already defined!
   #endif
   #define TFM_ASM
#endif
#ifdef TFM_SSE2
   #ifdef TFM_ASM
      #error TFM_ASM already defined!
   #endif
   #define TFM_ASM
#endif
#ifdef TFM_ARM
   #ifdef TFM_ASM
      #error TFM_ASM already defined!
   #endif
   #define TFM_ASM
#endif

/* we want no asm? */
#ifdef TFM_NO_ASM
   #undef TFM_X86
   #undef TFM_X86_64
   #undef TFM_SSE2
   #undef TFM_ARM
   #undef TFM_ASM   
#endif

/* some default configurations.
 */
#if defined(FP_64BIT)
   /* for GCC only on supported platforms */
#ifndef CRYPT
   typedef unsigned long ulong64;
#endif
   typedef ulong64            fp_digit;
   typedef unsigned long      fp_word __attribute__ ((mode(TI)));
#else
   /* this is to make porting into LibTomCrypt easier :-) */
#ifndef CRYPT
   #if defined(_MSC_VER) || defined(__BORLANDC__) 
      typedef unsigned __int64   ulong64;
      typedef signed __int64     long64;
   #else
      typedef unsigned long long ulong64;
      typedef signed long long   long64;
   #endif
#endif
   typedef unsigned long      fp_digit;
   typedef ulong64            fp_word;
#endif

/* # of digits this is */
#define DIGIT_BIT  (int)((CHAR_BIT) * sizeof(fp_digit))
#define FP_MASK    (fp_digit)(-1)
#define FP_SIZE    (FP_MAX_SIZE/DIGIT_BIT)

/* signs */
#define FP_ZPOS     0
#define FP_NEG      1

/* return codes */
#define FP_OKAY     0
#define FP_VAL      1
#define FP_MEM      2

/* equalities */
#define FP_LT        -1   /* less than */
#define FP_EQ         0   /* equal to */
#define FP_GT         1   /* greater than */

/* replies */
#define FP_YES        1   /* yes response */
#define FP_NO         0   /* no response */

/* a FP type */
typedef struct {
    fp_digit dp[FP_SIZE];
    int      used, 
             sign;
} fp_int;

/* functions */

/* returns a TFM ident string useful for debugging... */
const char *fp_ident(void);

/* initialize [or zero] an fp int */
#define fp_init(a)  (void)memset((a), 0, sizeof(fp_int))
#define fp_zero(a)  fp_init(a)

/* zero/even/odd ? */
#define fp_iszero(a) (((a)->used == 0) ? FP_YES : FP_NO)
#define fp_iseven(a) (((a)->used > 0 && (((a)->dp[0] & 1) == 0)) ? FP_YES : FP_NO)
#define fp_isodd(a)  (((a)->used > 0 && (((a)->dp[0] & 1) == 1)) ? FP_YES : FP_NO)

/* set to a small digit */
void fp_set(fp_int *a, fp_digit b);

/* copy from a to b */
#define fp_copy(a, b)      (void)(((a) != (b)) && memcpy((b), (a), sizeof(fp_int)))
#define fp_init_copy(a, b) fp_copy(b, a)

/* negate and absolute */
#define fp_neg(a, b)  { fp_copy(a, b); (b)->sign ^= 1; }
#define fp_abs(a, b)  { fp_copy(a, b); (b)->sign  = 0; }

/* clamp digits */
#define fp_clamp(a)   { while ((a)->used && (a)->dp[(a)->used-1] == 0) --((a)->used); (a)->sign = (a)->used ? (a)->sign : FP_ZPOS; }

/* right shift x digits */
void fp_rshd(fp_int *a, int x);

/* left shift x digits */
void fp_lshd(fp_int *a, int x);

/* signed comparison */
int fp_cmp(fp_int *a, fp_int *b);

/* unsigned comparison */
int fp_cmp_mag(fp_int *a, fp_int *b);

/* power of 2 operations */
void fp_div_2d(fp_int *a, int b, fp_int *c, fp_int *d);
void fp_mod_2d(fp_int *a, int b, fp_int *c);
void fp_mul_2d(fp_int *a, int b, fp_int *c);
void fp_2expt (fp_int *a, int b);
void fp_mul_2(fp_int *a, fp_int *c);
void fp_div_2(fp_int *a, fp_int *c);

/* Counts the number of lsbs which are zero before the first zero bit */
int fp_cnt_lsb(fp_int *a);

/* c = a + b */
void fp_add(fp_int *a, fp_int *b, fp_int *c);

/* c = a - b */
void fp_sub(fp_int *a, fp_int *b, fp_int *c);

/* c = a * b */
void fp_mul(fp_int *a, fp_int *b, fp_int *c);

/* b = a*a  */
void fp_sqr(fp_int *a, fp_int *b);

/* a/b => cb + d == a */
int fp_div(fp_int *a, fp_int *b, fp_int *c, fp_int *d);

/* c = a mod b, 0 <= c < b  */
int fp_mod(fp_int *a, fp_int *b, fp_int *c);

/* compare against a single digit */
int fp_cmp_d(fp_int *a, fp_digit b);

/* c = a + b */
void fp_add_d(fp_int *a, fp_digit b, fp_int *c);

/* c = a - b */
void fp_sub_d(fp_int *a, fp_digit b, fp_int *c);

/* c = a * b */
void fp_mul_d(fp_int *a, fp_digit b, fp_int *c);

/* a/b => cb + d == a */
int fp_div_d(fp_int *a, fp_digit b, fp_int *c, fp_digit *d);

/* c = a mod b, 0 <= c < b  */
int fp_mod_d(fp_int *a, fp_digit b, fp_digit *c);

/* ---> number theory <--- */
/* d = a + b (mod c) */
int fp_addmod(fp_int *a, fp_int *b, fp_int *c, fp_int *d);

/* d = a - b (mod c) */
int fp_submod(fp_int *a, fp_int *b, fp_int *c, fp_int *d);

/* d = a * b (mod c) */
int fp_mulmod(fp_int *a, fp_int *b, fp_int *c, fp_int *d);

/* c = a * a (mod b) */
int fp_sqrmod(fp_int *a, fp_int *b, fp_int *c);

/* c = 1/a (mod b) */
int fp_invmod(fp_int *a, fp_int *b, fp_int *c);

/* c = (a, b) */
void fp_gcd(fp_int *a, fp_int *b, fp_int *c);

/* c = [a, b] */
void fp_lcm(fp_int *a, fp_int *b, fp_int *c);

/* setups the montgomery reduction */
int fp_montgomery_setup(fp_int *a, fp_digit *mp);

/* computes a = B**n mod b without division or multiplication useful for
 * normalizing numbers in a Montgomery system.
 */
void fp_montgomery_calc_normalization(fp_int *a, fp_int *b);

/* computes x/R == x (mod N) via Montgomery Reduction */
void fp_montgomery_reduce(fp_int *a, fp_int *m, fp_digit mp);

/* d = a**b (mod c) */
int fp_exptmod(fp_int *a, fp_int *b, fp_int *c, fp_int *d);

/* primality stuff */

/* perform a Miller-Rabin test of a to the base b and store result in "result" */
void fp_prime_miller_rabin (fp_int * a, fp_int * b, int *result);

/* 256 trial divisions + 8 Miller-Rabins, returns FP_YES if probable prime  */
int fp_isprime(fp_int *a);

/* Primality generation flags */
#define TFM_PRIME_BBS      0x0001 /* BBS style prime */
#define TFM_PRIME_SAFE     0x0002 /* Safe prime (p-1)/2 == prime */
#define TFM_PRIME_2MSB_OFF 0x0004 /* force 2nd MSB to 0 */
#define TFM_PRIME_2MSB_ON  0x0008 /* force 2nd MSB to 1 */

/* callback for fp_prime_random, should fill dst with random bytes and return how many read [upto len] */
typedef int tfm_prime_callback(unsigned char *dst, int len, void *dat);

#define fp_prime_random(a, t, size, bbs, cb, dat) fp_prime_random_ex(a, t, ((size) * 8) + 1, (bbs==1)?TFM_PRIME_BBS:0, cb, dat)

int fp_prime_random_ex(fp_int *a, int t, int size, int flags, tfm_prime_callback cb, void *dat);

/* radix conersions */
int fp_count_bits(fp_int *a);

int fp_unsigned_bin_size(fp_int *a);
void fp_read_unsigned_bin(fp_int *a, unsigned char *b, int c);
void fp_to_unsigned_bin(fp_int *a, unsigned char *b);

int fp_signed_bin_size(fp_int *a);
void fp_read_signed_bin(fp_int *a, unsigned char *b, int c);
void fp_to_signed_bin(fp_int *a, unsigned char *b);

int fp_read_radix(fp_int *a, char *str, int radix);
int fp_toradix(fp_int *a, char *str, int radix);
int fp_toradix_n(fp_int * a, char *str, int radix, int maxlen);


/* VARIOUS LOW LEVEL STUFFS */
void s_fp_add(fp_int *a, fp_int *b, fp_int *c);
void s_fp_sub(fp_int *a, fp_int *b, fp_int *c);
void bn_reverse(unsigned char *s, int len);
void fp_mul_comba(fp_int *A, fp_int *B, fp_int *C);
#ifdef TFM_HUGE
void fp_mul_comba32(fp_int *A, fp_int *B, fp_int *C);
#endif
#ifdef TFM_LARGE
void fp_mul_comba16(fp_int *A, fp_int *B, fp_int *C);
#endif
void fp_mul_comba8(fp_int *A, fp_int *B, fp_int *C);
void fp_mul_comba4(fp_int *A, fp_int *B, fp_int *C);

void fp_sqr_comba(fp_int *A, fp_int *B);
void fp_sqr_comba4(fp_int *A, fp_int *B);
void fp_sqr_comba8(fp_int *A, fp_int *B);
#ifdef TFM_LARGE
void fp_sqr_comba16(fp_int *A, fp_int *B);
#endif
#ifdef TFM_HUGE
void fp_sqr_comba32(fp_int *A, fp_int *B);
void fp_sqr_comba64(fp_int *A, fp_int *B);
#endif
extern const char *fp_s_rmap;

#endif

