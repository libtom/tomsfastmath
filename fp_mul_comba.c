/* TomsFastMath, a fast ISO C bignum library.
 * 
 * This project is meant to fill in where LibTomMath
 * falls short.  That is speed ;-)
 *
 * This project is public domain and free for all purposes.
 * 
 * Tom St Denis, tomstdenis@iahu.ca
 */

/* About this file...

*/

#include <tfm.h>

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
asm (                                                     \
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
   #define MULADD(i, j)                                      \
   asm volatile (                                            \
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
"  ADC    %2, %2, #0            \n\t"                         \
:"=r"(c0), "=r"(c1), "=r"(c2) : "0"(c0), "1"(c1), "2"(c2), "r"(i), "r"(j) : "r0", "r1", "%cc");

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
  /* store final carry */
  COMBA_STORE2(dst->dp[ix]);
  COMBA_FINI;

  dst->used = pa;
  fp_clamp(dst);
  dst->sign = dst->used ? A->sign ^ B->sign : FP_ZPOS;
  fp_copy(dst, C);
}

void fp_mul_comba4(fp_int *A, fp_int *B, fp_int *C)
{
   fp_digit c0, c1, c2, at[8];

   memcpy(at, A->dp, 4 * sizeof(fp_digit));
   memcpy(at+4, B->dp, 4 * sizeof(fp_digit));
   COMBA_START;

   COMBA_CLEAR;
   /* 0 */
   MULADD(at[0], at[4]); 
   COMBA_STORE(C->dp[0]);
   /* 1 */
   COMBA_FORWARD;
   MULADD(at[0], at[5]);    MULADD(at[1], at[4]); 
   COMBA_STORE(C->dp[1]);
   /* 2 */
   COMBA_FORWARD;
   MULADD(at[0], at[6]);    MULADD(at[1], at[5]);    MULADD(at[2], at[4]); 
   COMBA_STORE(C->dp[2]);
   /* 3 */
   COMBA_FORWARD;
   MULADD(at[0], at[7]);    MULADD(at[1], at[6]);    MULADD(at[2], at[5]);    MULADD(at[3], at[4]); 
   COMBA_STORE(C->dp[3]);
   /* 4 */
   COMBA_FORWARD;
   MULADD(at[1], at[7]);    MULADD(at[2], at[6]);    MULADD(at[3], at[5]); 
   COMBA_STORE(C->dp[4]);
   /* 5 */
   COMBA_FORWARD;
   MULADD(at[2], at[7]);    MULADD(at[3], at[6]); 
   COMBA_STORE(C->dp[5]);
   /* 6 */
   COMBA_FORWARD;
   MULADD(at[3], at[7]); 
   COMBA_STORE(C->dp[6]);
   COMBA_STORE2(C->dp[7]);
   C->used = 8;
   C->sign = A->sign ^ B->sign;
   fp_clamp(C);
   COMBA_FINI;
}


void fp_mul_comba8(fp_int *A, fp_int *B, fp_int *C)
{
   fp_digit c0, c1, c2, at[16];

   memcpy(at, A->dp, 8 * sizeof(fp_digit));
   memcpy(at+8, B->dp, 8 * sizeof(fp_digit));
   COMBA_START;

   COMBA_CLEAR;
   /* 0 */
   MULADD(at[0], at[8]); 
   COMBA_STORE(C->dp[0]);
   /* 1 */
   COMBA_FORWARD;
   MULADD(at[0], at[9]);    MULADD(at[1], at[8]); 
   COMBA_STORE(C->dp[1]);
   /* 2 */
   COMBA_FORWARD;
   MULADD(at[0], at[10]);    MULADD(at[1], at[9]);    MULADD(at[2], at[8]); 
   COMBA_STORE(C->dp[2]);
   /* 3 */
   COMBA_FORWARD;
   MULADD(at[0], at[11]);    MULADD(at[1], at[10]);    MULADD(at[2], at[9]);    MULADD(at[3], at[8]); 
   COMBA_STORE(C->dp[3]);
   /* 4 */
   COMBA_FORWARD;
   MULADD(at[0], at[12]);    MULADD(at[1], at[11]);    MULADD(at[2], at[10]);    MULADD(at[3], at[9]);    MULADD(at[4], at[8]); 
   COMBA_STORE(C->dp[4]);
   /* 5 */
   COMBA_FORWARD;
   MULADD(at[0], at[13]);    MULADD(at[1], at[12]);    MULADD(at[2], at[11]);    MULADD(at[3], at[10]);    MULADD(at[4], at[9]);    MULADD(at[5], at[8]); 
   COMBA_STORE(C->dp[5]);
   /* 6 */
   COMBA_FORWARD;
   MULADD(at[0], at[14]);    MULADD(at[1], at[13]);    MULADD(at[2], at[12]);    MULADD(at[3], at[11]);    MULADD(at[4], at[10]);    MULADD(at[5], at[9]);    MULADD(at[6], at[8]); 
   COMBA_STORE(C->dp[6]);
   /* 7 */
   COMBA_FORWARD;
   MULADD(at[0], at[15]);    MULADD(at[1], at[14]);    MULADD(at[2], at[13]);    MULADD(at[3], at[12]);    MULADD(at[4], at[11]);    MULADD(at[5], at[10]);    MULADD(at[6], at[9]);    MULADD(at[7], at[8]); 
   COMBA_STORE(C->dp[7]);
   /* 8 */
   COMBA_FORWARD;
   MULADD(at[1], at[15]);    MULADD(at[2], at[14]);    MULADD(at[3], at[13]);    MULADD(at[4], at[12]);    MULADD(at[5], at[11]);    MULADD(at[6], at[10]);    MULADD(at[7], at[9]); 
   COMBA_STORE(C->dp[8]);
   /* 9 */
   COMBA_FORWARD;
   MULADD(at[2], at[15]);    MULADD(at[3], at[14]);    MULADD(at[4], at[13]);    MULADD(at[5], at[12]);    MULADD(at[6], at[11]);    MULADD(at[7], at[10]); 
   COMBA_STORE(C->dp[9]);
   /* 10 */
   COMBA_FORWARD;
   MULADD(at[3], at[15]);    MULADD(at[4], at[14]);    MULADD(at[5], at[13]);    MULADD(at[6], at[12]);    MULADD(at[7], at[11]); 
   COMBA_STORE(C->dp[10]);
   /* 11 */
   COMBA_FORWARD;
   MULADD(at[4], at[15]);    MULADD(at[5], at[14]);    MULADD(at[6], at[13]);    MULADD(at[7], at[12]); 
   COMBA_STORE(C->dp[11]);
   /* 12 */
   COMBA_FORWARD;
   MULADD(at[5], at[15]);    MULADD(at[6], at[14]);    MULADD(at[7], at[13]); 
   COMBA_STORE(C->dp[12]);
   /* 13 */
   COMBA_FORWARD;
   MULADD(at[6], at[15]);    MULADD(at[7], at[14]); 
   COMBA_STORE(C->dp[13]);
   /* 14 */
   COMBA_FORWARD;
   MULADD(at[7], at[15]); 
   COMBA_STORE(C->dp[14]);
   COMBA_STORE2(C->dp[15]);
   C->used = 16;
   C->sign = A->sign ^ B->sign;
   fp_clamp(C);
   COMBA_FINI;
}

#if defined(TFM_LARGE)

void fp_mul_comba16(fp_int *A, fp_int *B, fp_int *C)
{
   fp_digit c0, c1, c2, at[32];

   memcpy(at, A->dp, 16 * sizeof(fp_digit));
   memcpy(at+16, B->dp, 16 * sizeof(fp_digit));
   COMBA_START;

   COMBA_CLEAR;
   /* 0 */
   MULADD(at[0], at[16]); 
   COMBA_STORE(C->dp[0]);
   /* 1 */
   COMBA_FORWARD;
   MULADD(at[0], at[17]);    MULADD(at[1], at[16]); 
   COMBA_STORE(C->dp[1]);
   /* 2 */
   COMBA_FORWARD;
   MULADD(at[0], at[18]);    MULADD(at[1], at[17]);    MULADD(at[2], at[16]); 
   COMBA_STORE(C->dp[2]);
   /* 3 */
   COMBA_FORWARD;
   MULADD(at[0], at[19]);    MULADD(at[1], at[18]);    MULADD(at[2], at[17]);    MULADD(at[3], at[16]); 
   COMBA_STORE(C->dp[3]);
   /* 4 */
   COMBA_FORWARD;
   MULADD(at[0], at[20]);    MULADD(at[1], at[19]);    MULADD(at[2], at[18]);    MULADD(at[3], at[17]);    MULADD(at[4], at[16]); 
   COMBA_STORE(C->dp[4]);
   /* 5 */
   COMBA_FORWARD;
   MULADD(at[0], at[21]);    MULADD(at[1], at[20]);    MULADD(at[2], at[19]);    MULADD(at[3], at[18]);    MULADD(at[4], at[17]);    MULADD(at[5], at[16]); 
   COMBA_STORE(C->dp[5]);
   /* 6 */
   COMBA_FORWARD;
   MULADD(at[0], at[22]);    MULADD(at[1], at[21]);    MULADD(at[2], at[20]);    MULADD(at[3], at[19]);    MULADD(at[4], at[18]);    MULADD(at[5], at[17]);    MULADD(at[6], at[16]); 
   COMBA_STORE(C->dp[6]);
   /* 7 */
   COMBA_FORWARD;
   MULADD(at[0], at[23]);    MULADD(at[1], at[22]);    MULADD(at[2], at[21]);    MULADD(at[3], at[20]);    MULADD(at[4], at[19]);    MULADD(at[5], at[18]);    MULADD(at[6], at[17]);    MULADD(at[7], at[16]); 
   COMBA_STORE(C->dp[7]);
   /* 8 */
   COMBA_FORWARD;
   MULADD(at[0], at[24]);    MULADD(at[1], at[23]);    MULADD(at[2], at[22]);    MULADD(at[3], at[21]);    MULADD(at[4], at[20]);    MULADD(at[5], at[19]);    MULADD(at[6], at[18]);    MULADD(at[7], at[17]);    MULADD(at[8], at[16]); 
   COMBA_STORE(C->dp[8]);
   /* 9 */
   COMBA_FORWARD;
   MULADD(at[0], at[25]);    MULADD(at[1], at[24]);    MULADD(at[2], at[23]);    MULADD(at[3], at[22]);    MULADD(at[4], at[21]);    MULADD(at[5], at[20]);    MULADD(at[6], at[19]);    MULADD(at[7], at[18]);    MULADD(at[8], at[17]);    MULADD(at[9], at[16]); 
   COMBA_STORE(C->dp[9]);
   /* 10 */
   COMBA_FORWARD;
   MULADD(at[0], at[26]);    MULADD(at[1], at[25]);    MULADD(at[2], at[24]);    MULADD(at[3], at[23]);    MULADD(at[4], at[22]);    MULADD(at[5], at[21]);    MULADD(at[6], at[20]);    MULADD(at[7], at[19]);    MULADD(at[8], at[18]);    MULADD(at[9], at[17]);    MULADD(at[10], at[16]); 
   COMBA_STORE(C->dp[10]);
   /* 11 */
   COMBA_FORWARD;
   MULADD(at[0], at[27]);    MULADD(at[1], at[26]);    MULADD(at[2], at[25]);    MULADD(at[3], at[24]);    MULADD(at[4], at[23]);    MULADD(at[5], at[22]);    MULADD(at[6], at[21]);    MULADD(at[7], at[20]);    MULADD(at[8], at[19]);    MULADD(at[9], at[18]);    MULADD(at[10], at[17]);    MULADD(at[11], at[16]); 
   COMBA_STORE(C->dp[11]);
   /* 12 */
   COMBA_FORWARD;
   MULADD(at[0], at[28]);    MULADD(at[1], at[27]);    MULADD(at[2], at[26]);    MULADD(at[3], at[25]);    MULADD(at[4], at[24]);    MULADD(at[5], at[23]);    MULADD(at[6], at[22]);    MULADD(at[7], at[21]);    MULADD(at[8], at[20]);    MULADD(at[9], at[19]);    MULADD(at[10], at[18]);    MULADD(at[11], at[17]);    MULADD(at[12], at[16]); 
   COMBA_STORE(C->dp[12]);
   /* 13 */
   COMBA_FORWARD;
   MULADD(at[0], at[29]);    MULADD(at[1], at[28]);    MULADD(at[2], at[27]);    MULADD(at[3], at[26]);    MULADD(at[4], at[25]);    MULADD(at[5], at[24]);    MULADD(at[6], at[23]);    MULADD(at[7], at[22]);    MULADD(at[8], at[21]);    MULADD(at[9], at[20]);    MULADD(at[10], at[19]);    MULADD(at[11], at[18]);    MULADD(at[12], at[17]);    MULADD(at[13], at[16]); 
   COMBA_STORE(C->dp[13]);
   /* 14 */
   COMBA_FORWARD;
   MULADD(at[0], at[30]);    MULADD(at[1], at[29]);    MULADD(at[2], at[28]);    MULADD(at[3], at[27]);    MULADD(at[4], at[26]);    MULADD(at[5], at[25]);    MULADD(at[6], at[24]);    MULADD(at[7], at[23]);    MULADD(at[8], at[22]);    MULADD(at[9], at[21]);    MULADD(at[10], at[20]);    MULADD(at[11], at[19]);    MULADD(at[12], at[18]);    MULADD(at[13], at[17]);    MULADD(at[14], at[16]); 
   COMBA_STORE(C->dp[14]);
   /* 15 */
   COMBA_FORWARD;
   MULADD(at[0], at[31]);    MULADD(at[1], at[30]);    MULADD(at[2], at[29]);    MULADD(at[3], at[28]);    MULADD(at[4], at[27]);    MULADD(at[5], at[26]);    MULADD(at[6], at[25]);    MULADD(at[7], at[24]);    MULADD(at[8], at[23]);    MULADD(at[9], at[22]);    MULADD(at[10], at[21]);    MULADD(at[11], at[20]);    MULADD(at[12], at[19]);    MULADD(at[13], at[18]);    MULADD(at[14], at[17]);    MULADD(at[15], at[16]); 
   COMBA_STORE(C->dp[15]);
   /* 16 */
   COMBA_FORWARD;
   MULADD(at[1], at[31]);    MULADD(at[2], at[30]);    MULADD(at[3], at[29]);    MULADD(at[4], at[28]);    MULADD(at[5], at[27]);    MULADD(at[6], at[26]);    MULADD(at[7], at[25]);    MULADD(at[8], at[24]);    MULADD(at[9], at[23]);    MULADD(at[10], at[22]);    MULADD(at[11], at[21]);    MULADD(at[12], at[20]);    MULADD(at[13], at[19]);    MULADD(at[14], at[18]);    MULADD(at[15], at[17]); 
   COMBA_STORE(C->dp[16]);
   /* 17 */
   COMBA_FORWARD;
   MULADD(at[2], at[31]);    MULADD(at[3], at[30]);    MULADD(at[4], at[29]);    MULADD(at[5], at[28]);    MULADD(at[6], at[27]);    MULADD(at[7], at[26]);    MULADD(at[8], at[25]);    MULADD(at[9], at[24]);    MULADD(at[10], at[23]);    MULADD(at[11], at[22]);    MULADD(at[12], at[21]);    MULADD(at[13], at[20]);    MULADD(at[14], at[19]);    MULADD(at[15], at[18]); 
   COMBA_STORE(C->dp[17]);
   /* 18 */
   COMBA_FORWARD;
   MULADD(at[3], at[31]);    MULADD(at[4], at[30]);    MULADD(at[5], at[29]);    MULADD(at[6], at[28]);    MULADD(at[7], at[27]);    MULADD(at[8], at[26]);    MULADD(at[9], at[25]);    MULADD(at[10], at[24]);    MULADD(at[11], at[23]);    MULADD(at[12], at[22]);    MULADD(at[13], at[21]);    MULADD(at[14], at[20]);    MULADD(at[15], at[19]); 
   COMBA_STORE(C->dp[18]);
   /* 19 */
   COMBA_FORWARD;
   MULADD(at[4], at[31]);    MULADD(at[5], at[30]);    MULADD(at[6], at[29]);    MULADD(at[7], at[28]);    MULADD(at[8], at[27]);    MULADD(at[9], at[26]);    MULADD(at[10], at[25]);    MULADD(at[11], at[24]);    MULADD(at[12], at[23]);    MULADD(at[13], at[22]);    MULADD(at[14], at[21]);    MULADD(at[15], at[20]); 
   COMBA_STORE(C->dp[19]);
   /* 20 */
   COMBA_FORWARD;
   MULADD(at[5], at[31]);    MULADD(at[6], at[30]);    MULADD(at[7], at[29]);    MULADD(at[8], at[28]);    MULADD(at[9], at[27]);    MULADD(at[10], at[26]);    MULADD(at[11], at[25]);    MULADD(at[12], at[24]);    MULADD(at[13], at[23]);    MULADD(at[14], at[22]);    MULADD(at[15], at[21]); 
   COMBA_STORE(C->dp[20]);
   /* 21 */
   COMBA_FORWARD;
   MULADD(at[6], at[31]);    MULADD(at[7], at[30]);    MULADD(at[8], at[29]);    MULADD(at[9], at[28]);    MULADD(at[10], at[27]);    MULADD(at[11], at[26]);    MULADD(at[12], at[25]);    MULADD(at[13], at[24]);    MULADD(at[14], at[23]);    MULADD(at[15], at[22]); 
   COMBA_STORE(C->dp[21]);
   /* 22 */
   COMBA_FORWARD;
   MULADD(at[7], at[31]);    MULADD(at[8], at[30]);    MULADD(at[9], at[29]);    MULADD(at[10], at[28]);    MULADD(at[11], at[27]);    MULADD(at[12], at[26]);    MULADD(at[13], at[25]);    MULADD(at[14], at[24]);    MULADD(at[15], at[23]); 
   COMBA_STORE(C->dp[22]);
   /* 23 */
   COMBA_FORWARD;
   MULADD(at[8], at[31]);    MULADD(at[9], at[30]);    MULADD(at[10], at[29]);    MULADD(at[11], at[28]);    MULADD(at[12], at[27]);    MULADD(at[13], at[26]);    MULADD(at[14], at[25]);    MULADD(at[15], at[24]); 
   COMBA_STORE(C->dp[23]);
   /* 24 */
   COMBA_FORWARD;
   MULADD(at[9], at[31]);    MULADD(at[10], at[30]);    MULADD(at[11], at[29]);    MULADD(at[12], at[28]);    MULADD(at[13], at[27]);    MULADD(at[14], at[26]);    MULADD(at[15], at[25]); 
   COMBA_STORE(C->dp[24]);
   /* 25 */
   COMBA_FORWARD;
   MULADD(at[10], at[31]);    MULADD(at[11], at[30]);    MULADD(at[12], at[29]);    MULADD(at[13], at[28]);    MULADD(at[14], at[27]);    MULADD(at[15], at[26]); 
   COMBA_STORE(C->dp[25]);
   /* 26 */
   COMBA_FORWARD;
   MULADD(at[11], at[31]);    MULADD(at[12], at[30]);    MULADD(at[13], at[29]);    MULADD(at[14], at[28]);    MULADD(at[15], at[27]); 
   COMBA_STORE(C->dp[26]);
   /* 27 */
   COMBA_FORWARD;
   MULADD(at[12], at[31]);    MULADD(at[13], at[30]);    MULADD(at[14], at[29]);    MULADD(at[15], at[28]); 
   COMBA_STORE(C->dp[27]);
   /* 28 */
   COMBA_FORWARD;
   MULADD(at[13], at[31]);    MULADD(at[14], at[30]);    MULADD(at[15], at[29]); 
   COMBA_STORE(C->dp[28]);
   /* 29 */
   COMBA_FORWARD;
   MULADD(at[14], at[31]);    MULADD(at[15], at[30]); 
   COMBA_STORE(C->dp[29]);
   /* 30 */
   COMBA_FORWARD;
   MULADD(at[15], at[31]); 
   COMBA_STORE(C->dp[30]);
   COMBA_STORE2(C->dp[31]);
   C->used = 32;
   C->sign = A->sign ^ B->sign;
   fp_clamp(C);
   COMBA_FINI;
}

#endif /* TFM_LARGE */

#ifdef TFM_HUGE

void fp_mul_comba32(fp_int *A, fp_int *B, fp_int *C)
{
   fp_digit c0, c1, c2, at[64];

   memcpy(at, A->dp, 32 * sizeof(fp_digit));
   memcpy(at+32, B->dp, 32 * sizeof(fp_digit));
   COMBA_START;

   COMBA_CLEAR;
   /* 0 */
   MULADD(at[0], at[32]); 
   COMBA_STORE(C->dp[0]);
   /* 1 */
   COMBA_FORWARD;
   MULADD(at[0], at[33]);    MULADD(at[1], at[32]); 
   COMBA_STORE(C->dp[1]);
   /* 2 */
   COMBA_FORWARD;
   MULADD(at[0], at[34]);    MULADD(at[1], at[33]);    MULADD(at[2], at[32]); 
   COMBA_STORE(C->dp[2]);
   /* 3 */
   COMBA_FORWARD;
   MULADD(at[0], at[35]);    MULADD(at[1], at[34]);    MULADD(at[2], at[33]);    MULADD(at[3], at[32]); 
   COMBA_STORE(C->dp[3]);
   /* 4 */
   COMBA_FORWARD;
   MULADD(at[0], at[36]);    MULADD(at[1], at[35]);    MULADD(at[2], at[34]);    MULADD(at[3], at[33]);    MULADD(at[4], at[32]); 
   COMBA_STORE(C->dp[4]);
   /* 5 */
   COMBA_FORWARD;
   MULADD(at[0], at[37]);    MULADD(at[1], at[36]);    MULADD(at[2], at[35]);    MULADD(at[3], at[34]);    MULADD(at[4], at[33]);    MULADD(at[5], at[32]); 
   COMBA_STORE(C->dp[5]);
   /* 6 */
   COMBA_FORWARD;
   MULADD(at[0], at[38]);    MULADD(at[1], at[37]);    MULADD(at[2], at[36]);    MULADD(at[3], at[35]);    MULADD(at[4], at[34]);    MULADD(at[5], at[33]);    MULADD(at[6], at[32]); 
   COMBA_STORE(C->dp[6]);
   /* 7 */
   COMBA_FORWARD;
   MULADD(at[0], at[39]);    MULADD(at[1], at[38]);    MULADD(at[2], at[37]);    MULADD(at[3], at[36]);    MULADD(at[4], at[35]);    MULADD(at[5], at[34]);    MULADD(at[6], at[33]);    MULADD(at[7], at[32]); 
   COMBA_STORE(C->dp[7]);
   /* 8 */
   COMBA_FORWARD;
   MULADD(at[0], at[40]);    MULADD(at[1], at[39]);    MULADD(at[2], at[38]);    MULADD(at[3], at[37]);    MULADD(at[4], at[36]);    MULADD(at[5], at[35]);    MULADD(at[6], at[34]);    MULADD(at[7], at[33]);    MULADD(at[8], at[32]); 
   COMBA_STORE(C->dp[8]);
   /* 9 */
   COMBA_FORWARD;
   MULADD(at[0], at[41]);    MULADD(at[1], at[40]);    MULADD(at[2], at[39]);    MULADD(at[3], at[38]);    MULADD(at[4], at[37]);    MULADD(at[5], at[36]);    MULADD(at[6], at[35]);    MULADD(at[7], at[34]);    MULADD(at[8], at[33]);    MULADD(at[9], at[32]); 
   COMBA_STORE(C->dp[9]);
   /* 10 */
   COMBA_FORWARD;
   MULADD(at[0], at[42]);    MULADD(at[1], at[41]);    MULADD(at[2], at[40]);    MULADD(at[3], at[39]);    MULADD(at[4], at[38]);    MULADD(at[5], at[37]);    MULADD(at[6], at[36]);    MULADD(at[7], at[35]);    MULADD(at[8], at[34]);    MULADD(at[9], at[33]);    MULADD(at[10], at[32]); 
   COMBA_STORE(C->dp[10]);
   /* 11 */
   COMBA_FORWARD;
   MULADD(at[0], at[43]);    MULADD(at[1], at[42]);    MULADD(at[2], at[41]);    MULADD(at[3], at[40]);    MULADD(at[4], at[39]);    MULADD(at[5], at[38]);    MULADD(at[6], at[37]);    MULADD(at[7], at[36]);    MULADD(at[8], at[35]);    MULADD(at[9], at[34]);    MULADD(at[10], at[33]);    MULADD(at[11], at[32]); 
   COMBA_STORE(C->dp[11]);
   /* 12 */
   COMBA_FORWARD;
   MULADD(at[0], at[44]);    MULADD(at[1], at[43]);    MULADD(at[2], at[42]);    MULADD(at[3], at[41]);    MULADD(at[4], at[40]);    MULADD(at[5], at[39]);    MULADD(at[6], at[38]);    MULADD(at[7], at[37]);    MULADD(at[8], at[36]);    MULADD(at[9], at[35]);    MULADD(at[10], at[34]);    MULADD(at[11], at[33]);    MULADD(at[12], at[32]); 
   COMBA_STORE(C->dp[12]);
   /* 13 */
   COMBA_FORWARD;
   MULADD(at[0], at[45]);    MULADD(at[1], at[44]);    MULADD(at[2], at[43]);    MULADD(at[3], at[42]);    MULADD(at[4], at[41]);    MULADD(at[5], at[40]);    MULADD(at[6], at[39]);    MULADD(at[7], at[38]);    MULADD(at[8], at[37]);    MULADD(at[9], at[36]);    MULADD(at[10], at[35]);    MULADD(at[11], at[34]);    MULADD(at[12], at[33]);    MULADD(at[13], at[32]); 
   COMBA_STORE(C->dp[13]);
   /* 14 */
   COMBA_FORWARD;
   MULADD(at[0], at[46]);    MULADD(at[1], at[45]);    MULADD(at[2], at[44]);    MULADD(at[3], at[43]);    MULADD(at[4], at[42]);    MULADD(at[5], at[41]);    MULADD(at[6], at[40]);    MULADD(at[7], at[39]);    MULADD(at[8], at[38]);    MULADD(at[9], at[37]);    MULADD(at[10], at[36]);    MULADD(at[11], at[35]);    MULADD(at[12], at[34]);    MULADD(at[13], at[33]);    MULADD(at[14], at[32]); 
   COMBA_STORE(C->dp[14]);
   /* 15 */
   COMBA_FORWARD;
   MULADD(at[0], at[47]);    MULADD(at[1], at[46]);    MULADD(at[2], at[45]);    MULADD(at[3], at[44]);    MULADD(at[4], at[43]);    MULADD(at[5], at[42]);    MULADD(at[6], at[41]);    MULADD(at[7], at[40]);    MULADD(at[8], at[39]);    MULADD(at[9], at[38]);    MULADD(at[10], at[37]);    MULADD(at[11], at[36]);    MULADD(at[12], at[35]);    MULADD(at[13], at[34]);    MULADD(at[14], at[33]);    MULADD(at[15], at[32]); 
   COMBA_STORE(C->dp[15]);
   /* 16 */
   COMBA_FORWARD;
   MULADD(at[0], at[48]);    MULADD(at[1], at[47]);    MULADD(at[2], at[46]);    MULADD(at[3], at[45]);    MULADD(at[4], at[44]);    MULADD(at[5], at[43]);    MULADD(at[6], at[42]);    MULADD(at[7], at[41]);    MULADD(at[8], at[40]);    MULADD(at[9], at[39]);    MULADD(at[10], at[38]);    MULADD(at[11], at[37]);    MULADD(at[12], at[36]);    MULADD(at[13], at[35]);    MULADD(at[14], at[34]);    MULADD(at[15], at[33]);    MULADD(at[16], at[32]); 
   COMBA_STORE(C->dp[16]);
   /* 17 */
   COMBA_FORWARD;
   MULADD(at[0], at[49]);    MULADD(at[1], at[48]);    MULADD(at[2], at[47]);    MULADD(at[3], at[46]);    MULADD(at[4], at[45]);    MULADD(at[5], at[44]);    MULADD(at[6], at[43]);    MULADD(at[7], at[42]);    MULADD(at[8], at[41]);    MULADD(at[9], at[40]);    MULADD(at[10], at[39]);    MULADD(at[11], at[38]);    MULADD(at[12], at[37]);    MULADD(at[13], at[36]);    MULADD(at[14], at[35]);    MULADD(at[15], at[34]);    MULADD(at[16], at[33]);    MULADD(at[17], at[32]); 
   COMBA_STORE(C->dp[17]);
   /* 18 */
   COMBA_FORWARD;
   MULADD(at[0], at[50]);    MULADD(at[1], at[49]);    MULADD(at[2], at[48]);    MULADD(at[3], at[47]);    MULADD(at[4], at[46]);    MULADD(at[5], at[45]);    MULADD(at[6], at[44]);    MULADD(at[7], at[43]);    MULADD(at[8], at[42]);    MULADD(at[9], at[41]);    MULADD(at[10], at[40]);    MULADD(at[11], at[39]);    MULADD(at[12], at[38]);    MULADD(at[13], at[37]);    MULADD(at[14], at[36]);    MULADD(at[15], at[35]);    MULADD(at[16], at[34]);    MULADD(at[17], at[33]);    MULADD(at[18], at[32]); 
   COMBA_STORE(C->dp[18]);
   /* 19 */
   COMBA_FORWARD;
   MULADD(at[0], at[51]);    MULADD(at[1], at[50]);    MULADD(at[2], at[49]);    MULADD(at[3], at[48]);    MULADD(at[4], at[47]);    MULADD(at[5], at[46]);    MULADD(at[6], at[45]);    MULADD(at[7], at[44]);    MULADD(at[8], at[43]);    MULADD(at[9], at[42]);    MULADD(at[10], at[41]);    MULADD(at[11], at[40]);    MULADD(at[12], at[39]);    MULADD(at[13], at[38]);    MULADD(at[14], at[37]);    MULADD(at[15], at[36]);    MULADD(at[16], at[35]);    MULADD(at[17], at[34]);    MULADD(at[18], at[33]);    MULADD(at[19], at[32]); 
   COMBA_STORE(C->dp[19]);
   /* 20 */
   COMBA_FORWARD;
   MULADD(at[0], at[52]);    MULADD(at[1], at[51]);    MULADD(at[2], at[50]);    MULADD(at[3], at[49]);    MULADD(at[4], at[48]);    MULADD(at[5], at[47]);    MULADD(at[6], at[46]);    MULADD(at[7], at[45]);    MULADD(at[8], at[44]);    MULADD(at[9], at[43]);    MULADD(at[10], at[42]);    MULADD(at[11], at[41]);    MULADD(at[12], at[40]);    MULADD(at[13], at[39]);    MULADD(at[14], at[38]);    MULADD(at[15], at[37]);    MULADD(at[16], at[36]);    MULADD(at[17], at[35]);    MULADD(at[18], at[34]);    MULADD(at[19], at[33]);    MULADD(at[20], at[32]); 
   COMBA_STORE(C->dp[20]);
   /* 21 */
   COMBA_FORWARD;
   MULADD(at[0], at[53]);    MULADD(at[1], at[52]);    MULADD(at[2], at[51]);    MULADD(at[3], at[50]);    MULADD(at[4], at[49]);    MULADD(at[5], at[48]);    MULADD(at[6], at[47]);    MULADD(at[7], at[46]);    MULADD(at[8], at[45]);    MULADD(at[9], at[44]);    MULADD(at[10], at[43]);    MULADD(at[11], at[42]);    MULADD(at[12], at[41]);    MULADD(at[13], at[40]);    MULADD(at[14], at[39]);    MULADD(at[15], at[38]);    MULADD(at[16], at[37]);    MULADD(at[17], at[36]);    MULADD(at[18], at[35]);    MULADD(at[19], at[34]);    MULADD(at[20], at[33]);    MULADD(at[21], at[32]); 
   COMBA_STORE(C->dp[21]);
   /* 22 */
   COMBA_FORWARD;
   MULADD(at[0], at[54]);    MULADD(at[1], at[53]);    MULADD(at[2], at[52]);    MULADD(at[3], at[51]);    MULADD(at[4], at[50]);    MULADD(at[5], at[49]);    MULADD(at[6], at[48]);    MULADD(at[7], at[47]);    MULADD(at[8], at[46]);    MULADD(at[9], at[45]);    MULADD(at[10], at[44]);    MULADD(at[11], at[43]);    MULADD(at[12], at[42]);    MULADD(at[13], at[41]);    MULADD(at[14], at[40]);    MULADD(at[15], at[39]);    MULADD(at[16], at[38]);    MULADD(at[17], at[37]);    MULADD(at[18], at[36]);    MULADD(at[19], at[35]);    MULADD(at[20], at[34]);    MULADD(at[21], at[33]);    MULADD(at[22], at[32]); 
   COMBA_STORE(C->dp[22]);
   /* 23 */
   COMBA_FORWARD;
   MULADD(at[0], at[55]);    MULADD(at[1], at[54]);    MULADD(at[2], at[53]);    MULADD(at[3], at[52]);    MULADD(at[4], at[51]);    MULADD(at[5], at[50]);    MULADD(at[6], at[49]);    MULADD(at[7], at[48]);    MULADD(at[8], at[47]);    MULADD(at[9], at[46]);    MULADD(at[10], at[45]);    MULADD(at[11], at[44]);    MULADD(at[12], at[43]);    MULADD(at[13], at[42]);    MULADD(at[14], at[41]);    MULADD(at[15], at[40]);    MULADD(at[16], at[39]);    MULADD(at[17], at[38]);    MULADD(at[18], at[37]);    MULADD(at[19], at[36]);    MULADD(at[20], at[35]);    MULADD(at[21], at[34]);    MULADD(at[22], at[33]);    MULADD(at[23], at[32]); 
   COMBA_STORE(C->dp[23]);
   /* 24 */
   COMBA_FORWARD;
   MULADD(at[0], at[56]);    MULADD(at[1], at[55]);    MULADD(at[2], at[54]);    MULADD(at[3], at[53]);    MULADD(at[4], at[52]);    MULADD(at[5], at[51]);    MULADD(at[6], at[50]);    MULADD(at[7], at[49]);    MULADD(at[8], at[48]);    MULADD(at[9], at[47]);    MULADD(at[10], at[46]);    MULADD(at[11], at[45]);    MULADD(at[12], at[44]);    MULADD(at[13], at[43]);    MULADD(at[14], at[42]);    MULADD(at[15], at[41]);    MULADD(at[16], at[40]);    MULADD(at[17], at[39]);    MULADD(at[18], at[38]);    MULADD(at[19], at[37]);    MULADD(at[20], at[36]);    MULADD(at[21], at[35]);    MULADD(at[22], at[34]);    MULADD(at[23], at[33]);    MULADD(at[24], at[32]); 
   COMBA_STORE(C->dp[24]);
   /* 25 */
   COMBA_FORWARD;
   MULADD(at[0], at[57]);    MULADD(at[1], at[56]);    MULADD(at[2], at[55]);    MULADD(at[3], at[54]);    MULADD(at[4], at[53]);    MULADD(at[5], at[52]);    MULADD(at[6], at[51]);    MULADD(at[7], at[50]);    MULADD(at[8], at[49]);    MULADD(at[9], at[48]);    MULADD(at[10], at[47]);    MULADD(at[11], at[46]);    MULADD(at[12], at[45]);    MULADD(at[13], at[44]);    MULADD(at[14], at[43]);    MULADD(at[15], at[42]);    MULADD(at[16], at[41]);    MULADD(at[17], at[40]);    MULADD(at[18], at[39]);    MULADD(at[19], at[38]);    MULADD(at[20], at[37]);    MULADD(at[21], at[36]);    MULADD(at[22], at[35]);    MULADD(at[23], at[34]);    MULADD(at[24], at[33]);    MULADD(at[25], at[32]); 
   COMBA_STORE(C->dp[25]);
   /* 26 */
   COMBA_FORWARD;
   MULADD(at[0], at[58]);    MULADD(at[1], at[57]);    MULADD(at[2], at[56]);    MULADD(at[3], at[55]);    MULADD(at[4], at[54]);    MULADD(at[5], at[53]);    MULADD(at[6], at[52]);    MULADD(at[7], at[51]);    MULADD(at[8], at[50]);    MULADD(at[9], at[49]);    MULADD(at[10], at[48]);    MULADD(at[11], at[47]);    MULADD(at[12], at[46]);    MULADD(at[13], at[45]);    MULADD(at[14], at[44]);    MULADD(at[15], at[43]);    MULADD(at[16], at[42]);    MULADD(at[17], at[41]);    MULADD(at[18], at[40]);    MULADD(at[19], at[39]);    MULADD(at[20], at[38]);    MULADD(at[21], at[37]);    MULADD(at[22], at[36]);    MULADD(at[23], at[35]);    MULADD(at[24], at[34]);    MULADD(at[25], at[33]);    MULADD(at[26], at[32]); 
   COMBA_STORE(C->dp[26]);
   /* 27 */
   COMBA_FORWARD;
   MULADD(at[0], at[59]);    MULADD(at[1], at[58]);    MULADD(at[2], at[57]);    MULADD(at[3], at[56]);    MULADD(at[4], at[55]);    MULADD(at[5], at[54]);    MULADD(at[6], at[53]);    MULADD(at[7], at[52]);    MULADD(at[8], at[51]);    MULADD(at[9], at[50]);    MULADD(at[10], at[49]);    MULADD(at[11], at[48]);    MULADD(at[12], at[47]);    MULADD(at[13], at[46]);    MULADD(at[14], at[45]);    MULADD(at[15], at[44]);    MULADD(at[16], at[43]);    MULADD(at[17], at[42]);    MULADD(at[18], at[41]);    MULADD(at[19], at[40]);    MULADD(at[20], at[39]);    MULADD(at[21], at[38]);    MULADD(at[22], at[37]);    MULADD(at[23], at[36]);    MULADD(at[24], at[35]);    MULADD(at[25], at[34]);    MULADD(at[26], at[33]);    MULADD(at[27], at[32]); 
   COMBA_STORE(C->dp[27]);
   /* 28 */
   COMBA_FORWARD;
   MULADD(at[0], at[60]);    MULADD(at[1], at[59]);    MULADD(at[2], at[58]);    MULADD(at[3], at[57]);    MULADD(at[4], at[56]);    MULADD(at[5], at[55]);    MULADD(at[6], at[54]);    MULADD(at[7], at[53]);    MULADD(at[8], at[52]);    MULADD(at[9], at[51]);    MULADD(at[10], at[50]);    MULADD(at[11], at[49]);    MULADD(at[12], at[48]);    MULADD(at[13], at[47]);    MULADD(at[14], at[46]);    MULADD(at[15], at[45]);    MULADD(at[16], at[44]);    MULADD(at[17], at[43]);    MULADD(at[18], at[42]);    MULADD(at[19], at[41]);    MULADD(at[20], at[40]);    MULADD(at[21], at[39]);    MULADD(at[22], at[38]);    MULADD(at[23], at[37]);    MULADD(at[24], at[36]);    MULADD(at[25], at[35]);    MULADD(at[26], at[34]);    MULADD(at[27], at[33]);    MULADD(at[28], at[32]); 
   COMBA_STORE(C->dp[28]);
   /* 29 */
   COMBA_FORWARD;
   MULADD(at[0], at[61]);    MULADD(at[1], at[60]);    MULADD(at[2], at[59]);    MULADD(at[3], at[58]);    MULADD(at[4], at[57]);    MULADD(at[5], at[56]);    MULADD(at[6], at[55]);    MULADD(at[7], at[54]);    MULADD(at[8], at[53]);    MULADD(at[9], at[52]);    MULADD(at[10], at[51]);    MULADD(at[11], at[50]);    MULADD(at[12], at[49]);    MULADD(at[13], at[48]);    MULADD(at[14], at[47]);    MULADD(at[15], at[46]);    MULADD(at[16], at[45]);    MULADD(at[17], at[44]);    MULADD(at[18], at[43]);    MULADD(at[19], at[42]);    MULADD(at[20], at[41]);    MULADD(at[21], at[40]);    MULADD(at[22], at[39]);    MULADD(at[23], at[38]);    MULADD(at[24], at[37]);    MULADD(at[25], at[36]);    MULADD(at[26], at[35]);    MULADD(at[27], at[34]);    MULADD(at[28], at[33]);    MULADD(at[29], at[32]); 
   COMBA_STORE(C->dp[29]);
   /* 30 */
   COMBA_FORWARD;
   MULADD(at[0], at[62]);    MULADD(at[1], at[61]);    MULADD(at[2], at[60]);    MULADD(at[3], at[59]);    MULADD(at[4], at[58]);    MULADD(at[5], at[57]);    MULADD(at[6], at[56]);    MULADD(at[7], at[55]);    MULADD(at[8], at[54]);    MULADD(at[9], at[53]);    MULADD(at[10], at[52]);    MULADD(at[11], at[51]);    MULADD(at[12], at[50]);    MULADD(at[13], at[49]);    MULADD(at[14], at[48]);    MULADD(at[15], at[47]);    MULADD(at[16], at[46]);    MULADD(at[17], at[45]);    MULADD(at[18], at[44]);    MULADD(at[19], at[43]);    MULADD(at[20], at[42]);    MULADD(at[21], at[41]);    MULADD(at[22], at[40]);    MULADD(at[23], at[39]);    MULADD(at[24], at[38]);    MULADD(at[25], at[37]);    MULADD(at[26], at[36]);    MULADD(at[27], at[35]);    MULADD(at[28], at[34]);    MULADD(at[29], at[33]);    MULADD(at[30], at[32]); 
   COMBA_STORE(C->dp[30]);
   /* 31 */
   COMBA_FORWARD;
   MULADD(at[0], at[63]);    MULADD(at[1], at[62]);    MULADD(at[2], at[61]);    MULADD(at[3], at[60]);    MULADD(at[4], at[59]);    MULADD(at[5], at[58]);    MULADD(at[6], at[57]);    MULADD(at[7], at[56]);    MULADD(at[8], at[55]);    MULADD(at[9], at[54]);    MULADD(at[10], at[53]);    MULADD(at[11], at[52]);    MULADD(at[12], at[51]);    MULADD(at[13], at[50]);    MULADD(at[14], at[49]);    MULADD(at[15], at[48]);    MULADD(at[16], at[47]);    MULADD(at[17], at[46]);    MULADD(at[18], at[45]);    MULADD(at[19], at[44]);    MULADD(at[20], at[43]);    MULADD(at[21], at[42]);    MULADD(at[22], at[41]);    MULADD(at[23], at[40]);    MULADD(at[24], at[39]);    MULADD(at[25], at[38]);    MULADD(at[26], at[37]);    MULADD(at[27], at[36]);    MULADD(at[28], at[35]);    MULADD(at[29], at[34]);    MULADD(at[30], at[33]);    MULADD(at[31], at[32]); 
   COMBA_STORE(C->dp[31]);
   /* 32 */
   COMBA_FORWARD;
   MULADD(at[1], at[63]);    MULADD(at[2], at[62]);    MULADD(at[3], at[61]);    MULADD(at[4], at[60]);    MULADD(at[5], at[59]);    MULADD(at[6], at[58]);    MULADD(at[7], at[57]);    MULADD(at[8], at[56]);    MULADD(at[9], at[55]);    MULADD(at[10], at[54]);    MULADD(at[11], at[53]);    MULADD(at[12], at[52]);    MULADD(at[13], at[51]);    MULADD(at[14], at[50]);    MULADD(at[15], at[49]);    MULADD(at[16], at[48]);    MULADD(at[17], at[47]);    MULADD(at[18], at[46]);    MULADD(at[19], at[45]);    MULADD(at[20], at[44]);    MULADD(at[21], at[43]);    MULADD(at[22], at[42]);    MULADD(at[23], at[41]);    MULADD(at[24], at[40]);    MULADD(at[25], at[39]);    MULADD(at[26], at[38]);    MULADD(at[27], at[37]);    MULADD(at[28], at[36]);    MULADD(at[29], at[35]);    MULADD(at[30], at[34]);    MULADD(at[31], at[33]); 
   COMBA_STORE(C->dp[32]);
   /* 33 */
   COMBA_FORWARD;
   MULADD(at[2], at[63]);    MULADD(at[3], at[62]);    MULADD(at[4], at[61]);    MULADD(at[5], at[60]);    MULADD(at[6], at[59]);    MULADD(at[7], at[58]);    MULADD(at[8], at[57]);    MULADD(at[9], at[56]);    MULADD(at[10], at[55]);    MULADD(at[11], at[54]);    MULADD(at[12], at[53]);    MULADD(at[13], at[52]);    MULADD(at[14], at[51]);    MULADD(at[15], at[50]);    MULADD(at[16], at[49]);    MULADD(at[17], at[48]);    MULADD(at[18], at[47]);    MULADD(at[19], at[46]);    MULADD(at[20], at[45]);    MULADD(at[21], at[44]);    MULADD(at[22], at[43]);    MULADD(at[23], at[42]);    MULADD(at[24], at[41]);    MULADD(at[25], at[40]);    MULADD(at[26], at[39]);    MULADD(at[27], at[38]);    MULADD(at[28], at[37]);    MULADD(at[29], at[36]);    MULADD(at[30], at[35]);    MULADD(at[31], at[34]); 
   COMBA_STORE(C->dp[33]);
   /* 34 */
   COMBA_FORWARD;
   MULADD(at[3], at[63]);    MULADD(at[4], at[62]);    MULADD(at[5], at[61]);    MULADD(at[6], at[60]);    MULADD(at[7], at[59]);    MULADD(at[8], at[58]);    MULADD(at[9], at[57]);    MULADD(at[10], at[56]);    MULADD(at[11], at[55]);    MULADD(at[12], at[54]);    MULADD(at[13], at[53]);    MULADD(at[14], at[52]);    MULADD(at[15], at[51]);    MULADD(at[16], at[50]);    MULADD(at[17], at[49]);    MULADD(at[18], at[48]);    MULADD(at[19], at[47]);    MULADD(at[20], at[46]);    MULADD(at[21], at[45]);    MULADD(at[22], at[44]);    MULADD(at[23], at[43]);    MULADD(at[24], at[42]);    MULADD(at[25], at[41]);    MULADD(at[26], at[40]);    MULADD(at[27], at[39]);    MULADD(at[28], at[38]);    MULADD(at[29], at[37]);    MULADD(at[30], at[36]);    MULADD(at[31], at[35]); 
   COMBA_STORE(C->dp[34]);
   /* 35 */
   COMBA_FORWARD;
   MULADD(at[4], at[63]);    MULADD(at[5], at[62]);    MULADD(at[6], at[61]);    MULADD(at[7], at[60]);    MULADD(at[8], at[59]);    MULADD(at[9], at[58]);    MULADD(at[10], at[57]);    MULADD(at[11], at[56]);    MULADD(at[12], at[55]);    MULADD(at[13], at[54]);    MULADD(at[14], at[53]);    MULADD(at[15], at[52]);    MULADD(at[16], at[51]);    MULADD(at[17], at[50]);    MULADD(at[18], at[49]);    MULADD(at[19], at[48]);    MULADD(at[20], at[47]);    MULADD(at[21], at[46]);    MULADD(at[22], at[45]);    MULADD(at[23], at[44]);    MULADD(at[24], at[43]);    MULADD(at[25], at[42]);    MULADD(at[26], at[41]);    MULADD(at[27], at[40]);    MULADD(at[28], at[39]);    MULADD(at[29], at[38]);    MULADD(at[30], at[37]);    MULADD(at[31], at[36]); 
   COMBA_STORE(C->dp[35]);
   /* 36 */
   COMBA_FORWARD;
   MULADD(at[5], at[63]);    MULADD(at[6], at[62]);    MULADD(at[7], at[61]);    MULADD(at[8], at[60]);    MULADD(at[9], at[59]);    MULADD(at[10], at[58]);    MULADD(at[11], at[57]);    MULADD(at[12], at[56]);    MULADD(at[13], at[55]);    MULADD(at[14], at[54]);    MULADD(at[15], at[53]);    MULADD(at[16], at[52]);    MULADD(at[17], at[51]);    MULADD(at[18], at[50]);    MULADD(at[19], at[49]);    MULADD(at[20], at[48]);    MULADD(at[21], at[47]);    MULADD(at[22], at[46]);    MULADD(at[23], at[45]);    MULADD(at[24], at[44]);    MULADD(at[25], at[43]);    MULADD(at[26], at[42]);    MULADD(at[27], at[41]);    MULADD(at[28], at[40]);    MULADD(at[29], at[39]);    MULADD(at[30], at[38]);    MULADD(at[31], at[37]); 
   COMBA_STORE(C->dp[36]);
   /* 37 */
   COMBA_FORWARD;
   MULADD(at[6], at[63]);    MULADD(at[7], at[62]);    MULADD(at[8], at[61]);    MULADD(at[9], at[60]);    MULADD(at[10], at[59]);    MULADD(at[11], at[58]);    MULADD(at[12], at[57]);    MULADD(at[13], at[56]);    MULADD(at[14], at[55]);    MULADD(at[15], at[54]);    MULADD(at[16], at[53]);    MULADD(at[17], at[52]);    MULADD(at[18], at[51]);    MULADD(at[19], at[50]);    MULADD(at[20], at[49]);    MULADD(at[21], at[48]);    MULADD(at[22], at[47]);    MULADD(at[23], at[46]);    MULADD(at[24], at[45]);    MULADD(at[25], at[44]);    MULADD(at[26], at[43]);    MULADD(at[27], at[42]);    MULADD(at[28], at[41]);    MULADD(at[29], at[40]);    MULADD(at[30], at[39]);    MULADD(at[31], at[38]); 
   COMBA_STORE(C->dp[37]);
   /* 38 */
   COMBA_FORWARD;
   MULADD(at[7], at[63]);    MULADD(at[8], at[62]);    MULADD(at[9], at[61]);    MULADD(at[10], at[60]);    MULADD(at[11], at[59]);    MULADD(at[12], at[58]);    MULADD(at[13], at[57]);    MULADD(at[14], at[56]);    MULADD(at[15], at[55]);    MULADD(at[16], at[54]);    MULADD(at[17], at[53]);    MULADD(at[18], at[52]);    MULADD(at[19], at[51]);    MULADD(at[20], at[50]);    MULADD(at[21], at[49]);    MULADD(at[22], at[48]);    MULADD(at[23], at[47]);    MULADD(at[24], at[46]);    MULADD(at[25], at[45]);    MULADD(at[26], at[44]);    MULADD(at[27], at[43]);    MULADD(at[28], at[42]);    MULADD(at[29], at[41]);    MULADD(at[30], at[40]);    MULADD(at[31], at[39]); 
   COMBA_STORE(C->dp[38]);
   /* 39 */
   COMBA_FORWARD;
   MULADD(at[8], at[63]);    MULADD(at[9], at[62]);    MULADD(at[10], at[61]);    MULADD(at[11], at[60]);    MULADD(at[12], at[59]);    MULADD(at[13], at[58]);    MULADD(at[14], at[57]);    MULADD(at[15], at[56]);    MULADD(at[16], at[55]);    MULADD(at[17], at[54]);    MULADD(at[18], at[53]);    MULADD(at[19], at[52]);    MULADD(at[20], at[51]);    MULADD(at[21], at[50]);    MULADD(at[22], at[49]);    MULADD(at[23], at[48]);    MULADD(at[24], at[47]);    MULADD(at[25], at[46]);    MULADD(at[26], at[45]);    MULADD(at[27], at[44]);    MULADD(at[28], at[43]);    MULADD(at[29], at[42]);    MULADD(at[30], at[41]);    MULADD(at[31], at[40]); 
   COMBA_STORE(C->dp[39]);
   /* 40 */
   COMBA_FORWARD;
   MULADD(at[9], at[63]);    MULADD(at[10], at[62]);    MULADD(at[11], at[61]);    MULADD(at[12], at[60]);    MULADD(at[13], at[59]);    MULADD(at[14], at[58]);    MULADD(at[15], at[57]);    MULADD(at[16], at[56]);    MULADD(at[17], at[55]);    MULADD(at[18], at[54]);    MULADD(at[19], at[53]);    MULADD(at[20], at[52]);    MULADD(at[21], at[51]);    MULADD(at[22], at[50]);    MULADD(at[23], at[49]);    MULADD(at[24], at[48]);    MULADD(at[25], at[47]);    MULADD(at[26], at[46]);    MULADD(at[27], at[45]);    MULADD(at[28], at[44]);    MULADD(at[29], at[43]);    MULADD(at[30], at[42]);    MULADD(at[31], at[41]); 
   COMBA_STORE(C->dp[40]);
   /* 41 */
   COMBA_FORWARD;
   MULADD(at[10], at[63]);    MULADD(at[11], at[62]);    MULADD(at[12], at[61]);    MULADD(at[13], at[60]);    MULADD(at[14], at[59]);    MULADD(at[15], at[58]);    MULADD(at[16], at[57]);    MULADD(at[17], at[56]);    MULADD(at[18], at[55]);    MULADD(at[19], at[54]);    MULADD(at[20], at[53]);    MULADD(at[21], at[52]);    MULADD(at[22], at[51]);    MULADD(at[23], at[50]);    MULADD(at[24], at[49]);    MULADD(at[25], at[48]);    MULADD(at[26], at[47]);    MULADD(at[27], at[46]);    MULADD(at[28], at[45]);    MULADD(at[29], at[44]);    MULADD(at[30], at[43]);    MULADD(at[31], at[42]); 
   COMBA_STORE(C->dp[41]);
   /* 42 */
   COMBA_FORWARD;
   MULADD(at[11], at[63]);    MULADD(at[12], at[62]);    MULADD(at[13], at[61]);    MULADD(at[14], at[60]);    MULADD(at[15], at[59]);    MULADD(at[16], at[58]);    MULADD(at[17], at[57]);    MULADD(at[18], at[56]);    MULADD(at[19], at[55]);    MULADD(at[20], at[54]);    MULADD(at[21], at[53]);    MULADD(at[22], at[52]);    MULADD(at[23], at[51]);    MULADD(at[24], at[50]);    MULADD(at[25], at[49]);    MULADD(at[26], at[48]);    MULADD(at[27], at[47]);    MULADD(at[28], at[46]);    MULADD(at[29], at[45]);    MULADD(at[30], at[44]);    MULADD(at[31], at[43]); 
   COMBA_STORE(C->dp[42]);
   /* 43 */
   COMBA_FORWARD;
   MULADD(at[12], at[63]);    MULADD(at[13], at[62]);    MULADD(at[14], at[61]);    MULADD(at[15], at[60]);    MULADD(at[16], at[59]);    MULADD(at[17], at[58]);    MULADD(at[18], at[57]);    MULADD(at[19], at[56]);    MULADD(at[20], at[55]);    MULADD(at[21], at[54]);    MULADD(at[22], at[53]);    MULADD(at[23], at[52]);    MULADD(at[24], at[51]);    MULADD(at[25], at[50]);    MULADD(at[26], at[49]);    MULADD(at[27], at[48]);    MULADD(at[28], at[47]);    MULADD(at[29], at[46]);    MULADD(at[30], at[45]);    MULADD(at[31], at[44]); 
   COMBA_STORE(C->dp[43]);
   /* 44 */
   COMBA_FORWARD;
   MULADD(at[13], at[63]);    MULADD(at[14], at[62]);    MULADD(at[15], at[61]);    MULADD(at[16], at[60]);    MULADD(at[17], at[59]);    MULADD(at[18], at[58]);    MULADD(at[19], at[57]);    MULADD(at[20], at[56]);    MULADD(at[21], at[55]);    MULADD(at[22], at[54]);    MULADD(at[23], at[53]);    MULADD(at[24], at[52]);    MULADD(at[25], at[51]);    MULADD(at[26], at[50]);    MULADD(at[27], at[49]);    MULADD(at[28], at[48]);    MULADD(at[29], at[47]);    MULADD(at[30], at[46]);    MULADD(at[31], at[45]); 
   COMBA_STORE(C->dp[44]);
   /* 45 */
   COMBA_FORWARD;
   MULADD(at[14], at[63]);    MULADD(at[15], at[62]);    MULADD(at[16], at[61]);    MULADD(at[17], at[60]);    MULADD(at[18], at[59]);    MULADD(at[19], at[58]);    MULADD(at[20], at[57]);    MULADD(at[21], at[56]);    MULADD(at[22], at[55]);    MULADD(at[23], at[54]);    MULADD(at[24], at[53]);    MULADD(at[25], at[52]);    MULADD(at[26], at[51]);    MULADD(at[27], at[50]);    MULADD(at[28], at[49]);    MULADD(at[29], at[48]);    MULADD(at[30], at[47]);    MULADD(at[31], at[46]); 
   COMBA_STORE(C->dp[45]);
   /* 46 */
   COMBA_FORWARD;
   MULADD(at[15], at[63]);    MULADD(at[16], at[62]);    MULADD(at[17], at[61]);    MULADD(at[18], at[60]);    MULADD(at[19], at[59]);    MULADD(at[20], at[58]);    MULADD(at[21], at[57]);    MULADD(at[22], at[56]);    MULADD(at[23], at[55]);    MULADD(at[24], at[54]);    MULADD(at[25], at[53]);    MULADD(at[26], at[52]);    MULADD(at[27], at[51]);    MULADD(at[28], at[50]);    MULADD(at[29], at[49]);    MULADD(at[30], at[48]);    MULADD(at[31], at[47]); 
   COMBA_STORE(C->dp[46]);
   /* 47 */
   COMBA_FORWARD;
   MULADD(at[16], at[63]);    MULADD(at[17], at[62]);    MULADD(at[18], at[61]);    MULADD(at[19], at[60]);    MULADD(at[20], at[59]);    MULADD(at[21], at[58]);    MULADD(at[22], at[57]);    MULADD(at[23], at[56]);    MULADD(at[24], at[55]);    MULADD(at[25], at[54]);    MULADD(at[26], at[53]);    MULADD(at[27], at[52]);    MULADD(at[28], at[51]);    MULADD(at[29], at[50]);    MULADD(at[30], at[49]);    MULADD(at[31], at[48]); 
   COMBA_STORE(C->dp[47]);
   /* 48 */
   COMBA_FORWARD;
   MULADD(at[17], at[63]);    MULADD(at[18], at[62]);    MULADD(at[19], at[61]);    MULADD(at[20], at[60]);    MULADD(at[21], at[59]);    MULADD(at[22], at[58]);    MULADD(at[23], at[57]);    MULADD(at[24], at[56]);    MULADD(at[25], at[55]);    MULADD(at[26], at[54]);    MULADD(at[27], at[53]);    MULADD(at[28], at[52]);    MULADD(at[29], at[51]);    MULADD(at[30], at[50]);    MULADD(at[31], at[49]); 
   COMBA_STORE(C->dp[48]);
   /* 49 */
   COMBA_FORWARD;
   MULADD(at[18], at[63]);    MULADD(at[19], at[62]);    MULADD(at[20], at[61]);    MULADD(at[21], at[60]);    MULADD(at[22], at[59]);    MULADD(at[23], at[58]);    MULADD(at[24], at[57]);    MULADD(at[25], at[56]);    MULADD(at[26], at[55]);    MULADD(at[27], at[54]);    MULADD(at[28], at[53]);    MULADD(at[29], at[52]);    MULADD(at[30], at[51]);    MULADD(at[31], at[50]); 
   COMBA_STORE(C->dp[49]);
   /* 50 */
   COMBA_FORWARD;
   MULADD(at[19], at[63]);    MULADD(at[20], at[62]);    MULADD(at[21], at[61]);    MULADD(at[22], at[60]);    MULADD(at[23], at[59]);    MULADD(at[24], at[58]);    MULADD(at[25], at[57]);    MULADD(at[26], at[56]);    MULADD(at[27], at[55]);    MULADD(at[28], at[54]);    MULADD(at[29], at[53]);    MULADD(at[30], at[52]);    MULADD(at[31], at[51]); 
   COMBA_STORE(C->dp[50]);
   /* 51 */
   COMBA_FORWARD;
   MULADD(at[20], at[63]);    MULADD(at[21], at[62]);    MULADD(at[22], at[61]);    MULADD(at[23], at[60]);    MULADD(at[24], at[59]);    MULADD(at[25], at[58]);    MULADD(at[26], at[57]);    MULADD(at[27], at[56]);    MULADD(at[28], at[55]);    MULADD(at[29], at[54]);    MULADD(at[30], at[53]);    MULADD(at[31], at[52]); 
   COMBA_STORE(C->dp[51]);
   /* 52 */
   COMBA_FORWARD;
   MULADD(at[21], at[63]);    MULADD(at[22], at[62]);    MULADD(at[23], at[61]);    MULADD(at[24], at[60]);    MULADD(at[25], at[59]);    MULADD(at[26], at[58]);    MULADD(at[27], at[57]);    MULADD(at[28], at[56]);    MULADD(at[29], at[55]);    MULADD(at[30], at[54]);    MULADD(at[31], at[53]); 
   COMBA_STORE(C->dp[52]);
   /* 53 */
   COMBA_FORWARD;
   MULADD(at[22], at[63]);    MULADD(at[23], at[62]);    MULADD(at[24], at[61]);    MULADD(at[25], at[60]);    MULADD(at[26], at[59]);    MULADD(at[27], at[58]);    MULADD(at[28], at[57]);    MULADD(at[29], at[56]);    MULADD(at[30], at[55]);    MULADD(at[31], at[54]); 
   COMBA_STORE(C->dp[53]);
   /* 54 */
   COMBA_FORWARD;
   MULADD(at[23], at[63]);    MULADD(at[24], at[62]);    MULADD(at[25], at[61]);    MULADD(at[26], at[60]);    MULADD(at[27], at[59]);    MULADD(at[28], at[58]);    MULADD(at[29], at[57]);    MULADD(at[30], at[56]);    MULADD(at[31], at[55]); 
   COMBA_STORE(C->dp[54]);
   /* 55 */
   COMBA_FORWARD;
   MULADD(at[24], at[63]);    MULADD(at[25], at[62]);    MULADD(at[26], at[61]);    MULADD(at[27], at[60]);    MULADD(at[28], at[59]);    MULADD(at[29], at[58]);    MULADD(at[30], at[57]);    MULADD(at[31], at[56]); 
   COMBA_STORE(C->dp[55]);
   /* 56 */
   COMBA_FORWARD;
   MULADD(at[25], at[63]);    MULADD(at[26], at[62]);    MULADD(at[27], at[61]);    MULADD(at[28], at[60]);    MULADD(at[29], at[59]);    MULADD(at[30], at[58]);    MULADD(at[31], at[57]); 
   COMBA_STORE(C->dp[56]);
   /* 57 */
   COMBA_FORWARD;
   MULADD(at[26], at[63]);    MULADD(at[27], at[62]);    MULADD(at[28], at[61]);    MULADD(at[29], at[60]);    MULADD(at[30], at[59]);    MULADD(at[31], at[58]); 
   COMBA_STORE(C->dp[57]);
   /* 58 */
   COMBA_FORWARD;
   MULADD(at[27], at[63]);    MULADD(at[28], at[62]);    MULADD(at[29], at[61]);    MULADD(at[30], at[60]);    MULADD(at[31], at[59]); 
   COMBA_STORE(C->dp[58]);
   /* 59 */
   COMBA_FORWARD;
   MULADD(at[28], at[63]);    MULADD(at[29], at[62]);    MULADD(at[30], at[61]);    MULADD(at[31], at[60]); 
   COMBA_STORE(C->dp[59]);
   /* 60 */
   COMBA_FORWARD;
   MULADD(at[29], at[63]);    MULADD(at[30], at[62]);    MULADD(at[31], at[61]); 
   COMBA_STORE(C->dp[60]);
   /* 61 */
   COMBA_FORWARD;
   MULADD(at[30], at[63]);    MULADD(at[31], at[62]); 
   COMBA_STORE(C->dp[61]);
   /* 62 */
   COMBA_FORWARD;
   MULADD(at[31], at[63]); 
   COMBA_STORE(C->dp[62]);
   COMBA_STORE2(C->dp[63]);
   C->used = 64;
   C->sign = A->sign ^ B->sign;
   fp_clamp(C);
   COMBA_FINI;
}

#endif
