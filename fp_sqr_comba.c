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

/* About this file...
*/

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
   c0 = c1; c1 = c2; c2 = 0;

#define COMBA_FINI

#define SQRADD(i, j)                                      \
asm volatile (                                            \
     "movl  %6,%%eax     \n\t"                            \
     "mull  %%eax        \n\t"                            \
     "addl  %%eax,%0     \n\t"                            \
     "adcl  %%edx,%1     \n\t"                            \
     "adcl  $0,%2        \n\t"                            \
     :"=r"(c0), "=r"(c1), "=r"(c2): "0"(c0), "1"(c1), "2"(c2), "m"(i) :"%eax","%edx","%cc");

#define SQRADD2(i, j)                                     \
asm volatile (                                            \
     "movl  %6,%%eax     \n\t"                            \
     "mull  %7           \n\t"                            \
     "addl  %%eax,%0     \n\t"                            \
     "adcl  %%edx,%1     \n\t"                            \
     "adcl  $0,%2        \n\t"                            \
     "addl  %%eax,%0     \n\t"                            \
     "adcl  %%edx,%1     \n\t"                            \
     "adcl  $0,%2        \n\t"                            \
     :"=r"(c0), "=r"(c1), "=r"(c2): "0"(c0), "1"(c1), "2"(c2), "m"(i), "m"(j)  :"%eax","%edx","%cc");

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
   c0 = c1; c1 = c2; c2 = 0;

#define COMBA_FINI

#define SQRADD(i, j)                                      \
asm volatile (                                            \
     "movq  %6,%%rax     \n\t"                            \
     "mulq  %%rax        \n\t"                            \
     "addq  %%rax,%0     \n\t"                            \
     "adcq  %%rdx,%1     \n\t"                            \
     "adcq  $0,%2        \n\t"                            \
     :"=r"(c0), "=r"(c1), "=r"(c2): "0"(c0), "1"(c1), "2"(c2), "m"(i) :"%rax","%rdx","%cc");

#define SQRADD2(i, j)                                     \
asm volatile (                                            \
     "movq  %6,%%rax     \n\t"                            \
     "mulq  %7           \n\t"                            \
     "addq  %%rax,%0     \n\t"                            \
     "adcq  %%rdx,%1     \n\t"                            \
     "adcq  $0,%2        \n\t"                            \
     "addq  %%rax,%0     \n\t"                            \
     "adcq  %%rdx,%1     \n\t"                            \
     "adcq  $0,%2        \n\t"                            \
     :"=r"(c0), "=r"(c1), "=r"(c2): "0"(c0), "1"(c1), "2"(c2), "m"(i), "m"(j)  :"%rax","%rdx","%cc");


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
   c0 = c1; c1 = c2; c2 = 0;

#define COMBA_FINI \
   asm("emms");

#define SQRADD(i, j)                                      \
asm volatile (                                            \
     "movd  %6,%%mm0     \n\t"                            \
     "pmuludq %%mm0,%%mm0\n\t"                            \
     "movd  %%mm0,%%eax  \n\t"                            \
     "psrlq $32,%%mm0    \n\t"                            \
     "movd  %%mm0,%%edx  \n\t"                            \
     "addl  %%eax,%0     \n\t"                            \
     "adcl  %%edx,%1     \n\t"                            \
     "adcl  $0,%2        \n\t"                            \
     :"=r"(c0), "=r"(c1), "=r"(c2): "0"(c0), "1"(c1), "2"(c2), "m"(i) :"%eax","%edx","%cc");

#define SQRADD2(i, j)                                     \
asm volatile (                                            \
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
   c0 = c1; c1 = c2; c2 = 0;

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

#else

/* ISO C portable code */

#define COMBA_START

#define CLEAR_CARRY \
   c0 = c1 = c2 = 0;

#define COMBA_STORE(x) \
   x = c0;

#define COMBA_STORE2(x) \
   x = c1;

#define CARRY_FORWARD \
   c0 = c1; c1 = c2; c2 = 0;

#define COMBA_FINI

/* multiplies point i and j, updates carry "c1" and digit c2 */
#define SQRADD(i, j)                       \
   t  = ((fp_word)i) * ((fp_word)j);       \
   c0 = (c0 + t);              if (c0 < ((fp_digit)t))  ++c1; \
   c1 = (c1 + (t>>DIGIT_BIT)); if (c1 < (t>>DIGIT_BIT)) ++c2; 

/* for squaring some of the terms are doubled... */
#define SQRADD2(i, j)                       \
   t  = ((fp_word)i) * ((fp_word)j);       \
   c0 = (c0 + t);              if (c0 < ((fp_digit)t))  ++c1; \
   c1 = (c1 + (t>>DIGIT_BIT)); if (c1 < (t>>DIGIT_BIT)) ++c2; \
   c0 = (c0 + t);              if (c0 < ((fp_digit)t))  ++c1; \
   c1 = (c1 + (t>>DIGIT_BIT)); if (c1 < (t>>DIGIT_BIT)) ++c2; 

#endif

/* generic comba squarer */
void fp_sqr_comba(fp_int *A, fp_int *B)
{
  int       pa, ix, iz;
  fp_digit  c0, c1, c2;
  fp_int    tmp, *dst;
  fp_word   t;

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

      /* this is the number of times the loop will iterrate, essentially its 
         while (tx++ < a->used && ty-- >= 0) { ... }
       */
      iy = MIN(A->used-tx, ty+1);

      /* now for squaring tx can never equal ty 
       * we halve the distance since they approach at a rate of 2x
       * and we have to round because odd cases need to be executed
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
  COMBA_STORE2(dst->dp[ix]);

  COMBA_FINI;

  /* setup dest */
  dst->used = pa;
  fp_clamp (dst);
  if (dst != B) {
     fp_copy(dst, B);
  }
}

void fp_sqr_comba4(fp_int *A, fp_int *B)
{
   fp_word t;
   fp_digit *a, b[8], c0, c1, c2;

   a = A->dp;
   COMBA_START; 

   /* clear carries */
   CLEAR_CARRY;

   /* output 0 */
   SQRADD(a[0],a[0]);
   COMBA_STORE(b[0]);

   /* output 1 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[1]); 
   COMBA_STORE(b[1]);

   /* output 2 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[2]); SQRADD(a[1], a[1]); 
   COMBA_STORE(b[2]);

   /* output 3 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[3]); SQRADD2(a[1], a[2]); 
   COMBA_STORE(b[3]);

   /* output 4 */
   CARRY_FORWARD;
   SQRADD2(a[1], a[3]); SQRADD(a[2], a[2]); 
   COMBA_STORE(b[4]);

   /* output 5 */
   CARRY_FORWARD;
   SQRADD2(a[2], a[3]); 
   COMBA_STORE(b[5]);

   /* output 6 */
   CARRY_FORWARD;
   SQRADD(a[3], a[3]); 
   COMBA_STORE(b[6]);
   COMBA_STORE2(b[7]);
   COMBA_FINI;

   B->used = 8;
   B->sign = FP_ZPOS;
   memcpy(B->dp, b, 8 * sizeof(fp_digit));
   fp_clamp(B);
}


void fp_sqr_comba8(fp_int *A, fp_int *B)
{
   fp_word t;
   fp_digit *a, b[16], c0, c1, c2;

   a = A->dp;
   COMBA_START; 

   /* clear carries */
   CLEAR_CARRY;

   /* output 0 */
   SQRADD(a[0],a[0]);
   COMBA_STORE(b[0]);

   /* output 1 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[1]); 
   COMBA_STORE(b[1]);

   /* output 2 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[2]); SQRADD(a[1], a[1]); 
   COMBA_STORE(b[2]);

   /* output 3 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[3]); SQRADD2(a[1], a[2]); 
   COMBA_STORE(b[3]);

   /* output 4 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[4]); SQRADD2(a[1], a[3]); SQRADD(a[2], a[2]); 
   COMBA_STORE(b[4]);

   /* output 5 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[5]); SQRADD2(a[1], a[4]); SQRADD2(a[2], a[3]); 
   COMBA_STORE(b[5]);

   /* output 6 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[6]); SQRADD2(a[1], a[5]); SQRADD2(a[2], a[4]); SQRADD(a[3], a[3]); 
   COMBA_STORE(b[6]);

   /* output 7 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[7]); SQRADD2(a[1], a[6]); SQRADD2(a[2], a[5]); SQRADD2(a[3], a[4]); 
   COMBA_STORE(b[7]);

   /* output 8 */
   CARRY_FORWARD;
   SQRADD2(a[1], a[7]); SQRADD2(a[2], a[6]); SQRADD2(a[3], a[5]); SQRADD(a[4], a[4]); 
   COMBA_STORE(b[8]);

   /* output 9 */
   CARRY_FORWARD;
   SQRADD2(a[2], a[7]); SQRADD2(a[3], a[6]); SQRADD2(a[4], a[5]); 
   COMBA_STORE(b[9]);

   /* output 10 */
   CARRY_FORWARD;
   SQRADD2(a[3], a[7]); SQRADD2(a[4], a[6]); SQRADD(a[5], a[5]); 
   COMBA_STORE(b[10]);

   /* output 11 */
   CARRY_FORWARD;
   SQRADD2(a[4], a[7]); SQRADD2(a[5], a[6]); 
   COMBA_STORE(b[11]);

   /* output 12 */
   CARRY_FORWARD;
   SQRADD2(a[5], a[7]); SQRADD(a[6], a[6]); 
   COMBA_STORE(b[12]);

   /* output 13 */
   CARRY_FORWARD;
   SQRADD2(a[6], a[7]); 
   COMBA_STORE(b[13]);

   /* output 14 */
   CARRY_FORWARD;
   SQRADD(a[7], a[7]); 
   COMBA_STORE(b[14]);
   COMBA_STORE2(b[15]);
   COMBA_FINI;

   B->used = 16;
   B->sign = FP_ZPOS;
   memcpy(B->dp, b, 16 * sizeof(fp_digit));
   fp_clamp(B);
}


void fp_sqr_comba16(fp_int *A, fp_int *B)
{
   fp_word t;
   fp_digit *a, b[32], c0, c1, c2;

   a = A->dp;
   COMBA_START; 

   /* clear carries */
   CLEAR_CARRY;

   /* output 0 */
   SQRADD(a[0],a[0]);
   COMBA_STORE(b[0]);

   /* output 1 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[1]); 
   COMBA_STORE(b[1]);

   /* output 2 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[2]); SQRADD(a[1], a[1]); 
   COMBA_STORE(b[2]);

   /* output 3 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[3]); SQRADD2(a[1], a[2]); 
   COMBA_STORE(b[3]);

   /* output 4 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[4]); SQRADD2(a[1], a[3]); SQRADD(a[2], a[2]); 
   COMBA_STORE(b[4]);

   /* output 5 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[5]); SQRADD2(a[1], a[4]); SQRADD2(a[2], a[3]); 
   COMBA_STORE(b[5]);

   /* output 6 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[6]); SQRADD2(a[1], a[5]); SQRADD2(a[2], a[4]); SQRADD(a[3], a[3]); 
   COMBA_STORE(b[6]);

   /* output 7 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[7]); SQRADD2(a[1], a[6]); SQRADD2(a[2], a[5]); SQRADD2(a[3], a[4]); 
   COMBA_STORE(b[7]);

   /* output 8 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[8]); SQRADD2(a[1], a[7]); SQRADD2(a[2], a[6]); SQRADD2(a[3], a[5]); SQRADD(a[4], a[4]); 
   COMBA_STORE(b[8]);

   /* output 9 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[9]); SQRADD2(a[1], a[8]); SQRADD2(a[2], a[7]); SQRADD2(a[3], a[6]); SQRADD2(a[4], a[5]); 
   COMBA_STORE(b[9]);

   /* output 10 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[10]); SQRADD2(a[1], a[9]); SQRADD2(a[2], a[8]); SQRADD2(a[3], a[7]); SQRADD2(a[4], a[6]); SQRADD(a[5], a[5]); 
   COMBA_STORE(b[10]);

   /* output 11 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[11]); SQRADD2(a[1], a[10]); SQRADD2(a[2], a[9]); SQRADD2(a[3], a[8]); SQRADD2(a[4], a[7]); SQRADD2(a[5], a[6]); 
   COMBA_STORE(b[11]);

   /* output 12 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[12]); SQRADD2(a[1], a[11]); SQRADD2(a[2], a[10]); SQRADD2(a[3], a[9]); SQRADD2(a[4], a[8]); SQRADD2(a[5], a[7]); SQRADD(a[6], a[6]); 
   COMBA_STORE(b[12]);

   /* output 13 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[13]); SQRADD2(a[1], a[12]); SQRADD2(a[2], a[11]); SQRADD2(a[3], a[10]); SQRADD2(a[4], a[9]); SQRADD2(a[5], a[8]); SQRADD2(a[6], a[7]); 
   COMBA_STORE(b[13]);

   /* output 14 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[14]); SQRADD2(a[1], a[13]); SQRADD2(a[2], a[12]); SQRADD2(a[3], a[11]); SQRADD2(a[4], a[10]); SQRADD2(a[5], a[9]); SQRADD2(a[6], a[8]); SQRADD(a[7], a[7]); 
   COMBA_STORE(b[14]);

   /* output 15 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[15]); SQRADD2(a[1], a[14]); SQRADD2(a[2], a[13]); SQRADD2(a[3], a[12]); SQRADD2(a[4], a[11]); SQRADD2(a[5], a[10]); SQRADD2(a[6], a[9]); SQRADD2(a[7], a[8]); 
   COMBA_STORE(b[15]);

   /* output 16 */
   CARRY_FORWARD;
   SQRADD2(a[1], a[15]); SQRADD2(a[2], a[14]); SQRADD2(a[3], a[13]); SQRADD2(a[4], a[12]); SQRADD2(a[5], a[11]); SQRADD2(a[6], a[10]); SQRADD2(a[7], a[9]); SQRADD(a[8], a[8]); 
   COMBA_STORE(b[16]);

   /* output 17 */
   CARRY_FORWARD;
   SQRADD2(a[2], a[15]); SQRADD2(a[3], a[14]); SQRADD2(a[4], a[13]); SQRADD2(a[5], a[12]); SQRADD2(a[6], a[11]); SQRADD2(a[7], a[10]); SQRADD2(a[8], a[9]); 
   COMBA_STORE(b[17]);

   /* output 18 */
   CARRY_FORWARD;
   SQRADD2(a[3], a[15]); SQRADD2(a[4], a[14]); SQRADD2(a[5], a[13]); SQRADD2(a[6], a[12]); SQRADD2(a[7], a[11]); SQRADD2(a[8], a[10]); SQRADD(a[9], a[9]); 
   COMBA_STORE(b[18]);

   /* output 19 */
   CARRY_FORWARD;
   SQRADD2(a[4], a[15]); SQRADD2(a[5], a[14]); SQRADD2(a[6], a[13]); SQRADD2(a[7], a[12]); SQRADD2(a[8], a[11]); SQRADD2(a[9], a[10]); 
   COMBA_STORE(b[19]);

   /* output 20 */
   CARRY_FORWARD;
   SQRADD2(a[5], a[15]); SQRADD2(a[6], a[14]); SQRADD2(a[7], a[13]); SQRADD2(a[8], a[12]); SQRADD2(a[9], a[11]); SQRADD(a[10], a[10]); 
   COMBA_STORE(b[20]);

   /* output 21 */
   CARRY_FORWARD;
   SQRADD2(a[6], a[15]); SQRADD2(a[7], a[14]); SQRADD2(a[8], a[13]); SQRADD2(a[9], a[12]); SQRADD2(a[10], a[11]); 
   COMBA_STORE(b[21]);

   /* output 22 */
   CARRY_FORWARD;
   SQRADD2(a[7], a[15]); SQRADD2(a[8], a[14]); SQRADD2(a[9], a[13]); SQRADD2(a[10], a[12]); SQRADD(a[11], a[11]); 
   COMBA_STORE(b[22]);

   /* output 23 */
   CARRY_FORWARD;
   SQRADD2(a[8], a[15]); SQRADD2(a[9], a[14]); SQRADD2(a[10], a[13]); SQRADD2(a[11], a[12]); 
   COMBA_STORE(b[23]);

   /* output 24 */
   CARRY_FORWARD;
   SQRADD2(a[9], a[15]); SQRADD2(a[10], a[14]); SQRADD2(a[11], a[13]); SQRADD(a[12], a[12]); 
   COMBA_STORE(b[24]);

   /* output 25 */
   CARRY_FORWARD;
   SQRADD2(a[10], a[15]); SQRADD2(a[11], a[14]); SQRADD2(a[12], a[13]); 
   COMBA_STORE(b[25]);

   /* output 26 */
   CARRY_FORWARD;
   SQRADD2(a[11], a[15]); SQRADD2(a[12], a[14]); SQRADD(a[13], a[13]); 
   COMBA_STORE(b[26]);

   /* output 27 */
   CARRY_FORWARD;
   SQRADD2(a[12], a[15]); SQRADD2(a[13], a[14]); 
   COMBA_STORE(b[27]);

   /* output 28 */
   CARRY_FORWARD;
   SQRADD2(a[13], a[15]); SQRADD(a[14], a[14]); 
   COMBA_STORE(b[28]);

   /* output 29 */
   CARRY_FORWARD;
   SQRADD2(a[14], a[15]); 
   COMBA_STORE(b[29]);

   /* output 30 */
   CARRY_FORWARD;
   SQRADD(a[15], a[15]); 
   COMBA_STORE(b[30]);
   COMBA_STORE2(b[31]);
   COMBA_FINI;

   B->used = 32;
   B->sign = FP_ZPOS;
   memcpy(B->dp, b, 32 * sizeof(fp_digit));
   fp_clamp(B);
}

#ifdef TFM_HUGE

void fp_sqr_comba32(fp_int *A, fp_int *B)
{
   fp_word t;
   fp_digit *a, b[64], c0, c1, c2;

   a = A->dp;
   COMBA_START; 

   /* clear carries */
   CLEAR_CARRY;

   /* output 0 */
   SQRADD(a[0],a[0]);
   COMBA_STORE(b[0]);

   /* output 1 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[1]); 
   COMBA_STORE(b[1]);

   /* output 2 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[2]); SQRADD(a[1], a[1]); 
   COMBA_STORE(b[2]);

   /* output 3 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[3]); SQRADD2(a[1], a[2]); 
   COMBA_STORE(b[3]);

   /* output 4 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[4]); SQRADD2(a[1], a[3]); SQRADD(a[2], a[2]); 
   COMBA_STORE(b[4]);

   /* output 5 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[5]); SQRADD2(a[1], a[4]); SQRADD2(a[2], a[3]); 
   COMBA_STORE(b[5]);

   /* output 6 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[6]); SQRADD2(a[1], a[5]); SQRADD2(a[2], a[4]); SQRADD(a[3], a[3]); 
   COMBA_STORE(b[6]);

   /* output 7 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[7]); SQRADD2(a[1], a[6]); SQRADD2(a[2], a[5]); SQRADD2(a[3], a[4]); 
   COMBA_STORE(b[7]);

   /* output 8 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[8]); SQRADD2(a[1], a[7]); SQRADD2(a[2], a[6]); SQRADD2(a[3], a[5]); SQRADD(a[4], a[4]); 
   COMBA_STORE(b[8]);

   /* output 9 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[9]); SQRADD2(a[1], a[8]); SQRADD2(a[2], a[7]); SQRADD2(a[3], a[6]); SQRADD2(a[4], a[5]); 
   COMBA_STORE(b[9]);

   /* output 10 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[10]); SQRADD2(a[1], a[9]); SQRADD2(a[2], a[8]); SQRADD2(a[3], a[7]); SQRADD2(a[4], a[6]); SQRADD(a[5], a[5]); 
   COMBA_STORE(b[10]);

   /* output 11 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[11]); SQRADD2(a[1], a[10]); SQRADD2(a[2], a[9]); SQRADD2(a[3], a[8]); SQRADD2(a[4], a[7]); SQRADD2(a[5], a[6]); 
   COMBA_STORE(b[11]);

   /* output 12 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[12]); SQRADD2(a[1], a[11]); SQRADD2(a[2], a[10]); SQRADD2(a[3], a[9]); SQRADD2(a[4], a[8]); SQRADD2(a[5], a[7]); SQRADD(a[6], a[6]); 
   COMBA_STORE(b[12]);

   /* output 13 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[13]); SQRADD2(a[1], a[12]); SQRADD2(a[2], a[11]); SQRADD2(a[3], a[10]); SQRADD2(a[4], a[9]); SQRADD2(a[5], a[8]); SQRADD2(a[6], a[7]); 
   COMBA_STORE(b[13]);

   /* output 14 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[14]); SQRADD2(a[1], a[13]); SQRADD2(a[2], a[12]); SQRADD2(a[3], a[11]); SQRADD2(a[4], a[10]); SQRADD2(a[5], a[9]); SQRADD2(a[6], a[8]); SQRADD(a[7], a[7]); 
   COMBA_STORE(b[14]);

   /* output 15 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[15]); SQRADD2(a[1], a[14]); SQRADD2(a[2], a[13]); SQRADD2(a[3], a[12]); SQRADD2(a[4], a[11]); SQRADD2(a[5], a[10]); SQRADD2(a[6], a[9]); SQRADD2(a[7], a[8]); 
   COMBA_STORE(b[15]);

   /* output 16 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[16]); SQRADD2(a[1], a[15]); SQRADD2(a[2], a[14]); SQRADD2(a[3], a[13]); SQRADD2(a[4], a[12]); SQRADD2(a[5], a[11]); SQRADD2(a[6], a[10]); SQRADD2(a[7], a[9]); SQRADD(a[8], a[8]); 
   COMBA_STORE(b[16]);

   /* output 17 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[17]); SQRADD2(a[1], a[16]); SQRADD2(a[2], a[15]); SQRADD2(a[3], a[14]); SQRADD2(a[4], a[13]); SQRADD2(a[5], a[12]); SQRADD2(a[6], a[11]); SQRADD2(a[7], a[10]); SQRADD2(a[8], a[9]); 
   COMBA_STORE(b[17]);

   /* output 18 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[18]); SQRADD2(a[1], a[17]); SQRADD2(a[2], a[16]); SQRADD2(a[3], a[15]); SQRADD2(a[4], a[14]); SQRADD2(a[5], a[13]); SQRADD2(a[6], a[12]); SQRADD2(a[7], a[11]); SQRADD2(a[8], a[10]); SQRADD(a[9], a[9]); 
   COMBA_STORE(b[18]);

   /* output 19 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[19]); SQRADD2(a[1], a[18]); SQRADD2(a[2], a[17]); SQRADD2(a[3], a[16]); SQRADD2(a[4], a[15]); SQRADD2(a[5], a[14]); SQRADD2(a[6], a[13]); SQRADD2(a[7], a[12]); SQRADD2(a[8], a[11]); SQRADD2(a[9], a[10]); 
   COMBA_STORE(b[19]);

   /* output 20 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[20]); SQRADD2(a[1], a[19]); SQRADD2(a[2], a[18]); SQRADD2(a[3], a[17]); SQRADD2(a[4], a[16]); SQRADD2(a[5], a[15]); SQRADD2(a[6], a[14]); SQRADD2(a[7], a[13]); SQRADD2(a[8], a[12]); SQRADD2(a[9], a[11]); SQRADD(a[10], a[10]); 
   COMBA_STORE(b[20]);

   /* output 21 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[21]); SQRADD2(a[1], a[20]); SQRADD2(a[2], a[19]); SQRADD2(a[3], a[18]); SQRADD2(a[4], a[17]); SQRADD2(a[5], a[16]); SQRADD2(a[6], a[15]); SQRADD2(a[7], a[14]); SQRADD2(a[8], a[13]); SQRADD2(a[9], a[12]); SQRADD2(a[10], a[11]); 
   COMBA_STORE(b[21]);

   /* output 22 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[22]); SQRADD2(a[1], a[21]); SQRADD2(a[2], a[20]); SQRADD2(a[3], a[19]); SQRADD2(a[4], a[18]); SQRADD2(a[5], a[17]); SQRADD2(a[6], a[16]); SQRADD2(a[7], a[15]); SQRADD2(a[8], a[14]); SQRADD2(a[9], a[13]); SQRADD2(a[10], a[12]); SQRADD(a[11], a[11]); 
   COMBA_STORE(b[22]);

   /* output 23 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[23]); SQRADD2(a[1], a[22]); SQRADD2(a[2], a[21]); SQRADD2(a[3], a[20]); SQRADD2(a[4], a[19]); SQRADD2(a[5], a[18]); SQRADD2(a[6], a[17]); SQRADD2(a[7], a[16]); SQRADD2(a[8], a[15]); SQRADD2(a[9], a[14]); SQRADD2(a[10], a[13]); SQRADD2(a[11], a[12]); 
   COMBA_STORE(b[23]);

   /* output 24 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[24]); SQRADD2(a[1], a[23]); SQRADD2(a[2], a[22]); SQRADD2(a[3], a[21]); SQRADD2(a[4], a[20]); SQRADD2(a[5], a[19]); SQRADD2(a[6], a[18]); SQRADD2(a[7], a[17]); SQRADD2(a[8], a[16]); SQRADD2(a[9], a[15]); SQRADD2(a[10], a[14]); SQRADD2(a[11], a[13]); SQRADD(a[12], a[12]); 
   COMBA_STORE(b[24]);

   /* output 25 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[25]); SQRADD2(a[1], a[24]); SQRADD2(a[2], a[23]); SQRADD2(a[3], a[22]); SQRADD2(a[4], a[21]); SQRADD2(a[5], a[20]); SQRADD2(a[6], a[19]); SQRADD2(a[7], a[18]); SQRADD2(a[8], a[17]); SQRADD2(a[9], a[16]); SQRADD2(a[10], a[15]); SQRADD2(a[11], a[14]); SQRADD2(a[12], a[13]); 
   COMBA_STORE(b[25]);

   /* output 26 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[26]); SQRADD2(a[1], a[25]); SQRADD2(a[2], a[24]); SQRADD2(a[3], a[23]); SQRADD2(a[4], a[22]); SQRADD2(a[5], a[21]); SQRADD2(a[6], a[20]); SQRADD2(a[7], a[19]); SQRADD2(a[8], a[18]); SQRADD2(a[9], a[17]); SQRADD2(a[10], a[16]); SQRADD2(a[11], a[15]); SQRADD2(a[12], a[14]); SQRADD(a[13], a[13]); 
   COMBA_STORE(b[26]);

   /* output 27 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[27]); SQRADD2(a[1], a[26]); SQRADD2(a[2], a[25]); SQRADD2(a[3], a[24]); SQRADD2(a[4], a[23]); SQRADD2(a[5], a[22]); SQRADD2(a[6], a[21]); SQRADD2(a[7], a[20]); SQRADD2(a[8], a[19]); SQRADD2(a[9], a[18]); SQRADD2(a[10], a[17]); SQRADD2(a[11], a[16]); SQRADD2(a[12], a[15]); SQRADD2(a[13], a[14]); 
   COMBA_STORE(b[27]);

   /* output 28 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[28]); SQRADD2(a[1], a[27]); SQRADD2(a[2], a[26]); SQRADD2(a[3], a[25]); SQRADD2(a[4], a[24]); SQRADD2(a[5], a[23]); SQRADD2(a[6], a[22]); SQRADD2(a[7], a[21]); SQRADD2(a[8], a[20]); SQRADD2(a[9], a[19]); SQRADD2(a[10], a[18]); SQRADD2(a[11], a[17]); SQRADD2(a[12], a[16]); SQRADD2(a[13], a[15]); SQRADD(a[14], a[14]); 
   COMBA_STORE(b[28]);

   /* output 29 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[29]); SQRADD2(a[1], a[28]); SQRADD2(a[2], a[27]); SQRADD2(a[3], a[26]); SQRADD2(a[4], a[25]); SQRADD2(a[5], a[24]); SQRADD2(a[6], a[23]); SQRADD2(a[7], a[22]); SQRADD2(a[8], a[21]); SQRADD2(a[9], a[20]); SQRADD2(a[10], a[19]); SQRADD2(a[11], a[18]); SQRADD2(a[12], a[17]); SQRADD2(a[13], a[16]); SQRADD2(a[14], a[15]); 
   COMBA_STORE(b[29]);

   /* output 30 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[30]); SQRADD2(a[1], a[29]); SQRADD2(a[2], a[28]); SQRADD2(a[3], a[27]); SQRADD2(a[4], a[26]); SQRADD2(a[5], a[25]); SQRADD2(a[6], a[24]); SQRADD2(a[7], a[23]); SQRADD2(a[8], a[22]); SQRADD2(a[9], a[21]); SQRADD2(a[10], a[20]); SQRADD2(a[11], a[19]); SQRADD2(a[12], a[18]); SQRADD2(a[13], a[17]); SQRADD2(a[14], a[16]); SQRADD(a[15], a[15]); 
   COMBA_STORE(b[30]);

   /* output 31 */
   CARRY_FORWARD;
   SQRADD2(a[0], a[31]); SQRADD2(a[1], a[30]); SQRADD2(a[2], a[29]); SQRADD2(a[3], a[28]); SQRADD2(a[4], a[27]); SQRADD2(a[5], a[26]); SQRADD2(a[6], a[25]); SQRADD2(a[7], a[24]); SQRADD2(a[8], a[23]); SQRADD2(a[9], a[22]); SQRADD2(a[10], a[21]); SQRADD2(a[11], a[20]); SQRADD2(a[12], a[19]); SQRADD2(a[13], a[18]); SQRADD2(a[14], a[17]); SQRADD2(a[15], a[16]); 
   COMBA_STORE(b[31]);

   /* output 32 */
   CARRY_FORWARD;
   SQRADD2(a[1], a[31]); SQRADD2(a[2], a[30]); SQRADD2(a[3], a[29]); SQRADD2(a[4], a[28]); SQRADD2(a[5], a[27]); SQRADD2(a[6], a[26]); SQRADD2(a[7], a[25]); SQRADD2(a[8], a[24]); SQRADD2(a[9], a[23]); SQRADD2(a[10], a[22]); SQRADD2(a[11], a[21]); SQRADD2(a[12], a[20]); SQRADD2(a[13], a[19]); SQRADD2(a[14], a[18]); SQRADD2(a[15], a[17]); SQRADD(a[16], a[16]); 
   COMBA_STORE(b[32]);

   /* output 33 */
   CARRY_FORWARD;
   SQRADD2(a[2], a[31]); SQRADD2(a[3], a[30]); SQRADD2(a[4], a[29]); SQRADD2(a[5], a[28]); SQRADD2(a[6], a[27]); SQRADD2(a[7], a[26]); SQRADD2(a[8], a[25]); SQRADD2(a[9], a[24]); SQRADD2(a[10], a[23]); SQRADD2(a[11], a[22]); SQRADD2(a[12], a[21]); SQRADD2(a[13], a[20]); SQRADD2(a[14], a[19]); SQRADD2(a[15], a[18]); SQRADD2(a[16], a[17]); 
   COMBA_STORE(b[33]);

   /* output 34 */
   CARRY_FORWARD;
   SQRADD2(a[3], a[31]); SQRADD2(a[4], a[30]); SQRADD2(a[5], a[29]); SQRADD2(a[6], a[28]); SQRADD2(a[7], a[27]); SQRADD2(a[8], a[26]); SQRADD2(a[9], a[25]); SQRADD2(a[10], a[24]); SQRADD2(a[11], a[23]); SQRADD2(a[12], a[22]); SQRADD2(a[13], a[21]); SQRADD2(a[14], a[20]); SQRADD2(a[15], a[19]); SQRADD2(a[16], a[18]); SQRADD(a[17], a[17]); 
   COMBA_STORE(b[34]);

   /* output 35 */
   CARRY_FORWARD;
   SQRADD2(a[4], a[31]); SQRADD2(a[5], a[30]); SQRADD2(a[6], a[29]); SQRADD2(a[7], a[28]); SQRADD2(a[8], a[27]); SQRADD2(a[9], a[26]); SQRADD2(a[10], a[25]); SQRADD2(a[11], a[24]); SQRADD2(a[12], a[23]); SQRADD2(a[13], a[22]); SQRADD2(a[14], a[21]); SQRADD2(a[15], a[20]); SQRADD2(a[16], a[19]); SQRADD2(a[17], a[18]); 
   COMBA_STORE(b[35]);

   /* output 36 */
   CARRY_FORWARD;
   SQRADD2(a[5], a[31]); SQRADD2(a[6], a[30]); SQRADD2(a[7], a[29]); SQRADD2(a[8], a[28]); SQRADD2(a[9], a[27]); SQRADD2(a[10], a[26]); SQRADD2(a[11], a[25]); SQRADD2(a[12], a[24]); SQRADD2(a[13], a[23]); SQRADD2(a[14], a[22]); SQRADD2(a[15], a[21]); SQRADD2(a[16], a[20]); SQRADD2(a[17], a[19]); SQRADD(a[18], a[18]); 
   COMBA_STORE(b[36]);

   /* output 37 */
   CARRY_FORWARD;
   SQRADD2(a[6], a[31]); SQRADD2(a[7], a[30]); SQRADD2(a[8], a[29]); SQRADD2(a[9], a[28]); SQRADD2(a[10], a[27]); SQRADD2(a[11], a[26]); SQRADD2(a[12], a[25]); SQRADD2(a[13], a[24]); SQRADD2(a[14], a[23]); SQRADD2(a[15], a[22]); SQRADD2(a[16], a[21]); SQRADD2(a[17], a[20]); SQRADD2(a[18], a[19]); 
   COMBA_STORE(b[37]);

   /* output 38 */
   CARRY_FORWARD;
   SQRADD2(a[7], a[31]); SQRADD2(a[8], a[30]); SQRADD2(a[9], a[29]); SQRADD2(a[10], a[28]); SQRADD2(a[11], a[27]); SQRADD2(a[12], a[26]); SQRADD2(a[13], a[25]); SQRADD2(a[14], a[24]); SQRADD2(a[15], a[23]); SQRADD2(a[16], a[22]); SQRADD2(a[17], a[21]); SQRADD2(a[18], a[20]); SQRADD(a[19], a[19]); 
   COMBA_STORE(b[38]);

   /* output 39 */
   CARRY_FORWARD;
   SQRADD2(a[8], a[31]); SQRADD2(a[9], a[30]); SQRADD2(a[10], a[29]); SQRADD2(a[11], a[28]); SQRADD2(a[12], a[27]); SQRADD2(a[13], a[26]); SQRADD2(a[14], a[25]); SQRADD2(a[15], a[24]); SQRADD2(a[16], a[23]); SQRADD2(a[17], a[22]); SQRADD2(a[18], a[21]); SQRADD2(a[19], a[20]); 
   COMBA_STORE(b[39]);

   /* output 40 */
   CARRY_FORWARD;
   SQRADD2(a[9], a[31]); SQRADD2(a[10], a[30]); SQRADD2(a[11], a[29]); SQRADD2(a[12], a[28]); SQRADD2(a[13], a[27]); SQRADD2(a[14], a[26]); SQRADD2(a[15], a[25]); SQRADD2(a[16], a[24]); SQRADD2(a[17], a[23]); SQRADD2(a[18], a[22]); SQRADD2(a[19], a[21]); SQRADD(a[20], a[20]); 
   COMBA_STORE(b[40]);

   /* output 41 */
   CARRY_FORWARD;
   SQRADD2(a[10], a[31]); SQRADD2(a[11], a[30]); SQRADD2(a[12], a[29]); SQRADD2(a[13], a[28]); SQRADD2(a[14], a[27]); SQRADD2(a[15], a[26]); SQRADD2(a[16], a[25]); SQRADD2(a[17], a[24]); SQRADD2(a[18], a[23]); SQRADD2(a[19], a[22]); SQRADD2(a[20], a[21]); 
   COMBA_STORE(b[41]);

   /* output 42 */
   CARRY_FORWARD;
   SQRADD2(a[11], a[31]); SQRADD2(a[12], a[30]); SQRADD2(a[13], a[29]); SQRADD2(a[14], a[28]); SQRADD2(a[15], a[27]); SQRADD2(a[16], a[26]); SQRADD2(a[17], a[25]); SQRADD2(a[18], a[24]); SQRADD2(a[19], a[23]); SQRADD2(a[20], a[22]); SQRADD(a[21], a[21]); 
   COMBA_STORE(b[42]);

   /* output 43 */
   CARRY_FORWARD;
   SQRADD2(a[12], a[31]); SQRADD2(a[13], a[30]); SQRADD2(a[14], a[29]); SQRADD2(a[15], a[28]); SQRADD2(a[16], a[27]); SQRADD2(a[17], a[26]); SQRADD2(a[18], a[25]); SQRADD2(a[19], a[24]); SQRADD2(a[20], a[23]); SQRADD2(a[21], a[22]); 
   COMBA_STORE(b[43]);

   /* output 44 */
   CARRY_FORWARD;
   SQRADD2(a[13], a[31]); SQRADD2(a[14], a[30]); SQRADD2(a[15], a[29]); SQRADD2(a[16], a[28]); SQRADD2(a[17], a[27]); SQRADD2(a[18], a[26]); SQRADD2(a[19], a[25]); SQRADD2(a[20], a[24]); SQRADD2(a[21], a[23]); SQRADD(a[22], a[22]); 
   COMBA_STORE(b[44]);

   /* output 45 */
   CARRY_FORWARD;
   SQRADD2(a[14], a[31]); SQRADD2(a[15], a[30]); SQRADD2(a[16], a[29]); SQRADD2(a[17], a[28]); SQRADD2(a[18], a[27]); SQRADD2(a[19], a[26]); SQRADD2(a[20], a[25]); SQRADD2(a[21], a[24]); SQRADD2(a[22], a[23]); 
   COMBA_STORE(b[45]);

   /* output 46 */
   CARRY_FORWARD;
   SQRADD2(a[15], a[31]); SQRADD2(a[16], a[30]); SQRADD2(a[17], a[29]); SQRADD2(a[18], a[28]); SQRADD2(a[19], a[27]); SQRADD2(a[20], a[26]); SQRADD2(a[21], a[25]); SQRADD2(a[22], a[24]); SQRADD(a[23], a[23]); 
   COMBA_STORE(b[46]);

   /* output 47 */
   CARRY_FORWARD;
   SQRADD2(a[16], a[31]); SQRADD2(a[17], a[30]); SQRADD2(a[18], a[29]); SQRADD2(a[19], a[28]); SQRADD2(a[20], a[27]); SQRADD2(a[21], a[26]); SQRADD2(a[22], a[25]); SQRADD2(a[23], a[24]); 
   COMBA_STORE(b[47]);

   /* output 48 */
   CARRY_FORWARD;
   SQRADD2(a[17], a[31]); SQRADD2(a[18], a[30]); SQRADD2(a[19], a[29]); SQRADD2(a[20], a[28]); SQRADD2(a[21], a[27]); SQRADD2(a[22], a[26]); SQRADD2(a[23], a[25]); SQRADD(a[24], a[24]); 
   COMBA_STORE(b[48]);

   /* output 49 */
   CARRY_FORWARD;
   SQRADD2(a[18], a[31]); SQRADD2(a[19], a[30]); SQRADD2(a[20], a[29]); SQRADD2(a[21], a[28]); SQRADD2(a[22], a[27]); SQRADD2(a[23], a[26]); SQRADD2(a[24], a[25]); 
   COMBA_STORE(b[49]);

   /* output 50 */
   CARRY_FORWARD;
   SQRADD2(a[19], a[31]); SQRADD2(a[20], a[30]); SQRADD2(a[21], a[29]); SQRADD2(a[22], a[28]); SQRADD2(a[23], a[27]); SQRADD2(a[24], a[26]); SQRADD(a[25], a[25]); 
   COMBA_STORE(b[50]);

   /* output 51 */
   CARRY_FORWARD;
   SQRADD2(a[20], a[31]); SQRADD2(a[21], a[30]); SQRADD2(a[22], a[29]); SQRADD2(a[23], a[28]); SQRADD2(a[24], a[27]); SQRADD2(a[25], a[26]); 
   COMBA_STORE(b[51]);

   /* output 52 */
   CARRY_FORWARD;
   SQRADD2(a[21], a[31]); SQRADD2(a[22], a[30]); SQRADD2(a[23], a[29]); SQRADD2(a[24], a[28]); SQRADD2(a[25], a[27]); SQRADD(a[26], a[26]); 
   COMBA_STORE(b[52]);

   /* output 53 */
   CARRY_FORWARD;
   SQRADD2(a[22], a[31]); SQRADD2(a[23], a[30]); SQRADD2(a[24], a[29]); SQRADD2(a[25], a[28]); SQRADD2(a[26], a[27]); 
   COMBA_STORE(b[53]);

   /* output 54 */
   CARRY_FORWARD;
   SQRADD2(a[23], a[31]); SQRADD2(a[24], a[30]); SQRADD2(a[25], a[29]); SQRADD2(a[26], a[28]); SQRADD(a[27], a[27]); 
   COMBA_STORE(b[54]);

   /* output 55 */
   CARRY_FORWARD;
   SQRADD2(a[24], a[31]); SQRADD2(a[25], a[30]); SQRADD2(a[26], a[29]); SQRADD2(a[27], a[28]); 
   COMBA_STORE(b[55]);

   /* output 56 */
   CARRY_FORWARD;
   SQRADD2(a[25], a[31]); SQRADD2(a[26], a[30]); SQRADD2(a[27], a[29]); SQRADD(a[28], a[28]); 
   COMBA_STORE(b[56]);

   /* output 57 */
   CARRY_FORWARD;
   SQRADD2(a[26], a[31]); SQRADD2(a[27], a[30]); SQRADD2(a[28], a[29]); 
   COMBA_STORE(b[57]);

   /* output 58 */
   CARRY_FORWARD;
   SQRADD2(a[27], a[31]); SQRADD2(a[28], a[30]); SQRADD(a[29], a[29]); 
   COMBA_STORE(b[58]);

   /* output 59 */
   CARRY_FORWARD;
   SQRADD2(a[28], a[31]); SQRADD2(a[29], a[30]); 
   COMBA_STORE(b[59]);

   /* output 60 */
   CARRY_FORWARD;
   SQRADD2(a[29], a[31]); SQRADD(a[30], a[30]); 
   COMBA_STORE(b[60]);

   /* output 61 */
   CARRY_FORWARD;
   SQRADD2(a[30], a[31]); 
   COMBA_STORE(b[61]);

   /* output 62 */
   CARRY_FORWARD;
   SQRADD(a[31], a[31]); 
   COMBA_STORE(b[62]);
   COMBA_STORE2(b[63]);
   COMBA_FINI;

   B->used = 64;
   B->sign = FP_ZPOS;
   memcpy(B->dp, b, 64 * sizeof(fp_digit));
   fp_clamp(B);
}

#endif

