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

#if defined(TFM_X86) 

/* x86-32 code */

#define MONT_START 

#define MONT_FINI

#define LOOP_START \
   mu = c[x] * mp;

#define INNERMUL \
asm(                                                                                          \
"movl %7,%%eax                \n\t"                                                           \
"mull %6                      \n\t"                                                           \
"addl %%eax,%0                \n\t"                                                           \
"adcl %%edx,%1                \n\t"                                                           \
"adcl $0,%2                   \n\t"                                                           \
:"=g"(_c[OFF0]), "=g"(_c[OFF1]), "=g"(_c[OFF2]):"0"(_c[OFF0]), "1"(_c[OFF1]), "2"(_c[OFF2]),  \
                                                "g"(mu), "g"(*tmpm++)                          \
                                               : "%eax", "%edx", "%cc");

#define PROPCARRY \
asm(                                                                                               \
"movl %1,%%eax                \n\t"                                                                \
"addl  %%eax,%6               \n\t"                                                                \
"movl %2,%%eax                \n\t"                                                                \
"adcl  %%eax,%7               \n\t"                                                                \
"adcl $0,%8                   \n\t"                                                                \
:"=g"(_c[OFF0]), "=g"(_c[OFF1]), "=g"(_c[OFF2]):"0"(_c[OFF0]), "1"(_c[OFF1]), "2"(_c[OFF2]),       \
                                                "m"(_c[OFF0+1]), "m"(_c[OFF1+1]), "m"(_c[OFF2+1])  \
: "%eax", "%cc");

#elif defined(TFM_X86_64)
/* x86-64 code */

#define MONT_START 

#define MONT_FINI

#define LOOP_START \
   mu = c[x] * mp;

#define INNERMUL \
asm(                                                                                          \
"movq %7,%%rax                \n\t"                                                           \
"mulq %6                      \n\t"                                                           \
"addq %%rax,%0                \n\t"                                                           \
"adcq %%rdx,%1                \n\t"                                                           \
"adcq $0,%2                   \n\t"                                                           \
:"=g"(_c[OFF0]), "=g"(_c[OFF1]), "=g"(_c[OFF2]):"0"(_c[OFF0]), "1"(_c[OFF1]), "2"(_c[OFF2]),  \
                                                "g"(mu), "g"(*tmpm++)                          \
                                               : "%rax", "%rdx", "%cc");

#define PROPCARRY \
asm(                                                                                               \
"movq %1,%%rax                \n\t"                                                                \
"movq %2,%%rbx                \n\t"                                                                \
"addq  %%rax,%6               \n\t"                                                                \
"adcq  %%rbx,%7               \n\t"                                                                \
"adcq $0,%8                   \n\t"                                                                \
:"=g"(_c[OFF0]), "=g"(_c[OFF1]), "=g"(_c[OFF2]):"0"(_c[OFF0]), "1"(_c[OFF1]), "2"(_c[OFF2]),       \
                                                "m"(_c[OFF0+1]), "m"(_c[OFF1+1]), "m"(_c[OFF2+1])  \
: "%rax", "%rbx", "%cc");

#elif defined(TFM_SSE2)

/* SSE2 code */

#define MONT_START \
asm("movd %0,%%mm2"::"g"(mp));

#define MONT_FINI \
asm("emms");

#define LOOP_START \
asm(\
"movd %0,%%mm1                \n\t" \
"pmuludq %%mm2,%%mm1          \n\t" \
:: "g"(c[x]));

#define INNERMUL \
asm(                                                                                          \
"movd %6,%%mm0                \n\t"                                                           \
"pmuludq %%mm1,%%mm0          \n\t"                                                           \
"movd %%mm0,%%eax             \n\t"                                                           \
"psrlq $32, %%mm0             \n\t"                                                           \
"addl %%eax,%0                \n\t"                                                           \
"movd %%mm0,%%eax             \n\t"                                                           \
"adcl %%eax,%1                \n\t"                                                           \
"adcl $0,%2                   \n\t"                                                           \
:"=g"(_c[OFF0]), "=g"(_c[OFF1]), "=g"(_c[OFF2]):"0"(_c[OFF0]), "1"(_c[OFF1]), "2"(_c[OFF2]),  \
                                                "g"(*tmpm++)                                  \
                                               : "%eax", "%cc");

#define PROPCARRY \
asm(                                                                                               \
"movl %1,%%eax                \n\t"                                                                \
"addl  %%eax,%6               \n\t"                                                                \
"movl %2,%%eax                \n\t"                                                                \
"adcl  %%eax,%7               \n\t"                                                                \
"adcl $0,%8                   \n\t"                                                                \
:"=g"(_c[OFF0]), "=g"(_c[OFF1]), "=g"(_c[OFF2]):"0"(_c[OFF0]), "1"(_c[OFF1]), "2"(_c[OFF2]),       \
                                                "g"(_c[OFF0+1]), "g"(_c[OFF1+1]), "g"(_c[OFF2+1])  \
: "%eax", "%cc");

#elif defined(TFM_ARM)

/* ISO C code */
#define MONT_START 

#define MONT_FINI

#define LOOP_START \
   mu = c[x] * mp;

/* NOTE: later write it using two regs instead of three for _c + ... */
#define INNERMUL \
asm(                                             \
"UMULL r0,r1,%0,%1                \n\t"          \
"LDR   r2,[%2]                    \n\t"          \
"ADDS  r2,r2,r0                   \n\t"          \
"STR   r2,[%2]                    \n\t"          \
"LDR   r2,[%3]                    \n\t"          \
"ADCS  r2,r2,r1                   \n\t"          \
"STR   r2,[%3]                    \n\t"          \
"LDR   r2,[%4]                    \n\t"          \
"ADC   r2,r2,#0                   \n\t"          \
"STR   r2,[%4]                    \n\t"          \
::"r"(mu),"r"(*tmpm++),"r"(_c + OFF0),"r"(_c + OFF1),"r"(_c + OFF2):"r0", "r1", "r2", "%cc");

#define PROPCARRY \
asm(                                             \
"LDR   r0,[%1]                    \n\t"          \
"LDR   r1,[%0,#4]                 \n\t"          \
"ADDS  r0,r0,r1                   \n\t"          \
"STR   r0,[%0,#4]                 \n\t"          \
"LDR   r0,[%2]                    \n\t"          \
"LDR   r1,[%1,#4]                 \n\t"          \
"ADCS  r0,r0,r1                   \n\t"          \
"STR   r0,[%1,#4]                 \n\t"          \
"LDR   r0,[%2,#4]                 \n\t"          \
"ADC   r0,r0,#0                   \n\t"          \
"STR   r0,[%2,#4]                 \n\t"          \
::"r"(_c + OFF0),"r"(_c + OFF1),"r"(_c + OFF2):"r0", "r1", "%cc");

#else

/* ISO C code */
#define MONT_START 

#define MONT_FINI

#define LOOP_START \
   mu = c[x] * mp;

#define INNERMUL \
   do { fp_word t;                                                           \
   t = (fp_word)_c[OFF0] + ((fp_word)mu) * ((fp_word)*tmpm++); _c[OFF0] = t; \
   t = (fp_word)_c[OFF1] + (t >> DIGIT_BIT);                   _c[OFF1] = t; \
   _c[OFF2] += (t >> DIGIT_BIT);                                             \
   } while (0);

#define PROPCARRY \
   do { fp_word t;                                                           \
   t = (fp_word)_c[OFF0+1] + (fp_word)_c[OFF1];                    _c[OFF0+1] = t; \
   t = (fp_word)_c[OFF1+1] + (t >> DIGIT_BIT) + (fp_word)_c[OFF2]; _c[OFF1+1] = t; \
   _c[OFF2+1] += (t >> DIGIT_BIT);                                           \
   } while (0);

#endif


#define OFF0  (0)
#define OFF1  (FP_SIZE)
#define OFF2  (FP_SIZE+FP_SIZE)

/* computes x/R == x (mod N) via Montgomery Reduction */
void fp_montgomery_reduce(fp_int *a, fp_int *m, fp_digit mp)
{
   fp_digit c[3*FP_SIZE], *_c, *tmpm, mu;
   int      oldused, x, y, pa;

   /* now zero the buff */
   pa = m->used;
   memset(c, 0, sizeof(c));

   /* copy the input */
   oldused = a->used;
   for (x = 0; x < oldused; x++) {
       c[x] = a->dp[x];
   }

   MONT_START;

   /* now let's get bizz-sy! */
   for (x = 0; x < pa; x++) {
       /* get Mu for this round */
       LOOP_START;

       /* our friendly neighbourhood alias */
       _c   = c + x;
       tmpm = m->dp;

       for (y = 0; y < pa; y++) {
          INNERMUL;
          ++_c;
       }
       /* send carry up man... */
       _c = c + x;
       PROPCARRY;
  }         

  /* fix the rest of the carries */
  _c = c + pa;
  for (x = pa; x < pa * 2 + 2; x++) {
     PROPCARRY;
     ++_c;
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
