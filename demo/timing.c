/* TFM timing analysis */
#define _GNU_SOURCE
#include <tfm.h>
#include <time.h>
#include <unistd.h>

/* RDTSC from Scott Duplichan */
static ulong64 TIMFUNC(void)
{
#if defined __GNUC__
   #if defined(INTEL_CC)
 ulong64 a;
      asm ("rdtsc":"=A"(a));
      return a;
   #elif defined(__i386__) || defined(__x86_64__)
      /* version from http://www.mcs.anl.gov/~kazutomo/rdtsc.html
       * the old code always got a warning issued by gcc, clang did not complain...
       */
      unsigned hi, lo;
      __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
      return ((ulong64)lo)|( ((ulong64)hi)<<32);
   #elif defined(TFM_PPC32)
      unsigned long a, b;
      __asm__ __volatile__ ("mftbu %1 \nmftb %0\n":"=r"(a), "=r"(b));
      return (((ulong64)b) << 32ULL) | ((ulong64)a);
   #elif defined(TFM_AVR32)
 FILE *in;
      char buf[20];
 in = fopen("/sys/devices/system/cpu/cpu0/pccycles", "r");
 fgets(buf, 20, in);
 fclose(in);
 return strtoul(buf, NULL, 10);
   #else /* gcc-IA64 version */
      unsigned long result;
      __asm__ __volatile__("mov %0=ar.itc" : "=r"(result) :: "memory");
      while (__builtin_expect ((int) result == -1, 0))
      __asm__ __volatile__("mov %0=ar.itc" : "=r"(result) :: "memory");
      return result;
   #endif

// Microsoft and Intel Windows compilers
#elif defined _M_IX86
  __asm rdtsc
#elif defined _M_AMD64
  return __rdtsc ();
#elif defined _M_IA64
  #if defined __INTEL_COMPILER
    #include <ia64intrin.h>
  #endif
   return __getReg (3116);
#else
  #error need rdtsc function for this build
#endif
}

static ulong64 ticks;
static const char* p_str;

static void print_start(const char* s)
{
   p_str = s;
}

static void print_line(ulong64 b, ulong64 t)
{
   printf("%llu;%s;%llu;%llu\n", ticks, p_str, b, t);
}

static void Addition(ulong64 t1)
{
   fp_int a,b,c;
   ulong64 t2;
   unsigned long t, ix;
   for (t = 2; t <= FP_SIZE / 2; t += 2) {
      fp_zero(&a);
      fp_zero(&b);
      fp_zero(&c);
      for (ix = 0; ix < t; ix++) {
         a.dp[ix] = ix;
         b.dp[ix] = ix;
      }
      a.used = t;
      b.used = t;
      t2 = -1;
      for (ix = 0; ix < 25000; ++ix) {
         t1 = TIMFUNC();
         fp_add(&a, &b, &c);
         fp_add(&a, &b, &c);
         fp_add(&a, &b, &c);
         fp_add(&a, &b, &c);
         fp_add(&a, &b, &c);
         fp_add(&a, &b, &c);
         fp_add(&a, &b, &c);
         fp_add(&a, &b, &c);
         t2 = (TIMFUNC() - t1) >> 3;
         if (t1 < t2) {
            --ix;
            t2 = t1;
         }
      }
      print_line(t * DIGIT_BIT, t2);
   }
}

static void Multiplication(ulong64 t1)
{
   fp_int a,b,c;
   ulong64 t2;
   unsigned long t, ix;
   for (t = 2; t < FP_SIZE / 2; t += 2) {
      fp_zero(&a);
      fp_zero(&b);
      fp_zero(&c);
      for (ix = 0; ix < t; ix++) {
         a.dp[ix] = ix;
         b.dp[ix] = ix;
      }
      a.used = t;
      b.used = t;
      t2 = -1;
      for (ix = 0; ix < 100; ++ix) {
         t1 = TIMFUNC();
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         fp_mul(&a, &b, &c);
         t2 = (TIMFUNC() - t1) >> 7;
         if (t1 < t2) {
            --ix;
            t2 = t1;
         }
      }
      print_line(t * DIGIT_BIT, t2);
   }
}

static void Squaring(ulong64 t1)
{
   fp_int a,b;
   ulong64 t2;
   unsigned long t, ix;
   for (t = 2; t < FP_SIZE / 2; t += 2) {
      fp_zero(&a);
      fp_zero(&b);
      for (ix = 0; ix < t; ix++) {
         a.dp[ix] = ix;
      }
      a.used = t;
      t2 = -1;
      for (ix = 0; ix < 100; ++ix) {
         t1 = TIMFUNC();
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         fp_sqr(&a, &b);
         t2 = (TIMFUNC() - t1) >> 7;
         if (t1 < t2) {
            --ix;
            t2 = t1;
         }
      }
      print_line(t * DIGIT_BIT, t2);
   }
}

static void Invmod(ulong64 t1)
{
   fp_int a,b,c;
   ulong64 t2;
   unsigned long t, ix;
   for (t = 2; t < FP_SIZE / 2; t += 2) {
      fp_zero(&a);
      for (ix = 0; ix < t; ix++) {
         a.dp[ix] = ix | 1;
      }
      a.used = t;
      fp_zero(&b);
      for (ix = 0; ix < t; ix++) {
         b.dp[ix] = rand();
      }
      b.used = t;
      fp_clamp(&b);
      fp_zero(&c);
      t2 = -1;
      for (ix = 0; ix < 100; ++ix) {
         t1 = TIMFUNC();
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         fp_invmod(&b, &a, &c);
         t2 = (TIMFUNC() - t1) >> 6;
         if (t1 < t2) {
            --ix;
            t2 = t1;
         }
      }
      print_line(t * DIGIT_BIT, t2);
   }
}

static void Montgomery(ulong64 t1)
{
   fp_int a,b,c,d;
   ulong64 t2;
   fp_digit fp;
   unsigned long t, ix;
   for (t = 2; t <= (FP_SIZE / 2) - 4; t += 2) {
      //      printf("%5lu-bit: %9llu\n", t * DIGIT_BIT, t2);
      fp_zero(&a);
      for (ix = 0; ix < t; ix++) {
         a.dp[ix] = ix | 1;
      }
      a.used = t;

      fp_montgomery_setup(&a, &fp);
      fp_sub_d(&a, 3, &b);
      fp_sqr(&b, &b);
      fp_copy(&b, &c);
      fp_copy(&b, &d);

      t2 = -1;
      for (ix = 0; ix < 100; ++ix) {
         t1 = TIMFUNC();
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         fp_montgomery_reduce(&c, &a, fp);
         fp_montgomery_reduce(&d, &a, fp);
         t2 = (TIMFUNC() - t1) >> 6;
         fp_copy(&b, &c);
         fp_copy(&b, &d);
         if (t1 < t2) {
            --ix;
            t2 = t1;
         }
      }
      print_line(t * DIGIT_BIT, t2);
   }
}

static void Exptmod(ulong64 t1)
{
   fp_int a,b,c,d;
   ulong64 t2;
   unsigned long t, ix;
   for (t = 512 / DIGIT_BIT; t <= (FP_SIZE / 2) - 2; t += 256 / DIGIT_BIT) {
      fp_zero(&a);
      fp_zero(&b);
      fp_zero(&c);
      for (ix = 0; ix < t; ix++) {
         a.dp[ix] = ix + 1;
         b.dp[ix] = (fp_digit) rand() * (fp_digit) rand();
         c.dp[ix] = ix;
      }
      a.used = t;
      b.used = t;
      c.used = t;

      t2 = -1;
      for (ix = 0; ix < 500; ++ix) {
         t1 = TIMFUNC();
         fp_exptmod(&c, &b, &a, &d);
         fp_exptmod(&c, &b, &a, &d);
         t2 = (TIMFUNC() - t1) >> 1;
         fp_copy(&b, &c);
         fp_copy(&b, &d);
         if (t1 < t2) {
            t2 = t1;
            --ix;
         }
      }
      print_line(t * DIGIT_BIT, t2);
   }
}

#define FN(n) { n, #n }
struct {
   void (*fn)(ulong64 t1);
   const char* name;
} funcs[] = {
             FN(Addition),
             FN(Multiplication),
             FN(Squaring),
             FN(Invmod),
             FN(Montgomery),
             FN(Exptmod),
};

int main(int argc, char **argv)
{
   ulong64 t1;
   unsigned int t;
   char* arg = NULL;

   if (argc > 1) arg = argv[1];

   t1 = TIMFUNC();
   sleep(1);
   ticks = TIMFUNC() - t1;
   fprintf(stderr, "Ticks per second: %llu\n", ticks);

   printf("Ticks/sec;Algorithm;bits;time\n");
   /* do some timings... */
   for (t = 0; t < sizeof(funcs)/sizeof(funcs[0]); ++t) {
      if(arg != NULL && strcasestr(funcs[t].name, arg) == NULL) continue;
      print_start(funcs[t].name);
      funcs[t].fn(t1);
   }

   return 0;
}
