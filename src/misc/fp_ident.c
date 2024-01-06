/* TomsFastMath, a fast ISO C bignum library. -- Tom St Denis */
/* SPDX-License-Identifier: Unlicense */
#include <tfm_private.h>

#ifndef GIT_VERSION
#define GIT_VERSION TFM_VERSION_S
#endif

#define dnstrcon(D, N, S) dnmemcpy((D), (N), S, sizeof(S)-1)
static signed long dnmemcpy(char **d, size_t *n, const char *s, size_t len) {
   if (len >= *n) return -1;
   memcpy(*d, s, len);
   *n -= len;
   *d += len;
   **d = '\0';
   return len;
}

/* log(2)/log(10) ~= 0.30102999... ~= 30103 / 100000
 * we need to add one byte because this rounds to zero, and one for sign
 * these provide exact answer for integers up to 4720 bytes wide... */
#define U_DIGITS(T) (1 + ((sizeof(T) * 8UL)     * 30103UL) / 100000UL)
#define S_DIGITS(T) (2 + ((sizeof(T) * 8UL - 1) * 30103UL) / 100000UL)

static signed long dnstrul(char **d, size_t *n, unsigned long value) {
   char digits[U_DIGITS(unsigned long)+1]; /* fit digits plus null byte */
   char *digit = digits + (sizeof(digits) - 1);
   size_t len = 0;
   *digit = '\0';
   do {
      *--digit = '0' + (value % 10);
      value /= 10;
      ++len;
      if (digit < digits) return -1;
   } while (value);
   if (len >= *n) return -1;
   return dnmemcpy(d, n, digit, len);
}

const char *fp_ident(void)
{
   static char buf[1024];
   char *d = buf;
   size_t n = sizeof(buf);

   dnstrcon(&d, &n,
"TomsFastMath " GIT_VERSION "\n"
#if defined(TFM_IDENT_BUILD_DATE)
"Built on " __DATE__ " at " __TIME__ "\n"
#endif
"\n"
"Sizeofs\n"
"\tfp_digit = "
   );
   dnstrul(&d, &n, sizeof(fp_digit));
   dnstrcon(&d, &n,
"\n"
"\tfp_word  = "
   );
   dnstrul(&d, &n, sizeof(fp_word));
   dnstrcon(&d, &n,
"\n\n"
"FP_MAX_SIZE = "
   );
   dnstrul(&d, &n, FP_MAX_SIZE);
   dnstrcon(&d, &n,
"\n\n"
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
"\n"
   );

   if (sizeof(fp_digit) == sizeof(fp_word)) {
      dnstrcon(&d, &n,
         "WARNING: sizeof(fp_digit) == sizeof(fp_word),"
         " this build is likely to not work properly.\n"
      );
   }

   memset(d, 0, n);
   return buf;
}

#ifdef STANDALONE

int main(void)
{
   printf("%s\n", fp_ident());
   return 0;
}

#endif
