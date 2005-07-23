/* generate montgomery reductions for m->used = 1...16 */

#include <stdio.h>

int main(void)
{
   int N;
   
   for (N = 1; N <= 16; N++) {
       
printf("void fp_montgomery_reduce_%d(fp_int *a, fp_int *m, fp_digit mp)\n", N);
printf(
"{\n"
"   fp_digit c[3*FP_SIZE], *_c, *tmpm, mu;\n"
"   int      oldused, x, y;\n"
"\n"
"   /* now zero the buff */\n"
"   memset(c, 0, sizeof(c));\n"
"\n"
"   /* copy the input */\n"
"   oldused = a->used;\n"
"   for (x = 0; x < oldused; x++) {\n"
"       c[x] = a->dp[x];\n"
"   }\n"
"\n"
"   MONT_START;\n"
"\n"
"   /* now let's get bizz-sy! */\n"
"   for (x = 0; x < %d; x++) {\n"
"       /* get Mu for this round */\n"
"       LOOP_START;\n"
"\n"
"       /* our friendly neighbourhood alias */\n"
"       _c   = c + x;\n"
"       tmpm = m->dp;\n"
"\n"
"       for (y = 0; y < %d; y++) {\n"
"          INNERMUL;\n"
"          ++_c;\n"
"       }\n"
"       /* send carry up man... */\n"
"       _c = c + x;\n"
"       PROPCARRY;\n"
"   }         \n"
"\n"
"  /* fix the rest of the carries */\n"
"  _c = c + %d;\n"
"  for (x = %d; x < %d * 2 + 2; x++) {\n"
"     PROPCARRY;\n"
"     ++_c;\n"
"  }\n"
"\n"
"  /* now copy out */\n"
"  _c   = c + %d;\n"
"  tmpm = a->dp;\n"
"  for (x = 0; x < %d+1; x++) {\n"
"     *tmpm++ = *_c++;\n"
"  }\n"
"\n"
"  for (; x < oldused; x++)   {\n"
"     *tmpm++ = 0;\n"
"  }\n"
"\n"
"  MONT_FINI;\n"
"\n"
"  a->used = %d+1;\n"
"  fp_clamp(a);\n"
"\n"  
"  /* if A >= m then A = A - m */\n"
"  if (fp_cmp_mag (a, m) != FP_LT) {\n"
"    s_fp_sub (a, m, a);\n"
"  }\n"
"}\n", N,N,N,N,N,N,N,N);
}

return 0;
}




