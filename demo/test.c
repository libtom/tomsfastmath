/* TFM demo program */
#include <tfm.h>
#include <time.h>
#include <unistd.h>


#ifndef TFM_DEMO_TEST_VS_MTEST
#define TFM_DEMO_TEST_VS_MTEST 1
#endif

void draw(fp_int *a)
{
  int x;
  printf("%d, %d, ", a->used, a->sign);
  for (x = a->used - 1; x >= 0; x--) {
#if SIZEOF_FP_DIGIT == 4
      printf("%08lx ", a->dp[x]);
#else
      printf("%016llx ", a->dp[x]);
#endif
  }
  printf("\n");
}

int myrng(unsigned char *dst, int len, void *dat)
{
   int x;
   (void)dat;
   for (x = 0; x < len; x++) dst[x] = rand() & 0xFF;
   return len;
}

   char cmd[4096], buf[4096];

int main(void)
{
  fp_int a,b,c,d,e,f;
  unsigned long ix;
#if TFM_DEMO_TEST_VS_MTEST
  unsigned long expt_n, add_n, sub_n, mul_n, div_n, sqr_n, mul2d_n, div2d_n, gcd_n, lcm_n, inv_n,
                 div2_n, mul2_n, add_d_n, sub_d_n, mul_d_n, cnt, rr;
#else
  fp_digit fp;
  int n, err;
#endif

  srand(time(NULL));
  printf("TFM Ident string:\n%s\n\n", fp_ident());
  fp_zero(&b); fp_zero(&c); fp_zero(&d); fp_zero(&e); fp_zero(&f);
  fp_zero(&a);

#if TFM_DEMO_TEST_VS_MTEST == 0

  draw(&a);

  /* test set and simple shifts */
  printf("Testing mul/div 2\n");
  fp_set(&a, 1); draw(&a);
  for (n = 0; n <= DIGIT_BIT; n++) {
      fp_mul_2(&a, &a); printf("(%d) ", fp_count_bits(&a));
      draw(&a);

  }
  for (n = 0; n <= (DIGIT_BIT + 1); n++) {
      fp_div_2(&a, &a);
      draw(&a);
  }
  fp_set(&a, 1);

  /* test lshd/rshd */
  printf("testing lshd/rshd\n");
  fp_lshd(&a, 3); draw(&a);
  fp_rshd(&a, 3); draw(&a);

  /* test more complicated shifts */
  printf("Testing mul/div 2d\n");
  fp_mul_2d(&a, DIGIT_BIT/2, &a); draw(&a);
  fp_div_2d(&a, DIGIT_BIT/2, &a, NULL); draw(&a);

  fp_mul_2d(&a, DIGIT_BIT + DIGIT_BIT/2, &a); draw(&a);
  fp_div_2d(&a, DIGIT_BIT + DIGIT_BIT/2, &a, NULL); draw(&a);

  /* test neg/abs  */
  printf("testing neg/abs\n");
  fp_neg(&a, &a); draw(&a);
  fp_neg(&a, &a); draw(&a);
  fp_neg(&a, &a); draw(&a);
  fp_abs(&a, &a); draw(&a);

  /* test comparisons */
  fp_set(&b, 3);
  fp_set(&c, 4); fp_neg(&c, &c);
  fp_set(&d, 1);
  printf("Testing compares\n%d, %d, %d, %d\n", fp_cmp(&a, &b), fp_cmp(&a, &c), fp_cmp(&a, &d), fp_cmp(&b, &c));

  /* test add/sub */
  printf("Testing add/sub \n");
  fp_set(&a, ((fp_digit)1)<<(DIGIT_BIT-1)); draw(&a);
  fp_set(&b, ((fp_digit)1)<<(DIGIT_BIT-2));
  fp_add(&a, &b, &a); draw(&a);
  fp_add(&a, &b, &a); draw(&a);
  fp_add(&a, &b, &a); draw(&a);
  printf("sub...\n");
  printf("cmp returns: %d, ", fp_cmp(&a, &b)); fp_sub(&a, &b, &a); draw(&a);
  printf("cmp returns: %d, ", fp_cmp(&a, &b)); fp_sub(&a, &b, &a); draw(&a);
  printf("cmp returns: %d, ", fp_cmp(&a, &b)); fp_sub(&a, &b, &a); draw(&a);
  printf("cmp returns: %d, ", fp_cmp(&a, &b)); fp_sub(&a, &b, &a); draw(&a);
  printf("cmp returns: %d, ", fp_cmp(&a, &b)); fp_sub(&a, &b, &a); draw(&a);
  printf("cmp returns: %d, ", fp_cmp(&a, &b)); fp_sub(&a, &b, &a); draw(&a);
  fp_read_radix(&a, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001", 16); draw(&a);
  fp_sub_d(&a, 3, &b); draw(&b);
  fp_read_radix(&a, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFE", 16);
  printf("cmp returns: %d, ", fp_cmp(&a, &b)); fp_sub(&a, &b, &a); draw(&a);

  /* test mul_d */
  printf("Testing mul_d and div_d\n");
  fp_set(&a, 1);
  fp_mul_d(&a, ((fp_digit)1)<<(DIGIT_BIT/2), &a); draw(&a);
  fp_mul_d(&a, ((fp_digit)1)<<(DIGIT_BIT/2), &a); draw(&a);
  fp_mul_d(&a, ((fp_digit)1)<<(DIGIT_BIT/2), &a); draw(&a);
  printf("div_d\n");
  fp_div_d(&a, ((fp_digit)1)<<(DIGIT_BIT/2), &a, NULL); draw(&a);
  fp_div_d(&a, ((fp_digit)1)<<(DIGIT_BIT/2), &a, NULL); draw(&a);
  fp_div_d(&a, ((fp_digit)1)<<(DIGIT_BIT/2), &a, NULL); draw(&a);

  /* testing read radix */
  printf("Testing read_radix\n");
  fp_read_radix(&a, "123456789012345678901234567890", 16); draw(&a);

  /* test mont */
  printf("Montgomery test #1\n");
  fp_set(&a, 0x1234567ULL);
  fp_montgomery_setup(&a, &fp);
  fp_montgomery_calc_normalization(&b, &a);

  fp_read_radix(&d, "123456789123", 16);
  for (n = 0; n < 1000000; n++) {
      fp_add_d(&d, 1, &d); fp_sqrmod(&d, &a, &d);
      fp_mul(&d, &b, &c);
      fp_montgomery_reduce(&c, &a, fp);
      if (fp_cmp(&c, &d) != FP_EQ) {
         printf("Failed mont %d\n", n);
         draw(&a);
         draw(&d);
         draw(&c);
         return EXIT_FAILURE;
      }
  }
  printf("Passed.\n");

  printf("Montgomery test #2\n");
  fp_set(&a, 0x1234567ULL);
  fp_lshd(&a, 4);
  fp_add_d(&a, 1, &a);
  fp_montgomery_setup(&a, &fp);
  fp_montgomery_calc_normalization(&b, &a);

  fp_read_radix(&d, "123456789123", 16);
  for (n = 0; n < 1000000; n++) {
      fp_add_d(&d, 1, &d); fp_sqrmod(&d, &a, &d);
      fp_mul(&d, &b, &c);
      fp_montgomery_reduce(&c, &a, fp);
      if (fp_cmp(&c, &d) != FP_EQ) {
         printf("Failed mont %d\n", n);
         draw(&a);
         draw(&d);
         draw(&c);
         return EXIT_FAILURE;
      }
  }
  printf("Passed.\n");

   /* test for size */
   for (ix = 8*DIGIT_BIT; ix < 10*DIGIT_BIT; ix++) {
       printf("Testing (not safe-prime): %9lu bits    \r", ix); fflush(stdout);
       err = fp_prime_random_ex(&a, 8, ix, (rand()&1)?TFM_PRIME_2MSB_OFF:TFM_PRIME_2MSB_ON, myrng, NULL);
       if (err != FP_OKAY) {
          printf("failed with err code %d\n", err);
          return EXIT_FAILURE;
       }
       if ((unsigned long)fp_count_bits(&a) != ix) {
          printf("Prime is %d not %lu bits!!!\n", fp_count_bits(&a), ix);
          return EXIT_FAILURE;
       }
   }
   printf("\n\n");

return 0;

#else

  fp_zero(&b); fp_zero(&c); fp_zero(&d); fp_zero(&e); fp_zero(&f); fp_zero(&a);


   div2_n = mul2_n = inv_n = expt_n = lcm_n = gcd_n = add_n =
   sub_n = mul_n = div_n = sqr_n = mul2d_n = div2d_n = cnt = add_d_n = sub_d_n= mul_d_n = 0;

   for (;;) {
       printf("%4lu/%4lu/%4lu/%4lu/%4lu/%4lu/%4lu/%4lu/%4lu/%4lu/%4lu/%4lu/%4lu/%4lu/%4lu/%4lu ", add_n, sub_n, mul_n, div_n, sqr_n, mul2d_n, div2d_n, gcd_n, lcm_n, expt_n, inv_n, div2_n, mul2_n, add_d_n, sub_d_n, mul_d_n);
       fgets(cmd, 4095, stdin);
       cmd[strlen(cmd)-1] = 0;
       printf("%s  ]\r",cmd); fflush(stdout);
       if (!strcmp(cmd, "mul2d")) { ++mul2d_n;
          fgets(buf, 4095, stdin); fp_read_radix(&a, buf, 64);
          fgets(buf, 4095, stdin); sscanf(buf, "%lu", &rr);
          fgets(buf, 4095, stdin); fp_read_radix(&b, buf, 64);

          fp_mul_2d(&a, rr, &a);
          a.sign = b.sign;
          if (fp_cmp(&a, &b) != FP_EQ) {
             printf("\nmul2d failed, rr == %lu\n",rr);
             draw(&a);
             draw(&b);
             return 0;
          }
       } else if (!strcmp(cmd, "div2d")) { ++div2d_n;
          fgets(buf, 4095, stdin); fp_read_radix(&a, buf, 64);
          fgets(buf, 4095, stdin); sscanf(buf, "%lu", &rr);
          fgets(buf, 4095, stdin); fp_read_radix(&b, buf, 64);

          fp_div_2d(&a, rr, &a, &e);
          a.sign = b.sign;
          if (a.used == b.used && a.used == 0) { a.sign = b.sign = FP_ZPOS; }
          if (fp_cmp(&a, &b) != FP_EQ) {
             printf("\ndiv2d failed, rr == %lu\n",rr);
             draw(&a);
             draw(&b);
             return 0;
          }
       } else if (!strcmp(cmd, "add")) { ++add_n;
          fgets(buf, 4095, stdin);  fp_read_radix(&a, buf, 64);
          fgets(buf, 4095, stdin);  fp_read_radix(&b, buf, 64);
          fgets(buf, 4095, stdin);  fp_read_radix(&c, buf, 64);
          fp_copy(&a, &d);
          fp_add(&d, &b, &d);
          if (fp_cmp(&c, &d) != FP_EQ) {
             printf("\nadd %lu failure!\n", add_n);
draw(&a);draw(&b);draw(&c);draw(&d);
             return 0;
          }

          /* test the sign/unsigned storage functions */

          rr = fp_signed_bin_size(&c);
          fp_to_signed_bin(&c, (unsigned char *)cmd);
          memset(cmd+rr, rand()&255, sizeof(cmd)-rr);
          fp_read_signed_bin(&d, (unsigned char *)cmd, rr);
          if (fp_cmp(&c, &d) != FP_EQ) {
             printf("f\np_signed_bin failure!\n");
             draw(&c);
             draw(&d);
             return 0;
          }

          rr = fp_unsigned_bin_size(&c);
          fp_to_unsigned_bin(&c, (unsigned char *)cmd);
          memset(cmd+rr, rand()&255, sizeof(cmd)-rr);
          fp_read_unsigned_bin(&d, (unsigned char *)cmd, rr);
          if (fp_cmp_mag(&c, &d) != FP_EQ) {
             printf("\nfp_unsigned_bin failure!\n");
             draw(&c);
             draw(&d);
             return 0;
          }

       } else if (!strcmp(cmd, "sub")) { ++sub_n;
          fgets(buf, 4095, stdin);  fp_read_radix(&a, buf, 64);
          fgets(buf, 4095, stdin);  fp_read_radix(&b, buf, 64);
          fgets(buf, 4095, stdin);  fp_read_radix(&c, buf, 64);
          fp_copy(&a, &d);
          fp_sub(&d, &b, &d);
          if (fp_cmp(&c, &d) != FP_EQ) {
             printf("\nsub %lu failure!\n", sub_n);
draw(&a);draw(&b);draw(&c);draw(&d);
             return 0;
          }
       } else if (!strcmp(cmd, "mul")) { ++mul_n;
          fgets(buf, 4095, stdin);  fp_read_radix(&a, buf, 64);
          fgets(buf, 4095, stdin);  fp_read_radix(&b, buf, 64);
          fgets(buf, 4095, stdin);  fp_read_radix(&c, buf, 64);
//continue;
          fp_copy(&a, &d);
          fp_mul(&d, &b, &d);
          if (fp_cmp(&c, &d) != FP_EQ) {
             printf("\nmul %lu failure!\n", mul_n);
draw(&a);draw(&b);draw(&c);draw(&d);
             return 0;
          }
       } else if (!strcmp(cmd, "div")) { ++div_n;
          fgets(buf, 4095, stdin); fp_read_radix(&a, buf, 64);
          fgets(buf, 4095, stdin); fp_read_radix(&b, buf, 64);
          fgets(buf, 4095, stdin); fp_read_radix(&c, buf, 64);
          fgets(buf, 4095, stdin); fp_read_radix(&d, buf, 64);
// continue;
          fp_div(&a, &b, &e, &f);
          if (fp_cmp(&c, &e) != FP_EQ || fp_cmp(&d, &f) != FP_EQ) {
             printf("\ndiv %lu failure!\n", div_n);
draw(&a);draw(&b);draw(&c);draw(&d); draw(&e); draw(&f);
             return 0;
          }

       } else if (!strcmp(cmd, "sqr")) { ++sqr_n;
          fgets(buf, 4095, stdin);  fp_read_radix(&a, buf, 64);
          fgets(buf, 4095, stdin);  fp_read_radix(&b, buf, 64);
// continue;
          fp_copy(&a, &c);
          fp_sqr(&c, &c);
          if (fp_cmp(&b, &c) != FP_EQ) {
             printf("\nsqr %lu failure!\n", sqr_n);
draw(&a);draw(&b);draw(&c);
             return 0;
          }
       } else if (!strcmp(cmd, "gcd")) { ++gcd_n;
          fgets(buf, 4095, stdin);  fp_read_radix(&a, buf, 64);
          fgets(buf, 4095, stdin);  fp_read_radix(&b, buf, 64);
          fgets(buf, 4095, stdin);  fp_read_radix(&c, buf, 64);
// continue;
          fp_copy(&a, &d);
          fp_gcd(&d, &b, &d);
          d.sign = c.sign;
          if (fp_cmp(&c, &d) != FP_EQ) {
             printf("\ngcd %lu failure!\n", gcd_n);
draw(&a);draw(&b);draw(&c);draw(&d);
             return 0;
          }
       } else if (!strcmp(cmd, "lcm")) { ++lcm_n;
             fgets(buf, 4095, stdin);  fp_read_radix(&a, buf, 64);
             fgets(buf, 4095, stdin);  fp_read_radix(&b, buf, 64);
             fgets(buf, 4095, stdin);  fp_read_radix(&c, buf, 64);
//continue;
             fp_copy(&a, &d);
             fp_lcm(&d, &b, &d);
             d.sign = c.sign;
             if (fp_cmp(&c, &d) != FP_EQ) {
                printf("\nlcm %lu failure!\n", lcm_n);
   draw(&a);draw(&b);draw(&c);draw(&d);
                return 0;
             }
       } else if (!strcmp(cmd, "expt")) { ++expt_n;
             fgets(buf, 4095, stdin);  fp_read_radix(&a, buf, 64);
             fgets(buf, 4095, stdin);  fp_read_radix(&b, buf, 64);
             fgets(buf, 4095, stdin);  fp_read_radix(&c, buf, 64);
             fgets(buf, 4095, stdin);  fp_read_radix(&d, buf, 64);
// continue;
             fp_copy(&a, &e);
             fp_exptmod(&e, &b, &c, &e);
             if (fp_cmp(&d, &e) != FP_EQ) {
                printf("\nexpt %lu failure!\n", expt_n);
   draw(&a);draw(&b);draw(&c);draw(&d); draw(&e);
                return 0;
             }
       } else if (!strcmp(cmd, "invmod")) { ++inv_n;
             fgets(buf, 4095, stdin);  fp_read_radix(&a, buf, 64);
             fgets(buf, 4095, stdin);  fp_read_radix(&b, buf, 64);
             fgets(buf, 4095, stdin);  fp_read_radix(&c, buf, 64);
//continue;
             fp_invmod(&a, &b, &d);
#if 1
             fp_mulmod(&d,&a,&b,&e);
             if (fp_cmp_d(&e, 1) != FP_EQ) {
#else
             if (fp_cmp(&d, &c) != FP_EQ) {
#endif
                printf("\ninv [wrong value from MPI?!] failure\n");
                draw(&a);draw(&b);draw(&c);draw(&d);
                return 0;
             }

       } else if (!strcmp(cmd, "div2")) { ++div2_n;
             fgets(buf, 4095, stdin);  fp_read_radix(&a, buf, 64);
             fgets(buf, 4095, stdin);  fp_read_radix(&b, buf, 64);
             fp_div_2(&a, &c);
             if (fp_cmp(&c, &b) != FP_EQ) {
                 printf("\ndiv_2 %lu failure\n", div2_n);
                 draw(&a);
                 draw(&b);
                 draw(&c);
                 return 0;
             }
       } else if (!strcmp(cmd, "mul2")) { ++mul2_n;
             fgets(buf, 4095, stdin);  fp_read_radix(&a, buf, 64);
             fgets(buf, 4095, stdin);  fp_read_radix(&b, buf, 64);
             fp_mul_2(&a, &c);
             if (fp_cmp(&c, &b) != FP_EQ) {
                 printf("\nmul_2 %lu failure\n", mul2_n);
                 draw(&a);
                 draw(&b);
                 draw(&c);
                 return 0;
             }
       } else if (!strcmp(cmd, "add_d")) { ++add_d_n;
              fgets(buf, 4095, stdin); fp_read_radix(&a, buf, 64);
              fgets(buf, 4095, stdin); sscanf(buf, "%lu", &ix);
              fgets(buf, 4095, stdin); fp_read_radix(&b, buf, 64);
              fp_add_d(&a, ix, &c);
              if (fp_cmp(&b, &c) != FP_EQ) {
                 printf("\nadd_d %lu failure\n", add_d_n);
                 draw(&a);
                 draw(&b);
                 draw(&c);
                 printf("d == %lu\n", ix);
                 return 0;
              }
       } else if (!strcmp(cmd, "sub_d")) { ++sub_d_n;
              fgets(buf, 4095, stdin); fp_read_radix(&a, buf, 64);
              fgets(buf, 4095, stdin); sscanf(buf, "%lu", &ix);
              fgets(buf, 4095, stdin); fp_read_radix(&b, buf, 64);
              fp_sub_d(&a, ix, &c);
              if (fp_cmp(&b, &c) != FP_EQ) {
                 printf("\nsub_d %lu failure\n", sub_d_n);
                 draw(&a);
                 draw(&b);
                 draw(&c);
                 printf("d == %lu\n", ix);
                 return 0;
              }
       } else if (!strcmp(cmd, "mul_d")) { ++mul_d_n;
              fgets(buf, 4095, stdin); fp_read_radix(&a, buf, 64);
              fgets(buf, 4095, stdin); sscanf(buf, "%lu", &ix);
              fgets(buf, 4095, stdin); fp_read_radix(&b, buf, 64);
              fp_mul_d(&a, ix, &c);
              if (fp_cmp(&b, &c) != FP_EQ) {
                 printf("\nmul_d %lu failure\n", mul_d_n);
                 draw(&a);
                 draw(&b);
                 draw(&c);
                 printf("d == %lu\n", ix);
                 return 0;
              }
       } else if (!strcmp(cmd, "exit")) {
              printf("\nokay, exiting now\n");
              return 0;
       }

   }
#endif
}



/* $Source$ */
/* $Revision$ */
/* $Date$ */
