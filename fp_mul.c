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

/* c = a * b */
void fp_mul(fp_int *A, fp_int *B, fp_int *C)
{
    int    r, y, yy, s;
    fp_int ac, bd, comp, amb, cmd, t1, t2;

    /* call generic if we're out of range */
    if (A->used + B->used > FP_SIZE) {
       fp_mul_comba(A, B, C);
       return ;
    }

     y  = MAX(A->used, B->used);
     yy = MIN(A->used, B->used);
     if (yy <= 8 || y <= 64) {

    /* pick a comba (unrolled 4/8/16/32 x or rolled) based on the size
       of the largest input.  We also want to avoid doing excess mults if the 
       inputs are not close to the next power of two.  That is, for example,
       if say y=17 then we would do (32-17)^2 = 225 unneeded multiplications 
    */

#ifdef TFM_SMALL_SET
        if (y <= 16) {
           fp_mul_comba_small(A,B,C);
#elif defined(TFM_HUGE)
        if (0) { 1;
#endif
#if defined(TFM_MUL32)
        } else if (y <= 32) {
           fp_mul_comba32(A,B,C);
#endif
#if defined(TFM_MUL48)
        } else if (y <= 48) {
           fp_mul_comba48(A,B,C);
#endif
#if defined(TFM_MUL64)
        } else if (y <= 64) {
           fp_mul_comba64(A,B,C);
#endif
#if !defined(TFM_HUGE) && !defined(TFM_SMALL_SET)
        {
#else
        } else {
#endif
           fp_mul_comba(A,B,C);
        }
    } else {
        /* do the karatsuba action 

           if A = ab and B = cd for ||a|| = r we need to solve 

           ac*r^2 + ((a+b)(c+d) - (ac + bd))*r + bd

           So we solve for the three products then we form the final result with careful shifting 
           and addition.

Obvious points of optimization

- "ac" parts can be memcpy'ed with an offset [all you have to do is zero upto the next 8 digits]
- Similarly the "bd" parts can be memcpy'ed and zeroed to 8
- 

        */
        /* get our value of r */
        r = yy >> 1;

        /* now solve for ac */
//        fp_copy(A, &t1); fp_rshd(&t1, r); 
        for (s = 0; s < A->used - r; s++) {
            t1.dp[s] = A->dp[s+r];
        }
        for (; s < FP_SIZE; s++) {
            t1.dp[s] = 0; 
        }
        if (A->used >= r) {
           t1.used = A->used - r;
        } else {
           t1.used = 0;
        }
        t1.sign = 0;

//        fp_copy(B, &t2); fp_rshd(&t2, r); 
        for (s = 0; s < B->used - r; s++) {
            t2.dp[s] = B->dp[s+r];
        }
        for (; s < FP_SIZE; s++) {
            t2.dp[s] = 0; 
        }
        if (B->used >= r) {
           t2.used = B->used - r;
        } else {
           t2.used = 0;
        }
        t2.sign = 0;

        fp_copy(&t1, &amb); fp_copy(&t2, &cmd);
        fp_zero(&ac);
        fp_mul(&t1, &t2, &ac);

        /* now solve for bd */
//        fp_mod_2d(A, r * DIGIT_BIT, &t1);
//        fp_mod_2d(B, r * DIGIT_BIT, &t2);
        for (s = 0; s < r; s++) {
            t1.dp[s] = A->dp[s];
            t2.dp[s] = B->dp[s];
        }
        for (; s < FP_SIZE; s++) {
            t1.dp[s] = 0; 
            t2.dp[s] = 0; 
        }
        t1.used = r;
        t2.used = r;
        fp_clamp(&t1);
        fp_clamp(&t2);
        
        s_fp_add(&amb, &t1, &amb); s_fp_add(&cmd, &t2, &cmd);
        fp_zero(&bd);
        fp_mul(&t1, &t2, &bd);

        /* now get the (a+b)(c+d) term */
        fp_zero(&comp);
        fp_mul(&amb, &cmd, &comp);

        /* now solve the system, do the middle term first */
        s_fp_sub(&comp, &ac, &comp);
        s_fp_sub(&comp, &bd, &comp);
        fp_lshd(&comp, r);
  
        /* leading term */
        fp_lshd(&ac, r+r);

        /* now sum them together */
        s = A->sign ^ B->sign;
        fp_zero(C);
        fp_add(&ac, &comp, C);
        fp_add(&bd, C, C);    
        C->sign = C->used ? s : FP_ZPOS;
    }
}


/* $Source$ */
/* $Revision$ */
/* $Date$ */
