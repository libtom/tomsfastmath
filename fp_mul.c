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

/* c = a * b */
void fp_mul(fp_int *A, fp_int *B, fp_int *C)
{
    int    r, y, yy, s;
    fp_int ac, bd, comp, amb, cmd, t1, t2;

     y  = MAX(A->used, B->used);
     yy = MIN(A->used, B->used);
     if (yy <= 8 || y <= 64) {

    /* pick a comba (unrolled 4/8/16/32 x or rolled) based on the size
       of the largest input.  We also want to avoid doing excess mults if the 
       inputs are not close to the next power of two.  That is, for example,
       if say y=17 then we would do (32-17)^2 = 225 unneeded multiplications 
    */
        if (y <= 4) {
           fp_mul_comba4(A,B,C);
        } else if (y <= 8) {
           fp_mul_comba8(A,B,C);
#if defined(TFM_LARGE)
        } else if (y <= 16 && y >= 10) {
           fp_mul_comba16(A,B,C);
#endif
#if defined(TFM_HUGE)
        } else if (y <= 32 && y >= 24) {
           fp_mul_comba32(A,B,C);
#endif
        } else {
           fp_mul_comba(A,B,C);
        }
    } else {
        /* do the karatsuba action 

           if A = ab and B = cd for ||a|| = r we need to solve 

           ac*r^2 + (-(a-b)(c-d) + ac + bd)*r + bd

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
        t1.sign = A->sign;

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
        t2.sign = B->sign;

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
            t1.dp[s]   = 0; 
            t2.dp[s] = 0; 
        }
        t1.used = r;
        t2.used = r;
        fp_clamp(&t1);
        fp_clamp(&t2);
        
        fp_sub(&amb, &t1, &amb); fp_sub(&cmd, &t2, &cmd);
        fp_zero(&bd);
        fp_mul(&t1, &t2, &bd);

        /* now get the (a-b)(c-d) term */
        fp_zero(&comp);
        fp_mul(&amb, &cmd, &comp);

        /* now solve the system, do the middle term first */
        comp.sign ^= 1;
        fp_add(&comp, &ac, &comp);
        fp_add(&comp, &bd, &comp);
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

