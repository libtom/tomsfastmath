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

/* b = a*a  */
void fp_sqr(fp_int *A, fp_int *B)
{
    int    r, y, s;
    fp_int aa, bb, comp, amb, t1;

    y = A->used;
    if (y <= 64) { 
        if (y <= 4) {
           fp_sqr_comba4(A,B);
        } else if (y <= 8) {
           fp_sqr_comba8(A,B);
#if defined(TFM_LARGE)
        } else if (y <= 16 && y >= 12) {
           fp_sqr_comba16(A,B);
#endif
#if defined(TFM_HUGE)
        } else if (y <= 32 && y >= 20) {
           fp_sqr_comba32(A,B);
        } else if (y <= 64 && y >= 48) {
           fp_sqr_comba64(A,B);
#endif
        } else {
           fp_sqr_comba(A, B);
        }
       
    } else {
        /* do the karatsuba action 

           if A = ab ||a|| = r we need to solve 

           a^2*r^2 + (-(a-b)^2 + a^2 + b^2)*r + b^2

           So we solve for the three products then we form the final result with careful shifting 
           and addition.

Obvious points of optimization

- "ac" parts can be memcpy'ed with an offset [all you have to do is zero upto the next 8 digits]
- Similarly the "bd" parts can be memcpy'ed and zeroed to 8
- 

        */
        /* get our value of r */
        r = y >> 1;

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
        fp_copy(&t1, &amb); 
        fp_zero(&aa);
        fp_sqr(&t1, &aa);

        /* now solve for bd */
//        fp_mod_2d(A, r * DIGIT_BIT, &t1);
        for (s = 0; s < r; s++) {
            t1.dp[s] = A->dp[s];
        }
        for (; s < FP_SIZE; s++) {
            t1.dp[s]   = 0; 
        }
        t1.used = r;
        fp_clamp(&t1);

        fp_sub(&amb, &t1, &amb); 
        fp_zero(&bb);
        fp_sqr(&t1, &bb);

        /* now get the (a-b) term */
        fp_zero(&comp);
        fp_sqr(&amb, &comp);

        /* now solve the system, do the middle term first */
        comp.sign ^= 1;
        fp_add(&comp, &aa, &comp);
        fp_add(&comp, &bb, &comp);
        fp_lshd(&comp, r);
  
        /* leading term */
        fp_lshd(&aa, r+r);

        /* now sum them together */
        fp_zero(B);
        fp_add(&aa, &comp, B);
        fp_add(&bb, B, B);    
        B->sign = FP_ZPOS;
    }
}

