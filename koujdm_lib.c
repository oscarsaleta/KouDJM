#include "koujdm_lib.h"

#include <inttypes.h>

/****
 * GENERIC FUNCTIONS
 ****/
// binomial coefficient computed recursively
int64_t binomial(int64_t n, int64_t k) {
    if (n<k)
        return -1;
    if (k==0)
        return 1;
    if (k > n/2)
        return binomial(n,n-k);
    return n*binomial(n-1,k-1)/k;
}

// coefs. of the polynomial of degree 4
Quartic polynomial(long double complex alpha, Parameters prm) {
    Quartic p;
    long double sigma2 = prm.sigma*prm.sigma;
    p.a = 1;
    p.b = -prm.eta1 + prm.eta2 + 2.*prm.mu/sigma2;
    p.c = -prm.eta1*prm.eta2 + 2./sigma2*(-prm.mu*prm.eta1 - prm.lambda + prm.mu*prm.eta2 - alpha);
    p.d = 2./sigma2*(-prm.eta1*prm.eta2*prm.mu - prm.eta1*prm.lambda*(prm.p-1) - prm.eta2*prm.lambda*prm.p + alpha*(prm.eta1-prm.eta2));
    p.e = 2.*prm.eta1*prm.eta2*alpha/sigma2;

    return p;
}

// direct computation of the roots of the polynomial of degree 4
//http://math.stackexchange.com/a/786
QuarticSolutions quartic_solve(Quartic p) {
    long double complex p1 = 2*p.c*p.c*p.c-9*p.b*p.c*p.d+27*p.a*p.d*p.d+27*p.b*p.b*p.e-72*p.a*p.c*p.e;
    long double complex p2 = p1+csqrt(-4*cpow((p.c*p.c-3*p.b*p.d+12*p.a*p.e),3.)+p1*p1);
    long double complex p3 = (p.c*p.c-3*p.b*p.d+12*p.a*p.e)/(3*p.a*cpow(p2/2.,1/3.))+cpow(p2/2.,1/3.)/(3*p.a);
    long double complex p4 = csqrt(p.b*p.b/(4*p.a*p.a)-2*p.c/(3*p.a)+p3);
    long double complex p5 = p.b*p.b/(2*p.a*p.a)-4*p.c/(3*p.a)-p3;
    long double complex p6 = (-cpow(p.b/p.a,3.)+4*p.b*p.c/(p.a*p.a)-8*p.d/p.a)/(4*p4);
    QuarticSolutions result;
    result.x1 = -p.b/(4*p.a)-p4/2-csqrt(p5-p6)/2;
    result.x2 = -p.b/(4*p.a)-p4/2+csqrt(p5-p6)/2;
    result.x3 = -p.b/(4*p.a)+p4/2-csqrt(p5+p6)/2;
    result.x4 = -p.b/(4*p.a)+p4/2+csqrt(p5+p6)/2;
    return result;
}

/*
 * FIRST BLOCK (f1, f2)
 */
// used in f1.c and f2.c
long double suma(long double complex (*f)(long double complex, Parameters, void*), int64_t t, int64_t A, int64_t n, Parameters prm, void* prm2) {
    long double sum_a_k = 0;
    int64_t sign;
    for (int_fast64_t k=1;k<=n;k++) {
        sign = k%2 ? -1 : 1;
        sum_a_k += sign*creal(f((A+2.*k*M_PI*I)/(2.*t),prm,prm2));
    }
    long double s_n = exp(A/2.)/(2.*t)*creal(f(A/(2.*t),prm,prm2))+exp(A/2.)/t*sum_a_k;
    return s_n;
}

// used in f1.c and f2.c
long double euler(long double complex (*f)(long double complex, Parameters, void*), int64_t t, int64_t A, int64_t n, int64_t B, Parameters prm, void* prm2) {
    long double E = pow(2.,(long double)(-n))*suma(f,t,A,B,prm,prm2);
    for (int_fast64_t k=1;k<=n;k++) {
        E += pow(2.,(long double)(-n))*binomial(n,k)*suma(f,t,A,B+k,prm,prm2);
    }
    return E;
}

// used in f1.c and f3.c
long double complex hat_f1(long double complex alpha, Parameters prm, void *prm2) {
    Parametersf1 *prm2_ = (Parametersf1*)(prm2);

    Quartic p = polynomial(alpha, prm);
    QuarticSolutions beta = quartic_solve(p);

    long double complex max1,min1,max2,min2;
    
    if (creal(beta.x1)>creal(beta.x2)) {
        max1 = beta.x1;
        min1 = beta.x2;
    } else {
        max1 = beta.x2;
        min1 = beta.x1;
    }
    if (creal(beta.x3)>creal(beta.x4)) {
        max2 = beta.x3;
        min2 = beta.x4;
    } else {
        max2 = beta.x4;
        min2 = beta.x3;
    }

    long double complex beta1,beta2;
    if (creal(max1)<creal(min2)) {
        beta1 = min2;
        beta2 = max2;
    } else if (creal(max2)<creal(min1)) {
        beta1 = min1;
        beta2 = max1;
    } else if (creal(max1)<creal(max2)) {
        beta1 = max1;
        beta2 = max2;
    } else {
        beta1 = max2;
        beta2 = max1;
    }

    long double complex denom = alpha*prm.eta1*(beta2-beta1);
    long double complex num = beta2*(prm.eta1-beta1)*cexp(-beta1*prm2_->b)+beta1*(beta2-prm.eta1)*cexp(-beta2*prm2_->b);
    return num/denom;
}

// used in f2.c
long double complex aux(long double complex beta1, long double complex beta2, long double eta1, long double b) {
    return (eta1-beta1)*cexp(-beta1*b)/(beta2-beta1);
}

// used in f2.c
long double complex dG(long double complex z, Parameters prm) {
    return prm.mu+prm.sigma*prm.sigma*z+prm.lambda*(prm.p*prm.eta1/((prm.eta1-z)*(prm.eta1-z))-(1-prm.p)*prm.eta2/((prm.eta2+z)*(prm.eta2+z)));
}

// used in f2.c
long double complex hat_f2(long double complex alpha, Parameters prm, void *prm2) {
    Parametersf2 *prm2_ = (Parametersf2*)(prm2);

    Quartic p = polynomial(alpha, prm);
    QuarticSolutions beta = quartic_solve(p);

    long double complex max1,min1,max2,min2;
    
    if (creal(beta.x1)>creal(beta.x2)) {
        max1 = beta.x1;
        min1 = beta.x2;
    } else {
        max1 = beta.x2;
        min1 = beta.x1;
    }
    if (creal(beta.x3)>creal(beta.x4)) {
        max2 = beta.x3;
        min2 = beta.x4;
    } else {
        max2 = beta.x4;
        min2 = beta.x3;
    }

    long double complex beta1,beta2,beta3,beta4;
    if (creal(max1)<creal(min2)) {
        beta1 = min2;
        beta2 = max2;
        beta3 = -max1;
        beta4 = -min1;
    } else if (creal(max2)<creal(min1)) {
        beta1 = min1;
        beta2 = max1;
        beta3 = -max2;
        beta4 = -min2;
    } else if (creal(max1)<creal(max2)) {
        beta1 = max1;
        beta2 = max2;
        if (creal(-min1)<creal(-min2)) {
            beta3 = -min2;
            beta4 = -min1;    
        } else {
            beta3 = -min1;
            beta4 = -min2;
        }
    } else {
        beta1 = max2;
        beta2 = max1;
        if (creal(-min1)<creal(-min2)) {
            beta3 = -min2;
            beta4 = -min1;    
        } else {
            beta3 = -min1;
            beta4 = -min2;
        }
    }

    complex long double Aalpha, Balpha;
    Aalpha = aux(beta1,beta2,prm.eta1,prm2_->b)+aux(beta2,beta1,prm.eta1,prm2_->b);
    Balpha = aux(beta1,beta2,prm.eta1,prm2_->b)*(beta2-prm.eta1)/prm.eta1
        -aux(beta2,beta1,prm.eta1,prm2_->b)*(prm.eta1-beta1)/prm.eta1;

    long double complex C3,C4,D3,D4;
    C3 = 1./(beta3*dG(-beta3,prm));
    C4 = 1./(beta4*dG(-beta4,prm));
    D3 = prm.eta1/((prm.eta1+beta3)*beta3*dG(-beta3,prm));
    D4 = prm.eta2/((prm.eta2+beta4)*beta4*dG(-beta4,prm));

    long double c = prm2_->a-prm2_->b;

    return (Aalpha+Balpha)/alpha + (C3*Aalpha+D3*Balpha)*cexp(c*beta3)+(C4*Aalpha+D4*Balpha)*cexp(c*beta4);
}






/*
 * SECOND BLOCK (MULTIPRECISION)
 */

// factorial, uses GNU GMP to avoid integer overflow, returns result in GNU MPC format
void factorial(mpc_t result, int64_t n) {
    if (n==0) {
        mpc_set_si(result,1,MPCRND);
        return;
    }
    mpz_t multi_precision;
    mpz_init(multi_precision);
    mpz_set_ui(multi_precision,1);
    for (int_fast64_t k=n;k>0;k--) {
        mpz_mul_si(multi_precision,multi_precision,k);
    }
    mpc_set_z(result,multi_precision,MPCRND);
    mpz_clear(multi_precision);
}

// coefs. of the polynomial of degree 4 (multiprecision version)
Quartic_mp polynomial_mp(mpc_t alpha, Parameters_mp prm) {
    // aux vars
    mpfr_t aux, sigma2;
    mpfr_inits2(prm.bits,aux,sigma2,NULL);
    mpfr_sqr(sigma2,prm.sigma,MPFRND); // sigma2 = sigma^2

    Quartic_mp p;
    // p.a
    mpc_init2(p.a,prm.bits);
    mpc_set_si(p.a,1,MPCRND);                   // p.a = 1
    // p.b
    mpc_init2(p.b,prm.bits);
    mpc_set_fr(p.b,prm.mu,MPCRND);              // p.b = mu
    mpc_mul_ui(p.b,p.b,2,MPCRND);               //          * 2
    mpc_div_fr(p.b,p.b,sigma2,MPCRND);          //          / sigma^2
    mpc_add_fr(p.b,p.b,prm.eta2,MPCRND);        //          + eta2
    mpc_sub_fr(p.b,p.b,prm.eta1,MPCRND);        //          - eta1
    // p.c
    mpc_init2(p.c,prm.bits);
    mpc_set_fr(p.c,prm.eta1,MPCRND);            // p.c = eta1
    mpc_mul_fr(p.c,p.c,prm.mu,MPCRND);          //          * mu
    mpc_neg(p.c,p.c,MPCRND);                    //          * -1
    mpc_sub_fr(p.c,p.c,prm.lambda,MPCRND);      //          - lambda
    mpc_sub(p.c,p.c,alpha,MPCRND);              //          - alpha
    mpfr_set(aux,prm.mu,MPFRND);                // aux = mu
    mpfr_mul(aux,aux,prm.eta2,MPFRND);          //          * eta2
    mpc_add_fr(p.c,p.c,aux,MPCRND);             // p.c += aux
    mpc_div_fr(p.c,p.c,sigma2,MPCRND);          //          / sigma^2
    mpc_mul_ui(p.c,p.c,2,MPCRND);               //          * 2
    mpfr_set(aux,prm.eta2,MPFRND);              // aux = eta2
    mpfr_mul(aux,aux,prm.eta1,MPFRND);          //          * eta1
    mpfr_neg(aux,aux,MPFRND);                   //          * -1
    mpc_add_fr(p.c,p.c,aux,MPCRND);             // p.c += aux
    // p.d
    mpc_init2(p.d,prm.bits);
    mpc_set_fr(p.d,prm.eta1,MPCRND);            // p.d = eta1
    mpc_sub_fr(p.d,p.d,prm.eta2,MPCRND);            //          - eta2
    mpc_mul(p.d,p.d,alpha,MPCRND);              //          * alpha
    mpfr_set(aux,prm.p,MPFRND);                 // aux = p
    mpfr_mul(aux,aux,prm.lambda,MPFRND);        //          * lambda
    mpfr_mul(aux,aux,prm.eta2,MPFRND);          //          * eta2
    mpc_sub_fr(p.d,p.d,aux,MPCRND);             // p.d -= aux
    mpfr_set(aux,prm.p,MPFRND);                 // aux = p
    mpfr_sub_ui(aux,aux,1,MPFRND);              //          - 1
    mpfr_mul(aux,aux,prm.lambda,MPFRND);        //          * lambda
    mpfr_mul(aux,aux,prm.eta1,MPFRND);          //          * eta1
    mpc_sub_fr(p.d,p.d,aux,MPCRND);             // p.d -= aux
    mpfr_set(aux,prm.mu,MPFRND);                // aux = mu
    mpfr_mul(aux,aux,prm.eta2,MPFRND);          //          * eta2
    mpfr_mul(aux,aux,prm.eta1,MPFRND);          //          * eta1
    mpc_sub_fr(p.d,p.d,aux,MPCRND);             // p.d -= aux
    mpc_div_fr(p.d,p.d,sigma2,MPCRND);          //          / sigma^2
    mpc_mul_ui(p.d,p.d,2,MPCRND);               //          * 2
    // p.e
    mpc_init2(p.e,prm.bits);
    mpc_set(p.e,alpha,MPCRND);                  // p.e = alpha
    mpc_div_fr(p.e,p.e,sigma2,MPCRND);          //          / sigma^2
    mpc_mul_fr(p.e,p.e,prm.eta2,MPCRND);        //          * eta2
    mpc_mul_fr(p.e,p.e,prm.eta1,MPCRND);        //          * eta1
    mpc_mul_ui(p.e,p.e,2,MPCRND);               //          * 2

    mpfr_clears(aux,sigma2,NULL);
    return p;
}

// direct computation of the roots of the polynomial of degree 4 (multiprecision version)
//http://math.stackexchange.com/a/786
QuarticSolutions_mp quartic_solve_mp(Quartic_mp p, mpfr_prec_t bits) {
    // aux
    mpc_t aux, aux2;
    mpc_init2(aux,bits);
    mpc_init2(aux2,bits);
    // p1
    mpc_t p1;
    mpc_init2(p1,bits);
    mpc_pow_ui(p1,p.c,3,MPCRND);                // p1 = (p.c)^3
    mpc_mul_ui(p1,p1,2,MPCRND);                 //          * 2
    mpc_set(aux,p.b,MPCRND);                    // aux = p.b
    mpc_mul(aux,aux,p.c,MPCRND);                //          * p.c
    mpc_mul(aux,aux,p.d,MPCRND);                //          * p.d
    mpc_mul_ui(aux,aux,9,MPCRND);               //          * 9
    mpc_sub(p1,p1,aux,MPCRND);                  // p1 -= aux
    mpc_sqr(aux,p.d,MPCRND);                    // aux = (p.d)^2
    mpc_mul(aux,aux,p.a,MPCRND);                //          * p.a
    mpc_mul_ui(aux,aux,27,MPCRND);              //          * 27
    mpc_add(p1,p1,aux,MPCRND);                  // p1 += aux
    mpc_sqr(aux,p.b,MPCRND);                    // aux = (p.b)^2
    mpc_mul(aux,aux,p.e,MPCRND);                //          * p.e
    mpc_mul_ui(aux,aux,27,MPCRND);              //          * 27
    mpc_add(p1,p1,aux,MPCRND);                  // p1 += aux
    mpc_set(aux,p.a,MPCRND);                    // aux = p.a
    mpc_mul(aux,aux,p.c,MPCRND);                //          * p.c
    mpc_mul(aux,aux,p.e,MPCRND);                //          * p.e
    mpc_mul_ui(aux,aux,72,MPCRND);              //          * 72
    mpc_sub(p1,p1,aux,MPCRND);                  // p1 -= aux
    // p2
    mpc_t p2;
    mpc_init2(p2,bits);
    mpc_sqr(p2,p.c,MPCRND);                     // p2 = (p.c)^2
    mpc_set(aux,p.b,MPCRND);                    // aux = p.b
    mpc_mul(aux,aux,p.d,MPCRND);                //          * p.d
    mpc_mul_ui(aux,aux,3,MPCRND);               //          * 3
    mpc_sub(p2,p2,aux,MPCRND);                  // p2 -= aux
    mpc_set(aux,p.a,MPCRND);                    // aux = p.a
    mpc_mul(aux,aux,p.e,MPCRND);                //          * p.e
    mpc_mul_ui(aux,aux,12,MPCRND);              //          * 12
    mpc_add(p2,p2,aux,MPCRND);                  // p2 += aux
    mpc_pow_ui(p2,p2,3,MPCRND);                 //          ^ 3
    mpc_mul_si(p2,p2,-4,MPCRND);                //          * -4
    mpc_sqr(aux,p1,MPCRND);                     // aux = p1^2
    mpc_add(p2,p2,aux,MPCRND);                  // p2 += aux
    mpc_sqrt(p2,p2,MPCRND);                     //          ^ 1/2
    mpc_add(p2,p2,p1,MPCRND);
    // p3
    mpc_t p3;
    mpc_init2(p3,bits);
    mpc_sqr(p3,p.c,MPCRND);                     // p3 = (p.c)^2
    mpc_set(aux,p.b,MPCRND);                    // aux = p.b
    mpc_mul(aux,aux,p.d,MPCRND);                //          * p.b
    mpc_mul_ui(aux,aux,3,MPCRND);               //          * 3
    mpc_sub(p3,p3,aux,MPCRND);                  // p3 -= aux
    mpc_set(aux,p.a,MPCRND);                    // aux = p.a
    mpc_mul(aux,aux,p.e,MPCRND);                //          * p.e
    mpc_mul_ui(aux,aux,12,MPCRND);              //          * 12
    mpc_add(p3,p3,aux,MPCRND);                  // p3 += aux
    mpc_set(aux,p2,MPCRND);                     // aux = p2
    mpc_div_ui(aux,aux,2,MPCRND);               //          / 2
    mpc_set_ui(aux2,1,MPCRND);                  // aux2 = 1
    mpc_div_ui(aux2,aux2,3,MPCRND);             //          / 3
    mpc_pow(aux,aux,aux2,MPCRND);               // aux ^= aux2
    mpc_mul(aux,aux,p.a,MPCRND);                //          * p.a
    mpc_mul_ui(aux,aux,3,MPCRND);               //          * 3
    mpc_div(p3,p3,aux,MPCRND);                  // p3 /= aux
    mpc_set(aux,p2,MPCRND);                     // aux = p2
    mpc_div_ui(aux,aux,2,MPCRND);               //          / 2
    mpc_pow(aux,aux,aux2,MPCRND);               // aux ^= aux2
    mpc_div(aux,aux,p.a,MPCRND);                //          / p.a
    mpc_div_ui(aux,aux,3,MPCRND);               //          / 3
    mpc_add(p3,p3,aux,MPCRND);                  // p3 += aux
    // p4
    mpc_t p4;
    mpc_init2(p4,bits);
    mpc_sqr(p4,p.b,MPCRND);                     // p4 = (p.b)^2
    mpc_sqr(aux,p.a,MPCRND);                    // aux = (p.a)^2
    mpc_mul_ui(aux,aux,4,MPCRND);               //          * 4
    mpc_div(p4,p4,aux,MPCRND);                  // p4 /= aux
    mpc_set(aux,p.c,MPCRND);                    // aux = p.c
    mpc_mul_ui(aux,aux,2,MPCRND);               //          * 2
    mpc_div_ui(aux,aux,3,MPCRND);               //          / 3
    mpc_div(aux,aux,p.a,MPCRND);                //          / p.a
    mpc_sub(p4,p4,aux,MPCRND);                  // p4 -= aux
    mpc_add(p4,p4,p3,MPCRND);                   //          + p3
    mpc_sqrt(p4,p4,MPCRND);                     //          sqrt
    // p5
    mpc_t p5;
    mpc_init2(p5,bits);
    mpc_sqr(p5,p.b,MPCRND);                     // p5 = (p.b)^2
    mpc_sqr(aux,p.a,MPCRND);                    // aux = (p.a)^2
    mpc_mul_ui(aux,aux,2,MPCRND);               //          * 2
    mpc_div(p5,p5,aux,MPCRND);                  // p5 /= aux
    mpc_set(aux,p.c,MPCRND);                    // aux = p.c
    mpc_mul_ui(aux,aux,4,MPCRND);               //          * 4
    mpc_div_ui(aux,aux,3,MPCRND);               //          / 3
    mpc_div(aux,aux,p.a,MPCRND);                //          / p.a
    mpc_sub(p5,p5,aux,MPCRND);                  // p5 -= aux
    mpc_sub(p5,p5,p3,MPCRND);                   //          - p3
    // p6
    mpc_t p6;
    mpc_init2(p6,bits);
    mpc_div(p6,p.b,p.a,MPCRND);                 // p6 = p.b / p.a
    mpc_pow_ui(p6,p6,3,MPCRND);                 //          ^ 3
    mpc_neg(p6,p6,MPCRND);                      //          * -1
    mpc_mul(aux,p.b,p.c,MPCRND);                // aux = p.b * p.c
    mpc_mul_ui(aux,aux,4,MPCRND);               //          * 4
    mpc_sqr(aux2,p.a,MPCRND);                   // aux2 = (p.a)^2
    mpc_div(aux,aux,aux2,MPCRND);               // aux /= aux2
    mpc_add(p6,p6,aux,MPCRND);                  // p6 += aux
    mpc_set(aux,p.d,MPCRND);                    // aux = p.d
    mpc_mul_ui(aux,aux,8,MPCRND);               //          * 8
    mpc_div(aux,aux,p.a,MPCRND);                //          / p.a
    mpc_sub(p6,p6,aux,MPCRND);                  // p6 -= aux
    mpc_div_ui(p6,p6,4,MPCRND);                 //          / 4
    mpc_div(p6,p6,p4,MPCRND);                   //          / p4

    QuarticSolutions_mp rs;
    // rs.x1
    mpc_init2(rs.x1,bits);
    mpc_set(rs.x1,p.b,MPCRND);                  // rs.x1 = p.b
    mpc_div_ui(rs.x1,rs.x1,4,MPCRND);           //          / 4
    mpc_div(rs.x1,rs.x1,p.a,MPCRND);            //          / p.a
    mpc_set(aux,p4,MPCRND);                     // aux = p4
    mpc_div_ui(aux,aux,2,MPCRND);               //          / 2
    mpc_add(rs.x1,rs.x1,aux,MPCRND);            // rs.x1 += aux
    mpc_set(aux,p5,MPCRND);                     // aux = p5
    mpc_sub(aux,aux,p6,MPCRND);                 //          - p6
    mpc_sqrt(aux,aux,MPCRND);                   //          ^ 1/2
    mpc_div_ui(aux,aux,2,MPCRND);               //          / 2
    mpc_add(rs.x1,rs.x1,aux,MPCRND);            // rs.x1 += aux
    mpc_neg(rs.x1,rs.x1,MPCRND);                //          * -1
    // rs.x2
    mpc_init2(rs.x2,bits);
    mpc_set(rs.x2,p.b,MPCRND);                  // rs.x2 = p.b
    mpc_div_ui(rs.x2,rs.x2,4,MPCRND);           //          / 4
    mpc_div(rs.x2,rs.x2,p.a,MPCRND);            //          / p.a
    mpc_set(aux,p4,MPCRND);                     // aux = p4
    mpc_div_ui(aux,aux,2,MPCRND);               //          / 2
    mpc_add(rs.x2,rs.x2,aux,MPCRND);            // rs.x2 += aux
    mpc_set(aux,p5,MPCRND);                     // aux = p5
    mpc_sub(aux,aux,p6,MPCRND);                 //          - p6
    mpc_sqrt(aux,aux,MPCRND);                   //          ^ 1/2
    mpc_div_ui(aux,aux,2,MPCRND);               //          / 2
    mpc_sub(rs.x2,rs.x2,aux,MPCRND);            // rs.x2 -= aux
    mpc_neg(rs.x2,rs.x2,MPCRND);                //          * -1
    // rs.x3
    mpc_init2(rs.x3,bits);
    mpc_set(rs.x3,p.b,MPCRND);                  // rs.x3 = p.b
    mpc_div_ui(rs.x3,rs.x3,4,MPCRND);           //          / 4
    mpc_div(rs.x3,rs.x3,p.a,MPCRND);            //          / p.a
    mpc_set(aux,p4,MPCRND);                     // aux = p4
    mpc_div_ui(aux,aux,2,MPCRND);               //          / 2
    mpc_sub(rs.x3,rs.x3,aux,MPCRND);            // rs.x3 -= aux
    mpc_set(aux,p5,MPCRND);                     // aux = p5
    mpc_add(aux,aux,p6,MPCRND);                 //          + p6
    mpc_sqrt(aux,aux,MPCRND);                   //          ^ 1/2
    mpc_div_ui(aux,aux,2,MPCRND);               //          / 2
    mpc_add(rs.x3,rs.x3,aux,MPCRND);            // rs.x3 += aux
    mpc_neg(rs.x3,rs.x3,MPCRND);                //          * -1
    // rs.x4
    mpc_init2(rs.x4,bits);
    mpc_set(rs.x4,p.b,MPCRND);                  // rs.x4 = p.b
    mpc_div_ui(rs.x4,rs.x4,4,MPCRND);           //          / 4
    mpc_div(rs.x4,rs.x4,p.a,MPCRND);            //          / p.a
    mpc_set(aux,p4,MPCRND);                     // aux = p4
    mpc_div_ui(aux,aux,2,MPCRND);               //          / 2
    mpc_sub(rs.x4,rs.x4,aux,MPCRND);            // rs.x4 -= aux
    mpc_set(aux,p5,MPCRND);                     // aux = p5
    mpc_add(aux,aux,p6,MPCRND);                 //          + p6
    mpc_sqrt(aux,aux,MPCRND);                   //          ^ 1/2
    mpc_div_ui(aux,aux,2,MPCRND);               //          / 2
    mpc_sub(rs.x4,rs.x4,aux,MPCRND);            // rs.x4 -= aux
    mpc_neg(rs.x4,rs.x4,MPCRND);                //          * -1

    mpc_clear(aux); mpc_clear(aux2);

    return rs;
}

// used in f3.c
void hat_f1_mp(mpc_t result, mpc_t alpha, Parameters_mp prm, void *prm2) {
    Parametersf3_mp *prm2_ = (Parametersf3_mp*)(prm2);

    Quartic_mp p = polynomial_mp(alpha, prm);
    QuarticSolutions_mp beta = quartic_solve_mp(p, prm.bits);

    mpc_t max1,min1,max2,min2;
    mpc_init2(max1,prm.bits);
    mpc_init2(min1,prm.bits);
    mpc_init2(max2,prm.bits);
    mpc_init2(min2,prm.bits);
    
    if (mpfr_greater_p(mpc_realref(beta.x1),mpc_realref(beta.x2))) { // Re(beta.x1)>Re(beta.x2) ?
        mpc_set(max1,beta.x1,MPCRND);           // max1 = beta.x1
        mpc_set(min1,beta.x2,MPCRND);           // min1 = beta.x2
    } else {
        mpc_set(max1,beta.x2,MPCRND);           // max1 = beta.x2
        mpc_set(min1,beta.x1,MPCRND);           // min1 = beta.x1
    }
    if (mpfr_greater_p(mpc_realref(beta.x3),mpc_realref(beta.x4))) { // Re(beta.x3)>Re(beta.x4) ?
        mpc_set(max2,beta.x3,MPCRND);           // max2 = beta.x3
        mpc_set(min2,beta.x4,MPCRND);           // min1 = beta.x4
    } else {
        mpc_set(max2,beta.x4,MPCRND);           // max2 = beta.x4
        mpc_set(min2,beta.x3,MPCRND);           // min2 = beta.x3
    }

    mpc_t beta1,beta2;
    mpc_init2(beta1,prm.bits);
    mpc_init2(beta2,prm.bits);

    if (mpfr_less_p(mpc_realref(max1),mpc_realref(min2))) { // Re(max1)<Re(min2) ?
        mpc_set(beta1,min2,MPCRND);             // beta1 = min2
        mpc_set(beta2,max2,MPCRND);             // beta2 = max2
    } else if (mpfr_less_p(mpc_realref(max2),mpc_realref(min1))) { // Re(max2)<Re(min1) ?
        mpc_set(beta1,min1,MPCRND);             // beta1 = min1
        mpc_set(beta2,max1,MPCRND);             // beta2 = max1
    } else if (mpfr_less_p(mpc_realref(max1),mpc_realref(max2))) { // Re(max1)<Re(max2) ?
        mpc_set(beta1,max1,MPCRND);             // beta1 = max1
        mpc_set(beta2,max2,MPCRND);             // beta2 = max2
    } else {
        mpc_set(beta1,max2,MPCRND);             // beta1 = max2
        mpc_set(beta2,max1,MPCRND);             // beta2 = max1
    }

    mpc_clear(max1); mpc_clear(min1);
    mpc_clear(max2); mpc_clear(min2);

    mpc_t num, denom, aux, aux2;
    mpc_init2(num,prm.bits);
    mpc_init2(denom,prm.bits);
    mpc_init2(aux,prm.bits);
    mpc_init2(aux2,prm.bits);

    mpc_set(denom,beta2,MPCRND);                // denom = beta2
    mpc_sub(denom,denom,beta1,MPCRND);          //          - beta1
    mpc_mul_fr(denom,denom,prm.eta1,MPCRND);    //          * eta1
    mpc_mul(denom,denom,alpha,MPCRND);          //          * alpha
    mpc_set(num,beta1,MPCRND);                  // num = beta1
    mpc_neg(num,num,MPCRND);                    //          * -1
    mpc_add_fr(num,num,prm.eta1,MPCRND);        //          + eta1
    mpc_mul(num,num,beta2,MPCRND);              //          * beta2
    mpc_set(aux,beta1,MPCRND);                  // aux = beta1
    mpc_neg(aux,aux,MPCRND);                    //          * -1
    mpc_mul_fr(aux,aux,prm2_->b,MPCRND);        //          * b
    mpc_exp(aux,aux,MPCRND);                    //          exp
    mpc_mul(num,num,aux,MPCRND);                // num *= aux
    mpc_set(aux,beta2,MPCRND);                  // aux = beta2
    mpc_sub_fr(aux,aux,prm.eta1,MPCRND);        //          - eta1
    mpc_mul(aux,aux,beta1,MPCRND);              //          * beta1
    mpc_set(aux2,beta2,MPCRND);                 // aux2 = beta2
    mpc_neg(aux2,aux2,MPCRND);                  //          * -1
    mpc_mul_fr(aux2,aux2,prm2_->b,MPCRND);      //          * b
    mpc_exp(aux2,aux2,MPCRND);                  //          exp
    mpc_mul(aux,aux,aux2,MPCRND);               // aux *= aux2
    mpc_add(num,num,aux,MPCRND);                // num += aux
    mpc_div(result,num,denom,MPCRND);           // result = num/denom

    mpc_clear(beta1); mpc_clear(beta2);
    mpc_clear(aux); mpc_clear(aux2);
    mpc_clear(num); mpc_clear(denom);
}



// used in f3.c
void fn(mpc_t suma, void (*f)(mpc_t, mpc_t, Parameters_mp, void*), int64_t t, int64_t n, Parameters_mp prm, void *prm2) {
    int64_t sign;
    
    // variables for calling f
    mpc_t resultf, alpha;
    mpc_init2(resultf,prm.bits);
    mpc_init2(alpha,prm.bits);

    // log(2)
    mpc_t log2;
    mpc_init2(log2,prm.bits);
    mpc_set_ui(log2,2,MPCRND);                  // log2 = 2
    mpc_log(log2,log2,MPCRND);                  //          log
    
    // factorials
    mpc_t f2n;
    mpc_init2(f2n,prm.bits);
    factorial(f2n,2*n);                         // f2n = (2n)!

    mpc_t fnm1;
    mpc_init2(fnm1,prm.bits);
    factorial(fnm1,n-1);                        // fnm1 = (n-1)!
    mpc_t fk,fnmk;
    mpc_init2(fk,prm.bits);
    mpc_init2(fnmk,prm.bits);

    // loop
    mpc_set_si(suma,0,MPCRND);                  // suma = 0
    for (int64_t k=0; k<=n; k++) {
        sign = k%2 ? -1 : 1;
        factorial(fk,k);                        // fk = k!
        factorial(fnmk,n-k);                    // fnmk = (n-k)!
        mpc_mul(fk,fk,fnmk,MPCRND);             // fk *= fnmk
        mpc_mul(fk,fk,fnm1,MPCRND);             //          * fnm1
        mpc_mul_si(fk,fk,t,MPCRND);             //          * t
        mpc_div(fk,f2n,fk,MPCRND);              // fk = f2n / fk
        mpc_mul(fk,log2,fk,MPCRND);             //          * log2
        mpc_mul_si(fk,fk,sign,MPCRND);          //          * sign
        mpc_set_ui(alpha,n+k,MPCRND);           // alpha = n
        mpc_mul(alpha,alpha,log2,MPCRND);       //          * log2
        mpc_div_ui(alpha,alpha,t,MPCRND);       //          / t
        f(resultf,alpha,prm,prm2);              // resultf = f(alpha)
        mpc_mul(fk,fk,resultf,MPCRND);          // fk *= resultf
        mpc_add(suma,suma,fk,MPCRND);           // suma += fk        
    }

    mpc_clear(log2);
    mpc_clear(f2n); mpc_clear(fnm1); mpc_clear(fk); mpc_clear(fnmk);
    mpc_clear(resultf); mpc_clear(alpha);
}

// used in f3.c
void gaver(mpc_t result, void (*f)(mpc_t, mpc_t, Parameters_mp, void*), int64_t t, int64_t n, int64_t B, Parameters_mp prm, void *prm2) {
    int64_t sign;

    // aux variables
    mpc_t aux;
    mpc_init2(aux,prm.bits);

    // variables for calling fn
    mpc_t suma_fn;
    mpc_init2(suma_fn,prm.bits);

    // factorials
    mpc_t fk, fnmk;
    mpc_init2(fk,prm.bits);
    mpc_init2(fnmk,prm.bits);

    mpc_set_ui(result,0,MPCRND);
    for (int_fast64_t k=1; k<=n; k++) {
        sign = (n-k)%2 ? -1 : 1;
        factorial(fk,k);                        // fk = k!
        factorial(fnmk,n-k);                    // fnmk = (n-k)!
        mpc_set_ui(aux,k,MPCRND);               // aux = k
        mpc_pow_ui(aux,aux,n,MPCRND);           //          ^ n
        mpc_mul_si(aux,aux,sign,MPCRND);        //          * sign
        mpc_div(aux,aux,fk,MPCRND);             //          / fk
        mpc_div(aux,aux,fnmk,MPCRND);           //          / fnmk
        fn(suma_fn,f,t,k+B,prm,prm2);           // suma_fn = fn(...)
        mpc_mul(aux,aux,suma_fn,MPCRND);        // aux *= suma_fn
        mpc_add(result,result,aux,MPCRND);      // result += aux

    }

    mpc_clear(aux);
    mpc_clear(suma_fn);
    mpc_clear(fk); mpc_clear(fnmk);
}





// resolvent (resultant of the polynomial P_alpha, numerator
// of G(z) in [5], and its derivative)
Resolvent resolvent(Parameters prm) {
    long double a,b,c,d,e,f;
    
    long double mu,sigma,eta1,eta2,lambda,p;
    mu = prm.mu; sigma = prm.sigma; eta1 = prm.eta1; eta2 = prm.eta2; lambda = prm.lambda; p = prm.p;
    
    long double aux1,aux2,aux3,aux4,aux5,aux6;

    long double sigma2 = sigma*sigma;
    long double sigma3 = sigma*sigma2;
    long double sigma4 = sigma*sigma3;
    long double sigma5 = sigma*sigma4;
    long double sigma6 = sigma*sigma5;
    long double sigma8 = sigma2*sigma6;
    long double sigma10 = sigma2*sigma8;
    long double sigma12 = sigma2*sigma10;

    long double eta12 = eta1*eta1;
    long double eta13 = eta1*eta12;
    long double eta14 = eta1*eta13;
    long double eta15 = eta1*eta14;
    long double eta16 = eta1*eta15;

    long double eta22 = eta2*eta2;
    long double eta23 = eta2*eta22;
    long double eta24 = eta2*eta23;
    long double eta25 = eta2*eta24;
    long double eta26 = eta2*eta25;

    long double mu2 = mu*mu;
    long double mu3 = mu*mu2;
    long double mu4 = mu*mu3;
    long double mu5 = mu*mu4;
    long double mu6 = mu*mu5;

    long double lambda2 = lambda*lambda;
    long double lambda3 = lambda*lambda2;
    long double lambda4 = lambda*lambda3;
    long double lambda5 = lambda*lambda4;

    long double p2 = p*p;
    long double p3 = p*p2;
    long double p4 = p*p3;

    long double eta1p2  = eta1+eta2;

    // resolvent is a*alpha^5+...+e*alpha+f
    // a
    a = -128.*sigma4*(eta1p2*eta1p2);
    // b
    aux1 = 128*(eta12+eta22)*(eta1p2*eta1p2)*sigma6;
    aux2 = 128*eta1p2*(2*eta12*mu+2*eta1*lambda*p-2*eta22*mu-2*eta2*lambda*p-5*eta1*lambda
        -3*eta2*lambda)*sigma4;
    aux3 = -64*mu2*(eta1p2*eta1p2)*sigma2;
    b = aux1+aux2+aux3;
    // c
    aux1 = -32*(eta14+4*eta12*eta22+eta24)*(eta1p2*eta1p2)*sigma8;
    aux2 = 128*mu2*eta1p2*(eta12*mu+eta1*lambda*p-eta22*mu-eta2*lambda*p-2*eta1*lambda-eta2*lambda)*sigma2;
    aux3 = (-64*eta14*mu2+384*eta13*eta2*mu2-1088*eta13*lambda*mu*p+896*eta12*eta22*mu2
        -1728*eta12*eta2*lambda*mu*p-128*eta12*lambda2*p2+384*eta1*eta23*mu2-1728*eta1*eta22*lambda*mu*p
        -256*eta1*eta2*lambda2*p2-64*eta24*mu2-1088*eta23*lambda*mu*p-128*eta22*lambda2*p2
        +1024*eta13*lambda*mu+1344*eta12*eta2*lambda*mu+1024*eta12*lambda2*p+384*eta1*eta22*lambda*mu
        +256*eta1*eta2*lambda2*p+64*eta23*lambda*mu-768*eta22*lambda2*p-1280*eta12*lambda2
        -1536*eta1*eta2*lambda2-384*eta22*lambda2)*sigma4;
    aux4 = -(64*eta1p2)*(2*eta14*mu-2*eta13*eta2*mu+13*eta13*lambda*p+7*eta12*eta2*lambda*p
        +2*eta1*eta23*mu-7*eta1*eta22*lambda*p-2*eta24*mu-13*eta23*lambda*p-8*eta13*lambda
        -6*eta12*eta2*lambda+eta1*eta22*lambda+5*eta23*lambda)*sigma6;
    c = aux1+aux2+aux3+aux4;
    // d
    aux1 = 32*eta22*eta12*(eta1p2*eta1p2)*(eta12+eta22)*sigma10;
    aux2 = -(32*eta1p2)*(2*eta15*eta2*mu-3*eta15*lambda*p-2*eta14*eta22*mu-eta14*eta2*lambda*p
        -4*eta13*eta22*lambda*p+2*eta12*eta24*mu+4*eta12*eta23*lambda*p-2*eta1*eta25*mu
        +eta1*eta24*lambda*p+3*eta25*lambda*p+3*eta15*lambda+3*eta14*eta2*lambda-6*eta13*eta22*lambda
        -10*eta12*eta23*lambda+2*eta1*eta24*lambda)*sigma8;
    aux3 = -64*mu2*(eta14*mu2-2*eta13*eta2*mu2+8*eta13*lambda*mu*p-6*eta12*eta22*mu2
        +12*eta12*eta2*lambda*mu*p+eta12*lambda2*p2-2*eta1*eta23*mu2+12*eta1*eta22*lambda*mu*p
        +2*eta1*eta2*lambda2*p2+eta24*mu2+8*eta23*lambda*mu*p+eta22*lambda2*p2-6*eta13*lambda*mu
        -8*eta12*eta2*lambda*mu-6*eta12*lambda2*p-4*eta1*eta22*lambda*mu-2*eta1*eta2*lambda2*p
        -2*eta23*lambda*mu+4*eta22*lambda2*p+6*eta12*lambda2+6*eta1*eta2*lambda2+eta22*lambda2)*sigma2;
    aux4 = (-64*eta15*mu3-256*eta14*eta2*mu3-64*eta14*lambda*mu2*p-192*eta13*eta22*mu3
        -704*eta13*eta2*lambda*mu2*p+1408*eta13*lambda2*mu*p2+192*eta12*eta23*mu3
        +1408*eta12*eta2*lambda2*mu*p2+256*eta1*eta24*mu3+704*eta1*eta23*lambda*mu2*p
        -1408*eta1*eta22*lambda2*mu*p2+64*eta25*mu3+64*eta24*lambda*mu2*p-1408*eta23*lambda2*mu*p2
        -192*eta14*lambda*mu2-64*eta13*eta2*lambda*mu2-3264*eta13*lambda2*mu*p
        -384*eta12*eta22*lambda*mu2-5440*eta12*eta2*lambda2*mu*p-384*eta12*lambda3*p2
        -768*eta1*eta23*lambda*mu2-2624*eta1*eta22*lambda2*mu*p-768*eta1*eta2*lambda3*p2
        -256*eta24*lambda*mu2-448*eta23*lambda2*mu*p-384*eta22*lambda3*p2+1536*eta13*lambda2*mu
        +2496*eta12*eta2*lambda2*mu+1536*eta12*lambda3*p+1536*eta1*eta22*lambda2*mu
        +768*eta1*eta2*lambda3*p+320*eta23*lambda2*mu-768*eta22*lambda3*p-1280*eta12*lambda3
        -1024*eta1*eta2*lambda3-128*eta22*lambda3)*sigma4;
    aux5 = (-16*eta16*mu2-288*eta15*eta2*mu2+352*eta15*lambda*mu*p-336*eta14*eta22*mu2
        -384*eta14*eta2*lambda*mu*p+1712*eta14*lambda2*p2-128*eta13*eta23*mu2
        -1376*eta13*eta22*lambda*mu*p+2112*eta13*eta2*lambda2*p2-336*eta12*eta24*mu2
        -1376*eta12*eta23*lambda*mu*p+800*eta12*eta22*lambda2*p2-288*eta1*eta25*mu2
        -384*eta1*eta24*lambda*mu*p+2112*eta1*eta23*lambda2*p2-16*eta26*mu2+352*eta25*lambda*mu*p
        +1712*eta24*lambda2*p2-384*eta15*lambda*mu-576*eta14*eta2*lambda*mu-2496*eta14*lambda2*p
        +128*eta13*eta22*lambda*mu-3936*eta13*eta2*lambda2*p+1248*eta12*eta23*lambda*mu
        -800*eta12*eta22*lambda2*p+960*eta1*eta24*lambda*mu-288*eta1*eta23*lambda2*p+32*eta25*lambda*mu
        -928*eta24*lambda2*p+768*eta14*lambda2+1152*eta13*eta2*lambda2-144*eta12*eta22*lambda2
        -672*eta1*eta23*lambda2-16*eta24*lambda2)*sigma6;
    d = aux1+aux2+aux3+aux4+aux5;
    // e
    aux1 = (32*eta16*eta23*mu+80*eta16*eta22*lambda*p+32*eta15*eta24*mu+192*eta15*eta23*lambda*p
        -32*eta14*eta25*mu-32*eta13*eta26*mu-192*eta13*eta25*lambda*p-80*eta12*eta26*lambda*p
        -80*eta16*eta22*lambda-160*eta15*eta23*lambda-48*eta14*eta24*lambda+32*eta13*eta25*lambda)*sigma10;
    aux2 = (-16*eta16*eta22*mu2-16*eta16*eta2*lambda*mu*p-96*eta16*lambda2*p2+96*eta15*eta23*mu2
        +256*eta15*eta22*lambda*mu*p-208*eta15*eta2*lambda2*p2+224*eta14*eta24*mu2
        +1232*eta14*eta23*lambda*mu*p+1152*eta14*eta22*lambda2*p2+96*eta13*eta25*mu2
        +1232*eta13*eta24*lambda*mu*p+2528*eta13*eta23*lambda2*p2-16*eta12*eta26*mu2
        +256*eta12*eta25*lambda*mu*p+1152*eta12*eta24*lambda2*p2-16*eta1*eta26*lambda*mu*p
        -208*eta1*eta25*lambda2*p2-96*eta26*lambda2*p2+16*eta16*eta2*lambda*mu+192*eta16*lambda2*p
        -288*eta15*eta22*lambda*mu+400*eta15*eta2*lambda2*p-784*eta14*eta23*lambda*mu
        -1408*eta14*eta22*lambda2*p-448*eta13*eta24*lambda*mu-2528*eta13*eta23*lambda2*p
        +32*eta12*eta25*lambda*mu-896*eta12*eta24*lambda2*p+16*eta1*eta25*lambda2*p-96*eta16*lambda2
        -192*eta15*eta2*lambda2+240*eta14*eta22*lambda2+352*eta13*eta23*lambda2
        -16*eta12*eta24*lambda2)*sigma8;
    aux3 = (-8*eta16*eta24-16*eta15*eta25-8*eta14*eta26)*sigma12;
    aux4 = (-128*eta14*eta2*mu5+128*eta14*lambda*mu4*p-128*eta13*eta22*mu5-128*eta13*eta2*lambda*mu4*p
        +640*eta13*lambda2*mu3*p2+128*eta12*eta23*mu5+640*eta12*eta2*lambda2*mu3*p2+128*eta1*eta24*mu5
        +128*eta1*eta23*lambda*mu4*p-640*eta1*eta22*lambda2*mu3*p2-128*eta24*lambda*mu4*p
        -640*eta23*lambda2*mu3*p2-128*eta14*lambda*mu4-256*eta13*eta2*lambda*mu4
        -1024*eta13*lambda2*mu3*p-512*eta12*eta22*lambda*mu4-1664*eta12*eta2*lambda2*mu3*p
        -128*eta12*lambda3*mu2*p2-384*eta1*eta23*lambda*mu4-384*eta1*eta22*lambda2*mu3*p
        -256*eta1*eta2*lambda3*mu2*p2+256*eta23*lambda2*mu3*p-128*eta22*lambda3*mu2*p2
        +384*eta13*lambda2*mu3+640*eta12*eta2*lambda2*mu3+384*eta12*lambda3*mu2*p
        +384*eta1*eta22*lambda2*mu3+256*eta1*eta2*lambda3*mu2*p-128*eta22*lambda3*mu2*p
        -256*eta12*lambda3*mu2-128*eta1*eta2*lambda3*mu2)*sigma2;
    aux5 = (-128*eta15*eta2*mu4+128*eta15*lambda*mu3*p-256*eta14*eta22*mu4-448*eta14*eta2*lambda*mu3*p
        +320*eta14*lambda2*mu2*p2-256*eta13*eta23*mu4-1728*eta13*eta22*lambda*mu3*p
        -1600*eta13*eta2*lambda2*mu2*p2-576*eta13*lambda3*mu*p3-256*eta12*eta24*mu4
        -1728*eta12*eta23*lambda*mu3*p-3840*eta12*eta22*lambda2*mu2*p2-1728*eta12*eta2*lambda3*mu*p3
        -128*eta1*eta25*mu4-448*eta1*eta24*lambda*mu3*p-1600*eta1*eta23*lambda2*mu2*p2
        -1728*eta1*eta22*lambda3*mu*p3+128*eta25*lambda*mu3*p+320*eta24*lambda2*mu2*p2
        -576*eta23*lambda3*mu*p3-128*eta15*lambda*mu3-192*eta14*eta2*lambda*mu3-128*eta14*lambda2*mu2*p
        +512*eta13*eta22*lambda*mu3+1728*eta13*eta2*lambda2*mu2*p+2816*eta13*lambda3*mu*p2
        +1216*eta12*eta23*lambda*mu3+3840*eta12*eta22*lambda2*mu2*p+4544*eta12*eta2*lambda3*mu*p2
        +640*eta1*eta24*lambda*mu3+1472*eta1*eta23*lambda2*mu2*p+640*eta1*eta22*lambda3*mu*p2
        -512*eta24*lambda2*mu2*p-1088*eta23*lambda3*mu*p2-192*eta14*lambda2*mu2
        -1280*eta13*eta2*lambda2*mu2-3264*eta13*lambda3*mu*p-2304*eta12*eta22*lambda2*mu2
        -5696*eta12*eta2*lambda3*mu*p-384*eta12*lambda4*p2-1152*eta1*eta23*lambda2*mu2
        -1792*eta1*eta22*lambda3*mu*p-768*eta1*eta2*lambda4*p2+640*eta23*lambda3*mu*p
        -384*eta22*lambda4*p2+1024*eta13*lambda3*mu+1984*eta12*eta2*lambda3*mu+1024*eta12*lambda4*p
        +896*eta1*eta22*lambda3*mu+768*eta1*eta2*lambda4*p-256*eta22*lambda4*p-640*eta12*lambda4
        -256*eta1*eta2*lambda4)*sigma4;
    aux6 = (-32*eta16*eta2*mu3+32*eta16*lambda*mu2*p-128*eta15*eta22*mu3-96*eta15*eta2*lambda*mu2*p
        -320*eta15*lambda2*mu*p2-96*eta14*eta23*mu3-576*eta14*eta22*lambda*mu2*p
        -1344*eta14*eta2*lambda2*mu*p2-1440*eta14*lambda3*p3+96*eta13*eta24*mu3
        -1024*eta13*eta22*lambda2*mu*p2-2880*eta13*eta2*lambda3*p3+128*eta12*eta25*mu3
        +576*eta12*eta24*lambda*mu2*p+1024*eta12*eta23*lambda2*mu*p2+32*eta1*eta26*mu3
        +96*eta1*eta25*lambda*mu2*p+1344*eta1*eta24*lambda2*mu*p2+2880*eta1*eta23*lambda3*p3
        -32*eta26*lambda*mu2*p+320*eta25*lambda2*mu*p2+1440*eta24*lambda3*p3-32*eta16*lambda*mu2
        +704*eta15*lambda2*mu*p+96*eta14*eta22*lambda*mu2+2400*eta14*eta2*lambda2*mu*p
        +3424*eta14*lambda3*p2-320*eta13*eta23*lambda*mu2+864*eta13*eta22*lambda2*mu*p
        +6240*eta13*eta2*lambda3*p2-480*eta12*eta24*lambda*mu2-1184*eta12*eta23*lambda2*mu*p
        +1312*eta12*eta22*lambda3*p2-96*eta1*eta25*lambda*mu2-288*eta1*eta24*lambda2*mu*p
        -2400*eta1*eta23*lambda3*p2+64*eta25*lambda2*mu*p-896*eta24*lambda3*p2-384*eta15*lambda2*mu
        -1152*eta14*eta2*lambda2*mu-2496*eta14*lambda3*p-416*eta13*eta22*lambda2*mu
        -4032*eta13*eta2*lambda3*p+576*eta12*eta23*lambda2*mu-1312*eta12*eta22*lambda3*p
        +96*eta1*eta24*lambda2*mu+192*eta1*eta23*lambda3*p-32*eta24*lambda3*p+512*eta14*lambda3
        +640*eta13*eta2*lambda3-224*eta12*eta22*lambda3-32*eta1*eta23*lambda3)*sigma6;
    e = aux1+aux2+aux3+aux4+aux5+aux6;
    // f
    aux1 = (-4*eta16*eta24*mu2-8*eta16*eta23*lambda*mu*p-4*eta16*eta22*lambda2*p2-8*eta15*eta25*mu2
        -24*eta15*eta24*lambda*mu*p-16*eta15*eta23*lambda2*p2-4*eta14*eta26*mu2
        -24*eta14*eta25*lambda*mu*p-24*eta14*eta24*lambda2*p2-8*eta13*eta26*lambda*mu*p
        -16*eta13*eta25*lambda2*p2-4*eta12*eta26*lambda2*p2+8*eta16*eta23*lambda*mu
        +8*eta16*eta22*lambda2*p+16*eta15*eta24*lambda*mu+24*eta15*eta23*lambda2*p
        +8*eta14*eta25*lambda*mu+24*eta14*eta24*lambda2*p+8*eta13*eta25*lambda2*p
        -4*eta16*eta22*lambda2-8*eta15*eta23*lambda2-4*eta14*eta24*lambda2)*sigma10;
    aux2 = (16*eta16*eta23*mu3+64*eta16*eta22*lambda*mu2*p+80*eta16*eta2*lambda2*mu*p2
        +32*eta16*lambda3*p3+16*eta15*eta24*mu3+144*eta15*eta23*lambda*mu2*p+272*eta15*eta22*lambda2*mu*p2+144*eta15*eta2*lambda3*p3-16*eta14*eta25*mu3+192*eta14*eta23*lambda2*mu*p2
        +192*eta14*eta22*lambda3*p3-16*eta13*eta26*mu3-144*eta13*eta25*lambda*mu2*p
        -192*eta13*eta24*lambda2*mu*p2-64*eta12*eta26*lambda*mu2*p-272*eta12*eta25*lambda2*mu*p2
        -192*eta12*eta24*lambda3*p3-80*eta1*eta26*lambda2*mu*p2-144*eta1*eta25*lambda3*p3
        -32*eta26*lambda3*p3-64*eta16*eta22*lambda*mu2-160*eta16*eta2*lambda2*mu*p-96*eta16*lambda3*p2
        -96*eta15*eta23*lambda*mu2-416*eta15*eta22*lambda2*mu*p-352*eta15*eta2*lambda3*p2
        +16*eta14*eta24*lambda*mu2-160*eta14*eta23*lambda2*mu*p-336*eta14*eta22*lambda3*p2
        +48*eta13*eta25*lambda*mu2+224*eta13*eta24*lambda2*mu*p+80*eta13*eta23*lambda3*p2
        +128*eta12*eta25*lambda2*mu*p+240*eta12*eta24*lambda3*p2+80*eta1*eta25*lambda3*p2
        +80*eta16*eta2*lambda2*mu+96*eta16*lambda3*p+144*eta15*eta22*lambda2*mu+272*eta15*eta2*lambda3*p
        +16*eta14*eta23*lambda2*mu+160*eta14*eta22*lambda3*p-48*eta13*eta24*lambda2*mu
        -80*eta13*eta23*lambda3*p-64*eta12*eta24*lambda3*p-32*eta16*lambda3-64*eta15*eta2*lambda3
        -16*eta14*eta22*lambda3+16*eta13*eta23*lambda3)*sigma8;
    aux3 = (-64*eta14*eta22*mu6-128*eta14*eta2*lambda*mu5*p-64*eta14*lambda2*mu4*p2-128*eta13*eta23*mu6
        -640*eta13*eta22*lambda*mu5*p-768*eta13*eta2*lambda2*mu4*p2-256*eta13*lambda3*mu3*p3
        -64*eta12*eta24*mu6-640*eta12*eta23*lambda*mu5*p-1408*eta12*eta22*lambda2*mu4*p2
        -768*eta12*eta2*lambda3*mu3*p3-128*eta1*eta24*lambda*mu5*p-768*eta1*eta23*lambda2*mu4*p2
        -768*eta1*eta22*lambda3*mu3*p3-64*eta24*lambda2*mu4*p2-256*eta23*lambda3*mu3*p3
        +128*eta14*eta2*lambda*mu5+128*eta14*lambda2*mu4*p+384*eta13*eta22*lambda*mu5
        +1152*eta13*eta2*lambda2*mu4*p+640*eta13*lambda3*mu3*p2+256*eta12*eta23*lambda*mu5
        +1408*eta12*eta22*lambda2*mu4*p+1408*eta12*eta2*lambda3*mu3*p2+384*eta1*eta23*lambda2*mu4*p
        +896*eta1*eta22*lambda3*mu3*p2+128*eta23*lambda3*mu3*p2-64*eta14*lambda2*mu4
        -384*eta13*eta2*lambda2*mu4-512*eta13*lambda3*mu3*p-384*eta12*eta22*lambda2*mu4
        -896*eta12*eta2*lambda3*mu3*p-64*eta12*lambda4*mu2*p2-384*eta1*eta22*lambda3*mu3*p
        -128*eta1*eta2*lambda4*mu2*p2-64*eta22*lambda4*mu2*p2+128*eta13*lambda3*mu3
        +256*eta12*eta2*lambda3*mu3+128*eta12*lambda4*mu2*p+128*eta1*eta2*lambda4*mu2*p
        -64*eta12*lambda4*mu2)*sigma2;
    aux4 = (-64*eta15*eta22*mu5-128*eta15*eta2*lambda*mu4*p-64*eta15*lambda2*mu3*p2-64*eta14*eta23*mu5
        -448*eta14*eta22*lambda*mu4*p-576*eta14*eta2*lambda2*mu3*p2-192*eta14*lambda3*mu2*p3
        +64*eta13*eta24*mu5-512*eta13*eta22*lambda2*mu3*p2-384*eta13*eta2*lambda3*mu2*p3
        +64*eta12*eta25*mu5+448*eta12*eta24*lambda*mu4*p+512*eta12*eta23*lambda2*mu3*p2
        +128*eta1*eta25*lambda*mu4*p+576*eta1*eta24*lambda2*mu3*p2+384*eta1*eta23*lambda3*mu2*p3
        +64*eta25*lambda2*mu3*p2+192*eta24*lambda3*mu2*p3+128*eta15*eta2*lambda*mu4
        +128*eta15*lambda2*mu3*p+64*eta14*eta22*lambda*mu4+512*eta14*eta2*lambda2*mu3*p
        +320*eta14*lambda3*mu2*p2-448*eta13*eta23*lambda*mu4-1344*eta13*eta22*lambda2*mu3*p
        -1216*eta13*eta2*lambda3*mu2*p2-576*eta13*lambda4*mu*p3-384*eta12*eta24*lambda*mu4
        -2368*eta12*eta23*lambda2*mu3*p-3648*eta12*eta22*lambda3*mu2*p2-1728*eta12*eta2*lambda4*mu*p3
        -640*eta1*eta24*lambda2*mu3*p-2368*eta1*eta23*lambda3*mu2*p2-1728*eta1*eta22*lambda4*mu*p3
        -256*eta24*lambda3*mu2*p2-576*eta23*lambda4*mu*p3-64*eta15*lambda2*mu3+64*eta14*eta2*lambda2*mu3
        -64*eta14*lambda3*mu2*p+960*eta13*eta22*lambda2*mu3+2432*eta13*eta2*lambda3*mu2*p
        +1408*eta13*lambda4*mu*p2+896*eta12*eta23*lambda2*mu3+3648*eta12*eta22*lambda3*mu2*p
        +3136*eta12*eta2*lambda4*mu*p2+1152*eta1*eta23*lambda3*mu2*p+2048*eta1*eta22*lambda4*mu*p2
        +320*eta23*lambda4*mu*p2-64*eta14*lambda3*mu2-832*eta13*eta2*lambda3*mu2
        -1088*eta13*lambda4*mu*p-1024*eta12*eta22*lambda3*mu2-1984*eta12*eta2*lambda4*mu*p
        -128*eta12*lambda5*p2-896*eta1*eta22*lambda4*mu*p-256*eta1*eta2*lambda5*p2-128*eta22*lambda5*p2
        +256*eta13*lambda4*mu+576*eta12*eta2*lambda4*mu+256*eta12*lambda5*p+256*eta1*eta2*lambda5*p
        -128*eta12*lambda5)*sigma4;
    aux5 = (-16*eta16*eta22*mu4-32*eta16*eta2*lambda*mu3*p-16*eta16*lambda2*mu2*p2+32*eta15*eta23*mu4
        +128*eta15*eta22*lambda*mu3*p+192*eta15*eta2*lambda2*mu2*p2+96*eta15*lambda3*mu*p3
        +96*eta14*eta24*mu4+736*eta14*eta23*lambda*mu3*p+1728*eta14*eta22*lambda2*mu2*p2
        +1536*eta14*eta2*lambda3*mu*p3+432*eta14*lambda4*p4+32*eta13*eta25*mu4
        +736*eta13*eta24*lambda*mu3*p+3040*eta13*eta23*lambda2*mu2*p2+4128*eta13*eta22*lambda3*mu*p3
        +1728*eta13*eta2*lambda4*p4-16*eta12*eta26*mu4+128*eta12*eta25*lambda*mu3*p
        +1728*eta12*eta24*lambda2*mu2*p2+4128*eta12*eta23*lambda3*mu*p3+2592*eta12*eta22*lambda4*p4
        -32*eta1*eta26*lambda*mu3*p+192*eta1*eta25*lambda2*mu2*p2+1536*eta1*eta24*lambda3*mu*p3
        +1728*eta1*eta23*lambda4*p4-16*eta26*lambda2*mu2*p2+96*eta25*lambda3*mu*p3+432*eta24*lambda4*p4
        +32*eta16*eta2*lambda*mu3+32*eta16*lambda2*mu2*p-192*eta15*eta22*lambda*mu3
        -480*eta15*eta2*lambda2*mu2*p-320*eta15*lambda3*mu*p2-512*eta14*eta23*lambda*mu3
        -2688*eta14*eta22*lambda2*mu2*p-3744*eta14*eta2*lambda3*mu*p2-1440*eta14*lambda4*p3
        -224*eta13*eta24*lambda*mu3-3040*eta13*eta23*lambda2*mu2*p-7456*eta13*eta22*lambda3*mu*p2
        -4608*eta13*eta2*lambda4*p3+64*eta12*eta25*lambda*mu3-768*eta12*eta24*lambda2*mu2*p
        -4928*eta12*eta23*lambda3*mu*p2-5184*eta12*eta22*lambda4*p3+96*eta1*eta25*lambda2*mu2*p
        -864*eta1*eta24*lambda3*mu*p2-2304*eta1*eta23*lambda4*p3+32*eta25*lambda3*mu*p2
        -288*eta24*lambda4*p3-16*eta16*lambda2*mu2+288*eta15*eta2*lambda2*mu2+352*eta15*lambda3*mu*p
        +864*eta14*eta22*lambda2*mu2+2784*eta14*eta2*lambda3*mu*p+1712*eta14*lambda4*p2
        +480*eta13*eta23*lambda2*mu2+3680*eta13*eta22*lambda3*mu*p+4128*eta13*eta2*lambda4*p2
        -96*eta12*eta24*lambda2*mu2+1152*eta12*eta23*lambda3*mu*p+3104*eta12*eta22*lambda4*p2
        -96*eta1*eta24*lambda3*mu*p+672*eta1*eta23*lambda4*p2-16*eta24*lambda4*p2-128*eta15*lambda3*mu
        -576*eta14*eta2*lambda3*mu-832*eta14*lambda4*p-416*eta13*eta22*lambda3*mu-1376*eta13*eta2*lambda4*p
        +64*eta12*eta23*lambda3*mu-512*eta12*eta22*lambda4*p+32*eta1*eta23*lambda4*p+128*eta14*lambda4
        +128*eta13*eta2*lambda4-16*eta12*eta22*lambda4)*sigma6;
    f = aux1+aux2+aux3+aux4+aux5;

    Resolvent r={a,b,c,d,e,f};
    return r;
}
