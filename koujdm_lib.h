#pragma once

#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>

#define MPCRND MPC_RNDNN
#define MPFRND MPFR_RNDN


#ifndef M_PI
#define M_PI acos(-1.0)
#endif

/* Structures */

// quartic polynomial
typedef struct {
    long double complex a,b,c,d,e;
} Quartic;

typedef struct {
    long double complex x1,x2,x3,x4;
} QuarticSolutions;

// resultant
typedef struct {
    long double a,b,c,d,e,f;
} Resultant;

// generic parameter structure
typedef struct {
    long double mu,sigma,lambda,eta1,eta2,p;
} Parameters;

// specific parameter structures
typedef struct {
    long double b;
    int64_t t,A,n,B;
} Parametersf1;

typedef struct {
    long double a,b;
    int64_t t,A,n,B;
} Parametersf2;


/* Multiprecision structures */

// quartic polynomial
typedef struct {
    mpc_t a,b,c,d,e;
} Quartic_mp;

typedef struct {
    mpc_t x1,x2,x3,x4;
} QuarticSolutions_mp;


// generic parameter structure
typedef struct {
    mpfr_t mu,sigma,lambda,eta1,eta2,p;
    mpfr_prec_t bits;
} Parameters_mp;

// specific parameter structures
typedef struct {
    mpfr_t b;
    int64_t t,n,B;
} Parametersf3_mp;






/* Generic functions */
// binomial
int64_t binomial(int64_t, int64_t);
// quartic polynomial
Quartic polynomial(long double complex, Parameters);
// roots of quartic polynomial
QuarticSolutions quartic_solve(Quartic);
// resultant of a quartic polynomial
Resultant resultant(Parameters);

// first block functions
long double suma(long double complex (*f)(long double complex, Parameters, void*), int64_t, int64_t, int64_t, Parameters, void*);
long double euler(long double complex (*f)(long double complex, Parameters, void*), int64_t, int64_t, int64_t, int64_t, Parameters, void*);
// f1.c
long double complex hat_f1(long double complex, Parameters, void*);
// f2.c
long double complex aux(long double complex, long double complex, long double, long double);
long double complex dG(long double complex, Parameters);
long double complex hat_f2(long double complex, Parameters, void*);

// second block functions (arbitrary precision with mpc)
Quartic_mp polynomial_mp(mpc_t, Parameters_mp);
QuarticSolutions_mp quartic_solve_mp(Quartic_mp, mpfr_prec_t);
void factorial(mpc_t, int64_t);
void fn(mpc_t, void (*f)(mpc_t, mpc_t, Parameters_mp, void *), int64_t t, int64_t n, Parameters_mp, void *);
void gaver(mpc_t, void (*f)(mpc_t, mpc_t, Parameters_mp, void *), int64_t t, int64_t n, int64_t B, Parameters_mp, void *);
void hat_f1_mp(mpc_t result, mpc_t alpha, Parameters_mp, void *);

