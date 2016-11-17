#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>

#include "koujdm_lib.h"

// function for reading fractions from a string and computing them with arbitrary precision
int readfraction(mpfr_t x, char * from);

// main program: executed as ./f3 digits mu sigma lambda eta1 eta2 p t b n B
int main(int argc, char *argv[]) {
    uint8_t k; // counter for reading arguments from command line

    int64_t digits; // number of exact decimal places for the results
    Parameters_mp prm;
    Parametersf3_mp prm2;

    // read digits and compute number of necessary bits to hold the results
    if (argc != 12 || sscanf(argv[1],"%"PRIi64,&digits)!=1) {
        fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
        return -1;
    }
    prm.bits = (digits+10)*log(10.)/log(2.)+16; // add 2 extra bytes to be sure

    // initialise multiprecision structures
    mpfr_inits2(prm.bits,prm.mu,prm.sigma,prm.lambda,prm.eta1,prm.eta2,prm.p,NULL);
    mpfr_init2(prm2.b,prm.bits);

    // read input arguments
    k=2; // mu
    if (strchr(argv[k],'/')==NULL) {
        if (mpfr_inp_str(prm.mu,fmemopen(argv[k],strlen(argv[k]),"r"),10,MPFRND)==0) {
            fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
            return -1;
        }
    } else {
        if (readfraction(prm.mu,argv[k])!=0) {
            fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
            return -1;
        }
    }
    k=3; // sigma
    if (strchr(argv[k],'/')==NULL) {
        if (mpfr_inp_str(prm.sigma,fmemopen(argv[k],strlen(argv[k]),"r"),10,MPFRND)==0) {
            fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
            return -1;
        }
    } else {
        if (readfraction(prm.sigma,argv[k])!=0) {
            fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
            return -1;
        }
    }
    k=4; // lambda
    if (strchr(argv[k],'/')==NULL) {
        if (mpfr_inp_str(prm.lambda,fmemopen(argv[k],strlen(argv[k]),"r"),10,MPFRND)==0) {
            fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
            return -1;
        }
    } else {
        if (readfraction(prm.lambda,argv[k])!=0) {
            fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
            return -1;
        }
    }
    k=5; // eta1
    if (strchr(argv[k],'/')==NULL) {
        if (mpfr_inp_str(prm.eta1,fmemopen(argv[k],strlen(argv[k]),"r"),10,MPFRND)==0) {
            fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
            return -1;
        }
    } else {
        if (readfraction(prm.eta1,argv[k])!=0) {
            fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
            return -1;
        }
    }
    k=6; // eta2
    if (strchr(argv[k],'/')==NULL) {
        if (mpfr_inp_str(prm.eta2,fmemopen(argv[k],strlen(argv[k]),"r"),10,MPFRND)==0) {
            fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
            return -1;
        }
    } else {
        if (readfraction(prm.eta2,argv[k])!=0) {
            fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
            return -1;
        }
    }
    k=7; // p
    if (strchr(argv[k],'/')==NULL) {
        if (mpfr_inp_str(prm.p,fmemopen(argv[k],strlen(argv[k]),"r"),10,MPFRND)==0) {
            fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
            return -1;
        }
    } else {
        if (readfraction(prm.p,argv[k])!=0) {
            fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
            return -1;
        }
    }
    k=8; // t
    if (sscanf(argv[k],"%"PRIi64,&prm2.t)==0) {
        fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
        return -1;
    }
    k=9; // b
    if (strchr(argv[k],'/')==NULL) {
        if (mpfr_inp_str(prm2.b,fmemopen(argv[k],strlen(argv[k]),"r"),10,MPFRND)==0) {
            fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
            return -1;
        }
    } else {
        if (readfraction(prm2.b,argv[k])!=0) {
            fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
            return -1;
        }
    }
    k=10; // n
    if (sscanf(argv[k],"%"PRIi64,&prm2.n)==0) {
        fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
        return -1;
    }
    k=11; // B
    if (sscanf(argv[k],"%"PRIi64,&prm2.B)==0) {
        fprintf(stderr,"%s digits mu sigma lambda eta1 eta2 p t b n B\n",argv[0]);
        return -1;
    }

    // declare and initialise result variable (complex)
    mpc_t result;
    mpc_init2(result,prm.bits);
    // run gaver method
    gaver(result,&hat_f1_mp,prm2.t,prm2.n,prm2.B,prm,(void*)&prm2);
    // print result in form (a b) where a is the real part and b is the imaginary part
    mpfr_out_str(stdout,10,digits,mpc_realref(result),MPCRND);
    printf("\n");

    return 0;
}


int readfraction(mpfr_t x, char * from) {
    // make a copy of the input string
    char *cpy;
    cpy = (char*)malloc(strlen(from)+1);
    strcpy(cpy,from);
    // tokenize it with '/' as separator
    int64_t num,denom;
    char *token = strtok(cpy,"/");
    // read the first token as numerator and the second as denominator
    if (sscanf(token,"%"PRIi64,&num)==0)
        return -1;
    token = strtok(NULL,"/");
    if (sscanf(token,"%"PRIi64,&denom)==0)
        return -1;
    // compute the multiprecision division
    mpfr_set_si(x,num,MPFRND);
    mpfr_div_si(x,x,denom,MPFRND);

    return 0;
}
