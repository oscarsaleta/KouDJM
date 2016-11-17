#include "koujdm_lib.h"

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

int main(int argc, char *argv[]) {
    Parameters prm;
    Parametersf2 *prm2 = (Parametersf2 *)malloc(sizeof(Parametersf2));

    if (argc != 13
        || sscanf(argv[1],"%Lf",&prm.mu)!=1
        || sscanf(argv[2],"%Lf",&prm.sigma)!=1
        || sscanf(argv[3],"%Lf",&prm.lambda)!=1
        || sscanf(argv[4],"%Lf",&prm.eta1)!=1
        || sscanf(argv[5],"%Lf",&prm.eta2)!=1
        || sscanf(argv[6],"%Lf",&prm.p)!=1
        || sscanf(argv[7],"%"PRIi64,&prm2->t)!=1
        || sscanf(argv[8],"%Lf",&prm2->a)!=1
        || sscanf(argv[9],"%Lf",&prm2->b)!=1
        || sscanf(argv[10],"%"PRIi64,&prm2->A)!=1
        || sscanf(argv[11],"%"PRIi64,&prm2->n)!=1
        || sscanf(argv[12],"%"PRIi64,&prm2->B)!=1) {
        fprintf(stderr,"%s mu sigma lambda eta1 eta2 p t a b A n B\n",argv[0]);
        return -1;
    }

    long double result = euler(&hat_f2,prm2->t,prm2->A,prm2->n,prm2->B,prm,prm2);
    printf("result=%30.29Lg\n",result);

    return 0;
}