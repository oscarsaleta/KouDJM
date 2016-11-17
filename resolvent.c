#include "koujdm_lib.h"

#include <stdio.h>

int main(int argc, char *argv[]) {
    Parameters prm;

    if (argc != 7
        || sscanf(argv[1],"%Lf",&prm.mu)!=1
        || sscanf(argv[2],"%Lf",&prm.sigma)!=1
        || sscanf(argv[3],"%Lf",&prm.lambda)!=1
        || sscanf(argv[4],"%Lf",&prm.eta1)!=1
        || sscanf(argv[5],"%Lf",&prm.eta2)!=1
        || sscanf(argv[6],"%Lf",&prm.p)!=1) {
        fprintf(stderr,"%s mu sigma lambda eta1 eta2 p\n",argv[0]);
        return -1;
    }

    Resolvent r = resolvent(prm);
    printf("Resolvent = (%10.9Le) * alpha^5 +\n\t (%10.9Le) * alpha^4 +\n\t (%10.9Le) * alpha^3 +\n\t (%10.9Le) * alpha^2 +\n\t (%10.9Le) * alpha +\n\t (%10.9Le)\n",r.a,r.b,r.c,r.d,r.e,r.f);


    return 0;
}