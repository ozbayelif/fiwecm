/**
 * \file ecm.c
 * \brief Implementation of ecm.h library.
 */ 

#include <gmp.h>
#include <stdlib.h>
#include <time.h>
#include "montgomery.h"
#include "mplib.h"
#include "ecm.h"

int primes[39] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167};
int prime_powers[4] = {1,2,3,4};

void generate_n(ui n, ui_t nl) {
    int i, j, k, prime, power, exceed = 0;
    ui_t n_[2 * nl], factor[nl];

    n_[0] = 1L;
    for(i = 1; i < nl; i++) {
        n_[i] = 0L;
        factor[i] = 0L;
    }
    while(exceed != 1) {
        big_cpy(n, n_, 0, nl);
        prime = primes[rand() % PRIMESL];
        power = prime_powers[rand() % POWERSL];
        factor[0] = 1;
        for(j = 0; j < power; j++) {
            factor[0] *= prime;
        }
        big_mul(n_, n, nl, factor, nl);
        for(k = nl; k < 2 * nl; k++) {
            if(n_[k] != 0L) {
                k = 2 * nl;
                exceed = 1;
            }
        }
    }
}

int ecm(ui d, ui n, ui_t nl) {
    ui A24 = (ui)malloc(sizeof(ui_t) * nl);
    MONTG_CURVE c = (MONTG_CURVE)malloc(sizeof(MONTG_CURVE_t) * 1);
    PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    ui_t mu[nl + 1], B[1];
    int i, j, is_one, is_n, is_zero, flag;

    big_get_mu(mu, n, nl);
    for(i = 0; i <nl; i++) {
        d[i] = 0L;
    }
    for(i = 0; i < CRV_THRESHOLD; i++) {
        do {
            pro_curve_point(d, c, p1, n, nl, mu, nl + 1, &flag);
        } while(flag == -1);
        if(flag == 0) {
            return 1;
        } else {
            big_get_A24(A24, c->A, n, nl, mu, nl + 1, &flag);
            if(flag == 1) {
                B[0] = 1;
                for(j = 1; j < 11; j++) {
                    B[0] *= j;
                }
                // TODO: How to choose B?
                // TODO: How big is B?
                for(j = 0; j < B_THRESHOLD; j++) {
                    pro_ladder(p, p1, A24, B, 1, n, nl, mu, nl + 1);
                    // big_print(stdout, p->Z, nl, "Z", NULL);
                    big_is_equal_ui(&is_zero, p->Z, nl, 0L);
                    if(is_zero == 1) {
                        printf("zero\n");
                        big_gcd(d, nl, B, 1, n, nl);
                        big_is_equal_ui(&is_one, d, nl, 1L);
                        big_is_equal(&is_n, d, n, nl);
                        if(!is_one && !is_n) {
                            return 2;
                        }
                    }
                    B[0]++;
                    // TODO: How to increment B?
                }
            } else {
                return 3;
                // FIXME: Does not find in (A+2)/4!
            }
        }
    }
    return 0;
}