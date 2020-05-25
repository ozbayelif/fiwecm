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

int primes[21] = {5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83};
int prime_powers[5] = {1,2,3,4,5};

void generate_B_smooth(ui z, ui_t l) {
    int i, j, k, prime, power, exceed = 0;
    ui_t z_[2 * l], factor[l];

    z_[0] = 1L;
    for(i = 1; i < l; i++) {
        z_[i] = 0L;
        factor[i] = 0L;
    }
    while(exceed != 1) {
        big_cpy(z, z_, 0, l);
        prime = primes[rand() % PRIMESL];
        power = prime_powers[rand() % POWERSL];
        factor[0] = 1L;
        for(j = 0; j < power; j++) {
            factor[0] *= prime;
        }
        big_mul(z_, z, l, factor, l);
        for(k = l; k < 2 * l; k++) {
            if(z_[k] != 0L) {
                k = 2 * l;
                exceed = 1;
            }
        }
    }
}

int ecm(ui d, ui n, ui_t nl) {
    MONTG_CURVE_t c; PRO_POINT_t p, p1;
    ui_t A24[nl], mu[nl + 1], k[1];
    int i, j, is_one, is_n, is_zero, flag;

    big_get_mu(mu, n, nl);
    for(i = 0; i < nl; i++) {
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
                generate_B_smooth(k, 1);
                for(j = 0; j < k_THRESHOLD; j++) {
                    pro_ladder(p, p1, A24, k, 1, n, nl, mu, nl + 1);
                    big_is_equal_ui(&is_zero, p->Z, nl, 0L);
                    if(is_zero == 1) {
                        big_gcd(d, nl, k, 1, n, nl);
                        big_is_equal_ui(&is_one, d, nl, 1L);
                        big_is_equal(&is_n, d, n, nl);
                        if(!is_one && !is_n) {
                            return 2;
                        }
                    }
                    generate_B_smooth(k, 1);
                }
            } else {
                return 3;
            }
        }
    }
    return 0;
}