#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "mplib.h"
#include "montgomery.h"

void pro_curve_point_test() {
    MONTG_CURVE c = (MONTG_CURVE)malloc(sizeof(MONTG_CURVE_t) * 1);
    PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    mpz_t mp_n, mp_a, mp_b, mp_x, mp_y, mp_z, mp_y2, mp_x2, mp_x3, mp_by2, mp_by2z, mp_ax2, mp_ax2z, mp_z2, mp_xz2, mp_d, mp_mod, mp_dif;
    ui_t nl;
    int i, j, res = 0, flag = 0, trues = 0, falses = 0;

    mpz_init(mp_n);
    mpz_init(mp_a);
    mpz_init(mp_b);
    mpz_init(mp_x);
    mpz_init(mp_y);
    mpz_init(mp_z);
    mpz_init(mp_y2);
    mpz_init(mp_x2);
    mpz_init(mp_x3);
    mpz_init(mp_by2);
    mpz_init(mp_by2z);
    mpz_init(mp_ax2);
    mpz_init(mp_ax2z);
    mpz_init(mp_z2);
    mpz_init(mp_xz2);
    mpz_init(mp_d);
    mpz_init(mp_mod);
    mpz_init(mp_dif);
    mpz_set_ui(mp_mod, 0L);
    mpz_set_ui(mp_dif, 0L);

    for (i = 0; i < 1000000; i++) {
        nl = (ui_t)(rand() % 100 + 1);
        ui_t n[nl], mu[nl + 1], d[nl];
        
        flag = 0;
        big_rand(n, nl);
        big_get_mu(mu, n, nl);
        for (j = 0; j < nl; j++) {
            d[j] = 0;
        }

        pro_curve_point(d, c, p, n, nl, mu, nl + 1, &flag);

        mpz_import(mp_n, nl, -1, 4, 0, 0, n);
        if(flag == 0) {
            mpz_import(mp_d, nl, -1, 4, 0, 0, d);
            mpz_mod(mp_mod, mp_n, mp_d);
            res = mpz_cmp_ui(mp_mod, 0);
            if(res == 0) {
                trues++;
            } else {
                falses++;
                printf("i = %d\n", i);
            }
        } else {
            mpz_import(mp_a, nl, -1, 4, 0, 0, c->A);
            mpz_import(mp_b, nl, -1, 4, 0, 0, c->B);
            mpz_import(mp_x, nl, -1, 4, 0, 0, p->X);
            mpz_import(mp_y, nl, -1, 4, 0, 0, p->Y);
            mpz_import(mp_z, nl, -1, 4, 0, 0, p->Z);
            mpz_mul(mp_y2, mp_y, mp_y);
            mpz_mul(mp_by2, mp_b, mp_y2);
            mpz_mul(mp_by2z, mp_by2, mp_z);
            mpz_mul(mp_x2, mp_x, mp_x);
            mpz_mul(mp_x3, mp_x2, mp_x);
            mpz_mul(mp_ax2, mp_a, mp_x2);
            mpz_mul(mp_ax2z, mp_ax2, mp_z);
            mpz_mul(mp_z2, mp_z, mp_z);
            mpz_mul(mp_xz2, mp_x, mp_z2);
            mpz_sub(mp_dif, mp_by2z, mp_x3);
            mpz_sub(mp_dif, mp_dif, mp_ax2z);
            mpz_sub(mp_dif, mp_dif, mp_xz2);
            mpz_mod(mp_dif, mp_dif, mp_n);
            res = mpz_cmp_ui(mp_dif, 0);
            if(res == 0) {
                trues++;
            } else {
                falses++;
                printf("i = %d\n", i);
            }
        }
    }
    printf("TRUE: %d\n", trues);
    printf("FALSE: %d\n", falses);
}

void aff_curve_point_test() {
    MONTG_CURVE c = (MONTG_CURVE)malloc(sizeof(MONTG_CURVE_t) * 1);
    AFF_POINT p = (AFF_POINT)malloc(sizeof(AFF_POINT_t) * 1);
    mpz_t mp_n, mp_a, mp_b, mp_x, mp_y, mp_y2, mp_x2, mp_x3, mp_by2, mp_ax2, mp_d, mp_mod, mp_dif;
    ui_t nl;
    int i, j, res = 0, flag = 0, trues = 0, falses = 0;

    mpz_init(mp_n);
    mpz_init(mp_a);
    mpz_init(mp_b);
    mpz_init(mp_x);
    mpz_init(mp_y);
    mpz_init(mp_y2);
    mpz_init(mp_x2);
    mpz_init(mp_x3);
    mpz_init(mp_by2);
    mpz_init(mp_ax2);
    mpz_init(mp_d);
    mpz_init(mp_mod);
    mpz_init(mp_dif);
    mpz_set_ui(mp_mod, 0L);
    mpz_set_ui(mp_dif, 0L);

    for (i = 0; i < 1000000; i++) {
        nl = (ui_t)(rand() % 100 + 1);
        ui_t n[nl], mu[nl + 1], d[nl];
        
        flag = 0;
        big_rand(n, nl);
        big_get_mu(mu, n, nl);
        for (j = 0; j < nl; j++) {
            d[j] = 0;
        }

        aff_curve_point(d, c, p, n, nl, mu, nl + 1, &flag);

        mpz_import(mp_n, nl, -1, 4, 0, 0, n);
        if(flag == 0) {
            mpz_import(mp_d, nl, -1, 4, 0, 0, d);
            mpz_mod(mp_mod, mp_n, mp_d);
            res = mpz_cmp_ui(mp_mod, 0);
            if(res == 0) {
                trues++;
            } else {
                falses++;
            }
        } else {
            mpz_import(mp_a, nl, -1, 4, 0, 0, c->A);
            mpz_import(mp_b, nl, -1, 4, 0, 0, c->B);
            mpz_import(mp_x, nl, -1, 4, 0, 0, p->x);
            mpz_import(mp_y, nl, -1, 4, 0, 0, p->y);
            mpz_mul(mp_y2, mp_y, mp_y);
            mpz_mul(mp_by2, mp_b, mp_y2);
            mpz_mul(mp_x2, mp_x, mp_x);
            mpz_mul(mp_x3, mp_x2, mp_x);
            mpz_mul(mp_ax2, mp_a, mp_x2);
            mpz_sub(mp_dif, mp_by2, mp_x3);
            mpz_sub(mp_dif, mp_dif, mp_ax2);
            mpz_sub(mp_dif, mp_dif, mp_x);
            mpz_mod(mp_dif, mp_dif, mp_n);
            res = mpz_cmp_ui(mp_dif, 0);
            if(res == 0) {
                trues++;
            } else {
                falses++;
            }
        }
    }
    printf("TRUE: %d\n", trues);
    printf("FALSE: %d\n", falses);
}

void pro_add_test() {
    MONTG_CURVE c = (MONTG_CURVE)malloc(sizeof(MONTG_CURVE_t) * 1);
    PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p2 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    mpz_t  mp_n, mp_x1, mp_x2, mp_z1, mp_z2, mp_a, mp_a2; 
    ui_t nl;
    int i, trues = 0, falses = 0;

    mpz_init(mp_n);
    mpz_init(mp_x1);
    mpz_init(mp_x2);
    mpz_init(mp_z1);
    mpz_init(mp_z2);
    mpz_init(mp_a);
    mpz_init(mp_a2);
    
    for (i = 0; i < 1000000; i++) {
        nl = (ui_t)(rand() % 100 + 1);
        ui_t n[nl], mu[nl + 1], X1[nl], X2[nl], Z1[nl], Z2[nl];
        ui_t a1[nl]; // temp

        mpz_set_ui(mp_a, 0L);
        mpz_set_ui(mp_a2, 0L);

        big_rand(n, nl);
        big_get_mu(mu, n, nl);
        big_rand(X1, nl);
        big_rand(Z1, nl);
        big_rand(X2, nl);
        big_rand(Z2, nl);
    
        mpz_import(mp_n, nl, -1, 4, 0, 0, n);
        mpz_import(mp_x1, nl, -1, 4, 0, 0, X1);
        mpz_import(mp_x2, nl, -1, 4, 0, 0, X2);
        mpz_import(mp_z1, nl, -1, 4, 0, 0, Z1);
        mpz_import(mp_z2, nl, -1, 4, 0, 0, Z2);

        p1->X = X1;
        p2->X = X2;
        p1->Z = Z1;
        p2->Z = Z2;

        pro_add(p, p1, p2, p, n, nl, mu, nl + 1);

        mpz_sub(mp_a, mp_x2, mp_z2);
        while(mpz_sgn(mp_a) == -1) {
            mpz_add(mp_a, mp_a, mp_n);
        }
        mpz_import(mp_a2, nl, -1, 4, 0, 0, p->X); // p->X keeps a temp
        int res = mpz_cmp(mp_a, mp_a2);
        if(res == 0) {
            trues++;
        } else {
            falses++;
        }
    }
    printf("TRUES: %d\n", trues);
    printf("FALSES: %d\n", falses);
}

int main() {
    // pro_curve_point_test();
    // aff_curve_point_test();
    pro_add_test();
}