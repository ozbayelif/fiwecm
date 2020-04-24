/**
 * \file test.c
 * \brief Implementation of test.h library.
 */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "mplib.h"
#include "montgomery.h"
#include "ecm.h"
#include "test.h"

void pro_curve_point_gmp_test() {
    MONTG_CURVE c = (MONTG_CURVE)malloc(sizeof(MONTG_CURVE_t) * 1);
    PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    mpz_t mp_n, mp_a, mp_b, mp_x, mp_y, mp_z, mp_y2, mp_x2, mp_x3, mp_by2, mp_by2z, mp_ax2, mp_ax2z, mp_z2, mp_xz2, mp_d, mp_mod, mp_dif;
    ui_t nl;
    int i, j, res = 0, flag = 0, trues = 0, falses = 0, singulars = 0;

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

    for (i = 0; i < 10000; i++) {
        nl = (ui_t)(rand() % 100 + 1);
        ui_t n[nl], mu[nl + 1], d[nl];
        
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
                printf("False at index: %d\n", i);
            }
        } else if(flag == 1){
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
                printf("False at index: %d\n", i);
            }
        } else {
            singulars++;
            printf("Singularity at index: %d\n", i);
        }
    }
    printf("TRUE: %d\n", trues);
    printf("FALSE: %d\n", falses);
    printf("SINGULAR: %d\n", singulars);
}

void aff_curve_point_gmp_test() {
    MONTG_CURVE c = (MONTG_CURVE)malloc(sizeof(MONTG_CURVE_t) * 1);
    AFF_POINT p = (AFF_POINT)malloc(sizeof(AFF_POINT_t) * 1);
    mpz_t mp_n, mp_a, mp_b, mp_x, mp_y, mp_y2, mp_x2, mp_x3, mp_by2, mp_ax2, mp_d, mp_mod, mp_dif;
    ui_t nl;
    int i, j, res = 0, flag = 0, trues = 0, falses = 0, singulars = 0;

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

    for (i = 0; i < 10000; i++) {
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
                printf("False at index: %d\n", i);
            }
        } else if(flag == 1){
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
                printf("False at index: %d\n", i);
            }
        } else {
            singulars++;
            printf("Singularity at index: %d\n", i);
        }
    }
    printf("TRUE: %d\n", trues);
    printf("FALSE: %d\n", falses);
    printf("SINGULAR: %d\n", singulars);
}

void pro_add_gmp_test() {
    MONTG_CURVE c = (MONTG_CURVE)malloc(sizeof(MONTG_CURVE_t) * 1);
    PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p2 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT pd = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    mpz_t  mp_n, mp_x1, mp_x2, mp_xd, mp_z1, mp_z2, mp_zd, mp_a, mp_b, mp_c, mp_d, mp_da, mp_cb, mp_e, mp_f, mp_e2, mp_f2, mp_g, mp_h, mp_X, mp_Z; 
    ui_t nl;
    int i, trues = 0, falses = 0;

    mpz_init(mp_n);
    mpz_init(mp_x1);
    mpz_init(mp_x2);
    mpz_init(mp_xd);
    mpz_init(mp_z1);
    mpz_init(mp_z2);
    mpz_init(mp_zd);
    mpz_init(mp_a);
    mpz_init(mp_b);
    mpz_init(mp_c);
    mpz_init(mp_d);
    mpz_init(mp_da);
    mpz_init(mp_cb);
    mpz_init(mp_e);
    mpz_init(mp_f);
    mpz_init(mp_e2);
    mpz_init(mp_f2);
    mpz_init(mp_g);
    mpz_init(mp_h);
    mpz_init(mp_X);
    mpz_init(mp_Z);
    
    for (i = 0; i < 10000; i++) {
        nl = (ui_t)(rand() % 100 + 1);
        ui_t n[nl], mu[nl + 1], X1[nl], X2[nl], Xd[nl], Z1[nl], Z2[nl], Zd[nl];

        mpz_set_ui(mp_X, 0L);
        mpz_set_ui(mp_Z, 0L);

        big_rand(n, nl);
        big_get_mu(mu, n, nl);
        big_mod_rand(X1, nl, n, nl, mu, nl + 1);
        big_mod_rand(X2, nl, n, nl, mu, nl + 1);
        big_mod_rand(Xd, nl, n, nl, mu, nl + 1);
        big_mod_rand(Z1, nl, n, nl, mu, nl + 1);
        big_mod_rand(Z2, nl, n, nl, mu, nl + 1);
        big_mod_rand(Zd, nl, n, nl, mu, nl + 1);
    
        mpz_import(mp_n, nl, -1, 4, 0, 0, n);
        mpz_import(mp_x1, nl, -1, 4, 0, 0, X1);
        mpz_import(mp_x2, nl, -1, 4, 0, 0, X2);
        mpz_import(mp_xd, nl, -1, 4, 0, 0, Xd);
        mpz_import(mp_z1, nl, -1, 4, 0, 0, Z1);
        mpz_import(mp_z2, nl, -1, 4, 0, 0, Z2);
        mpz_import(mp_zd, nl, -1, 4, 0, 0, Zd);

        p1->X = X1;
        p2->X = X2;
        pd->X = Xd;
        p1->Z = Z1;
        p2->Z = Z2;
        pd->Z = Zd;

        pro_add(p, p1, p2, pd, n, nl, mu, nl + 1);
        mpz_import(mp_X, nl, -1, 4, 0, 0, p->X);
        mpz_import(mp_Z, nl, -1, 4, 0, 0, p->Z);

        mpz_add(mp_a, mp_x2, mp_z2);
        mpz_mod(mp_a, mp_a, mp_n);
        mpz_sub(mp_b, mp_x2, mp_z2);
        while(mpz_sgn(mp_b) == -1) {
            mpz_add(mp_b, mp_b, mp_n);
        }
        mpz_add(mp_c, mp_x1, mp_z1);
        mpz_mod(mp_c, mp_c, mp_n);
        mpz_sub(mp_d, mp_x1, mp_z1);
        while(mpz_sgn(mp_d) == -1) {
            mpz_add(mp_d, mp_d, mp_n);
        }
        mpz_mul(mp_da, mp_d, mp_a);
        mpz_mod(mp_da, mp_da, mp_n);
        mpz_mul(mp_cb, mp_c, mp_b);
        mpz_mod(mp_cb, mp_cb, mp_n);
        mpz_add(mp_e, mp_da, mp_cb);
        mpz_mod(mp_e, mp_e, mp_n);
        mpz_sub(mp_f, mp_da, mp_cb);
        while(mpz_sgn(mp_f) == -1) {
            mpz_add(mp_f, mp_f, mp_n);
        }
        mpz_mul(mp_e2, mp_e, mp_e);
        mpz_mod(mp_e2, mp_e2, mp_n);
        mpz_mul(mp_f2, mp_f, mp_f);
        mpz_mod(mp_f2, mp_f2, mp_n);
        mpz_mul(mp_g, mp_zd, mp_e2);
        mpz_mod(mp_g, mp_g, mp_n);
        mpz_mul(mp_h, mp_xd, mp_f2);
        mpz_mod(mp_h, mp_h, mp_n);

        int res = mpz_cmp(mp_g, mp_X) + mpz_cmp(mp_h, mp_Z);
        if(res == 0) {
            trues++;
        } else {
            falses++;
        }
    }
    printf("TRUE: %d\n", trues);
    printf("FALSE: %d\n", falses);
}

void pro_add_magma_test() {
    FILE *fp = fopen("/home/ozbayelif/Development/FIWE/ecm/pro_add_test.magma", "a");
    PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p2 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT pd = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    ui_t nl;
    int i, trues = 0, falses = 0;

    fprintf(fp, "clear;\n");
    fprintf(fp, "/****************************************************************************/\n");
    fprintf(fp, "ADDM:=function(X1, Z1, X2, Z2, Xd, Zd, n)\n");
    fprintf(fp, "return (Zd * (((X1 - Z1) * (X2 + Z2)) + ((X1 + Z1) * (X2 - Z2)))^2) mod n, (Xd * (((X1 - Z1) * (X2 + Z2)) - ((X1 + Z1) * (X2 - Z2)))^2) mod n;\n");
    fprintf(fp, "end function;\n");
    fprintf(fp, "trues := 0;\n");
    fprintf(fp, "falses := 0;\n");

    for (i = 0; i < 10000; i++) {
        nl = (ui_t)(rand() % 100 + 1);
        ui_t n[nl], mu[nl + 1], X1[nl], X2[nl], Xd[nl], Z1[nl], Z2[nl], Zd[nl];

        big_rand(n, nl);
        big_get_mu(mu, n, nl);
        big_mod_rand(X1, nl, n, nl, mu, nl + 1);
        big_mod_rand(X2, nl, n, nl, mu, nl + 1);
        big_mod_rand(Xd, nl, n, nl, mu, nl + 1);
        big_mod_rand(Z1, nl, n, nl, mu, nl + 1);
        big_mod_rand(Z2, nl, n, nl, mu, nl + 1);
        big_mod_rand(Zd, nl, n, nl, mu, nl + 1);

        p1->X = X1;
        p2->X = X2;
        pd->X = Xd;
        p1->Z = Z1;
        p2->Z = Z2;
        pd->Z = Zd;

        pro_add(p, p1, p2, pd, n, nl, mu, nl + 1);

        big_print(fp, n, nl, "n", NULL);
        big_print(fp, p1->X, nl, "X1", NULL);
        big_print(fp, p1->Z, nl, "Z1", NULL);
        big_print(fp, p2->X, nl, "X2", NULL);
        big_print(fp, p2->Z, nl, "Z2", NULL);
        big_print(fp, pd->X, nl, "Xd", NULL);
        big_print(fp, pd->Z, nl, "Zd", NULL);

        big_print(fp, p->X, nl, "pX", NULL);
        big_print(fp, p->Z, nl, "pZ", NULL);

        fprintf(fp, "qX, qZ := ADDM(X1, Z1, X2, Z2, Xd, Zd, n);\n");
        fprintf(fp, "if qX eq pX then\n");
        fprintf(fp, "trues := trues + 1;\n");
        fprintf(fp, "else\n");
        fprintf(fp, "falses := falses + 1;\n");
        fprintf(fp, "end if;\n");
        
    }
    fprintf(fp, "Write(\"/home/ozbayelif/Development/FIWE/ecm/pro_add_results.magma\", trues);\n");
    fprintf(fp, "Write(\"/home/ozbayelif/Development/FIWE/ecm/pro_add_results.magma\", falses);\n");

    fclose(fp);
}

void pro_dbl_magma_test() {
    FILE *fp = fopen("/home/ozbayelif/Development/FIWE/ecm/pro_dbl_test.magma", "a");
    PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    ui_t nl;
    int i, trues = 0, falses = 0;

    fprintf(fp, "clear;\n");
    fprintf(fp, "/****************************************************************************/\n");
    fprintf(fp, "DBLM:=function(X1, Z1, A24)\n");
    fprintf(fp, "return ((X1+Z1)*(X1+Z1))*((X1-Z1)*(X1-Z1)),(((X1+Z1)*(X1+Z1))-((X1-Z1)*(X1-Z1)))*(((X1-Z1)*(X1-Z1))+A24*(((X1+Z1)*(X1+Z1)) - ((X1-Z1)*(X1-Z1))));\n");
    fprintf(fp, "end function;\n");
    fprintf(fp, "trues := 0;\n");
    fprintf(fp, "falses := 0;\n");

    for (i = 0; i < 10000; i++) {
        nl = (ui_t)(rand() % 100 + 1);
        ui_t n[nl], mu[nl + 1], X1[nl], Z1[nl], A24[nl];

        big_rand(n, nl);
        big_get_mu(mu, n, nl);
        big_mod_rand(X1, nl, n, nl, mu, nl + 1);
        big_mod_rand(Z1, nl, n, nl, mu, nl + 1);
        big_mod_rand(A24, nl, n, nl, mu, nl + 1);

        p1->X = X1;
        p1->Z = Z1;

        pro_dbl(p, p1, A24, n, nl, mu, nl + 1);

        big_print(fp, n, nl, "n", NULL);
        fprintf(fp, "R:=Integers(n);\n");
        big_print(fp, p1->X, nl, "X1", "R");
        big_print(fp, p1->Z, nl, "Z1", "R");
        big_print(fp, A24, nl, "A24", "R");

        big_print(fp, p->X, nl, "pX", "R");
        big_print(fp, p->Z, nl, "pZ", "R");

        fprintf(fp, "qX, qZ := DBLM(X1, Z1, A24);\n");
        fprintf(fp, "if qX eq pX then\n");
        fprintf(fp, "trues := trues + 1;\n");
        fprintf(fp, "else\n");
        fprintf(fp, "falses := falses + 1;\n");
        fprintf(fp, "end if;\n");
        
    }
    fprintf(fp, "Write(\"/home/ozbayelif/Development/FIWE/ecm/pro_dbl_results.magma\", trues);\n");
    fprintf(fp, "Write(\"/home/ozbayelif/Development/FIWE/ecm/pro_dbl_results.magma\", falses);\n");

    fclose(fp);
}

void pro_ladder_gmp_test() {
    MONTG_CURVE c = (MONTG_CURVE)malloc(sizeof(MONTG_CURVE_t) * 1);
    PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p2 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p3 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p4 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p5 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    int i, j, nl, kl, ll, res1 ,res2, flag, true = 0, false = 0;
    mpz_t mp_p3X, mp_p3Z, mp_p5X, mp_p5Z;
    // nl = (ui_t)(rand() % 1 + 1);
    nl = 1;

    mpz_init(mp_p3X);
    mpz_init(mp_p3Z);
    mpz_init(mp_p5X);
    mpz_init(mp_p5Z);

    mpz_set_ui(mp_p3X, 0L);
    mpz_set_ui(mp_p3Z, 0L);
    mpz_set_ui(mp_p5X, 0L);
    mpz_set_ui(mp_p5Z, 0L);

    for (i = 0; i < 10; i++) {
        kl = 1, ll = 1;
        ui_t n[nl], mu[nl + 1], A24[nl];
        ui d = (ui)malloc(sizeof(ui_t) * nl);
        ui k = (ui)malloc(sizeof(ui_t) * kl);
        ui l = (ui)malloc(sizeof(ui_t) * ll);

        // TODO: Work with prime n
        big_rand(n, nl);
        big_get_mu(mu, n, nl);
        big_rand(k, kl);
        big_rand(l, ll);

        do {
            pro_curve_point(d, c, p1, n, nl, mu, nl + 1, &flag);
        } while(flag != 1);
        do {
            big_get_A24(A24, c->A, n, nl, mu, nl + 1, &flag);
        } while(flag != 1);

        pro_ladder(p2, p1, A24, k, kl, n, nl, mu, nl + 1);  // p2 = k*P
        pro_ladder(p3, p2, A24, l, ll, n, nl, mu, nl + 1);  // p3 = l*(k*P) 

        pro_ladder(p4, p1, A24, l, ll, n, nl, mu, nl + 1);  // p4 = l*P
        pro_ladder(p5, p4, A24, k, kl, n, nl, mu, nl + 1);  // p5 = k*(l*P)

        mpz_import(mp_p3X, nl, -1, 4, 0, 0, p3->X);
        mpz_import(mp_p5X, nl, -1, 4, 0, 0, p5->X);
        mpz_import(mp_p3Z, nl, -1, 4, 0, 0, p3->Z);
        mpz_import(mp_p5Z, nl, -1, 4, 0, 0, p5->Z);

        if((mpz_cmp(mp_p3X, mp_p5X) == 0) && (mpz_cmp(mp_p3Z, mp_p5Z) == 0)) {
            true++;
        } else {
            false++;
        }
    }
    printf("TRUE: %d\n", true);
    printf("FALSE: %d\n", false);
}

void pro_ladder_magma_test() {
    // FILE *fp = fopen("/home/ozbayelif/Development/FIWE/ecm/pro_ladder_test.magma", "a");
    FILE *fp = stdout;
    MONTG_CURVE c = (MONTG_CURVE)malloc(sizeof(MONTG_CURVE_t) * 1);
    PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    AFF_POINT ap1 = (AFF_POINT)malloc(sizeof(AFF_POINT_t) * 1);
    PRO_POINT pp1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    ui_t nl, kl;
    int i, j, flag, trues = 0, falses = 0;

    fprintf(fp, "trues := 0;\n");
    fprintf(fp, "falses := 0;\n");

    fprintf(fp, "Fp:=GF(2^255-19);\n");

    for (i = 0; i < 1; i++) {
        nl = (ui_t)(rand() % 1 + 1);
        kl = 1;
        ui_t n[nl], mu[nl + 1], X1[nl], Z1[nl], k[kl], c_1[nl];
        ui d = (ui)malloc(sizeof(ui_t) * nl);
        ui A24 = (ui)malloc(sizeof(ui_t) * nl);
        flag = 0;

        c_1[0] = 1L;
        for(j = 1; j < nl; j++) {
            c_1[j] = 0L;
        }
        big_rand(n, nl);
        n[0]--;
        big_get_mu(mu, n, nl);

        do {
            aff_curve_point(d, c, ap1, n, nl, mu, nl + 1, &flag);
        } while(flag != 1);
        pp1->X = ap1->x;
        pp1->Y = ap1->y;
        pp1->Z = c_1;

        do {
            big_get_A24(A24, c->A, n, nl, mu, nl + 1, &flag);
        } while(flag != 1);

        big_print(fp, n, nl, "n", NULL);
        fprintf(fp, "R:=RingOfIntegers(n);\n");
        big_print(fp, k, kl, "k", NULL);
        big_print(fp, c->A, nl, "A", NULL);
        big_print(fp, c->B, nl, "B", NULL);
        big_print(fp, pp1->X, nl, "X1", "R");
        big_print(fp, pp1->Y, nl, "Y1", "R");
        big_print(fp, pp1->Z, nl, "Z1", "R");

        fprintf(fp, "E:=EllipticCurve([A,B]);\n");
        fprintf(fp, "P:=E![X1,Y1];\n");

        pro_ladder(p, pp1, A24, k, kl, n, nl, mu, nl + 1);

        big_print(fp, p->X, nl, "X", "R");
        big_print(fp, p->Z, nl, "Z", "R");

        // fprintf(fp, "P2 := k*P1;\n");
        // fprintf(fp, "if (P2[1] eq X and P2[3] eq Z) then\n");
        // fprintf(fp, "trues := trues + 1;\n");
        // fprintf(fp, "else\n");
        // fprintf(fp, "falses := falses + 1;\n");
        // fprintf(fp, "end if;\n\n");
    }
    // fprintf(fp, "Write(\"/home/ozbayelif/Development/FIWE/ecm/pro_ladder_results.magma\", trues);\n");
    // fprintf(fp, "Write(\"/home/ozbayelif/Development/FIWE/ecm/pro_ladder_results.magma\", falses);\n");

    fclose(fp);
}

void ecm_test() {
    int i, nl, res, true = 0, false = 0, success = 0, fail = 0, success_type[3] = {0};
    mpz_t mp_n, mp_d, mp_mod;

    mpz_init(mp_n);
    mpz_init(mp_mod);
    mpz_init(mp_d);
    mpz_set_ui(mp_d, 0L);

    for(i = 0; i < 50; i++) {
        nl = (ui_t)(rand() % 2 + 1);
        ui_t n[nl], d[nl];
        big_rand(n, nl);
        n[0]--; // To make n most probably odd
        int ret = ecm(d, n, nl);
        if(ret) {
            success++;
            success_type[ret - 1]++;
            mpz_import(mp_n, nl, -1, 4, 0, 0, n);
            mpz_import(mp_d, nl, -1, 4, 0, 0, d);
            mpz_mod(mp_mod, mp_n, mp_d);
            res = mpz_cmp_ui(mp_mod, 0);
            if(res == 0) {
                gmp_printf("n:=%Zd;\n", mp_n);
                gmp_printf("d:=%Zd;\n", mp_d);
                printf("n mod d eq 0;\n");
                true++;
            } else {
                false++;
                printf("False at index: %d\n", i);
            }
        } else {
            fail++;
        }
    }
    printf("TRUE: %d\n", true);
    printf("FALSE: %d\n", false);
    printf("SUCCESS: %d\n", success);
    printf("FAIL: %d\n", fail);
    printf("FOUND IN pro_curve_point: %d\n", success_type[0]);
    printf("FOUND IN ladder: %d\n", success_type[1]);
    printf("FOUND IN A24: %d\n", success_type[2]);
}