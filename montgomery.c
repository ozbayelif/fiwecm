#include <gmp.h>
#include <stdlib.h>
#include <time.h>
#include "montgomery.h"
#include "mplib.h"

void pro_curve_point(ui d, MONTG_CURVE c, PRO_POINT p, ui n, ui_t nl, ui mu, ui_t mul, int *flag) {
    ui_t A2_[2 * nl], A2[nl], B_[2 * nl], X2_[2 * nl], X2[nl], AX2_[2 * nl], AX2[nl], AX2Z_[2 * nl], AX2Z[nl], X3_[2 * nl], X3[nl], Z2_[2 * nl], Z2[nl], XZ2_[2 * nl], XZ2[nl], Y2_[2 * nl], Y2[nl], Y2Z_[2 * nl], Y2Z[nl], iY2Z[nl], RHS1_[2 * nl], RHS1[nl], RHS_[2 * nl], RHS[nl];
    ui A = (ui)malloc(sizeof(ui_t) * nl);
    ui B = (ui)malloc(sizeof(ui_t) * nl);
    ui X = (ui)malloc(sizeof(ui_t) * nl);
    ui Y = (ui)malloc(sizeof(ui_t) * nl);
    ui Z = (ui)malloc(sizeof(ui_t) * nl);
    mpz_t mp_n, mp_a, mp_a2, mp_b, mp_Y2Z, mp_iY2Z, mp_d;
    int i, s = 1;

    mpz_init(mp_n);
    mpz_init(mp_b);
    mpz_init(mp_Y2Z);
    mpz_init(mp_iY2Z);
    mpz_init(mp_d);
    mpz_set_ui(mp_iY2Z, 0L);
    mpz_set_ui(mp_d, 0L);

    mpz_import(mp_n, nl, -1, 4, 0, 0, n);
    big_rand(A, nl);
    big_rand(X, nl);
    big_rand(Y, nl);
    big_rand(Z, nl);
    
    big_mul(A2_, A, nl, A, nl);
    barret_reduction(A2, A2_, 2 * nl, n, nl, mu, mul);
    if(A2[0] == 4L) {
        s = 0;
    }
    for(i = 1; i < nl; i++) {
        if(A2[i] != 0L) {
            s = 1;
        }
    }
    if(s != 0) {
        // AX^2Z
        big_mul(X2_, X, nl, X, nl);
        barret_reduction(X2, X2_, 2 * nl, n, nl, mu, mul);
        big_mul(AX2_, A, nl, X2, nl);
        barret_reduction(AX2, AX2_, 2 * nl, n, nl, mu, mul);
        big_mul(AX2Z_, AX2, nl, Z, nl);
        barret_reduction(AX2Z, AX2Z_, 2 * nl, n, nl, mu, mul);

        // X^3
        big_mul(X3_, X2, nl, X, nl);
        barret_reduction(X3, X3_, 2 * nl, n, nl, mu, mul);

        // XZ^2
        big_mul(Z2_, Z, nl, Z, nl);
        barret_reduction(Z2, Z2_, 2 * nl, n, nl, mu, mul);
        big_mul(XZ2_, X, nl, Z2, nl);
        barret_reduction(XZ2, XZ2_, 2 * nl, n, nl, mu, mul);
        
        // Y^2Z
        big_mul(Y2_, Y, nl, Y, nl);
        barret_reduction(Y2, Y2_, 2 * nl, n, nl, mu, mul);
        big_mul(Y2Z_, Y2, nl, Z, nl);
        barret_reduction(Y2Z, Y2Z_, 2 * nl, n, nl, mu, mul);
        
        // RHS = X^3 + AX^2Z + XZ^2
        for(i = nl + 1; i < 2 * nl; i++) {
            RHS1_[i] = 0L;
            RHS_[i] = 0L;
        }
        big_add(RHS1_, X3, nl, AX2Z, nl);
        barret_reduction(RHS1, RHS1_, 2 * nl, n, nl, mu, mul);
        big_add(RHS_, RHS1, nl, XZ2, nl);
        barret_reduction(RHS, RHS_, 2 * nl, n, nl, mu, mul);

        // B or factorize
        mpz_import(mp_Y2Z, nl, -1, 4, 0, 0, Y2Z);
        mpz_gcd(mp_d, mp_Y2Z, mp_n);
        mpz_export(d, NULL, -1, 4, 0, 0, mp_d);

        if(mpz_cmp_ui(mp_d, 1) == 0 || mpz_cmp(mp_d, mp_n) == 0) {
            for(i = 0; i < nl; i++) { // mpz_export does not write anything
                iY2Z[i] = 0L;           // if the number is 0
            }
            mpz_invert(mp_iY2Z, mp_Y2Z, mp_n);
            mpz_export(iY2Z, NULL, -1, 4, 0, 0, mp_iY2Z);
            big_mul(B_, RHS, nl, iY2Z, nl);
            barret_reduction(B, B_, 2 * nl, n, nl, mu, mul);
            mpz_import(mp_b, nl, -1, 4, 0, 0, B);
            s = mpz_cmp_ui(mp_b, 0L);
            if(s != 0) {
                c->A = A;
                c->B = B;
                c->n = n;
                p->X = X;
                p->Y = Y;
                p->Z = Z;
                *flag = 1;
            }
        }
    } else {
        *flag = -1;
    }
}

void aff_curve_point(ui d, MONTG_CURVE c, AFF_POINT p, ui n, ui_t nl, ui mu, ui_t mul, int *flag) {
    ui_t A2_[2 * nl], A2[nl], B_[2 * nl], x2_[2 * nl], x2[nl], x3_[2 * nl], x3[nl], Ax2_[2 * nl], Ax2[nl], y2_[2 * nl], y2[nl], iy2[nl], rhs1_[2 * nl], rhs1[nl], rhs_[2 * nl], rhs[nl];
    ui A = (ui)malloc(sizeof(ui_t) * nl);
    ui B = (ui)malloc(sizeof(ui_t) * nl);
    ui x = (ui)malloc(sizeof(ui_t) * nl);
    ui y = (ui)malloc(sizeof(ui_t) * nl);
    mpz_t mp_n, mp_b, mp_y2, mp_iy2, mp_d;
    int i, s = 1;

    mpz_init(mp_n);
    mpz_init(mp_b);
    mpz_init(mp_y2);
    mpz_init(mp_iy2);
    mpz_init(mp_d);
    mpz_set_ui(mp_iy2, 0L);
    mpz_set_ui(mp_d, 0L);

    mpz_import(mp_n, nl, -1, 4, 0, 0, n);
    big_rand(A, nl);
    big_rand(x, nl);
    big_rand(y, nl);

    big_mul(A2_, A, nl, A, nl);
    barret_reduction(A2, A2_, 2 * nl, n, nl, mu, mul);
    if(A2[0] == 4L) {
        s = 0;
    }
    for(i = 1; i < nl; i++) {
        if(A2[i] != 0L) {
            s = 1;
        }
    }
    if(s != 0) {
        // Ax^2
        big_mul(x2_, x, nl, x, nl);
        barret_reduction(x2, x2_, 2 * nl, n, nl, mu, mul);
        big_mul(Ax2_, A, nl, x2, nl);
        barret_reduction(Ax2, Ax2_, 2 * nl, n, nl, mu, mul);

        // x^3
        big_mul(x3_, x2, nl, x, nl);
        barret_reduction(x3, x3_, 2 * nl, n, nl, mu, mul);

        // y^2
        big_mul(y2_, y, nl, y, nl);
        barret_reduction(y2, y2_, 2 * nl, n, nl, mu, mul);
        
        // rhs = x^3 + Ax^2 + X
        for(i = nl + 1; i < 2 * nl; i++) {
            rhs1_[i] = 0L;
            rhs_[i] = 0L;
        }
        big_add(rhs1_, x3, nl, Ax2, nl);
        barret_reduction(rhs1, rhs1_, 2 * nl, n, nl, mu, mul);
        big_add(rhs_, rhs1, nl, x, nl);
        barret_reduction(rhs, rhs_, 2 * nl, n, nl, mu, mul);

        // B or factorize
        mpz_import(mp_y2, nl, -1, 4, 0, 0, y2);
        mpz_gcd(mp_d, mp_y2, mp_n);
        mpz_export(d, NULL, -1, 4, 0, 0, mp_d);

        if(mpz_cmp_ui(mp_d, 1) == 0 || mpz_cmp(mp_d, mp_n) == 0) {
            for(i = 0; i < nl; i++) { // mpz_export does not write anything
                iy2[i] = 0L;           // if the number is 0
            }
            mpz_invert(mp_iy2, mp_y2, mp_n);
            mpz_export(iy2, NULL, -1, 4, 0, 0, mp_iy2);
            big_mul(B_, rhs, nl, iy2, nl);
            barret_reduction(B, B_, 2 * nl, n, nl, mu, mul);
            mpz_import(mp_b, nl, -1, 4, 0, 0, B);
            s = mpz_cmp_ui(mp_b, 0L);
            if(s != 0) {
                c->A = A;
                c->B = B;
                c->n = n;
                p->x = x;
                p->y = y;
                *flag = 1;
            }
        }
    } else {
        *flag = -1;
    }
}

void pro_add(PRO_POINT p, PRO_POINT p1, PRO_POINT p2, PRO_POINT pd, ui n, ui_t nl, ui mu, ui_t mul) {
    ui_t a_[2 * nl], a[nl], b[nl], bs[nl + 1], c_[2 * nl], c[nl], d[nl], d2[nl + 1], da_[2 * nl], da[nl], cb_[2 * nl], cb[nl], e_[2 * nl], e[nl], f[nl], fs[nl + 1], e2_[2 * nl], e2[nl], f2_[2 * nl], f2[nl], X_[2 * nl], Z_[2 * nl];
    ui X = (ui)malloc(sizeof(ui_t) * nl);
    ui Z = (ui)malloc(sizeof(ui_t) * nl);
    int i, s;

    for(i = nl + 1; i < 2 * nl; i++) {
        a_[i] = 0L;
        c_[i] = 0L;
        da_[i] = 0L;
        cb_[i] = 0L;
        e_[i] = 0L;
        e2_[i] = 0L;
        f2_[i] = 0L;
        X_[i] = 0L;
        Z_[i] = 0L;
    }

    big_add(a_, p2->X, nl, p2->Z, nl);
    barret_reduction(a, a_, 2 * nl, n, nl, mu, mul);        // a = X2 + Z2
    s = big_sub(b, p2->X, nl, p2->Z, nl);
    if(s) {
        do {
            big_add(bs, b, nl, n, nl);
            big_cpy(b, bs, 0, nl);
        } while((bs[nl] == 0));
    }                                                       // b = X2 - Z2
    big_add(c_, p1->X, nl, p1->Z, nl);
    barret_reduction(c, c_, 2 * nl, n, nl, mu, mul);        // c = X1 + Z1
    s = big_sub(d, p1->X, nl, p1->Z, nl);
    if(s) {
        do {
            big_add(d2, d, nl, n, nl);
            big_cpy(d, d2, 0, nl);
        } while((d2[nl] == 0));
    }                                                       // d = X1 - Z1
    big_mul(da_, d, nl, a, nl);
    barret_reduction(da, da_, 2 * nl, n, nl, mu, mul);      // da = d * a
    big_mul(cb_, c, nl, b, nl);
    barret_reduction(cb, cb_, 2 * nl, n, nl, mu, mul);      // cb = c * b
    big_add(e_, da, nl, cb, nl);
    barret_reduction(e, e_, 2 * nl, n, nl, mu, mul);        // e = da + cb
    s = big_sub(f, da, nl, cb, nl);
    if(s) {
        do {
            big_add(fs, f, nl, n, nl);
            big_cpy(f, fs, 0, nl);
        } while((fs[nl] == 0));
    }                                                       // f = da - cb
    big_mul(e2_, e, nl, e, nl);
    barret_reduction(e2, e2_, 2 * nl, n, nl, mu, mul);      // e2 = e^2
    big_mul(f2_, f, nl, f, nl);
    barret_reduction(f2, f2_, 2 * nl, n, nl, mu, mul);      // f2 = f^2
    big_mul(X_, pd->Z, nl, e2, nl);
    barret_reduction(X, X_, 2 * nl, n, nl, mu, mul);        // X = Zd * e2
    big_mul(Z_, pd->X, nl, f2, nl);
    barret_reduction(Z, Z_, 2 * nl, n, nl, mu, mul);        // Z = Xd * f2

    p->X = X;
    p->Z = Z;
}

// void aff_add(ui x, ui y, ui x1, ui y1, ui x2, ui y2, ui ..)

void pro_dbl(ui X, ui Z, ui X1, ui Z1, ui A24, ui n) {
}

void aff_dbl(ui x, ui z, ui x1, ui y1, ui A, ui B, ui n) {
}

void pro_ladder(ui X, ui Z, ui X1, ui Z1, ui k, ui n) {
}

void aff_ladder(ui x, ui y, ui x1, ui y1, ui k, ui n) {
}

int pro_is_on_curve(ui A, ui B, ui X, ui Y, ui Z, ui n) {
}

int aff_is_on_curve(ui A, ui B, ui x, ui y, ui n) {
}