#include <gmp.h>
#include <stdlib.h>
#include <time.h>
#include "montgomery.h"
#include "mplib.h"

void pro_curve_point(ui d, MONTG_CURVE c, PRO_POINT p, ui n, ui_t nl, ui mu, ui_t mul, int *flag) {
    ui A = (ui)malloc(sizeof(ui_t) * nl);
    ui B = (ui)malloc(sizeof(ui_t) * nl);
    ui X = (ui)malloc(sizeof(ui_t) * nl);
    ui Y = (ui)malloc(sizeof(ui_t) * nl);
    ui Z = (ui)malloc(sizeof(ui_t) * nl);
    ui_t B_[2 * nl];
    ui_t X2_[2 * nl], X2[nl];
    ui_t AX2_[2 * nl], AX2[nl];
    ui_t AX2Z_[2 * nl], AX2Z[nl];
    ui_t X3_[2 * nl], X3[nl];
    ui_t Z2_[2 * nl], Z2[nl];
    ui_t XZ2_[2 * nl], XZ2[nl];
    ui_t Y2_[2 * nl], Y2[nl];
    ui_t Y2Z_[2 * nl], Y2Z[nl];
    ui_t iY2Z[nl];
    ui_t RHS1_[2 * nl], RHS1[nl];
    ui_t RHS_[2 * nl], RHS[nl];
    mpz_t mp_n, mp_Y2Z, mp_iY2Z, mp_d;
    int i;

    mpz_init(mp_n);
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
    
    // R = X^3 + AX^2Z + XZ^2
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
        mpz_invert(mp_iY2Z, mp_Y2Z, mp_n);
        mpz_export(iY2Z, NULL, -1, 4, 0, 0, mp_iY2Z);
        big_mul(B_, RHS, nl, iY2Z, nl);
        barret_reduction(B, B_, 2 * nl, n, nl, mu, mul);

        c->A = A;
        c->B = B;
        c->n = n;
        p->X = X;
        p->Y = Y;
        p->Z = Z;
        *flag = 1;
    }
}

void aff_curve_point(MONTG_CURVE c, AFF_POINT p, ui n) {
}

void pro_add(ui X, ui Z, ui X1, ui Z1, ui X2, ui Z2, ui Xd, ui Zd, ui n) {
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