#include <gmp.h>
#include <stdlib.h>
#include "montgomery.h"
#include "mplib.h"
// B*Y^2*Z eq X^3 + A*X^2*Z + X*Z^2;

void curve_pro_point_init_rand(ui d, MONTG_CURVE c, PRO_POINT p, ui n, ui_t nl, ui mu, ui_t mul, int *flag) {
    // ui_t A[nl], X[nl], Y[nl], Z[nl], mu[nl + 1];
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

    mpz_init(mp_n);
    mpz_init(mp_Y2Z);
    mpz_init(mp_iY2Z);
    mpz_init(mp_d);
    // mpz_set_ui(mp_Y2Z, 0L);
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
    // big_print_R(AX2Z, nl, "AX2Z");
    // printf("assert AX2Z eq A*X^2*Z;\n\n");

    // X^3
    big_mul(X3_, X2, nl, X, nl);
    barret_reduction(X3, X3_, 2 * nl, n, nl, mu, mul);
    // big_print_R(X3, nl, "X3");
    // printf("assert X3 eq X^3;\n\n");

    // XZ^2
    big_mul(Z2_, Z, nl, Z, nl);
    barret_reduction(Z2, Z2_, 2 * nl, n, nl, mu, mul);
    big_mul(XZ2_, X, nl, Z2, nl);
    barret_reduction(XZ2, XZ2_, 2 * nl, n, nl, mu, mul);
    // big_print_R(XZ2, nl, "XZ2");
    // printf("assert XZ2 eq X*Z^2;\n\n");
    
    // Y^2Z
    big_mul(Y2_, Y, nl, Y, nl);
    barret_reduction(Y2, Y2_, 2 * nl, n, nl, mu, mul);
    big_mul(Y2Z_, Y2, nl, Z, nl);
    barret_reduction(Y2Z, Y2Z_, 2 * nl, n, nl, mu, mul);
    // big_print_R(Y2Z, nl, "Y2Z");
    // printf("assert Y2Z eq Y^2*Z;\n\n");
    
    // R = X^3 + AX^2Z + XZ^2
    for(int i = nl + 1; i < 2 * nl; i++) {
        RHS1_[i] = 0L;
        RHS_[i] = 0L;
    }
    big_add(RHS1_, X3, nl, AX2Z, nl);
    barret_reduction(RHS1, RHS1_, 2 * nl, n, nl, mu, mul);
    big_add(RHS_, RHS1, nl, XZ2, nl);
    barret_reduction(RHS, RHS_, 2 * nl, n, nl, mu, mul);
    // big_print_R(RHS, nl, "RHS");
    // printf("assert RHS eq X^3 + A*X^2*Z + X*Z^2;\n\n");

    // B or factorize
    mpz_import(mp_Y2Z, nl, -1, 4, 0, 0, Y2Z);
    mpz_gcd(mp_d, mp_Y2Z, mp_n);
    mpz_export(d, NULL, -1, 4, 0, 0, mp_d);
    // big_print(stdout, d, nl, "GCD", NULL);

    if(mpz_cmp_ui(mp_d, 1) == 0 || mpz_cmp(mp_d, mp_n) == 0) {
        mpz_invert(mp_iY2Z, mp_Y2Z, mp_n);
        mpz_export(iY2Z, NULL, -1, 4, 0, 0, mp_iY2Z);
        // big_print_R(iY2Z, nl, "iY2Z");
        // printf("assert iY2Z eq (Y^2*Z)^-1;\n\n");
        
        big_mul(B_, RHS, nl, iY2Z, nl);
        barret_reduction(B, B_, 2 * nl, n, nl, mu, mul);
        // big_print_R(B, nl, "B");
        // printf("assert B*Y^2*Z eq X^3 + A*X^2*Z + X*Z^2;\n\n");

        c->A = A;
        c->B = B;
        c->n = n;
        p->X = X;
        p->Y = Y;
        p->Z = Z;
        *flag = 1;
    }
}

void curve_aff_point_init_rand(MONTG_CURVE c, AFF_POINT p, ui n) {
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

int main() {
    MONTG_CURVE c = (MONTG_CURVE)malloc(sizeof(MONTG_CURVE_t) * 1);
    PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    ui_t nl = 5;
    ui_t n[nl], mu[nl + 1], d[nl];
    int flag = 0;
    // FILE *fp = fopen("/home/ozbayelif/Development/FIWE/ecm/input.magma", "a");

    for (int i = 0; i < nl; i++) {
        d[i] = 0;
    }

    for (int i = 0; i < 10000; i++) {
        flag = 0;
        
        big_rand(n, nl);
        big_get_mu(mu, n, nl);

        big_print(stdout, n, nl, "n", NULL);
        fprintf(stdout, "R := Integers(n);\n\n");

        curve_pro_point_init_rand(d, c, p, n, nl, mu, nl + 1, &flag);
        if(flag == 0) {
            big_print(stdout, d, nl, "d", NULL);
            fprintf(stdout, "assert n mod d eq 0;\n\n");
        } else {
            big_print(stdout, c->A, nl, "A", "R");
            big_print(stdout, c->B, nl, "B", "R");
            big_print(stdout, p->X, nl, "X", "R");
            big_print(stdout, p->Y, nl, "Y", "R");
            big_print(stdout, p->Z, nl, "Z", "R");
            fprintf(stdout, "assert B*Y^2*Z eq X^3 + A*X^2*Z + X*Z^2;\n\n");
        }
    }
    // fclose(fp);

    return 0;
}