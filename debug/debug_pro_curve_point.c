#include <gmp.h>
#include <stdlib.h>
#include <time.h>
#include "montgomery.h"
#include "mplib.h"

void curve_pro_point_init_rand(ui d, MONTG_CURVE c, PRO_POINT p, ui n, ui_t nl, ui mu, ui_t mul, int *flag) {
    FILE *fp = stdout;
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

    fprintf(fp, "nl:=%lu;\n\n", nl);
    big_print(fp, A, nl, "A", "R");
    big_print(fp, X, nl, "X", "R");
    big_print(fp, Y, nl, "Y", "R");
    big_print(fp, Z, nl, "Z", "R");

    // AX^2Z
    big_mul(X2_, X, nl, X, nl);
    barret_reduction(X2, X2_, 2 * nl, n, nl, mu, mul);
    big_mul(AX2_, A, nl, X2, nl);
    barret_reduction(AX2, AX2_, 2 * nl, n, nl, mu, mul);
    big_mul(AX2Z_, AX2, nl, Z, nl);
    barret_reduction(AX2Z, AX2Z_, 2 * nl, n, nl, mu, mul);
    big_print(fp, AX2Z, nl, "AX2Z", "R");
    printf("assert AX2Z eq A*X^2*Z;\n\n");

    // X^3
    big_mul(X3_, X2, nl, X, nl);
    barret_reduction(X3, X3_, 2 * nl, n, nl, mu, mul);
    big_print(fp, X3, nl, "X3", "R");
    printf("assert X3 eq X^3;\n\n");

    // XZ^2
    big_mul(Z2_, Z, nl, Z, nl);
    barret_reduction(Z2, Z2_, 2 * nl, n, nl, mu, mul);
    big_mul(XZ2_, X, nl, Z2, nl);
    barret_reduction(XZ2, XZ2_, 2 * nl, n, nl, mu, mul);
    big_print(fp, XZ2, nl, "XZ2", "R");
    printf("assert XZ2 eq X*Z^2;\n\n");
    
    // Y^2Z
    big_mul(Y2_, Y, nl, Y, nl);
    barret_reduction(Y2, Y2_, 2 * nl, n, nl, mu, mul);
    big_mul(Y2Z_, Y2, nl, Z, nl);
    barret_reduction(Y2Z, Y2Z_, 2 * nl, n, nl, mu, mul);
    big_print(fp, Y2Z, nl, "Y2Z", "R");
    printf("assert Y2Z eq Y^2*Z;\n\n");
    
    // R = X^3 + AX^2Z + XZ^2
    for(i = nl + 1; i < 2 * nl; i++) {
        RHS1_[i] = 0L;
        RHS_[i] = 0L;
    }
    big_add(RHS1_, X3, nl, AX2Z, nl);
    barret_reduction(RHS1, RHS1_, 2 * nl, n, nl, mu, mul);
    big_add(RHS_, RHS1, nl, XZ2, nl);
    barret_reduction(RHS, RHS_, 2 * nl, n, nl, mu, mul);
    big_print(fp, RHS, nl, "RHS", "R");
    printf("assert RHS eq X^3 + A*X^2*Z + X*Z^2;\n\n");

    // B or factorize
    mpz_import(mp_Y2Z, nl, -1, 4, 0, 0, Y2Z);
    mpz_gcd(mp_d, mp_Y2Z, mp_n);
    mpz_export(d, NULL, -1, 4, 0, 0, mp_d);
    big_print(fp, d, nl, "GCD", NULL);

    if(mpz_cmp_ui(mp_d, 1) == 0 || mpz_cmp(mp_d, mp_n) == 0) {
        mpz_invert(mp_iY2Z, mp_Y2Z, mp_n);
        mpz_export(iY2Z, NULL, -1, 4, 0, 0, mp_iY2Z);
        big_print(fp, iY2Z, nl, "iY2Z", "R");
        printf("assert iY2Z eq (Y^2*Z)^-1;\n\n");
        
        big_mul(B_, RHS, nl, iY2Z, nl);
        barret_reduction(B, B_, 2 * nl, n, nl, mu, mul);
        big_print(fp, B, nl, "B", "R");
        printf("assert B*Y^2*Z eq X^3 + A*X^2*Z + X*Z^2;\n\n");

        c->A = A;
        c->B = B;
        c->n = n;
        p->X = X;
        p->Y = Y;
        p->Z = Z;
        *flag = 1;
    }
}
