#include <gmp.h>
#include <stdlib.h>
#include <time.h>
#include "montgomery.h"
#include "mplib.h"

void pro_curve_point(ui d, MONTG_CURVE c, PRO_POINT p, ui n, ui_t nl, ui mu, ui_t mul, int *flag) {
    ui_t A__[2 * nl], X__[2 * nl], Y__[2 * nl], Z__[2 * nl], A2[nl], X2[nl], AX2[nl], AX2Z[nl], X3[nl], Z2[nl], XZ2[nl], Y2[nl], Y2Z[nl], iY2Z[nl], RHS1[nl], RHS[nl];
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
    big_mod_rand(A, nl, n, nl, mu, mul);
    big_mod_rand(X, nl, n, nl, mu, mul);
    big_mod_rand(Y, nl, n, nl, mu, mul);
    big_mod_rand(Z, nl, n, nl, mu, mul);
    
    big_mod_mul(A2, A, nl, A, nl, n, nl, mu, mul);
    if(A2[0] == 4L) {
        s = 0;
    }
    for(i = 1; i < nl; i++) {
        if(A2[i] != 0L) {
            s = 1;
        }
    }                                                           // Check singularity
    if(s != 0) {
        big_mod_mul(X2, X, nl, X, nl, n, nl, mu, mul);          // X2 = X^2
        big_mod_mul(AX2, A, nl, X2, nl, n, nl, mu, mul);        // AX = AX^2
        big_mod_mul(AX2Z, AX2, nl, Z, nl, n, nl, mu, mul);      // AX2Z = AX^2Z
        big_mod_mul(X3, X2, nl, X, nl, n, nl, mu, mul);         // X3 = X^3
        big_mod_mul(Z2, Z, nl, Z, nl, n, nl, mu, mul);          // Z2 = Z^2
        big_mod_mul(XZ2, X, nl, Z2, nl, n, nl, mu, mul);        // XZ2 = XZ^2
        big_mod_mul(Y2, Y, nl, Y, nl, n, nl, mu, mul);          // Y2 = Y^2
        big_mod_mul(Y2Z, Y2, nl, Z, nl, n, nl, mu, mul);        // Y2Z = Y^2Z
        big_mod_add(RHS1, X3, nl, AX2Z, nl, n, nl);                 // RHS1 = X^3 + AX^2Z
        big_mod_add(RHS, RHS1, nl, XZ2, nl, n, nl);                 // RHS = X^3 + AX^2Z + XZ^2
        mpz_import(mp_Y2Z, nl, -1, 4, 0, 0, Y2Z);
        mpz_gcd(mp_d, mp_Y2Z, mp_n);
        mpz_export(d, NULL, -1, 4, 0, 0, mp_d);                 // d = GCD(Y^2Z, n)

        if(mpz_cmp_ui(mp_d, 1) == 0 || mpz_cmp(mp_d, mp_n) == 0) {
            for(i = 0; i < nl; i++) {
                iY2Z[i] = 0L;
            }
            mpz_invert(mp_iY2Z, mp_Y2Z, mp_n);
            mpz_export(iY2Z, NULL, -1, 4, 0, 0, mp_iY2Z);       // iY2Z = Inv(Y^2Z)
            big_mod_mul(B, RHS, nl, iY2Z, nl, n, nl, mu, mul);     // B = RHS * Inv(Y^2Z)
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
    ui_t A2[nl], x2[nl], x3[nl], Ax2[nl], y2[nl], iy2[nl], rhs1[nl], rhs[nl];
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
    big_mod_rand(A, nl, n, nl, mu, mul);
    big_mod_rand(x, nl, n, nl, mu, mul);
    big_mod_rand(y, nl, n, nl, mu, mul);

    big_mod_mul(A2, A, nl, A, nl, n, nl, mu, mul);
    if(A2[0] == 4L) {
        s = 0;
    }
    for(i = 1; i < nl; i++) {
        if(A2[i] != 0L) {
            s = 1;
        }
    }                                                           // Check singularity
    if(s != 0) {
        big_mod_mul(x2, x, nl, x, nl, n, nl, mu, mul);          // x2 = x^2
        big_mod_mul(Ax2, A, nl, x2, nl, n, nl, mu, mul);        // Ax2 = Ax^2
        big_mod_mul(x3, x2, nl, x, nl, n, nl, mu, mul);         // x3 = x^3
        big_mod_mul(y2, y, nl, y, nl, n, nl, mu, mul);          // y2 = y^2
        big_mod_add(rhs1, x3, nl, Ax2, nl, n, nl);                  // rhs1 = x^3 + Ax^2
        big_mod_add(rhs, rhs1, nl, x, nl, n, nl);                   // rhs = x^3 + Ax^2 + x
        mpz_import(mp_y2, nl, -1, 4, 0, 0, y2);
        mpz_gcd(mp_d, mp_y2, mp_n);
        mpz_export(d, NULL, -1, 4, 0, 0, mp_d);                 // d = GCD(y^2, n)

        if(mpz_cmp_ui(mp_d, 1) == 0 || mpz_cmp(mp_d, mp_n) == 0) {
            for(i = 0; i < nl; i++) {
                iy2[i] = 0L;
            }
            mpz_invert(mp_iy2, mp_y2, mp_n);
            mpz_export(iy2, NULL, -1, 4, 0, 0, mp_iy2);         // iy2 = Inv(y^2)
            big_mod_mul(B, rhs, nl, iy2, nl, n, nl, mu, mul);       // B = rhs * Inv(y^2)
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
    ui_t a[nl], b[nl], bs[nl + 1], c[nl], d[nl], ds[nl + 1], da[nl], cb[nl], e[nl], f[nl], fs[nl + 1], e2[nl], f2[nl];
    ui X = (ui)malloc(sizeof(ui_t) * nl);
    ui Z = (ui)malloc(sizeof(ui_t) * nl);
    int i, s;

    big_mod_add(a, p2->X, nl, p2->Z, nl, n, nl);                 // a = X2 + Z2
    big_mod_sub(b, p2->X, nl, p2->Z, nl, n, nl);
    big_mod_add(c, p1->X, nl, p1->Z, nl, n, nl);                // c = X1 + Z1
    big_mod_sub(d, p1->X, nl, p1->Z, nl, n, nl);
    big_mod_mul(da, d, nl, a, nl, n, nl, mu, mul);          // da = d * a
    big_mod_mul(cb, c, nl, b, nl, n, nl, mu, mul);          // cb = c * b
    big_mod_add(e, da, nl, cb, nl, n, nl);                      // e = da + cb
    big_mod_sub(f, da, nl, cb, nl, n, nl);
    big_mod_mul(e2, e, nl, e, nl, n, nl, mu, mul);          // e2 = e^2
    big_mod_mul(f2, f, nl, f, nl, n, nl, mu, mul);          // f2 = f^2
    big_mod_mul(X, pd->Z, nl, e2, nl, n, nl, mu, mul);      // X = Zd * e2
    big_mod_mul(Z, pd->X, nl, f2, nl, n, nl, mu, mul);      // Z = Xd * f2

    p->X = X;
    p->Z = Z;
}

// Reuqires division
void aff_add(AFF_POINT p, AFF_POINT p1, AFF_POINT p2, ui A, ui B, ui n, ui_t nl, ui mu, ui_t mul) {
    // ui_t c[nl], cs[nl + 1], c2_[2 * nl], c2[nl], c3_[2 * nl], c3[nl], d[nl], ds[nl + 1], d2_[2 * nl], d2[nl], d3_[2 * nl], d3[nl];
    // ui x = (ui)malloc(sizeof(ui_t) * nl);
    // ui y = (ui)malloc(sizeof(ui_t) * nl);
    // int i, s;

    // for(i = nl + 1; i < 2 * nl; i++) {
    //     c2_[i] = 0L;
    //     d2_[i] = 0L;
    //     c3_[i] = 0L;
    // }

    // s = big_sub(c, p2->y, nl, p1->y, nl);
    // if(s) {
    //     do {
    //         big_mod_add(cs, c, nl, n, nl, n, nl);
    //         big_cpy(c, cs, 0, nl);
    //     } while((cs[nl] == 0));
    // }                                                       // c = y2 - y1
    // s = big_sub(d, p2->x, nl, p1->x, nl);
    // if(s) {
    //     do {
    //         big_mod_add(ds, d, nl, n, nl, n, nl);
    //         big_cpy(d, ds, 0, nl);
    //     } while((ds[nl] == 0));
    // }                                                       // d = x2 - x1
    // big_mul(c2_, c, nl, c, nl);
    // barret_reduction(c2, c2_, 2 * nl, n, nl, mu, mul);      // c2 = c^2
    // big_mul(d2_, d, nl, d, nl);
    // barret_reduction(d2, d2_, 2 * nl, n, nl, mu, mul);      // d2 = d^2
    // big_mul(c3_, c2, nl, c, nl);
    // barret_reduction(c3, c3_, 2 * nl, n, nl, mu, mul);      // c3 = c^3
    // big_mul(d3_, d2, nl, d, nl);
    // barret_reduction(d3, d3_, 2 * nl, n, nl, mu, mul);      // d3 = d^3
}

void pro_dbl(PRO_POINT p, PRO_POINT p1, ui A24, ui n, ui_t nl, ui mu, ui_t mul) {
    ui_t a [nl], a2[nl], bs[nl + 1], b[nl], b2[nl], cs[nl + 1], c[nl], d[nl], e[nl];
    ui X = (ui)malloc(sizeof(ui_t) * nl);
    ui Z = (ui)malloc(sizeof(ui_t) * nl);
    int i, s;

    big_mod_add(a, p1->X, nl, p1->Z, nl, n, nl);                    // a = X + Z
    big_mod_mul(a2, a, nl, a, nl, n, nl, mu, mul);              // a2 = a^2
    s = big_sub(b, p1->X, nl, p1->Z, nl);
    if(s) {
        do {
            big_mod_add(bs, b, nl, n, nl, n, nl);
            big_cpy(b, bs, 0, nl);
        } while((bs[nl] == 0));
    }                                                           // b = X - Z
    big_mod_mul(b2, b, nl, b, nl, n, nl, mu, mul);              // b2 = b^2
    s = big_sub(c, a2, nl, b2, nl);
    if(s) {
        do {
            big_mod_add(cs, c, nl, n, nl, n, nl);
            big_cpy(c, cs, 0, nl);
        } while((cs[nl] == 0));
    }                                                           // c = a2 - b2
    big_mod_mul(X, a2, nl, b2, nl, n, nl, mu, mul);             // X = a2 * b2
    big_mod_mul(d, A24, nl, c, nl, n, nl, mu, mul);             // d = a24 * c
    big_mod_add(e, b2 ,nl, d, nl, n, nl);                           // e = b2 + d
    big_mod_mul(Z, c, nl, e, nl, n, nl, mu, mul);               // Z = c * e

    p->X = X;
    p->Z = Z;
}

// // // Reuqires division
// // void aff_dbl(ui x, ui z, ui x1, ui y1, ui A, ui B, ui n) {
// // }

// // void pro_ladder(PRO_POINT p, PRO_POINT p1, ui A24, ui k, ui_t kl, ui n, ui_t nl, ui mu, ui_t mul) {
// //     ui_t a, x;
// //     PRO_POINT R0 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
// //     PRO_POINT R1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
// //     PRO_POINT_t R0_, R1_;
// //     int i, j;

// //     R0->X = p1->X;
// //     R0->Z = p1->Z;
// //     pro_dbl(R1, p1, A24, n, nl, mu, mul);

// //     a = k[kl - 1];
// //     j = 0;
// //     while(a > 0) {
// //         a = a >> 1;
// //         j++;
// //     }                                                           // Find the index of the first 1
// //     for(i = j - 2; i >= 0; i--) {
// //         x = 1;
// //         x <<= i;
// //         R0_->X = R0->X;
// //         R0_->Z = R0->Z;
// //         R1_->X = R1->X;
// //         R1_->Z = R1->Z;
// //         if(!(k[kl - 1] & x)) {
// //             pro_dbl(R0, R0_, A24, n, nl, mu, mul);
// //             pro_add(R1, R0_, R1_, p1, n, nl, mu, mul);
// //         } else {
// //             pro_add(R0, R0_, R1_, p1, n, nl, mu, mul);
// //             pro_dbl(R1, R1_, A24, n, nl, mu, mul);
// //         }
// //     }
// //     for (i = kl - 2; i >= 0; i--) {
// //         for(j = W - 1; j >= 0; j--) {
// //             x = 1;
// //             x <<= i;
// //             R0_->X = R0->X;
// //             R0_->Z = R0->Z;
// //             R1_->X = R1->X;
// //             R1_->Z = R1->Z;
// //             if(!(k[kl - 1] & x)) {
// //                 pro_dbl(R0, R0_, A24, n, nl, mu, mul);
// //                 pro_add(R1, R0_, R1_, p1, n, nl, mu, mul);
// //             } else {
// //                 pro_add(R0, R0_, R1_, p1, n, nl, mu, mul);
// //                 pro_dbl(R1, R1_, A24, n, nl, mu, mul);
// //             }
// //         }
// //     }
// //     p->X = R0->X;
// //     p->Z = R0->Z;
// // }

// // void aff_ladder(ui x, ui y, ui x1, ui y1, ui k, ui n) {
// // }

// // int pro_is_on_curve(ui A, ui B, ui X, ui Y, ui Z, ui n) {
// // }

// // int aff_is_on_curve(ui A, ui B, ui x, ui y, ui n) {
// // }