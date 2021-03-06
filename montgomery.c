/**
 * \file montgomery.c
 * \brief Implementation of montgomery.h library.
 */ 

#include <gmp.h>
#include <stdlib.h>
#include <time.h>
#include "montgomery.h"
#include "mplib.h"

void pro_curve_point(ui d, MONTG_CURVE c, PRO_POINT p, ui n, ui_t nl, ui mu, ui_t mul, int *flag) {
    ui_t A2[nl], X2[nl], AX2[nl], AX2Z[nl], X3[nl], Z2[nl], XZ2[nl], Y2[nl], Y2Z[nl], iY2Z[nl], RHS1[nl], RHS[nl];
    ui A = (ui)malloc(sizeof(ui_t) * nl);
    ui B = (ui)malloc(sizeof(ui_t) * nl);
    ui X = (ui)malloc(sizeof(ui_t) * nl);
    ui Y = (ui)malloc(sizeof(ui_t) * nl);
    ui Z = (ui)malloc(sizeof(ui_t) * nl);
    int is_four = 1, is_zero, is_one, is_n;

    fiwe_mod_rand(A, nl, n, nl, mu, mul);
    fiwe_mod_rand(X, nl, n, nl, mu, mul);
    fiwe_mod_rand(Y, nl, n, nl, mu, mul);
    fiwe_mod_rand(Z, nl, n, nl, mu, mul);
    
    fiwe_mod_mul(A2, A, nl, A, nl, n, nl, mu, mul);
    fiwe_is_equal_ui(&is_four, A2, nl, 4L);
    if(is_four != 1) {
        fiwe_mod_mul(X2, X, nl, X, nl, n, nl, mu, mul);          // X2 = X^2
        fiwe_mod_mul(AX2, A, nl, X2, nl, n, nl, mu, mul);        // AX = AX^2
        fiwe_mod_mul(AX2Z, AX2, nl, Z, nl, n, nl, mu, mul);      // AX2Z = AX^2Z
        fiwe_mod_mul(X3, X2, nl, X, nl, n, nl, mu, mul);         // X3 = X^3
        fiwe_mod_mul(Z2, Z, nl, Z, nl, n, nl, mu, mul);          // Z2 = Z^2
        fiwe_mod_mul(XZ2, X, nl, Z2, nl, n, nl, mu, mul);        // XZ2 = XZ^2
        fiwe_mod_mul(Y2, Y, nl, Y, nl, n, nl, mu, mul);          // Y2 = Y^2
        fiwe_mod_mul(Y2Z, Y2, nl, Z, nl, n, nl, mu, mul);        // Y2Z = Y^2Z
        fiwe_mod_add(RHS1, X3, nl, AX2Z, nl, n, nl, mu, mul);    // RHS1 = X^3 + AX^2Z
        fiwe_mod_add(RHS, RHS1, nl, XZ2, nl, n, nl, mu, mul);    // RHS = X^3 + AX^2Z + XZ^2
        fiwe_gcd(d, nl, Y2Z, nl, n, nl);                         // d = GCD(Y^2Z, n)
        fiwe_is_equal_ui(&is_one, d, nl, 1L);
        fiwe_is_equal(&is_n, d, n, nl);
        if(is_one || is_n) {
            fiwe_invert(iY2Z, Y2Z, nl, n, nl);                   // iY2Z = Inv(Y2Z)
            fiwe_mod_mul(B, RHS, nl, iY2Z, nl, n, nl, mu, mul);  // B = RHS * Inv(Y^2Z)
            fiwe_is_equal_ui(&is_zero, B, nl, 0L);               // If B mod n = 0, singular
            if(!is_zero) {
                c->A = A;
                c->B = B;
                c->n = n;
                p->X = X;
                p->Y = Y;
                p->Z = Z;
                *flag = 1;
            } else {
                *flag = -1;
            }
        } else {
            *flag = 0;
        }
    } else {
        *flag = -1;
    }
}

void pro_add(PRO_POINT p, PRO_POINT p1, PRO_POINT p2, PRO_POINT pd, ui A24, ui n, ui_t nl, ui mu, ui_t mul) {
    ui_t a[nl], b[nl],  c[nl], d[nl], da[nl], cb[nl], e[nl], f[nl], e2[nl], f2[nl];
    ui X = (ui)malloc(sizeof(ui_t) * nl);
    ui Z = (ui)malloc(sizeof(ui_t) * nl);
    int is_zero1, is_zero2, is_dbl1, is_dbl2;

    fiwe_is_equal_ui(&is_zero1, p1->Z, nl, 0L);
    fiwe_is_equal_ui(&is_zero2, p2->Z, nl, 0L);
    fiwe_is_equal(&is_dbl1, p1->X, p2->X, nl);
    fiwe_is_equal(&is_dbl2, p1->Z, p2->Z, nl);

    if(is_zero1 == 1) {
        p->X = p2->X;
        p->Z = p2->Z;
    } else if(is_zero2 == 1) {
        p->X = p1->X;
        p->Z = p1->Z;
    } else if(is_dbl1 == 1 && is_dbl2 == 1) {
        pro_dbl(p, p1, A24, n, nl, mu, mul);
    } else {
        fiwe_mod_add(a, p2->X, nl, p2->Z, nl, n, nl, mu, mul);   // a = X2 + Z2
        fiwe_mod_sub(b, p2->X, nl, p2->Z, nl, n, nl);            // b = X2 - Z2
        fiwe_mod_add(c, p1->X, nl, p1->Z, nl, n, nl, mu, mul);   // c = X1 + Z1
        fiwe_mod_sub(d, p1->X, nl, p1->Z, nl, n, nl);            // d = X1 - Z1
        fiwe_mod_mul(da, d, nl, a, nl, n, nl, mu, mul);          // da = d * a
        fiwe_mod_mul(cb, c, nl, b, nl, n, nl, mu, mul);          // cb = c * b
        fiwe_mod_add(e, da, nl, cb, nl, n, nl, mu, mul);         // e = da + cb
        fiwe_mod_sub(f, da, nl, cb, nl, n, nl);                  // f = da - cb
        fiwe_mod_mul(e2, e, nl, e, nl, n, nl, mu, mul);          // e2 = e^2
        fiwe_mod_mul(f2, f, nl, f, nl, n, nl, mu, mul);          // f2 = f^2
        fiwe_mod_mul(X, pd->Z, nl, e2, nl, n, nl, mu, mul);      // X = Zd * e2
        fiwe_mod_mul(Z, pd->X, nl, f2, nl, n, nl, mu, mul);      // Z = Xd * f2
        p->X = X;
        p->Z = Z;
    }

}

void pro_dbl(PRO_POINT p, PRO_POINT p1, ui A24, ui n, ui_t nl, ui mu, ui_t mul) {
    ui_t a [nl], a2[nl], b[nl], b2[nl], c[nl], d[nl], e[nl];
    ui X = (ui)malloc(sizeof(ui_t) * nl);
    ui Z = (ui)malloc(sizeof(ui_t) * nl);
    int is_zero1;

    fiwe_is_equal_ui(&is_zero1, p1->Z, nl, 0L);

    if(is_zero1 == 1) {
        p->X = p1->X;
        p->Z = p1->Z;
    } else {
        fiwe_mod_add(a, p1->X, nl, p1->Z, nl, n, nl, mu, mul);   // a = X + Z
        fiwe_mod_mul(a2, a, nl, a, nl, n, nl, mu, mul);          // a2 = a^2
        fiwe_mod_sub(b, p1->X, nl, p1->Z, nl, n, nl);            // b = X - Z
        fiwe_mod_mul(b2, b, nl, b, nl, n, nl, mu, mul);          // b2 = b^2
        fiwe_mod_sub(c, a2, nl, b2, nl, n, nl);                  // c = a2 - b2
        fiwe_mod_mul(X, a2, nl, b2, nl, n, nl, mu, mul);         // X = a2 * b2
        fiwe_mod_mul(d, A24, nl, c, nl, n, nl, mu, mul);         // d = a24 * c
        fiwe_mod_add(e, b2 ,nl, d, nl, n, nl, mu, mul);          // e = b2 + d
        fiwe_mod_mul(Z, c, nl, e, nl, n, nl, mu, mul);           // Z = c * e
        p->X = X;
        p->Z = Z;
    }
}

void pro_ladder(PRO_POINT p, PRO_POINT p1, ui A24, ui k, ui_t kl, ui n, ui_t nl, ui mu, ui_t mul) {
    ui_t a, x;
    PRO_POINT_t R0, R1, R0_, R1_;
    ui_t R0X[nl], R0Y[nl], R0Z[nl], R1X[nl], R1Y[nl], R1Z[nl], R0_X[nl], R0_Y[nl], R0_Z[nl], R1_X[nl], R1_Y[nl], R1_Z[nl];
    R0->X= R0X; R0->Y = R0Y; R0->Z = R0Z; R0_->X = R0_X; R0_->Y = R0_Y; R0_->Z = R0_Z; R1->X= R1X; R1->Y = R1Y; R1->Z = R1Z; R1_->X = R1_X; R1_->Y = R1_Y; R1_->Z = R1_Z;
    p->X = (ui)malloc(sizeof(ui_t) * nl);
    p->Z = (ui)malloc(sizeof(ui_t) * nl);
    int m, l, is_zero;

    fiwe_cpy(R0->X, p1->X, 0, nl);
    fiwe_cpy(R0->Z, p1->Z, 0, nl);
    pro_dbl(R1, p1, A24, n, nl, mu, mul);

    a = k[kl - 1];
    m = 0;
    while(a > 0) {
        a = a >> 1;
        m++;
    }                                                       // Find the index of the first 1
    for(l = m - 2; l >= 0; l--) {
        x = 1;
        x <<= l;
        fiwe_cpy(R0_->X, R0->X, 0, nl);
        fiwe_cpy(R0_->Z, R0->Z, 0, nl);
        fiwe_cpy(R1_->X, R1->X, 0, nl);
        fiwe_cpy(R1_->Z, R1->Z, 0, nl);
        if(!(k[0] & x)) {
            pro_dbl(R0, R0_, A24, n, nl, mu, mul);
            pro_add(R1, R0_, R1_, p1, A24, n, nl, mu, mul);
        } else {
            pro_add(R0, R0_, R1_, p1, A24, n, nl, mu, mul);
            pro_dbl(R1, R1_, A24, n, nl, mu, mul);
        }
        fiwe_is_equal_ui(&is_zero, R0->Z, nl, 0L);
        if(is_zero == 1) {
            fiwe_cpy(p->X, R0->X, 0, nl);
            fiwe_cpy(p->Z, R0->Z, 0, nl);
            k[0] >>= l;
            return;
        }
    }
    fiwe_cpy(p->X, R0->X, 0, nl);
    fiwe_cpy(p->Z, R0->Z, 0, nl);
}