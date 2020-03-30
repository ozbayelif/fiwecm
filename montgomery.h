#include "mplib.h"

typedef struct Montg_Curve_s {
    ui A;
    ui B; 
    ui n;
} MONTG_CURVE_t[1], *MONTG_CURVE;

typedef struct Pro_Point_s {
    ui X;
    ui Y; // May be ignored
    ui Z;
} PRO_POINT_t[1], *PRO_POINT;

typedef struct Aff_Point_s {
    ui x;
    ui y;
} AFF_POINT_t[1], *AFF_POINT;

// Random montgomery curve and projective point initialization
void pro_curve_point(ui d, MONTG_CURVE c, PRO_POINT p, ui n, ui_t nl, ui mu, ui_t mul, int *flag);

// Random montgomery curve and affine point initialization
void aff_curve_point(ui d, MONTG_CURVE c, AFF_POINT p, ui n, ui_t nl, ui mu, ui_t mul, int *flag);

// Differential addition on projective coordinates
void pro_add(PRO_POINT p, PRO_POINT p1, PRO_POINT p2, PRO_POINT pd, ui n, ui_t nl, ui mu, ui_t mul);

// Addition on affine coordinates
void aff_add(AFF_POINT p, AFF_POINT p1, AFF_POINT p2, ui A, ui B, ui n, ui_t nl, ui mu, ui_t mul);

// Doubling on projective coordinates
void pro_dbl(PRO_POINT p, PRO_POINT p1, ui A24, ui n, ui_t nl, ui mu, ui_t mul);

// Doubling on affine coordinates
void aff_dbl(ui x, ui z, ui x1, ui y1, ui A, ui B, ui n);

// Ladder on projective coordinates
void pro_ladder(PRO_POINT p, PRO_POINT p1, ui A24, ui k, ui_t kl, ui n, ui_t nl, ui mu, ui_t mul);

// Ladder on affine coordinates
void aff_ladder(ui x, ui y, ui x1, ui y1, ui k, ui n);

// Check if the point with projective coordinates is on the curve
int pro_is_on_curve(ui A, ui B, ui X, ui Y, ui Z, ui n);

// Check if the point with affine coordinates is on the curve
int aff_is_on_curve(ui A, ui B, ui x, ui y, ui n);