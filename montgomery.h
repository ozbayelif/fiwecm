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
void curve_pro_point_init_rand(ui d, MONTG_CURVE c, PRO_POINT p, ui n, ui_t nl, ui mu, ui_t mul, int *flag);

// Random montgomery curve and affine point initialization
void curve_aff_point_init_rand(MONTG_CURVE c, AFF_POINT p, ui n);

// Differential addition on projective coordinates
void pro_add(ui X, ui Z, ui X1, ui Z1, ui X2, ui Z2, ui Xd, ui Zd, ui n);

// Differential addition on affine coordinates
// void aff_add(ui x, ui y, ui x1, ui y1, ui x2, ui y2, ui ..)

// Doubling on projective coordinates
void pro_dbl(ui X, ui Z, ui X1, ui Z1, ui A24, ui n);

// Doubling on affine coordinates
void aff_dbl(ui x, ui z, ui x1, ui y1, ui A, ui B, ui n);

// Ladder on projective coordinates
void pro_ladder(ui X, ui Z, ui X1, ui Z1, ui k, ui n);

// Ladder on affine coordinates
void aff_ladder(ui x, ui y, ui x1, ui y1, ui k, ui n);

// Check if the point with projective coordinates is on the curve
int pro_is_on_curve(ui A, ui B, ui X, ui Y, ui Z, ui n);

// Check if the point with affine coordinates is on the curve
int aff_is_on_curve(ui A, ui B, ui x, ui y, ui n);