/**
 * \file montgomery.h
 * \brief A library to represent Montgomery curves and
 *      points on them, and to do arithmetic with that
 *      points.
 */ 

#include "mplib.h"

/**
 * \brief Structure to represent a Montgomery curve
 * 
 * Curve Equations:
 * - Affine Coordinates: B*y^2 = x^3 + A*x^2 + x
 * - Projective Coordinates: B*y^2*Z = x^3 + A*x^2*Z + x*Z^2
 */
typedef struct Montg_Curve_s {
    ui A; /**< Coefficient A in the equation*/
    ui B; /**< Coefficient B in the equation*/
    ui n; /**< Modular base of the curve*/
} MONTG_CURVE_t[1], *MONTG_CURVE;

/**
 * \brief Structure to represent a projective point
 */
typedef struct Pro_Point_s {
    ui X; /**< X coordinate of the point*/
    ui Y; /**< Y coordinate of the point*/
    ui Z; /**< Z coordinate of the point*/
} PRO_POINT_t[1], *PRO_POINT;

/**
 * \brief Initializes a randomly generated Montgomery curve
 *      and a projective point on the curve
 * @param[out] d factor of n, that may be found while generating the curve
 * @param[out] c curve to be initialized
 * @param[out] p projective point to be initialized
 * @param[in] n modular base of the curve that is going to be generated
 * @param[in] nl number of digits of n in base \f$2^W\f$
 * @param[in] mu precalculated value of \f$(2^{W})^{2*nl} / n\f$
 * @param[in] mul number of digits of mu in base \f$2^W\f$
 * @param[in] flag 0 when factor found, 1 when curve and point generated, -1 when function failed due to singular curve generation
 */
void pro_curve_point(ui d, MONTG_CURVE c, PRO_POINT p, ui n, ui_t nl, ui mu, ui_t mul, int *flag);

/**
 * \brief Adds two projective points on a curve
 * @param[out] p resulted projective point
 * @param[in] p1 first operand of the addition
 * @param[in] p2 second operand of the addition
 * @param[in] pd differences of the first and second operands
 * @param[in] A24 \f$(A+2)/4\f$ where A is the coefficient of the curve
 * @param[in] n modular base of the curve that is going to be generated
 * @param[in] nl number of digits of n in base \f$2^W\f$
 * @param[in] mu precalculated value of \f$(2^{W})^{2*nl} / n\f$
 * @param[in] mul number of digits of mu in base \f$2^W\f$
 */
void pro_add(PRO_POINT p, PRO_POINT p1, PRO_POINT p2, PRO_POINT pd, ui A24, ui n, ui_t nl, ui mu, ui_t mul);

/**
 * \brief Doubles a projective point on a curve
 * @param[out] p resulted projective point
 * @param[in] p1 point to be doubled
 * @param[in] A24 \f$(A+2)/4\f$ where A is the coefficient of the curve
 * @param[in] n modular base of the curve that is going to be generated
 * @param[in] nl number of digits of n in base \f$2^W\f$
 * @param[in] mu precalculated value of \f$(2^{W})^{2*nl} / n\f$
 * @param[in] mul number of digits of mu in base \f$2^W\f$
 */
void pro_dbl(PRO_POINT p, PRO_POINT p1, ui A24, ui n, ui_t nl, ui mu, ui_t mul);

/**
 * \brief Multiplies a projective point with a constant
 * @param[out] p resulted projective point
 * @param[in] p1 point to be multiplied
 * @param[in] A24 \f$(A+2)/4\f$ where A is the coefficient of the curve
 * @param[in] k constant to multiply with p1
 * @param[in] kl number of digits of k in base \f$2^W\f$
 * @param[in] n modular base of the curve that is going to be generated
 * @param[in] nl number of digits of n in base \f$2^W\f$
 * @param[in] mu precalculated value of \f$(2^{W})^{2*nl} / n\f$
 * @param[in] mul number of digits of mu in base \f$2^W\f$
 */
void pro_ladder(PRO_POINT p, PRO_POINT p1, ui A24, ui k, ui_t kl, ui n, ui_t nl, ui mu, ui_t mul);