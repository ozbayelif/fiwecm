/**
 * \file mplib.h
 * \brief A library to represent multi-precision
 *      integers and do arithmetic on them.
 */ 

#include <stdio.h>
#include <stdlib.h>

/**
 * \def W
 * \brief Word size of the computer
 * 
 * The numbers are represented in base \f$2^{W}\f$
 * such that a digit of the number cannot exceed
 * \f$2^{W}\f$.
 */
#define W 32

/**
 * \brief Type definition for unsigned long
 *      pointer
 */
typedef unsigned long *uni;
/**
 * \brief Type definition for unsigned long
 */
typedef unsigned long uni_t;
/**
 * \brief Type definition for unsigned integer
 *      pointer
 */
typedef unsigned int *ui;
/**
 * \brief Type definition for unsigned integer
 */
typedef unsigned int ui_t;

/**
 * \brief Copies \f$end - start\f$ elements
 *      from a to z.
 * @param[out] z destination of the copy operation
 * @param[in] a source of the copy operation
 * @param[in] start starting index for the copy operation
 * @param[in] end ending index for the copy operation
 */
#define fiwe_cpy(z, a, start, end) if(1) { \
    int i, j; \
    for(i = 0, j = (start); i < (end); i++, j++) { \
        z[i] = a[j]; \
    } \
};

/**
 * \brief Initializes z with a random multi-precision
 *      number
 * @param[out] z multi-precision number to be initialized
 * @param[in] l number of digits of z in base \f$2^W\f$
 */
void fiwe_rand(ui z, ui_t l);

/**
 * \brief Initializes z with a random multi-precision
 *      number mod n
 * @param[out] z multi-precision number to be initialized
 * @param[in] l number of digits of z in base \f$2^W\f$
 * @param[in] n modular base for z
 * @param[in] nl number of digits of n in base \f$2^W\f$
 * @param[in] mu precalculated value of  \f$(2^{W})^{2*nl} / n\f$
 * @param[in] mul number of digits of mu in base \f$2^W\f$
 */
void fiwe_mod_rand(ui z, ui_t l, ui n, ui_t nl, ui mu, ui_t mul);

/**
 * \brief Prints the given multi-precision number in
 *      Magma assignment format
 * @param[in] fp pointer to the file to print
 * @param[in] a multi-precision number to be printed
 * @param[in] al number of digits of a in base \f$2^W\f$
 * @param[in] s name of the variable going to be assigned to a
 * @param[in] R name of the ring that a is going to be defined in (optional)
 */
void fiwe_print(FILE *fp, ui a, ui_t al, char *s, char *R);

/**
 * \brief Checks if two multi-precision numbers are equal
 * @param[out] z 1 if equal, 0 otw
 * @param[in] a first number
 * @param[in] b second number
 * @param[in] l number of digits of a and b in base \f$2^W\f$
 */
void fiwe_is_equal(int *z, ui a, ui b, ui_t l);

/**
 * \brief Checks if a multi-precision number is equal to given
 *      unsigned integer
 * @param[out] z 1 if equal, 0 otw
 * @param[in] a first number
 * @param[in] al number of digits of a in base \f$2^W\f$
 * @param[in] b unsigned int to be compared
 */
void fiwe_is_equal_ui(int *z, ui a, ui_t al, ui_t b);

/**
 * \brief Adds two multi-precision numbers
 * @param[out] z result of the addition
 * @param[in] a first number
 * @param[in] al number of digits of a in base \f$2^W\f$
 * @param[in] b second number
 * @param[in] bl number of digits of b in base \f$2^W\f$
 */
void fiwe_add(ui z, ui a, ui_t al, ui b, ui_t bl);

/**
 * \brief Adds two multi-precision numbers in mod n
 * @param[out] z result of the addition
 * @param[in] a first number
 * @param[in] al number of digits of a in base \f$2^W\f$
 * @param[in] b second number
 * @param[in] bl number of digits of b in base \f$2^W\f$
 * @param[in] n modular base for the addition
 * @param[in] nl number of digits of n in base \f$2^W\f$
 * @param[in] mu precalculated value of \f$(2^{W})^{2*nl} / n\f$
 * @param[in] mul number of digits of mu in base \f$2^W\f$
 */
void fiwe_mod_add(ui z, ui a, ui_t al, ui b, ui_t bl, ui n, ui_t nl, ui mu, ui_t mul);

/**
 * \brief Subtracts two multi-precision numbers
 * @param[out] z result of the subtraction
 * @param[in] a first number
 * @param[in] al number of digits of a in base \f$2^W\f$
 * @param[in] b second number
 * @param[in] bl number of digits of b in base \f$2^W\f$
 */
void fiwe_sub(ui z, int *d, ui a, ui_t al, ui b, ui_t bl);

/**
 * \brief Subtracts two multi-precision numbers in mod n
 * @param[out] z result of the subtraction
 * @param[in] a first number
 * @param[in] al number of digits of a in base \f$2^W\f$
 * @param[in] b second number
 * @param[in] bl number of digits of b in base \f$2^W\f$
 * @param[in] n modular base for the subtraction
 * @param[in] nl number of digits of n in base \f$2^W\f$
 */
void fiwe_mod_sub(ui z, ui a, ui_t al, ui b, ui_t bl, ui n, ui_t nl);

/**
 * \brief Multiplies two multi-precision numbers
 * @param[out] z result of the multiplication
 * @param[in] a first number
 * @param[in] al number of digits of a in base \f$2^W\f$
 * @param[in] b second number
 * @param[in] bl number of digits of b in base \f$2^W\f$
 */
void fiwe_mul(ui z, ui a, ui_t al, ui b, ui_t bl);

/**
 * \brief Multiplies two multi-precision numbers in mod n
 * @param[out] z result of the multiplication
 * @param[in] a first number
 * @param[in] al number of digits of a in base \f$2^W\f$
 * @param[in] b second number
 * @param[in] bl number of digits of b in base \f$2^W\f$
 * @param[in] n modular base for the multiplication
 * @param[in] nl number of digits of n in base \f$2^W\f$
 * @param[in] mu precalculated value of \f$(2^{W})^{2*nl} / n\f$
 * @param[in] mul number of digits of mu in base \f$2^W\f$
 */
void fiwe_mod_mul(ui z, ui a, ui_t al, ui b, ui_t bl, ui n, ui_t nl, ui mu, ui_t mul);

/**
 * \brief Calculates  \f$(2^{W})^{2*nl} / n\f$
 * @param[out] z result of the calculation
 * @param[in] n n in the equation
 * @param[in] nl number of digits of n in base \f$2^W\f$
 */
void fiwe_get_mu(ui z, ui n, ui_t nl);

/**
 * \brief Calculates \f$(A + 2) / 4\f$
 * \param[out] A24 result of \f$(A + 2) / 4\f$ or a factor of n
 * @param[in] A A in the equation
 * @param[in] n n modular base for the calculation
 * @param[in] number of digits of n in base \f$2^W\f$
 * @param[in] mu precalculated value of \f$(2^{W})^{2*nl} / n\f$
 * @param[in] mul number of digits of mu in base \f$2^W\f$
 * @param[in] flag 1 when calculation succeeds, 0 when factor gets found
 */
void fiwe_get_A24(ui A24, ui A, ui n, ui_t nl, ui mu, ui_t mul, int *flag);

/**
 * \brief Calculates m mod n
 * @param[out] z result of the reduction
 * @param[in] m number to be reduced
 * @param[in] ml number of digits of m in base \f$2^W\f$
 * @param[in] n modular base for the reduction
 * @param[in] nl number of digits of n in base \f$2^W\f$
 * @param[in] mu precalculated value of \f$(2^{W})^{2*nl} / n\f$
 * @param[in] mul number of digits of mu in base \f$2^W\f$
 */
void barret_reduction(ui z, ui m, ui_t ml, ui n, ui_t nl, ui mu, ui_t mul);

/**
 * \brief Calculates GCD(a,b)
 * @param[out] d result of the GCD
 * @param[in] dl number of digits of d in base \f$2^W\f$
 * @param[in] a first operand
 * @param[in] al number of digits of a in base \f$2^W\f$
 * @param[in] b first operand
 * @param[in] bl number of digits of b in base \f$2^W\f$
 */
void fiwe_gcd(ui d, ui_t dl, ui a, ui_t al, ui b, ui_t bl);

/**
 * \brief Calculates \f$a^{-1} \pmod b\f$
 * @param[out] z result of the inversion, \f$a^{-1}\f$
 * @param[in] a first operand
 * @param[in] al number of digits of a in base \f$2^W\f$
 * @param[in] b first operand
 * @param[in] bl number of digits of b in base \f$2^W\f$
 */
int fiwe_invert(ui z, ui a, ui_t al, ui b, ui_t bl);