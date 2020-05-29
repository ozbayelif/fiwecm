/**
 * \file ecm.h
 * \brief A library to implement Elliptic Curve Method
 *      to factorize composite numbers.
 */ 

#include "mplib.h"

/**
 * \brief Number of different k values to try
 *      during the factorization
 * 
 * If a k value does not produce a factor,
 * another one is tried until the number of
 * trials exceed k_THRESHOLD.
 */
#define k_THRESHOLD 1000

/**
 * \brief Number of different curves to try
 *      during the factorization
 * 
 * If a curve does not produce a factor,
 * another one is tried until the number of
 * trials exceed CRV_THRESHOLD.
 */
#define CRV_THRESHOLD 25

/**
 * \brief Total number of primes that are
 *      put on prime list
 */
#define PRIMESL 21

/**
 * \brief Total number of power values that
 *      are put on power list
 */
#define POWERSL 5

/**
 * \brief A list of small primes
 * 
 * Used by the generate_B_smooth
 * function
 */
extern int primes[PRIMESL];

/**
 * \brief A list of small powers for
 *      primes
 * 
 * Used by the generate_B_smooth
 * function
 */
extern int prime_powers[POWERSL];

/**
 * @brief Generates a B-smooth number
 * @param[out] z the number generated
 * @param[in] l number of digits of z in base \f$2^W\f$
 */
void generate_B_smooth(ui z, ui_t l);

/**
 * \brief Factorizes the given composite using Elliptic Curve Method
 * @param[out] d factor of n
 * @param[in] n number to be factorized
 * @param[in] nl number of digits of n in base \f$2^W\f$
 * @return an integer, positive when ECM succeeds, 0 otw
 */
int ecm(ui d, ui n, ui_t nl);