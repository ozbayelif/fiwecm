/**
 * \file main.c
 * \brief Main of the software
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "test.h"

int main() {
    time_t start, end;
    start = time(NULL);
    // pro_curve_point_gmp_test(10000);
    // aff_curve_point_gmp_test(10000);
    // pro_add_gmp_test(10000);
    // pro_add_magma_test(10000);
    // pro_dbl_magma_test(10000);
    // pro_ladder_gmp_test(10000);
    pro_ladder_magma_test(100);
    // ecm_test(10000);
    end = time(NULL);
    printf("Time spent: %ld\n", end - start);
}