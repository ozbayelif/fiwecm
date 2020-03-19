#define W 32
#include <stdio.h>
#include <stdlib.h>

typedef unsigned long *uni;
typedef unsigned long uni_t;
typedef unsigned int *ui;
typedef unsigned int ui_t;

#define big_cpy(z, a, start, end) if(1) { \
    int i, j; \
    for(i = 0, j = (start); i < (end); i++, j++) { \
        z[i] = a[j]; \
    } \
};
void big_rand(ui z, ui_t l);
void big_print(FILE *fp, ui a, ui_t al, char *s, char *R);
void big_add(ui z, ui a, ui_t al, ui b, ui_t bl);
int big_sub(ui z, ui a, ui_t al, ui b, ui_t bl);
void big_mul(ui z, ui a, ui_t al, ui b, ui_t bl);
void big_get_mu(ui z, ui n, ui_t nl);
uni_t barret_reduction_UL(uni_t p, uni_t b, uni_t k, uni_t z, uni_t m, uni_t L);
void barret_reduction(ui z, ui m, ui_t ml, ui n, ui_t nl, ui mu, ui_t mul);