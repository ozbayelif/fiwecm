#include <stdlib.h>
#include <stdio.h>
#include "mplib.h"

void big_rand(ui z, ui_t l) {
	int i;
	for(i=0; i<l; i++) {
		z[i] = ((ui_t)rand()) * ((ui_t)rand()) * ((ui_t)rand()) * ((ui_t)rand());
	}
}

void big_print(ui a, ui_t al, char *s) {
    printf("%s := ", s);
    printf("%u", a[0]);
    for(int i = 1; i < al; i++) {
        printf(" + %u * (2^%d)^%d", a[i], W, i);
    }
    printf(";\n\n");
}

void big_print_F(ui a, ui_t al, char *s) {
	printf("%s := F!(", s);
    printf("%u", a[0]);
    for(int i = 1; i < al; i++) {
        printf(" + %u * (2^%d)^%d", a[i], W, i);
    }
    printf(");\n\n");
}

void big_fprint(FILE *fp, ui a, ui_t al, char *s) {
    fprintf(fp, "%s := ", s);
    fprintf(fp, "%u", a[0]);
    for(int i = 1; i < al; i++) {
        fprintf(fp, " + %u * 2^(%d * %d)", a[i], W, i);
    }
    fprintf(fp, ";\n");
}

void big_add(ui z, ui a, ui_t al, ui b, ui_t bl) {
	int j;
	ui_t carry_bit = 0;
	
	for(j = 0; j < al; j++) {
		z[j] = a[j] + b[j] + carry_bit;
		if(z[j] < a[j]) {
			carry_bit = 1;
		} else if(z[j] > a[j]) {
			carry_bit = 0;
		}
	}
	z[al] = carry_bit;
}

int big_sub(ui z, ui a, ui_t al, ui b, ui_t bl) {
	int j;
	ui_t borrow_bit = 0;

	for(j = 0; j < al; j++) {
		z[j] = a[j] - b[j] - borrow_bit;
		if(z[j] < a[j]) {
			borrow_bit = 0;
		} else if(z[j] > a[j]) {
			borrow_bit = 1;
		}
	}
    return borrow_bit;
}

void big_mul(ui z, ui a, ui_t al, ui b, ui_t bl) {
	int i, j;
	ui_t u, v;
	uni_t uv;

	for(int i = 0; i <= al; i++) {
		z[i] = 0;
	}
 	for(i = 0; i < al; i++) {
	 	u = 0;	
		for(j = 0; j < bl; j++) {
			uv = (uni_t)z[i + j] + (uni_t)a[i] * (uni_t)b[j] + (uni_t)u;
			u = uv >> 32;
			v = uv & 0xFFFFFFFF;
			z[i + j] = v; 
		}
		z[i + bl] = u;
	}		
}

uni_t barret_reduction_UL(uni_t p, uni_t b, uni_t k, uni_t z, uni_t m, uni_t L) { // Calculate z mod p where z < 2^W and p < 2^W
    uni_t bkpp = (k + 1) * L;
    uni_t bkmp = (k - 1) * L;
    uni_t bkp = 1 << bkpp;
    uni_t bkm = 1 << bkmp;
    
    uni_t q = ((z >> bkmp) * m) >> bkpp;
    uni_t r = (z & (bkp - 1)) - ((q * p) & (bkp - 1));
    if(r < 0) {
        r = r + bkp;
    }
    while(r >= p) {
        r -= p;
    }
    
    return r;
}

// ml = 2 * nl
void barret_reduction(ui z, ui m, ui_t ml, ui n, ui_t nl, ui mu, ui_t mul) { // Calculate m mod n
    ui_t k = nl, md[k + 1], mdmu[mul + k + 1], q[mul], mm[k + 1], qn[mul + nl], qnm[k + 1], r2[k + 1];
    int i, j, b;
    
    big_cpy(md, m, k - 1, k + 1); // md = m / b^(k - 1) 
    big_mul(mdmu, md, k + 1, mu, mul); // mdmu = md * mu 
    big_cpy(q, mdmu, k + 1, mul); // q = (m / b^(k - 1) * mu) / b^(k + 1) 
    big_cpy(mm, m, 0, k + 1); // mm = m mod b^(k + 1) 
    big_mul(qn, q, mul, n, nl); // qn = q * n 
    big_cpy(qnm, qn, 0, k + 1); // qnm = qn mod b^(k + 1) 
    big_sub(z, mm, k + 1, qnm, k + 1); // r = mm - qnm
    b = big_sub(r2, z, nl + 1, n, nl + 1); // while r >= n do: r <- r - n
    while(!b) {
        big_cpy(z, r2, 0, k + 1);
        b = big_sub(r2, z, nl + 1, n, nl + 1);
    }
}