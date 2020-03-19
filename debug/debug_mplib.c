#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "mplib.h"

void big_rand(ui z, ui_t l) {
	int i;
	for(i=0; i<l; i++) {
		z[i] = ((ui_t)rand()) * ((ui_t)rand()) * ((ui_t)rand()) * ((ui_t)rand());
	}
}

void big_print(FILE *fp, ui a, ui_t al, char *s, char *R) {
    if(R != NULL) {
	    fprintf(fp, "%s := %s!(", s, R);
    } else {
        fprintf(fp, "%s := (", s);
    }
    fprintf(fp, "%u", a[0]);
    for(int i = 1; i < al; i++) {
        fprintf(fp, " + %u * (2^%d)^%d", a[i], W, i);
    }
    fprintf(fp, ");\n\n");
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

void big_get_mu(ui mu, ui n, ui_t nl) {
	mpz_t mp_n, mp_b2k, mp_mu;

    mpz_init(mp_n);
    mpz_init(mp_b2k);
    mpz_init(mp_mu);

	mpz_set_ui(mp_b2k, 0L);
	mpz_set_ui(mp_mu, 0L);
	
	mpz_import(mp_n, nl, -1, 4, 0, 0, n);
	mpz_add_ui(mp_b2k, mp_b2k, 1);
	mpz_mul_2exp(mp_b2k, mp_b2k, W * (2 * nl));
	mpz_fdiv_q(mp_mu, mp_b2k, mp_n);
	mpz_export(mu, NULL, -1, 4, 0, 0, mp_mu);
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

	FILE *fp = stdout;

	fprintf(fp, "b:=2^%d;\n\n", W);
    fprintf(fp, "k:=nl;\n\n");

	big_print(fp, m, ml, "m", NULL);
	big_print(fp, n, nl, "n", NULL);
	big_print(fp, mu, mul, "mu", NULL);

    big_cpy(md, m, k - 1, k + 1); // md = m / b^(k - 1) 
	big_print(fp, md, k + 1, "md", NULL);
	fprintf(fp, "assert md eq m div b^(k - 1);\n\n");

    big_mul(mdmu, md, k + 1, mu, mul); // mdmu = md * mu 
	big_print(fp, mdmu, mul + k + 1, "mdmu", NULL);
	fprintf(fp, "assert mdmu eq md * mu;\n\n");

    big_cpy(q, mdmu, k + 1, mul); // q = (m / b^(k - 1) * mu) / b^(k + 1) 
	big_print(fp, q, mul, "q", NULL);
	fprintf(fp, "assert q eq mdmu div b^(k + 1);\n\n");

    big_cpy(mm, m, 0, k + 1); // mm = m mod b^(k + 1) 
	big_print(fp, mm, k + 1, "mm", NULL);
	fprintf(fp, "assert mm eq m mod b^(k + 1);\n\n");

    big_mul(qn, q, mul, n, nl); // qn = q * n 
	big_print(fp, qn, mul + nl, "qn", NULL);
	fprintf(fp, "assert qn eq q * n;\n\n");

    big_cpy(qnm, qn, 0, k + 1); // qnm = qn mod b^(k + 1) 
	big_print(fp, qnm, k + 1, "qnm", NULL);
	fprintf(fp, "assert qnm eq qn mod b^(k + 1);\n\n");

    big_sub(z, mm, k + 1, qnm, k + 1); // r = mm - qnm
	big_print(fp, z, k + 1, "z", NULL);
	fprintf(fp, "assert z eq mm - qnm;\n\n");

    b = big_sub(r2, z, nl, n, nl); // while r >= n do: r <- r - n
	r2[nl] = z[nl] - b;
	big_print(fp, r2, nl + 1, "r2", NULL);
    while((r2[nl] >> (W - 1)) == 0) {
        big_cpy(z, r2, 0, k + 1);
        b = big_sub(r2, z, nl, n, nl);
		r2[nl] = z[nl] - b;
    }
}