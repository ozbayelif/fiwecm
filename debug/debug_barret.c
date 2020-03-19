#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "mplib.h"

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