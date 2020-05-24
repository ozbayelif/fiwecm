/**
 * \file main.c
 * \brief Main of the software
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gmp.h>
#include <string.h>
#include "mplib.h"
#include "montgomery.h"
#include "test.h"
#include "ecm.h"

int SIZE;

int get_size(char *nstr, int len) {
    int i, j = 0, max = 0;

    for(i = 0; i < len; i++) {
        if(nstr[i] != '\n') {
            j++;
        } else {
            if(j > max) {
                max = j;
            }
            j = 0;
        }
    }
    if(j > max) {
        max = j;
    }

    return max;
}

ui str2arr(char *nstr, int len, int nl_dec, int nl) {
    int i, j, k;
    mpz_t mp_n;
    int firstIndexN = 0, lastIndexN = nl;

    ui_t n_[nl];
    char *num = (char *)malloc(sizeof(char) * nl_dec);
    ui n = (ui)malloc(sizeof(ui_t) * SIZE * nl);

    mpz_init(mp_n);

    for(j = 0; j < nl_dec; j++) {
        num[j] = '0';
    }
    j = 0;
    for(i = 0; i < len; i++) {
        if(nstr[i] != '\n') {
            num[j] = nstr[i];
            j++;
        } else {
            num[j] = '\0';
            num = realloc(num, sizeof(char) * j);

            mpz_set_str(mp_n, num, 10);
            for(k = 0; k < nl; k++) {
                n_[k] = 0L;
            }
            mpz_export(n_, NULL, -1, 4, 0, 0, mp_n);

            j = 0;
            for(k = firstIndexN; k < lastIndexN; k++) {
                n[k] = n_[j];
                j++;
            }
            firstIndexN += nl;
            lastIndexN += nl;

            num = realloc(num, sizeof(char) * nl_dec);
            for(j = 0; j < nl_dec; j++) {
                num[j] = '0';
            }
            j = 0;
        }
    }
    return n;
}

void str2txt(char *dstr, int nl, char *path) {
    FILE *fp = fopen(path, "w");
    fwrite(dstr, nl, SIZE, fp);
    fclose(fp);
}

char *txt2str(char *path) {
    FILE *fp = fopen(path, "r");

    fseek(fp, 0, SEEK_END);
    long len = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char *nstr = (char *)malloc((len + 1) * sizeof(char));
    fread(nstr, 1, len, fp);

    fclose(fp);

    return nstr;
}

char *factorize(char *nstr1, int len, FILE *fp) {
    int i, j, dIndex = 0;
    mpz_t mp_d;

    char *res = (char *)malloc(sizeof(char) * len * SIZE + SIZE);
    mpz_init(mp_d);

    char *nstr = (char *)malloc(sizeof(char) * (len + 1));
    strcpy(nstr, nstr1);
    nstr[len] = '\n';
    len++;

    int nl_dec = get_size(nstr, len);
    int nl = (nl_dec * 4) / W + 1;

    ui d = (ui)malloc(sizeof(ui_t) * nl * SIZE);
    ui n = str2arr(nstr, len, nl_dec, nl);
    ui_t d_[nl];

    ecm(d, n, nl);

    for(i = 0; i < SIZE; i++) {
        for(j = 0; j < nl; j++) {
            d_[j] = d[dIndex];
            dIndex++;
        }
        mpz_import(mp_d, nl, -1, 4, 0, 0, d_);
        res = mpz_get_str(res, 10, mp_d);

        fprintf(fp, "%s\n", res);
    }
    free(n);
    free(d);

    return res;
}


int main(int argc, char *argv[]) {
    clock_t start, end;
    start = clock();

    // pro_curve_point_gmp_test(10000);
    // pro_add_gmp_test(10000);
    // pro_add_magma_test(10000);
    // pro_dbl_magma_test(10000);
    // pro_ladder_gmp_test(10000);
    // pro_ladder_magma_test(100);
    // ecm_gmp_test(100);

    SIZE = 1;

    if(argc == 1) {                                             // fiwecm
        printf("Invalid arguments. See fiwecm --help.\n");
        return 1;
    }
    if(argc == 2 && !strcmp(argv[1], "--help")) {               // fiwecm --help
        printf("\n /5asena55 /555555 /55      /55 /5furkan5     \n"
                 "| 55_____/|_  55_/| 55  /5 | 55| 55_____/     \n"
                 "| 55        | 55  | 55 /555| 55| 55           \n"     
                 "| 55555     | 55  | 55/55 55 55| 55555        \n"   
                 "| 55__/     | 55  | ozan_  5555| 55__/        \n"   
                 "| 55        | 55  | 555/ \\  555| 55          \n"      
                 "| 55       /5aslı5| 55/   \\  55| 55elif55    \n"
                 "|__/      |______/|__/     \\__/|________/    \n\n");
        printf("Usage: fiwecm [n] [options]\n\n");
        printf("Parameters:\n");
        printf("%3c  n %11c number to be factorized in base 10\n\n", ' ', ' ');
        printf("Options:\n");
        printf("%3c -f filename %1c read numbers to be factorized from file \"filename\"\n", ' ', ' ');
        printf("%3c -o filename %1c output the factors to file \"filename\"\n\n", ' ', ' ');
        printf("Examples:\n");
        printf("%3c  fiwecm 126378132 %13c factorize 126378132 and print the factor to console\n", ' ', ' ');
        printf("%3c  fiwecm 126378132 -o out.txt %2c factorize 126378132 and print the factor to out.txt\n", ' ', ' ');
        printf("%3c  fiwecm -f in.txt %13c factorize the numbers in in.txt and print the factors to console\n", ' ', ' ');
        printf("%3c  fiwecm -f in.txt -o out.txt %2c factorize the numbers in in.txt and print the factors to out.txt\n\n", ' ', ' ');
        return 1;
    }

    FILE *fp;
    int i = 0, len;

    if(strcmp(argv[1], "-f")) {                                 // fiwecm num
        do {
            if(((int)(argv[1][i]) < 48) || ((int)(argv[1][i]) > 57)) {
                printf("Invalid number.\n");
                return 0;
            }
            i++;
        } while(argv[1][i] != '\0');
        len = i;

        if(argc >= 3 && !strcmp(argv[2], "-o")) {               // fiwecm num -o out.txt
            fp = fopen(argv[3], "w");
            factorize(argv[1], len, fp);
            fclose(fp);
            return 1;
        } else {   
            factorize(argv[1], len, stdout);
            return 1;
        }
    } else {    
        if(argc <= 2) {
            printf("Please insert input file name.\n");
            return 0;
        } else {
            char *nstr = txt2str(argv[2]);
            do {
                if(((int)(nstr[i]) < 48) || ((int)(nstr[i]) > 57)) {
                    if(nstr[i] != '\n') {
                        printf("Invalid number\n");
                        return 0;
                    }
                }
                i++;
            } while(nstr[i] != '\0');
            len = i;

            if(argc > 3 && !strcmp(argv[3], "-o")) {          // fiwecm -f in.txt -o out.txt
                fp = fopen(argv[4], "w");
                factorize(nstr, len, fp);
                fclose(fp);
                return 1;
            } else {
                factorize(nstr, len, stdout);
                return 1;
            }
        }
    }
   
    end = clock();
    printf("// Time spent: %.4fms\n", (float)((end - start) / CLOCKS_PER_SEC));

    return 0;
}