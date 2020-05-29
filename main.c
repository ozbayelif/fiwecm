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
    int i, j, k, was_newline;
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
            was_newline = 0;
        } else {
            if(was_newline == 0) {
                num = realloc(num, sizeof(char) * (j + 1));
                num[j] = '\0';

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
            was_newline = 1;
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
    if(fp == NULL) {
        printf("Please enter a valid input file name.\n");
        exit(1);
    }
    fseek(fp, 0, SEEK_END);
    long len = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char *nstr = (char *)malloc(sizeof(char) * (len + 1));
    fread(nstr, 1, len, fp);
    fclose(fp);

    return nstr;
}

void factorize(char *nstr1, int len, FILE *fp) {
    int i, j, dIndex = 0;
    mpz_t mp_d;

    char *nstr = (char *)malloc(sizeof(char) * (len + 1)); // for num + '\n'
    char *res = (char *)malloc(sizeof(char) * ((len + 1) * SIZE)); // for factors
    mpz_init(mp_d);

    strcpy(nstr, nstr1);
    nstr1[len] = '\n';
    len++;

    int nl_dec = get_size(nstr1, len);
    int nl = (nl_dec * 4) / W + 1;

    ui d = (ui)malloc(sizeof(ui_t) * nl * SIZE);
    ui n = str2arr(nstr1, len, nl_dec, nl);
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
}

int bash_main(int argc, char *argv[]) {
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
    char *nstr;
    int i = 0, len, out_file;

    if(!strcmp(argv[1], "-f")) {                                                // file input
        if(argc >= 3) {
            nstr = txt2str(argv[2]);
        } else {
            printf("Please enter a file name.\n");
            return 0;
        }
    } else {                                                                    // normal input
        nstr = argv[1];
    }

    do {
        if(((nstr[i]) < 48) || ((nstr[i]) > 57)) {
            printf("Please enter a valid number.\n");
            return 0;
        }
        if(nstr[i] == '\n') {
            SIZE++;
            continue;
        }
        i++;
    } while(nstr[i] != '\0');
    len = i;

    if(argc == 4 && !strcmp(argv[2], "-o")) {
        fp = fopen(argv[3], "w");
        out_file = 1;
    } else if(argc == 5 && !strcmp(argv[3], "-o")) {                            // output to file
        fp = fopen(argv[4], "w");
        out_file = 1;
    } else if(argc == 3 && !strcmp(argv[2], "-o")) {
        printf("Please enter a file name.\n");
        return 0;
    } else if(argc == 4 && !strcmp(argv[3], "-o")) {
        printf("Please enter a file name.\n");
        return 0;
    } else {                                                                    // output to console
        fp = stdout;
        out_file = 0;
    }
    printf("Factoring %s...\n\n", nstr);
    if(out_file == 0) {
        printf("Factor found: ");
    } else {
        if(fp == NULL) {
            printf("Please enter a valid output file name.\n");
            return 0;
        }
    }

    factorize(nstr, len, fp);

    if(out_file == 1) {
        fclose(fp);
    }

    return 1;
}


int main(int argc, char *argv[]) {
    clock_t start, end;
    start = clock();
    srand(time(NULL));

    // pro_curve_point_gmp_test(10000);
    // pro_add_gmp_test(10000);
    // pro_add_magma_test(10000);
    // pro_dbl_magma_test(10000);
    // pro_ladder_gmp_test(10000);
    // pro_ladder_magma_test(100);
    // ecm_gmp_test(1000);

    bash_main(argc, argv);

    end = clock();
    printf("\nTime spent: %.4fms\n\n", (float)((end - start) / CLOCKS_PER_SEC));

    return 0;
}