#include <stdio.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <zconf.h>
#include "include/bam.h"

void print_matrix(const gsl_matrix *M, FILE *file);

void print_help() {
    printf("Usage: BAM c <N> <M> <P> <file>\n"
                   "  or:  BAM a|d <fileA> <fileB> <file>\n"
                   "  or:  BAM l|r <fileA> <fileB> <file>\n"
                   "  or:  BAM --help\n\n"
                   "  c\tcreate new BAM in <file>, where\n"
                   "\t  N - height of weight matrix\n"
                   "\t  M - width of matrix A and height of matrix B\n"
                   "\t  P - width of weight matrix\n"
                   "  a\tadd a pair from <fileA> and <fileB> in BAM from <infile>\n"
                   "  d\tdelete a pair from <fileA> and <fileB> from BAM from <infile>\n"
                   "  l\tget association of image X from <file>\n"
                   "\t  from BAM in <infile>\n"
                   "  r\tget association of image Y from <file>\n"
                   "\t  from BAM in <infile>\n"
    );
}

int main(int argc, char **argv) {
    const char *usage = "Usage: %s c <N> <M> <P> <file>\n"
            "  or:  %s a|d <fileA> <fileB> <file>\n"
            "  or:  %s l|r <fileA> <fileB> <file>\n"
            "  or:  %s --help\n";
    if(!(argc == 2 || argc == 5 || argc == 6)) {
        printf(usage, argv[0], argv[0], argv[0], argv[0]);
        return 1;
    }

    if (argc == 2) {
        if (!strcmp(argv[1], "--help")) {
            print_help();
            return 0;
        } else {
            fprintf(stderr, usage, argv[0], argv[0]);
            return 1;
        }
    }

    char option = 0;

    if(strlen(argv[1]) > 1) {
        fprintf(stderr, "wrong option!\n");
        return 1;
    }
    option = *argv[1];

    BAM_model *bam_model = NULL;
    FILE *fileA = NULL;
    FILE *fileB = NULL;
    gsl_matrix *A = NULL;
    gsl_matrix *B = NULL;
    size_t N = 0;
    size_t M = 0;
    size_t P = 0;
    size_t nA = 0;
    size_t mA = 0;
    size_t nB = 0;
    size_t mB = 0;
    int ret_val = 0;

    bam_model = malloc(sizeof(BAM_model));

    switch (option){
        case 'c':
            if(argc != 6) {
                printf(usage, argv[0], argv[0], argv[0], argv[0]);
                return 1;
            }
            N = strtoull(argv[2], NULL, 10);
            M = strtoull(argv[3], NULL, 10);
            P = strtoull(argv[4], NULL, 10);

            ret_val = BAM_create(argv[5], N, M, P);
            if(FILE_ERR == ret_val) {
                fprintf(stderr, "error reading of file!\n");
                return 1;
            }

            printf("useful capacity: %d\n",
                   BAM_max_vectors_count((int) (MIN(N * M, M * P))));

            break;

        case 'a':
            if(argc != 5) {
                printf(usage, argv[0], argv[0], argv[0], argv[0]);
                return 1;
            }
            if(NULL == (fileA = fopen(argv[2], "r"))) {
                fprintf(stderr, "could not open: %s\n", argv[2]);
                return 1;
            };
            if(NULL == (fileB = fopen(argv[3], "r"))) {
                fprintf(stderr, "could not open: %s\n", argv[2]);
                return 1;
            };

            ret_val = BAM_load(bam_model, argv[4]);
            if(FILE_ERR == ret_val){
                fprintf(stderr, "could not open: %s\n", argv[4]);
                return 1;
            }

            fscanf(fileA, "%zu %zu\n", &nA, &mA);
            fscanf(fileB, "%zu %zu\n", &nB, &mB);

            A = gsl_matrix_alloc(nA, mA);
            B = gsl_matrix_alloc(nB, mB);

            gsl_matrix_fscanf(fileA, A);
            gsl_matrix_fscanf(fileB, B);

            ret_val = BAM_add(bam_model, A, B);
            if(ALG_ERR == ret_val) {
                fprintf(stderr, "matrices have wrong sizes\n");
                return 1;
            }

            printf("useful capacity left: %d\n", bam_model->useful_capacity);

            BAM_destroy(bam_model);
            gsl_matrix_free(A);
            gsl_matrix_free(B);
            fclose(fileA);
            fclose(fileB);
            break;

        case 'd':
            if(argc != 5) {
                printf(usage, argv[0], argv[0], argv[0], argv[0]);
                return 1;
            }
            if(NULL == (fileA = fopen(argv[2], "r"))) {
                fprintf(stderr, "could not open: %s\n", argv[2]);
                return 1;
            };
            if(NULL == (fileB = fopen(argv[3], "r"))) {
                fprintf(stderr, "could not open: %s\n", argv[2]);
                return 1;
            };

            ret_val = BAM_load(bam_model, argv[4]);
            if(FILE_ERR == ret_val){
                fprintf(stderr, "could not open: %s\n", argv[4]);
                return 1;
            }

            fscanf(fileA, "%zu %zu\n", &nA, &mA);
            fscanf(fileB, "%zu %zu\n", &nB, &mB);

            A = gsl_matrix_alloc(nA, mA);
            B = gsl_matrix_alloc(nB, mB);

            gsl_matrix_fscanf(fileA, A);
            gsl_matrix_fscanf(fileB, B);

            ret_val = BAM_delete(bam_model, A, B);
            if(ALG_ERR == ret_val) {
                fprintf(stderr, "matrices have wrong sizes\n");
                return 1;
            }

            printf("useful capacity left: %d\n", bam_model->useful_capacity);

            gsl_matrix_free(A);
            gsl_matrix_free(B);
            fclose(fileA);
            fclose(fileB);
            break;

        case 'l':
            if(argc != 5) {
                printf(usage, argv[0], argv[0], argv[0], argv[0]);
                return 1;
            }
            if(NULL == (fileA = fopen(argv[2], "r"))) {
                fprintf(stderr, "could not open: %s\n", argv[2]);
                return 1;
            };

            printf("--associating\n");

            BAM_load(bam_model, argv[4]);

            fscanf(fileA, "%zu %zu", &nA, &mA);
            if(nA * mA != bam_model->W->size1) {
                fprintf(stderr, "sizes of matrix A do not fit to the BAM\n");
                return 1;
            }

            A = gsl_matrix_alloc(nA, mA);
            gsl_matrix_fscanf(fileA, A);

            B = gsl_matrix_alloc(bam_model->M, bam_model->P);
            if(NULL == B) {
                fprintf(stderr, "internal memory error\n");
                return 1;
            }

            ret_val = BAM_associate_left(bam_model, A, B);
            if(ALG_ERR == ret_val) {
                fprintf(stderr, "linear algebra error\n");
                return 1;
            }

            if(NULL == (fileB = fopen(argv[3], "w"))) {
                fprintf(stderr, "could not open: %s\n", argv[3]);
                return 1;
            };

            print_matrix(B, fileB);

            BAM_destroy(bam_model);
            gsl_matrix_free(A);
            gsl_matrix_free(B);
            fclose(fileA);
            fclose(fileB);

            break;

        case 'r':
            if(argc != 5) {
                printf(usage, argv[0], argv[0], argv[0], argv[0]);
                return 1;
            }
            if(NULL == (fileB = fopen(argv[2], "r"))) {
                fprintf(stderr, "could not open: %s\n", argv[3]);
                return 1;
            };

            printf("--associating\n");

            BAM_load(bam_model, argv[4]);

            fscanf(fileB, "%zu %zu", &nB, &mB);
            if(nB * mB != bam_model->W->size2) {
                fprintf(stderr, "sizes of matrix B do not fit to the BAM\n");
                return 1;
            }

            B = gsl_matrix_alloc(nB, mB);
            gsl_matrix_fscanf(fileB, B);

            A = gsl_matrix_alloc(bam_model->N, bam_model->M);
            if(NULL == A) {
                fprintf(stderr, "internal memory error\n");
                return 1;
            }
            ret_val = BAM_associate_right(bam_model, A, B);
            if(ALG_ERR == ret_val) {
                fprintf(stderr, "linear algebra error\n");
                return 1;
            }

            if(NULL == (fileA = fopen(argv[3], "w"))) {
                fprintf(stderr, "could not open: %s\n", argv[2]);
                return 1;
            };

            print_matrix(A, fileA);

            BAM_destroy(bam_model);
            gsl_matrix_free(A);
            gsl_matrix_free(B);
            fclose(fileA);
            fclose(fileB);

            break;
        default:
            fprintf(stderr, "wrong option!\n");
            return 1;
    }

    free(bam_model);

    return 0;
}

void print_matrix(const gsl_matrix *M, FILE *file) {
    int number = 0;

    fprintf(file, "%zu %zu\n", M->size1, M->size2);

    for(size_t i = 0; i < M->size1; ++i) {
        for (size_t j = 0; j < M->size2; ++j) {
            number = (int) gsl_matrix_get(M, i, j);
            fprintf(file, "%d ", number == 1 ? 1 : 0);
        }
        fprintf(file, "\n");
    }
}
