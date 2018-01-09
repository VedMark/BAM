#ifndef ICOMPRESSOR_COMPRESSOR_H
#define ICOMPRESSOR_COMPRESSOR_H


#include <stdbool.h>
#include <gsl/gsl_blas.h>

#define PRINT_MATRIX(M) {                               \
    printf("MxN: %lu %lu\n", (M)->size1, (M)->size2);   \
    for(int i = 0; i < (M)->size1; ++i) {               \
        for(int j = 0; j < (M)->size2; ++j) {           \
            printf("%lf ", gsl_matrix_get((M), i, j));  \
        }                                               \
    printf("\n");                                       \
    }                                                   \
}

#define PRINT_VECTOR(V) {                           \
    printf("M: %lu\n", (V)->size);                  \
    for(int i = 0; i < (V)->size; ++i) {            \
        printf("%lf ", gsl_vector_get((V), i));     \
    }                                               \
    printf("\n");                                   \
}

#define MIN(a, b) (a) < (b) ? (a) : (b)

#define MEM_ERR (1)
#define ALG_ERR (-1)
#define FILE_ERR (-2)
#define SUCCESS (0)

typedef struct BAM_model {
    gsl_matrix *X;
    gsl_matrix *Y;
    gsl_matrix *W;
    char *file_name;
    FILE *matrix_file;
    size_t N;
    size_t M;
    size_t P;
    int useful_capacity;
} BAM_model;

int BAM_max_vectors_count(int n);
int BAM_create(char *file_name, size_t N, size_t M, size_t P);
int BAM_load(BAM_model *model, char *file_name);
void BAM_destroy(BAM_model *model);
int BAM_add(BAM_model *model, gsl_matrix *A, gsl_matrix *B);
int BAM_delete(BAM_model *model, gsl_matrix *A, gsl_matrix *B);
int BAM_associate_left(BAM_model *model, const gsl_matrix *A, gsl_matrix *B);
int BAM_associate_right(BAM_model *model, gsl_matrix *A, const gsl_matrix const *B);

#endif //ICOMPRESSOR_COMPRESSOR_H
