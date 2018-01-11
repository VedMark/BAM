#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include "include/bam.h"

#define LAMBDA (5)

void activate_m(gsl_matrix *Y);
double F(double S);
int read_matrix(gsl_matrix **matrix, FILE *file, size_t n, size_t m);
void transform_to_bipolar_vector(const gsl_matrix const *A, gsl_matrix *X);
void transform_to_matrix(const gsl_matrix const *X, gsl_matrix *A);
void write_matrix(const gsl_matrix *matrix, FILE *file);

double F(double S) {
    return 1 / (1 + exp(-LAMBDA * S)) >= 0.5 ? 1 : -1;
}

int BAM_max_vectors_count(int n) {
    int res = (int) (n / (2. * log2(n)));
    return res < 0 ? 0 : res;
}

int BAM_create(char *file_name, size_t N, size_t M, size_t P) {
    FILE *file = NULL;
    int useful_capacity = 0;
    int ret_val = 0;

    if(NULL == (file = fopen(file_name, "w"))) return FILE_ERR;

    useful_capacity = BAM_max_vectors_count((int) (MIN(N * M, M * P)));

    ret_val = fprintf(file, "%d\n%zu %zu %zu\n", useful_capacity, N, M, P);
    if(0 == ret_val) return FILE_ERR;

    for(size_t i = 0; i < N * M; ++i) {
        for(size_t j = 0; j < M * P; ++j) {
            ret_val = fprintf(file, "%d ", 0);
            if(0 == ret_val) return FILE_ERR;
        }
        ret_val = fprintf(file, "\n");
        if(0 == ret_val) return FILE_ERR;
    }

    fclose(file);

    return SUCCESS;
}

int BAM_load(BAM_model *model, char *file_name)
{
    int ret_val = 0;

    model->matrix_file = fopen(file_name, "r");
    if(NULL == model->matrix_file) return FILE_ERR;

    fscanf(model->matrix_file, "%d\n", &model->useful_capacity);
    fscanf(model->matrix_file, "%zu %zu %zu\n", &model->N, &model->M, &model->P);

    ret_val = read_matrix(&model->W, model->matrix_file,
                          model->N * model->M, model->M * model->P);
    if(FILE_ERR == ret_val) return FILE_ERR;
    if(MEM_ERR == ret_val) return MEM_ERR;

    model->X = gsl_matrix_alloc(1, model->W->size1);
    if(NULL == model->X) return MEM_ERR;

    model->Y = gsl_matrix_alloc(1, model->W->size2);
    if(NULL == model->Y) return MEM_ERR;

    return SUCCESS;
}

void BAM_destroy(BAM_model *model) {
    gsl_matrix_free(model->X);
    gsl_matrix_free(model->Y);
    gsl_matrix_free(model->W);
    fclose(model->matrix_file);
}

int BAM_add(BAM_model *model, gsl_matrix *A, gsl_matrix *B) {
    int ret_val = 0;

    if(A->size1 * A->size2 != model->W->size1) return ALG_ERR;
    if(B->size1 * B->size2 != model->W->size2) return ALG_ERR;

    transform_to_bipolar_vector(A, model->X);
    transform_to_bipolar_vector(B, model->Y);

    PRINT_MATRIX(model->X);
    PRINT_MATRIX(model->Y);

    ret_val = gsl_blas_dgemm(CblasTrans, CblasNoTrans,
                             1, model->X, model->Y,
                             1, model->W);
    if(0 != ret_val) return ALG_ERR;

    PRINT_MATRIX(model->W);

    model->useful_capacity--;

    freopen(model->file_name, "w", model->matrix_file);
    fprintf(model->matrix_file, "%d\n", model->useful_capacity);
    fprintf(model->matrix_file, "%zu %zu %zu\n", model->N, model->M, model->P);
    write_matrix(model->W, model->matrix_file);

    return SUCCESS;
}

void transform_to_bipolar_vector(const gsl_matrix const *A, gsl_matrix *X){
    for(size_t i = 0; i < A->size1; ++i) {
        for(size_t j = 0; j < A->size2; ++j) {
            gsl_matrix_set(X, 0, i * A->size2 + j,
                           gsl_matrix_get(A, i, j) ? 1 : -1);
        }
    }
}

void transform_to_matrix(const gsl_matrix const *X, gsl_matrix *A) {
    for(size_t i = 0; i < A->size1; ++i) {
        for(size_t j = 0; j < A->size2; ++j) {
            gsl_matrix_set(A, i, j, gsl_matrix_get(X, 0,  i * A->size2 + j));
        }
    }
}

int read_matrix(gsl_matrix **matrix, FILE *file, size_t n, size_t m) {
    int number = 0;

    *matrix = gsl_matrix_alloc(n, m);
    if(NULL == *matrix) return MEM_ERR;

    for(size_t i = 0; i < (*matrix)->size1; ++i) {
        for(size_t j = 0; j < (*matrix)->size2; ++j) {
            fscanf(file, "%d ", &number);
            gsl_matrix_set(*matrix, i, j, number);
        }
        fscanf(file, "\n");
    }

    return SUCCESS;
}

void write_matrix(const gsl_matrix const *matrix, FILE *file) {
    int number = 0;

    for(size_t i = 0; i < matrix->size1; ++i) {
        for(size_t j = 0; j < matrix->size2; ++j) {
            number = (int) gsl_matrix_get(matrix, i, j);
            fprintf(file, "%d ", number);
        }
        fprintf(file, "\n");
    }
}

int BAM_delete(BAM_model *model, gsl_matrix *A, gsl_matrix *B) {
    int ret_val = 0;

    if(A->size1 * A->size2 != model->W->size1) return ALG_ERR;
    if(B->size1 * B->size2 != model->W->size2) return ALG_ERR;

    transform_to_bipolar_vector(A, model->X);
    transform_to_bipolar_vector(B, model->Y);

    ret_val = gsl_blas_dgemm(CblasTrans, CblasNoTrans,
                             -1, model->X, model->Y,
                             1, model->W);
    if(0 != ret_val) return ALG_ERR;

    model->useful_capacity++;

    freopen(model->file_name, "w", model->matrix_file);
    fprintf(model->matrix_file, "%d\n", model->useful_capacity);
    fprintf(model->matrix_file, "%zu %zu %zu\n", model->N, model->M, model->P);
    write_matrix(model->W, model->matrix_file);

    return SUCCESS;
}

int BAM_associate_left(BAM_model *model, const gsl_matrix const*A, gsl_matrix *B) {
    int ret_val = 0;
    if(A->size1 * A->size2 != model->W->size1) return ALG_ERR;
    transform_to_bipolar_vector(A, model->X);

    gsl_matrix *X_prev = NULL;
    gsl_matrix *X_prev_prev = NULL;
    int iter = 0;

    X_prev = gsl_matrix_alloc(model->X->size1, model->X->size2);
    if(NULL == X_prev) return MEM_ERR;

    X_prev_prev = gsl_matrix_alloc(model->X->size1, model->X->size2);
    if(NULL == X_prev_prev) return MEM_ERR;

    do {
        ret_val = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                                 1, model->X, model->W,
                                 0, model->Y);
        if (0 != ret_val) return ALG_ERR;

        gsl_matrix_memcpy(X_prev_prev, X_prev);
        gsl_matrix_memcpy(X_prev, model->X);

        ret_val = gsl_blas_dgemm(CblasNoTrans, CblasTrans,
                                 1, model->Y, model->W,
                                 0, model->X);
        if (0 != ret_val) return ALG_ERR;

        activate_m(model->X);

        iter++;
    } while(!(gsl_matrix_equal(X_prev, model->X)
              || gsl_matrix_equal(X_prev_prev, model->X)));

    transform_to_matrix(model->Y, B);

    printf("iterations: %d\n", iter);

    gsl_matrix_free(X_prev);
    gsl_matrix_free(X_prev_prev);

    return SUCCESS;
}

void activate_m(gsl_matrix *Y) {
    for(size_t i = 0; i < Y->size1; ++i) {
        for(size_t j = 0; j < Y->size2; ++j) {
            gsl_matrix_set(Y, i, j, F(gsl_matrix_get(Y, i, j)));
        }
    }
}

int BAM_associate_right(BAM_model *model, gsl_matrix *A, const gsl_matrix const *B) {
    int ret_val = 0;
    if(B->size1 * B->size2 != model->W->size2) return ALG_ERR;
    transform_to_bipolar_vector(B, model->Y);

    gsl_matrix *Y_prev = NULL;
    gsl_matrix *Y_prev_prev = NULL;
    int iter = 0;

    Y_prev = gsl_matrix_alloc(model->Y->size1, model->Y->size2);
    if(NULL == Y_prev) return MEM_ERR;

    Y_prev_prev = gsl_matrix_alloc(model->Y->size1, model->Y->size2);
    if(NULL == Y_prev_prev) return MEM_ERR;

    do {
        ret_val = gsl_blas_dgemm(CblasNoTrans, CblasTrans,
                                 1, model->Y, model->W,
                                 0, model->X);
        if (0 != ret_val) return ALG_ERR;

        activate_m(model->X);

        gsl_matrix_memcpy(Y_prev_prev, Y_prev);
        gsl_matrix_memcpy(Y_prev, model->Y);

        ret_val = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                                 1, model->X, model->W,
                                 0, model->Y);
        if (0 != ret_val) return ALG_ERR;

        activate_m(model->Y);

        iter++;

    } while(!(gsl_matrix_equal(Y_prev, model->Y)
              || gsl_matrix_equal(Y_prev_prev, model->Y)));

    printf("iterations: %d\n", iter);

    transform_to_matrix(model->X, A);

    gsl_matrix_free(Y_prev);
    gsl_matrix_free(Y_prev_prev);

    return SUCCESS;
}
