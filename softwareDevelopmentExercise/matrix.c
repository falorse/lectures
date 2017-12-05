#include <stdio.h>
#include "matrix.h"

struct Matrix* create_mat(int rows, int cols) {
    struct Matrix* mat = NULL;

    mat = malloc(sizeof (struct Matrix));
    if (mat == NULL) {
        return NULL;
    }

    mat->rows = rows;
    mat->cols = cols;
    mat->bufs = NULL;

    double* buf = (double*) malloc(sizeof (double) * rows * cols);
    if (buf == NULL) {
        destroy_mat(mat);
        return NULL;
    }
    mat->bufs = buf;

    return mat;
};

void set_mat_value(struct Matrix* mat, int row, int col, double value) {
    int index = row * mat->cols + col;
    mat->bufs[index] = value;
};

double get_mat_value(struct Matrix* mat, int row, int col) {
    int index = row * mat->cols + col;
    return mat->bufs[index];
}

void destroy_mat(struct Matrix* mat) {
    if (mat != NULL) {
        if (mat->bufs != NULL) {
            free(mat->bufs);
        }
        mat->bufs = NULL;

        free(mat);
        mat = NULL;
    }
};

void show_mat(struct Matrix* mat) {
    int i, j;
    for (i = 0; i < mat->rows; i++) {
        for (j = 0; j < mat->cols; j++) {
            printf("%f ", mat->bufs[i * mat->cols + j]);
        }
        printf("\n");
    }
};

struct Matrix* transepose_mat(const struct Matrix* original_mat) {

    struct Matrix* transeposed_mat = create_mat(original_mat->cols, original_mat->rows);
    int i, j;
    for (i = 0; i < original_mat->rows; i++) {
        for (j = 0; j < original_mat->cols; j++) {
            double value = get_mat_value(original_mat, i, j);
            set_mat_value(transeposed_mat, j, i, value);
        }
    }

    return transeposed_mat;
};
