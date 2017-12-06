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

void set_mat_value(struct Matrix* mat, double value, int row, int col) {
    int index = row * mat->cols + col;
    mat->bufs[index] = value;
};

void set_mat_value_by_index(struct Matrix* mat, double value, int index){
    mat->bufs[index] = value;
}

double get_mat_value(const struct Matrix* mat, int row, int col) {
    int index = row * mat->cols + col;
    return mat->bufs[index];
}

double get_mat_value_by_index(const struct Matrix* mat, int index){
    return mat->bufs[index];
}

double* get_mat_bufs(const struct Matrix* mat){
    if(mat == NULL)
        return NULL;
    else
        return mat->bufs;
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
    printf("--------------------\n");
    for (i = 0; i < mat->rows; i++) {
        for (j = 0; j < mat->cols; j++) {
            printf("%f ", mat->bufs[i * mat->cols + j]);
        }
        printf("\n");
    }
    printf("--------------------\n");
};

struct Matrix* transpose_mat(const struct Matrix* original_mat) {

    struct Matrix* transposed_mat = create_mat(original_mat->cols, original_mat->rows);
    int i, j;
    for (i = 0; i < original_mat->rows; i++) {
        for (j = 0; j < original_mat->cols; j++) {
            double value = get_mat_value(original_mat, i, j);
            set_mat_value(transposed_mat, value, j, i);
        }
    }

    return transposed_mat;
};


struct Matrix* read_mat_file(const char* file_path) {
    FILE* fp;
    int rows, cols;
    double value;
    struct Matrix* mat = NULL;
    int i, j;
    
#ifdef NO_MAT_FILE
    int mat_rows = 2;
    matrix = (double *) malloc(mat_rows * mat_rows * sizeof (double));
    for (i = 0; i < mat_rows * mat_rows; i++) {
        matrix[i] = i;
    }
    *rows_num = mat_rows;
    *cols_num = mat_rows;
    return matrix;
#else

    fp = fopen(file_path, "rb");
    if (fp != NULL) {
        fread(&rows, sizeof (int), 1, fp);
        fread(&cols, sizeof (int), 1, fp);
        
        mat = create_mat(rows, cols);
        mat->rows = rows;
        mat->cols = cols;

        for (i = 0; i < rows; ++i) {
            for (j = 0; j < cols; ++j) {
                fread(&value, sizeof (double), 1, fp);
                set_mat_value(mat, value, i, j);
            }
        }
    } else {
        return NULL;
    }
    fclose(fp);

    return mat;
#endif
}


void write_mat_file(const char* file_path, struct Matrix* mat) {
    FILE* fp;
    int i, j;
    double value;

    int rows = mat->rows;
    int cols = mat->cols;
    fp = fopen(file_path, "wb");
    if (fp != NULL) {

        fwrite(&rows, sizeof (int), 1, fp);
        fwrite(&cols, sizeof (int), 1, fp);
        for (i = 0; i < rows; ++i) {
            for (j = 0; j < cols; ++j) {
                value = get_mat_value(mat, i, j);
                fwrite(&value, sizeof (double), 1, fp);
            }
        }
    }
    fclose(fp);
}