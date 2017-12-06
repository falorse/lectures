#include <stdio.h>
#include "matrix.h"

struct Matrix* create_mat(int rows, int cols) {
    struct Matrix* mat = NULL;

    mat = malloc(sizeof (struct Matrix));

    mat->rows = rows;
    mat->cols = cols;

    double* buf = (double*) malloc(sizeof (double) * rows * cols);
    mat->bufs = buf;

    return mat;
};

struct Matrix* create_mat_without_bufs(int rows, int cols){
    struct Matrix* mat = NULL;
    
    mat = malloc(sizeof (struct Matrix));
    
    mat->rows = rows;
    mat->cols = cols;
    double* buf = (double*) malloc(sizeof (double) * 0);
    
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

void add_mat_value(struct Matrix* mat, double value, int row, int col){
    double pre_value = get_mat_value(mat, row, col);
    set_mat_value(mat, pre_value + value, row, col);
}

double* get_mat_bufs(const struct Matrix* mat){
    if(mat == NULL)
        return NULL;
    else
        return mat->bufs;
}

int get_mat_size(const struct Matrix* mat){
    return mat->rows * mat->cols;
};

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
            printf("%f ", get_mat_value(mat, i, j));
        }
        printf("\n");
    }
    printf("--------------------\n");
};

struct Matrix* transepose_mat(const struct Matrix* original_mat) {

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