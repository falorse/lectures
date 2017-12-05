#ifndef MATRIX_H
#define MATRIX_H

struct Matrix{
    int rows;
    int cols;
    double* bufs;
};

struct Matrix* create_mat(int rows, int cols);

void destroy_mat(struct Matrix* mat);

void set_mat_value(struct Matrix*, int row, int cols, double value);

double get_mat_value(struct Matrix* mat, int row, int cols);

void show_mat(struct Matrix* mat);

struct Matrix* transepose_mat();

#endif /* MATRIX_H */

