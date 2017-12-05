#ifndef MATRIX_H
#define MATRIX_H

struct Matrix{
    int rows_;
    int cols_;
    double* bufs_;
};

struct Matrix* create_matrix(int rows, int cols);

void destroy_matrix(struct Matrix* mat);

void set_value(int row, int cols, double value);

double get_value(int row, int cols);

struct Matrix* transepose();

double* get_bufs();
#endif /* MATRIX_H */

