#ifndef MATRIX_H
#define MATRIX_H

struct Matrix{
    int rows;
    int cols;
    double* bufs;
};

struct Matrix* create_mat(int rows, int cols);

void destroy_mat(struct Matrix* mat);

void set_mat_value(struct Matrix*, double value, int row, int cols);

void set_mat_value_by_index(struct Matrix*, double value, int index);

double get_mat_value(const struct Matrix* mat, int row, int cols);

double get_mat_value_by_index(const struct Matrix* mat, int index);

double* get_mat_bufs(const struct Matrix* mat);

void show_mat(struct Matrix* mat);

struct Matrix* transpose_mat(const struct Matrix* original_mat);

struct Matrix* read_mat_file(const char* file_path);

void write_mat_file(const char* file_path, struct Matrix* mat);

#endif /* MATRIX_H */

