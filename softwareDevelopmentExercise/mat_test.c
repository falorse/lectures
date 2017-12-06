#include "matrix.h"
#include <stdio.h>

int main(int argc, char* argv) {
    struct Matrix* A = create_mat(32, 32);
    struct Matrix* B = create_mat(2, 3);

    int i, j;
    for (i = 0; i < 32; i++) {
        for (j = 0; j < 32; j++) {
            set_mat_value(A, i * A->cols + j, i, j);
        }
    }
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 3; j++) {
            set_mat_value(B, i * B->cols + j, i, j);
            printf("%d,%f\n", i * B->cols + j, get_mat_value(B, i, j));
        }
    }

    show_mat(A);
    show_mat(B);

    struct Matrix* t_B = transpose_mat(B);
    show_mat(t_B);
    
    write_mat_file("./matrix_files/Mat32", A);
}
