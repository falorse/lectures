#include "matrix.h"
#include <stdio.h>

int main(int argc, char* argv) {
    struct Matrix* A = create_mat(2, 2);
    struct Matrix* B = create_mat(2, 3);

    int i,j;
    for(i = 0; i < 2; i++){
        for(j = 0; j < 2; j++){
            set_mat_value(A, i, j, i * A->cols + j);
        }
    }
    for(i = 0; i < 2; i++){
        for(j = 0; j < 3; j++){
            set_mat_value(B, i, j, i * B->cols + j);
            printf("%d,%f\n", i * B->cols + j, get_mat_value(B,i,j));
        }
    }

    show_mat(A);
    show_mat(B);

    struct Matrix* t_B = transepose_mat(B);
    show_mat(t_B);
}
