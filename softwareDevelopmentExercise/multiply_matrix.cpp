#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

#define MAXLEN 1024
#define OUT 0

double* read_mat_file(const char* path, int* rows_num, int* cols_num);
//double** readMatrix(const char* path, int* rows_num, int* cols_num);

void write_mat_file(const char* path, double* matrix, int rows_num, int cols_num);
double* mulutipy_mat(const double* leftMat, const double* rightMat);
double* transposeMat(const double* originalMat);
void distributeRightMat();

int main(int argc, char* argv[]) {

    // 引数チェックとファイルパスの取得
    char A_path[MAXLEN];
    char B_path[MAXLEN];
    char C_path[MAXLEN];
    if (argc == 4) {
        strcpy(A_path, argv[1]);
        strcpy(B_path, argv[2]);
        strcpy(C_path, argv[3]);
    } else {
        printf("ERROR USAGE: mul_matrix leftMat_path rightMat_path ansMat_path\n");
        return 0;
    }

    double start, start2, end, end2;

    // MPI Initialize
    int rank, num_procs, alloc_A_elems_num, alloc_A_rows_num, alloc_C_elems_num;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // 行列の読み取り
    int A_rows, B_cols, A_cols;

    if (rank == 0) {
        double* A = read_mat_file(A_path, &A_rows, &A_cols);
        double* B = read_mat_file(B_path, &A_cols, &B_cols);
        double* transeposedB = transposeMat(B);
    }
}