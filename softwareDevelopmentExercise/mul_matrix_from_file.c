
#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>
#include <omp.h>
#include "matrix.h"

#define MAXLEN 1024
#define OUT 0
//#define NO_MAT_FILE

double* read_mat_file(const char* path, int* rows_num, int* cols_num);
//double** readMatrix(const char* path, int* rows_num, int* cols_num);

void write_mat_file(const char* path, double* matrix, int rows_num, int cols_num);
double* mulutipy_mat(const double* leftMat, const double* rightMat);
void broadcast_mat_size(int* A_rows_num,int* B_cols_num,int* A_cols_num,int* alloc_C_elems_num,int* alloc_A_rows_num,int* alloc_A_elems_num, MPI_Comm comm);
void broadcast_mat(double* partialA, double* partialC, double* transposedB, MPI_Comm comm);
int main2(int argc, char* argv[]);


int main(int argc, char *argv[]) {

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

    int i, j, k;

    // matrix sizes ,Aの行&Cの行、Bの列&Cの列、Aの列&Bの行
    int A_rows_num, B_cols_num, A_cols_num;

    double *A, *B, *loacl_A;
    double *transposed_B, *C, *local_C;

    double start, start2, end, end2;

    /* MPI Initialize */
    int myid, num_procs, local_A_elems_num, local_A_rows_num, local_C_elems_num;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    start = MPI_Wtime();

    /*rootプロセスの処理内容*/
    // TODO 次ここからstructを使って書き直す
    if (myid == 0) {
        A = read_mat_file(A_path, &A_rows_num, &A_cols_num);
        B = read_mat_file(B_path, &A_cols_num, &B_cols_num);
        show_mat(B);
//        transposedB = (double *) malloc(B_cols_num * A_cols_num * sizeof (double));
        transposed_B = transepose_mat(B);
        show_mat(transposed_B);
        if (A == NULL || B == NULL) {
            printf("ERROR filepath is wrong");
        }
        free(B);

        local_C_elems_num = A_rows_num * B_cols_num / num_procs;
        local_A_elems_num = A_rows_num * A_cols_num / num_procs;
        local_A_rows_num = A_rows_num / num_procs;
    }
    
    broadcast_mat_size(&A_rows_num, &B_cols_num, &A_cols_num, &local_C_elems_num, 
            &local_A_rows_num, &local_A_elems_num, MPI_COMM_WORLD);

    loacl_A = (double *) malloc(local_A_elems_num * sizeof (double));
    local_C = (double *) malloc(local_C_elems_num * sizeof (double));

//    broadcast_mat(partialA, partialA, transposedB, MPI_COMM_WORLD);
    if (myid == 0) {
        C = (double *) malloc(A_rows_num * B_cols_num * sizeof (double));

        //データを送信
        for (i = 1; i < num_procs; i++) {
            int a_first = i * local_A_elems_num;
            for (k = 0; k < local_A_elems_num; k++) {
                loacl_A[k] = A[a_first + k];
            }
            MPI_Send(loacl_A, local_C_elems_num, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }

        for (i = 0; i < local_A_elems_num; i++) {
            loacl_A[i] = A[i];
        }

        free(A);
        printf("calc start\n");
    } else {
        MPI_Status status;
        MPI_Recv(loacl_A, local_A_elems_num, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    }

    MPI_Bcast(transposed_B, B_cols_num * A_cols_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*各プロセスの処理内容*/
    double temp;
    start2 = MPI_Wtime();

    // 計算部分 tempをレジスタに置くこととアンローリングで高速化
#pragma omp parallel for private(i,k,temp)
    for (j = 0; j < B_cols_num; j += 8) {
        for (i = 0; i < local_A_rows_num; i++) {
            for (k = 0; k < A_cols_num; k++) {
                temp = loacl_A[i * A_cols_num + k];
                local_C[i * B_cols_num + j] += temp * transposed_B[j * A_cols_num + k];
                local_C[i * B_cols_num + j + 1] += temp * transposed_B[(j + 1) * A_cols_num + k];
                local_C[i * B_cols_num + j + 2] += temp * transposed_B[(j + 2) * A_cols_num + k];
                local_C[i * B_cols_num + j + 3] += temp * transposed_B[(j + 3) * A_cols_num + k];
                local_C[i * B_cols_num + j + 4] += temp * transposed_B[(j + 4) * A_cols_num + k];
                local_C[i * B_cols_num + j + 5] += temp * transposed_B[(j + 5) * A_cols_num + k];
                local_C[i * B_cols_num + j + 6] += temp * transposed_B[(j + 6) * A_cols_num + k];
                local_C[i * B_cols_num + j + 7] += temp * transposed_B[(j + 7) * A_cols_num + k];
            }
        }
    }

    MPI_Gather(local_C, local_C_elems_num, MPI_DOUBLE, C, local_C_elems_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /*各プロセスの処理内容終わり*/

    //rootに戻って出力
    if (myid == 0) {

        end2 = MPI_Wtime();

        show_mat(C);
        write_mat_file(C_path, C, A_rows_num, B_cols_num);
        end = MPI_Wtime();
        free(C);

        printf("N:%d, all_time:%f,calc_time:%f,processes:%d(alloc_c_num:%d).threads:%d\n", A_rows_num * B_cols_num, (end - start), (end2 - start2), num_procs, local_C_elems_num, omp_get_max_threads());
    }

    //終了処理
    printf("finalize");
    free(loacl_A);
    free(transposed_B);
    free(local_C);
    MPI_Finalize();

    return 0;
}

double* read_mat_file(const char* path, int* rows_num, int* cols_num) {
    FILE* fp;
    int row, col;
    double value;
    double *matrix;
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

    fp = fopen(path, "rb");
    if (fp != NULL) {
        fread(&row, sizeof (int), 1, fp);
        fread(&col, sizeof (int), 1, fp);
        printf("%d,%d\n", row, col);
        *rows_num = row;
        *cols_num = col;

        matrix = (double *) malloc(row * col * sizeof (double));

        for (i = 0; i < row; ++i) {
            for (j = 0; j < col; ++j) {
                fread(&value, sizeof (double), 1, fp);
                matrix[i * col + j] = value;
            }
        }
    } else {
        return NULL;
    }
    fclose(fp);

    return matrix;
#endif
}


void write_mat_file(const char* path, double* matrix, int rows_num, int cols_num) {
    FILE* fp;
    int i, j;
    double value;

    fp = fopen(path, "wb");
    if (fp != NULL) {

        fwrite(&rows_num, sizeof (int), 1, fp);
        fwrite(&cols_num, sizeof (int), 1, fp);
        for (i = 0; i < rows_num; ++i) {
            for (j = 0; j < cols_num; ++j) {
                value = matrix[i * cols_num + j];
                fwrite(&value, sizeof (double), 1, fp);
            }
        }
    }
    fclose(fp);
}

void broadcast_mat_size(int* A_rows_num,int* B_cols_num,int* A_cols_num,int* alloc_C_elems_num,int* alloc_A_rows_num,int* alloc_A_elems_num, MPI_Comm comm){
    MPI_Bcast(A_rows_num, 1, MPI_INT, 0, comm);
    MPI_Bcast(B_cols_num, 1, MPI_INT, 0, comm);
    MPI_Bcast(A_cols_num, 1, MPI_INT, 0, comm);
    MPI_Bcast(alloc_C_elems_num, 1, MPI_INT, 0, comm);
    MPI_Bcast(alloc_A_rows_num, 1, MPI_INT, 0, comm);
    MPI_Bcast(alloc_A_elems_num, 1, MPI_INT, 0, comm);    
}
