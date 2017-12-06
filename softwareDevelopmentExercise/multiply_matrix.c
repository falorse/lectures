
#include <stdio.h>
#include <stdlib.h>
#include <mpi/mpi.h>
#include <omp.h>
#include "matrix.h"

#define MAXLEN 1024
#define OUT 0
//#define NO_MAT_FILE

double* mulutipy_mat(const double* leftMat, const double* rightMat);
void broadcast_mat_size(int* A_rows_num,int* B_cols_num,int* A_cols_num,int* alloc_C_elems_num,int* alloc_A_rows_num,int* alloc_A_elems_num, MPI_Comm comm);
void broadcast_mat(double* partialA, double* partialC, double* transposedB, MPI_Comm comm);

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

    struct Matrix *A, *B, *local_A;
    struct Matrix *transposed_B, *C, *local_C;

    double start, start2, end, end2;

    /* MPI Initialize */
    int myid, num_procs, local_A_elems_num, local_A_rows_num, local_C_elems_num;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    start = MPI_Wtime();

    /*rootプロセスの処理内容*/
    if (myid == 0) {
        A = read_mat_file(A_path);
        A_rows_num = A->rows;
        A_cols_num = A->cols;
        B = read_mat_file(B_path);
        B_cols_num = B->cols;
        C = create_mat(A_rows_num, B_cols_num);
        
        transposed_B = transpose_mat(B);
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

    local_A = create_mat(local_A_rows_num, A_cols_num);
    local_C = create_mat(local_A_rows_num, A_cols_num);
    
    if (myid == 0) {

        double set_value;

        //データを送信
        for (i = 1; i < num_procs; i++) {
            int a_first = i * local_A_elems_num;
            for (k = 0; k < local_A_elems_num; k++) {
                set_value = get_mat_value_by_index(A, a_first + k);
                set_mat_value_by_index(local_A, set_value, k);
            }
            MPI_Send(local_A->bufs, local_A_elems_num, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }

        for (i = 0; i < local_A_elems_num; i++) {
            set_value = get_mat_value_by_index(local_A, i);
            set_mat_value_by_index(local_A, set_value, i);
        }
        

        free(A);
    } else {
        MPI_Status status;
        MPI_Recv(local_A->bufs, local_A_elems_num, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        transposed_B = create_mat(B_cols_num, A_cols_num);
    }
    
    MPI_Bcast(transposed_B->bufs, transposed_B->rows * transposed_B->cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*各プロセスの処理内容*/
    double temp;
    start2 = MPI_Wtime();

    // 計算部分 tempをレジスタに置くこととアンローリングで高速化
#pragma omp parallel for private(i,k,temp)
    for (j = 0; j < B_cols_num; j += 8) {
        for (i = 0; i < local_A_rows_num; i++) {
            for (k = 0; k < A_cols_num; k++) {
                temp = get_mat_value(local_A, i, k);
                set_mat_value(local_C, temp * get_mat_value(transposed_B, j, k), i, j);
                set_mat_value(local_C, temp * get_mat_value(transposed_B, j + 1, k), i, j + 1);
                set_mat_value(local_C, temp * get_mat_value(transposed_B, j + 2, k), i, j + 2);
                set_mat_value(local_C, temp * get_mat_value(transposed_B, j + 3, k), i, j + 3);
                set_mat_value(local_C, temp * get_mat_value(transposed_B, j + 4, k), i, j + 4);
                set_mat_value(local_C, temp * get_mat_value(transposed_B, j + 5, k), i, j + 5);
                set_mat_value(local_C, temp * get_mat_value(transposed_B, j + 6, k), i, j + 6);
                set_mat_value(local_C, temp * get_mat_value(transposed_B, j + 7, k), i, j + 7);            }
        }
    }
    
    MPI_Gather(local_C->bufs, local_C_elems_num, MPI_DOUBLE, get_mat_bufs(C), local_C_elems_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /*各プロセスの処理内容終わり*/

    printf("calc end\n");
    //rootに戻って出力
    if (myid == 0) {

        end2 = MPI_Wtime();

        show_mat(C);
        write_mat_file(C_path, C);
        end = MPI_Wtime();
        free(C);

        printf("N:%d, all_time:%f,calc_time:%f,processes:%d(alloc_c_num:%d).threads:%d\n", A_rows_num * B_cols_num, (end - start), (end2 - start2), num_procs, local_C_elems_num, omp_get_max_threads());
    }

    //終了処理
    printf("finalize");
    free(local_A);
    free(transposed_B);
    free(local_C);
    MPI_Finalize();

    return 0;
}


void broadcast_mat_size(int* A_rows_num,int* B_cols_num,int* A_cols_num,int* alloc_C_elems_num,int* alloc_A_rows_num,int* alloc_A_elems_num, MPI_Comm comm){
    MPI_Bcast(A_rows_num, 1, MPI_INT, 0, comm);
    MPI_Bcast(B_cols_num, 1, MPI_INT, 0, comm);
    MPI_Bcast(A_cols_num, 1, MPI_INT, 0, comm);
    MPI_Bcast(alloc_C_elems_num, 1, MPI_INT, 0, comm);
    MPI_Bcast(alloc_A_rows_num, 1, MPI_INT, 0, comm);
    MPI_Bcast(alloc_A_elems_num, 1, MPI_INT, 0, comm);    
}
