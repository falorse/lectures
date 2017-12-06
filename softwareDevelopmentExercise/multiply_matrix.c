
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi/mpi.h>
#include <omp.h>
#include "matrix.h"

#define MAXLEN 1024
#define OUT 0
//#define NO_MAT_FILE

double* mulutipy_mat(const double* leftMat, const double* rightMat);
void broadcast_mat_sizes(int* A_rows, int* B_cols, int* A_cols_num,
        int* local_C_size, int* local_A_rows, int* local_A_size, MPI_Comm comm);
void broadcast_mat(double* partialA, double* partialC, double* transposedB, MPI_Comm comm);
struct Matrix* read_input_file(char* file_path, int myid);

int main(int argc, char *argv[]) {

    // 引数チェックとファイルパスの取得
    char A_path[MAXLEN], B_path[MAXLEN], C_path[MAXLEN];
    if (argc == 4) {
        strcpy(A_path, argv[1]);
        strcpy(B_path, argv[2]);
        strcpy(C_path, argv[3]);
    } else {
        printf("ERROR USAGE: mul_matrix leftMat_path rightMat_path ansMat_path\n");
        return 0;
    }

    int i, j, k;
    struct Matrix *A, *B, *local_A, *local_C, *transeposed_B, *C;
    double start, start2, end, end2;

    /* MPI Initialize */
    int myid, procs_num;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procs_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    start = MPI_Wtime();

    // make matrixes
    // in id != 0, create A, B, C without bufs
    A = read_input_file(A_path, myid);
    B = read_input_file(B_path, myid);
    local_A = create_mat(A->rows / procs_num, A->cols);
    local_C = create_mat(A->rows / procs_num, B->cols);
    if (myid == 0) {
        C = create_mat(A->rows, B->cols);
        transeposed_B = transepose_mat(B);
    } else {
        C = create_mat_without_bufs(A->rows, B->cols);
        transeposed_B = create_mat(B->cols, B->rows);
    }

    destroy_mat(B);

    printf("A: %d, B: %d, t_B: %d, C: %d, l_A: %d, l_C: %d id: %d\n",
            get_mat_size(A), get_mat_size(B), get_mat_size(transeposed_B),
            get_mat_size(C), get_mat_size(local_A), get_mat_size(local_C), myid);

    if (myid == 0) {

        double value;

        //データを送信
        for (i = 1; i < procs_num; i++) {
            int local_A_size = get_mat_size(local_A);
            int local_A_first_index = i * local_A_size;
            for (k = 0; k < local_A_size; k++) {
                value = get_mat_value_by_index(A, local_A_first_index + k);
                set_mat_value_by_index(local_A, value, k);
            }
            MPI_Send(local_A->bufs, local_A_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }

        for (i = 0; i < get_mat_size(local_A); i++) {
            value = get_mat_value_by_index(A, i);
            set_mat_value_by_index(local_A, value, i);
        }

    } else {
        MPI_Recv(local_A->bufs, get_mat_size(local_A), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, 0);
        //        transeposed_B = create_mat(B->cols, B->rows);
    }
    destroy_mat(A);

    MPI_Bcast(transeposed_B->bufs, get_mat_size(transeposed_B),
            MPI_DOUBLE, 0, MPI_COMM_WORLD);

    printf("A: %f, B: %f, t_B: %f, C: %f, l_A: %f, l_C: %f id: %d\n",
            A->bufs[1], B->bufs[1], transeposed_B->bufs[1],
            C->bufs[1], local_A->bufs[1], local_C->bufs[1], myid);

    printf("calc start\n");
    /*各プロセスの処理内容*/
    double temp;
    start2 = MPI_Wtime();

    // 計算部分 tempをレジスタに置くこととアンローリングで高速化
#pragma omp parallel for private(i,k,temp)
    for (j = 0; j < transeposed_B->rows; j += 8) {
        for (i = 0; i < local_A->rows; i++) {
            for (k = 0; k < local_A->cols; k++) {
                temp = get_mat_value(local_A, i, k);
                add_mat_value(local_C, temp * get_mat_value(transeposed_B, j, k), i, j);
                add_mat_value(local_C, temp * get_mat_value(transeposed_B, j + 1, k), i, j + 1);
                add_mat_value(local_C, temp * get_mat_value(transeposed_B, j + 2, k), i, j + 2);
                add_mat_value(local_C, temp * get_mat_value(transeposed_B, j + 3, k), i, j + 3);
                add_mat_value(local_C, temp * get_mat_value(transeposed_B, j + 4, k), i, j + 4);
                add_mat_value(local_C, temp * get_mat_value(transeposed_B, j + 5, k), i, j + 5);
                add_mat_value(local_C, temp * get_mat_value(transeposed_B, j + 6, k), i, j + 6);
                add_mat_value(local_C, temp * get_mat_value(transeposed_B, j + 7, k), i, j + 7);
            }
        }
    }

    printf("A: %f, B: %f, t_B: %f, C: %f, l_A: %f, l_C: %f id: %d\n",
            A->bufs[1], B->bufs[1], transeposed_B->bufs[1],
            C->bufs[1], local_A->bufs[1], local_C->bufs[1], myid);
    //        show_mat(C);

    MPI_Gather(local_C->bufs, get_mat_size(local_C), MPI_DOUBLE, get_mat_bufs(C), get_mat_size(local_C), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /*各プロセスの処理内容終わり*/
    //    printf("%f,%f,%f,%f\n", local_C->bufs[1], local_C->bufs[3], local_C->bufs[4], local_C->bufs[10]);

    printf("calc end\n");
    //rootに戻って出力
    if (myid == 0) {

        end2 = MPI_Wtime();

        show_mat(C);
        write_mat_file(C_path, C);
        end = MPI_Wtime();

        printf("N:%d, all_time:%f,calc_time:%f,processes:%d(alloc_c_num:%d).threads:%d\n", get_mat_size(C), (end - start), (end2 - start2), procs_num, get_mat_size(local_C), omp_get_max_threads());
        free(C);
    }

    //終了処理
    printf("finalize");
    free(local_A);
    free(transeposed_B);
    free(local_C);
    MPI_Finalize();

    return 0;
}

struct Matrix* read_input_file(char* file_path, int myid) {
    FILE* fp;
    int i, j, rows, cols;

    fp = fopen(file_path, "rb");
    if (fp == NULL) {
        return NULL;
    } else {
        fread(&rows, sizeof (int), 1, fp);
        fread(&cols, sizeof (int), 1, fp);
    }

    struct Matrix* mat = NULL;
    mat = malloc(sizeof (struct Matrix));

    if (myid == 0) {
        mat = create_mat(rows, cols);
        for (i = 0; i < rows; ++i) {
            for (j = 0; j < cols; ++j) {
                double value;
                fread(&value, sizeof (double), 1, fp);
                set_mat_value(mat, value, i, j);
            }
        }
    } else {
        mat = create_mat_without_bufs(rows, cols);
    }
    fclose(fp);

    return mat;
};
