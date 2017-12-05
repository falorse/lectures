
#include <stdio.h>
#include <stdlib.h>
#include </usr/lib/openmpi/include/mpi.h>
#include <omp.h>
#include "matrix.h"

#define MAXLEN 1024
#define OUT 0
//#define NO_MAT_FILE

double* read_mat_file(const char* path, int* rows_num, int* cols_num);
//double** readMatrix(const char* path, int* rows_num, int* cols_num);

void write_mat_file(const char* path, double* matrix, int rows_num, int cols_num);
double* mulutipy_mat(const double* leftMat, const double* rightMat);
//double* transpose_mat(const double* originalMat, int original_rows, int original_cols);
void broadcastMatSize(int* A_rows_num,int* B_cols_num,int* A_cols_num,int* alloc_C_elems_num,int* alloc_A_rows_num,int* alloc_A_elems_num, MPI_Comm comm);
void broadcastMat(double* partialA, double* partialC, double* transposedB, MPI_Comm comm);
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

    double *A, *B, *partialA;
    double *transposedB, *C, *partialC;

    double start, start2, end, end2;

    /* MPI Initialize */
    int rank, num_procs, alloc_A_elems_num, alloc_A_rows_num, alloc_C_elems_num;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    start = MPI_Wtime();

    /*rootプロセスの処理内容*/
    if (rank == 0) {
        A = read_mat_file(A_path, &A_rows_num, &A_cols_num);
        B = read_mat_file(B_path, &A_cols_num, &B_cols_num);
        show_mat(B);
//        transposedB = (double *) malloc(B_cols_num * A_cols_num * sizeof (double));
        transposedB = transepose_mat(B);
        show_mat(transposedB);
        if (A == NULL || B == NULL) {
            printf("ERROR filepath is wrong");
        }
        free(B);

        alloc_C_elems_num = A_rows_num * B_cols_num / num_procs;
        alloc_A_elems_num = A_rows_num * A_cols_num / num_procs;
        alloc_A_rows_num = A_rows_num / num_procs;
    }
    
    broadcastMatSize(&A_rows_num, &B_cols_num, &A_cols_num, &alloc_C_elems_num, 
            &alloc_A_rows_num, &alloc_A_elems_num, MPI_COMM_WORLD);

    partialA = (double *) malloc(alloc_A_elems_num * sizeof (double));
    partialC = (double *) malloc(alloc_C_elems_num * sizeof (double));

//    broadcastMat(partialA, partialA, transposedB, MPI_COMM_WORLD);
//    if (rank == 0) {
//        C = (double *) malloc(A_rows_num * B_cols_num * sizeof (double));
//
//        //データを送信
//        for (i = 1; i < num_procs; i++) {
//            int a_first = i * alloc_A_elems_num;
//            for (k = 0; k < alloc_A_elems_num; k++) {
//                partialA[k] = A[a_first + k];
//            }
//            MPI_Send(partialA, alloc_C_elems_num, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
//        }
//
//        for (i = 0; i < alloc_A_elems_num; i++) {
//            partialA[i] = A[i];
//        }
//
//        free(A);
//        printf("calc start\n");
//    } else {
//        MPI_Status status;
//        MPI_Recv(partialA, alloc_A_elems_num, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
//    }
//
//    MPI_Bcast(transposedB, B_cols_num * A_cols_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*各プロセスの処理内容*/
    double temp;
    start2 = MPI_Wtime();

    // 計算部分 tempをレジスタに置くこととアンローリングで高速化
#pragma omp parallel for private(i,k,temp)
    for (j = 0; j < B_cols_num; j += 8) {
        for (i = 0; i < alloc_A_rows_num; i++) {
            for (k = 0; k < A_cols_num; k++) {
                temp = partialA[i * A_cols_num + k];
                partialC[i * B_cols_num + j] += temp * transposedB[j * A_cols_num + k];
                partialC[i * B_cols_num + j + 1] += temp * transposedB[(j + 1) * A_cols_num + k];
                partialC[i * B_cols_num + j + 2] += temp * transposedB[(j + 2) * A_cols_num + k];
                partialC[i * B_cols_num + j + 3] += temp * transposedB[(j + 3) * A_cols_num + k];
                partialC[i * B_cols_num + j + 4] += temp * transposedB[(j + 4) * A_cols_num + k];
                partialC[i * B_cols_num + j + 5] += temp * transposedB[(j + 5) * A_cols_num + k];
                partialC[i * B_cols_num + j + 6] += temp * transposedB[(j + 6) * A_cols_num + k];
                partialC[i * B_cols_num + j + 7] += temp * transposedB[(j + 7) * A_cols_num + k];
            }
        }
    }

    MPI_Gather(partialC, alloc_C_elems_num, MPI_DOUBLE, C, alloc_C_elems_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /*各プロセスの処理内容終わり*/

    //rootに戻って出力
    if (rank == 0) {

        end2 = MPI_Wtime();

        show_mat(C);
        write_mat_file(C_path, C, A_rows_num, B_cols_num);
        end = MPI_Wtime();
        free(C);

        printf("N:%d, all_time:%f,calc_time:%f,processes:%d(alloc_c_num:%d).threads:%d\n", A_rows_num * B_cols_num, (end - start), (end2 - start2), num_procs, alloc_C_elems_num, omp_get_max_threads());
    }

    //終了処理
    printf("finalize");
    free(partialA);
    free(transposedB);
    free(partialC);
    MPI_Finalize();

    return 0;
}

//double* read_mat_file(const char* path, int* rows_num, int* cols_num) {
//    FILE* fp;
//    int row, col;
//    double value;
//    double *matrix;
//    int i, j;
//
//#ifdef NO_MAT_FILE
//    int mat_rows = 2;
//    matrix = (double *) malloc(mat_rows * mat_rows * sizeof (double));
//    for (i = 0; i < mat_rows * mat_rows; i++) {
//        matrix[i] = i;
//    }
//    *rows_num = mat_rows;
//    *cols_num = mat_rows;
//    return matrix;
//#else
//
//    fp = fopen(path, "rb");
//    if (fp != NULL) {
//        fread(&row, sizeof (int), 1, fp);
//        fread(&col, sizeof (int), 1, fp);
//        printf("%d,%d\n", row, col);
//        *rows_num = row;
//        *cols_num = col;
//
//        matrix = (double *) malloc(row * col * sizeof (double));
//
//        for (i = 0; i < row; ++i) {
//            for (j = 0; j < col; ++j) {
//                fread(&value, sizeof (double), 1, fp);
//                matrix[i * col + j] = value;
//            }
//        }
//    } else {
//        return NULL;
//    }
//    fclose(fp);
//
//    return matrix;
//#endif
//}

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

//double* transpose_mat(const double* originalMat, int original_rows, int original_cols) {
//    double* transposedMat = (double *) malloc(original_rows * original_cols * sizeof (double));
//
//    int i,j;
//    for (i = 0; i < original_rows; i++) {
//        for (j = 0; j < original_cols; j++) {
//            transposedMat[i * original_rows + j] = originalMat[j * original_cols + i];
//        }
//    }
//    
//    return transposedMat;
//}


void broadcastMatSize(int* A_rows_num,int* B_cols_num,int* A_cols_num,int* alloc_C_elems_num,int* alloc_A_rows_num,int* alloc_A_elems_num, MPI_Comm comm){
    MPI_Bcast(A_rows_num, 1, MPI_INT, 0, comm);
    MPI_Bcast(B_cols_num, 1, MPI_INT, 0, comm);
    MPI_Bcast(A_cols_num, 1, MPI_INT, 0, comm);
    MPI_Bcast(alloc_C_elems_num, 1, MPI_INT, 0, comm);
    MPI_Bcast(alloc_A_rows_num, 1, MPI_INT, 0, comm);
    MPI_Bcast(alloc_A_elems_num, 1, MPI_INT, 0, comm);    
}

//void broadcastMat(double* partialA, double* partialA, double* transposedB, MPI_Comm comm ){
//        if (rank == 0) {
//        C = (double *) malloc(A_rows_num * B_cols_num * sizeof (double));
//
//        //データを送信
//        for (i = 1; i < num_procs; i++) {
//            int a_first = i * alloc_A_elems_num;
//            for (k = 0; k < alloc_A_elems_num; k++) {
//                partialA[k] = A[a_first + k];
//            }
//            MPI_Send(partialA, alloc_C_elems_num, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
//        }
//
//        for (i = 0; i < alloc_A_elems_num; i++) {
//            partialA[i] = A[i];
//        }
//
//        free(A);
//        printf("calc start\n");
//    } else {
//        MPI_Status status;
//        MPI_Recv(partialA, alloc_A_elems_num, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
//    }
//
//    MPI_Bcast(transposedB, B_cols_num * A_cols_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//}
