
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

#define MAXLEN 1024
#define OUT 0

double* readMatFile(const char* path, int* rows_num, int* cols_num);
//double** readMatrix(const char* path, int* rows_num, int* cols_num);

void writeMat(const char* path, double* matrix,int rows_num, int cols_num);
double* mulutipyMat(const double* leftMat, const double* rightMat);
double* transposeMat(const double* originalMat);
void distributeRightMat();

int main2(int argc, char* argv[]);

int main(int argc, char *argv[]) {
    
    // 行列A,Bに対して、A * B = C の計算を行い、Cをファイルに出力する
    int i, j, k;
	
    // matrix sizes ,Aの行&Cの行、Bの列&Cの列、Aの列&Bの行
    int A_rows, B_cols, A_cols;

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
        // MatFileの読み込みと行列の行数が正しいかどうか確認
        A = readMatFile(A_path, &A_rows, &A_cols);
        int tmp = A_cols;

        B = readMatFile(B_path, &A_cols, &B_cols);
        transposedB = transposeMat(B);
        
        if (tmp != A_cols) {
            printf("ERROR A cols_num != B rows_num");
            return 0;
        }

        alloc_C_elems_num = A_rows * B_cols / num_procs;
        alloc_A_elems_num = A_rows * A_cols / num_procs;
        alloc_A_rows_num = A_rows / num_procs;
    }

    MPI_Bcast(&A_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&B_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&A_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alloc_C_elems_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alloc_A_rows_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alloc_A_elems_num, 1, MPI_INT, 0, MPI_COMM_WORLD);

    partialA = (double *) malloc(alloc_A_elems_num * sizeof (double));
    partialC = (double *) malloc(alloc_C_elems_num * sizeof (double));

    transposedB = (double *) malloc(B_cols * A_cols * sizeof (double));

    if (rank == 0) {
        C = (double *) malloc(A_rows * B_cols * sizeof (double));
        transposedB = (double *) malloc(B_cols * A_cols * sizeof (double));

        //bを転置してb2に格納
        for (k = 0; k < A_cols; k++) {
            for (j = 0; j < B_cols; j++) {
                transposedB[k * B_cols + j] = B[j * A_cols + k];
            }
        }

        //データを送信
        for (i = 1; i < num_procs; i++) {
            int a_first = i*alloc_A_elems_num;
            for (k = 0; k < alloc_A_elems_num; k++) {
                partialA[k] = A[a_first + k];
            }
            MPI_Send(partialA, alloc_C_elems_num, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }

        for (i = 0; i < alloc_A_elems_num; i++) {
            partialA[i] = A[i];
        }

        free(A);
        free(B);
        printf("calc start\n");
    } else {
        MPI_Status status;
        MPI_Recv(partialA, alloc_A_elems_num, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    }

    MPI_Bcast(transposedB, B_cols * A_cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*各プロセスの処理内容*/
    double temp;
    start2 = MPI_Wtime();

    // 計算部分 tempをレジスタに置くこととアンローリングで高速化
#pragma omp parallel for private(i,k,temp)
    for (j = 0; j < B_cols; j += 8) {
        for (i = 0; i < alloc_A_rows_num; i++) {
            for (k = 0; k < A_cols; k++) {
                temp = partialA[i * A_cols + k];
                partialC[i * B_cols + j] += temp * transposedB[j * A_cols + k];
                partialC[i * B_cols + j + 1] += temp * transposedB[(j + 1) * A_cols + k];
                partialC[i * B_cols + j + 2] += temp * transposedB[(j + 2) * A_cols + k];
                partialC[i * B_cols + j + 3] += temp * transposedB[(j + 3) * A_cols + k];
                partialC[i * B_cols + j + 4] += temp * transposedB[(j + 4) * A_cols + k];
                partialC[i * B_cols + j + 5] += temp * transposedB[(j + 5) * A_cols + k];
                partialC[i * B_cols + j + 6] += temp * transposedB[(j + 6) * A_cols + k];
                partialC[i * B_cols + j + 7] += temp * transposedB[(j + 7) * A_cols + k];
            }
        }
    }

    MPI_Gather(partialC, alloc_C_elems_num, MPI_DOUBLE, C, alloc_C_elems_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*各プロセスの処理内容終わり*/

    //rootに戻って出力
    if (rank == 0) {

        end2 = MPI_Wtime();

        writeMat(C_path, C, A_rows, B_cols);
        end = MPI_Wtime();
        free(C);

        printf("N:%d, all_time:%f,calc_time:%f,processes:%d(alloc_c_num:%d).threads:%d\n", A_rows*B_cols, (end - start), (end2 - start2), num_procs, alloc_C_elems_num, omp_get_max_threads());
    }

    //終了処理
    free(partialA);
    free(transposedB);
    free(partialC);
    MPI_Finalize();

    return 0;
}

double* readMatFile(const char* path, int* rows_num, int* cols_num) {
    FILE* fp;
    int row, col;
    double value;
    double *matrix;

    int i, j;
    printf("read matrix %s\n", path);
    fp = fopen(path, "rb");
    if (fp != NULL) {
        fread(&row, sizeof (int), 1, fp);
        fread(&col, sizeof (int), 1, fp);
        *rows_num = row;
        *cols_num = col;

        matrix = (double *) malloc(row * col * sizeof (double));

        for (i = 0; i < row; ++i) {
            for (j = 0; j < col; ++j) {
                fread(&value, sizeof (double), 1, fp);
                matrix[i * col + j] = value;
            }
        }
    }
    fclose(fp);

    return matrix;
}

void writeMat(const char* path, double* matrix, int rows_num, int cols_num) {
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

