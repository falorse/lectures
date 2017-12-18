
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi/mpi.h>
#include <omp.h>
#include "matrix.h"

#define MAXLEN 1024
#define OUT 0
//#define NO_MAT_FILE

void get_file_paths(char* A_path, char* B_path, char* C_path, char** argv);
bool is_input_matrixes_invalid(struct Matrix* A, struct Matrix* B);
void create_global_matrixes(struct Matrix* A, char* A_path, struct Matrix* B,
        char* B_path, struct Matrix* transeposed_B, struct Matrix* C, int myid);
void create_local_matrixes(struct Matrix* local_A, struct Matrix* local_C,
        int a_rows, int a_cols, int b_cols, int myid, int procs_num);
void communicate_by_mpi(struct Matrix* A, struct Matrix* local_A,
        struct Matrix* B, struct Matrix* transeposed_B);
void calculate(const struct Matrix* local_A,
        const struct Matrix* transeposed_B, struct Matrix* local_C);
struct Matrix* read_input_file(char* file_path, int myid);

int main(int argc, char *argv[]) {

    // check argc and get file paths
    if (argc != 4) {
        printf("ERROR USAGE: ./mul_matrix left_mat_path right_mat_path ans_mat_path\n");
        return 0;
    }
    char A_path[MAXLEN], B_path[MAXLEN], C_path[MAXLEN];
    get_file_paths(A_path, B_path, C_path, argv);

    // MPI Initialize
    int myid, procs_num;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procs_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    double start, end;
    start = MPI_Wtime();

    // create matrixes
    // create A, B, C without bufs in myid != 0
    struct Matrix *A, *B, *C, *transeposed_B, *local_A, *local_C;
    A = read_input_file(A_path, myid);
    B = read_input_file(B_path, myid);
    if (myid == 0) {
        C = create_mat(A->rows, B->cols);
        transeposed_B = transepose_mat(B);
    } else {
        C = create_mat_without_bufs(A->rows, B->cols);
        transeposed_B = create_mat(B->cols, B->rows);
    }
    local_A = create_mat(A->rows / procs_num, A->cols);
    local_C = create_mat(A->rows / procs_num, B->cols);

    if(is_input_matrixes_invalid(A, B)){
        printf("ERROR: an input matrix file is invalid\n");
        return 0;
    }

    communicate_by_mpi(A, local_A, B, transeposed_B);

    calculate(local_A, transeposed_B, local_C);

    MPI_Gather(local_C->bufs, get_mat_size(local_C), MPI_DOUBLE,
            C->bufs, get_mat_size(local_C), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // write ans matrix to file and output run time
    if (myid == 0) {
        write_mat_file(C_path, C);
        end = MPI_Wtime();
        printf("N:%d, time:%f, processes:%d, threads:%d\n",
                get_mat_size(C), (end - start), procs_num, omp_get_max_threads());
    }

    // finalize
    destroy_mat(C);
    destroy_mat(local_A);
    destroy_mat(transeposed_B);
    destroy_mat(local_C);
    MPI_Finalize();
    return 0;
}

void get_file_paths(char* A_path, char* B_path, char* C_path, char** argv) {
    strcpy(A_path, argv[1]);
    strcpy(B_path, argv[2]);
    strcpy(C_path, argv[3]);
};

bool is_input_matrixes_invalid(struct Matrix* A, struct Matrix* B){
  if(A == NULL || B == NULL || A->cols != B->rows) {
      return true;
  }
  return false;
};

struct Matrix* read_input_file(char* file_path, int myid) {
    int i, j, rows, cols;

    FILE* fp = fopen(file_path, "rb");
    if (fp == NULL) {
        return NULL;
    } else {
        fread(&rows, sizeof (int), 1, fp);
        fread(&cols, sizeof (int), 1, fp);
    }

    struct Matrix* mat = malloc(sizeof (struct Matrix));

    if (myid == 0) {
        mat = create_mat(rows, cols);
        for (i = 0; i < rows; ++i) {
            for (j = 0; j < cols; ++j) {
                double value;
                fread(&value, sizeof (double), 1, fp);
                set_value(mat, value, i, j);
            }
        }
    } else {
        mat = create_mat_without_bufs(rows, cols);
    }
    fclose(fp);

    return mat;
};

void communicate_by_mpi(struct Matrix* A, struct Matrix* local_A,
        struct Matrix* B, struct Matrix* transeposed_B) {
    MPI_Scatter(A->bufs, get_mat_size(local_A), MPI_DOUBLE, local_A->bufs,
            get_mat_size(local_A), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    destroy_mat(A);

    MPI_Bcast(transeposed_B->bufs, get_mat_size(transeposed_B),
            MPI_DOUBLE, 0, MPI_COMM_WORLD);
    destroy_mat(B);
};

void create_global_matrixes(struct Matrix* A, char* A_path, struct Matrix* B,
        char* B_path, struct Matrix* transeposed_B, struct Matrix* C, int myid){
    A = read_input_file(A_path, myid);
    B = read_input_file(B_path, myid);
    if (myid == 0) {
        C = create_mat(A->rows, B->cols);
        transeposed_B = transepose_mat(B);
    } else {
        C = create_mat_without_bufs(A->rows, B->cols);
        transeposed_B = create_mat(B->cols, B->rows);
    }
}

void calculate(const struct Matrix* local_A,
        const struct Matrix* transeposed_B, struct Matrix* local_C) {
    int i, j, k;
    double tmp;

    // turning by loop unwinging and putting tmp on register
#pragma omp parallel for private(i,k,tmp)
    for (j = 0; j < transeposed_B->rows; j += 8) {
        for (i = 0; i < local_A->rows; i++) {
            for (k = 0; k < local_A->cols; k++) {
                tmp = get_value(local_A, i, k);
                add_value(local_C, tmp * get_value(transeposed_B, j, k), i, j);
                add_value(local_C, tmp * get_value(transeposed_B, j + 1, k), i, j + 1);
                add_value(local_C, tmp * get_value(transeposed_B, j + 2, k), i, j + 2);
                add_value(local_C, tmp * get_value(transeposed_B, j + 3, k), i, j + 3);
                add_value(local_C, tmp * get_value(transeposed_B, j + 4, k), i, j + 4);
                add_value(local_C, tmp * get_value(transeposed_B, j + 5, k), i, j + 5);
                add_value(local_C, tmp * get_value(transeposed_B, j + 6, k), i, j + 6);
                add_value(local_C, tmp * get_value(transeposed_B, j + 7, k), i, j + 7);
            }
        }
    }
};
