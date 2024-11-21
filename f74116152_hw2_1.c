#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

void readData(char *filename, int *t, int *w, int *h, int ***arr, int *D, int ***filter) {
    FILE *fp = fopen(filename, "r");
    fscanf(fp, "%d", t);
    fscanf(fp, "%d %d", h, w);
    *arr = (int **)malloc(sizeof(int *) * *h);
    for (int i = 0; i < *h; i++) {
        (*arr)[i] = (int *)malloc(sizeof(int) * *w);
        for (int j = 0; j < *w; j++) {
            fscanf(fp, "%d", &(*arr)[i][j]);
        }
    }
    fscanf(fp, "%d", D);
    *filter = (int **)malloc(sizeof(int *) * *D);
    for (int i = 0; i < *D; i++) {
        (*filter)[i] = (int *)malloc(sizeof(int) * *D);
        for (int j = 0; j < *D; j++) {
            fscanf(fp, "%d", &(*filter)[i][j]);
        }
    }
    fclose(fp);
}

// calculate dot product of a and b with start point (startX, startY) and size (w, h)
int dotProduct2(int **a, int startX, int startY, int aW, int aH, int **b, int w, int h) {
    int sum = 0;
    int tx, ty;
    ty = startY;
    if (ty < 0) ty += aH;
    for (int y = 0; y < h; ++y) {
        tx = startX;
        if (tx < 0) tx += aW;
        for (int x = 0; x < w; ++x) {
            sum += a[ty][tx] * b[y][x];
            ++tx;
            if (tx >= aW) tx -= aW;
        }
        ++ty;
        if (ty >= aH) ty -= aH;
    }
    return sum;
}

static inline void copyMatrix(int **src, int **dest, int h, int w) {
    for (int i = 0; i < h; i++) {
        memcpy(dest[i], src[i], sizeof(int) * w);
    }
}

static inline int** mallocMatrix(int h, int w) {
    int **arr = (int **)malloc(sizeof(int *) * h);
    for (int i = 0; i < h; i++) {
        arr[i] = (int *)malloc(sizeof(int) * w);
    }
    return arr;
}

void printAnswer(int **arr, int h, int w) {
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            printf("%d ", arr[i][j]);
        }
    }
}

int main(int argc, char **argv) {

    int rank, procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    int t, w, h, D;
    int **arr, **filter;

    if (rank == 0) {
        char filename[100];
        scanf("%s", filename);
        readData(filename, &t, &w, &h, &arr, &D, &filter);
    }

    MPI_Bcast(&t, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&w, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&h, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&D, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        arr = mallocMatrix(h, w);
        filter = mallocMatrix(D, D);
    }

    for (int i = 0; i < h; ++i) {
        MPI_Bcast(arr[i], w, MPI_INT, 0, MPI_COMM_WORLD);
    }
    for (int i = 0; i < D; ++i) {
        MPI_Bcast(filter[i], D, MPI_INT, 0, MPI_COMM_WORLD);
    }

    int dSquare = D * D;
    int start = -(D - 1) / 2;

    int **next = mallocMatrix(h, w);

    for (int i = 0; i < t; ++i) {
        // calculate A_t
        for (int y = rank; y < h; y += procs) {
            for (int x = 0; x < w; x++) {
                next[y][x] = dotProduct2(arr, x+start, y+start, w, h, filter, D, D) / dSquare;
            }
        }
        // broadcast result
        for (int y = 0; y < h; y++) {
            MPI_Bcast(next[y], w, MPI_INT, y % procs, MPI_COMM_WORLD); 
        }
        copyMatrix(next, arr, h, w);
    }
    
    if (rank == 0) {
        printAnswer(arr, h, w);
    }

    // 2.7-2.9 sec (time judge)
    // judge limit: 1.5 sec
    MPI_Finalize();
    return 0;
}