
#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>


void dfs(uint32_t *test, int *testCounts, int *testCosts, int cur, bool *curTests, int *result, int m, int n) {
    if (cur >= m) {
        uint32_t testBuffer = 0;

        for (int i = 0; i < m; ++i) {
            if (curTests[i]) {
                testBuffer |= test[i];
            }
        }
        if (testBuffer == (1 << n) - 1) {
            *result += 1;
        }
        return;
    }
    dfs(test, testCounts, testCosts, cur+1, curTests, result, m, n);
    curTests[cur] = true;
    dfs(test, testCounts, testCosts, cur+1, curTests, result, m, n);
    curTests[cur] = false;
}

// generate the initial curTests array, separate to 64 parts
void getStartTest(int index, bool *curTests) {
    curTests[5] = (index & 1);
    curTests[4] = (index & 2) >> 1;
    curTests[3] = (index & 4) >> 2;
    curTests[2] = (index & 8) >> 3;
    curTests[1] = (index & 16) >> 4;
    curTests[0] = (index & 32) >> 5;
}

int bruteForceForMPI(int processId, int totalProcess, uint32_t *test, int *testCounts, int *testCosts, int m, int n) {
    bool *curTests = (bool*)calloc(m, sizeof(bool));
    int result = 0;

    if (totalProcess == 1 || m < 7) {
        if (processId != 0) {
            free(curTests);
            return 0;
        }
        dfs(test, testCounts, testCosts, 0, curTests, &result, m, n);
        free(curTests);
        return result;
    }

    int idx = processId;

    while (idx < 64) {
        getStartTest(idx, curTests);
        dfs(test, testCounts, testCosts, 6, curTests, &result, m, n);
        idx += totalProcess;
    }
    free(curTests);
    return result;
}

int readFile(int *mPtr, int *nPtr, uint32_t **testPtr, int **testCountsPtr, int **testCostsPtr) {
    char filename[50];
    scanf("%s", filename);

    FILE *inFile = fopen(filename, "r");
    if (inFile == NULL) {
        printf("Could not open the file %s\n", filename);
        return 1;
    }

    int m, n;
    fscanf(inFile, "%d %d", &n, &m);
    *mPtr = m;
    *nPtr = n;

    uint32_t *test = (uint32_t*)calloc(m, sizeof(uint32_t));
    int *testCounts = (int*)malloc(m * sizeof(int));
    int *testCosts = (int*)malloc(m * sizeof(int));

    *testPtr = test;
    *testCountsPtr = testCounts;
    *testCostsPtr = testCosts;

    int testCount, testCost;

    for(int i = 0; i < m; ++i) {
        fscanf(inFile, "%d %d", &testCount, &testCost);
        testCounts[i] = testCount;
        testCosts[i] = testCost;
        int part;
        for (int j = 0; j < testCount; ++j) {
            fscanf(inFile, "%d", &part);
            // test[i][j] = part - 1; // Note: change to 0-based index
            test[i] |= (1 << (part - 1));
        }
    }

    fclose(inFile);
    return 0;
}


int main(int argc,char *argv[]) {

    int rank, numprocs;

    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // initialize and read file
    int m, n;
    int *testCounts = NULL, *testCosts = NULL;
    uint32_t *test = NULL;
    int result = 0;

    if (rank == 0) {
        if (readFile(&m, &n, &test, &testCounts, &testCosts)) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // broadcast
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank != 0) {
        testCounts = (int*)malloc(m * sizeof(int));
        testCosts = (int*)malloc(m * sizeof(int));
    }

    MPI_Bcast(testCounts, m, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(testCosts, m, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        test = (uint32_t*)malloc(m * sizeof(uint32_t));
    }

    MPI_Bcast(test, m, MPI_INT, 0, MPI_COMM_WORLD);

    int localResult = bruteForceForMPI(rank, numprocs, test, testCounts, testCosts, m, n);

    // sum up the result and print
    MPI_Reduce(&localResult, &result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("%d", result);
    }

    free(test);
    free(testCounts);
    free(testCosts);
    
    MPI_Finalize();

    return 0;

}

