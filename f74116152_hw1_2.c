#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include <math.h>

typedef struct Point{
    int x, y;
    int ori_idx; // original index
} Point;

int compare(const void *a, const void *b) {
    const Point *p1 = (const Point *)a;
    const Point *p2 = (const Point *)b;
    if (p1->x == p2->x) {
        return p1->y - p2->y;
    }
    return p1->x - p2->x;
}

static inline bool equal(Point p1, Point p2) {
    return p1.x == p2.x && p1.y == p2.y;
}

static inline void pcpy(Point *p1, Point *p2) {
    p1->x = p2->x;
    p1->y = p2->y; 
}

// return > 0 if p1 to p2 from p0 is counterclockwise 
long long cross(Point p0, Point p1, Point p2) {
    return (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x);
}
void Andrew_monotone_chain(Point *p, int size, Point *resultHull, int *resultSize, int rank) {
    
    int result = 0;
    
    if (rank == 0) {
        // upper hull
        for (int i = size - 1; i >= 0; --i) {
            while(result > 1 && cross(resultHull[result-2], resultHull[result-1], p[i]) <= 0) {
                --result;
            }
            resultHull[result] = p[i];
            ++result;
        }
    } else if (rank == 1) {
            // lower hull
        for (size_t i = 0; i < size; ++i) {
            while(result > 1 && cross(resultHull[result-2], resultHull[result-1], p[i]) <= 0) {
                --result;
            }
            resultHull[result] = p[i];
            ++result;
        }
    } else {
        *resultSize = 0;
        return;
    }
    *resultSize = result;
}

int main(int argc, char *argv[]) {

    int rank, size;

    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Datatype MPI_POINT;
    MPI_Type_contiguous(3, MPI_INT, &MPI_POINT);
    MPI_Type_commit(&MPI_POINT);

    int n;
    Point *p = NULL;

    if (rank == 0) {
        char filename[50];
        scanf("%s", filename);
        FILE *fp = fopen(filename, "r");
        fscanf(fp, "%d", &n);
        p = malloc(n * sizeof(Point));
        for (int i = 0; i < n; ++i) {
            fscanf(fp, "%d %d", &p[i].x, &p[i].y);
            p[i].ori_idx = i + 1;
        }
        qsort(p, n, sizeof(Point), compare);
        fclose(fp);
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        p = malloc(n * sizeof(Point));
    }

    MPI_Bcast(p, n, MPI_POINT, 0, MPI_COMM_WORLD);

    Point *resultHull = malloc(n * sizeof(Point));
    int resultSize = 0;
    Andrew_monotone_chain(p, n, resultHull, &resultSize, rank);
    
    int receivedSize;
    if (rank == 0) {
        // print upper hull
        for (int i = resultSize - 1; i >= 0; --i) {
            printf("%d ", resultHull[i].ori_idx);
        }
        for (int _rank = 1; _rank < size; ++_rank) {
            MPI_Recv(&receivedSize, 1, MPI_INT, _rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            Point *receivedHull = malloc(receivedSize * sizeof(Point));
            MPI_Recv(receivedHull, receivedSize, MPI_POINT, _rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // print lower hull
            for (int i = receivedSize - 2; i > 0; --i) {
                printf("%d ", receivedHull[i].ori_idx);
            }
        }
    } else {
        MPI_Send(&resultSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(resultHull, resultSize, MPI_POINT, 0, 0, MPI_COMM_WORLD);
    }

    free(p);
    free(resultHull);
    MPI_Finalize();
    return 0;

}