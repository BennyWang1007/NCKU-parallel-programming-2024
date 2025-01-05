#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <math.h>

typedef struct Point{
    int x, y;
} Point;

// comapre x then y, return true if p1 smaller
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

void Andrew_monotone_chain(Point *p, size_t size, Point *resultHull, int *resultSize) {
    int upper = 0;
    int lower = 0;
    Point *upperHull = malloc((size + 1) * sizeof(Point));
    Point *lowerHull = malloc((size + 1) * sizeof(Point));

    // upper hull
    for (int i = size - 1; i >= 0; --i) {
        while(upper > 1 && cross(upperHull[upper-2], upperHull[upper-1], p[i]) <= 0) {
            --upper;
        }
        upperHull[upper] = p[i];
        ++upper;
    }

    // lower hull
    for (int i = 0; i < size; ++i) {
        while(lower > 1 && cross(lowerHull[lower-2], lowerHull[lower-1], p[i]) <= 0) {
            --lower;
        }
        lowerHull[lower] = p[i];
        ++lower;
    }

    // conbine
    int result = 0;
    for (int i = upper - 1; i >= 0; --i) {
        resultHull[result++] = upperHull[i];
    }
    for (int i = lower - 2; i > 0; --i) { // discard start and end
        resultHull[result++] = lowerHull[i];
    }
    *resultSize = result;
    free(upperHull);
    free(lowerHull);
}

void readData(char *filename, Point **points, int *n) {
    FILE *fp = fopen(filename, "r");
    fscanf(fp, "%d", n);
    *points = malloc(*n * sizeof(Point));
    for (int i = 0; i < *n; ++i) {
        fscanf(fp, "%d %d", &(*points)[i].x, &(*points)[i].y);
    }
    fclose(fp);
}

// minimum spanning tree
typedef struct MST_Node {
    int parent;
    int rank;
} MST_Node;

typedef struct Edge {
    int u, v;
    double w;
} Edge;

typedef struct MST {
    MST_Node root;
} MST;

int find(struct MST_Node *nodes, int i) {
    if (nodes[i].parent != i) {
        nodes[i].parent = find(nodes, nodes[i].parent);
    }
    return nodes[i].parent;
}

void Union(struct MST_Node *nodes, int x, int y) {
    int xroot = find(nodes, x);
    int yroot = find(nodes, y);

    if (nodes[xroot].rank < nodes[yroot].rank) {
        nodes[xroot].parent = yroot;
    } else if (nodes[xroot].rank > nodes[yroot].rank) {
        nodes[yroot].parent = xroot;
    } else {
        nodes[yroot].parent = xroot;
        nodes[xroot].rank++;
    }
}

int compareEdge(const void *a, const void *b) {
    const Edge *e1 = (const Edge *)a;
    const Edge *e2 = (const Edge *)b;
    if (e1->w < e2->w) {
        return -1;
    } else if (e1->w > e2->w) {
        return 1;
    } else {
        return 0;
    }
}

Edge* getEdges(Point *points, int n, int *edgeSizePtr) {
    Edge *edges = malloc(n * n * sizeof(Edge));
    int edgeSize = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (i == j) continue;
            edges[edgeSize].u = i;
            edges[edgeSize].v = j;
            double w = sqrt((points[i].x - points[j].x) * (points[i].x - points[j].x) + (points[i].y - points[j].y) * (points[i].y - points[j].y));
            w = floor(w * 10000) / 10000;
            edges[edgeSize].w = w;
            ++edgeSize;
        }
    }
    *edgeSizePtr = edgeSize;
    return edges;
}

void prim(Edge *edges, int edgeCount, Edge **resultPtr, int *resultSizePtr) {
    qsort(edges, edgeCount, sizeof(Edge), compareEdge);
    struct MST_Node *nodes = malloc(edgeCount * sizeof(struct MST_Node));
    for (int i = 0; i < edgeCount; ++i) {
        nodes[i].parent = i;
        nodes[i].rank = 0;
    }
    Edge *result = malloc(edgeCount * sizeof(Edge));
    int resultSize = 0;
    for (int i = 0; i < edgeCount; ++i) {
        int x = find(nodes, edges[i].u);
        int y = find(nodes, edges[i].v);
        if (x != y) {
            result[resultSize++] = edges[i];
            Union(nodes, x, y);
        }
    }
    *resultSizePtr = resultSize;
    *resultPtr = result;
    free(nodes);

}

double prim_mod(Edge *edges, int edgeCount, MST_Node *nodes, int pointSize, Edge *result, int *resultSizePtr, bool *inMST) {
    for (int i = 0; i < pointSize; ++i) {
        nodes[i].parent = i;
        nodes[i].rank = 0;
    }
    int resultSize = 0;
    for (int i = 0; i < edgeCount; ++i) {
        if (!inMST[edges[i].u] || !inMST[edges[i].v]) {
            continue;
        }
        int x = find(nodes, edges[i].u);
        int y = find(nodes, edges[i].v);
        if (x != y) {
            result[resultSize++] = edges[i];
            Union(nodes, x, y);
        }
    }
    *resultSizePtr = resultSize;
    double sum = 0;
    for (int i = 0; i < resultSize; ++i) {
        sum += result[i].w;
    }
    return sum;

}

double bruteForce(bool *inMST, MST_Node *nodes, const int pointSize, const int curIdx, const int *hullIdx, const int hullSize, Edge *edges, const int edgeSize, Edge *resultMST, int *MSTSize, double min) {
    if (curIdx == pointSize) {
        for (int i = 0; i < hullSize; ++i) {
            if (!inMST[hullIdx[i]]) {
                return min;
            }
        }
        double sum = prim_mod(edges, edgeSize, nodes, pointSize, resultMST, MSTSize, inMST);
        if (sum < min) return sum;
        return min;
    }
    inMST[curIdx] = true;
    min = bruteForce(inMST, nodes, pointSize, curIdx + 1, hullIdx, hullSize, edges, edgeSize, resultMST, MSTSize, min);
    for (int i = 0; i < hullSize; ++i) {
        if(hullIdx[i] == curIdx) return min;
    }
    inMST[curIdx] = false;
    min = bruteForce(inMST, nodes, pointSize, curIdx + 1, hullIdx, hullSize, edges, edgeSize, resultMST, MSTSize, min);
    return min;
}


int main(int argc, char **argv) {
    int thread_num = atoi(argv[1]);

    char filename[100];
    scanf("%s", filename);

    Point *points;
    int n;
    readData(filename, &points, &n);
    qsort(points, n, sizeof(Point), compare);

    // find convex hull
    Point *hull = malloc(n * sizeof(Point));
    int hullSize;
    Andrew_monotone_chain(points, n, hull, &hullSize);
    qsort(hull, hullSize, sizeof(Point), compare);

    // find index of hull points in points
    int *hullIdx = malloc(hullSize * sizeof(int));
    {
        int curHullIdx = 0;
        for (int i = 0; i < n; ++i) {
            if (equal(points[i], hull[curHullIdx])) {
                hullIdx[curHullIdx] = i;
                ++curHullIdx;
            }
        }
    }

    int edgeSize = 0;
    Edge *edges = getEdges(points, n, &edgeSize);
    qsort(edges, edgeSize, sizeof(Edge), compareEdge);
    
    double min = INT_MAX;
    int work_count = n > 4 ? 16 : 1;

    // array of works
    bool **inMSTs = malloc(work_count * sizeof(bool *));

    for (int i = 0; i < work_count; ++i) {
        inMSTs[i] = malloc(n * sizeof(bool));
        for (int j = 0; j < n; ++j) {
            inMSTs[i][j] = false;
        }
    }

    for (int j = 0; j < work_count; ++j) {
        inMSTs[j][0] = j & 1;
        inMSTs[j][1] = (j >> 1) & 1;
        inMSTs[j][2] = (j >> 2) & 1;
        inMSTs[j][3] = (j >> 3) & 1;
    }
    
#pragma omp parallel num_threads(thread_num)
    {
    Edge *resultMST = malloc(edgeSize * sizeof(Edge));
    int MSTSize;
    int start_idx = n > 4 ? 4 : 0;
#pragma omp for 
        for (int i = 0; i < work_count; ++i) {
            MST_Node *local_nodes = malloc(n * sizeof(MST_Node));
            double my_min = bruteForce(inMSTs[i], local_nodes, n, start_idx, hullIdx, hullSize, edges, edgeSize, resultMST, &MSTSize, min);
            #pragma omp critical
            {
                if (my_min < min) {
                    min = my_min;
                }
            }
        }
    }

    printf("%.4f", min);

    return 0;
}