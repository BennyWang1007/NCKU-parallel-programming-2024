#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

int rank, procs;

typedef struct {
    unsigned short from;
    unsigned short dest;
    unsigned short cost;
} MoveCost;

typedef struct {
    int id;
    int dist;
} Node;

typedef struct {
    int size;
    int capacity;
    Node *nodes;
} Heap;

void heapInit(Heap *heap, int capacity, int **id_to_idx_ptr, int n) {
    heap->size = 0;
    heap->capacity = capacity;
    heap->nodes = (Node *)malloc(sizeof(Node) * capacity);
    int *id_to_idx = (int *)malloc(sizeof(int) * n);
    for (int i = 0; i < n; i++) id_to_idx[i] = -1; 
    *id_to_idx_ptr = id_to_idx;
}

Node heapTop(Heap *heap) {
    if (heap->size == 0) return (Node){-1, -1};
    return heap->nodes[0];
}

void minHeapify(Heap *heap, int idx, int *id_to_idx) {
    int left = idx * 2 + 1;
    int right = idx * 2 + 2;
    int smallest = idx;
    if (left < heap->size && heap->nodes[left].dist < heap->nodes[smallest].dist) {
        smallest = left;
    }
    if (right < heap->size && heap->nodes[right].dist < heap->nodes[smallest].dist) {
        smallest = right;
    }
    if (smallest != idx) {
        id_to_idx[heap->nodes[smallest].id] = idx;
        id_to_idx[heap->nodes[idx].id] = smallest;
        Node temp = heap->nodes[smallest];
        heap->nodes[smallest] = heap->nodes[idx];
        heap->nodes[idx] = temp;
        minHeapify(heap, smallest, id_to_idx);
    }
}

Node heapPop(Heap *heap, int *id_to_idx) {
    if (heap->size == 0) {
        return (Node){-1, -1};
    }
    Node ret = heap->nodes[0];
    if (heap->size == 1) {
        id_to_idx[ret.id] = -1;
        heap->size--;
        return ret;
    }
    heap->nodes[0] = heap->nodes[heap->size - 1];
    heap->size--;
    minHeapify(heap, 0, id_to_idx);
    return ret;
}

int heapSearch(Heap *heap, int id, int *id_to_idx) {
    return id_to_idx[id];
}

void heapDecreaseKey(Heap *heap, int id, int dist, int *id_to_idx) {
    int idx = heapSearch(heap, id, id_to_idx);
    if (idx == -1) {
        printf("(proc %d)Error: id %d not found\n", rank, id);
        return;
    }
    Node node = heap->nodes[idx];
    node.dist = dist;
    while (idx > 0) {
        int parent = (idx - 1) / 2;
        if (heap->nodes[parent].dist <= node.dist) break;
        id_to_idx[heap->nodes[parent].id] = idx;
        heap->nodes[idx] = heap->nodes[parent];
        idx = parent;
    }
    heap->nodes[idx] = node;
    id_to_idx[node.id] = idx;
}

int compareMoveCostFrom(const void *a, const void *b) {
    MoveCost *ma = (MoveCost *)a;
    MoveCost *mb = (MoveCost *)b;
    if (ma->from == mb->from) {
        return ma->dest - mb->dest;
    }
    return ma->from - mb->from;
}

int compareMoveCostDest(const void *a, const void *b) {
    MoveCost *ma = (MoveCost *)a;
    MoveCost *mb = (MoveCost *)b;
    if (ma->dest == mb->dest) {
        return ma->from - mb->from;
    }
    return ma->dest - mb->dest;
}

void readData(char *filename, int *n, MoveCost *moveCosts) {
    FILE *fp = fopen(filename, "r");
    fscanf(fp, "%d", n);
    int idx = 0;
    int from, dest, cost;
    while (!feof(fp)) {
        fscanf(fp, "%d %d %d", &from, &dest, &cost);
        moveCosts[idx] = (MoveCost){from, dest, cost};
        idx++;
    }
    fclose(fp);
}

int getPathCount(char *filename) {
    FILE *fp = fopen(filename, "r");
    int n;
    fscanf(fp, "%d", &n);
    int from, dest, cost;
    int count = 0;
    while (!feof(fp)) {
        fscanf(fp, "%d %d %d", &from, &dest, &cost);
        count++;
    }
    fclose(fp);
    return count;
}

void relex(int *distTable, MoveCost *moveCosts, int *moveCostStartIdx, int node, int node_cnt, int pathCount, Heap *heap, int *id_to_idx) {
    int startIdx = moveCostStartIdx[node];
    int from, dest, cost;
    for (int k = startIdx; k < pathCount; k++) {
        from = moveCosts[k].from;
        dest = moveCosts[k].dest;
        cost = moveCosts[k].cost;
        if (from != node) break;
        if (distTable[dest] > distTable[from] + cost) {
            distTable[dest] = distTable[from] + cost;
            heapDecreaseKey(heap, dest, distTable[dest], id_to_idx);
        }
    }
}

void printHeap(Heap *heap) {
    for (int i = 0; i < heap->size; i++) {
        printf("%d %d\n", heap->nodes[i].id, heap->nodes[i].dist);
    }
}

bool checkHeap(Heap *heap) {
    for (int i = 0; i < heap->size; i++) {
        int left = i * 2 + 1;
        int right = i * 2 + 2;
        if (left < heap->size && heap->nodes[i].dist > heap->nodes[left].dist) {
            return false;
        }
        if (right < heap->size && heap->nodes[i].dist > heap->nodes[right].dist) {
            return false;
        }
    }
    return true;
}

int compareNode(const void *a, const void *b) {
    Node *na = (Node *)a;
    Node *nb = (Node *)b;
    return na->dist - nb->dist;
}

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    int n, pathCount;
    MoveCost *sparsedMoveCost;
    int *moveCostStartIdx;

    if (rank == 0) {
        char filename[100];
        scanf("%s", filename);
        pathCount = getPathCount(filename);
        sparsedMoveCost = (MoveCost *)malloc(sizeof(MoveCost) * pathCount);
        readData(filename, &n, sparsedMoveCost);
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int color = (rank < n) ? 0 : MPI_UNDEFINED;
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &comm);
    
    // make procs <= n
    if (rank >= n) {
        MPI_Finalize();
        return 0;
    }
    procs = procs < n ? procs : n;

    // send and preprocess edges
    MPI_Bcast(&pathCount, 1, MPI_INT, 0, comm);
    moveCostStartIdx = (int *)malloc(sizeof(int) * n);
    if (rank != 0) {
        sparsedMoveCost = (MoveCost *)malloc(sizeof(MoveCost) * pathCount);
    }

    MPI_Datatype MPI_MoveCost;
    MPI_Type_contiguous(3, MPI_UNSIGNED_SHORT, &MPI_MoveCost);
    MPI_Type_commit(&MPI_MoveCost);

    MPI_Bcast(sparsedMoveCost, pathCount, MPI_MoveCost, 0, comm);
    int localStartId = rank * n / procs, localEndId = (rank + 1) * n / procs;
    int localPathCount = 0;
    MoveCost *localMoveCosts = (MoveCost *)malloc(sizeof(MoveCost) * pathCount);
    for (int i = 0; i < pathCount; i++) {
        if (sparsedMoveCost[i].dest >= localStartId && sparsedMoveCost[i].dest < localEndId) {
            localMoveCosts[localPathCount] = sparsedMoveCost[i];
            localPathCount++;
        }
    }
    qsort(localMoveCosts, localPathCount, sizeof(MoveCost), compareMoveCostFrom);

    for (int i = 0; i < n; i++) {
        moveCostStartIdx[i] = -1;
    }
    for (int i = 0; i < localPathCount; i++) {
        if (moveCostStartIdx[localMoveCosts[i].from] == -1) {
            moveCostStartIdx[localMoveCosts[i].from] = i;
        }
    }

    // set distTableFrom0 to infinity except distTableFrom0[0] = 0
    int *distTableFrom0 = (int *)calloc(n, sizeof(int));
    for (int i = 1; i < n; i++) {
        distTableFrom0[i] = INT32_MAX;
    }
    // make heap
    Heap heap;
    int *id_to_idx;
    int heapSize = (rank + 1) * n / procs - rank * n / procs;
    heapInit(&heap, heapSize, &id_to_idx, n);
    heap.size = heapSize;
    for (int i = 0; i < heapSize; i++) {
        Node node = {rank * n / procs + i, INT32_MAX};
        heap.nodes[i] = node;
        id_to_idx[node.id] = i;
    }
    // set dist of node 0 to 0
    if (rank == 0) heap.nodes[0].dist = 0;

    for (int i = 0; i < n; i++) {

        Node node = heapTop(&heap);
        int minCost2 = INT32_MAX;
        int minId2 = -1;

        if (node.dist != -1 && node.dist < minCost2) {
            minCost2 = node.dist;
            minId2 = node.id;
        }

        int globalMinCost2, globalMinId2=-1;
        MPI_Allreduce(&minCost2, &globalMinCost2, 1, MPI_INT, MPI_MIN, comm);
        if (globalMinCost2 != minCost2) minId2 = -1;
        MPI_Allreduce(&minId2, &globalMinId2, 1, MPI_INT, MPI_MAX, comm);

        if (globalMinId2 == -1) break;

        if (globalMinId2 == minId2) {
            heapPop(&heap, id_to_idx);
        }

        MPI_Allreduce(MPI_IN_PLACE, &distTableFrom0[globalMinId2], 1, MPI_INT, MPI_MIN, comm);
        // relex(distTableFrom0, localMoveCosts, moveCostStartIdx, globalMinId, n, localPathCount);
        relex(distTableFrom0, localMoveCosts, moveCostStartIdx, globalMinId2, n, localPathCount, &heap, id_to_idx);
    }
    MPI_Allreduce(MPI_IN_PLACE, distTableFrom0, n, MPI_INT, MPI_MIN, comm);

    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            printf("%d ", distTableFrom0[i]);
        }
    }

    MPI_Type_free(&MPI_MoveCost);

    free(sparsedMoveCost);
    free(moveCostStartIdx);
    free(distTableFrom0);
    free(localMoveCosts);
    free(id_to_idx);
    free(heap.nodes);

    // 0.5-0.6 sec in judge.conf
    MPI_Finalize();
    return 0;
}