#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <pthread.h>


int **arr, **filter, **next;
int t, w, h, D1, D2;

int thread_count;

void readData(char *filename) {
    FILE *fp = fopen(filename, "r");
    fscanf(fp, "%d", &t);
    fscanf(fp, "%d %d", &h, &w);
    arr = (int **)malloc(sizeof(int *) * h);
    for (int i = 0; i < h; i++) {
        arr[i] = (int *)malloc(sizeof(int) * w);
        for (int j = 0; j < w; j++) {
            fscanf(fp, "%d", &arr[i][j]);
        }
    }
    fscanf(fp, "%d %d", &D1, &D2);
    filter = (int **)malloc(sizeof(int *) * D1);
    for (int i = 0; i < D1; i++) {
        filter[i] = (int *)malloc(sizeof(int) * D2);
        for (int j = 0; j < D2; j++) {
            fscanf(fp, "%d", &filter[i][j]);
        }
    }
    fclose(fp);
}

// calculate dot product of a and b with start point (startX, startY) and size (w, h)
int dotProduct2(int **a, int startX, int startY, int aW, int aH, int **b, int w, int h) {
    long long sum = 0;
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

int barrier_counter = 0;
pthread_mutex_t barrier_mutex;
pthread_cond_t barrier_cond;

void barrier() {
    pthread_mutex_lock(&barrier_mutex);
    barrier_counter++;
    if (barrier_counter == thread_count) {
        barrier_counter = 0;
        pthread_cond_broadcast(&barrier_cond);
    } else {
        pthread_cond_wait(&barrier_cond, &barrier_mutex);
    }
    pthread_mutex_unlock(&barrier_mutex);
}

void barrier_with_function(void (*func)()) {
    pthread_mutex_lock(&barrier_mutex);
    barrier_counter++;
    if (barrier_counter == thread_count) {
        barrier_counter = 0;
        func();
        pthread_cond_broadcast(&barrier_cond);
    } else {
        pthread_cond_wait(&barrier_cond, &barrier_mutex);
    }
    pthread_mutex_unlock(&barrier_mutex);
}

void swap_next() {
    int **tmp = arr;
    arr = next;
    next = tmp;
}

void* thread_runner(void *args) {

    int D1D2 = D1 * D2;
    int startH = -(D1 - 1) / 2;
    int startW= -(D2 - 1) / 2;
    int rank = *((int *)args);

    for (int i = 0; i < t; ++i) {
        // calculate next
        for (int y = rank; y < h; y += thread_count) {
            for (int x = 0; x < w; x++) {
                next[y][x] = dotProduct2(arr, x+startW, y+startH, w, h, filter, D2, D1) / D1D2;
            }
        }
        // barrier();
        // if (rank == 0) {
        //     swap_next()
        // }
        // barrier();
        barrier_with_function(swap_next);
    }
    return NULL;
}

int main(int argc, char **argv) {

    long thread;
    pthread_t *thread_handles;

    thread_count = strtol(argv[1], NULL, 10);
    thread_handles = malloc(thread_count * sizeof(pthread_t));

    char filename[100];
    scanf("%s", filename);
    readData(filename);
    next = mallocMatrix(h, w);

    pthread_mutex_init(&barrier_mutex, NULL);
    pthread_cond_init(&barrier_cond, NULL);

    int *rank = (int *)malloc(sizeof(int) * thread_count);
    for (int i = 0; i < thread_count; i++) {
        rank[i] = i;
        pthread_create(&thread_handles[i], NULL, thread_runner, &rank[i]);
    }

    for (int i = 0; i < thread_count; i++) {
        pthread_join(thread_handles[i], NULL);
    }

    printAnswer(arr, h, w);
    
    // judge limit: 0.73-0.74s
    return 0;
}