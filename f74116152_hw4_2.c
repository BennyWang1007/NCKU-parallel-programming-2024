#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <pthread.h>


int n;
int *neurons;
int **dp_result;

int thread_count;

void readData(char *filename) {
    FILE *fp = fopen(filename, "r");
    fscanf(fp, "%d", &n);
    n += 1;
    neurons = (int *)malloc(sizeof(int) * n);
    for (int i = 1; i < n; i++) {
        fscanf(fp, "%d", &neurons[i]);
    }
    neurons[0] = 1;
}

void init_result() {
    dp_result = (int **)malloc(sizeof(int *) * n);
    for (int i = 0; i < n; i++) {
        dp_result[i] = (int *)malloc(sizeof(int) * n);
        dp_result[i][i] = 0;
        for (int j = i + 2; j < n; j++) {
            dp_result[i][j] = INT_MAX;
        }
    }
    // for (int l = 2; l < n; l++) {
    //     for (int i = 0; i < n - l; i++) {
    //         dp_result[i][i + l] = INT_MAX;
    //     }
    // }

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
pthread_mutex_t barrier_mutex, result_mutex;
pthread_cond_t barrier_cond;

void barrier() {
    // implement your barrier here
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
    // implement your barrier here
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

int **s;

void* thread_runner(void *args) {

    int rank = *(int *)args;

    // // single: about 0.0012s
    // // s = mallocMatrix(n, n);
    // if (rank != 0) {
    //     return NULL;
    // }
    // for (int l = 2; l < n; l++) {
    //     for (int i = 0; i < n - l; i++) {
    //         int j = i + l;
    //         // dp_result[i][j] = INT_MAX;
    //         for (int k = i+1; k < j; k++) {
    //             int q = dp_result[i][k] + dp_result[k][j] + neurons[i] * neurons[k] * neurons[j];
    //             if (q < dp_result[i][j]) {
    //                 dp_result[i][j] = q;
    //                 // s[i][j] = k;
    //             }
    //         }
    //     }
    // }

    // multi: about 0.0025s
    for (int l = 2; l < n; l++) {
        for (int i = 0; i < n - l; i++) {
            int j = i + l;
            int local_min = INT_MAX;
            for (int k = i + 1 + rank; k < j; k += thread_count) {
                int q = dp_result[i][k] + dp_result[k][j] + neurons[i] * neurons[k] * neurons[j];
                if (q < local_min) {
                    local_min = q;
                }
            }
            // pthread_mutex_lock(&result_mutex);
            // if (local_min < dp_result[i][j]) {
            //     dp_result[i][j] = local_min;
            // }
            // pthread_mutex_unlock(&result_mutex);
            // barrier();
            
            pthread_mutex_lock(&barrier_mutex);
            barrier_counter++;
            if (local_min < dp_result[i][j]) {
                dp_result[i][j] = local_min;
            }
            if (barrier_counter == thread_count) {
                barrier_counter = 0;
                pthread_cond_broadcast(&barrier_cond);
            } else {
                pthread_cond_wait(&barrier_cond, &barrier_mutex);
            }
            pthread_mutex_unlock(&barrier_mutex);
        }
    }
    
    return NULL;
}

void printMatrix(int **arr, int h, int w) {
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            printf("%d ", arr[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {

    long thread;
    pthread_t *thread_handles;

    thread_count = strtol(argv[1], NULL, 10);
    thread_handles = malloc(thread_count * sizeof(pthread_t));

    char filename[100];
    scanf("%s", filename);
    readData(filename);
    init_result();

    // printMatrix(dp_result, n, n);

    pthread_mutex_init(&barrier_mutex, NULL);
    pthread_mutex_init(&result_mutex, NULL);
    pthread_cond_init(&barrier_cond, NULL);

    int *rank = (int *)malloc(sizeof(int) * thread_count);
    for (int i = 0; i < thread_count; i++) {
        rank[i] = i;
        pthread_create(&thread_handles[i], NULL, thread_runner, &rank[i]);
    }

    for (int i = 0; i < thread_count; i++) {
        pthread_join(thread_handles[i], NULL);
    }

    printf("%d", dp_result[0][n-1]);
    // printMatrix(dp_result, n, n);
    
    return 0;
}