
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int n;
int *neurons;
int result = 0;
int thread_count;

// read neurons from file to neurons[0..n]
void readData(char *filename) {
    FILE *fp = fopen(filename, "r");
    fscanf(fp, "%d", &n);
    n += 1;
    neurons = (int *)malloc(sizeof(int) * n);
    neurons[0] = 1;
    for (int i = 1; i < n; i++) {
        fscanf(fp, "%d", &neurons[i]);
    }
    fclose(fp);
}

pthread_mutex_t result_mutex;

void* thread_runner(void *args) {

    int rank = *(int *)args;
    int local_result = 0;

    // Divide work [1, n-1] among threads
    int local_start = rank * (n - 2) / thread_count + 1;
    int local_end = (rank + 1) * (n - 2) / thread_count + 1;

    for (int i = local_start; i < local_end; ++i) {
        local_result += neurons[i] * neurons[i + 1];
    }

    pthread_mutex_lock(&result_mutex);
    result += local_result;
    pthread_mutex_unlock(&result_mutex);

    return NULL;
}

int main(int argc, char **argv) {

    thread_count = strtol(argv[1], NULL, 10);

    char filename[100];
    scanf("%s", filename);
    readData(filename);

    pthread_t *threads;
    threads = malloc(thread_count * sizeof(pthread_t));
    pthread_mutex_init(&result_mutex, NULL);

    int *ranks = (int *)malloc(sizeof(int) * thread_count);
    for (int i = 0; i < thread_count; i++) {
        ranks[i] = i;
        pthread_create(&threads[i], NULL, thread_runner, &ranks[i]);
    }

    for (int i = 0; i < thread_count; i++) {
        pthread_join(threads[i], NULL);
    }

    printf("%d", result);

    pthread_mutex_destroy(&result_mutex);
    free(neurons);
    free(ranks);
    free(threads);
    
    return 0;
}