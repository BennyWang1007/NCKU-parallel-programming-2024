#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

static inline void decompose_gate(long gate, long *gate_and, long *gate_dcare) {
    long gand = 0, gdcare = 0;
    for (int i = 0; i < 20; ++i) {
        gand <<= 1;
        gdcare <<= 1;
        int r = gate % 3;
        gate /= 3;
        if (r == 0) {
            ++gdcare;
        } else if (r == 1) {
            ++gand;
        }
    }
    *gate_and = gand;
    *gate_dcare = gdcare;
    return;
}

int read_file(int *n_ptr, long **gate_and_ptr, long **gate_dcare_ptr) {
    char filename[50];
    scanf("%s", filename);

    FILE *inFile = fopen(filename, "r");
    if (inFile == NULL) {
        printf("Could not open the file %s\n", filename);
        return 1;
    }

    int n;
    fscanf(inFile, "%d", &n);

    long gate = 0;
    long *gate_and = (long*)malloc(n * sizeof(long));
    long *gate_dcare = (long*)malloc(n * sizeof(long));

    for(int i = 0; i < n; ++i) {
        fscanf(inFile, "%ld", &gate);
        decompose_gate(gate, &gate_and[i], &gate_dcare[i]);
    }

    *n_ptr = n;
    *gate_and_ptr = gate_and;
    *gate_dcare_ptr = gate_dcare;

    fclose(inFile);
    return 0;
}

// void print_binary(long n) {
//     for (int i = 31; i >= 0; --i) {
//         printf("%ld", (n >> i) & 1);
//     }
// }

long find_accept_count(long *gate_and, long *gate_dcare, int n, long start, long end) {
    long count = 0;
    // long accept = -1;
    for (long in = start; in < end; ++in) {
        for (int i = 0; i < n; ++i) {
            if ((~(in ^ gate_and[i]) | gate_dcare[i]) == -1) {
                ++count;
                break;
            }
        }
    }
    return count;
}

pthread_mutex_t count_mutex;

typedef struct {
    long *gate_and;
    long *gate_dcare;
    int n;
    long start;
    long end;
    int thread_id;
    int thread_count;
} ThreadArgs;

long count = 0;

void *find_accept_count_thread(void *args) {
    ThreadArgs *targs = (ThreadArgs*)args;
    long *gate_and = targs->gate_and;
    long *gate_dcare = targs->gate_dcare;
    int n = targs->n;
    long start = targs->start;
    long end = targs->end;
    int thread_id = targs->thread_id;
    int thread_count = targs->thread_count;

    long my_start = start + thread_id * (end - start) / thread_count;
    long my_end = start + (thread_id + 1) * (end - start) / thread_count;
    long my_count = find_accept_count(gate_and, gate_dcare, n, my_start, my_end);

    pthread_mutex_lock(&count_mutex);
    count += my_count;
    pthread_mutex_unlock(&count_mutex);
    
    return NULL;
}

int main(int argc, char *argv[]) {

    int rank, numprocs;
    int thread_count = atoi(argv[1]);

    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // initialize and read file
    int n;
    long *gate_and, *gate_dcare;

    if (rank == 0) {
        if (read_file(&n, &gate_and, &gate_dcare) != 0) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // broadcast
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank != 0) {
        gate_and = malloc(n * sizeof(long));
        gate_dcare = malloc(n * sizeof(long));
    }

    MPI_Bcast(gate_and, n, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(gate_dcare, n, MPI_LONG, 0, MPI_COMM_WORLD);

    long end = (1 << 20);
    long start_x = rank * end / numprocs;
    long end_x = (rank+1) * end / numprocs;
    
    pthread_mutex_init(&count_mutex, NULL);
    pthread_t *threads = malloc(thread_count * sizeof(pthread_t));
    for (int i = 0; i < thread_count; ++i) {
        ThreadArgs *args = malloc(sizeof(ThreadArgs));
        args->gate_and = gate_and;
        args->gate_dcare = gate_dcare;
        args->n = n;
        args->start = start_x;
        args->end = end_x;
        args->thread_id = i;
        args->thread_count = thread_count;
        pthread_create(&threads[i], NULL, find_accept_count_thread, args);
    }
    for (int i = 0; i < thread_count; ++i) {
        pthread_join(threads[i], NULL);
    }

    long total_count;

    MPI_Reduce(&count, &total_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("%ld", total_count);
    }

    pthread_mutex_destroy(&count_mutex);
    free(threads);
    free(gate_and);
    free(gate_dcare);
    
    MPI_Finalize();

    return 0;

}