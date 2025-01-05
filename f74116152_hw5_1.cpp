#include <iostream>
#include <limits.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <time.h>
#include <algorithm>
#include <cstring> // memcpy

using namespace std;

// overloading cout << for vector
// for printing out all permutation case at line 40-42
// template <typename T>
// ostream &operator<<(ostream &os, const vector<T> &v) {
//     os << "[";
//     for (int i = 0; i < v.size(); ++i) {
//         os << v[i];
//         if (i != v.size() - 1)
//             os << ", ";
//     }
//     os << "]";
//     return os;
// }

typedef struct Task {
    int cost;
    int release_time;
    int weight;
} Task;

int *gantt_chart_template;

// int permutation_idx = 0;
// void gen_all_permute(int *permutation, int n, int depth, vector<int*> &all_permutate) {
//     if (depth == n) {
//         int *tmp = new int[n];
//         memcpy(tmp, permutation, n * sizeof(int));
//         all_permutate.emplace_back(tmp);
//         permutation_idx++;
//         return;
//     }
//     for (int i = depth; i < n; i++) {
//         swap(permutation[depth], permutation[i]);
//         gen_all_permute(permutation, n, depth + 1, all_permutate);
//         swap(permutation[depth], permutation[i]);
//     }
// }

inline int read_data(string input_file_name, Task **task_list_ptr, int &task_num, int &max_ri, int &sum_of_all_pi) {
    std::ifstream myfile;
    myfile.open(input_file_name);
    myfile >> task_num;

    max_ri = 0;
    sum_of_all_pi = 0;

    Task *task_list = new Task[task_num];

    int cost, rel, w; // cost, release time, weight

    for (int i = 0; i < task_num; i++) {
        myfile >> cost >> rel >> w;
        task_list[i].cost = cost;
        task_list[i].release_time = rel;
        task_list[i].weight = w;
        if (rel > max_ri) max_ri = rel;
        sum_of_all_pi += cost;
    }
    *task_list_ptr = task_list;
    myfile.close();
    return 0;
}

inline void gen_permute_start(int *permutation, int n, int **all_permutate_start) {
    for (int i = 0; i < n; i++) {
        swap(permutation[0], permutation[i]);
        int *tmp = new int[n];
        memcpy(tmp, permutation, n * sizeof(int));
        all_permutate_start[i] = tmp;
        swap(permutation[0], permutation[i]);
    }
}

inline bool test_permute(int *permutation, int n, Task *task_list, int &min, int *schedule_gantt_chart) {
    int cur_total_weighted_cost = 0;
    for (int j = 0; j < n; j++) {
        int cur_task_id = permutation[j];
        int cur_pi = task_list[cur_task_id].cost;
        int cur_ri = task_list[cur_task_id].release_time;

        while (cur_pi) {
            if (schedule_gantt_chart[cur_ri] == -1) {
                schedule_gantt_chart[cur_ri] = cur_task_id;
                cur_pi--;
            }
            cur_ri++;
        }
        // after fill current taks into schedule_gantt_chart,
        // calculate and sum to current total weighted cost
        if (cur_pi == 0){
            cur_total_weighted_cost += cur_ri * task_list[cur_task_id].weight;
        }
        if (cur_total_weighted_cost > min) {
            return false;
        }
    }

    // update min
    min = cur_total_weighted_cost;
    return true;
}

// return the weighted cost of the task at idx
inline bool test_permute_with_index(int *permutation, int idx, int n, Task *task_list, const int min, int *gantt_chart_buffer, int gantt_chart_length) {
    int cur_wc = 0;
    memcpy(gantt_chart_buffer, gantt_chart_template, gantt_chart_length * sizeof(int));

    for (int j = 0; j < n; j++) {
        int cur_task_id = permutation[j];
        int cur_pi = task_list[cur_task_id].cost;
        int cur_ri = task_list[cur_task_id].release_time;

        while (cur_pi) {
            if (gantt_chart_buffer[cur_ri] == -1) {
                gantt_chart_buffer[cur_ri] = cur_task_id;
                cur_pi--;
            }
            cur_ri++;
        }
        // after fill current taks into schedule_gantt_chart,
        // calculate and sum to current total weighted cost
        if (cur_pi == 0)
            cur_wc += cur_ri * task_list[cur_task_id].weight;

        if (cur_wc > min)
            return false;
            
        if (j == idx)
            return true;
    }
    return true;
}

inline void swap(int &a, int &b) {
    int tmp = a;
    a = b;
    b = tmp;
}


void test_all_permute_with_start(int *permute_start, const int n, int depth, Task *task_list, int &min_weighted_cost, int *gantt_chart_buffer, const int gantt_chart_length) {
    if (depth == n) {
        memcpy(gantt_chart_buffer, gantt_chart_template, gantt_chart_length * sizeof(int));
        test_permute(permute_start, n, task_list, min_weighted_cost, gantt_chart_buffer);
        return;
    }

    for (int i = depth; i < n; i++) {
        swap(permute_start[depth], permute_start[i]);
        // pruning
        if (!test_permute_with_index(permute_start, depth, n, task_list, min_weighted_cost, gantt_chart_buffer, gantt_chart_length)) {
            swap(permute_start[depth], permute_start[i]);
            break;
        }
        test_all_permute_with_start(permute_start, n, depth + 1, task_list, min_weighted_cost, gantt_chart_buffer, gantt_chart_length);
        swap(permute_start[depth], permute_start[i]);
    }
}

int main(int argc, char **argv) {

    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);

    int thread_count = atoi(argv[1]);

    string input_file_name = "";
    cin >> input_file_name;

    std::ifstream myfile;
    myfile.open(input_file_name);

    int task_num = 0;
    myfile >> task_num;

    int max_ri = 0;
    int sum_of_all_pi = 0;
    Task *task_list;
    read_data(input_file_name, &task_list, task_num, max_ri, sum_of_all_pi);

    int global_min_weighted_cost = INT_MAX;
    int schedule_gantt_chart_length = max_ri + sum_of_all_pi;
    
    gantt_chart_template = new int[schedule_gantt_chart_length];
    for (int i = 0; i < schedule_gantt_chart_length; i++) {
        gantt_chart_template[i] = -1;
    }

    // split the task_list into task_num parts
    int **all_permute_start;
    all_permute_start = new int*[task_num];
    int permutation_case[task_num];

    for (int i = 0; i < task_num; i++) {
        permutation_case[i] = i;
    }
    gen_permute_start(permutation_case, task_num, all_permute_start);

#pragma omp parallel num_threads(thread_count)
    {
        int *gantt_chart_buffer = new int[schedule_gantt_chart_length];
#pragma omp for
        for (int i = 0; i < task_num; i++) {
            int local_min = INT_MAX;
            test_all_permute_with_start(all_permute_start[i], task_num, 1, task_list, local_min, gantt_chart_buffer, schedule_gantt_chart_length);
#pragma omp critical
            {
                if (global_min_weighted_cost > local_min)
                    global_min_weighted_cost = local_min;
            }
        }
    }
    cout << global_min_weighted_cost;
    return 0;
}