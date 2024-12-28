#include <algorithm>
#include <cstdio>
#include <list>
#include <vector>
#include <stdint.h>
#include <fstream>
#include <iostream>

using std::ifstream, std::vector, std::cin, std::cout, std::endl, std::list, std::string;

struct Fraction {
    long a, b;
    Fraction (long aa, long bb = 1) : a(aa), b(bb) {}
    bool operator < (Fraction rhs) const {
        return (__int128) a * rhs.b < (__int128) b * rhs.a;
    }
    bool operator > (Fraction rhs) const { return rhs < *this; }
    bool operator <= (Fraction rhs) const { return !(rhs < *this); }
    bool operator >= (Fraction rhs) const { return !(*this < rhs); }
};

struct Arc {
    size_t first, second;
    list<Arc *> ceiling;
    list<Arc *>::iterator hm_it;
    long numerator, denominator;
    Fraction support() const {
        return Fraction(numerator, denominator);
    }
    void init() {
        hm_it = std::max_element(ceiling.begin(), ceiling.end(),
                [](const Arc *lhs, const Arc *rhs) {
                    return lhs->support() < rhs->support();
                });
    }
    const Arc *get_hm() const {
        return *hm_it;
    }
    Arc *get_hm() {
        return const_cast<Arc *>(const_cast<const Arc *>(this)->get_hm());
    }
    void merge_hm() {
        auto it = hm_it;
        Arc *hm = *it;
        it = ceiling.erase(it);
        ceiling.splice(it, hm->ceiling);
        hm_it = std::max_element(ceiling.begin(), ceiling.end(),
                [](const Arc *lhs, const Arc *rhs) {
                    return lhs->support() < rhs->support();
                });
    }
};

void get_arcs(const vector<long> &a, vector<Arc> &arcs) {
    size_t n = a.size() - 1;
    vector<long> v;
    vector<Arc *> w;
    for (size_t i = 0, j = 0; i <= n && j < n - 3; i++) {
        while (j < n - 3 && v.size() >= 2 && a[i] <= a[v.back()]) {
            arcs[j].first = v[v.size() - 2];
            arcs[j].second = i;
            while (!w.empty() &&
                    arcs[j].first <= w.back()->first &&
                    w.back()->second <= arcs[j].second) {
                arcs[j].ceiling.push_front(w.back());
                w.pop_back();
            }
            w.push_back(&arcs[j]);
            j++;
            v.pop_back();
        }
        v.push_back(i);
    }

    arcs[n - 3].first = 0;
    arcs[n - 3].second = n;
    arcs[n - 3].ceiling.assign(w.begin(), w.end());
}


int n;
int thread_count;
long result = 0;

vector<long> neurons;

void* solve(void *arg) {
    long rank = (long)arg;
    // cout << "thread " << (long)rank << " start" << endl;
    if (rank != 0) {
        return NULL;
    }
    size_t n = neurons.size();
    if (n <= 2)
        return NULL;
    if (n == 3) {
        result = neurons[0] * neurons[1] * neurons[2];
        return NULL;
    }
    std::rotate(neurons.begin(), std::min_element(neurons.begin(), neurons.end()), neurons.end());
    neurons.push_back(neurons[0]);

    vector<long> accum(n + 1);
    for (size_t i = 1; i <= n; i++)
        accum[i] = accum[i - 1] + neurons[i] * neurons[i - 1];

    vector<Arc> arcs(n - 2);
    get_arcs(neurons, arcs);

    long ans = 0;
    for (Arc &arc : arcs) {
        if (arc.first + 2 == arc.second) {
            // leaf nodes
            arc.numerator = neurons[arc.first] * neurons[arc.first + 1] * neurons[arc.second];
            arc.denominator = accum[arc.second] - accum[arc.first] -
                neurons[arc.first] * neurons[arc.second];
            ans += arc.numerator;
            continue;
        }

        arc.init();

        arc.denominator = accum[arc.second] - accum[arc.first] -
            neurons[arc.first] * neurons[arc.second];
        for (Arc *jp : arc.ceiling) {
            size_t jf = jp->first, js = jp->second;
            arc.denominator -= accum[js] - accum[jf] - neurons[jf] * neurons[js];
        }

        // step 1
        while (!arc.ceiling.empty() && arc.get_hm()->support() >=
                std::min(neurons[arc.first], neurons[arc.second])) {
            Arc &hm = *arc.get_hm();
            ans -= hm.numerator;
            arc.denominator += hm.denominator;
            arc.merge_hm();
        }

        // calculate support
        size_t c1 = neurons[arc.first] <= neurons[arc.second] ? arc.first : arc.second;
        size_t c2 = c1 == 0 ? n : c1;
        arc.numerator = arc.denominator + neurons[arc.first] * neurons[arc.second];
        if (arc.first == c1)
            arc.numerator -= neurons[c1] * neurons[c1 + 1];
        if (arc.second == c2)
            arc.numerator -= neurons[c2] * neurons[c2 - 1];

        if (!arc.ceiling.empty()) {
            Arc *jp = arc.ceiling.front();
            if (jp->first == c1) {
                arc.numerator += neurons[c1] * neurons[c1 + 1];
                arc.numerator -= neurons[jp->first] * neurons[jp->second];
            }
            Arc *jq = arc.ceiling.back();
            if (jq->second == c2) {
                arc.numerator += neurons[c2] * neurons[c2 - 1];
                arc.numerator -= neurons[jq->first] * neurons[jq->second];
            }
        }
        arc.numerator *= neurons[c1];
        ans += arc.numerator;

        // step 2
        while (!arc.ceiling.empty() && arc.support() <= arc.get_hm()->support()) {
            Arc &hm = *arc.get_hm();
            arc.numerator += hm.numerator;
            arc.denominator += hm.denominator;
            arc.merge_hm();
        }
    }

    result = ans;
    return NULL;
}


void readData(string filename) {
    ifstream file(filename);
    file >> n;
    n += 1;
    neurons = vector<long>(n);
    neurons[0] = 1;
    for (int i = 1; i < n; i++) {
        file >> neurons[i];
    }
    file.close();
}

int main(int argc, char **argv) {

    thread_count = strtol(argv[1], NULL, 10);
    vector<pthread_t> thread_handles(thread_count);

    string filename;
    cin >> filename;
    
    readData(filename);

    for (long rank = 0; rank < thread_count; rank++) {
        pthread_create(&thread_handles[rank], NULL, solve, (void *)rank);
    }
    for (auto &thread_handle : thread_handles) {
        pthread_join(thread_handle, NULL);
    }

    cout << result;
}