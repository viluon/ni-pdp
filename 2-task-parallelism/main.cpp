#include <cstdint>
#include <vector>
#include <iostream>
#include <algorithm>

// typedefs for speed, compactness, and clarity
typedef uint_fast64_t u64;
typedef  int_fast64_t i64;
typedef uint_fast32_t u32;
typedef  int_fast32_t i32;
typedef uint_fast16_t u16;
typedef  int_fast16_t i16;
typedef uint_fast8_t   u8;
typedef  int_fast8_t   i8;

using namespace std;


vector<vector<u32>> load_input() {
    vector<vector<u32>> graph;
    u32 n;

    cin >> n;
    graph.reserve(n);
    for (u32 i = 0; i < n; ++i) {
        graph.emplace_back();
        graph[i].reserve(n);
        u32 w;
        for (u32 j = 0; j < n; j++) {
            cin >> w;
            graph[i].emplace_back(w);
        }
    }

    return graph;
}

struct Sln {
    const vector<vector<u32>>& graph;
    vector<bool> assignment;
    bool valid;
    u32 weight;
    u32 n_set;

    explicit Sln(const vector<vector<u32>>& graph)
    : graph(graph)
    , assignment(graph.size(), false)
    , valid(false)
    , weight(0)
    , n_set(0)
    {}

    Sln(const Sln& other)
    : graph(other.graph)
    , assignment(other.assignment)
    , valid(other.valid)
    , weight(other.weight)
    , n_set(other.n_set)
    {}

    Sln(Sln&& other)
    : graph(other.graph)
    , assignment(move(other.assignment))
    , valid(other.valid)
    , weight(other.weight)
    , n_set(other.n_set)
    {}

    Sln& operator=(const Sln& other) {
        assignment = other.assignment;
        valid = other.valid;
        weight = other.weight;
        n_set = other.n_set;
        return *this;
    }

    Sln& operator=(Sln&& other) {
        assignment = move(other.assignment);
        valid = other.valid;
        weight = other.weight;
        n_set = other.n_set;
        return *this;
    }

    Sln include(u32 node) const { return this->update(node, true); }
    Sln exclude(u32 node) const { return this->update(node, false); }

    Sln update(u32 node, bool included) const {
        Sln result = *this;
        result.n_set += included - assignment[node];
        result.assignment[node] = included;

        for (u32 i = 0; i < node; i++) {
            if (result.assignment[node] != result.assignment[i]) {
                result.weight += graph[node][i];
            }
        }
        return result;
    }

    // relies on invariants:
    // - the assignment vector is decomposable into an assigned prefix and an unassigned suffix
    // - the about-to-be-assigned node is passed as the sole argument
    // - the graph matrix is symmetric (the graph is undirected)
    // [ 0 1 2 3 4 5 6 7 8 9 ]
    // [ 1 0 1 0 0 0 0 0 0 0 ]
    //    node ^
    u32 lower_bound(u32 node) const {
        u32 result = this->weight;
        for (u32 i = node; i < graph.size(); i++) {
            u32 included = 0;
            u32 excluded = 0;

            for (u32 j = 0; j < node; j++) {
                if (assignment[j]) {
                    included += graph[i][j];
                } else {
                    excluded += graph[i][j];
                }
            }
            result += min(included, excluded);
        }

        return result;
    }
};

const Sln& pick(const Sln& sln_x, const Sln& sln_y) {
    if (!sln_x.valid) {
        return sln_y;
    }
    return (!sln_y.valid || sln_x.weight < sln_y.weight) ? sln_x : sln_y;
}

Sln optimum = Sln(vector<vector<u32>>());

void solve(u16 depth, Sln* result, u32 a, const vector<vector<u32>>& graph, const Sln& sln, Sln best, u32 node) {
    #pragma omp critical
    {
        best = pick(best, optimum);
    }
    if (false
    || (node >= graph.size())
    || (min(sln.n_set, graph.size() - sln.n_set) > a)
    || (best.valid && sln.weight > best.weight)
    || (best.valid && sln.lower_bound(node) > best.weight)
    ) {
        *result = sln;
        if (node >= graph.size()) {
            result->valid = (result->n_set == a) || (result->n_set == graph.size() - a);
        }
        return;
    }

    Sln incl_sln = best;
    Sln excl_sln = best;

    if (depth < 3) {
        #pragma omp task shared(graph, sln, best, node, incl_sln) firstprivate(depth, a)
        solve(depth + 1, &incl_sln, a, graph, sln.include(node), best, node + 1);

        #pragma omp task shared(graph, sln, best, node, excl_sln) firstprivate(depth, a)
        solve(depth + 1, &excl_sln, a, graph, sln.exclude(node), best, node + 1);

        #pragma omp taskwait
        *result = pick(excl_sln, pick(incl_sln, best));
        #pragma omp critical
        {
            if (result->valid && result->weight < optimum.weight) {
                optimum = *result;
            }
        }
    } else {
        solve(depth + 1, &incl_sln, a, graph, sln.include(node), best, node + 1);
        const auto new_best = pick(incl_sln, best);
        solve(depth + 1, &excl_sln, a, graph, sln.exclude(node), new_best, node + 1);
        *result = pick(excl_sln, new_best);
        #pragma omp critical
        {
            if (result->valid && result->weight < optimum.weight) {
                optimum = *result;
            }
        }
    }
}

int main() {
    auto graph = load_input();
    auto init = Sln(graph);
    auto sln = init;
    #pragma omp parallel num_threads(4)
    {
        u32 a = 15;
        #pragma omp single
        solve(0, &sln, a, graph, init, init, 0);
    }

    cout << sln.weight << " (" << sln.n_set << ")" << endl;
    for (u32 i = 0; i < graph.size(); i++) {
        cout << sln.assignment[i] << " ";
    }
    cout << endl;

    return 0;
}
