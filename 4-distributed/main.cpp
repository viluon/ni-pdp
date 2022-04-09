#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <cstdint>
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <deque>

namespace mpi = boost::mpi;

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

// TODO: use Boost MPI for serialization of C++ objects
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

void solve(Sln* result, u32 a, const vector<vector<u32>>& graph, const Sln& sln, Sln best, u32 node) {
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

    solve(&incl_sln, a, graph, sln.include(node), best, node + 1);
    const auto new_best = pick(incl_sln, best);
    solve(&excl_sln, a, graph, sln.exclude(node), new_best, node + 1);
    *result = pick(excl_sln, new_best);
    #pragma omp critical
    {
        if (result->valid && (!optimum.valid || result->weight < optimum.weight)) {
            optimum = *result;
        }
    }
}

deque<pair<u16, Sln>> gen_initial_configurations(const vector<vector<u32>>& graph) {
    constexpr u8 p = 4;
    constexpr u8 z = 15; // found optimal with Hyperfine
    auto q = deque<pair<u16, Sln>>();
    q.emplace_back(make_pair(0, Sln(graph)));

    while (q.size() < z * p) {
        auto pair = q.front();
        q.pop_front();
        auto node = pair.first;
        auto sln = pair.second;
        q.emplace_back(make_pair(node + 1, sln.include(node)));
        q.emplace_back(make_pair(node + 1, sln.exclude(node)));
    }

    // for (const auto& pair : q) {
    //     cout << pair.first << " ";
    //     for (auto b : pair.second.assignment) {
    //         cout << (b ? "1" : "0");
    //     }
    //     cout << endl;
    // }

    return q;
}

void solve(u32 a, const vector<vector<u32>>& graph) {
    // generate a queue of configurations
    // each thread will process one of these configurations
    auto queue = gen_initial_configurations(graph);

    #pragma omp parallel for schedule(dynamic)
    for (auto pair : queue) {
        auto node = pair.first;
        auto sln = pair.second;
        solve(&sln, a, graph, sln, sln, node);
    }
}

void master(mpi::communicator& world) {
    auto graph = load_input();
    u32 a = 10;
    solve(a, graph);

    auto sln = Sln(graph);
    #pragma omp critical
    {
        sln = optimum;
    }

    cout << sln.weight << " (" << sln.n_set << ")" << endl;
    for (u32 i = 0; i < graph.size(); i++) {
        cout << sln.assignment[i] << " ";
    }
    cout << endl;
}

void slave(mpi::communicator& world) {
    // TODO
}

int main() {
    mpi::environment env;
    mpi::communicator world;
    cout
        << "I am process " << world.rank()
        << " of " << world.size()
        << "." << endl;

    if (world.rank() == 0) {
        master(world);
    } else {
        slave(world);
    }

    return 0;
}
