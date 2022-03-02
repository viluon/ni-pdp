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

struct State {
    u32 a;
    u32 n;
    u32 n_set;
    u32 best = (u32)-1;
    u32 weight = 0;
    vector<bool> assignment;
    vector<vector<u32>>& graph;

    State(u32 n, vector<vector<u32>> &graph, const vector<bool>& assignment, u32 n_set, u32 a)
        : a(a)
        , n(n)
        , n_set(n_set)
        , assignment(assignment)
        , graph(graph)
        {}

    State(const State& other)
        : a(other.a)
        , n(other.n)
        , n_set(other.n_set)
        , assignment(other.assignment)
        , graph(other.graph)
        {}

    // assignment operators
    State& operator=(const State& other) {
        this->a = other.a;
        this->n = other.n;
        this->n_set = other.n_set;
        this->best = other.best;
        this->assignment = other.assignment;
        this->graph = other.graph;
        return *this;
    }

    // comparison operators
    bool operator==(const State& other) const {
        return this->a == other.a
            && this->n == other.n
            && this->best == other.best
            && this->n_set == other.n_set
            && this->assignment == other.assignment;
    }

    bool operator!=(const State& other) const {
        return !(*this == other);
    }

    bool operator<(const State& other) const {
        return best < other.best;
    }

    bool operator>(const State& other) const {
        return best > other.best;
    }
};

State solve(State state, u32 node) {
    if (node >= state.n) { return state; }

    for (bool set : {false, true}) {
        // FIXME complete bs, assignment of the node determines which edges are cut

        state.assignment[node] = set;
        state.n_set += set ? 1 : 0;

        auto sln = solve(state, node + 1);
        // sum up the weights of edges that lead to nodes outside our set
        u32 weight_of_cut_adjacent_edges = 0;
        for (u32 i = 0; i < sln.n; i++) {
            if (sln.assignment[i] != set) {
                weight_of_cut_adjacent_edges += sln.graph[node][i];
            }
        }

        sln.weight += weight_of_cut_adjacent_edges;

        if ((sln.n_set == state.a || sln.n_set == state.n - state.a) && sln < state) {
            state.best = sln.weight;
        }
    }

    return state;
}

int main() {
    auto graph = load_input();
    u32 a = 5;
    auto initial = State(graph.size(), graph, vector<bool>(graph.size(), false), 8, a);
    cout << solve(initial, 0).best << endl;

    return 0;
}
