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
    u32 weight;
    u32 n_set;

    explicit Sln(const vector<vector<u32>>& graph)
    : graph(graph)
    , assignment(graph.size(), false)
    , weight(0)
    , n_set(0)
    {}

    Sln(const Sln& other)
    : graph(other.graph)
    , assignment(other.assignment)
    , weight(other.weight)
    , n_set(other.n_set)
    {}

    Sln(Sln&& other)
    : graph(other.graph)
    , assignment(move(other.assignment))
    , weight(other.weight)
    , n_set(other.n_set)
    {}

    Sln& operator=(const Sln& other) {
        assignment = other.assignment;
        weight = other.weight;
        n_set = other.n_set;
        return *this;
    }

    Sln& operator=(Sln&& other) {
        assignment = move(other.assignment);
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

    u32 lower_bound(u32 node) const {
        u32 result = this->weight;
        for (u32 i = node; i < graph.size(); i++) {
            u32 included = 0;
            u32 excluded = 0;
            for (u32 j = 0; j < graph.size(); j++) {
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

    bool valid(u32 a) const {
        return (n_set == a) || (n_set == graph.size() - a);
    }
};

const Sln& pick(u32 a, const Sln& sln_x, const Sln& sln_y) {
    if (!sln_x.valid(a)) {
        return sln_y;
    }
    return (!sln_y.valid(a) || sln_x.weight < sln_y.weight) ? sln_x : sln_y;
}

Sln solve(u32 a, const vector<vector<u32>>& graph, const Sln& sln, const Sln& best, u32 node) {
    if (false
    || (node >= graph.size())
    || (min(sln.n_set, graph.size() - sln.n_set) > a)
    || (best.valid(a) && sln.weight > best.weight)
    ) { return sln; }

    const auto incl_sln = solve(a, graph, sln.include(node), best, node + 1);
    const auto new_best = pick(a, incl_sln, best);
    const auto excl_sln = solve(a, graph, sln.exclude(node), new_best, node + 1);
    return pick(a, excl_sln, new_best);
}

int main() {
    auto graph = load_input();
    u32 a = 10;
    auto init = Sln(graph);
    auto sln = solve(a, graph, init, init, 0);
    cout << sln.weight << " (" << sln.n_set << ")" << endl;
    for (u32 i = 0; i < graph.size(); i++) {
        if (sln.assignment[i]) {
            cout << sln.assignment[i] << " ";
        }
    }
    cout << endl;

    return 0;
}
