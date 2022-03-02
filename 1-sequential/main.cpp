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

u32 cut_weight(const vector<vector<u32>>& graph, const vector<bool>& assignment) {
    const u32 a = 10;
    u32 w = 0;
    u32 n_set = 0;
    for (u32 i = 0; i < graph.size(); i++) {
        n_set += assignment[i];
        for (u32 j = 0; j < i; j++) {
            if (assignment[i] != assignment[j]) {
                w += graph[i][j];
            }
        }
    }

    if ((n_set == a) || (n_set == graph.size() - a)) {
        return w;
    }
    return (u32)-1;
}

vector<bool> solve(const vector<vector<u32>>& graph, const vector<bool>& assignment, u32 node) {
    if (node > graph.size()) {
        return assignment;
    }

    auto included = assignment;
    included[node] = true;
    auto excluded = assignment;
    excluded[node] = false;

    auto inc_sln = solve(graph, included, node + 1);
    auto exc_sln = solve(graph, excluded, node + 1);

    if (cut_weight(graph, inc_sln) < cut_weight(graph, exc_sln)) {
        return inc_sln;
    }
    return exc_sln;
}

int main() {
    auto graph = load_input();

    auto sln = solve(graph, vector<bool>(graph.size(), false), 0);
    cout << cut_weight(graph, sln) << endl;

    return 0;
}
