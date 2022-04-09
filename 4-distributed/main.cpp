#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/vector.hpp>
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

auto dummy = vector<vector<u32>>();
constexpr u32 a = 10;

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
    vector<vector<u32>>& graph;
    vector<bool> assignment;
    bool valid;
    u32 weight;
    u32 n_set;

    explicit Sln(vector<vector<u32>>& graph)
    : graph(graph)
    , assignment(graph.size(), false)
    , valid(false)
    , weight(0)
    , n_set(0)
    {}

    explicit Sln() : Sln(dummy) {}

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
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int __attribute__((__unused__)) version)
    {
        ar & assignment;
        ar & valid;
        ar & weight;
        ar & n_set;
    }
};

Sln optimum = Sln(dummy);

const Sln& pick(const Sln& sln_x, const Sln& sln_y) {
    if (!sln_x.valid) {
        return sln_y;
    }
    return (!sln_y.valid || sln_x.weight < sln_y.weight) ? sln_x : sln_y;
}

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

deque<pair<u16, Sln>> gen_initial_configurations(vector<vector<u32>>& graph) {
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

void solve(u32 a, vector<vector<u32>>& graph) {
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

constexpr u32 TAG_WORK_NODE_VECTOR = 0xbeef;
constexpr u32 TAG_WORK_SLN_VECTOR = 0xbecca;
constexpr u32 TAG_TERMINATE = 0xdead;
constexpr u32 TAG_GRAPH = 0xface;
constexpr u32 TAG_DONE = 0xcafe;
// 64 configurations at a time to each slave
constexpr u32 batch_size = 64;

// invariant: only called when the queue is non-empty
pair<mpi::request, mpi::request> send_configurations(mpi::communicator& world, i32 dest, deque<pair<u16, Sln>>& queue) {
    auto batch = make_pair(vector<u16>(), vector<Sln>());
    for (u32 i = 0; i < batch_size; i++) {
        if (queue.empty()) {
            break;
        }
        batch.first.emplace_back(queue.front().first);
        batch.second.emplace_back(queue.front().second);
        queue.pop_front();
    }
    return make_pair(
        world.isend(dest, TAG_WORK_NODE_VECTOR, batch.first),
        world.isend(dest, TAG_WORK_SLN_VECTOR, batch.second)
    );
}

deque<pair<u16, Sln>> init_slaves(mpi::communicator& world, vector<vector<u32>>& graph) {
    auto queue = gen_initial_configurations(graph);
    vector<mpi::request> requests;
    requests.reserve(3 * world.size());
    for (i32 dest = 1; dest < world.size(); dest++) {
        // send the graph and the initial configurations to the workers
        requests.emplace_back(world.isend(dest, TAG_GRAPH, graph));
        auto p = send_configurations(world, dest, queue);
        requests.emplace_back(p.first);
        requests.emplace_back(p.second);
    }
    return queue;
}

void master(mpi::communicator& world) {
    auto graph = load_input();
    auto queue = init_slaves(world, graph);

    Sln best = Sln(graph);
    u32 working_slaves = world.size() - 1;
    while (working_slaves > 0) {
        // listen for TAG_DONE from a slave
        Sln sln = best;
        auto status = world.recv(mpi::any_source, TAG_DONE, sln);
        if (queue.empty()) {
            // all work is done, terminate the slave
            world.send(status.source(), TAG_TERMINATE);
            working_slaves--;
        } else {
            // send the next batch of configurations to the slave
            send_configurations(world, status.source(), queue);
        }
        best = pick(best, sln);
    }

    cout << best.weight << " (" << best.n_set << ")" << endl;
    for (u32 i = 0; i < graph.size(); i++) {
        cout << best.assignment[i] << " ";
    }
    cout << endl;
}

void slave(mpi::communicator& world) {
    // receive the graph
    vector<vector<u32>> graph;
    world.recv(0, TAG_GRAPH, graph);

    while (true) {
        cout << "slave " << world.rank() << ": waiting for work" << endl;
        vector<u16> node_vector;
        vector<Sln> sln_vector;

        // listen for TAG_WORK_NODE_VECTOR, TAG_WORK_SLN_VECTOR, and TAG_TERMINATE
        // asynchronously. Only move on if both NODE_VECTOR and SLN_VECTOR are received, or
        // if only TAG_TERMINATE is received.
        mpi::request reqs[] = {
            world.irecv(0, TAG_TERMINATE),
            world.irecv(0, TAG_WORK_NODE_VECTOR, node_vector),
            world.irecv(0, TAG_WORK_SLN_VECTOR, sln_vector)
        };
        auto [status, req_ptr] = mpi::wait_any(reqs, reqs + 3);

        if (req_ptr == &reqs[0]) {
            cout << "slave " << world.rank() << ": terminating" << endl;
            reqs[1].cancel();
            reqs[2].cancel();
            break;
        } else {
            cout << "slave " << world.rank() << ": waiting for the remaining vector" << endl;
            mpi::wait_all(reqs + 1, reqs + 3);
            reqs[0].cancel();
        }

        // process the batch
        #pragma omp parallel for schedule(dynamic)
        for (u32 i = 0; i < node_vector.size(); i++) {
            auto node = node_vector[i];
            auto sln = sln_vector[i];
            // FIXME shitty semantics lead to a copy
            sln.graph = graph;
            Sln result(sln); // we read the actual result from the optimum
            solve(&result, a, graph, sln, sln, node);
        }

        // send TAG_DONE to the master
        world.send(0, TAG_DONE, optimum);
    }
}

int main() {
    mpi::environment env;
    mpi::communicator world;

    if (world.rank() == 0) {
        master(world);
        cout << "master " << world.rank() << " finished" << endl;
    } else {
        slave(world);
        cout << "slave " << world.rank() << " finished" << endl;
    }

    return 0;
}
