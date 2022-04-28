#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/vector.hpp>
#include <cstdint>
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <deque>
#include <chrono>

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
    vector<vector<u32>>* graph;
    vector<u8> assignment;
    u32 lower_bound;
    bool valid;
    u32 weight;
    u32 n_set;

    explicit Sln(vector<vector<u32>>* graph)
    : graph(graph)
    , assignment(graph->size(), false)
    , lower_bound(0)
    , valid(false)
    , weight(0)
    , n_set(0)
    {}

    explicit Sln() : Sln(&dummy) {}

    Sln(const Sln& other)
    : graph(other.graph)
    , assignment(other.assignment)
    , lower_bound(other.lower_bound)
    , valid(other.valid)
    , weight(other.weight)
    , n_set(other.n_set)
    {}

    Sln(Sln&& other)
    : graph(other.graph)
    , assignment(move(other.assignment))
    , lower_bound(other.lower_bound)
    , valid(other.valid)
    , weight(other.weight)
    , n_set(other.n_set)
    {}

    Sln& operator=(const Sln& other) {
        graph = other.graph;
        assignment = other.assignment;
        lower_bound = other.lower_bound;
        valid = other.valid;
        weight = other.weight;
        n_set = other.n_set;
        return *this;
    }

    Sln& operator=(Sln&& other) {
        graph = other.graph;
        assignment = move(other.assignment);
        lower_bound = other.lower_bound;
        valid = other.valid;
        weight = other.weight;
        n_set = other.n_set;
        return *this;
    }

    Sln include(u32 node) const { return this->update(node, true); }
    Sln exclude(u32 node) const { return this->update(node, false); }

    Sln update(u32 node, bool included) const {
        //cppcheck-suppress shadowVariable
        const auto& graph = *(this->graph);
        Sln result = *this;
        result.n_set += included - assignment[node];
        result.assignment[node] = included;

        for (u32 i = 0; i < node; i++) {
            if (result.assignment[node] != result.assignment[i]) {
                result.weight += graph[node][i];
            }
        }

        result.lower_bound += result.weight;
        result.lower_bound -= this->weight;
        result.update_lower_bound(node);
        return result;
    }

    // relies on invariants:
    // - the assignment vector is decomposable into an assigned prefix and an unassigned suffix
    // - the about-to-be-assigned node is passed as the sole argument
    // - the graph matrix is symmetric (the graph is undirected)
    // [ 0 1 2 3 4 5 6 7 8 9 ]
    // [ 1 0 1 0 0 0 0 0 0 0 ]
    //    node ^
    u32 compute_lower_bound(u32 node) const {
        //cppcheck-suppress shadowVariable
        const auto& graph = *(this->graph);
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

    // perform only one iteration of the above algorithm
    // for incremental updates
    void update_lower_bound(u32 node) {
        // TODO ensure that it's called with node + 1
        //cppcheck-suppress shadowVariable
        const auto& graph = *(this->graph);
        u32 included = 0;
        u32 excluded = 0;

        for (u32 j = 0; j < node; j++) {
            if (assignment[j]) {
                included += graph[node][j];
            } else {
                excluded += graph[node][j];
            }
        }

        lower_bound += min(included, excluded);
    }
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int __attribute__((__unused__)) version)
    {
        ar & assignment;
        ar & lower_bound;
        ar & valid;
        ar & weight;
        ar & n_set;
    }
};

Sln optimum = Sln(&dummy);

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
    || (best.valid && sln.compute_lower_bound(node) > best.weight)
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

deque<pair<u16, Sln>> gen_initial_configurations(u16 a, vector<vector<u32>>& graph) {
    const u16 p = 3 * a;
    constexpr u8 z = 15; // found optimal with Hyperfine
    auto q = deque<pair<u16, Sln>>();
    q.emplace_back(make_pair(0, Sln(&graph)));

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

constexpr u32 TAG_TERMINATE = 0xdead;
constexpr u32 TAG_GRAPH = 0xface;
constexpr u32 TAG_WORK = 0xbeef;
constexpr u32 TAG_DONE = 0xcafe;
// 64 configurations at a time to each slave
constexpr u32 batch_size = 64;

struct WorkData {
    vector<u16> node_ids;
    vector<Sln> configs;
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int __attribute__((__unused__)) version)
    {
        ar & node_ids;
        ar & configs;
    }
};

// invariant: only called when the queue is non-empty
mpi::request __attribute__((warn_unused_result)) send_configurations(
    mpi::communicator& world, i32 dest, deque<pair<u16, Sln>>& queue
) {
    WorkData batch;
    for (u32 i = 0; i < batch_size; i++) {
        if (queue.empty()) {
            break;
        }
        batch.node_ids.emplace_back(queue.front().first);
        batch.configs.emplace_back(queue.front().second);
        queue.pop_front();
    }
    return world.isend(dest, TAG_WORK, batch);
}

deque<pair<u16, Sln>> init_slaves(mpi::communicator& world, u16 a, vector<vector<u32>>& graph) {
    auto queue = gen_initial_configurations(a, graph);
    vector<mpi::request> requests;
    requests.reserve(2 * world.size());
    for (i32 dest = 1; dest < world.size(); dest++) {
        // send the graph and the initial configurations to the workers
        requests.emplace_back(world.isend(dest, TAG_GRAPH, graph));
        requests.emplace_back(send_configurations(world, dest, queue));
    }
    // wait for all the data to send
    mpi::wait_all(requests.begin(), requests.end());
    return queue;
}

void master(mpi::communicator& world, u16 a) {
    auto start = chrono::steady_clock::now();
    auto graph = load_input();
    chrono::duration<double> input_load_time = chrono::steady_clock::now() - start;
    cout << "loaded graph of n = " << graph.size() << " in " << input_load_time.count() << " s" << endl;

    start = chrono::steady_clock::now();
    auto queue = init_slaves(world, a, graph);

    Sln best = Sln(&graph);
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
            send_configurations(world, status.source(), queue).wait();
        }
        best = pick(best, sln);
    }

    cout << best.weight << " (" << best.n_set << ")" << endl;
    for (u32 i = 0; i < graph.size(); i++) {
        cout << (u32)best.assignment[i] << " ";
    }

    chrono::duration<double> total_time = chrono::steady_clock::now() - start;
    cout << endl << "took " << total_time.count() << " s" << endl;
}

// Listen for TAG_WORK_NODE_VECTOR, TAG_WORK_SLN_VECTOR, and TAG_TERMINATE
// asynchronously. Only move on if both NODE_VECTOR and SLN_VECTOR are received, or
// if only TAG_TERMINATE is received.
bool wait_for_work(mpi::communicator& world, WorkData& work) {
    mpi::request reqs[] = {
        world.irecv(0, TAG_TERMINATE),
        world.irecv(0, TAG_WORK, work),
    };
    cout << "slave " << world.rank() << ": waiting for work" << endl;
    auto [status, req_ptr] = mpi::wait_any(reqs, reqs + 2);

    if (req_ptr == &reqs[0]) {
        cout << "slave " << world.rank() << ": terminating" << endl;
        reqs[1].cancel();
        return false;
    } else {
        cout << "slave " << world.rank() << ": received work" << endl;
        reqs[0].cancel();
        return true;
    }
}

void slave(mpi::communicator& world, u16 a, u16 n_threads) {
    // receive the graph
    vector<vector<u32>> graph;
    world.recv(0, TAG_GRAPH, graph);

    for (WorkData work; wait_for_work(world, work);) {
        // process the batch
        #pragma omp parallel for schedule(dynamic) num_threads(n_threads)
        for (u32 i = 0; i < work.node_ids.size(); i++) {
            auto node = work.node_ids[i];
            auto sln = work.configs[i];
            sln.graph = &graph;
            Sln result(sln); // we read the actual result from the optimum
            solve(&result, a, graph, sln, sln, node);
        }

        // send TAG_DONE to the master
        world.send(0, TAG_DONE, optimum);
    }
}

int main(int argc, char** argv) {
    mpi::environment env;
    mpi::communicator world;

    if (argc < 3) {
        cout << "usage: " << argv[0] << " <a> <n_threads>" << endl;
        return 1;
    }

    u16 a = atoi(argv[1]);
    u16 n_threads = atoi(argv[2]);

    if (world.rank() == 0) {
        master(world, a);
        cout << "master " << world.rank() << " finished" << endl;
    } else {
        slave(world, a, n_threads);
        cout << "slave " << world.rank() << " finished" << endl;
    }

    return 0;
}
