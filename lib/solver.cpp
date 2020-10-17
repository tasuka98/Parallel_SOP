#include "solver.hpp"
#include "hash.hpp"
#include <iostream>
#include <fstream>
#include <stack>
#include <sstream>
#include <chrono>
#include <limits>
#include <algorithm>
#include <cstring>
#include <numeric>
#include <functional>
#include <vector>
#include <thread>
#include <queue>
#include <mutex>
#include <bits/stdc++.h>
#include <unordered_map>
#include <condition_variable>

using namespace std;

static Hash_Map history_table(1237935);
//static int recently_added = 0;

static int best_cost = 0;
static vector<int> best_solution;
static queue<solver> GPQ;

//Variable for locking;
static mutex GPQ_lock, Sol_lock, Split_lock;
static mutex Split_Call, asssign_mutex, unlock_mutex;
static mutex thread_load_mutex;
static condition_variable Idel;

//Protected by lock
atomic<int> active_thread (0);
atomic<int> idle_counter (0);
atomic<bool> time_out (false);
int* thread_load;

//Config variable (Reading Only)
static string enum_option;
static int global_pool_size = 0;
static int local_pool_size = 0;
static int t_limit = 0;
static int thread_total = 0;
static int local_depth = 0;
static string Assign_Opt;
static string Steal_Opt;

//Shared resources
std::chrono::time_point<std::chrono::system_clock> start_time_limit;
static vector<vector<int>> dependent_graph;
static vector<vector<edge>> in_degree;
static vector<vector<edge>> hung_graph;
static vector<vector<edge>> cost_graph;

void solver::assign_parameter(vector<string> setting) {
    enum_option = setting[0];
    //cout << enum_option << endl;
    t_limit = atoi(setting[1].c_str());
    //cout << t_limit << endl;
    global_pool_size = atoi(setting[2].c_str());
    //cout << global_pool_size << endl;
    local_pool_size = atoi(setting[3].c_str());
    //cout << local_pool_size << endl;
    thread_total = atoi(setting[4].c_str());
    //cout << thread_total << endl;
    local_depth = atoi(setting[5].c_str());
    //cout << local_depth << endl;
    Assign_Opt = setting[6];
    //cout << Assign_Opt << endl;
    Steal_Opt = setting[7];
    //cout << Steal_Opt << endl;
    return;
}

int solver::dynamic_hungarian(int src, int dest) {
    hungarian_solver.fix_row(src, dest);
	hungarian_solver.fix_column(dest, src);
	hungarian_solver.solve_dynamic();
    return hungarian_solver.get_matching_cost()/2;
}

bool solver::HistoryUtilization(int* lowerbound,bool* found) {
    string bit_string(node_count, '0');

    for (auto node : cur_solution) {
        bit_string[node] = '1';
    }   

    int last_element = cur_solution.back();
    auto key = make_pair(bit_string,last_element);

    HistoryNode history_node = history_table.retrieve(key);
    if (history_node.prefix_cost == -1) return true;

    *found = true;
    int history_prefix = history_node.prefix_cost;
    int history_lb = history_node.lower_bound;
    *lowerbound = history_lb;

    if (cur_cost >= history_prefix) return false;

    int imp = history_prefix - cur_cost;
    
    if (imp <= history_lb - best_cost) return false;
    history_node.prefix_cost = cur_cost;
    history_node.lower_bound = history_lb - imp;
    *lowerbound = history_node.lower_bound;
    
    history_table.insert(key,history_node);
    return true;
}

bool bound_sort(const node& src,const node& dest) {
    if (src.lb == dest.lb) return src.nc > dest.nc;
    return src.lb > dest.lb;
}

bool nearest_sort(const node& src,const node& dest) {
    if (src.nc == dest.nc) return src.lb > dest.lb;
    return src.nc > dest.nc;
}

//Get new subproblem state and push it back to the pool;
void solver::assign_workload(int taken_n, int lb) {
    int taken_node = taken_n;
    int cur_node = cur_solution.back();
    taken_arr[taken_node] = 1;
    for (int vertex : dependent_graph[taken_node]) depCnt[vertex]--;
    cur_cost += cost_graph[cur_node][taken_node].weight;
    cur_solution.push_back(taken_node);
    hungarian_solver.fix_row(cur_node, taken_node);
    hungarian_solver.fix_column(taken_node, cur_node);
    hungarian_solver.solve_dynamic();
    local_lb = lb;

    local_pool->push(*this);

    taken_arr[taken_node] = 0;
    for (int vertex : dependent_graph[taken_node]) depCnt[vertex]++;
    cur_cost -= cost_graph[cur_node][taken_node].weight;
    cur_solution.pop_back();
    hungarian_solver.undue_row(cur_node, taken_node);
    hungarian_solver.undue_column(taken_node, cur_node);
    local_lb = -1;

    return;
}

void solver::untake_node(int src, int dest, int taken_node) {
    for (int vertex : dependent_graph[taken_node]) depCnt[vertex]++;
    taken_arr[taken_node] = 0;
    cur_solution.pop_back();
    cur_cost -= cost_graph[src][dest].weight;
    if (cur_solution.size() >= 2) {
        hungarian_solver.undue_row(src,dest);
        hungarian_solver.undue_column(dest,src);
    }
    return;
}

void solver::take_node(int taken_node) {
    for (int vertex : dependent_graph[taken_node]) depCnt[vertex]--;
    cur_solution.push_back(taken_node);
    taken_arr[taken_node] = 1;
    return;
}

void solver::push_to_global_pool() {
    total_lb -= local_pool->front().local_lb;
    GPQ.push(local_pool->front());
    local_pool->pop();
    return;
}

void solver::notify_finished() {
    active_thread--;
    std::unique_lock<std::mutex> idel_lck(Split_lock);
    Idel.notify_all();
    idel_lck.unlock();
    return;
}

void solver::check_workload_request(int i) {
    if (idle_counter > 0) {
        asssign_mutex.lock();
        //cout << "visiting with thread id = " << std::this_thread::get_id() << "and i = " << i <<  endl;
        if (idle_counter > 0) {
            // Make sure we have at least two children in the ready list before splitting.
            bool push_to_global = take_from_local();
            if (push_to_global) {
                if (Assign_Opt == "SINGLE") {
                    idle_counter--;
                    std::unique_lock<std::mutex> idel_lck(Split_lock);
                    Idel.notify_one();
                    idel_lck.unlock();
                }
                else if (Assign_Opt == "FULL") {
                    GPQ_lock.lock();
                    for (unsigned i = 0; i < GPQ.size(); i++) {
                        idle_counter--;
                        std::unique_lock<std::mutex> idel_lck(Split_lock);
                        Idel.notify_one();
                        idel_lck.unlock();
                        if (idle_counter == 0) break;
                    }
                    GPQ_lock.unlock();
                }
            }
        }
        asssign_mutex.unlock();
    }
    return;
}

bool solver::check_tlimit() {
    auto cur_time = std::chrono::system_clock::now();
    if (std::chrono::duration<double>(cur_time - start_time_limit).count() > t_limit) {
        time_out = true;
        active_thread = 0;
        std::unique_lock<std::mutex> idel_lck(Split_lock);
        Idel.notify_all();
        idel_lck.unlock();
        return true;
    }
    return false;
}

bool solver::take_from_local() {
    if (local_pool->size() > 0) {
        //Find minimum load in the load array;
        int min = INT_MAX;
        bool steal = false;
        thread_load_mutex.lock();
        for (int k = 0; k < thread_total; k++) {
            if (thread_load[k] < min) min = thread_load[k];
        }
        if (min == total_lb / (int)local_pool->size()) steal = true;
        thread_load_mutex.unlock();

        if (Assign_Opt == "SINGLE" && steal) {
            push_to_global_pool();
            thread_load_mutex.lock();
            if (!local_pool->empty()) thread_load[thread_id] = total_lb / local_pool->size();
            else thread_load[thread_id] = INT_MAX;
            thread_load_mutex.unlock();
            return true;
        }
        else if (Assign_Opt == "FULL" && steal) {
            int p_count = 0;
            while (local_pool->size() > 1 && p_count < global_pool_size) {
                push_to_global_pool();
                p_count++;
            }
            thread_load_mutex.lock();
            if (!local_pool->empty()) thread_load[thread_id] = total_lb / local_pool->size();
            else thread_load[thread_id] = INT_MAX;
            thread_load_mutex.unlock();
            return true;
        }
    }

    return false;
}

void solver::enumerate(int i) {
    if (time_out) return;

    vector<node> ready_list;

    for (int i = node_count-1; i >= 0; i--) {
        if (!depCnt[i] && !taken_arr[i]) {
            //Push vertices with 0 depCnt into the ready list
            ready_list.push_back(node(i,-1));
        }
    }
    //update vertex dependent count for current node neighbors

    int taken_node = 0;
    int u = 0;
    int v = 0;
    int total_lb = 0;

    if (!cur_solution.empty()) {
        for (int i = 0; i < (int)ready_list.size(); i++) {
            node dest = ready_list[i];
            int src = cur_solution.back();
            cur_solution.push_back(dest.n);
            cur_cost += cost_graph[src][dest.n].weight;
            int temp_lb = -1;
            bool taken = false;
            
            //Backtrack
            if (cur_cost >= best_cost) {
                cur_solution.pop_back();
                cur_cost -= cost_graph[src][dest.n].weight;
                ready_list.erase(ready_list.begin()+i);
                i--;
                continue;
            }
            if (cur_solution.size() == (size_t)node_count) {
                if (cur_cost < best_cost) {
                    Sol_lock.lock();
                    best_solution = cur_solution;
                    best_cost = cur_cost;
                    Sol_lock.unlock();
                    cur_solution.pop_back();
                    cur_cost -= cost_graph[src][dest.n].weight;
                    ready_list.erase(ready_list.begin()+i);
                    i--;
                    continue;
                }
            }

            else {
                bool decision = HistoryUtilization(&temp_lb,&taken);
                if (!taken) {
                    string bit_string(node_count, '0');
                    for (auto node : cur_solution) bit_string[node] = '1';
                    int last_element = cur_solution.back();
                    auto key = make_pair(bit_string,last_element);
                    temp_lb = dynamic_hungarian(src,dest.n);

                    if (history_table.get_cur_size() < 0.8 * history_table.get_max_size()) {
                        history_table.insert(key,HistoryNode(cur_cost,temp_lb));
                    }
                    else if (i < int(0.5 * node_count)) {
                        history_table.insert(key,HistoryNode(cur_cost,temp_lb));
                    }

                    hungarian_solver.undue_row(src,dest.n);
                    hungarian_solver.undue_column(dest.n,src);
                }
                else if (taken && !decision) {
                    cur_solution.pop_back();
                    cur_cost -= cost_graph[src][dest.n].weight;
                    ready_list.erase(ready_list.begin()+i);
                    i--;
                    continue;
                }
                if (temp_lb >= best_cost) {
                    cur_solution.pop_back();
                    cur_cost -= cost_graph[src][dest.n].weight;
                    ready_list.erase(ready_list.begin()+i);
                    i--;
                    continue;
                }
                cur_solution.pop_back();
                cur_cost -= cost_graph[src][dest.n].weight;
                ready_list[i].nc = cost_graph[src][dest.n].weight;
                ready_list[i].lb = temp_lb;
                total_lb += temp_lb;
            }
        }
        if (enum_option == "DH") sort(ready_list.begin(),ready_list.end(),bound_sort);
        else if (enum_option == "NN") sort(ready_list.begin(),ready_list.end(),nearest_sort);
    }

    if (idle_counter > 0) {
        asssign_mutex.lock();
        //cout << "visiting with thread id = " << std::this_thread::get_id() << "and i = " << i <<  endl;
        if (idle_counter > 0) {
            // Make sure we have at least two children in the ready list before splitting.
            bool push_to_global = false;

            if (local_pool->size() > 0) {
                //Find minimum load in the load array;
                int min = INT_MAX;
                bool steal = false;
                thread_load_mutex.lock();
                for (int k = 0; k < thread_total; k++) {
                    if (thread_load[k] < min) min = thread_load[k];
                }
                if (min == total_lb / (int)local_pool->size()) steal = true;
                thread_load_mutex.unlock();

                if (Assign_Opt == "SINGLE" && steal) {
                    push_to_global_pool();
                    thread_load_mutex.lock();
                    if (!local_pool->empty()) thread_load[thread_id] = total_lb / local_pool->size();
                    else thread_load[thread_id] = INT_MAX;
                    thread_load_mutex.unlock();
                    push_to_global = true;
                }
                else if (Assign_Opt == "FULL" && steal) {
                    int p_count = 0;
                    while (local_pool->size() > 1 && p_count < global_pool_size) {
                        push_to_global_pool();
                        p_count++;
                    }
                    thread_load_mutex.lock();
                    if (!local_pool->empty()) thread_load[thread_id] = total_lb / local_pool->size();
                    else thread_load[thread_id] = INT_MAX;
                    thread_load_mutex.unlock();
                    push_to_global = true;
                }
            }

            if (push_to_global) {
                if (Assign_Opt == "SINGLE") {
                    idle_counter--;
                    std::unique_lock<std::mutex> idel_lck(Split_lock);
                    Idel.notify_one();
                    idel_lck.unlock();
                }
                else if (Assign_Opt == "FULL") {
                    GPQ_lock.lock();
                    for (unsigned i = 0; i < GPQ.size(); i++) {
                        idle_counter--;
                        std::unique_lock<std::mutex> idel_lck(Split_lock);
                        Idel.notify_one();
                        idel_lck.unlock();
                        if (idle_counter == 0) break;
                    }
                    GPQ_lock.unlock();
                }
            }
        }
        asssign_mutex.unlock();
    }

    while(!ready_list.empty()) {
        //Take the choosen node and back track if ready list is empty;
        taken_node = ready_list.back().n;
        ready_list.pop_back();
        if (!cur_solution.empty()) {
            u = cur_solution.back();
            v = taken_node;
            hungarian_solver.fix_row(u, v);
            hungarian_solver.fix_column(v, u);
            cur_cost += cost_graph[u][v].weight;
        }

        for (int vertex : dependent_graph[taken_node]) depCnt[vertex]--;
        cur_solution.push_back(taken_node);
        taken_arr[taken_node] = 1;
        
        //Push to local pool until it is full;
        if (local_pool->size() < (size_t)thread_total && ready_list.size() > 1 && i <= (float(local_depth) / float(100) * node_count)) {
            while (ready_list.size() > 1 && local_pool->size() < (size_t)local_pool_size) {
                assign_workload(ready_list.back().n,ready_list.back().lb);
                total_lb += ready_list.back().lb;
                ready_list.pop_back();
            }

            thread_load_mutex.lock();
            if (!local_pool->empty()) thread_load[thread_id] = total_lb / local_pool->size();
            else thread_load[thread_id] = INT_MAX;
            thread_load_mutex.unlock();
        }

        enumerate(i+1);
        
        for (int vertex : dependent_graph[taken_node]) depCnt[vertex]++;
        taken_arr[taken_node] = 0;
        cur_solution.pop_back();
        cur_cost -= cost_graph[u][v].weight;
        if (cur_solution.size() >= 2) {
            hungarian_solver.undue_row(u,v);
            hungarian_solver.undue_column(v,u);
        }

        auto cur_time = std::chrono::system_clock::now();
        if (std::chrono::duration<double>(cur_time - start_time_limit).count() > t_limit) {
            time_out = true;
            active_thread = 0;
            std::unique_lock<std::mutex> idel_lck(Split_lock);
            Idel.notify_all();
            idel_lck.unlock();
            return;
        }
    }
    

    bool terminate = true;

    if (i == initial_depth && active_thread >= 1) {
        notify_finished();
        if (!local_pool->empty()) {
            *this = local_pool->front();
            local_pool->pop();
            terminate = false;
        }
        else {
            GPQ_lock.lock();
            if (!GPQ.empty()) {
                GPQ.front().thread_id = thread_id;
                auto current_local_pool = local_pool;
                *this = GPQ.front();
                local_pool = current_local_pool;
                GPQ.pop();
                terminate = false;
            }
            GPQ_lock.unlock();

            if (active_thread > 0 && terminate) {
                idle_counter++;
                GPQ_lock.lock();
                int size = GPQ.size();
                GPQ_lock.unlock();

                std::unique_lock<std::mutex> idel_lck(Split_lock);
                while (size == 0 && active_thread > 0) {
                    Idel.wait(idel_lck);
                    GPQ_lock.lock();
                    size = GPQ.size();
                    GPQ_lock.unlock();
                }
                GPQ_lock.lock();
                size = GPQ.size();
                if (size != 0) {
                    GPQ.front().thread_id = thread_id;
                    auto current_local_pool = local_pool;
                    *this = GPQ.front();
                    local_pool = current_local_pool;
                    GPQ.pop();
                    terminate = false;
                }
                GPQ_lock.unlock();
            }
        }
    }

    if (!terminate) {
        active_thread++;
        total_lb = 0;
        enumerate(initial_depth);
    }

    return;
}

void solver::solve_parallel(int thread_num, int pool_size) {
    start_time_limit = std::chrono::system_clock::now();
    vector<solver> solvers(thread_num);
    vector<thread> Thread_manager(thread_num);
    vector<int> ready_thread;
    //Initially fill the GPQ with solvers:
    vector<node> ready_list;

    cur_solution.push_back(0);
    taken_arr[0] = 1;
    for (int vertex : dependent_graph[0]) depCnt[vertex]--;

    for (int i = node_count-1; i >= 0; i--) {
        if (!depCnt[i] && !taken_arr[i]) {
            //Push vertices with 0 depCnt into the ready list
            ready_list.push_back(node(i,-1));
        }
    }

    thread_load = new int [thread_num];

    //Initial filling of the GPQ
    for (auto node : ready_list) {
        solver target = *this;
        int taken_node = node.n;
        int cur_node = target.cur_solution.back();
        target.taken_arr[taken_node] = 1;
        for (int vertex : dependent_graph[taken_node]) target.depCnt[vertex]--;
        target.cur_cost += cost_graph[cur_node][taken_node].weight;
        target.cur_solution.push_back(taken_node);
        target.hungarian_solver.fix_row(cur_node, taken_node);
        target.hungarian_solver.fix_column(taken_node, cur_node);
        target.hungarian_solver.solve_dynamic();
        GPQ.push(target);
    }

    for (int i = thread_num - 1; i >= 0; i--) ready_thread.push_back(i);
    //While GPQ is not empty do split operation or assign threads with new node.

    while (true) {
        GPQ_lock.lock();
        if (GPQ.empty()) {
            GPQ_lock.unlock();
            break;
        }

        while (GPQ.size() < (size_t)pool_size) {
            auto target = GPQ.front();
            GPQ.pop();
            ready_list.clear();
            for (int i = node_count-1; i >= 0; i--) {
                if (!target.depCnt[i] && !target.taken_arr[i]) ready_list.push_back(node(i,-1));
            }

            if (!ready_list.empty()) {
                for (auto node : ready_list) {
                int vertex = node.n;
                //Push split node back into GPQ
                    if (!target.depCnt[vertex] && !target.taken_arr[vertex]) {
                        int taken_node = vertex;
                        int cur_node = target.cur_solution.back();
                        target.taken_arr[taken_node] = 1;
                        for (int vertex : dependent_graph[taken_node]) target.depCnt[vertex]--;
                        target.cur_cost += cost_graph[cur_node][taken_node].weight;
                        target.cur_solution.push_back(taken_node);
                        if (cur_solution.size() == (size_t)node_count && cur_cost < best_cost) {
                            Sol_lock.lock();
                            best_solution = cur_solution;
                            best_cost = cur_cost;
                            Sol_lock.unlock();
                        }
                        target.hungarian_solver.fix_row(cur_node, taken_node);
                        target.hungarian_solver.fix_column(taken_node, cur_node);
                        target.hungarian_solver.solve_dynamic();
                        
                        GPQ.push(target);

                        target.taken_arr[taken_node] = 0;
                        for (int vertex : dependent_graph[taken_node]) target.depCnt[vertex]++;
                        target.cur_cost -= cost_graph[cur_node][taken_node].weight;
                        target.cur_solution.pop_back();
                        target.hungarian_solver.undue_row(cur_node, taken_node);
                        target.hungarian_solver.undue_column(taken_node, cur_node);
                    }
                }
            }
            else break;
        }
        
        for (int i = 0; i < thread_num; i++) {
            if (!GPQ.empty()) {
                auto subproblem = GPQ.front();
                GPQ.pop();
                int k = subproblem.cur_solution.size() - 1;
                solvers[i] = subproblem;
                solvers[i].initial_depth = k;
                solvers[i].thread_id = active_thread;
                solvers[i].local_pool = new queue<solver>();
                active_thread++;
                Thread_manager[i] = thread(&solver::enumerate,move(solvers[i]),k);
            }
        }
        
        GPQ_lock.unlock();
        
        for (int i = 0; i < thread_num; i++) {
            if (Thread_manager[i].joinable()) {
                //cout << "waiting for thread to be joined" << endl;
                Thread_manager[i].join();
                delete solvers[i].local_pool;
                ready_thread.push_back(i);
            }
        }
        active_thread = 0;

        if (time_out) {
            return;
        }
    }
    
    delete thread_load;
    return;
}

void solver::solve(string filename) {
    retrieve_input(filename);
    //Remove redundant edges in the cost graph
    transitive_redundantcy();
    history_table.set_node_t(node_count);
    best_solution = nearest_neightbor();
    int max_edge_weight = get_maxedgeweight();
    hungarian_solver = Hungarian(node_count, max_edge_weight+1, get_cost_matrix(max_edge_weight+1));
    MMCP_static_lowerbound = hungarian_solver.start()/2;
    //picked_list = vector<bool>(node_count,false);
    //EGB_static_lowerbound = mmcp_lb();
    depCnt = vector<int>(node_count,0);
    taken_arr = vector<int>(node_count,0);

    for (int i = 0; i < node_count; i++) {
        for (unsigned k = 0; k < dependent_graph[i].size(); k++) {
            depCnt[dependent_graph[i][k]]++;
        }
    }
    
    /*
    cout << "best solution found using NN is " << best_cost << endl;
    cout << "the NN solution contains ";
    for (int i = 0; i < node_count; i++) {
        if (i != node_count - 1) cout << best_solution[i] << "-->";
        else cout << best_solution[i];
    }
    cout << endl;
    cout << "MMCP-based LB is " << MMCP_static_lowerbound << endl;
    */
    cur_cost = 0;
    auto start_time = chrono::high_resolution_clock::now();
    solve_parallel(thread_total,global_pool_size);
    auto end_time = chrono::high_resolution_clock::now();
    auto total_time = chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    cout << enum_option << ": " << best_cost << "," << setprecision(4) << total_time / (float)(1000000) << endl;

    /*
    cout << "the optimal solution contains ";
    for (int i = 0; i < node_count; i++) {
        if (i != node_count - 1) cout << best_solution[i] << "-->";
        else cout << best_solution[i];
    }
    cout << endl;
    */
    //cout << "Total calculated bounds are " << calculated_bounds << endl;
}

void solver::retrieve_input(string filename) {
    ifstream inFile;
    string line;
    inFile.open(filename);

    // Read input files and store it inside an array.
    vector<vector<int>> file_matrix;
    vector<int> matrix_temp;
    while (getline(inFile,line)) {  
        stringstream sstream;
        sstream << line;
        string weight;
        int weight_num;
        while (sstream >> weight) {
            stringstream(weight) >> weight_num;
            matrix_temp.push_back(weight_num);
        }
        file_matrix.push_back(matrix_temp);
        matrix_temp.clear();
    }

    unsigned size = file_matrix.size();
    
    cost_graph = vector<vector<edge>>(size);
    hung_graph = vector<vector<edge>>(size);
    dependent_graph = vector<vector<int>>(size);
    
    for (int i = 0; i < (int)file_matrix.size(); i++) {
        int j = 0;
        for (auto edge_weight: file_matrix[i]) {
            if (edge_weight < 0) {
                cost_graph[i].push_back(edge(i,j,file_matrix[j][i]));
                hung_graph[i].push_back(edge(i,j,-1));
                dependent_graph[j].push_back(i);
            }
            else { 
                cost_graph[i].push_back(edge(i,j,edge_weight));
                hung_graph[i].push_back(edge(i,j,edge_weight));
            }
            j++;
        }
    }
    node_count = cost_graph.size();

    //Trim redundant edges
    return;
}

void solver::transitive_redundantcy() {
    in_degree = std::vector<vector<edge>>(node_count);
    for (int i = 0; i < node_count; i++) {
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            int c = dependent_graph[i][k];
            in_degree[c].push_back(edge(i,c,hung_graph[i][c].weight));
        }
    }

    for(int i = 0; i < node_count; ++i) {
        vector<edge> preceding_nodes;
        for (int k = 0; k < (int)dependent_graph[i].size(); k++) preceding_nodes.push_back(edge(i,dependent_graph[i][k],-1));
        unordered_set<int> expanded_nodes;
        for(int j = 0; j < (int)preceding_nodes.size(); ++j){
            vector<edge> st;
            st.push_back(preceding_nodes[j]);
            while(!st.empty()){
                edge dependence_edge = st.back();
                st.pop_back();
                if(dependence_edge.src != i){
                    hung_graph[dependence_edge.dest][i].weight = -1;
                    hung_graph[i][dependence_edge.dest].weight = -1;
                    expanded_nodes.insert(dependence_edge.dest);
                }

                for(int dest : dependent_graph[dependence_edge.dest]){
                    if(expanded_nodes.find(dest) == expanded_nodes.end()){
                        st.push_back(edge(dependence_edge.dest,dest,-1));
                        expanded_nodes.insert(dest);
                    }
                }
            }
        } 
    }

    for(int i = 0; i < node_count; ++i) {
        const vector<edge> preceding_nodes = in_degree[i];
        unordered_set<int> expanded_nodes;
        for(int j = 0; j < (int)preceding_nodes.size(); ++j){
            vector<edge> st;
            st.push_back(preceding_nodes[j]);
            while(!st.empty()){
                edge dependence_edge = st.back();
                st.pop_back();
                if(dependence_edge.src != i){
                    hung_graph[i][dependence_edge.dest].weight = -1;
                    hung_graph[dependence_edge.dest][i].weight = -1;
                    expanded_nodes.insert(dependence_edge.dest);
                }
                for(const edge& e : in_degree[dependence_edge.dest]){
                    if(expanded_nodes.find(e.dest) == expanded_nodes.end()){
                        st.push_back(e);
                        expanded_nodes.insert(e.dest);
                    }
                }
            }
        } 
    }

    return;
}


bool compare (const edge src, const edge target) {
    int src_weight = src.weight;
    int dest_weight = target.weight;
    return (src_weight < dest_weight);
}

//Sort cost-graph weight by ascending order
void solver::sort_weight(vector<vector<edge>>& graph) {
    int size = graph.size();
    for (int i = 0; i < size; i++) {
        stable_sort(graph[i].begin(),graph[i].end(),compare);
    }

    return;
}

vector<int> solver::nearest_neightbor() {
    vector<int> solution;
    int current_node;
    bool visit_arr[node_count];
    bool selected = false;
    int depCnt_arr[node_count];
    vector<vector<edge>> sorted_costgraph = cost_graph; 
    sort_weight(sorted_costgraph);
    //sort input based on weight for NN heurestic
    memset(depCnt_arr,0,node_count*sizeof(int));

    for (int i = 0; i < node_count; i++) {
        visit_arr[i] = false;
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            depCnt_arr[dependent_graph[i][k]]++;
        }
    }

    for (int i = 0; i < node_count; i++) {
        if (depCnt_arr[i] == 0 && !selected) {
            current_node = i;
            visit_arr[current_node] = true;
            break;
        }
    }

    solution.push_back(current_node);
    
    int num = 1;
    int solution_cost = 0;
    
    while (num < node_count) {
        for (long unsigned int i = 0; i < dependent_graph[current_node].size(); i++) {
            depCnt_arr[dependent_graph[current_node][i]]--;
        }
        for (auto node: sorted_costgraph[current_node]) {
            if (!visit_arr[node.dest] && !depCnt_arr[node.dest]) {
                current_node = node.dest;
                solution_cost += node.weight;
                solution.push_back(current_node);
                num++;
                visit_arr[node.dest] = true;
                break;
            }
        }
    }

    best_cost = solution_cost;
    return solution;
}


int solver::get_maxedgeweight() {
    int max = 0;
    for (int i = 0; i < node_count; i++) {
        for (auto edge_weight : cost_graph[i]) {
            int weight = edge_weight.weight;
            if (weight > max) max = weight;
        }
    }
    return max;
}

vector<vector<int>> solver::get_cost_matrix(int max_edge_weight) {
    vector<vector<int>> matrix(node_count);
    for(int i = 0; i < node_count; ++i){
		matrix[i] = vector<int>(node_count, max_edge_weight*2);
	}

    for (vector<edge> edge_list : hung_graph) {
        for (auto edge : edge_list) {
            int i = edge.src;
            int k = edge.dest;
            int weight = edge.weight;
            if (weight != -1 && i != k) {
                matrix[i][k] = weight * 2;
            }
        }
    }

    return matrix;
}