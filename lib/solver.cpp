#include "solver.hpp"
#include "hash.hpp"
#include <iostream>
#include <fstream>
#include <stack>
#include <sstream>
#include <cmath>
#include <chrono>
#include <limits>
#include <algorithm>
#include <cstring>
#include <numeric>
#include <functional>
#include <vector>
#include <thread>
#include <queue>
#include <deque>
#include <mutex>
#include <bits/stdc++.h>
#include <unordered_map>
#include <condition_variable>

using namespace std;
#define HIS_MEM_LIMIT 0.7
#define TABLE_SIZE 12582917
#define TIME_FRAME 0.3

static Hash_Map history_table(TABLE_SIZE);
static int best_cost = 0;
static vector<int> best_solution;
static deque<sop_state> GPQ;

//Variable for locking;
static mutex GPQ_lock, Sol_lock, Split_lock;
static mutex asssign_mutex;
static mutex thread_load_mutex;
static mutex request_mutex;
static condition_variable Idel;
static vector<int> selected_orgin;
static mutex Info_Lock;
static mutex print_mutex;

//Protected by lock
//atomic<unsigned long> total_failed_nodes (1);
//atomic<unsigned long> total_failed_depth (15);
//atomic<unsigned long> number_of_history (0);
//atomic<unsigned long> number_of_lb (0);
atomic<int> selected_thread (-1);
atomic<unsigned> active_thread (0);
atomic<unsigned> idle_counter (0);
atomic<bool> time_out (false);
atomic<bool> limit_insert (false);
atomic<bool> freed(false);

//Data Collection
load_stats* thread_load;
//unsigned long* enumerated_bounds;
//int* steal_request;
//float* wait_time;

//Config variable (Reading Only)
static string enum_option;
static string assign_option;
static int global_pool_size = 0;
static int local_pool_size = 0;
static int t_limit = 0;
static int thread_total = 0;
static int local_depth = 0;

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
    local_depth = atoi(setting[4].c_str());
    return;
}

int solver::dynamic_hungarian(int src, int dest) {
    problem_state.hungarian_solver.fix_row(src, dest);
	problem_state.hungarian_solver.fix_column(dest, src);
	problem_state.hungarian_solver.solve_dynamic();
    //number_of_lb++;
    return problem_state.hungarian_solver.get_matching_cost()/2;
}

bool solver::HistoryUtilization(pair<vector<bool>,int>* key,int* lowerbound,bool* found,int cost) {
    size_t val = hash<vector<bool>>{}(key->first);
    int bucket_location = (val + key->second) % TABLE_SIZE;
    history_table.lock_table(bucket_location);
    HistoryNode* history_node = history_table.retrieve(key,bucket_location);

    if (history_node == NULL) {
        history_table.unlock_table(bucket_location);
        return true;
    }
    //number_of_history++;
    
    *found = true;
    int history_prefix = history_node->prefix_cost;
    int history_lb = history_node->lower_bound;
    *lowerbound = history_lb;

    if (cost >= history_prefix) {
        history_table.unlock_table(bucket_location);
        return false;
    }

    int imp = history_prefix - cost;
    
    if (imp <= history_lb - best_cost) {
        history_table.unlock_table(bucket_location);
        return false;
    }

    history_node->prefix_cost = cost;
    history_node->lower_bound = history_lb - imp;
    *lowerbound = history_node->lower_bound;

    history_table.unlock_table(bucket_location);
    return true;
}

bool Split_sort(const sop_state& src, const sop_state& dest) {
    if (src.cur_solution.size() == dest.cur_solution.size()) return src.load_info < dest.load_info;
    return src.cur_solution.size() < dest.cur_solution.size();
}

bool GPQ_sort(const sop_state& src, const sop_state& dest) {
    if (src.load_info == dest.load_info) return src.cur_solution.size()> dest.cur_solution.size();
    return src.load_info > dest.load_info;
}

bool local_sort(const sop_state& src, const sop_state& dest) {
    return src.load_info > dest.load_info;
}

bool bound_sort(const node& src,const node& dest) {
    if (src.lb == dest.lb) return src.nc > dest.nc;
    return src.lb > dest.lb;
}

bool nearest_sort(const node& src,const node& dest) {
    if (src.nc == dest.nc) return src.lb > dest.lb;
    return src.nc > dest.nc;
}

bool solver::Wlkload_Request(int i) {
    bool terminate = true;
    if (i == problem_state.initial_depth && active_thread >= 1) {
        enumerated_nodes = 0;
        active_thread--;
        thread_load_mutex.lock();
        thread_load[thread_id].out_of_work = true;
        thread_load_mutex.unlock();

        if (!local_pool->empty()) {
            active_thread++;
            problem_state = local_pool->back();
            local_pool->pop_back();
            terminate = false;
        }
        else {
            terminate = Grab_from_GPQ(false);
            if (active_thread > 1 && terminate) terminate = Steal_Workload();
            else if (active_thread == 1) notify_finished();
        }

        if (!terminate) {
            thread_load_mutex.lock();
            thread_load[thread_id].out_of_work = false;
            if (!local_pool->empty()) thread_load[thread_id].load = local_pool->back().load_info;
            else thread_load[thread_id].load = INT_MAX;
            thread_load_mutex.unlock();
        }
    }
    return terminate;
}

bool solver::Steal_Workload() {
    bool terminate = true;
    std::unique_lock<std::mutex> idel_lck(Split_lock);
                //auto start_time_wait = std::chrono::system_clock::now();
    GPQ_lock.lock();
    ////////////////
    //print_mutex.lock();
    //cout << std::chrono::duration<double>(std::chrono::system_clock::now() - start_time_limit).count() << ": " << "Thread " << thread_id << " started to wait idle count = " << idle_counter << " and active threads are " << active_thread << endl;
    //print_mutex.unlock();
    ///////////////
    int min = INT_MAX;
    int id = -1;
    
    //Select the thread with the lowest LB;
    if (selected_thread == -1) {
        thread_load_mutex.lock();
        for (int k = 0; k < thread_total; k++) {
            if (thread_load[k].load < min && thread_load[k].load > 0) {
                min = thread_load[k].load;
                id = k;
            }
        }
        selected_thread = id;
        thread_load_mutex.unlock();
    }
    //Else, we know a selection exists;
    while (GPQ.empty() && active_thread > 1) {
        if (idle_counter < (unsigned)thread_total) idle_counter++;
        GPQ_lock.unlock();
        Idel.wait(idel_lck);
        GPQ_lock.lock();
        idle_counter--;
        //steal_request[thread_id]++;
    }
    //Turn off thread selection if all of the threads are satisfied;
    if (idle_counter == 0) selected_thread = -1;
    //auto end_time_wait = std::chrono::system_clock::now();
    //auto wait_duration = (float) std::chrono::duration<double>(end_time_wait - start_time_wait).count();
    //wait_time[thread_id] += wait_duration;
    ////////////////
    //print_mutex.lock();
    //cout << std::chrono::duration<double>(std::chrono::system_clock::now() - start_time_limit).count() << ": " << "Thread " << thread_id << " released and finished waiting in " << wait_duration << " with GPQ size = " << GPQ.size() << endl;
    //print_mutex.unlock();
    ///////////////

    if (!GPQ.empty()) {
        active_thread++;
        problem_state = GPQ.back();
        GPQ.pop_back();
        terminate = false;
    }
    GPQ_lock.unlock();

    return terminate;
}

//Get new subproblem state and push it back to the pool;
void solver::assign_workload(int taken_n, int lb, const string dest) {
    //It is possible for other thread to update best solution and cause this subproblem to be pruned.
    sop_state target = problem_state;
    int taken_node = taken_n;
    int cur_node = target.cur_solution.back();
    target.taken_arr[taken_node] = 1;
    for (int vertex : dependent_graph[taken_node]) target.depCnt[vertex]--;
    target.cur_cost += cost_graph[cur_node][taken_node].weight;
    target.cur_solution.push_back(taken_node);
    target.hungarian_solver.fix_row(cur_node, taken_node);
    target.hungarian_solver.fix_column(taken_node, cur_node);
    target.initial_depth = target.cur_solution.size();
    target.load_info = lb;
    if (dest == "GPQ") {
        GPQ_lock.lock();
        GPQ.push_back(target);
        GPQ_lock.unlock();
    }
    else if (dest == "LOCAL") local_pool->push_back(target);
    return;
}

bool solver::push_to_global_pool(unsigned size) {
    bool insertion = false;

    while (size > 0) {
        if (local_pool->back().load_info >= best_cost) {
            local_pool->pop_back();
        }
        else {
            GPQ_lock.lock();
            GPQ.push_back(local_pool->back());
            GPQ_lock.unlock();
            local_pool->pop_back();
            insertion = true;
        }
        size--;
    }

    return insertion;
}

void solver::notify_finished() {
    if (active_thread == 1) {
        ////////////////
        //print_mutex.lock();
        //cout << std::chrono::duration<double>(std::chrono::system_clock::now() - start_time_limit).count() << ": " << "Thread " << thread_id << " called release all threads!" << endl;
        //print_mutex.unlock();
        ////////////////
        std::unique_lock<std::mutex> idel_lck(Split_lock);
        idel_lck.unlock();
        Idel.notify_all();
    }
    return;
}

void solver::push_to_historytable(pair<vector<bool>,int> key,int lower_bound) {
    if (problem_state.full_solution) {
        history_table.insert(key,problem_state.cur_cost,problem_state.suffix_cost,best_cost,thread_id);
    }
    else {
        history_table.insert(key,problem_state.cur_cost,lower_bound,best_cost,thread_id);
    }
    return;
}

bool solver::Split_level_check() {
    unsigned target_level = GPQ.front().cur_solution.size();
    for (auto node : GPQ) {
        if (node.cur_solution.size() != target_level) return true;
    }

    return false;
}

bool solver::Split_local_pool() {
    vector<int> ready_list;
    //Split until we have more than 1 node.
    while (local_pool->size() <= 1) {
        if (local_pool->empty()) return false;

        auto target = local_pool->front();
        local_pool->pop_front();

        //Create key for history table
        vector<bool> bit_vector(node_count, false);
        for (auto node : target.cur_solution) {
            bit_vector[node] = true;
        }
        int last_element = target.cur_solution.back();
        auto key = make_pair(bit_vector,last_element);

        //Get all children
        for (int i = node_count-1; i >= 0; i--) {
            if (!target.depCnt[i] && !target.taken_arr[i]) ready_list.push_back(i);
        }

        if (!ready_list.empty()) {
            for (auto node : ready_list) {
                int taken_node = node;
                int cur_node = target.cur_solution.back();
                target.taken_arr[taken_node] = 1;
                for (int vertex : dependent_graph[taken_node]) target.depCnt[vertex]--;
                target.cur_cost += cost_graph[cur_node][taken_node].weight;
                target.cur_solution.push_back(taken_node);
                target.hungarian_solver.fix_row(cur_node, taken_node);
                target.hungarian_solver.fix_column(taken_node, cur_node);

                if (target.cur_cost >= best_cost) {
                    //enumerated_bounds[thread_id]++;
                    key.first[taken_node] = false;
                    key.second = cur_node;
                    //total_failed_nodes++;
                    //total_failed_depth += target.cur_solution.size();
                    target.cur_solution.pop_back();
                    target.cur_cost -= cost_graph[cur_node][taken_node].weight;
                    target.taken_arr[taken_node] = 0;
                    for (int vertex : dependent_graph[taken_node]) target.depCnt[vertex]++;
                    target.hungarian_solver.undue_row(cur_node, taken_node);
                    target.hungarian_solver.undue_column(taken_node, cur_node);
                    continue;
                }
                if (target.cur_solution.size() == (size_t)node_count) {
                    //enumerated_bounds[thread_id]++;
                    if (target.cur_cost < best_cost) {
                        Sol_lock.lock();
                        best_solution = target.cur_solution;
                        best_cost = target.cur_cost;
                        Sol_lock.unlock();
                    }
                    key.first[taken_node] = false;
                    key.second = cur_node;
                    //total_failed_nodes++;
                    //total_failed_depth += target.cur_solution.size();
                    target.cur_solution.pop_back();
                    target.cur_cost -= cost_graph[cur_node][taken_node].weight;
                    target.taken_arr[taken_node] = 0;
                    for (int vertex : dependent_graph[taken_node]) target.depCnt[vertex]++;
                    target.hungarian_solver.undue_row(cur_node, taken_node);
                    target.hungarian_solver.undue_column(taken_node, cur_node);
                    continue;
                }
                else {
                    key.first[taken_node] = true;
                    key.second = taken_node;
                    int temp_lb = -1;
                    bool taken = false;
                    bool decision = HistoryUtilization(&key,&temp_lb,&taken,target.cur_cost);
                    if (!taken) {
                        //number_of_lb++;
                        target.hungarian_solver.solve_dynamic();
                        temp_lb = target.hungarian_solver.get_matching_cost()/2;
                        if (history_table.get_cur_size() < HIS_MEM_LIMIT * history_table.get_max_size()) {
                            if (target.full_solution) {
                                history_table.insert(key,target.cur_cost,target.suffix_cost,thread_id,target.cur_solution.size());
                            }
                            else history_table.insert(key,target.cur_cost,temp_lb,thread_id,target.cur_solution.size());
                        }
                    }
                    else if (taken && !decision) {
                        //enumerated_bounds[thread_id]++;
                        key.first[taken_node] = false;
                        key.second = cur_node;
                        //total_failed_nodes++;
                        //total_failed_depth += target.cur_solution.size();
                        target.cur_solution.pop_back();
                        target.cur_cost -= cost_graph[cur_node][taken_node].weight;
                        target.taken_arr[taken_node] = 0;
                        for (int vertex : dependent_graph[taken_node]) target.depCnt[vertex]++;
                        target.hungarian_solver.undue_row(cur_node, taken_node);
                        target.hungarian_solver.undue_column(taken_node, cur_node);
                        continue;
                    }
                    if (temp_lb >= best_cost) {
                        //enumerated_bounds[thread_id]++;
                        key.first[taken_node] = false;
                        key.second = cur_node;
                        //total_failed_nodes++;
                        //total_failed_depth += target.cur_solution.size();
                        target.cur_solution.pop_back();
                        target.cur_cost -= cost_graph[cur_node][taken_node].weight;
                        target.taken_arr[taken_node] = 0;
                        for (int vertex : dependent_graph[taken_node]) target.depCnt[vertex]++;
                        target.hungarian_solver.undue_row(cur_node, taken_node);
                        target.hungarian_solver.undue_column(taken_node, cur_node);
                        continue;
                    }
                    target.load_info = temp_lb;
                }
                //enumerated_bounds[thread_id]++;
                local_pool->push_back(target);
                key.first[taken_node] = false;
                key.second = cur_node;
                target.taken_arr[taken_node] = 0;
                for (int vertex : dependent_graph[taken_node]) target.depCnt[vertex]++;
                target.cur_cost -= cost_graph[cur_node][taken_node].weight;
                target.cur_solution.pop_back();
                target.hungarian_solver.undue_row(cur_node, taken_node);
                target.hungarian_solver.undue_column(taken_node, cur_node);
            }
        }
        ready_list.clear();
    }

    sort(local_pool->begin(),local_pool->end(),local_sort);
    return true;
}

void solver::Direct_Workload() {
    int min = INT_MAX;
    int id = -1;
    //Direct workload to new thread;
    thread_load_mutex.lock();
    for (int k = 0; k < thread_total; k++) {
        if (thread_load[k].load < min && thread_load[k].load > 0 && thread_load[k].out_of_work == false) {
            min = thread_load[k].load;
            id = k;
        }
    }
    selected_thread = id;
    thread_load_mutex.unlock();
    return;
}

void solver::Check_And_Distribute_Wlkload() {
    if (stolen) {
        auto time_frame_end = std::chrono::system_clock::now();
        if (std::chrono::duration<double>(time_frame_end - time_frame_start).count() > TIME_FRAME) {
            stolen = false;
            thread_load_mutex.lock();
            thread_load[thread_id].out_of_work = false;
            thread_load_mutex.unlock();
        }
    }

    bool insert = false;
    bool push_to_global = false;

    if (thread_id == selected_thread && idle_counter > 0) {
        if (!local_pool->empty()) {
            // Make sure we have at least one children in the ready list before splitting.
            unsigned release_size = 0;
            if (local_pool->size() <= 1) {
                if (Split_local_pool()) insert = true;
                else {
                    //print_mutex.lock();
                    //cout << "false insertion at thread " << thread_id << " and current total active threads are " << active_thread << endl;
                    //print_mutex.unlock();
                    insert = false;
                }
            }
            else insert = true;

            if (insert) {
                release_size = local_pool->size() / 2;
                if (push_to_global_pool(release_size)) push_to_global = true;
                else push_to_global = false;
            }

            thread_load_mutex.lock();
            if (!local_pool->empty()) thread_load[thread_id].load = local_pool->back().load_info;
            else thread_load[thread_id].load = INT_MAX;
            thread_load_mutex.unlock();

            if (push_to_global) {
                ///////////////
                //print_mutex.lock();
                //cout << std::chrono::duration<double>(std::chrono::system_clock::now() - start_time_limit).count() << ": " << "Thread " << thread_id << " just assigned and released potentially " << release_size << " threads\n";
                //print_mutex.unlock();
                //////////////

                if (release_size > idle_counter) {
                    std::unique_lock<std::mutex> idel_lck(Split_lock);
                    idel_lck.unlock();
                    Idel.notify_all();
                }
                else {
                    for (unsigned i = 0; i < release_size; i++) {
                        if (idle_counter > 0) {
                            std::unique_lock<std::mutex> idel_lck(Split_lock);
                            idel_lck.unlock();
                            Idel.notify_one();
                        }
                        else break;
                    }
                    Direct_Workload();
                }
                //Start time frame so this thread won't be appointed victim again for a certain period of time.
                stolen = true;
                time_frame_start = std::chrono::system_clock::now();
                thread_load_mutex.lock();
                thread_load[thread_id].out_of_work = true;
                thread_load_mutex.unlock();
            }
        }

        if (!push_to_global) Direct_Workload();
    }

    return;
}

int solver::enumerate(int i) {
    if (time_out) return -1;
    int next_level = i + 1;
    int max_depth_below = 0;

    if (abandon_work && next_level - 1 != problem_state.initial_depth) return -1;

    while (true) {
    deque<node> ready_list;

    for (int i = node_count-1; i >= 0; i--) {
        if (!problem_state.depCnt[i] && !problem_state.taken_arr[i]) {
            //Push vertices with 0 depCnt into the ready list
            ready_list.push_back(node(i,-1));
        }
    }
    //update vertex dependent count for current node neighbors
    int taken_node = 0;
    int u = 0;
    int v = 0;

    //Create key for history table;
    vector<bool> bit_vector(node_count, false);
    for (auto node : problem_state.cur_solution) {
        bit_vector[node] = true;
    }
    int last_element = problem_state.cur_solution.back();
    auto key = make_pair(bit_vector,last_element);

    if (!problem_state.cur_solution.empty()) {
        for (int i = 0; i < (int)ready_list.size(); i++) {
            node dest = ready_list[i];
            int src = problem_state.cur_solution.back();
            problem_state.cur_solution.push_back(dest.n);
            problem_state.cur_cost += cost_graph[src][dest.n].weight;
            int temp_lb = -1;
            bool taken = false;
            Check_And_Distribute_Wlkload();
            
            //Backtrack
            if (problem_state.cur_cost >= best_cost) {
                //enumerated_bounds[thread_id]++;
                //total_failed_nodes++;
                //total_failed_depth += problem_state.cur_solution.size();
                enumerated_nodes++;
                problem_state.cur_solution.pop_back();
                problem_state.cur_cost -= cost_graph[src][dest.n].weight;
                ready_list.erase(ready_list.begin()+i);
                i--;
                continue;
            }
            if (problem_state.cur_solution.size() == (size_t)node_count) {
                //enumerated_bounds[thread_id]++;
                enumerated_nodes++;
                if (problem_state.cur_cost < best_cost) {
                    Sol_lock.lock();
                    best_solution = problem_state.cur_solution;
                    best_cost = problem_state.cur_cost;
                    Sol_lock.unlock();
                }
                problem_state.full_solution = true;
                problem_state.suffix_cost = problem_state.cur_cost;
                //total_failed_nodes++;
                //total_failed_depth += problem_state.cur_solution.size();
                problem_state.cur_solution.pop_back();
                problem_state.cur_cost -= cost_graph[src][dest.n].weight;
                ready_list.erase(ready_list.begin()+i);
                i--;
                continue;
            }

            else {
                key.first[dest.n] = true;
                key.second = dest.n;
                bool decision = HistoryUtilization(&key,&temp_lb,&taken,problem_state.cur_cost);
                if (!taken) {
                    temp_lb = dynamic_hungarian(src,dest.n);
                    if (history_table.get_cur_size() < HIS_MEM_LIMIT * history_table.get_max_size()) push_to_historytable(key,temp_lb);
                    else if (!limit_insert) limit_insert = true;
                    problem_state.hungarian_solver.undue_row(src,dest.n);
                    problem_state.hungarian_solver.undue_column(dest.n,src);
                }
                else if (taken && !decision) {
                    enumerated_nodes++;
                    //enumerated_bounds[thread_id]++;
                    //total_failed_nodes++;
                    //total_failed_depth += problem_state.cur_solution.size();
                    problem_state.cur_solution.pop_back();
                    problem_state.cur_cost -= cost_graph[src][dest.n].weight;
                    key.first[dest.n] = false;
                    key.second = last_element;
                    ready_list.erase(ready_list.begin()+i);
                    i--;
                    continue;
                }
                if (temp_lb >= best_cost) {
                    enumerated_nodes++;
                    //if (!taken) push_to_historytable(key,temp_lb,i);
                    //enumerated_bounds[thread_id]++;
                    //total_failed_nodes++;
                    //total_failed_depth += problem_state.cur_solution.size();
                    problem_state.cur_solution.pop_back();
                    problem_state.cur_cost -= cost_graph[src][dest.n].weight;
                    key.first[dest.n] = false;
                    key.second = last_element;
                    ready_list.erase(ready_list.begin()+i);
                    i--;
                    continue;
                }
                enumerated_nodes++;
                //enumerated_bounds[thread_id]++;
                problem_state.cur_solution.pop_back();
                problem_state.cur_cost -= cost_graph[src][dest.n].weight;
                ready_list[i].nc = cost_graph[src][dest.n].weight;
                ready_list[i].pushed = taken;
                ready_list[i].lb = temp_lb;
                key.first[dest.n] = false;
                key.second = last_element;
            }
        }
        if (enum_option == "DH") sort(ready_list.begin(),ready_list.end(),bound_sort);
        else if (enum_option == "NN") sort(ready_list.begin(),ready_list.end(),nearest_sort);
    }

    int lb = -1;
    bool pushed_to_local = false;
    bool pushed_to_his = false;

    while(!ready_list.empty()) {
        //Take the choosen node and back track if ready list is empty
        pushed_to_his = ready_list.back().pushed;
        lb = ready_list.back().lb;
        taken_node = ready_list.back().n;
        key.first[taken_node] = true;
        key.second = taken_node;
        ready_list.pop_back();

        if (local_pool->size() < ((size_t)local_pool_size / 4) && !ready_list.empty()) {
            //int division = int(total_failed_depth/total_failed_nodes) - 10;
            if (i <= local_depth) {
                while (!ready_list.empty() && local_pool->size() < (size_t)local_pool_size) {
                    if (!pushed_to_local) pushed_to_local = true;
                    assign_workload(ready_list.back().n,ready_list.back().lb,"LOCAL");
                    ready_list.pop_back();
                }
            }
            else if (idle_counter > 0 && local_pool->size() <= 1) {
                while (!ready_list.empty() && local_pool->size() < (size_t)local_pool_size) {
                    if (!pushed_to_local) pushed_to_local = true;
                    assign_workload(ready_list.back().n,ready_list.back().lb,"LOCAL");
                    ready_list.pop_back();
                }
            }
            if (!local_pool->empty()) {
                sort(local_pool->begin(),local_pool->end(),local_sort);
                thread_load_mutex.lock();
                thread_load[thread_id].load = local_pool->back().load_info;
                thread_load_mutex.unlock();
            }
            else {
                thread_load_mutex.lock();
                thread_load[thread_id].load = INT_MAX;
                thread_load_mutex.unlock();
            }
        }

        if (!problem_state.cur_solution.empty()) {
            u = problem_state.cur_solution.back();
            v = taken_node;
            problem_state.hungarian_solver.fix_row(u, v);
            problem_state.hungarian_solver.fix_column(v, u);
            problem_state.cur_cost += cost_graph[u][v].weight;
        }

        for (int vertex : dependent_graph[taken_node]) problem_state.depCnt[vertex]--;
        problem_state.cur_solution.push_back(taken_node);
        problem_state.taken_arr[taken_node] = 1;
        problem_state.full_solution = false;
        problem_state.suffix_cost = 0;
        Check_And_Distribute_Wlkload();
        
        int children_depth = enumerate(next_level);
        if (children_depth > max_depth_below) max_depth_below = children_depth;

        if (limit_insert && !pushed_to_his) {
            //if (max_depth_below > 1) push_to_historytable(key,lb);
            if (max_depth_below > 1 && (history_table.get_cur_size() < (HIS_MEM_LIMIT+0.1) * history_table.get_max_size())) {
                push_to_historytable(key,lb);
            }
            else if (max_depth_below > 2 && (history_table.get_cur_size() < (HIS_MEM_LIMIT+0.2) * history_table.get_max_size())) {
                push_to_historytable(key,lb);
            }
            else if (max_depth_below > 4 && (history_table.get_cur_size() < history_table.get_max_size())) {
                push_to_historytable(key,lb);
            }
        }
        
        //////////////////////////
        //for (auto node : problem_state.cur_solution) {
        //    cout << node << ",";
        //}
        //cout << "---->" << max_depth_below;
        //cout << endl;
        //////////////////////////

        for (int vertex : dependent_graph[taken_node]) problem_state.depCnt[vertex]++;
        problem_state.taken_arr[taken_node] = 0;
        problem_state.cur_solution.pop_back();
        problem_state.cur_cost -= cost_graph[u][v].weight;
        if (problem_state.cur_solution.size() >= 2) {
            problem_state.hungarian_solver.undue_row(u,v);
            problem_state.hungarian_solver.undue_column(v,u);
        }
        key.first[taken_node] = false;
        key.second = problem_state.cur_solution.back();

        auto cur_time = std::chrono::system_clock::now();
        if (std::chrono::duration<double>(cur_time - start_time_limit).count() > t_limit) {
            time_out = true;
            active_thread = 0;
            std::unique_lock<std::mutex> idel_lck(Split_lock);
            idel_lck.unlock();
            Idel.notify_all();
            return -1;
        }

        if (pushed_to_local && !(local_pool->empty())) {
            pushed_to_local = false;
            unsigned starting_size = problem_state.cur_solution.size() + 1;
            vector<sop_state> selected_state;

            auto back_up_state = problem_state;

            for (int k = local_pool->size() - 1; k >= 0; k--) {
                if ((*local_pool)[k].cur_solution.size() >= starting_size) {
                    (*local_pool)[k].initial_depth = -1;
                    problem_state = (*local_pool)[k];
                    local_pool->erase(local_pool->begin() + k);
                    enumerate(problem_state.cur_solution.size()+1);
                    k = local_pool->size();
                }
            }

            problem_state = back_up_state;
        }
    }

    if (!Wlkload_Request(next_level-1)) {
        next_level = problem_state.initial_depth + 1;
    }
    else break;
    }

    return max_depth_below + 1;
}

void calculate_standard_deviation() {
    double mean = 0;
    for (auto node : GPQ) mean += node.load_info;
    mean = mean / GPQ.size();
    double summation = 0;

    for (auto node : GPQ) {
        summation += pow((node.load_info - mean),2);
    }
    summation = summation / GPQ.size();
    cout << "Average LB is " << mean << endl;
    summation = sqrt(summation);
    cout << "Stnadard deviation of the GPQ is " << summation << endl;
    cout << "Deviation percentage from the mean is " << summation / mean * 100 << endl;

    return;
}

bool solver::Grab_from_GPQ(bool reserve) {
    int assigned = false;
    int terminate = true;
    GPQ_lock.lock();
    if (!GPQ.empty()) {
        sort(GPQ.begin(),GPQ.end(),GPQ_sort);
        while (!assigned) {
            for (int k = GPQ.size() - 1; k >= 0; k--) {
                int origin_node = GPQ[k].originate;
                if (!count(selected_orgin.begin(),selected_orgin.end(),origin_node)) {
                    if (!reserve) {
                        active_thread++;
                        problem_state = GPQ[k];
                    }
                    else reserve_state = GPQ[k];
                    selected_orgin.push_back(origin_node);
                    GPQ.erase(GPQ.begin()+k);
                    terminate = false;
                    assigned = true;
                    break;
                }
            }
            if (!assigned) selected_orgin.clear();
        }
    }
    GPQ_lock.unlock();
    return terminate;
}

void solver::solve_parallel(int thread_num, int pool_size) {
    start_time_limit = std::chrono::system_clock::now();
    vector<solver> solvers(thread_num);
    vector<thread> Thread_manager(thread_num);
    vector<int> ready_thread;
    //Initially fill the GPQ with solvers:
    vector<node> ready_list;

    problem_state.cur_solution.push_back(0);
    problem_state.taken_arr[0] = 1;
    for (int vertex : dependent_graph[0]) problem_state.depCnt[vertex]--;

    for (int i = node_count-1; i >= 0; i--) {
        if (!problem_state.depCnt[i] && !problem_state.taken_arr[i]) {
            //Push vertices with 0 depCnt into the ready list
            ready_list.push_back(node(i,-1));
        }
    }

    //Initial filling of the GPQ
    for (auto node : ready_list) {
        sop_state target = problem_state;
        int taken_node = node.n;
        int cur_node = target.cur_solution.back();
        target.taken_arr[taken_node] = 1;
        for (int vertex : dependent_graph[taken_node]) target.depCnt[vertex]--;
        target.cur_cost += cost_graph[cur_node][taken_node].weight;
        target.cur_solution.push_back(taken_node);
        target.originate = target.cur_solution.back();
        target.hungarian_solver.fix_row(cur_node, taken_node);
        target.hungarian_solver.fix_column(taken_node, cur_node);
        target.hungarian_solver.solve_dynamic();
        target.load_info = target.hungarian_solver.get_matching_cost()/2;
        //number_of_lb++;
        GPQ.push_back(target);
    }

    for (int i = thread_num - 1; i >= 0; i--) ready_thread.push_back(i);

    //While GPQ is not empty do split operation or assign threads with new node.
    while (!GPQ.empty()) {
        while (GPQ.size() < (size_t)pool_size || Split_level_check()) {
            sort(GPQ.begin(),GPQ.end(),Split_sort);
            auto target = GPQ.front();
            if (target.cur_solution.size() == (unsigned) (node_count - 1)) break;
            GPQ.pop_front();
            ready_list.clear();
            for (int i = node_count-1; i >= 0; i--) {
                if (!target.depCnt[i] && !target.taken_arr[i]) ready_list.push_back(node(i,-1));
            }

            if (!ready_list.empty()) {
                for (auto node : ready_list) {
                    int vertex = node.n;
                    //Push split node back into GPQ
                    int taken_node = vertex;
                    int cur_node = target.cur_solution.back();
                    target.taken_arr[taken_node] = 1;
                    for (int vertex : dependent_graph[taken_node]) target.depCnt[vertex]--;
                    target.cur_cost += cost_graph[cur_node][taken_node].weight;
                    target.cur_solution.push_back(taken_node);
                    if (target.cur_solution.size() == (size_t)node_count && problem_state.cur_cost < best_cost) {
                        best_solution = target.cur_solution;
                        best_cost = target.cur_cost;
                    }
                    target.hungarian_solver.fix_row(cur_node, taken_node);
                    target.hungarian_solver.fix_column(taken_node, cur_node);
                    target.hungarian_solver.solve_dynamic();
                    target.load_info = target.hungarian_solver.get_matching_cost()/2;
                    //number_of_lb++;
                    GPQ.push_back(target);

                    target.taken_arr[taken_node] = 0;
                    for (int vertex : dependent_graph[taken_node]) target.depCnt[vertex]++;
                    target.cur_cost -= cost_graph[cur_node][taken_node].weight;
                    target.cur_solution.pop_back();
                    target.hungarian_solver.undue_row(cur_node, taken_node);
                    target.hungarian_solver.undue_column(taken_node, cur_node);
                }
            }
            else GPQ.push_back(target);
        }
        sort(GPQ.begin(),GPQ.end(),GPQ_sort);

        int thread_cnt = 0;
        //cout << "GPQ initial depth is " << GPQ.back().cur_solution.size() << endl;
        cout << "Initial GPQ size is " << GPQ.size() << endl;
        history_table.adjust_max_size(node_count,GPQ.size());

        //calculate_standard_deviation();
        
        //Assign GPQ subproblem into solver;
        while (thread_cnt < thread_num) {
            vector<int> selected_originator;
            for (int k = GPQ.size() - 1; k >= 0; k--) {
                if (thread_cnt >= thread_num) break;
                unsigned origin = GPQ[k].originate;
                if (!count(selected_originator.begin(),selected_originator.end(),origin)) {
                    solvers[thread_cnt].problem_state = GPQ[k];
                    GPQ.erase(GPQ.begin()+k);
                    solvers[thread_cnt].node_count = node_count;
                    solvers[thread_cnt].problem_state.initial_depth = solvers[thread_cnt].problem_state.cur_solution.size();
                    solvers[thread_cnt].thread_id = thread_cnt;
                    solvers[thread_cnt].local_pool = new deque<sop_state>();
                    selected_originator.push_back(origin);
                    thread_cnt++;
                }
            }
        }

        for (int i = 0; i < thread_num; i++) {
            int size = solvers[i].problem_state.initial_depth;
            Thread_manager[i] = thread(&solver::enumerate,move(solvers[i]),size);
            active_thread++;
        }
        
        for (int i = 0; i < thread_num; i++) {
            if (Thread_manager[i].joinable()) {
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

void solver::solve(string filename,int thread_num,string assignment_scheme) {
    if (thread_num == -1 || thread_num == 0) {
        cerr << "Incorrect thread number input" << endl;
        exit(-1);
    }
    assign_option = assignment_scheme;
    thread_total = thread_num;
    retrieve_input(filename);
    //Remove redundant edges in the cost graph
    transitive_redundantcy();
    history_table.set_node_t(node_count);
    //Only uses nearest neighbor as initial heuristic.
    ///////////////////////////////////////////////////////////
    //vector<int> temp_solution;
    //temp_solution.push_back(0);
    //best_solution = nearest_neightbor(&temp_solution,&best_cost);
    ////////////////////////////////////////////////////////////
    
    if (node_count <= 100) {
        cout << "roll out heuristic initialized" << endl;
        best_solution = roll_out();
    }
    else {
        cout << "nearest neighbor heuristic initialized" << endl;
        vector<int> temp_solution;
        temp_solution.push_back(0);
        best_solution = nearest_neightbor(&temp_solution,&best_cost);
    }

    //cout << "best initial cost is " << best_cost << endl;

    int max_edge_weight = get_maxedgeweight();
    problem_state.hungarian_solver = Hungarian(node_count, max_edge_weight+1, get_cost_matrix(max_edge_weight+1));
    problem_state.hungarian_solver.start()/2;
    //picked_list = vector<bool>(node_count,false);
    //EGB_static_lowerbound = mmcp_lb();
    problem_state.depCnt = vector<int>(node_count,0);
    problem_state.taken_arr = vector<int>(node_count,0);

    float average_depCnt = 0;

    for (int i = 0; i < node_count; i++) {
        for (unsigned k = 0; k < dependent_graph[i].size(); k++) {
            if (i != 0 && i != node_count - 1) average_depCnt++;
            problem_state.depCnt[dependent_graph[i][k]]++;
        }
    }

    average_depCnt = average_depCnt / (node_count - 2);
    cout << "Average dependance count is " << average_depCnt << endl;
    cout << "Instance size is " << node_count - 2 << endl;
    //cout << "best solution found using initial heuristic is " << best_cost << endl;
    /*
    cout << "the NN solution contains ";
    for (int i = 0; i < node_count; i++) {
        if (i != node_count - 1) cout << best_solution[i] << "-->";
        else cout << best_solution[i];
    }
    cout << endl;
    cout << "MMCP-based LB is " << MMCP_static_lowerbound << endl;
    */
    thread_load = new load_stats [thread_total];
    //enumerated_bounds = new unsigned long [thread_total];
    //steal_request = new int [thread_total];
    //wait_time = new float [thread_total];
    /*
    explored_nodes = new int [thread_total];
    explored_LB = new int [thread_total];
    LB_time = new float [thread_total]; 
    logic_time = new float [thread_total];
    node_time = new float [thread_total];
    history_wait = new float [thread_total];
    */
    //memset(enumerated_bounds,0,thread_total * sizeof(unsigned long));
    //memset(steal_request,0,thread_total * sizeof(int));
    //memset(wait_time,0,thread_total * sizeof(float));

    history_table.set_up_mem(thread_num,node_count);
    history_table.set_up_table(node_count);
    /*
    memset(explored_nodes,0,thread_total * sizeof(int));
    memset(explored_LB,0,thread_total * sizeof(int));
    memset(logic_time,0,thread_total * sizeof(float));
    memset(LB_time,0,thread_total * sizeof(float));
    memset(wait_time,0,thread_total * sizeof(float));
    memset(node_time,0,thread_total * sizeof(float));
    memset(history_wait,0,thread_total * sizeof(float));
    */
    auto start_time = chrono::high_resolution_clock::now();
    solve_parallel(thread_total,global_pool_size);
    auto end_time = chrono::high_resolution_clock::now();

    auto total_time = chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    cout << "------------------------" << thread_total << " thread" << "------------------------------" << endl;
    cout << enum_option << ": " << best_cost << "," << setprecision(4) << total_time / (float)(1000000) << endl;

    //long int sum = 0;
    //int thread_i = -1;
    //float max_wait = 0;
    
    /*
    for (int i = 1; i <= thread_total; i++) {
        if (i != thread_total) cout << enumerated_bounds[i-1] << ",";
        else cout << enumerated_bounds[i-1] << endl;
        sum += enumerated_bounds[i-1];
    }
    cout << "Total enumerated nodes = " << sum << endl;
    

    for (int i = 1; i <= thread_total; i++) {
        if (wait_time[i-1] > max_wait) {
            max_wait = wait_time[i-1];
            thread_i = i;
        }
        if (i != thread_total) cout << steal_request[i-1] << ":" << wait_time[i-1] << ",";
        else cout << steal_request[i-1] << ":" << wait_time[i-1] << endl;
    }
    */

    //cout << "......................Work Donation......................" << endl;
    //cout << "maximum wait time = " << max_wait << " happening in thread " << thread_i << " with total steal request = " << steal_request[thread_i-1] << endl;
    //cout << "......................Process Statistics.................." << endl;
    //cout << "Table look up: " << number_of_history << ", LB calculation: " << number_of_lb << endl;
    cout << "History table size is " << history_table.get_cur_size() / 1000000 << " megabyte" << endl;
    history_table.average_size();
    /*
    cout << "......................Global pool node info.............." << endl;
    cout << "Total nodes in the global pool are " << Info_arr.size() << endl;
    for (auto info : Info_arr) {
        info.Display_info();
    }
    */

    /*
    cout << "......................History Table......................" << endl;
    history_table.average_size();
    cout << "Total number of waits for the history table is " << history_table.num_of_waits << endl;
    cout << "...........................Work load info..........................." << endl;
    unsigned maximum = 0;
    unsigned minimum = INT_MAX;
    unsigned total_size = 0;
    unsigned size = 0;

    for (auto info : Info_arr) {
        if (info.Stolen_Load > maximum) maximum = info.Stolen_Load;
        if (info.Stolen_Load < minimum) minimum = info.Stolen_Load;
        total_size += info.Stolen_Load;
        size++;
    }

    cout << "number of steals are " << size << endl;
    cout << "maximum number of workload in a load is " << maximum << endl;
    cout << "minimum number of workload in a load is " << minimum << endl;
    cout << "average load is " << (float)total_size / (float)size << endl;
    
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
    if (inFile.fail()) {
        cerr << "Error: input file " << filename << " -> " << strerror(errno) << endl;
        exit(-1);
    }

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


vector<int> solver::roll_out() {
    vector<int> solution;
    int solution_cost = 0;
    int global_min = INT_MAX;
    int current_node;
    bool visit_arr[node_count];
    int depCnt_arr[node_count];
    //sort input based on weight for NN heurestic
    memset(depCnt_arr,0,node_count*sizeof(int));

    for (int i = 0; i < node_count; i++) {
        visit_arr[i] = false;
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            depCnt_arr[dependent_graph[i][k]]++;
        }
    }

    solution.push_back(0);
    current_node = 0;
    visit_arr[0] = true;
    
    int num = 1;

    start_time_limit = std::chrono::system_clock::now();

    while (num < node_count) {
        for (long unsigned int i = 0; i < dependent_graph[current_node].size(); i++) {
            depCnt_arr[dependent_graph[current_node][i]]--;
        }
        int local_min = INT_MAX;
        int taken = -1;
        int incre = -1;
        for (auto node: cost_graph[current_node]) {
            if (!visit_arr[node.dest] && !depCnt_arr[node.dest]) {
                solution.push_back(node.dest);
                solution_cost += node.weight;
                int initial_cost = INT_MAX;
                vector<int> initial_solution = nearest_neightbor(&solution,&initial_cost);
                int improved_cost = tour_improvement(initial_solution,initial_cost,num);
                if (improved_cost <= local_min) {
                    local_min = improved_cost;
                    taken = node.dest;
                    incre = node.weight;
                    if (improved_cost < global_min) global_min = improved_cost;
                }
                solution.pop_back();
                solution_cost -= node.weight;
                auto cur_time = std::chrono::system_clock::now();
                if (std::chrono::duration<double>(cur_time - start_time_limit).count() > 1) {
                    best_cost = global_min;
                    return solution;
                }
            }
        }
        solution.push_back(taken);
        solution_cost += incre;
        current_node = taken;
        visit_arr[taken] = true;
        num++;
    }
    best_cost = global_min;
    return solution;
}


vector<int> solver::nearest_neightbor(vector<int>* partial_solution,int* initial_cost) {
    vector<int> solution;
    int current_node;
    int solution_cost = 0;
    bool visit_arr[node_count];
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

    current_node = 0;
    for (auto node : *partial_solution) {
        visit_arr[node] = true;
        for (long unsigned int i = 0; i < dependent_graph[node].size(); i++) {
            depCnt_arr[dependent_graph[node][i]]--;
        }
        solution.push_back(node);
        solution_cost += cost_graph[current_node][node].weight;
        current_node = node;
    }

    int num = solution.size();
    
    while (num < node_count) {
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
        for (long unsigned int i = 0; i < dependent_graph[current_node].size(); i++) {
            depCnt_arr[dependent_graph[current_node][i]]--;
        }
    }

    *initial_cost = solution_cost;
    return solution;
}

int solver::tour_improvement(vector<int> initial_solution,int initial_cost,int initial_depth) {
    vector<int> initial_tour = initial_solution;
    vector<int> sub_optimal_solution = initial_solution;
    int cost = initial_cost;
    int n = node_count;
    bool improvement = true;
    //forward exchange
    auto run_time_start = std::chrono::system_clock::now();

    for (int h = initial_depth; h < n - 3; h ++) {
        for (int i = h + 1; i < n - 2; i ++) {
            for (int j = i + 1; j < n - 1; j ++) {
                int r1 = initial_tour[h + 1];
                int r2 = initial_tour[i + 1];
                int r3 = initial_tour[j + 1];
                initial_tour[h + 1] = r2;
                initial_tour[i + 1] = r3;
                initial_tour[j + 1] = r1;
                //Replacement
                initial_tour[i + 1] = r1;
                initial_tour[j + 1] = r2;
                initial_tour[h + 1] = r3;
                int local_cost = 0;

                auto cur_time = std::chrono::system_clock::now();
                if (std::chrono::duration<double>(cur_time - run_time_start).count() > 1) return cost;

                if (check_satisfiablity(&local_cost,&initial_tour) == false) {
                    improvement = false;
                    initial_tour[i + 1] = r2;
                    initial_tour[h + 1] = r1;
                    initial_tour[j + 1] = r3;
                }
                else if (local_cost < cost) {
                    improvement = true;
                    cost = local_cost;
                    sub_optimal_solution = initial_tour;
                }
                else {
                    improvement = false;
                    initial_tour[i + 1] = r2;
                    initial_tour[h + 1] = r1;
                    initial_tour[j + 1] = r3;
                }
            }
            if (improvement) break;
        }
    }

    //backward exchange
    for (int h = n - 1; h >= initial_depth + 2; h --) {
        for (int i = h - 1; i >= initial_depth + 1; i --) {
            for (int j = i - 1; j >= initial_depth; j --) {
                int r1 = initial_tour[h + 1];
                int r2 = initial_tour[i + 1];
                int r3 = initial_tour[j + 1];
                initial_tour[h + 1] = r2;
                initial_tour[i + 1] = r3;
                initial_tour[j + 1] = r1;
                //Replacement
                initial_tour[i + 1] = r1;
                initial_tour[j + 1] = r2;
                initial_tour[h + 1] = r3;
                int local_cost = 0;

                auto cur_time = std::chrono::system_clock::now();
                if (std::chrono::duration<double>(cur_time - run_time_start).count() > 1) return cost;

                if (check_satisfiablity(&local_cost,&initial_tour) == false) {
                    improvement = false;
                    initial_tour[i + 1] = r2;
                    initial_tour[h + 1] = r1;
                    initial_tour[j + 1] = r3;
                }
                else if (local_cost < cost) {
                    improvement = true;
                    cost = local_cost;
                    sub_optimal_solution = initial_tour;
                    break;
                }
                else {
                    improvement = false;
                    initial_tour[i + 1] = r2;
                    initial_tour[h + 1] = r1;
                    initial_tour[j + 1] = r3;
                }
            }
            if (improvement) break;
        }
    }

    return cost;
}

bool solver::check_satisfiablity(int* local_cost, vector<int>* tour) {
    int current_node;
    bool visit_arr[node_count];
    int depCnt_arr[node_count];
    //sort input based on weight for NN heurestic
    memset(depCnt_arr,0,node_count*sizeof(int));

    for (int i = 0; i < node_count; i++) {
        visit_arr[i] = false;
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            depCnt_arr[dependent_graph[i][k]]++;
        }
    }

    current_node = 0;
    visit_arr[0] = true;
    
    int num = 1;
    
    while (num < node_count) {
        for (long unsigned int i = 0; i < dependent_graph[current_node].size(); i++) {
            depCnt_arr[dependent_graph[current_node][i]]--;
        }
        bool valid = false;
        for (auto node: cost_graph[current_node]) {
            if (!visit_arr[node.dest] && !depCnt_arr[node.dest] && node.dest == (*tour)[num]) {
                current_node = node.dest;
                *local_cost += node.weight;
                valid = true;
                visit_arr[node.dest] = true;
                num++;
                break;
            }
        }
        if (valid == false) return false;
    }
    return true;
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