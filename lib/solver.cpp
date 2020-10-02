#include "solver.hpp"
#include "hash.hpp"
#include <iostream>
#include <fstream>
#include <stack>
#include <sstream>
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
static string enum_option;
//static int recently_added = 0;

static int best_cost = 0;
static vector<int> best_solution;
static queue<solver> GPQ;
//vector<bool> picked_list;

//Variable for locking;
static mutex GPQ_lock, Sol_lock, Split_lock, Boardcast_lock;
static mutex Split_Call, asssign_mutex, unlock_mutex;
static condition_variable Idel,thread_counter;


//Protected by lock
atomic<bool> Split (false);
atomic<int> active_thread (1);
atomic<int> steal_counter (1);

//Config variable (Reading Only)
static int p_size = 0;
static int t_limit = 0;
static int thread_total = 0;
static int split_depth = 0;
static string Split_Opt;
static string Steal_Opt;

void solver::process_solution() {
    if (cur_cost < best_cost) {
        best_solution = cur_solution;
        best_cost = cur_cost;
    }
    return;
}

int solver::dynamic_hungarian(int src, int dest) {
    hungarian_solver.fix_row(src, dest);
	hungarian_solver.fix_column(dest, src);
	hungarian_solver.solve_dynamic();
    return hungarian_solver.get_matching_cost()/2;
}

bool solver::HistoryUtilization(int* lowerbound,bool* found,bool suffix_exist) {
    string bit_string(node_count, '0');

    for (auto node : cur_solution) {
        bit_string[node] = '1';
    }   

    int last_element = cur_solution.back();
    auto key = make_pair(bit_string,last_element);

    HistoryNode history_node = history_table.retrieve(key);
    if (history_node.prefix_sched.empty()) return true;

    *found = true;
    int history_prefix = history_node.prefix_cost;
    int history_lb = history_node.lower_bound;
    *lowerbound = history_lb;

    if (cur_cost >= history_prefix) return false;

    int imp = history_prefix - cur_cost;
    
    if (imp <= history_lb - best_cost) return false;
    history_node.prefix_sched = cur_solution;
    history_node.prefix_cost = cur_cost;
    history_node.lower_bound = history_lb - imp;
    *lowerbound = history_node.lower_bound;
    
    if (!history_node.suffix_sched.empty()) {
        vector<int> temp = cur_solution;
        vector<int> suffix_sched;
        suffix_sched = history_node.suffix_sched;
        int suf_cost = history_node.suffix_cost;
        temp.insert(temp.end(),suffix_sched.begin(),suffix_sched.end());
        Sol_lock.lock();
        best_solution = temp;
        Sol_lock.unlock();
        best_cost = cur_cost + suf_cost;
        history_node.lower_bound = best_cost;
        *lowerbound = history_node.lower_bound;
        history_table.insert(key,history_node);
        return false;
    }
    else {
        if (suffix_exist) {
             for (int i = suffix.size() - 1; i >= 0; i--) {
                history_node.suffix_sched.push_back(suffix[i]);
            }
            history_node.suffix_cost = suffix_cost;
            history_node.lower_bound = cur_cost + suffix_cost;
            history_table.insert(key,history_node);
            return false;
        }
    }

    history_table.insert(key,history_node);
    return true;
}

void solver::assign_historytable(int prefix_cost,int lower_bound,int i) {
    bool taken = false;
    int lb = 0;

    if (suffix.empty()) HistoryUtilization(&lb,&taken,false);
    else HistoryUtilization(&lb,&taken,true);

    if (!taken) {
        string bit_string(node_count, '0');
        for (auto node : cur_solution) bit_string[node] = '1';
        int last_element = cur_solution.back();
        auto key = make_pair(bit_string,last_element);
    
        if (suffix.empty()) {
            history_table.insert(key,HistoryNode(prefix_cost,lower_bound,cur_solution));
        }

        else {
            std::vector<int> suffix_sched;
            int size = suffix_sched.size();
            for (int i = size - 1; i >= 0; i--) suffix_sched.push_back(suffix[i]);
            history_table.insert(key,HistoryNode(prefix_cost,suffix_cost,prefix_cost+suffix_cost,cur_solution,suffix_sched));
        }
    }
    
}   

bool bound_sort(const node& src,const node& dest) {
    if (src.lb == dest.lb) return src.nc > dest.nc;
    return src.lb > dest.lb;
}

bool nearest_sort(const node& src,const node& dest) {
    if (src.nc == dest.nc) return src.lb > dest.lb;
    return src.nc > dest.nc;
}


void solver::enumerate(int i) {
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

    int bound = INT_MAX;
    if (!cur_solution.empty()) {
            for (int i = 0; i < (int)ready_list.size(); i++) {
                node dest = ready_list[i];
                int src = cur_solution.back();
                cur_solution.push_back(dest.n);
                cur_cost += cost_graph[src][dest.n].weight;
                int temp_lb = -1;
                bool taken = false;
                
                if (cur_cost >= best_cost) {
                    cur_solution.pop_back();
                    cur_cost -= cost_graph[src][dest.n].weight;
                    ready_list.erase(ready_list.begin()+i);
                    i--;
                    continue;
                }
                if (cur_solution.size() == node_count) {
                    if (cur_cost < best_cost) {
                        //Lock the best solution found so far.
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
                    bool decision = HistoryUtilization(&temp_lb,&taken,false);
                    if (!taken) {
                        string bit_string(node_count, '0');
                        for (auto node : cur_solution) bit_string[node] = '1';
                        int last_element = cur_solution.back();
                        auto key = make_pair(bit_string,last_element);
                        temp_lb = dynamic_hungarian(src,dest.n);
                        history_table.insert(key,HistoryNode(cur_cost,temp_lb,cur_solution));
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
                }
            }
        if (enum_option == "DH") sort(ready_list.begin(),ready_list.end(),bound_sort);
        else if (enum_option == "NN") sort(ready_list.begin(),ready_list.end(),nearest_sort);
    }

    //cout << "running with thread id = "<< std::this_thread::get_id() << endl;
    while (Split && i <= split_depth) {
        asssign_mutex.lock();
        if (!Split) {
            asssign_mutex.unlock();
            break;
        }
        // Make sure we have at least two children in the ready list before splitting.
        if (ready_list.size() >= 2) {
            //Limit spliting to be 5 levels down the tree.
            if (Split_Opt == "SINGLE") {
                int taken_node = ready_list.back().n;
                int cur_node = cur_solution.back();
                taken_arr[taken_node] = 1;
                for (int vertex : dependent_graph[taken_node]) depCnt[vertex]--;
                cur_cost += cost_graph[cur_node][taken_node].weight;
                cur_solution.push_back(taken_node);
                hungarian_solver.fix_row(cur_node, taken_node);
                hungarian_solver.fix_column(taken_node, cur_node);
                hungarian_solver.solve_dynamic();
                ready_list.pop_back();

                GPQ_lock.lock();
                GPQ.push(*this);
                GPQ_lock.unlock();

                taken_arr[taken_node] = 0;
                for (int vertex : dependent_graph[taken_node]) depCnt[vertex]++;
                cur_cost -= cost_graph[cur_node][taken_node].weight;
                cur_solution.pop_back();
                hungarian_solver.undue_row(cur_node, taken_node);
                hungarian_solver.undue_column(taken_node, cur_node);
                Split = false;
                //cout << "finish taking from GPQ with thread id = " << std::this_thread::get_id() << endl;
            }
            else if (Split_Opt == "FULL") {
                int p_count = 0;
                while (ready_list.size() > 1 && p_count < p_size) {
                    int taken_node = ready_list.back().n;
                    int cur_node = cur_solution.back();
                    taken_arr[taken_node] = 1;
                    for (int vertex : dependent_graph[taken_node]) depCnt[vertex]--;
                    cur_cost += cost_graph[cur_node][taken_node].weight;
                    cur_solution.push_back(taken_node);
                    hungarian_solver.fix_row(cur_node, taken_node);
                    hungarian_solver.fix_column(taken_node, cur_node);
                    hungarian_solver.solve_dynamic();
                    ready_list.pop_back();

                    GPQ_lock.lock();
                    GPQ.push(*this);
                    GPQ_lock.unlock();

                    taken_arr[taken_node] = 0;
                    for (int vertex : dependent_graph[taken_node]) depCnt[vertex]++;
                    cur_cost -= cost_graph[cur_node][taken_node].weight;
                    cur_solution.pop_back();
                    hungarian_solver.undue_row(cur_node, taken_node);
                    hungarian_solver.undue_column(taken_node, cur_node);
                    p_count++;
                }
                Split = false;
                //cout << "finish taking from GPQ with thread id = " << std::this_thread::get_id() << endl;
            }
        }

        steal_counter++;

        if (Steal_Opt == "ASYNC") {
            if (!Split || steal_counter >= active_thread - 1) {
                Split = false;
                std::unique_lock<std::mutex> idel_lck(Split_lock);
                Idel.notify_all();
                idel_lck.unlock();
                //cout << "notified!\n";
            }
        }

        else if (Steal_Opt == "SYNC") {
            if (!Split || steal_counter >= active_thread - 1) {
                Split = false;
                std::unique_lock<std::mutex> unlock(unlock_mutex);
                std::unique_lock<std::mutex> idel_lck(Split_lock);
                thread_counter.notify_all();
                unlock.unlock();
                Idel.notify_all();
                idel_lck.unlock();
            }
            else {
                while (Split && steal_counter < active_thread - 1) {
                    asssign_mutex.unlock();
                    std::unique_lock<std::mutex> unlock(unlock_mutex);
                    thread_counter.wait(unlock);
                } 
                //cout << "unblocked!\n";
            }
            
        }
        asssign_mutex.unlock();
        
    }

    while(!ready_list.empty()) {
        //Take the choosen node;
        //start_time = chrono::high_resolution_clock::now();
       //Back track if ready list is empty;
        bound = ready_list.back().lb;
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
        full_solution = false;
        suffix.clear();
        suffix_cost = 0;
        previous_snode = 0;
        
        enumerate(i+1);

        for (int vertex : dependent_graph[taken_node]) depCnt[vertex]++;
        taken_arr[taken_node] = 0;
        
        if (bound != -1) assign_historytable(cur_cost,bound,i);
        else {
            bound = hungarian_solver.get_matching_cost()/2;
            assign_historytable(cur_cost,bound,i);
        }
        
        if ((int)cur_solution.size() == node_count) {
            full_solution = true;
            previous_snode = cur_solution.back();
        }
        if (full_solution) {
            suffix.push_back(cur_solution.back());
            suffix_cost += cost_graph[cur_solution.back()][previous_snode].weight;
            previous_snode = cur_solution.back();
        }
        cur_solution.pop_back();
        cur_cost -= cost_graph[u][v].weight;
        if (cur_solution.size() >= 2) {
            hungarian_solver.undue_row(u,v);
		    hungarian_solver.undue_column(v,u);
        }
    }
    

    int size = 0;

    if (i == initial_depth && active_thread > 1) {
        active_thread--;

        if (Steal_Opt == "SYNC") {
            std::unique_lock<std::mutex> lock(unlock_mutex);
            thread_counter.notify_all();
            lock.unlock();
        }
    
        if (active_thread == 1) {
            Split = false;
            std::unique_lock<std::mutex> idel_lck(Split_lock);
            Idel.notify_all();
            idel_lck.unlock();
        }

        Split_Call.lock();
        GPQ_lock.lock();
        size = GPQ.size();
        if (size != 0) {
            *this = GPQ.front();
            GPQ.pop();
        }
        GPQ_lock.unlock();

        if (size == 0 && active_thread > 1) {
            Split = true;
            steal_counter = 1;
            while (Split) {
                std::unique_lock<std::mutex> idel_lck(Split_lock);
                Idel.wait(idel_lck);
            }

            GPQ_lock.lock();
            size = GPQ.size();
            if (size != 0) {
                *this = GPQ.front();
                GPQ.pop();
            }
            GPQ_lock.unlock();
        }
        Split_Call.unlock();
    }


    if (size != 0) {
        active_thread++;
        enumerate(initial_depth);
    }

    return;
}


void solver::solve_parallel(int thread_num, int pool_size) {
    vector<solver> Temp;
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
    
    //Initial filling of the GPQ
    for (auto node : ready_list) {
        auto target = *this;
        int taken_node = node.n;
        int cur_node = target.cur_solution.back();
        target.taken_arr[taken_node] = 1;
        for (int vertex : dependent_graph[taken_node]) target.depCnt[vertex]--;
        target.cur_cost += cost_graph[cur_node][taken_node].weight;
        target.cur_solution.push_back(taken_node);
        target.hungarian_solver.fix_row(cur_node, taken_node);
        target.hungarian_solver.fix_column(taken_node, cur_node);
        target.hungarian_solver.solve_dynamic();
        if (GPQ.size() < pool_size) {
            GPQ.push(target);
        }
        else {
            Temp.push_back(target);
        }
    }

    for (int i = thread_num - 1; i >= 0; i--) ready_thread.push_back(i);

    //While GPQ is not empty do split operation or assign threads with new node.

    

    while (true) {
        GPQ_lock.lock();
        if (!GPQ.size()) {
            GPQ_lock.unlock();
            break;
        }

        if (GPQ.size() <= thread_num) {
            auto target = GPQ.front();
            GPQ.pop();
            ready_list.clear();
            for (int i = node_count-1; i >= 0; i--) {
                if (!target.depCnt[i] && !target.taken_arr[i]) ready_list.push_back(node(i,-1));
            }
            for (auto node : ready_list) {
                int vertex = node.n;
                //Push split node back into GPQ
                if (!target.depCnt[vertex] && !target.taken_arr[vertex]) {
                    int taken_node = vertex;
                    int cur_node = target.cur_solution.back();
                    target.taken_arr[taken_node] = 1;
                    for (int vertex : target.dependent_graph[taken_node]) target.depCnt[vertex]--;
                    target.cur_cost += cost_graph[cur_node][taken_node].weight;
                    target.cur_solution.push_back(taken_node);
                    target.hungarian_solver.fix_row(cur_node, taken_node);
                    target.hungarian_solver.fix_column(taken_node, cur_node);
                    target.hungarian_solver.solve_dynamic();
                    
                    GPQ.push(target);

                    target.taken_arr[taken_node] = 0;
                    for (int vertex : target.dependent_graph[taken_node]) target.depCnt[vertex]++;
                    target.cur_cost -= cost_graph[cur_node][taken_node].weight;
                    target.cur_solution.pop_back();
                    target.hungarian_solver.undue_row(cur_node, taken_node);
                    target.hungarian_solver.undue_column(taken_node, cur_node);
                }
            }
            if (GPQ.empty()) {
                GPQ.push(target);
            }
        }
        
        for (int i = 0; i < thread_num; i++) {
            if (!GPQ.empty()) {
                auto subproblem = GPQ.front();
                GPQ.pop();
                int k = subproblem.cur_solution.size() - 1;
                solvers[i] = subproblem;
                solvers[i].initial_depth = k;
                active_thread++;
                Thread_manager[i] = thread(&solver::enumerate,move(solvers[i]),k);
            }
        }
        
        GPQ_lock.unlock();
        
        for (int i = 0; i < thread_num; i++) {
            if (Thread_manager[i].joinable()) {
                //cout << "waiting for thread to be joined" << endl;
                Thread_manager[i].join();
                ready_thread.push_back(i);
            }
        }
        
        GPQ_lock.lock();
        if (GPQ.empty()) {
            while (!Temp.empty()) {
                if (GPQ.size() < pool_size) {
                    GPQ.push(Temp.back());
                    Temp.pop_back();
                }
                else break;
            }
        }
        GPQ_lock.unlock();
    }
    
    
    return;
}

void solver::solve(string filename,string enum_opt,long time_limit,int pool_size,int thread_num,int split_num,string split_option,string steal_option) {
    retrieve_input(filename);
    split_depth = split_num;
    enum_option = enum_opt;
    thread_total = thread_num;
    p_size = pool_size;
    Split_Opt = split_option;
    Steal_Opt = steal_option;
    t_limit = time_limit * 1000000;
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
    cout << "Edge-based LB is " << EGB_static_lowerbound << endl;
    cout << "MMCP-based LB is " << MMCP_static_lowerbound << endl;
    */
    cur_cost = 0;

    auto start_time = chrono::high_resolution_clock::now();
    solve_parallel(thread_num,pool_size);
    auto end_time = chrono::high_resolution_clock::now();
    auto total_time = chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    cout << enum_opt << ": " << best_cost << "," << setprecision(4) << total_time / (float)(1000000) << endl;

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
    transitive_redundantcy();
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

void solver::print_dep() {
    for (int i = 0; i < node_count; i++) {
        cout << "node " << i << "has children: ";
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            if (k != dependent_graph[i].size() - 1) cout << dependent_graph[i][k] << ",";
            else cout << dependent_graph[i][k];
        }
        cout << endl;
    }
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