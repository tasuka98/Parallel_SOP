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
#include <bits/stdc++.h>
#include <unordered_map>
#include<setjmp.h>

using namespace std;

Hash_Map history_table(12582917);
bool optimal_found = false;
int graph_index = 0;
int enumerated_nodes = 0;
int calculated_bounds = 0;
int recently_added = 0;
long eclipsed_time = 0;

HistoryNode visted;
bool full_solution = false;
float total_time_static = 0;
float ready_list_time = 0;
int lb = 0;
int suffix_cost = 0;
int previous_snode = 0;

std::chrono::time_point<std::chrono::system_clock> start_time_limit;
jmp_buf buf;
vector<int> raedy_list;
vector<int> suffix;
vector<bool> picked_list;

int* depCnt;
int* taken_arr;

void solver::push_to_historytable(pair<string,int> key,int lower_bound,int i) {
    HistoryNode node;
    node.prefix_cost = cur_cost;
    if (full_solution) node.lower_bound = suffix_cost;
    else node.lower_bound = lower_bound;

    if (history_table.get_cur_size() < 0.8 * history_table.get_max_size()) {
        history_table.insert(key,node);
    }
    else if (i < int(0.5 * node_count)) {
        history_table.insert(key,node);
    }
    return;
}

int solver::dynamic_hungarian(int src, int dest) {
    hungarian_solver.fix_row(src, dest);
	hungarian_solver.fix_column(dest, src);
	hungarian_solver.solve_dynamic();
    return hungarian_solver.get_matching_cost()/2;
}

bool solver::HistoryUtilization(pair<string,int>* key,int* lowerbound,bool* found) {
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

    history_table.insert(*key,history_node);
    return true;
}

int solver::dynamic_edb() {
    int picked_node = cur_solution.back();
    picked_list[picked_node] = true;
    recently_added = picked_node;
    return mmcp_lb();
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
    enumerated_nodes++;
    
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

    string bit_string(node_count, '0');
    for (auto node : cur_solution) {
        bit_string[node] = '1';
    }   
    int last_element = cur_solution.back();
    auto key = make_pair(bit_string,last_element);

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
            if (cur_solution.size() == (unsigned)node_count) {
                if (cur_cost < best_cost) {
                    best_solution = cur_solution;
                    best_cost = cur_cost;
                }
                cur_solution.pop_back();
                cur_cost -= cost_graph[src][dest.n].weight;
                ready_list.erase(ready_list.begin()+i);
                i--;
                continue;
            }
            
            else {
                key.first[dest.n] = '1';
                key.second = dest.n;
                bool decision = HistoryUtilization(&key,&temp_lb,&taken);

                if (!taken) {
                    temp_lb = dynamic_hungarian(src,dest.n);
                    push_to_historytable(key,temp_lb,i);
                    hungarian_solver.undue_row(src,dest.n);
                    hungarian_solver.undue_column(dest.n,src);
                }
                else if (taken && !decision) {
                    cur_solution.pop_back();
                    cur_cost -= cost_graph[src][dest.n].weight;
                    key.first[dest.n] = '0';
                    key.second = last_element;
                    ready_list.erase(ready_list.begin()+i);
                    i--;
                    continue;
                }
                if (temp_lb >= best_cost) {
                    cur_solution.pop_back();
                    cur_cost -= cost_graph[src][dest.n].weight;
                    key.first[dest.n] = '0';
                    key.second = last_element;
                    ready_list.erase(ready_list.begin()+i);
                    i--;
                    continue;
                }
                cur_solution.pop_back();
                cur_cost -= cost_graph[src][dest.n].weight;
                ready_list[i].nc = cost_graph[src][dest.n].weight;
                ready_list[i].lb = temp_lb;
                key.first[dest.n] = '0';
                key.second = last_element;
            }
        }
        if (enum_option == "DH") sort(ready_list.begin(),ready_list.end(),bound_sort);
        else if (enum_option == "NN") sort(ready_list.begin(),ready_list.end(),nearest_sort);
    }

    while(!ready_list.empty()) {
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
        suffix_cost = 0;
        ///////
        
        //cout << "current solution is";
        //for (auto node : cur_solution) cout << node << ",";
        //cout <<  endl;
        
        /////

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
            longjmp(buf, 1);
        }
        
       // cout << "assign to history table time: " << setprecision(4) << total_time / (float)(1000000) << endl;
    }
    return;
}


void solver::solve(string filename,string enum_opt,long time_limit) {
    retrieve_input(filename);
    history_table.set_node_t(node_count);
    enum_option = enum_opt;
    t_limit = time_limit;
    best_solution = nearest_neightbor();
    int max_edge_weight = get_maxedgeweight();
    hungarian_solver = Hungarian(node_count, max_edge_weight+1, get_cost_matrix(max_edge_weight+1));
    MMCP_static_lowerbound = hungarian_solver.start()/2;
    picked_list = vector<bool>(node_count,false);
    EGB_static_lowerbound = mmcp_lb();
    depCnt = new int[node_count];
    taken_arr = new int[node_count];
    memset(depCnt,0,node_count * sizeof(int));
    memset(taken_arr,0,node_count * sizeof(int));

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
    start_time_limit = std::chrono::system_clock::now();
    if (setjmp(buf)) {
        cout << "Timed out!" << endl;
    }
    else {
        cur_solution.push_back(0);
        taken_arr[0] = 1;
        for (int vertex : dependent_graph[0]) depCnt[vertex]--;
        enumerate(1);
    }
    auto end_time = chrono::high_resolution_clock::now();
    auto total_time = chrono::duration_cast<std::chrono::microseconds>( end_time - start_time ).count();

    cout << enum_opt << ": " << best_cost << "," << setprecision(4) << total_time / (float)(1000000) << endl;

    /*
    cout << "the optimal solution contains ";
    for (int i = 0; i < node_count; i++) {
        if (i != node_count - 1) cout << best_solution[i] << "-->";
        else cout << best_solution[i];
    }
    cout << endl;
    */

    //cout << "Total enumerated nodes are " << enumerated_nodes << endl;
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

    for(int i = 0; i < node_count; ++i){
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
        graph_index = i;
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


int solver::mmcp_lb() {
    int outsum = 0;
    int insum = 0;
    int out_array[node_count];
    int in_array[node_count];
    int depCnt_arr[node_count];
    bool trim_in = true;
    bool trim_out = true;
    
    memset(depCnt_arr,0,node_count*sizeof(int));
    for (int i = 0; i < node_count; i++) {
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            depCnt_arr[dependent_graph[i][k]]++;
        }
    }

    int out_max = 0;
    int in_max = 0;

    out_array[0] = 0;
    out_array[node_count-1] = 0;
    in_array[0] = 0;
    in_array[node_count-1] = 0;


    for (int i = 1; i < node_count - 1; i++) {
        int min_out = INT_MAX;
        int min_in = INT_MAX;
        bool calculate_in = true;
        bool calculate_out = true;


        if (picked_list[i]) {
            calculate_in = false;
            trim_in = false;
            if (recently_added != i) {
                calculate_out = false;
                trim_out = false;
            }
        }


        for (int k = 1; k < node_count - 1; k++) {
            int weight = cost_graph[i][k].weight;
            if (calculate_out) {
                if (find(in_degree[i].begin(),in_degree[i].end(),edge(k,i,-1)) == in_degree[i].end()) {
                    if (k != i && weight < min_out && !picked_list[k]) min_out = weight;
                }
            }
            
            if (calculate_in) {
                if (find(in_degree[k].begin(),in_degree[k].end(),edge(i,k,-1)) == in_degree[k].end()) {
                    if (recently_added == k && k != i && weight < min_in) {
                        min_in = weight;
                    }
                    else if (!picked_list[k] && k != i && weight < min_in) {
                        min_in = weight;
                    }
                }
            }
            
        }

        if (min_out == INT_MAX) {
            trim_out = false;
            min_out = 0;
        }

        if (min_in == INT_MAX) {
            trim_in = false;
            min_out = 0;
        }

        
        if (calculate_in) in_array[i] = min_in;
        else in_array[i] = 0;

        if (calculate_out) out_array[i] = min_out;
        else out_array[i] = 0;
        
    }

    for (int i = 1; i < node_count - 1; i++) {
        if (dependent_graph[i].size() == 1 && dependent_graph[i][0] == node_count-1 && out_array[i] > out_max) {
            out_max = out_array[i];
        }
        if (depCnt_arr[i] == 1 && in_degree[i][0].src == 0 && in_array[i] > in_max) {
            in_max = in_array[i];
        }
        outsum += out_array[i];
        insum += in_array[i];
    }

    if (trim_in == false && trim_out == false) return max(outsum,insum);
    else if (trim_in == false) return max(outsum-out_max,insum);
    else if (trim_out == false) return max(outsum,insum-in_max);
    
    return max(outsum-out_max,insum-in_max);
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