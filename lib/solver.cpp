#include "solver.hpp"
#include <iostream>
#include <map>
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

bool optimal_found = false;
int graph_index = 0;

vector<int> raedy_list;
int* depCnt;
int* taken_arr;

using namespace std;

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

bool solver::LB_Check(int src, int dest) {
    if (cur_cost > best_cost) return false;

    //Dyanmic hungarian lower bound
    int lower_bound = dynamic_hungarian(src,dest);

    if (lower_bound >= best_cost) {
        return false;
    }

    return true;
}


void solver::enumerate(int i) {
    if (optimal_found) return;
    if (i == node_count) {
        process_solution();
        return;
    }

    bool keep_explore = true;
    
    if (cur_solution.size() >= 2) {
        int u = cur_solution.end()[-2];
        int v = cur_solution.back();
        keep_explore = LB_Check(u,v);
    }
    

    if (!keep_explore) return;

    vector<int> ready_list;

    for (int i = 0; i < node_count; i++) {
        if (!depCnt[i] && !taken_arr[i]) {
            //Push vertices with 0 depCnt into the ready list
            ready_list.push_back(i);
        }
    }

    //update vertex dependent count for current node neighbors

    int taken_node = 0;
    int u = 0;
    int v = 0;

    while(!ready_list.empty()) {
        //Take the choosen node;
        taken_node = ready_list.back();
        ready_list.pop_back();
        if (cur_solution.empty()) u = -1;
        else {
            u = cur_solution.back();
            v = taken_node;
            cur_cost += cost_graph[u][v].retrieve_weight();
        }

        for (int vertex : dependent_graph[taken_node]) depCnt[vertex]--;
        cur_solution.push_back(taken_node);
        taken_arr[taken_node] = 1;

        enumerate(i+1);

        //Untake the choosen node;
        for (int vertex : dependent_graph[taken_node]) depCnt[vertex]++;
        taken_arr[taken_node] = 0;
        cur_solution.pop_back();
        if (u != -1) {
            cur_cost -= cost_graph[u][v].retrieve_weight();
            //Revert previous hungarian operations
            hungarian_solver.undue_row(u,v);
		    hungarian_solver.undue_column(v,u);
        }
    }
    
    return;
}


void solver::solve(string filename) {
    retrieve_input(filename);
    best_solution = nearest_neightbor();
    int max_edge_weight = get_maxedgeweight();
    hungarian_solver = Hungarian(node_count, max_edge_weight, get_cost_matrix(max_edge_weight));
    MMCP_static_lowerbound = hungarian_solver.start()/2;
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

    cout << "best solution found using NN is " << best_cost << endl;
    cout << "Edge-based LB is " << EGB_static_lowerbound << endl;
    cout << "MMCP-based LB is " << MMCP_static_lowerbound << endl;

    cur_cost = 0;
    enumerate(0);

    cout << "optimal solution using B&B is " << best_cost << endl;

    cout << "the optimal solution contains ";
    for (int i = 0; i < node_count; i++) {
        if (i != node_count - 1) cout << best_solution[i] << "-->";
        else cout << best_solution[i];
    }
    cout << endl;
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
                dependent_graph[j].push_back(i);
            }
            else { 
                cost_graph[i].push_back(edge(i,j,edge_weight));
                if (i != j) hung_graph[i].push_back(edge(i,j,edge_weight));
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
    stack<int> dfs_stack;
    stat visit_arr[node_count];

    for (int i = 0; i < node_count; i++) {
        visit_arr[i] = stat::access;
    }

    for (int i = 0; i < node_count; i++) {
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            dfs_stack.push(dependent_graph[i][k]);
            visit_arr[dependent_graph[i][k]] = stat::no_access;
        }

        //Use DFS on each and every node to fix transitive redundantcy.
        while (!dfs_stack.empty()) {
            int node_num = dfs_stack.top();
            dfs_stack.pop();

            if (visit_arr[node_num] == stat::no_access) visit_arr[node_num] = stat::access;
            // Remove redundent edges in the graph
            else if (find(dependent_graph[i].begin(),dependent_graph[i].end(),node_num) != dependent_graph[i].end()) {
                    for (auto k = dependent_graph[i].begin(); k != dependent_graph[i].end(); k++) {
                        if (*k == node_num) {
                            dependent_graph[i].erase(k);
                            break;
                        }
                    }
            }

            for (long unsigned int k = 0; k < dependent_graph[node_num].size(); k++) {
                dfs_stack.push(dependent_graph[node_num][k]);
            }
        }
        //Clear visit_array for the iteration of the dfs.
        for (int i = 0; i < node_count; i++) {
            visit_arr[i] = stat::access;
        }
    }

    in_degree = std::vector<vector<int>>(node_count);

    for (int i = 0; i < node_count; i++) {
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            int c = dependent_graph[i][k];
            in_degree[c].push_back(i);
        }
    }
}


bool compare (edge& src, edge& target) {
    int src_weight = src.retrieve_weight();
    int dest_weight = target.retrieve_weight();
    return (src_weight < dest_weight);
}

//Sort cost-graph weight by ascending order
void solver::sort_weight(vector<vector<edge>>& graph) {
    int size = graph.size();

    for (int i = 0; i < size; i++) {
        graph_index = i;
        sort(graph[i].begin(),graph[i].end(),compare);
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
             
            if (!visit_arr[node.retrieve_dest()] && !depCnt_arr[node.retrieve_dest()]) {
                current_node = node.retrieve_dest();
                solution_cost += node.retrieve_weight();
                solution.push_back(current_node);
                num++;
                visit_arr[node.retrieve_dest()] = true;
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
            int weight = edge_weight.retrieve_weight();
            if (weight > max) max = weight;
        }
    }
    return max;
}

int solver::get_nodecount() {
    return node_count;
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
    
    memset(depCnt_arr,0,node_count*sizeof(int));
    for (int i = 0; i < node_count; i++) {
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            depCnt_arr[dependent_graph[i][k]]++;
        }
    }

    int out_max = 0;
    int in_max = 0;

    for (int i = 0; i < node_count; i++) {
        int min = INT_MAX;


        for (auto dest : dependent_graph[i]) {
            int weight = cost_graph[i][dest].retrieve_weight();
            if (weight < min) min = weight;
        }

        for (int k = 0; k < node_count; k++) {
            int weight = cost_graph[i][k].retrieve_weight();
            if (find(dependent_graph[k].begin(),dependent_graph[k].end(),i) == dependent_graph[k].end()) {
                if (k != i && depCnt_arr[k] == 0 && weight < min) min = weight;
            }
        }

        out_array[i] = min;

        min = INT_MAX;
        for (auto dest : in_degree[i]) {
            int weight = cost_graph[dest][i].retrieve_weight();
            if (weight < min) min = weight;
        }

        for (int k = 0; k < node_count; k++) {
            int weight = cost_graph[k][i].retrieve_weight();
            if (find(in_degree[k].begin(),in_degree[k].end(),i) == dependent_graph[k].end()) {
                if (k != i && depCnt_arr[k] == 0 && weight < min) min = weight;
            }
        }

        in_array[i] = min;
    }

    for (int i = 0; i < node_count; i++) {
        if (!dependent_graph[i].size() && out_array[i] > out_max) out_max = out_array[i];
        if (!depCnt_arr[i] && in_array[i] > in_max) in_max = in_array[i];
        outsum += out_array[i];
        insum += in_array[i];
    }

    return max(outsum-out_max,insum-in_max);
}

vector<vector<int>> solver::get_cost_matrix(int max_edge_weight) {
    vector<vector<int>> matrix(cost_graph.size());
    for(int i = 0; i < node_count; ++i){
		matrix[i] = vector<int>(node_count, max_edge_weight*2);
	}

    int i = 0;
    int k = 0;
    for (vector<edge> edge_list : hung_graph) {
        k = 0;
        for (auto edge : edge_list) {
            matrix[edge.retrieve_src()][edge.retrieve_dest()] = edge.retrieve_weight()*2;
            k++;
        }
        i++;
    }

    /*
    for (int i = 0; i < size; i++) {
		std::cout << "matrix at row" << i << "is ";
		for (auto weight : matrix[i]) {
			std::cout << weight << ",";
		}
		std::cout << std::endl;
	}
    */
    return matrix;
}

int edge::retrieve_weight() {
    return weight;
}

int edge::retrieve_src() {
    return src;
}

int edge::retrieve_dest() {
    return dest;
}
    