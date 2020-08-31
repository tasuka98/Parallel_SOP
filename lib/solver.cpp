#include "solver.hpp"
#include <iostream>
#include <map>
#include <fstream>
#include <stack>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <functional>
#include <vector>

int graph_index = 0;

using namespace std;

void solver::solve(string filename) {
    retrieve_input(filename);
    nearest_neightbor();
    //hungarian_solver = Hungarian(get_nodecount(), get_maxedgeweight(), get_cost_matrix());
    //static_lowerbound = hungarian_solver.start()/2;
    cout << "best solution found using NN is " << best_solution << endl;
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
    dependent_graph = vector<vector<int>>(size);

    for (int i = 0; i < (int)file_matrix.size(); i++) {
        int j = 0;
        for (auto edge_weight: file_matrix[i]) {
            if (edge_weight < 0) {
                cost_graph[i].push_back(edge(i,j,file_matrix[j][i]));
                dependent_graph[j].push_back(i);
            }
            else cost_graph[i].push_back(edge(i,j,edge_weight));
            j++;
        }
    }

    
    //Trim redundant edges
    transitive_redundantcy();
    
    //sort input based on weight for NN heurestic
    sort_weight();

    

    return;
}

void solver::transitive_redundantcy() {
    stack<int> dfs_stack;
    stat visit_arr[dependent_graph.size()];

    for (long unsigned int i = 0; i < dependent_graph.size(); i++) {
        visit_arr[i] = stat::access;
    }

    for (long unsigned int i = 0; i < dependent_graph.size(); i++) {
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
        for (long unsigned int i = 0; i < dependent_graph.size(); i++) {
            visit_arr[i] = stat::access;
        }
    }

    
}

bool compare (edge& src, edge& target) {
    int src_weight = src.retrieve_weight();
    int dest_weight = target.retrieve_weight();
    return (src_weight < dest_weight);
}

void solver::sort_weight() {
    int size = cost_graph.size();

    for (int i = 0; i < size; i++) {
        graph_index = i;
        sort(cost_graph[i].begin(),cost_graph[i].end(),compare);
    }

    in_degree = std::vector<vector<int>>(dependent_graph.size());

    for (int i = 0; i < size; i++) {
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            int c = dependent_graph[i][k];
            in_degree[c].push_back(i);
        }
    }

    return;
}

void solver::nearest_neightbor() {
    
    vector<int> solution;
    int current_node;
    int dep_size = dependent_graph.size();
    bool visit_arr[dep_size];
    bool selected = false;
    int depCnt_arr[dep_size];
    memset(depCnt_arr,0,dep_size*sizeof(int));

    for (int i = 0; i < dep_size; i++) {
        visit_arr[i] = false;
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            depCnt_arr[dependent_graph[i][k]]++;
        }
    }

    for (int i = 0; i < dep_size; i++) {
        if (depCnt_arr[i] == 0 && !selected) {
            current_node = i;
            visit_arr[current_node] = true;
            break;
        }
    }

    solution.push_back(current_node);
    
    int num = 1;
    int solution_cost = 0;

    while (num < dep_size) {
        for (long unsigned int i = 0; i < dependent_graph[current_node].size(); i++) {
            depCnt_arr[dependent_graph[current_node][i]]--;
        }
        for (auto node: cost_graph[current_node]) {
             
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

    

    for (long unsigned int i = 0; i < solution.size(); i++) {
        cout << solution[i] << ",";
    }
    cout << endl;

    best_solution = solution_cost;
}


int solver::get_maxedgeweight() {
    int max = 0;
    for (long unsigned int i = 0; i < cost_graph.size(); i++) {
        for (auto edge_weight : cost_graph[i]) {
            int weight = edge_weight.retrieve_weight();
            if (weight > max) max = weight;
        }
    }
    return max;
}

int solver::get_nodecount() {
    return (int)cost_graph.size();
}

void solver::print_dep() {
    for (long unsigned int i = 0; i < dependent_graph.size(); i++) {
        cout << "node " << i << "has children: ";
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            if (k != dependent_graph[i].size() - 1) cout << dependent_graph[i][k] << ",";
            else cout << dependent_graph[i][k];
        }
        cout << endl;
    }
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
    