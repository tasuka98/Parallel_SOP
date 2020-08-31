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
    cout << "best solution found using NN is " << best_solution << endl;
}

void solver::retrieve_input(string filename) {
    ifstream inFile;
    string line;
    inFile.open(filename);
    int c = 0;

    // Read input files and store it inside an array.
    while (getline(inFile,line)) {  
        stringstream sstream;
        sstream << line;
        string weight;
        int weight_num;
        int k = 0;
        while (sstream >> weight) {
            stringstream(weight) >> weight_num;
            if (cost_graph.find(c) == cost_graph.end()) {
                dependent_graph.insert(pair<int,vector<int>>(c,{}));
                cost_graph.insert(pair<int,vector<int>>(c,{weight_num}));
            }
            else {
                cost_graph[c].push_back(weight_num);
            }
            if (weight_num == -1) {
                dependent_graph[c].push_back(k);
            }
            k++;
        }
        c++;
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

bool solver::compare (const int& src, const int& target) {
    return (cost_graph[src][graph_index] < cost_graph[target][graph_index]);
}

int solver::get_cost(int row, int col) {
    return cost_graph[col][row];
}

void solver::sort_weight() {
    int dep_size = dependent_graph.size();

    for (int i = 0; i < dep_size; i++) {
        graph_index = i;
        sort(dependent_graph[i].begin(),dependent_graph[i].end(),bind(&solver::compare, this, placeholders::_1, placeholders::_2));
    }

    for (int i = 0; i < dep_size; i++) {
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            int c = dependent_graph[i][k];
            if (in_degree.find(c) == in_degree.end()) {
                in_degree.insert(pair<int,vector<int>>(c,{i}));
            }
            else {
                in_degree[c].push_back(i);
            }
        }
    }

    return;
}

void solver::nearest_neightbor() {
    
    vector<int> solution;
    int current_node;
    int dep_size = dependent_graph.size();
    bool visit_arr[dep_size];
    int depCnt_arr[dep_size];
    memset(depCnt_arr,0,dep_size*sizeof(int));

    for (int i = 0; i < dep_size; i++) {
        visit_arr[i] = false;
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            depCnt_arr[dependent_graph[i][k]]++;
        }
        if (!depCnt_arr[i]) {
            current_node = i;
        }
    }
    
    solution.push_back(current_node);
    int num = 1;
    int solution_cost = 0;

    while (num < dep_size) {
        for (long unsigned int i = 0; i < dependent_graph[current_node].size(); i++) {
            depCnt_arr[dependent_graph[current_node][i]]--;
        }
        for (auto node: dependent_graph[current_node]) {
            if (!visit_arr[node] && !depCnt_arr[node]) {
                if (!in_degree[node].empty()) solution_cost += cost_graph[node][in_degree[node][0]];
                current_node = node;
                solution.push_back(current_node);
                num++;
                visit_arr[node] = true;
            }
        }
    }

    for (long unsigned int i = 0; i < solution.size(); i++) {
        cout << solution[i] << ",";
    }
    cout << endl;

    best_solution = solution_cost;
}


void solver::print_dep() {
    
    for (long unsigned int i = 0; i < dependent_graph.size(); i++) {
        cout << "node " << i << " has children " << ' ';
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            if (k != dependent_graph[i].size() - 1) cout << dependent_graph[i][k] << ",";
            else cout << dependent_graph[i][k];
        }
        cout << endl;
    }

}
    