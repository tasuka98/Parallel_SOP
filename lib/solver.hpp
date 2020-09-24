#include <map>
#include <string>
#include <vector>
#include <chrono>
#include "../lib/hungarian.h"

using namespace std;

class edge {
    public:
        int src;
        int dest;
        int weight;
        edge(int x, int y, int z): src{x},dest{y},weight{z} {}
        bool operator==(const edge &rhs) const {return this->src == rhs.src;}
};

class node {
    public:
        int n;
        int lb;
        int nc;
        node(int x, int y): n{x},lb{y} {}
};


class solver {
    private:
        string enum_option;
        vector<int> depCnt;
        vector<int> taken_arr;
        int node_count;
        int EGB_static_lowerbound;
        int MMCP_static_lowerbound;
        int cur_cost;
        long t_limit;
        vector<int> cur_solution;
        vector<vector<int>> dependent_graph;
        vector<vector<edge>> in_degree;
        vector<vector<edge>> hung_graph;
        vector<vector<edge>> cost_graph;
        Hungarian hungarian_solver;
    public:
        vector<vector<int>> get_cost_matrix(int max_edge_weight);
        vector<int> nearest_neightbor();
        bool HistoryUtilization(int* lowerbound,bool* found,bool suffix_exist);
        int get_maxedgeweight();
        int dynamic_edb();
        int dynamic_hungarian(int src, int dest);
        int mmcp_lb();
        int History_LB();
        void process_solution();
        void assign_historytable(int prefix_cost,int lower_bound,int i);
        void enumerate(int i);
        void solve(string filename,string enum_opt,long time_limit,int pool_size, int thread_num);
        void solve_parallel(int thread_num, int pool_size);
        void retrieve_input(string filename);
        void transitive_redundantcy();
        void sort_weight(vector<vector<edge>>& graph);
        void print_dep();
};