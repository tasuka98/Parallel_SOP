#include <map>
#include <string>
#include <vector>
#include <queue>
#include <chrono>
#include "../lib/hungarian.hpp"

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
        queue<solver> *local_pool = NULL;
        int total_lb = 0;
        int local_lb = -1;

        vector<int> depCnt;
        vector<int> taken_arr;
        bool full_solution = false;
        int node_count = 0;
        //int EGB_static_lowerbound = 0;
        int MMCP_static_lowerbound = 0;
        int cur_cost = 0;
        int initial_depth = 0;
        int thread_id = -1;

        vector<int> suffix;
        vector<int> cur_solution;
        Hungarian hungarian_solver;
    public:
        vector<vector<int>> get_cost_matrix(int max_edge_weight);
        vector<int> nearest_neightbor();
        bool HistoryUtilization(int* lowerbound,bool* found);
        bool check_tlimit();
        bool take_from_local();
        int get_maxedgeweight();
        int dynamic_hungarian(int src, int dest);
        void assign_parameter(vector<string> setting);
        void check_workload_request(int i);
        void notify_finished();
        void take_node(int taken_node);
        void untake_node(int src, int dest, int taken_node);
        void assign_workload(int taken_n, int lb);
        void push_to_global_pool();
        void assign_historytable(int prefix_cost,int lower_bound,int i);
        void enumerate(int i);
        void solve(string filename);
        void solve_parallel(int thread_num, int pool_size);
        void retrieve_input(string filename);
        void transitive_redundantcy();
        void sort_weight(vector<vector<edge>>& graph);
};