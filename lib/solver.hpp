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
        bool his_append;
        node(int x, int y): n{x},lb{y} {}
};

class solver {
    private:
        vector<solver> *local_pool = NULL;
        vector<int> depCnt;
        vector<int> taken_arr;
        int node_count = 0;
        //int EGB_static_lowerbound = 0;
        int MMCP_static_lowerbound = 0;
        int cur_cost = 0;
        int initial_depth = 0;
        int thread_id = -1;

        //Assignment Herustic Value
        bool full_solution = false;
        int suffix_cost = 0;
        int lower_depth = 0;
        
        vector<int> cur_solution;
        Hungarian hungarian_solver;
    public:
        int load_info = -1;
        vector<vector<int>> get_cost_matrix(int max_edge_weight);
        vector<int> nearest_neightbor(vector<int>* partial_solution,int* initial_cost);
        vector<int> tour_improvement(vector<int> initial_solution,int initial_cost,int initial_depth,int* final_cost);
        vector<int> roll_out();
        bool Wlkload_Request(int i);
        bool HistoryUtilization(pair<string,int>* key,int* lowerbound,bool* found);
        bool check_satisfiablity(int* local_cost, vector<int>* tour);
        int get_maxedgeweight();
        int dynamic_hungarian(int src, int dest);
        int enumerate(int i, int lb);
        int Two_Opt(vector<int> initial_tour);
        void Check_And_Distribute_Wlkload();
        void push_to_historytable(pair<string,int> key,int lower_bound,int i);
        void assign_parameter(vector<string> setting);
        void check_workload_request(int i);
        void notify_finished();
        void take_node(int taken_node);
        void untake_node(int src, int dest, int taken_node);
        void Check_Wkload_Request();
        void assign_workload(int taken_n, int lb);
        void push_to_global_pool();
        void append_suffix();
        void solve(string filename);
        void solve_parallel(int thread_num, int pool_size);
        void retrieve_input(string filename);
        void transitive_redundantcy();
        void sort_weight(vector<vector<edge>>& graph);
};
