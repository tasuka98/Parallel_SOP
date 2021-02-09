#ifndef SOLVER_H
#define SOLVER_H

#include <map>
#include <string>
#include <vector>
#include <queue>
#include <chrono>
#include <climits>
#include "../lib/hungarian.hpp"
#include "history.hpp"
#include "load_test.hpp"

using namespace std;

class load_stats {
    public:
        bool out_of_work;
        bool temp_disable;
        int padding[14];
        int load;
        load_stats() {
            out_of_work = false;
            temp_disable = false;
            load = 0;
        }
};

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
        int n = -1;
        int lb = -1;
        int nc = -1;
        bool pushed = false;;
        node(int x, int y): n{x},lb{y} {}
};

struct sop_state {
    vector<int> depCnt;
    vector<int> taken_arr;
    vector<int> cur_solution;
    Hungarian hungarian_solver;
    int cur_cost = 0;
    int initial_depth = 0;
    int suffix_cost = 0;
    int originate = -1;
    int load_info = -1;
    bool full_solution = false;
};

class solver {
    private:
        //////////////////////////////////////////////////////////////////////
        //                Used for workload distribution test
        //bool Stolen = false;
        //load_test Stolen_load_info;
        //////////////////////////////////////////////////////////////////////
        //Local memory block
        
        //Time frame to prevent double stealing;
        std::chrono::time_point<std::chrono::system_clock> time_frame_start;
        bool stolen = false;
        bool abandon_work = false;
        bool grabbed = false;
        int enumerated_nodes = 0;
        int node_count = 0;
        deque<sop_state> *local_pool = NULL;
        sop_state problem_state;
        sop_state reserve_state;
        //int EGB_static_lowerbound = 0;
        int thread_id = -1;
    public:
        HistoryNode* retrieve_his_node();
        vector<vector<int>> get_cost_matrix(int max_edge_weight);
        vector<int> nearest_neightbor(vector<int>* partial_solution,int* initial_cost);
        vector<int> roll_out();
        bool Wlkload_Request(int i);
        bool HistoryUtilization(pair<vector<bool>,int>* key,int* lowerbound,bool* found,int cost);
        //bool HistoryUtilization(pair<vector<bool>,int>* key,int* lowerbound,bool* found);
        bool check_satisfiablity(int* local_cost, vector<int>* tour);
        bool Split_local_pool();
        bool push_to_global_pool(unsigned size);
        bool Split_level_check();
        bool Grab_from_GPQ(bool reserve);
        bool Steal_Workload();
        int get_maxedgeweight();
        int tour_improvement(vector<int> initial_solution,int initial_cost,int initial_depth);
        int dynamic_hungarian(int src, int dest);
        int enumerate(int i);
        //void enumerate(int i);
        void Direct_Workload();
        void assign_thread_load();
        void Check_And_Distribute_Wlkload();
        void push_to_historytable(pair<vector<bool>,int> key,int lower_bound);
        void assign_parameter(vector<string> setting);
        void check_workload_request(int i);
        void notify_finished();
        void assign_workload(int taken_n, int lb, const string dest);
        void solve(string filename,int thread_num,string assignment_scheme);
        void solve_parallel(int thread_num, int pool_size);
        void retrieve_input(string filename);
        void transitive_redundantcy();
        void sort_weight(vector<vector<edge>>& graph);
};

#endif
