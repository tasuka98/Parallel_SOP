#ifndef SOLVER_H
#define SOLVER_H

#include <map>
#include <string>
#include <vector>
#include <gmpxx.h>
#include <gmp.h>
#include <queue>
#include <chrono>
#include "../lib/hungarian.hpp"
#include "history.hpp"
#include "load_test.hpp"

using namespace std;

class load_stats {
    public:
        bool out_of_work;
        int load;
        load_stats() {
            out_of_work = false;
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
        int n;
        int lb;
        int nc;
        node(int x, int y): n{x},lb{y} {}
};

class solver {
    private:
        //////////////////////////////////////////////////////////////////////
        //                Used for workload distribution test
        //bool Stolen = false;
        //load_test Stolen_load_info;
        //////////////////////////////////////////////////////////////////////

        //Local memory block
        HistoryNode* history_block = NULL;
        unsigned counter = 0;

        //Time frame to prevent double stealing;
        std::chrono::time_point<std::chrono::system_clock> time_frame_start;
        bool stolen = false;

        bool Just_Stole = false;

        vector<solver> *local_pool = NULL;
        vector<int> depCnt;
        vector<int> taken_arr;
        int node_count = 0;

        //int EGB_static_lowerbound = 0;
        int cur_cost = 0;
        int initial_depth = 0;
        int thread_id = -1;

        //Assignment Herustic Value
        bool full_solution = false;
        int suffix_cost = 0;
        int force_cnt = 0;

        //Assign originate node value
        int originate = -1;
        
        vector<int> cur_solution;
        Hungarian hungarian_solver;
    public:
        int load_info = -1;
        
        HistoryNode* retrieve_his_node();
        mpf_class retrieve_impact();
        mpf_class get_load_estimation() const;
        vector<vector<int>> get_cost_matrix(int max_edge_weight);
        vector<int> nearest_neightbor(vector<int>* partial_solution,int* initial_cost);
        vector<int> roll_out();
        bool Wlkload_Request(int i);
        bool HistoryUtilization(pair<vector<bool>,int>* key,int* lowerbound,bool* found,int cost);
        //bool HistoryUtilization(pair<vector<bool>,int>* key,int* lowerbound,bool* found);
        bool check_satisfiablity(int* local_cost, vector<int>* tour);
        bool Split_local_pool();
        bool push_to_global_pool();
        int get_maxedgeweight();
        int tour_improvement(vector<int> initial_solution,int initial_cost,int initial_depth);
        int dynamic_hungarian(int src, int dest);
        int Get_cur_depth() const;
        int Get_cur_cost() const;
        double get_impact_iter(int i, int current_node, vector<double> &impact_arr, vector<bool> &visited_arr);
        void Direct_Workload();
        void print_dep();
        void enumerate(int i);
        void Check_And_Distribute_Wlkload();
        void push_to_historytable(pair<vector<bool>,int> key,int lower_bound, int i);
        void assign_parameter(vector<string> setting);
        void check_workload_request(int i);
        void notify_finished();
        void Check_Wkload_Request();
        void assign_workload(int taken_n, int lb);
        void solve(string filename,int thread_num);
        void solve_parallel(int thread_num, int pool_size);
        void retrieve_input(string filename);
        void transitive_redundantcy();
        void sort_weight(vector<vector<edge>>& graph);
};

#endif
