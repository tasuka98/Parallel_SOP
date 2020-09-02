#include <map>
#include <string>
#include <vector>
#include "../lib/hungarian.h"

using namespace std;

class edge {
    private:
        int src;
        int dest;
        int weight;
    public:
        edge(int x, int y, int z): src{x},dest{y},weight{z} {}
        int retrieve_src();
        int retrieve_dest();
        int retrieve_weight();
};


class solver {
    private:
        int node_count;
        int EGB_static_lowerbound;
        int MMCP_static_lowerbound;
        int cur_cost;
        int best_cost;
        vector<int> best_solution;
        vector<int> cur_solution;
        enum stat {access,no_access};
        vector<vector<int>> in_degree;
        vector<vector<int>> dependent_graph;
        vector<vector<edge>> hung_graph;
        vector<vector<edge>> cost_graph;
        Hungarian hungarian_solver;
    public:
        vector<vector<int>> get_cost_matrix(int max_edge_weight);
        vector<int> nearest_neightbor();
        bool LB_Check(int src, int dest);
        int get_maxedgeweight();
        int mmcp_lb();
        int get_nodecount();
        int dynamic_hungarian(int src, int dest);
        void process_solution();
        void enumerate(int i);
        void solve(string filename);
        void retrieve_input(string filename);
        void transitive_redundantcy();
        void sort_weight(vector<vector<edge>>& graph);
        void print_dep();
        
};