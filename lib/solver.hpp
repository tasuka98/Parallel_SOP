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
        int static_lowerbound;
        int best_solution;
        enum stat {access,no_access};
        vector<vector<int>> in_degree;
        vector<vector<edge>> cost_graph;
        vector<vector<int>> dependent_graph;
        Hungarian hungarian_solver;
    public:
        vector<vector<int>> get_cost_matrix();
        int get_maxedgeweight();
        int get_nodecount();
        void solve(string filename);
        void retrieve_input(string filename);
        void transitive_redundantcy();
        void sort_weight();
        void nearest_neightbor();
        void print_dep();
};