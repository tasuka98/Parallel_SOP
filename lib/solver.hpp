#include <map>
#include <string>
#include <vector>

using namespace std;

class solver {
    private:
        int best_solution;
        enum stat {access,no_access};
        map<int,vector<int>> in_degree;
        map<int,vector<int>> cost_graph;
        map<int,vector<int>> dependent_graph;
    public:
        bool compare(const int& src, const int& target);
        int get_cost(int row, int col);
        void solve(string filename);
        void retrieve_input(string filename);
        void transitive_redundantcy();
        void sort_weight();
        void nearest_neightbor();
        void print_dep();
};