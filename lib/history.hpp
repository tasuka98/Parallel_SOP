#include <iostream>
#include <vector>

using namespace std;

class HistoryNode {
	public:
        HistoryNode(int prefix_cost, int suffix_cost, int bound, vector<int>prefix_sched, vector<int>suffix_sched);
        HistoryNode(int prefix_cost, int bound, vector<int>prefix_sched);
        HistoryNode();
		int prefix_cost;
        int suffix_cost;
		int lower_bound;
        vector<int> prefix_sched;
        vector<int> suffix_sched;
};