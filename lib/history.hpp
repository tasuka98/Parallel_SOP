#include <iostream>
#include <vector>

using namespace std;

#pragma once

class HistoryNode {
	public:
        HistoryNode(int prefix_cost, int bound);
        HistoryNode();
		int prefix_cost;
		int lower_bound;
};