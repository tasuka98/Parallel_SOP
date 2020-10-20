#include "history.hpp"
#include <limits>
#include <vector>

HistoryNode::HistoryNode(int prefix_cost, int bound){
	this->prefix_cost = prefix_cost;
	this->lower_bound = bound;
}

HistoryNode::HistoryNode(){
    prefix_cost = -1;
	lower_bound = -1;
}
