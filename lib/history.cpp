#include "history.hpp"
#include <limits>
#include <vector>

HistoryNode::HistoryNode(int prefix_cost, int suffix_cost, int bound, vector<int>suffix_sched){
	this->prefix_cost = prefix_cost;
	this->suffix_cost = suffix_cost;
	this->lower_bound = bound;
    this->suffix_sched = suffix_sched;
}


HistoryNode::HistoryNode(int prefix_cost, int bound){
	this->prefix_cost = prefix_cost;
	this->lower_bound = bound;
	this->suffix_cost = 0;
}

HistoryNode::HistoryNode(){
    prefix_cost = -1;
	suffix_cost = -1;
	lower_bound = -1;
}

