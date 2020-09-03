#include "history.hpp"
#include <limits>
#include <vector>

HistoryNode::HistoryNode(int prefix_cost, int suffix_cost, int bound, vector<int>prefix_sched, vector<int>suffix_sched){
	this->prefix_cost = prefix_cost;
	this->suffix_cost = suffix_cost;
	this->lower_bound = bound;
    this->prefix_sched = prefix_sched;
    this->suffix_sched = suffix_sched;
}


HistoryNode::HistoryNode(int prefix_cost, int bound, vector<int>prefix_sched){
	this->prefix_cost = prefix_cost;
	this->lower_bound = bound;
    this->prefix_sched = prefix_sched;
	this->suffix_cost = 0;
}

HistoryNode::HistoryNode(){
    prefix_cost = std::numeric_limits<int>::max();
	suffix_cost = std::numeric_limits<int>::max();
	lower_bound = std::numeric_limits<int>::max();
}

