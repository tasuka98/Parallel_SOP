#include <vector>
#include <unordered_map>
#include <iostream>
#include <list>
#include <mutex>
#include <cerrno>
#include <atomic>
#include <sys/sysinfo.h>
#include "history.hpp"

static vector<mutex> hash_lock;

class Hash_Map {
    private:
        atomic<size_t> cur_size;
        size_t size = 0;
        size_t max_size = 0;
        size_t node_size = 0;
        vector<list<pair<pair<vector<bool>,int>,HistoryNode*>>> History_table;
    public:
        atomic<long> num_of_waits;
        Hash_Map(int size);
        bool isEmpty();
        size_t get_max_size();
        size_t get_cur_size();
        uint32_t hash_func(pair<vector<bool>,int> item);
        void lock_table(int& key);
        void unlock_table(int& key);
        void increase_size(size_t size_incre);
        void set_node_t(int node_size);
        void insert(pair<vector<bool>,int>& item,HistoryNode* node);
        void average_size();
        HistoryNode* retrieve(pair<vector<bool>,int>* item,int key);
};

Hash_Map::Hash_Map(int size) {
    this->size = size;
    struct sysinfo info;
	if(sysinfo(&info) != 0){
        cout << "can't retrieve sys mem info\n";
		exit(1);
	}

    max_size = (double)info.freeram*0.85-size*sizeof(mutex);
    cur_size = 0;
    hash_lock = vector<mutex>(size);
    History_table.resize(size);
}

size_t Hash_Map::get_max_size() {
    return max_size;
}

size_t Hash_Map::get_cur_size() {
    return cur_size;
}

void Hash_Map::increase_size(size_t size_incre) {
    cur_size += size_incre;
    return;
}

void Hash_Map::set_node_t(int node_count) {
    node_size = 2*sizeof(HistoryNode*) + sizeof(HistoryNode) + sizeof(string) + sizeof(int) + (size_t)node_count + (size_t)max(node_count/8,32);
    return;
}

/*
uint32_t Hash_Map::hash_func(pair<string,int> item) {
    uint32_t hash, i;
    string key = item.first;
    for(hash = i = 0; i < key.size(); ++i)
    {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}
*/

void Hash_Map::average_size() {
    long unsigned max = 0;
    long unsigned sum = 0;
    long unsigned item = 0;
    for (unsigned i = 0; i < size; i++) {
        if (History_table[i].size() > 0) item++;
        if (History_table[i].size() > max) max = History_table[i].size();
        sum +=  History_table[i].size();
    }
    cout << "Total non empty entry is " << item << endl;
    cout << "Total history table size is " << sum << endl;
    cout << "Maximum bracket size is " << max << endl;
    cout << "Average bracket size is " << float(sum)/float(item) << endl;
}
    
void Hash_Map::insert(pair<vector<bool>,int>& item,HistoryNode* node) {
    size_t val = hash<vector<bool>>{}(item.first);
    int key = (val + item.second) % size;

    hash_lock[key].lock();
    History_table[key].push_back(make_pair(item,node));
    cur_size += node_size;
    hash_lock[key].unlock();
    return;
}

void Hash_Map::lock_table(int& key) {
    hash_lock[key].lock();
}

void Hash_Map::unlock_table(int& key) {
    hash_lock[key].unlock();
}

/*
HistoryNode* Hash_Map::retrieve(pair<string,int>* item, int* bucket_num) {
    size_t val = hash<string>{}(item->first);
    int key = (val + item->second) % size;
    *bucket_num = key;

    for (auto iter = begin(History_table[key]); iter != History_table[key].end(); iter++) {
        string permutation_src = item->first;
        string permutation_dest = iter->first.first;
        int ending_src = item->second;
        int ending_dest = iter->first.second;
        if (permutation_src == permutation_dest && ending_src == ending_dest) {
            return iter->second;
        }
    }
    
    return NULL;
}
*/

HistoryNode* Hash_Map::retrieve(pair<vector<bool>,int>* item,int key) {
    for (auto iter = History_table[key].begin(); iter != History_table[key].end(); iter++) {
        if (item->first == iter->first.first && item->second == iter->first.second) {
            return iter->second;
        }
    }
    return NULL;
}
