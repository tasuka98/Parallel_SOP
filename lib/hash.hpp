#include <vector>
#include <unordered_map>
#include <list>
#include "history.hpp"

class Hash_Map {
    private:
        int size = 0;
        vector<list<pair<pair<string,int>,HistoryNode>>> History_table;
    public:
        Hash_Map(int size);
        bool isEmpty();
        uint32_t hash_func(pair<string,int> item);
        void insert(pair<string,int>& item,HistoryNode node);
        HistoryNode retrieve(pair<string,int> item);
};



Hash_Map::Hash_Map(int size) {
    this->size = size;
    History_table.resize(size);
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


void Hash_Map::insert(pair<string,int>& item,HistoryNode node) {
    size_t val = hash<string>{}(item.first);
    int key = (val + item.second) % size;
    bool exist = false;
    if (!History_table[key].size()) {
        History_table[key].push_back(make_pair(item,node));
        return;
    }
    for (auto iter = History_table[key].begin(); iter != History_table[key].end(); iter++) {
        string permutation_src = item.first;
        string permutation_dest = iter->first.first;
        int ending_src = item.second;
        int ending_dest = iter->first.second;
        if (permutation_src == permutation_dest && ending_src == ending_dest) {
            exist = true;
            iter->second = node;
            return;
        }
    }

    if (!exist) {
        History_table[key].push_back(make_pair(item,node));
        return;
    }
    return;
}

HistoryNode Hash_Map::retrieve(pair<string,int> item) {
    size_t val = hash<string>{}(item.first);
    int key = (val + item.second) % size;

    for (auto iter = begin(History_table[key]); iter != History_table[key].end(); iter++) {
        string permutation_src = item.first;
        string permutation_dest = iter->first.first;
        int ending_src = item.second;
        int ending_dest = iter->first.second;
        if (permutation_src == permutation_dest && ending_src == ending_dest) {
            return iter->second;
        }
    }
    
    return HistoryNode();
}