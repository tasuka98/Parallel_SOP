#include <vector>
#include <unordered_map>
#include <iostream>
#include <list>
#include <mutex>
#include <cmath>
#include <cerrno>
#include <atomic>
#include <sys/sysinfo.h>
#include "history.hpp"

#define BUCKET_BLK_SIZE 81920
#define COVER_AREA 25
typedef list<pair<pair<vector<bool>,int>,HistoryNode*>> Bucket;

static vector<mutex> hash_lock;

class Mem_Allocator {
    private:
        Bucket* Bucket_blk;
        unsigned counter;
    public:
        Mem_Allocator();
        Bucket* Get_bucket();
};

Mem_Allocator::Mem_Allocator() {
    Bucket_blk = new Bucket[BUCKET_BLK_SIZE];
    counter = 0;
}

Bucket* Mem_Allocator::Get_bucket() {
    if (counter == BUCKET_BLK_SIZE) {
        Bucket_blk = new Bucket[BUCKET_BLK_SIZE];
        counter = 0;
    }
    Bucket* bucket = Bucket_blk + counter;
    counter++;
    return bucket;
}

class Hash_Map {
    private:
        atomic<size_t> cur_size;
        size_t size = 0;
        size_t max_size = 0;
        size_t node_size = 0;
        vector<Bucket*> History_table;
        //vector<vector<HistoryNode*>> Stratified_table;
        vector<Mem_Allocator> Mem_Manager;
    public:
        atomic<long> num_of_waits;
        Hash_Map(size_t size);
        bool isEmpty();
        size_t get_max_size();
        size_t get_cur_size();
        uint32_t hash_func(pair<vector<bool>,int> item);
        bool insert(pair<vector<bool>,int>& item,HistoryNode* node,unsigned thread_id,unsigned level);
        void calculate_SD();
        void lock_table(int& key);
        void unlock_table(int& key);
        void increase_size(size_t size_incre);
        void set_node_t(int node_size);
        void set_up_mem(int thread_num);
        void set_up_table(int ins_size);
        void average_size();
        HistoryNode* retrieve(pair<vector<bool>,int>* item,int key);
};

Hash_Map::Hash_Map(size_t size) {
    this->size = size;
    struct sysinfo info;
	if(sysinfo(&info) != 0){
        cout << "can't retrieve sys mem info\n";
		exit(1);
	}
    max_size = ((double)info.freeram * 0.85) - (size * 4) - ((size/COVER_AREA + 1) * sizeof(mutex));
    cur_size = 0;
    hash_lock = vector<mutex>(size/COVER_AREA + 1);
    History_table.resize(size);

    cout << "Max bucket size is " << max_size / 1000000 << " MB" << endl;
    
    for (unsigned i = 0; i < size; i++) History_table[i] = NULL;
}

void Hash_Map::set_up_mem(int thread_num) {
    Mem_Manager = vector<Mem_Allocator>(thread_num);
}

void Hash_Map::set_up_table(int ins_size) {
    //Stratified_table = vector<vector<HistoryNode*>>(node_size);
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
    node_size = 2*sizeof(HistoryNode*) + sizeof(HistoryNode) + sizeof(vector<bool>) + sizeof(int) + (size_t)max(node_count/8,48);
    return;
}

void Hash_Map::calculate_SD() {
    double mean = 0;
    int non_empty_buckets = 0;
    for (auto temp_bucket : History_table) {
        if (temp_bucket != NULL) {
            mean += temp_bucket->size();
            non_empty_buckets++;
        }
    }
    mean = mean / non_empty_buckets;
    double summation = 0;

    for (auto temp_bucket : History_table) {
        if (temp_bucket != NULL) {
            summation += pow((temp_bucket->size() - mean),2);
        }
    }
    summation = summation / non_empty_buckets;
    summation = sqrt(summation);
    cout << "Stnadard deviation of the history table is " << summation << endl;

    return;
}

void Hash_Map::average_size() {
    long unsigned max = 0;
    long unsigned min = INT_MAX;
    long unsigned sum = 0;
    long unsigned item = 0;
    for (unsigned i = 0; i < size; i++) {
        if (History_table[i] != NULL) {
            if (History_table[i]->size() > 0) item++;
            if (History_table[i]->size() < min) min = History_table[i]->size();
            if (History_table[i]->size() > max) max = History_table[i]->size();
            sum +=  History_table[i]->size();
        }
    }
    cout << "Total non empty entry is " << item << endl;
    cout << "Total history table size is " << sum << endl;
    cout << "Maximum bucket size is " << max << endl;
    cout << "Minimum bucket size is " << min << endl;
    cout << "Average bucket size is " << float(sum)/float(item) << endl;
    calculate_SD();
    return;
}
    
bool Hash_Map::insert(pair<vector<bool>,int>& item,HistoryNode* node,unsigned thread_id,unsigned level) {
    size_t val = hash<vector<bool>>{}(item.first);
    int key = (val + item.second) % size;

    hash_lock[key/COVER_AREA].lock();

    if (History_table[key] == NULL) {
        History_table[key] = Mem_Manager[thread_id].Get_bucket();
        History_table[key]->push_back(make_pair(item,node));
        cur_size += node_size;
        cur_size += sizeof(Bucket);
        hash_lock[key/COVER_AREA].unlock();
        return true;
    }

    for (auto iter = History_table[key]->begin(); iter != History_table[key]->end(); iter++) {
        if (item.first == iter->first.first && item.second == iter->first.second) {
            hash_lock[key/COVER_AREA].unlock();
            return false;
        }
    }

    History_table[key]->push_back(make_pair(item,node));
    cur_size += node_size;
    hash_lock[key/COVER_AREA].unlock();
    return true;
}

void Hash_Map::lock_table(int& key) {
    hash_lock[key/COVER_AREA].lock();
}

void Hash_Map::unlock_table(int& key) {
    hash_lock[key/COVER_AREA].unlock();
}

HistoryNode* Hash_Map::retrieve(pair<vector<bool>,int>* item,int key) {
    if (History_table[key] == NULL) return NULL;

    for (auto iter = History_table[key]->begin(); iter != History_table[key]->end(); iter++) {
        if (item->first == iter->first.first && item->second == iter->first.second) {
            return iter->second;
        }
    }
    return NULL;
}