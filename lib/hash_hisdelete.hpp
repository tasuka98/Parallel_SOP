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
#define HIS_BLK_SIZE 81920
#define FREED_SIZE 100000000
#define COVER_AREA 25
#define MEMORY_RESTRIC 0.85
#define DEPLETION_TIME_FRAME 10
typedef list<pair<pair<vector<bool>,int>,HistoryNode*>> Bucket;
typedef pair<Bucket*,vector<int>*> Memory_Record_Blk;

static vector<mutex> hash_lock;
atomic<bool> memory_limit_reached(false);

class Mem_Allocator {
    private:
        std::chrono::time_point<std::chrono::system_clock> start_timer;
        vector<vector<Memory_Record_Blk>> Bucket_blk;
        vector<Bucket*> Freed_blk;
        vector<HistoryNode*> Freed_node;
        vector<unsigned> counter;
        HistoryNode* history_block = NULL;
        unsigned node_counter = 0;
        bool depleted = false;
    public:
        Mem_Allocator(int node_count);
        Bucket* Get_bucket(int level,int key);
        HistoryNode* retrieve_his_node();
        HistoryNode* Retrieve_FreeNode();
        void Insert_FreeNode(HistoryNode* node);
        size_t Abandon_level(size_t node_size,vector<Bucket*>* History_table);
        bool Check_Depletion();
};

Mem_Allocator::Mem_Allocator(int node_count) {
    Bucket_blk.resize(node_count);
    for (int i = 0; i < node_count; i++) {
        Bucket* new_alloc = new Bucket[BUCKET_BLK_SIZE];
        vector<int>* Depth_record = new vector<int>();
        Bucket_blk[i].push_back(make_pair(new_alloc,Depth_record));
    }
    counter = vector<unsigned>(node_count,0);
}

HistoryNode* Mem_Allocator::retrieve_his_node() {
    if (node_counter >= HIS_BLK_SIZE || history_block == NULL) {
        history_block = new HistoryNode[HIS_BLK_SIZE];
        node_counter = 0;
    }
    HistoryNode* node = history_block + node_counter;
    node_counter++;
    return node;
}

size_t Mem_Allocator::Abandon_level(size_t node_size,vector<Bucket*>* History_table) {
    size_t freed_mem = 0;
    unsigned dest = 0;
    bool found = false;
    while (freed_mem < FREED_SIZE) {
        for (unsigned i = (Bucket_blk.size() * 0.5); i < Bucket_blk.size(); i++) {
            if (Bucket_blk[i].size() > 1 && i > dest) {
                found = true;
                dest = i; 
            }
        }
        
        if (!found) {
            start_timer = std::chrono::system_clock::now();
            depleted = true;
            break;
        }

        for (unsigned i = 0; i < Bucket_blk[dest].size(); i++) {
            for (auto key : *(Bucket_blk[dest][i].second)) {
                freed_mem += ((*History_table)[key]->size() * node_size + sizeof(Bucket));
                hash_lock[key/COVER_AREA].lock();
                (*History_table)[key] = NULL;
                hash_lock[key/COVER_AREA].unlock();
            }
            delete(Bucket_blk[dest][i].second);
            Freed_blk.push_back(Bucket_blk[dest][i].first);
        }
        Bucket_blk[dest].clear();
        counter[dest] = BUCKET_BLK_SIZE;
        dest = 0;
        found = false;
    }

    return freed_mem;
}

void Mem_Allocator::Insert_FreeNode(HistoryNode* node) {
    Freed_node.push_back(node);
}

HistoryNode* Mem_Allocator::Retrieve_FreeNode() {
    if (Freed_node.empty()) return NULL;

    HistoryNode* fn = Freed_node.back();
    Freed_node.pop_back();
    return fn;
}

bool Mem_Allocator::Check_Depletion() {
    if (depleted) {
        auto end_timer = std::chrono::system_clock::now();
        if (std::chrono::duration<double>(end_timer - start_timer).count() > DEPLETION_TIME_FRAME) {
            depleted = false;
        }
    }
    return depleted;
}


Bucket* Mem_Allocator::Get_bucket(int level,int key) {
    if (counter[level] == BUCKET_BLK_SIZE) {
        if (!Freed_blk.empty()) {
            vector<int>* Depth_record = new vector<int>();
            Bucket_blk[level].push_back(make_pair(Freed_blk.back(),Depth_record));
            counter[level] = 0;
            Freed_blk.pop_back();
        }
        else if (memory_limit_reached == false) {
            Bucket* new_alloc = new Bucket[BUCKET_BLK_SIZE];
            vector<int>* Depth_record = new vector<int>();
            Bucket_blk[level].push_back(make_pair(new_alloc,Depth_record));
            counter[level] = 0;
        }
        else {
            return NULL;
        }
    }
    Bucket* bucket = Bucket_blk[level].back().first + counter[level];
    Bucket_blk[level].back().second->push_back(key);
    counter[level]++;
    return bucket;
}

class Hash_Map {
    private:
        atomic<unsigned long> cur_size;
        unsigned long max_size = 0;
        unsigned size = 0;
        size_t node_size = 0;
        vector<Bucket*> History_table;
        vector<vector<vector<int>>> Depth_info;
        vector<Mem_Allocator> Mem_Manager;
    public:
        atomic<long> num_of_waits;
        Hash_Map(size_t size);
        bool isEmpty();
        size_t get_max_size();
        size_t get_cur_size();
        uint32_t hash_func(pair<vector<bool>,int> item);
        bool insert(pair<vector<bool>,int>& item,int prefix_cost,int lower_bound,unsigned thread_id,unsigned level);
        void calculate_SD();
        void adjust_max_size(int node_count, int GPQ_size);
        void lock_table(int& key);
        void unlock_table(int& key);
        void increase_size(size_t size_incre);
        void set_node_t(int node_size);
        void free_mem(int thread_id);
        void set_up_mem(int thread_num,int node_count);
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
    max_size = ((double)info.freeram * MEMORY_RESTRIC) - (size * 4) - ((size/COVER_AREA + 1) * sizeof(mutex));
    cur_size = 0;
    hash_lock = vector<mutex>(size/COVER_AREA + 1);
    History_table.resize(size);

    cout << "Max bucket size is " << max_size / 1000000 << " MB" << endl;
    
    for (unsigned i = 0; i < size; i++) History_table[i] = NULL;
}

void Hash_Map::free_mem(int thread_id) {
    if (!memory_limit_reached) {
        memory_limit_reached = true;
    }
    if (!Mem_Manager[thread_id].Check_Depletion()) {
        cur_size -= Mem_Manager[thread_id].Abandon_level(node_size,&History_table);
    }
}

void Hash_Map::adjust_max_size(int node_count,int GPQ_size) {
    max_size -= (GPQ_size * (12*node_count + 2*node_count*node_count) * 4);
}

void Hash_Map::set_up_mem(int thread_num,int node_count) {
    for (int i = 0; i < thread_num; i++) {
        Mem_Manager.push_back(Mem_Allocator(node_count));
    }
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
    
bool Hash_Map::insert(pair<vector<bool>,int>& item,int prefix_cost,int lower_bound,unsigned thread_id,unsigned level) {
    size_t val = hash<vector<bool>>{}(item.first);
    int key = (val + item.second) % size;

    hash_lock[key/COVER_AREA].lock();

    if (History_table[key] == NULL) {
        History_table[key] = Mem_Manager[thread_id].Get_bucket(level,key);
        if (History_table[key] == NULL) {
            hash_lock[key/COVER_AREA].unlock();
            return false;
        }

        if (!memory_limit_reached) {
            HistoryNode* node = Mem_Manager[thread_id].retrieve_his_node();
            node->prefix_cost = prefix_cost;
            node->lower_bound = lower_bound;
            History_table[key]->push_back(make_pair(item,node));
        }
        else {
            if (!History_table[key]->empty()) {
                while (History_table[key]->size() > 1) {
                    Mem_Manager[thread_id].Insert_FreeNode(History_table[key]->back().second);
                    History_table[key]->pop_back();
                }
                History_table[key]->back().first = item;
                History_table[key]->back().second->prefix_cost = prefix_cost;
                History_table[key]->back().second->lower_bound = lower_bound;
            }
            else {
                HistoryNode* fn = Mem_Manager[thread_id].Retrieve_FreeNode();
                if (fn == NULL) {
                    hash_lock[key/COVER_AREA].unlock();
                    return false;
                }
                else {
                    fn->lower_bound = lower_bound;
                    fn->prefix_cost = prefix_cost;
                    History_table[key]->push_back(make_pair(item,fn));
                }
            }
        }

        cur_size += (node_size + sizeof(int) + sizeof(Bucket));
        hash_lock[key/COVER_AREA].unlock();
        return true;
    }

    for (auto iter = History_table[key]->begin(); iter != History_table[key]->end(); iter++) {
        if (item.first == iter->first.first && item.second == iter->first.second) {
            hash_lock[key/COVER_AREA].unlock();
            return false;
        }
    }

    HistoryNode* node;
    if (!memory_limit_reached) {
        node = Mem_Manager[thread_id].retrieve_his_node();
        node->prefix_cost = prefix_cost;
        node->lower_bound = lower_bound;
        History_table[key]->push_back(make_pair(item,node));
    }
    else {
        HistoryNode* fn = Mem_Manager[thread_id].Retrieve_FreeNode();
        if (fn == NULL) {
            hash_lock[key/COVER_AREA].unlock();
            return false;
        }
        else {
            fn->prefix_cost = prefix_cost;
            fn->lower_bound = lower_bound;
            History_table[key]->push_back(make_pair(item,fn));
        }
    }

    cur_size += (node_size + sizeof(int));
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