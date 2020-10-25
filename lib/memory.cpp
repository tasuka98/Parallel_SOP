#include <vector>
#include "memory.hpp"
#include <mutex>

mutex manager_mutex;

memory_manager::memory_manager(int size) {
    block_size = size;
    counter = 0;
    HistoryNode* new_allocation = new HistoryNode[block_size];
    history_block.push_back(new_allocation);
    return;
}

HistoryNode* memory_manager::retrieve_history_table() {
    manager_mutex.lock();
    HistoryNode* his_node;
    if (counter >= block_size) {
        counter = 0;
        HistoryNode* new_allocation = new HistoryNode[block_size];
        history_block.push_back(new_allocation);
    }
    his_node = (history_block.back())+counter;
    counter++;
    manager_mutex.unlock();
    return his_node;
}

