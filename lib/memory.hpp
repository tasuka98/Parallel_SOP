#include <vector>
#include "history.hpp"

class memory_manager {
    private:
        vector<HistoryNode*> history_block;
        int counter = 0;
        int block_size = 0;
    public:
        memory_manager(int size);
        HistoryNode* retrieve_history_table();
};