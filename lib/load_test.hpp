#ifndef LOADTEST_H
#define LOADTEST_H

#include <iostream>

class load_test {
    public:
        int Stolen_Load = 0;
        int Stolen_Level = 0;
        int Num_of_Sol_Updates = 0;
        int Stolen_Lb = 0;
        void Display_info();
        void Clear_info();
};

#endif