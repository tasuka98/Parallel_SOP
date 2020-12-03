#include <iostream>
#include "load_test.hpp"

using namespace std;

void load_test::Display_info() {
    cout << ".....................Stolen Workload Info...................." << endl;
    cout << "Workload Level is " << Stolen_Level << endl;
    cout << "Initial LB is " << Stolen_Lb << endl;
    cout << "Total number of nodes stolen are " << Stolen_Load << endl;
    cout << "Number of Updated Solution during the interval is " << Num_of_Sol_Updates << endl;
    cout << ".....................Stolen Workload Info...................." << endl;
    cout << endl;
    return;
}

void load_test::Clear_info() {
    Stolen_Load = 0;
    Stolen_Level = 0;
    Num_of_Sol_Updates = 0;
    Stolen_Lb = 0;
    return;
}