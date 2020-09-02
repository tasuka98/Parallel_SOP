#include <iostream>
#include "../lib/solver.hpp"
#include "../lib/hungarian.h"
#include <stdio.h>
#include <chrono>

using namespace std;

int main(int argc, char*argv[]) {
    solver s;
    auto start_time = chrono::high_resolution_clock::now();
    s.solve(argv[1]);
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::microseconds>( end_time - start_time ).count();
    cout << "B&B solver run time is " << duration << " microseconds" << endl;
    return 0;
}