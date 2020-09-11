#include <iostream>
#include "../lib/solver.hpp"
#include "../lib/hungarian.h"
#include <stdio.h>
#include <chrono>

using namespace std;

int main(int argc, char*argv[]) {
    solver s;
    
    if (argc < 4) {
        cout << "Input Argument Incorrect";
        exit(-1);
    }
    s.solve(argv[1],argv[2],atoi(argv[3]));
    
    
    return 0;
}