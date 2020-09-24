#include <iostream>
#include "../lib/solver.hpp"
#include "../lib/hungarian.h"
#include <stdio.h>
#include <chrono>

using namespace std;

int main(int argc, char*argv[]) {
    solver s;
    if (argc < 6) {
        cout << "Incorrect format.....\n";
        cout << "<Instant Name> <Static Enumeration> <Time Limit> <Pool size> <Thread Size>\n";
        exit(-1);
    }
    s.solve(argv[1],argv[2],atoi(argv[3]),atoi(argv[4]),atoi(argv[5]));
    return 0;
}