#include <iostream>
#include "../lib/solver.hpp"
#include "../lib/hungarian.h"
#include <stdio.h>

using namespace std;

int main(int argc, char*argv[]) {
    solver s;
    s.solve(argv[1]);
    
    return 0;
}