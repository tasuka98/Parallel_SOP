#include <iostream>
#include "../lib/solver.hpp"
#include "../lib/hungarian.hpp"
#include <stdio.h>
#include <chrono>
#include <string>
#include <algorithm>
#include <fstream>

using namespace std;

int main(int argc, char*argv[]) {
    ifstream infile("config.txt");
    solver s;

    string line;
    vector<string> setting;
    bool execute = false;
    while (getline(infile, line))
    {
        string command;
        for (auto word : line) {
            if (execute && !isspace(word)) command += word;
            if (word == '=') execute = true;
        }
        if (execute) setting.push_back(command);
        execute = false;
        line.clear();
    }
    infile.close();

    if (argc < 2 || argc > 3) {
        cout << "Incorrect format.....\n";
        cout << "<Instant Name>\n";
        exit(-1);
    }
    else if (argc == 2) {
        s.assign_parameter(setting);
        s.solve(argv[1],-1);
    }
    else if (argc == 3) {
        s.assign_parameter(setting);
        s.solve(argv[1],atoi(argv[2]));
    }
    
    return 0;
}