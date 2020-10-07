#include <iostream>
#include "../lib/solver.hpp"
#include "../lib/hungarian.h"
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

    if (argc < 2) {
        cout << "Incorrect format.....\n";
        cout << "<Instant Name>\n";
        exit(-1);
    }
    s.solve(argv[1],setting[0],atoi(setting[1].c_str()),atoi(setting[2].c_str()),atoi(setting[3].c_str()),atoi(setting[4].c_str()),setting[5].c_str(),setting[6].c_str());
    return 0;
}