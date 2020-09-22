#include <iostream>
#include <stdio.h>
#include <chrono>
#include <string>
#include <stdio.h>

using namespace std;

int main() {
    freopen ("./result/benchmark_result.txt","w",stdout);

    cout << "------------------------TSPLIB------------------------------" << endl;
    string TSPLIB[] = {"rbg050c.sop","prob.42.sop"};

    for (int i = 0; i < 2; i++) {
        cout << TSPLIB[i] << " ";
        string command_first = "./sop_solver ./tsplib/" + TSPLIB[i]+" NN 100";
        string command_second = "./sop_solver ./tsplib/" + TSPLIB[i]+" DH 100";
        string command_third = "./main_searchorder ./tsplib/" + TSPLIB[i] + " 7200 1237935 1";
        string command_fourth = "./main_sequential ./tsplib/" + TSPLIB[i] + " 7200 1237935";
        cout << endl;
        system(command_first.c_str());
        system(command_second.c_str());
        system(command_third.c_str());
        system(command_fourth.c_str());
        cout << endl;
    }

    cout << "------------------------SOPLIB---------------------------" << endl;

    string SOPLIB[] = {"R.200.1000.15.sop","R.700.100.30.sop","R.700.1000.30.sop"};

    for (int i = 0; i < 3; i++) {
        cout << SOPLIB[i] << " ";
        string command_first = "./sop_solver ./soplib/" + SOPLIB[i]+" NN 100";
        string command_second = "./sop_solver ./soplib/" + SOPLIB[i]+" DH 100";
        string command_third = "./main_searchorder ./soplib/" + SOPLIB[i] + " 7200 1237935 1";
        string command_fourth = "./main_sequential ./soplib/" + SOPLIB[i] + " 7200 1237935";
        cout << endl;
        system(command_first.c_str());
        system(command_second.c_str());
        system(command_third.c_str());
        system(command_fourth.c_str());
        cout << endl;
    }

    fclose (stdout);
    return 0;
}