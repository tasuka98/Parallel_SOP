#include <iostream>
#include <stdio.h>
#include <chrono>
#include <string>
#include <stdio.h>

using namespace std;

int main() {
    freopen ("./result/benchmark_result.txt","w",stdout);

    cout << "------------------------TSPLIB------------------------------" << endl;
    string TSPLIB[] = {"rbg174a.sop"};

    for (int i = 0; i < 1; i++) {
            cout << TSPLIB[i] << " ";
            string command_first = "./sop_solver ./tsplib/" + TSPLIB[i]+" DH 7200";
            cout << endl;
            system(command_first.c_str());
            cout << endl;
    }

    cout << "------------------------SOPLIB---------------------------" << endl;

    string SOPLIB[] = {"R.200.100.1.sop","R.200.100.15.sop","R.200.1000.1.sop","R.300.100.1.sop","R.300.100.15.sop","R.300.1000.1.sop",
                        "R.300.1000.15.sop","R.400.100.1.sop","R.500.100.1.sop","R.600.100.1.sop","R.700.100.1.sop"};

    for (int i = 0; i < 11; i++) {
        cout << SOPLIB[i] << " ";
        string command_first = "./sop_solver ./soplib/" + SOPLIB[i]+" DH 7200";
        cout << endl;
        system(command_first.c_str());
        cout << endl;
    }

    fclose (stdout);
    return 0;
}