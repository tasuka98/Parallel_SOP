#include <iostream>
#include <stdio.h>
#include <chrono>
#include <string>
#include <stdio.h>

using namespace std;

int main() {
    freopen ("./result/benchmark_result.txt","w",stdout);

    cout << "------------------------TSPLIB------------------------------" << endl;
    string TSPLIB[] = {"ft53.4.sop","ft70.4.sop","rbg050c.sop","rbg109a.sop","rbg150a.sop","rbg174a.sop","prob.42.sop","ry48p.4.sop"};

    for (int i = 0; i < 8; i++) {
            cout << TSPLIB[i] << " ";
            string command_first = "./sop_solver ./tsplib/" + TSPLIB[i];
            cout << endl;
            system(command_first.c_str());
            cout << endl;
    }

    cout << "------------------------SOPLIB---------------------------" << endl;

    string SOPLIB[] = {"R.200.100.1.sop","R.200.100.15.sop","R.200.1000.1.sop","R.200.1000.15.sop","R.300.100.1.sop","R.300.100.15.sop","R.300.1000.1.sop",
                        "R.300.1000.15.sop","R.300.1000.30.sop","R.400.100.1.sop","R.400.100.30.sop","R.400.1000.30.sop","R.500.100.1.sop","R.500.100.30.sop",
                        "R.500.1000.30.sop","R.600.100.1.sop","R.600.100.30.sop","R.600.1000.30.sop","R.700.100.1.sop","R.700.100.30.sop","R.700.1000.30.sop"};

    for (int i = 0; i < 21; i++) {
        cout << SOPLIB[i] << " ";
        string command_first = "./sop_solver ./soplib/" + SOPLIB[i];
        cout << endl;
        system(command_first.c_str());
        cout << endl;
    }

    fclose (stdout);
    return 0;
}