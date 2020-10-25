#include <iostream>
#include <stdio.h>
#include <chrono>
#include <string>
#include <stdio.h>
#include <unistd.h>

using namespace std;

int main() {
    freopen ("./result/benchmark_result.txt","w",stdout);

    cout << "------------------------TSPLIB------------------------------" << endl;
    string TSPLIB[] = {"ft53.4.sop","ft70.4.sop","rbg050c.sop","rbg109a.sop","rbg150a.sop","rbg174a.sop","ry48p.4.sop"};

    for (int i = 0; i < 7; i++) {
            cout << TSPLIB[i] << " ";
            string command_first = "./sop_solver ./tsplib/" + TSPLIB[i];
            string command_second = "./sop_solver_full ./tsplib/" + TSPLIB[i];
            cout << endl;
            system(command_first.c_str());
            system(command_second.c_str());
            cout << endl;
            usleep(2000000);
    }

    cout << "------------------------SOPLIB---------------------------" << endl;

    string SOPLIB[] = {"R.200.100.1.sop","R.200.100.15.sop","R.200.1000.1.sop","R.200.1000.15.sop","R.300.100.1.sop",
                        "R.300.100.15.sop","R.300.1000.1.sop","R.300.1000.30.sop","R.400.100.1.sop","R.400.100.30.sop",
                        "R.400.1000.30.sop","R.500.100.1.sop","R.500.100.30.sop","R.500.1000.30.sop","R.600.100.30.sop",
                        "R.600.1000.30.sop","R.700.100.1.sop","R.700.100.30.sop","R.700.1000.30.sop"};

    for (int i = 0; i < 19; i++) {
        cout << SOPLIB[i] << " ";
        string command_first = "./sop_solver ./soplib/" + SOPLIB[i];
        string command_second = "./sop_solver_full ./soplib/" + SOPLIB[i];
        cout << endl;
        system(command_first.c_str());
        system(command_second.c_str());
        cout << endl;
        usleep(2000000);
    }

    fclose (stdout);
    return 0;
}