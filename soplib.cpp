#include <iostream>
#include <stdio.h>
#include <chrono>
#include <string>
#include <stdio.h>

using namespace std;

int main() {

    
    freopen ("./result/soplib_result.txt","w",stdout);

    cout << "------------------------SOPLIB---------------------------" << endl;

    //R.300.1000.30.sop
    //R.400.100.30.sop
    string SOPLIB[] = {"R.200.100.30.sop","R.200.100.60.sop","R.200.1000.30.sop","R.200.1000.60.sop",
                        "R.300.100.30.sop","R.300.100.60.sop","R.300.1000.60.sop","R.400.100.60.sop",
                        "R.400.1000.60.sop","R.500.100.60.sop","R.500.1000.60.sop","R.600.100.30.sop",
                        "R.600.100.60.sop","R.600.1000.60.sop","R.700.100.60.sop","R.700.1000.60.sop"};

    for (int i = 0; i < 16; i++) {
        cout << SOPLIB[i] << " ";
        string command_first = "./sop_solver ./soplib/" + SOPLIB[i]+" NN 100";
        string command_second = "./sop_solver ./soplib/" + SOPLIB[i]+" DH 100";
        string command_third = "./main ./soplib/" + SOPLIB[i] + " 7200 1237935";
        cout << endl;
        system(command_first.c_str());
        system(command_second.c_str());
        system(command_third.c_str());
        cout << endl;
    }

    fclose (stdout);
    return 0;
}