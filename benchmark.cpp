#include <iostream>
#include <stdio.h>
#include <chrono>
#include <string>
#include <stdio.h>

using namespace std;

int main() {
    freopen ("./result/benchmark_result.txt","w",stdout);

    cout << "------------------------Compiler-----------------------------" << endl;

    string Compiler[] = {"gsm.153.124.sop","gsm.444.350.sop","gsm.462.77.sop","jpeg.1483.25.sop",
                        "jpeg.3184.107.sop","jpeg.3195.85.sop","jpeg.3198.93.sop","jpeg.3203.135.sop",
                        "jpeg.3740.15.sop","jpeg.4154.36.sop","jpeg.4753.54.sop","susan.248.197.sop",
                        "susan.260.158.sop","susan.343.182.sop","typeset.10835.26.sop","typeset.12395.43.sop",
                        "typeset.15087.23.sop","typeset.15577.36.sop","typeset.16000.68.sop","typeset.1723.25.sop",
                        "typeset.19972.246.sop","typeset.4391.240.sop","typeset.4597.45.sop","typeset.4724.433.sop",
                        "typeset.5797.33.sop","typeset.5881.246.sop"};

    for (int i = 0; i < 26; i++) {
        cout << Compiler[i] << " ";
        string command_first = "./sop_solver ./mibench/" + Compiler[i]+" NN 100";
        string command_second = "./sop_solver ./mibench/" + Compiler[i]+" DH 100";
        string command_third = "./main_searchorder ./mibench/" + Compiler[i] + " 7200 1237935 1";
        string command_fourth = "./main_sequential ./mibench/" + Compiler[i] + " 7200 1237935";
        cout << endl;
        system(command_first.c_str());
        system(command_second.c_str());
        system(command_third.c_str());
        system(command_fourth.c_str());
        cout << endl;
    }

    cout << "------------------------TSPLIB------------------------------" << endl;
    
    string TSPLIB[] = {"br17.10.sop","br17.12.sop","ESC07_last.sop",
                        "ESC11_last.sop","ESC12.sop","ESC25.sop","ESC47.sop","ESC63.sop",
                        "ft53.4.sop","rbg109a.sop","rbg150a.sop","p43.4.sop","ry48p.4.sop"};

    for (int i = 0; i < 13; i++) {
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

    string SOPLIB[] = {"R.200.100.30.sop","R.200.100.60.sop","R.200.1000.30.sop","R.200.1000.60.sop",
                        "R.300.100.30.sop","R.300.100.60.sop","R.300.1000.30.sop","R.300.1000.60.sop",
                        "R.400.100.30.sop","R.400.100.60.sop","R.400.1000.30.sop","R.400.1000.60.sop",
                        "R.500.100.30.sop","R.500.100.60.sop","R.500.1000.30.sop","R.500.1000.60.sop",
                        "R.600.100.30.sop","R.600.100.60.sop","R.600.1000.30.sop","R.600.1000.60.sop",
                        "R.700.100.60.sop","R.700.1000.60.sop"};

    for (int i = 0; i < 22; i++) {
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