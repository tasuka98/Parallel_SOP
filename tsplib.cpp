#include <iostream>
#include <stdio.h>
#include <chrono>
#include <string>
#include <stdio.h>

using namespace std;

int main() {

    freopen ("./result/tsplib_result.txt","w",stdout);
    
    cout << "------------------------TSPLIB-----------------------------" << endl;
    
    string TSPLIB[] = {"br17.10.sop","br17.12.sop","ESC07_last.sop",
                        "ESC11_last.sop","ESC12.sop","ESC25.sop","ESC47.sop","ESC63.sop",
                        "ft53.4.sop","rbg109a.sop","rbg150a.sop","p43.4.sop","ry48p.4.sop"};

    for (int i = 0; i < 13; i++) {
        cout << TSPLIB[i] << " ";
        string command_first = "./sop_solver ./tsplib/" + TSPLIB[i]+" NN 100";
        string command_second = "./sop_solver ./tsplib/" + TSPLIB[i]+" DH 100";
        string command_third = "./main ./tsplib/" + TSPLIB[i] + " 7200 1237935";
        cout << endl;
        system(command_first.c_str());
        system(command_second.c_str());
        system(command_third.c_str());
        cout << endl;
    }

    fclose (stdout);
    return 0;
}