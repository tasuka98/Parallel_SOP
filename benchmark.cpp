#include <iostream>
#include <stdio.h>
#include <chrono>
#include <string>
#include <stdio.h>

using namespace std;

int main() {

    freopen ("./result/benchmark_result.txt","w",stdout);

    string TSPLIB[] = {"br17.10.sop","br17.12.sop","ESC07.sop",
                        "ESC11.sop","ESC12.sop","ESC25.sop","ESC47.sop","ESC63.sop",
                        "ft53.4.sop","ft53.4.sop","ft70.4.sop","rbg050c.sop"
                        ,"rbg109a.sop","rbg150a.sop","p43.4.sop","prob.42.sop","ry48p.4.sop"};
    
    string Compiler[] = {"gsm.153.124.sop","gsm.444.350.sop","gsm.462.77.sop","jpeg.1483.25.sop",
                        "jpeg.3184.107.sop","jpeg.3195.85.sop","jpeg.3198.93.sop","jpeg.3203.135.sop",
                        "jpeg.3740.15.sop","jpeg.4154.36.sop","jpeg.4753.54.sop","susan.248.197.sop",
                        "susan.260.158.sop","susan.343.182.sop","typeset.10835.26.sop","typeset.12395.43.sop"
                        ,"typeset.15087.23.sop","typeset.15577.36.sop","typeset.16000.68.sop","typeset.1723.25.sop"
                        ,"typeset.19972.246.sop","typeset.4391.240.sop","typeset.4597.45.sop","typeset.4724.433.sop"
                        ,"typeset.5797.33.sop","typeset.5881.246.sop"};

    for (int i = 0; i < 26; i++) {
        cout << Compiler[i] << " ";
        string command_first = "./sop_solver ./mibench/" + Compiler[i]+" NN 100";
        string command_second = "./sop_solver ./mibench/" + Compiler[i]+" DH 100";
        string command_third = "./main ./mibench/" + Compiler[i] + " 200 1237935";
        cout << endl;
        system(command_first.c_str());
        system(command_second.c_str());
        system(command_third.c_str());
        cout << endl;
    }
    
    fclose (stdout);

    return 0;
}