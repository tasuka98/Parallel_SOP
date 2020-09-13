#include <iostream>
#include "../lib/solver.hpp"
#include "../lib/hungarian.h"
#include <stdio.h>
#include <chrono>

using namespace std;

int main(int argc, char*argv[]) {
    solver s;
    
    string TSPLIB[] = {"br17.10.sop","br17.12.sop","ESC07.sop",
                        "ESC11.sop","ESC12.sop","ESC25.sop","ESC47.sop","ESC63.sop",
                        "ft53.4.sop","ft53.4.sop","ft70.4.sop","rbg050c.sop"
                        ,"rbg109a.sop","rbg150a.sop","p43.4.sop","prob.42.sop","ry48p.4.sop"};
    
    string Compiler[] = {"gsm.153.124.sop","gsm.444.350.sop","gsm.462.77.sop","jpeg.1483.25.sop"
                        "jpeg.3184.107.sop","jpeg.3195.85.sop","jpeg.3198.93.sop","jpeg.3203.135.sop",
                        "jpeg.3740.15.sop","jpeg.4154.36.sop","jpeg.4753.54.sop","susan.248.197.sop",
                        "susan.260.158.sop","susan.343.182.sop","typeset.10835.26.sop","typeset.12395.43.sop"
                        ,"typeset.15087.23.sop","typeset.15577.36.sop","typeset.16000.68.sop","typeset.1723.25.sop"
                        ,"typeset.19972.246.sop","typeset.4391.240.sop","typeset.4597.45.sop","typeset.4724.433.sop"
                        ,"typeset.5797.33.sop","typeset.5881.246.sop"};


    for (int i = 0; i < sizeof(Compiler)/sizeof(Compiler[0]); i++) {
        cout << Compiler[i] << " ";
        s.solve("./mibench/"+Compiler[i],"RD",7200);
    }
    
    
    
    return 0;
}