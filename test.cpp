#include <iostream>
#include <vector>
#include <string.h>
#include <map>

extern "C" {
    #include "glpk.h"
}

using namespace std;


int main(int argc, char const *argv[]) {
    map<int, bool> M;
    M[1] = true;
    cout << M[1] << endl;
    cout << M[2] << endl;
    return 0;
}
