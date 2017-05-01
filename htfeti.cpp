#include <iostream>
#include <string> 
#include "include/Options.hpp"
#include "include/Cluster.hpp"
//#include "Driver/DissectionSolver.hpp"


using namespace std;
int main(int argc, char *argv[]){


    cout << argv[0] << endl;
    Options options;


    int n_subdomOnCluster;
    string path2data;


    int var = 4;

    switch (var)
    {
        case (1):
            // small variant: 8 subdom, each 3x3x3 elements
            n_subdomOnCluster = 8;
            path2data = "../data/";
            break;
        case (2):
            // larger variant: 8 subdom, each 3x3x3 elements
            n_subdomOnCluster = 8;
            path2data = "../data2/";
            break;
        case (3):
            // 2nd larger variant: 8 subdom, each 3x3x3 elements
            n_subdomOnCluster = 27;
            path2data = "../data3/";
            break;
        case (4):
            // 2nd larger variant: 8 subdom, each 3x3x3 elements
            // 'kernel' gluing
            n_subdomOnCluster = 27;
            path2data = "../data4/";
            break;
    }


    options.set_values(n_subdomOnCluster, path2data);
//
    Cluster cluster(options);
    cout << "----------------- done -----------------\n" ;
    return 0;
}
