#include <iostream>
#include <string> 
#include "Options.cpp"
#include "Cluster.h"
//#include "Driver/DissectionSolver.hpp"

int main(){ 
    Options options;
    int n_subdomOnCluster = 2;
    std::string path2data = "data";
    options.set_values(n_subdomOnCluster, path2data);
    Cluster cluster(options);
    std::cout << "----------------- done -----------------\n" ;
    return 0;
}
