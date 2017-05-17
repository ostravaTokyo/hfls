#include <iostream>
#include <string> 
#include "src/Options.hpp"
#include "src/Cluster.hpp"



using namespace std;
int main(int argc, char *argv[]){


    cout << argv[0] << endl;
    Options options;


    int n_subdomOnCluster;
    string path2data;



    path2data = "../data/";
    options.set_values(n_subdomOnCluster, path2data);


    if (argc > 1) {
       options.meshSetting.N[0] =  atoi(argv[1]);
    }
    if (argc == 2) {
       options.meshSetting.N[1] =  atoi(argv[2]);
    }
    if (argc == 3) {
       options.meshSetting.N[2] =  atoi(argv[3]);
    }
    if (argc > 4) {
       options.meshSetting.n[0] =  atoi(argv[4]);
    }
    if (argc > 5) {
       options.meshSetting.n[1] =  atoi(argv[5]);
    }
    if (argc > 6) {
       options.meshSetting.n[2] =  atoi(argv[6]);
    }




//    if (argc == 1) {
//       var = atoi(argv[1]);
//    }
//
//    switch (var)
//    {
//        case (0):
//            // small variant: 8 subdom, each 3x3x3 elements
//            n_subdomOnCluster = 8;
//            path2data = "../data/";
//            break;
//        case (1):
//            // larger variant: 27 subdom, each 10x10x10 elements
//            n_subdomOnCluster = 27;
//            path2data = "../data1/";
//            break;
//        case (2):
//            // larger variant: 8 subdom, each 6x6x6 elements
//            n_subdomOnCluster = 8;
//            path2data = "../data2/";
//            break;
//        case (3):
//            // 2nd larger variant: 27 subdom, each 5x5x5 elements
//            n_subdomOnCluster = 27;
//            path2data = "../data3/";
//            break;
//        case (4):
//            // 2nd larger variant: 27 subdom, each 5x5x5 elements
//            // 'kernel' gluing
//            n_subdomOnCluster = 27;
//            path2data = "../data4/";
//            break;
//        case (5):
//            // data are created in application
//            // and 'path2data' is setup for dumped files
//            path2data = "../data5/";
//            break;
//    }

//
    Cluster cluster(options);
    cout << "----------------- done -----------------\n" ;
    return 0;
}
