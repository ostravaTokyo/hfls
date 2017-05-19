#include <iostream>
#include <string> 
#include "src/Options.hpp"
#include "src/Cluster.hpp"



using namespace std;
int main(int argc, char *argv[]){


    string path2data = "../data/";
/* norm of matrix can be changed by following 2 parameters */
    double young_modulus = 10000;
    double poissons_ratio = 0.3;


//    double young_modulus = 2.1e9;
//    double poissons_ratio = 0.4999;

/* ---------------------------------------------------------*/
/* ---------------------------------------------------------*/
/* ---------------------------------------------------------*/
/* ---------------------------------------------------------*/


    cout << argv[0] << endl;
    Options options;
    options.set_values(path2data, argc, argv, young_modulus, poissons_ratio);
    Cluster cluster(options);
    cout << "----------------- done -----------------\n" ;
    return 0;
}
