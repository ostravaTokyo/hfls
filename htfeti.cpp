#include <iostream>
#include <string> 
#include "src/Options.hpp"
#include "src/Cluster.hpp"



using namespace std;
int main(int argc, char *argv[]){

    /* data printed into */
    string path2data = "../data/";
    /* material constants */
    double young_modulus = 10000;
    double poissons_ratio = 0.3;
    /* linear solver */
    double pardiso_0_dissection_1 = 1;
    int print_matrices = 1;

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

    cout << argv[0] << endl;
    Options options;
    options.set_values(path2data, argc, argv,
                        young_modulus, poissons_ratio,
                        pardiso_0_dissection_1,print_matrices);
    Cluster cluster(options);
    cout << "----------------- done -----------------\n" ;
    return 0;
}
