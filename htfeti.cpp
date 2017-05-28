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
    int typeBc = 0;    // 0: corners, 1: zero and first approx. (ker), 2: all nodes on interf.)
    bool Ac_extended_by_kerGc = false;
    bool GcTGc_assembl_block_by_block= true;

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

    cout << argv[0] << endl;
    Options options;
    options.set_values(path2data, argc, argv,
                        young_modulus, poissons_ratio,
                        pardiso_0_dissection_1,print_matrices,
                       typeBc, Ac_extended_by_kerGc, GcTGc_assembl_block_by_block);
    Cluster cluster(options);
    cout << "----------------- done -----------------\n" ;
    return 0;
}
