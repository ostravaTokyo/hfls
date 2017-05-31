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

                /* type of Bc matrix
                 * 0: corners,
                 * 1: zero and first approx. (ker),
                 * 2: all nodes on interf.)*/
    int typeBc = 0;

                /* dumping matrices in MatrixMarket format */
    int print_matrices = 0;

                /* Ac matrix with or withour regularization */
    bool Ac_extended_by_kerGc = false;

                /* sparse-BLOCK or dense assembling of GcTGc matrix */
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
