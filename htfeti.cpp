#include <iostream>
#include <string>
#include "src/Options.hpp"
#include "src/Cluster.hpp"
#include <map>



using namespace std;
int main(int argc, char *argv[]){


    std::map<string,string > options2;


    options2["path2data"]                       =  "../data/";
    options2["young_modulus"]                   =  "1000";
    options2["poissons_ratio"]                  =  "0.3";
    options2["linear_solver"]                   =  "dissection";    /* dissection | pardiso */
    options2["print_matrices"]                  =  "0";             /* 0 | 1 | 2 | 3        */
    options2["typeBc"]                          =  "ker";           /* ker | cor | all      */
    options2["Ac_extended_by_kerGc"]            =  "false";          /* true | false         */
    options2["GcTGc_assembl_block_by_block"]    =  "true";          /* true | false         */
    options2["Bc_fullRank"]                     =  "true";         /* true | false         */
    options2["create_analytic_ker_K"]           =  "false";         /* true | false         */

    options2["Nx"]                              =  "2";
    options2["Ny"]                              =  "2";
    options2["Nz"]                              =  "2";
    options2["nx"]                              =  "5";
    options2["ny"]                              =  "5";
    options2["nz"]                              =  "5";



    /* data printed into */
    string path2data = "../data/";
                /* material constants */
    double young_modulus = 1000;
    double poissons_ratio = 0.3;
                /* linear solver */
    int pardiso_0_dissection_1 = 1;
                /* type of Bc matrix
                 * 0: corners,
                 * 1: zero and first approx. (ker),
                 * 2: all nodes on interf.)*/
    int typeBc = 1;

                /* if true, matrix Bc (corners) is full column rank */
                /*          matrix Bc (zero and first approx.) is full column rank matrix if
                 *          each subdomain consists of at least two or more elements */
    bool Bc_fullRank = false;
                /* dumping matrices in MatrixMarket format */
                /* print_matrices  = 0 : do nothing */
                /* print_matrices  = 1 : clust. objects (GcTGc, Fc, Ac)*/
                /* print_matrices  = 2 : print all previous + K, Bc, Fc_sub, Gc_sub ...  */
                /* print_matrices  = 3 : print all previous + BcT_dense, BcT_dense_new, ...  */
    int print_matrices = 0;
                /* Ac matrix with or withour regularization */
    bool Ac_extended_by_kerGc = true;

                /* sparse-BLOCK or dense assembling of GcTGc matrix */
    bool GcTGc_assembl_block_by_block= true;
    bool create_analytic_ker_K = false;

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
//

    if (pardiso_0_dissection_1 == 0 && !create_analytic_ker_K){
        cout <<"###########################################################" << endl;
        cout << "if pardiso, kernel must be provided by analytical formula" << endl;
        cout <<"###########################################################" << endl;
        create_analytic_ker_K = true;
    }




    if (argc > 1) {
       options2["Nx"]   =  (argv[1]);
    }
    if (argc > 2) {
       options2["Ny"]   =  (argv[2]);
    }
    if (argc > 3) {
       options2["Nz"]   =  (argv[3]);
    }
    if (argc > 4) {
       options2["nx"]   =  (argv[4]);
    }
    if (argc > 5) {
       options2["ny"]   =  (argv[5]);
    }
    if (argc > 6) {
       options2["nz"]   =  (argv[6]);
    }



    cout << argv[0] << endl;
    Options options;
//    options.set_values(path2data, argc, argv,
//                        young_modulus, poissons_ratio,
//                        pardiso_0_dissection_1,print_matrices,
//                       typeBc, Ac_extended_by_kerGc, GcTGc_assembl_block_by_block,
//                       Bc_fullRank,create_analytic_ker_K);
    Cluster cluster(options,options2);
    cout << "----------------- done -----------------\n" ;
    return 0;
}
