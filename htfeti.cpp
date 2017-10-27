#include <iostream>
#include <string>
#include "src/Options.hpp"
#include "src/Cluster.hpp"
#include <map>



using namespace std;
int main(int argc, char *argv[]){

    map<string,string > options2;

    options2["eps_iter"]                        = "1e-4";
    /* Dirichlet_explicit, Dirichlet_implicit, lumped */
    options2["preconditioner"]                  = "Dirichlet_implicit";
    //options2["preconditioner"]                  = "Dirichlet_explicit";
    //options2["preconditioner"]                  = "lumped";
    options2["path2data"]                       =  "data/";
    options2["young_modulus"]                   =  "1000";
    options2["poissons_ratio"]                  =  "0.3";

    /* dissection | pardiso */
    options2["linear_solver"]                   =  "dissection";

    /* {0, 1, 2, 3 }                                            */
    options2["print_matrices"]                  =  "0";

    /* ker, cor, all                                            */
    options2["typeBc"]                          =  "cor";

    /* true,  false                                             */
    options2["Ac_extended_by_kerGc"]            =  "false";

    /* true,  false                                             */
    options2["GcTGc_assembl_block_by_block"]    =  "true";

    /* true, false                                              */
    options2["Bc_fullRank"]                     =  "true";

    /* true, false                                              */
    options2["create_analytic_ker_K"]           =  "false";



    options2["Nx"]                              =  "2";
    options2["Ny"]                              =  "2";
    options2["Nz"]                              =  "2";
    options2["nx"]                              =  "5";
    options2["ny"]                              =  "5";
    options2["nz"]                              =  "5";


    printf("+++++++++++++++++++++++++++++++++++      %s\n", options2["path2data"].c_str());


    if (options2["linear_solver"] ==  "pardiso" &&
            options2["create_analytic_ker_K"] == "false"){
        cout <<"###########################################################" << endl;
        cout << "if pardiso, kernel must be provided by analytical formula" << endl;
        cout <<"###########################################################" << endl;
        options2["create_analytic_ker_K"] = "true";
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
    Cluster cluster(options,options2);
    cout << "----------------- done -----------------\n" ;


// HOW TO PRINT OPTIONS2

//    cout << "\t\t\toptions2.size() "  << options2.size() << endl;
//    for (map<string,string>::const_iterator it=options2.begin(); it!=options2.end(); ++it)
//        cout << "\t\t\t" << it->first << " => " << it->second << '\n';



    return 0;
}

