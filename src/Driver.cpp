#include "Driver.hpp"
#include <ctime>



using namespace std;

Driver::Driver()
{
   cout << "constructor of Diver" <<endl;
}


Driver::Driver(map <string,string> options2_)
{
    options2 = options2_;
}


void Driver::initialization(){
//    clock_t begin = clock();
    //options = options_;
//
//    bool reduceZeroRows, transpose, printCooOrDense, checkOrthogonality;
//    int printMat = atoi(options2["print_matrices"].c_str());
//    folder = options2["path2data"];
    printf("+++++++++++++++++++++++++++++++++++ cluster     %s\n", options2["path2data"].c_str());
//
    /* mesh belonging to i-th cluster (currently, i=0 only)*/
    mesh.createMesh(options2);
    mesh.ddm_metis(options2);
    cluster.initialization(options2,mesh);
//    solver.pcpg(cluster);
//    int nSubClst = cluster.get_nSubClst();



}

