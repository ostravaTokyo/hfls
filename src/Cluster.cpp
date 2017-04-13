#include "Cluster.h"

using namespace std;

Cluster::Cluster(Options options)
{

    K.resize(options.n_subdomOnCluster);
    RegMat.resize(options.n_subdomOnCluster);
    R.resize(options.n_subdomOnCluster);
    rhs.resize(options.n_subdomOnCluster);
    Bc.resize(options.n_subdomOnCluster);
    Bct.resize(options.n_subdomOnCluster);
    Bf.resize(options.n_subdomOnCluster);
    K_regularized.resize(options.n_subdomOnCluster);
    Fc_sub.resize(options.n_subdomOnCluster);
    Gc.resize(options.n_subdomOnCluster);

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*                         READING DATA (TXT)                              */
    /*                 matrices in coo format are one-based                    */
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    for (int i = 0 ; i < options.n_subdomOnCluster ; i++){
        std::cout << "\t\t\t====(" << i+1 << "/" <<
                     options.n_subdomOnCluster <<")====\n";
        /* K - stiffness matrix */
        std::string path2matrix = options.path2data+"/K"+to_string(i)+".txt";
        this->K[i].readCooFromFile(path2matrix);
        this->K[i].COO2CSR();

        /* R - kernel matrix */
        path2matrix = options.path2data+"/R1"+to_string(i)+".txt";
        this->R[i].readCooFromFile(path2matrix);
        this->R[i].COO2CSR();
        int two_means_not_lower_not_upper = 2;
        this->R[i].CSR2DNS(two_means_not_lower_not_upper);

        /* rhs - right-hand-side vector */
        path2matrix = options.path2data+"/R1"+to_string(i)+".txt";
        this->rhs[i].readCooFromFile(path2matrix,this->K[i].n_row);

        /* Bf - constraints matrix */
        path2matrix = options.path2data+"/B1"+to_string(i)+".txt";
        this->Bf[i].readCooFromFile(path2matrix);
        this->Bf[i].compresRows();
        this->Bf[i].COO2CSR();

        /* Bc - constraints matrix */
        path2matrix = options.path2data+"/B0"+to_string(i)+".txt";
        this->Bc[i].readCooFromFile(path2matrix);
        this->Bc[i].compresRows();

        /* create Bct as copy of Bc */
        this->Bct[i].CreateCopyFrom(this->Bc[i],0);
        this->Bc[i].COO2CSR();

//        this->Bct[i].COO2CSR();

        /* Regularization of sing. stiff. mat K*/
        path2matrix = options.path2data+"/RegMat"+to_string(i)+".txt";
        this->RegMat[i].readCooFromFile(path2matrix);
        this->RegMat[i].COO2CSR();
        SparseMatrix::dcsradd(this->K[i], this->RegMat[i], this->K_regularized[i]);
        this->K_regularized[i].RemoveLower();
    }


//    SparseMatrix::testPardiso();

    int  isUpperLowerOrBoth = 2;
//    for (int i = 0 ; i < options.n_subdomOnCluster ; i++){
//        SparseMatrix KplusBct;
//        this->Bct[i].COO2CSR();
//        this->Bct[i].printToFile("Bct_beforPardiso",i);
//        this->Bct[i].CSR2DNS(isUpperLowerOrBoth);
//        this->K_regularized[i].InitializeSolve();
//        this->K_regularized[i].solve(this->Bct[i], KplusBct);
//        this->Bct[i].DNS2CSR(isUpperLowerOrBoth);
//        this->Bct[i].printToFile("Bct_beforPardiso",i);
//        KplusBct.DNS2CSR(isUpperLowerOrBoth);
//        KplusBct.printToFile("KplusBct",i);
//        KplusBct.CSR2COO();
//        KplusBct.CSR2DNS(isUpperLowerOrBoth);
//
//        //      SparseMatrix _Fc_sub;
//        SparseMatrix::Acsr_mult_Bdns_is_Cdns(this->Bc[i],KplusBct,this->Fc_sub[i]);
//        this->Fc_sub[i].DNS2CSR(2);
//        SparseMatrix::Acsr_mult_Bdns_is_Cdns(this->Bc[i],this->R[i],this->Gc[i]);
//
//
//        std::cout << " ++++++     size  +++++     " << KplusBct.adns.size() << "   \n";
//        for (int ii = 0 ; ii < 20 ; ii ++ ){
//            std::cout << "  +++++++++++++++  "<< KplusBct.adns[ii] << "   \n";
//        }

 //   }

    /* dumping matrices */
    for (int i = 0; i < options.n_subdomOnCluster ; i++ ){
        this->K_regularized[i].printToFile("K_regularized",i);
        this->RegMat[i].printToFile("RegMat",i);
        this->K[i].printToFile("K",i);
        this->Bc[i].printToFile("Bc",i);
 //       this->Bct[i].printToFile("Bct",i);
        this->Bf[i].printToFile("Bf",i);
        this->R[i].printToFile("R",i);
    }



    for (int i = 0 ; i < options.n_subdomOnCluster ; i++){
        this->K_regularized[i].FinalizeSolve(i);
    }

}
