#include "../include/Cluster.hpp"

using namespace std;

Cluster::Cluster(Options options)
{


int symmetric, format, offset;
bool reduceZeroRows, transpose, printCooOrDense;

#if 0
    std::string path2matrix = options.path2data+"/testMat"+to_string(0)+".txt";
    Matrix A;

    symmetric = 0;  // 0-not, 1-lower tr., 2-upper tr.
    format = 1;     // 0-coo, 1-csr, 2-dense
    offset = 0;

    printCooOrDense = true;
    A.readCooFromFile(path2matrix,symmetric,format,offset);
    A.printToFile("modif",0,printCooOrDense);
    A.CSR2COO();
    A.printToFile("modif2",0,printCooOrDense);
    printCooOrDense = true;
    A.printToFile("modif3",0,printCooOrDense);

    return

#endif
//    Matrix B;
//    B = B.CreateCopyFrom(A);
//    B.COO2CSR();
//    B.printToFile("modif3",0,printCooOrDense);
//    bool reduceZeroRows = true;
//    bool transpose = true;
//    B.CSRorCOO2DNS(reduceZeroRows, transpose);
    //
    K.resize(options.n_subdomOnCluster);
//    RegMat.resize(options.n_subdomOnCluster);
    R.resize(options.n_subdomOnCluster);
    rhs.resize(options.n_subdomOnCluster);
    Bc.resize(options.n_subdomOnCluster);
    Bf.resize(options.n_subdomOnCluster);
    Bct.resize(options.n_subdomOnCluster);
//  K_regularized.resize(options.n_subdomOnCluster);
//    Fc_sub.resize(options.n_subdomOnCluster);
//    Gc.resize(options.n_subdomOnCluster);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                         READING DATA (TXT)                              */
/*                 matrices in coo format are one-based                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


    for (int i = 0 ; i < options.n_subdomOnCluster ; i++){
//        std::cout << "\t\t\t====(" << i+1 << "/" <<
//                   options.n_subdomOnCluster <<")====\n";

        /* K - stiffness matrix */
        string path2matrix = options.path2data+"/K"+to_string(i)+".txt";
        symmetric = 2; format = 1; offset = 1;
        this->K[i].readCooFromFile(path2matrix,symmetric,format,offset);

        /* R - kernel matrix */
        path2matrix = options.path2data+"/R1"+to_string(i)+".txt";
        symmetric = 0; format = 2; offset = 1;
        this->R[i].readCooFromFile(path2matrix,symmetric,format,offset);

//    if (_format == 2){
//    }
//    else{
        /* rhs - right-hand-side vector */
        path2matrix = options.path2data+"/f"+to_string(i)+".txt";
        this->rhs[i].readCooFromFile(path2matrix,this->K[i].n_row);
//
        /* Bf - constraints matrix */
        path2matrix = options.path2data+"/B1"+to_string(i)+".txt";
        symmetric = 0; format = 1; offset = 1;
        this->Bf[i].readCooFromFile(path2matrix,symmetric,format,offset);
//
        /* Bc - constraints matrix */
        path2matrix = options.path2data+"/B0"+to_string(i)+".txt";
        symmetric = 0; format = 1; offset = 1;
        this->Bc[i].readCooFromFile(path2matrix,symmetric,format,offset);


        reduceZeroRows = true;
        transpose = true;
        this->Bct[i] = this->Bc[i];
        this->Bct[i].CSRorCOO2DNS(reduceZeroRows,transpose);

//
//
//        /* Regularization of sing. stiff. mat K*/
//        path2matrix = options.path2data+"/RegMat"+to_string(i)+".txt";
//        this->RegMat[i].readCooFromFile(path2matrix,2,1,0);
//         Matrix::dcsradd(this->K[i], this->RegMat[i], this->K_regularized[i]);
//          this->K_regularized[i].RemoveLower();
  }
//
//
// //      Matrix::testPardiso();
//
//      int  isUpperLowerOrBoth = 2;
// //      for (int i = 0 ; i < options.n_subdomOnCluster ; i++){
// //          Matrix KplusBct;
// //          this->Bct[i].COO2CSR();
// //          this->Bct[i].printToFile("Bct_beforPardiso",i);
// //          this->Bct[i].CSR2DNS(isUpperLowerOrBoth);
// //          this->K_regularized[i].InitializeSolve();
// //          this->K_regularized[i].solve(this->Bct[i], KplusBct);
// //          this->Bct[i].DNS2CSR(isUpperLowerOrBoth);
// //          this->Bct[i].printToFile("Bct_beforPardiso",i);
// //          KplusBct.DNS2CSR(isUpperLowerOrBoth);
// //          KplusBct.printToFile("KplusBct",i);
// //          KplusBct.CSR2COO();
// //          KplusBct.CSR2DNS(isUpperLowerOrBoth);
// //
// //          //      Matrix _Fc_sub;
// //          Matrix::Acsr_mult_Bdns_is_Cdns(this->Bc[i],KplusBct,this->Fc_sub[i]);
// //          this->Fc_sub[i].DNS2CSR(2);
// //          Matrix::Acsr_mult_Bdns_is_Cdns(this->Bc[i],this->R[i],this->Gc[i]);
// //
// //
// //          std::cout << " ++++++     size  +++++     " << KplusBct.adns.size() << "   \n";
// //          for (int ii = 0 ; ii < 20 ; ii ++ ){
// //              std::cout << "  +++++++++++++++  "<< KplusBct.adns[ii] << "   \n";
// //          }
//
//             //   }
//
//                /* dumping matrices */



#if 0
    double *X = new double[this->K[0].n_row_cmprs];
    double *Y = new double[this->K[0].n_row_cmprs];
    for (int i = 0 ; i < this->K[0].n_row_cmprs; i++){
       X[i] = 1;
       Y[i] = 0;
    }
    this->K[0].mv_csr(X,Y,true,1);
    for (int i = 0 ; i < this->K[0].n_row_cmprs; i++){
        printf("%f %f \n", X[i], Y[i]);
    }
#else
    Matrix Y = this->Bct[0];
    printCooOrDense = true;
    Y.printToFile("K_Bct",0,printCooOrDense);



#endif


    for (int i = 0; i < options.n_subdomOnCluster ; i++ ){
        printCooOrDense = true;
        this->K[i].printToFile("K",i,printCooOrDense);
        printCooOrDense = true;
        this->R[i].printToFile("R",i,printCooOrDense);
//
//        printCooOrDense = true;
//        this->Bc[i].printToFile("Bc",i,printCooOrDense);
//        printCooOrDense = true;
//        this->Bf[i].printToFile("Bf",i,printCooOrDense);
//        printCooOrDense = true;
//        this->Bf[i].printToFile("Bct",i,printCooOrDense);

//       this->K_regularized[i].printToFile("K_regularized",i);
//        this->RegMat[i].printToFile("RegMat",i);
     }
//
//
//
//          //    for (int i = 0 ; i < options.n_subdomOnCluster ; i++){
//          //        this->K_regularized[i].FinalizeSolve(i);
//          //    }
//
}
