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
    const int nS = options.n_subdomOnCluster;
    K.resize(nS);
//    RegMat.resize(nS);
    R.resize(nS);
    rhs.resize(nS);
    Bc.resize(nS);
    Bf.resize(nS);
    Bct.resize(nS);
//  K_regularized.resize(nS);
//    Fc_sub.resize(nS);
//    Gc.resize(nS);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                         READING DATA (TXT)                              */
/*                 matrices in coo format are one-based                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


    for (int i = 0 ; i < nS ; i++){
//        std::cout << "\t\t\t====(" << i+1 << "/" <<
//                   nS <<")====\n";

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




        path2matrix = options.path2data+"/B0"+to_string(i)+".txt";
        symmetric = 0; format = 1; offset = 1;
//        this->Bct[i].readCooFromFile(path2matrix,symmetric,format,offset);
        this->Bct[i] = this->Bc[i];
        reduceZeroRows = true;
        transpose = true;
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
// //      for (int i = 0 ; i < nS ; i++){
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

//for (int i = 0 ; i < this->K.size(); i ++ ){
//    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
//    cout << "this->Bct[0].i_ptr[k].size() = " << this->Bct[i].n_row_cmprs << endl;
//    cout << "this->Bct[0].j_col[k].size() = " << this->Bct[i].n_col  << endl;
//    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
//    cout <<"\n";
//
//
//    for (int k = 0; k < this->Bct[i].dense.size(); k++){
//        if (this->Bct[i].dense[k] != 0)
//            cout << this->Bct[i].dense[k] << "/"<< endl;
//    }
//}




#if 1
    int i_sub = 3;

    for (int i_sub = 0; i_sub < nS; i_sub++){
    double *Y = new double[this->K[i_sub].n_row_cmprs * this->R[i_sub].n_col];
    for (int i = 0 ; i < this->K[i_sub].n_row_cmprs; i++){
       Y[i] = 0;
    }
//    int nn = this->K[0].n_row_cmprs;
//    this->K[0].mv_csr(&(this->Bct[0].dense[0]),Y,true,1);
    this->K[i_sub].mv_csr(&(this->R[i_sub].dense[0]),Y,true,6);
//    for (int i = 0 ; i < this->K[0].n_row_cmprs; i++){
//        printf("%d: %f \n",i,Y[i]);
//    }



    double normK = 0;
    for (int i = 0 ; i < this->K[i_sub].nnz; i++){
       normK += this->K[i_sub].val[i] * this->K[i_sub].val[i];
    }
//    normK /= this->K[0].nnz;
    double normR = 0;
    for (int i = 0 ; i < this->R[i_sub].dense.size(); i++){
       normR += this->R[i_sub].dense[i] * this->R[i_sub].dense[i];
    }

    double normKR = 0;
    for (int i = 0 ; i < this->K[i_sub].n_row_cmprs; i++){
        for (int j = 0; j < this->R[i_sub].n_col; j++){
            if (i_sub == -1){
                if (j==0){
                    printf("%d: ",i);
                }
                printf(" %f\t", Y[i +  j * this->K[i_sub].n_row_cmprs]);
            }
            normKR += Y[i +  j * this->K[i_sub].n_row_cmprs] *
                        Y[i +  j * this->K[i_sub].n_row_cmprs];
        }
        if (i_sub == -1){
            printf("\n");
        }
    }

    printf("||K||                 = %.16e ", sqrt(normK));
    printf("||K||                 = %.16e ", sqrt(normR));
    printf("||KR|| /(||K|| ||R||) = %.16e \n", sqrt(normKR) / (sqrt(normK*normR)));
    }


#else
//    Matrix Y = this->Bct[0];
//    printCooOrDense = true;
//    Y.printToFile("K_Bct",0,printCooOrDense);



#endif


    for (int i = 0; i < nS ; i++ ){
        printCooOrDense = true;
        this->K[i].printToFile("K",i,printCooOrDense);
        printCooOrDense = true;
        this->R[i].printToFile("R",i,printCooOrDense);
        printCooOrDense = true;
        this->Bc[i].printToFile("Bc",i,printCooOrDense);
        printCooOrDense = true;
        this->Bct[i].printToFile("Bct",i,printCooOrDense);
        printCooOrDense = true;
        this->Bf[i].printToFile("Bf",i,printCooOrDense);

//       this->K_regularized[i].printToFile("K_regularized",i);
//        this->RegMat[i].printToFile("RegMat",i);
     }
//
//
//
//          //    for (int i = 0 ; i < nS ; i++){
//          //        this->K_regularized[i].FinalizeSolve(i);
//          //    }
//
}
