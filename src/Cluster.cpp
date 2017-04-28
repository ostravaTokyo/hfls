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
    R.resize(nS);
    rhs.resize(nS);
    Bc.resize(nS);
    Bf.resize(nS);
    Bct.resize(nS);
//    Fc_sub.resize(nS);
//    Gc.resize(nS);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                         READING DATA (TXT)                              */
/*                 matrices in coo format are one-based                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
 
    for (int i = 0 ; i < nS ; i++){
        /* K - stiffness matrix */
        string path2matrix = options.path2data+"/K"+to_string(i)+".txt";
        symmetric = 2; format = 1; offset = 1;
        K[i].readCooFromFile(path2matrix,symmetric,format,offset);

        /* R - kernel matrix */
        path2matrix = options.path2data+"/R1"+to_string(i)+".txt";
        symmetric = 0; format = 2; offset = 1;
        R[i].readCooFromFile(path2matrix,symmetric,format,offset);

        /* rhs - right-hand-side vector */
        path2matrix = options.path2data+"/f"+to_string(i)+".txt";
        rhs[i].readCooFromFile(path2matrix,K[i].n_row);
//
        /* Bf - constraints matrix */
        path2matrix = options.path2data+"/B1"+to_string(i)+".txt";
        symmetric = 0; format = 1; offset = 1;
        Bf[i].readCooFromFile(path2matrix,symmetric,format,offset);
//
        /* Bc - constraints matrix */
        path2matrix = options.path2data+"/B0"+to_string(i)+".txt";
        symmetric = 0; format = 1; offset = 1;
        Bc[i].readCooFromFile(path2matrix,symmetric,format,offset);

        Bct[i] = Bc[i];
        reduceZeroRows = true;
        transpose = true;
        Bct[i].CSRorCOO2DNS(reduceZeroRows,transpose);

  }

// //      Matrix::testPardiso();


#if 1
    int i_sub = 3;

    double *Y = new double[K[i_sub].n_row_cmprs * R[i_sub].n_col];
    for (int i_sub = 0; i_sub < nS; i_sub++){
        for (int i = 0 ; i < K[i_sub].n_row_cmprs; i++){
           Y[i] = 0;
        }
        K[i_sub].mv_csr(&(R[i_sub].dense[0]),Y,true,6);

        double normK = 0;
        for (int i = 0 ; i < K[i_sub].nnz; i++){
           normK += K[i_sub].val[i] * K[i_sub].val[i];
        }
    //    normK /= K[0].nnz;
        double normR = 0;
        for (unsigned int i = 0 ; i < R[i_sub].dense.size(); i++){
           normR += R[i_sub].dense[i] * R[i_sub].dense[i];
        }

        double normKR = 0;
        for (int i = 0 ; i < K[i_sub].n_row_cmprs; i++){
            for (int j = 0; j < R[i_sub].n_col; j++){
                if (i_sub == -1){
                    if (j==0){
                        printf("%d: ",i);
                    }
                    printf(" %f\t", Y[i +  j * K[i_sub].n_row_cmprs]);
                }
                normKR += Y[i +  j * K[i_sub].n_row_cmprs] *
                            Y[i +  j * K[i_sub].n_row_cmprs];
            }
            if (i_sub == -1){
                printf("\n");
            }
        }

        printf("||K|| = %.6e, ", sqrt(normK));
        printf("||R|| = %.6e, ", sqrt(normR));
        printf("||KR|| = %.6e, ",sqrt(normKR) );
        printf("||KR|| /(||K|| ||R||) = %.6e \n", sqrt(normKR) / (sqrt(normK*normR)));
    }
    delete [] Y;
#endif


    i_sub = 0;
    Matrix Gc;
    int n_rowGc = Bc[i_sub].n_row_cmprs;
    int n_colGc = R[i_sub].n_col;
    Gc.zero_dense(n_rowGc, n_colGc );

    Bc[i_sub].mv_csr(&(R[i_sub].dense[0]),&(Gc.dense)[0],true ,R[i_sub].n_col);
//
    printCooOrDense = true;

    Gc.printToFile("Gc",0,printCooOrDense);


    for (int i = 0; i < nS ; i++ ){
        printCooOrDense = true;
        K[i].printToFile("K",i,printCooOrDense);
        printCooOrDense = true;
        R[i].printToFile("R",i,printCooOrDense);
        printCooOrDense = true;
        Bc[i].printToFile("Bc",i,printCooOrDense);
        printCooOrDense = true;
        Bct[i].printToFile("Bct",i,printCooOrDense);
        printCooOrDense = true;
        Bf[i].printToFile("Bf",i,printCooOrDense);
     }
}
