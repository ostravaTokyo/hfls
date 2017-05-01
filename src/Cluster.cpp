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
    K_reg.resize(nS);
    Lumped.resize(nS);
    R.resize(nS);
    rhs.resize(nS);
    Bc.resize(nS);
    Bf.resize(nS);
    Bc_dense.resize(nS);
    Fc.resize(nS);
    Gc.resize(nS);

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

//
        /* Bc_dense - constraints matrix */
        Bc_dense[i] = Bc[i];
        reduceZeroRows = true;
        transpose = true;
        Bc_dense[i].CSRorCOO2DNS(reduceZeroRows,transpose);

  }


    for (int i_sub = 0; i_sub < nS; i_sub++){
        Matrix Y;
        Y.zero_dense(K[i_sub].n_row_cmprs , R[i_sub].n_col);
        K[i_sub].mult(R[i_sub],Y,true);

        double normK = K[i_sub].norm2();
        double normR = R[i_sub].norm2();
        double normKR = Y.norm2();

        printf("|K| = %.6e, ", normK);
        printf("|R| = %.6e, ", normR);
        printf("|KR|/(|K|*|R|) = %.6e \n", normKR / (normK*normR));
    }

    printCooOrDense = true;
    for (int i = 0; i < nS ; i++ ){
        K_reg[i] = K[i];
        K_reg[i].factorization();

        Matrix BcK_dense;
        BcK_dense = Matrix::CreateCopyFrom(Bc_dense[i]);
        BcK_dense.setZero();


        K[i].mult(Bc_dense[i],BcK_dense,true);
        BcK_dense.printToFile("BcK_dense",i,printCooOrDense);
        Lumped[i].zero_dense(Bc[i].n_row_cmprs,  Bc[i].n_row_cmprs);
        Bc[i].mult(BcK_dense,Lumped[i],true);
        Lumped[i].printToFile("Lumped",i,printCooOrDense);

        Matrix BcKplus_dense;
        BcKplus_dense = Matrix::CreateCopyFrom(Bc_dense[i]);
        BcKplus_dense.setZero();

        K_reg[i].solve(Bc_dense[i],BcKplus_dense);
        BcKplus_dense.printToFile("BcKplus",i,printCooOrDense);
        Fc[i].zero_dense(Bc[i].n_row_cmprs,  Bc[i].n_row_cmprs);
        Bc[i].mult(BcKplus_dense,Fc[i],true);
        Fc[i].printToFile("Fc",i,printCooOrDense);

        /* Gc - constraints matrix */
        int n_rowGc = Bc[i].n_row_cmprs;
        int n_colGc = R[i].n_col;
        Gc[i].zero_dense(n_rowGc, n_colGc );
        Bc[i].mult(R[i],Gc[i],true);

        K[i].printToFile("K",i,printCooOrDense);
        K_reg[i].printToFile("K_reg",i,printCooOrDense);
        R[i].printToFile("R",i,printCooOrDense);
        Bc[i].printToFile("Bc",i,printCooOrDense);
        Bc_dense[i].printToFile("Bc_dense",i,printCooOrDense);
        Bf[i].printToFile("Bf",i,printCooOrDense);
        Gc[i].printToFile("Gc",i,printCooOrDense);


    }

//    Matrix::testPardiso();



    for (int i = 0 ; i < nS; i++){
       K_reg[i].FinalizeSolve(i);
    }
}
