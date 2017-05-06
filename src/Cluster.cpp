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
    A.printToFile("modif",folder,0,printCooOrDense);
    A.CSR2COO();
    A.printToFile("modif2",folder,0,printCooOrDense);
    printCooOrDense = true;
    A.printToFile("modif3",folder,0,printCooOrDense);

    return

#endif
    const int nS = options.n_subdomOnCluster;
    string folder = options.path2data;

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


#if 0
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
#endif
    printCooOrDense = true;
    for (int i = 0; i < nS ; i++ ){

        vector <int > nullPivots;
        R[i].getNullPivots(nullPivots);
//        for (int i = 0 ; i < nullPivots.size(); i++){
//            cout << nullPivots[i] << endl;
//        }

        K_reg[i] = K[i];
        K_reg[i].factorization(nullPivots);



        Matrix BcK_dense;
        //BcK_dense = Matrix::CreateCopyFrom(Bc_dense[i]);
        BcK_dense = Bc_dense[i];
        BcK_dense.setZero();


        K[i].mult(Bc_dense[i],BcK_dense,true);
        BcK_dense.printToFile("BcK_dense",folder,i,printCooOrDense);
        Lumped[i].zero_dense(Bc[i].n_row_cmprs,  Bc[i].n_row_cmprs);
        Bc[i].mult(BcK_dense,Lumped[i],true);
        Lumped[i].printToFile("Lumped",folder,i,printCooOrDense);

        Matrix BcKplus_dense;
//        BcKplus_dense = Matrix::CreateCopyFrom(Bc_dense[i]);
        BcKplus_dense = Bc_dense[i];
        BcKplus_dense.setZero();

        K_reg[i].solve(Bc_dense[i],BcKplus_dense);
        BcKplus_dense.printToFile("BcKplus",folder,i,printCooOrDense);
        Fc[i].zero_dense(Bc[i].n_row_cmprs,  Bc[i].n_row_cmprs);
        Bc[i].mult(BcKplus_dense,Fc[i],true);
        Fc[i].printToFile("Fc",folder,i,printCooOrDense);
        Fc[i].l2g_i_coo = Bc[i].l2g_i_coo;


        /* Gc - constraints matrix */
        int n_rowGc = Bc[i].n_row_cmprs;
        int n_colGc = R[i].n_col;
        Gc[i].zero_dense(n_rowGc, n_colGc );
        Bc[i].mult(R[i],Gc[i],true);


        /* print */
        /* K, R, etc.*/
        cout << "   ... printing matrices start  ...";
        K[i].printToFile("K",folder,i,printCooOrDense);
        K_reg[i].printToFile("K_reg",folder,i,printCooOrDense);
        R[i].printToFile("R",folder,i,printCooOrDense);
        Bc[i].printToFile("Bc",folder,i,printCooOrDense);
        Bc_dense[i].printToFile("Bc_dense",folder,i,printCooOrDense);
        Bf[i].printToFile("Bf",folder,i,printCooOrDense);
        Gc[i].printToFile("Gc",folder,i,printCooOrDense);
        Gc[i].l2g_i_coo = Bc[i].l2g_i_coo;
        cout << "  " << (i + 1) <<"/" << nS <<"  ...\n";
    }

//    Matrix::testPardiso();

    /* Fc_clust  */
    create_Fc_clust();
    Fc_clust.printToFile("Fc_clust",folder,0,printCooOrDense);

    create_Gc_clust();
    Gc_clust.printToFile("Gc_clust",folder,0,printCooOrDense);


    for (int i = 0 ; i < nS; i++){
       K_reg[i].FinalizeSolve(i);
    }
}

void Cluster::create_Fc_clust(){

    Matrix & A_clust = Fc_clust;
    vector <Matrix> & A_i = Fc;
    A_clust.symmetric = 2;
    bool remapCols = true;
    create_clust_object(A_clust, A_i, remapCols);
}

void Cluster::create_Gc_clust(){

    Matrix & A_clust = Gc_clust;
    vector <Matrix> & A_i = Gc;
    A_clust.symmetric = 0;
    bool remapCols = false;
    create_clust_object(A_clust, A_i, remapCols);
}


void Cluster::create_clust_object(Matrix &A_clust, vector <Matrix> & A_i, bool remapCols){

    const int nS = A_i.size();
    int init_nnz = 0;
    for (int i = 0; i < nS; i++){
        init_nnz += A_i[i].nnz;
    }

    A_clust.format = 0;

    vector < int_int_dbl > tmpVec;
    tmpVec.resize(init_nnz);

    A_clust.n_col = 0;


/* Sorting [I, J, V]  --------------------------------------------------------*/
    int cnt = 0;
    int _i, _j;
    for (int d = 0; d < nS ; d++){
        A_i[d].DNS2COO();
        for (int i = 0; i < A_i[d].i_coo_cmpr.size();i++){
            _i = A_i[d].l2g_i_coo[ A_i[d].i_coo_cmpr[i] ];
            _j = A_i[d].j_col[i];
            if (remapCols){
                _j = A_i[d].l2g_i_coo[ _j ];
            }
            else{
                _j += A_clust.n_col;
            }
            if (A_clust.symmetric == 0 ||
                    (A_clust.symmetric == 1 && _j <= _i) ||
                    (A_clust.symmetric == 2 && _i <= _j)    ){
                tmpVec[cnt].I = _i;
                tmpVec[cnt].J = _j;
                tmpVec[cnt].V = A_i[d].val[i];
                cnt++;
            }
        }
        A_clust.n_col += A_i[d].n_col;
    }
    A_clust.nnz = cnt;
    tmpVec.resize(A_clust.nnz);

    A_clust.sortAndUniqueCOO(tmpVec);

    A_clust.n_row = A_clust.n_row_cmprs;



    if (A_clust.symmetric > 0){
        A_clust.n_col = A_clust.n_row_cmprs;
    }

    A_clust.COO2CSR();

    tmpVec.clear();
    tmpVec.shrink_to_fit();
}
