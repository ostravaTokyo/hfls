#include "../include/Cluster.hpp"

using namespace std;

Cluster::Cluster(Options options)
{


int symmetric, format, offset;
bool reduceZeroRows, transpose, printCooOrDense;

#if 0
//    std::string path2matrix = options.path2data+"/testMat"+to_string(0)+".txt";
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

    mesh.createMesh();
//
//    for (int i = 0 ; i < mesh.nPoints; i++){
//       cout << mesh.points[i].x << "  " ;
//       cout << mesh.points[i].y << "  " ;
//       cout << mesh.points[i].z << " \n " ;
//    }





    nS = options.n_subdomOnCluster;
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
    Gf.resize(nS);
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



//    create_cluster_constraints();

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

        K_reg[i] = K[i];
        K_reg[i].factorization(nullPivots);



        Matrix BcK_dense;
        BcK_dense = Bc_dense[i];
        BcK_dense.setZero();


        K[i].mult(Bc_dense[i],BcK_dense,true);
        BcK_dense.printToFile("BcK_dense",folder,i,printCooOrDense);
        Lumped[i].zero_dense(Bc[i].n_row_cmprs,  Bc[i].n_row_cmprs);
        Bc[i].mult(BcK_dense,Lumped[i],true);
        Lumped[i].printToFile("Lumped",folder,i,printCooOrDense);

        Matrix BcKplus_dense;
        BcKplus_dense = Bc_dense[i];
        BcKplus_dense.setZero();

        K_reg[i].solve(Bc_dense[i],BcKplus_dense);
        BcKplus_dense.printToFile("BcKplus",folder,i,printCooOrDense);
        Fc[i].zero_dense(Bc[i].n_row_cmprs,  Bc[i].n_row_cmprs);
        Bc[i].mult(BcKplus_dense,Fc[i],true);
        Fc[i].printToFile("Fc",folder,i,printCooOrDense);
        Fc[i].l2g_i_coo = Bc[i].l2g_i_coo;

        /* Gf - constraints matrix */
        int n_rowGf = Bf[i].n_row_cmprs;
        int n_colGf = R[i].n_col;
        Gf[i].zero_dense(n_rowGf, n_colGf );
        Bf[i].mult(R[i],Gf[i],true);

        /* Gc - constraints matrix */
        int n_rowGc = Bc[i].n_row_cmprs;
        int n_colGc = R[i].n_col;
        Gc[i].zero_dense(n_rowGc, n_colGc );
        Bc[i].mult(R[i],Gc[i],true);


        /* print */
        /* K, R, etc.*/
#if VERBOSE_LEVER>3
        cout << "   ... printing matrices start  ...";
#endif
        K[i].printToFile("K",folder,i,printCooOrDense);
        K_reg[i].printToFile("K_reg",folder,i,printCooOrDense);
        R[i].printToFile("R",folder,i,printCooOrDense);
        Bc[i].printToFile("Bc",folder,i,printCooOrDense);
        Bc_dense[i].printToFile("Bc_dense",folder,i,printCooOrDense);
        Bf[i].printToFile("Bf",folder,i,printCooOrDense);
        Gc[i].printToFile("Gc",folder,i,printCooOrDense);
        Gc[i].l2g_i_coo = Bc[i].l2g_i_coo;
        Gf[i].printToFile("Gf",folder,i,printCooOrDense);
        Gf[i].l2g_i_coo = Bf[i].l2g_i_coo;
#if VERBOSE_LEVER>3
        cout << "  " << (i + 1) <<"/" << nS <<"  ...\n";
#endif
    }

//    Matrix::testPardiso();

    create_Gf_clust();
    Gf_clust.printToFile("Gf_clust",folder,0,printCooOrDense);

    /* Fc_clust  */
    create_Fc_clust();
    Fc_clust.printToFile("Fc_clust",folder,0,printCooOrDense);

    create_Gc_clust();
    Gc_clust.printToFile("Gc_clust",folder,0,printCooOrDense);

    create_Ac_clust();
    Ac_clust.getBasicMatrixInfo();

#if VERBOSE_LEVER>3
    cout << "0 Ac_clust.format = " << Ac_clust.format << endl;
#endif
    Ac_clust.printToFile("Ac_clust",folder,0,printCooOrDense);

    for (int i = 0 ; i < nS; i++){
       K_reg[i].FinalizeSolve(i);
    }
}


void Cluster::create_Gf_clust(){

    Matrix & A_clust = Gf_clust;
    vector <Matrix> & A_i = Gf;
    A_clust.symmetric = 0;
    bool remapCols = false;
    create_clust_object(A_clust, A_i, remapCols);
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
    /* updated on nnz */
    A_clust.nnz = cnt;
    tmpVec.resize(A_clust.nnz);
    A_clust.sortAndUniqueCOO(tmpVec);
    /* tmpVec released from the memory */
    tmpVec.clear();
    tmpVec.shrink_to_fit();

    /* both Fc or Gc do not have compressed rows,
     * therefore n_row = n_row_cmprs                        */
    A_clust.n_row = A_clust.n_row_cmprs;

    /* if Fc, n_col == n_row (or n_row_cmprs)               */
    if (A_clust.symmetric > 0)
        A_clust.n_col = A_clust.n_row_cmprs;
    /* transform to CSR format */
    A_clust.COO2CSR();
}

void Cluster::create_Ac_clust(){

    Ac_clust.symmetric = 2;
    Ac_clust.format= 0;


    Ac_clust.nnz = Fc_clust.nnz + Gc_clust.nnz;

    vector < int_int_dbl > tmpVec;
    tmpVec.resize(Ac_clust.nnz);

    int cnt = 0;
    /* Fc_clust */
    for (int i = 0; i < Fc_clust.n_row; i++) {
        for (int j = Fc_clust.i_ptr[i]; j < Fc_clust.i_ptr[i + 1]; j++) {
            tmpVec[cnt].I = i;
            tmpVec[cnt].J = Fc_clust.j_col[j];
            tmpVec[cnt].V = Fc_clust.val[j];
            cnt++;
        }
    }
    for (int i = 0; i < Gc_clust.n_row; i++) {
        for (int j = Gc_clust.i_ptr[i]; j < Gc_clust.i_ptr[i + 1]; j++) {
            tmpVec[cnt].I = i;
            tmpVec[cnt].J = Gc_clust.j_col[j] + Fc_clust.n_col;
            tmpVec[cnt].V = Gc_clust.val[j];
            cnt++;
        }
    }

    Ac_clust.sortAndUniqueCOO(tmpVec);
    tmpVec.clear();
    tmpVec.shrink_to_fit();

    Ac_clust.n_row = Ac_clust.n_row_cmprs;
    Ac_clust.n_col = Fc_clust.n_col + Gc_clust.n_col;

    Ac_clust.COO2CSR();

    /* only upper triangular part of symmetric matrix
     *      Ac = [Fc     Gc]    =   [ Ac[0,0] Ac[0,1] ]
     *           [Gc^T   O ]    =   [ Ac[1,0] Ac[1,1] ]
     * is kept in the memeory.
     * Due to dissection, zero diagonal matrix is placed
     * in the block (1,1)                                   */
    int init_nnz = Ac_clust.nnz;
    int new_nnz = init_nnz + Gc_clust.n_col;
    /* n_row = n_row_cmprs = n_col */
    Ac_clust.nnz = new_nnz;
    Ac_clust.n_row = Ac_clust.n_col;
    Ac_clust.n_row_cmprs = Ac_clust.n_col;
    /* update sizes of sparse structures*/
    Ac_clust.i_ptr.resize(Ac_clust.n_col + 1);
    Ac_clust.j_col.resize(new_nnz);
    Ac_clust.val.resize(new_nnz);
    Ac_clust.l2g_i_coo.resize(Ac_clust.n_col);


    for (int i = 0; i < Gc_clust.n_col; i++){
        Ac_clust.i_ptr[Fc_clust.n_col + 1 + i] = init_nnz + 1 + i;
        Ac_clust.j_col[init_nnz + i] = Fc_clust.n_row + i;
        Ac_clust.val[init_nnz + i] = 0;
        Ac_clust.l2g_i_coo[Fc_clust.n_row + i] = Fc_clust.n_row + i;
    }
    Ac_clust.i_ptr[Ac_clust.n_row] = new_nnz;

}






//  int first[] = {5,10,15,20,25};
//  int second[] = {50,40,30,20,10};
//  std::vector<int> v(10);                      // 0  0  0  0  0  0  0  0  0  0
//  std::vector<int>::iterator it;
//
//  std::sort (first,first+5);     //  5 10 15 20 25
//  std::sort (second,second+5);   // 10 20 30 40 50
//
//  it=std::set_intersection (first, first+5, second, second+5, v.begin());
//                                               // 10 20 0  0  0  0  0  0  0  0
//  v.resize(it-v.begin());                      // 10 20







void Cluster::create_cluster_constraints(){
    vector < vector < int > > subDOFset;


    vector<int>::iterator it;
#if 1

//     25----26----27----28----29
//      |           |           |
//      |           |           |
//     20----21----22----23----24
//      |           |           |
//      |           |           |
//     15    16    17    18    19
//      |           |           |
//      |           |           |
//     10----11----12----13----14
//      |           |           |
//      |           |           |
//      5     6     7     8     9
//      |           |           |
//      |           |           |
//      0-----1-----2-----3-----4

    int nS_ = 6;
    subDOFset.resize(nS_);

    int data0[] = {0, 1, 2, 5, 6, 7, 10, 11, 12};
    subDOFset[0].insert(subDOFset[0].begin(), data0, data0 + 9);

    int data1[] = {2, 3, 4, 7, 8, 9, 12, 13, 14};
    subDOFset[1].insert(subDOFset[1].begin(), data1, data1 + 9);

    int data2[] = {10, 11, 12, 15, 16, 17, 20, 21, 22};
    subDOFset[2].insert(subDOFset[2].begin(), data2, data2 + 9);

    int data3[] = {12, 13, 14, 17, 18, 19, 22, 23, 24};
    subDOFset[3].insert(subDOFset[3].begin(), data3, data3 + 9);

    int data4[] = {20, 21, 22, 25, 26, 27};
    subDOFset[4].insert(subDOFset[4].begin(), data4, data4 + 6);

    int data5[] = {22, 23, 24, 27, 28, 29};
    subDOFset[5].insert(subDOFset[5].begin(), data5, data5 + 6);
#endif


    /* filling subDOFset be DOF set over each subdomain */
    int maxSize=0;
    for (int i = 0 ; i < nS_ ; i++){
//        subDOFset[i] = Bf[i].l2g_i_coo;
        sort(subDOFset[i].begin(), subDOFset[i].end());
        it = unique (subDOFset[i].begin(), subDOFset[i].end());
        subDOFset[i].resize( distance(subDOFset[i].begin(),it));
        if (subDOFset[i].size() > maxSize)
            maxSize = subDOFset[i].size();
    }


    vector<int> v(2 * maxSize);

    for (int i = 0; i < nS_ - 1 ; i++){
        for (int j = i + 1; j < nS_ ; j++){
            cout << "(" << i << ":"  << j << ")\t";
            v.resize( 2 * maxSize );
            it=set_intersection (subDOFset[i].begin(), subDOFset[i].end(),
                                 subDOFset[j].begin(), subDOFset[j].end(), v.begin());
            v.resize(it-v.begin());
            for (int k = 0 ; k < v.size(); k++)
                cout << v[k] << " ";
            cout << endl;
        }
        cout << "\n";
    }



}
