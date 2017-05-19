#include "Cluster.hpp"
#include <unistd.h>
#include <cstdlib>
#ifdef DISSECTION
#include "Dissection.hpp"
#endif
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


/*_NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW_*/
    /* mesh belonging to i-th cluster (currently, i=0 only)*/
    mesh.createMesh(options);
    const int nSubClst = mesh.nSubClst;

    Bc_new.resize(nSubClst);
    Bc_dense_new.resize(nSubClst);
    K_new.resize(nSubClst);
    K_reg_new.resize(nSubClst);
    Fc_new.resize(nSubClst);
    Gc_new.resize(nSubClst);
    R_new.resize(nSubClst);
    Lumped_new.resize(nSubClst);



    cout << "assembling of K, f ... \n" ;
    data.fe_assemb_local_K_f(mesh);
    cout << "symbolic factorization etc. ... \n" ;
    data.feti_symbolic(mesh,K_new);
    data.feti_numeric(mesh,K_new);
    cout << "boolean matrix is being created. ... \n" ;
    create_cluster_constraints();
    cout << "ker(K) is being created ... \n" ;
    data.create_analytic_ker_K(mesh,R_new);



/*_NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW_*/



//    nS = options.n_subdomOnCluster;
    string folder = options.path2data;
//
//    K.resize(nS);
//    K_reg.resize(nS);
//    Lumped.resize(nS);
//    R.resize(nS);
//    rhs.resize(nS);
//    Bc.resize(nS);
//    Bf.resize(nS);
//    Bc_dense.resize(nS);
//    Fc.resize(nS);
//    Gf.resize(nS);
//    Gc.resize(nS);
    bool printMat = false;


    printCooOrDense = true;


/*_NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW_*/
    for (int d = 0; d < K_new.size(); d++){
        cout << "d = " << d << endl;
        if (printMat)
            K_new[d].printToFile("K_new",folder,d,printCooOrDense);
        //K_new[d].getBasicMatrixInfo();
        if (printMat)
            Bc_new[d].printToFile("Bc_new",folder,d,printCooOrDense);
        //Bc_new[d].getBasicMatrixInfo();
//
//
        Bc_dense_new[d] = Bc_new[d];
        reduceZeroRows = true;
        transpose = true;
        Bc_dense_new[d].CSRorCOO2DNS(reduceZeroRows,transpose);
        if (printMat)
            Bc_dense_new[d].printToFile("Bc_dense_new",folder,d,printCooOrDense);
//
//
//
        vector <int > nullPivots_new;
        R_new[d].getNullPivots(nullPivots_new);

//        for (int i = 0 ; i < nullPivots_new.size(); i++){
//            cout << nullPivots_new[i] << " ";
//        }
//        cout << endl;

        K_reg_new[d] = K_new[d];
        K_reg_new[d].factorization(nullPivots_new);
        if (printMat)
            K_reg_new[d].printToFile("K_reg_new",folder,d,printCooOrDense);

//
//
//
        Matrix BcK_dense_new;
        BcK_dense_new = Bc_dense_new[d];
        BcK_dense_new.setZero();
//
//
        K_new[d].mult(Bc_dense_new[d],BcK_dense_new,true);
        if (printMat)
            BcK_dense_new.printToFile("BcK_dense_new",folder,d,printCooOrDense);
        Lumped_new[d].zero_dense(Bc_new[d].n_row_cmprs,  Bc_new[d].n_row_cmprs);
        Bc_new[d].mult(BcK_dense_new,Lumped_new[d],true);
        if (printMat)
            Lumped_new[d].printToFile("Lumped_new",folder,d,printCooOrDense);
//

//      Matrix BcKplus_dense;
//      BcKplus_dense = Bc_dense[i];
//      BcKplus_dense.setZero();

//      K_reg[i].solve(Bc_dense[i],BcKplus_dense);
//      BcKplus_dense.printToFile("BcKplus",folder,i,printCooOrDense);
//      Fc[i].zero_dense(Bc[i].n_row_cmprs,  Bc[i].n_row_cmprs);
//      Bc[i].mult(BcKplus_dense,Fc[i],true);
//      Fc[i].printToFile("Fc",folder,i,printCooOrDense);
//      Fc[i].l2g_i_coo = Bc[i].l2g_i_coo;



        Matrix BcKplus_dense_new;
        BcKplus_dense_new = Bc_dense_new[d];
        BcKplus_dense_new.setZero();



        if (d == 0)
            cout << "Fc[i] (for each subdom) is being created ... \n" ;
        K_reg_new[d].solve(Bc_dense_new[d],BcKplus_dense_new);
        if (printMat)
            BcKplus_dense_new.printToFile("BcKplus_new",folder,d,printCooOrDense);
        Fc_new[d].zero_dense(Bc_new[d].n_row_cmprs,  Bc_new[d].n_row_cmprs);
        Bc_new[d].mult(BcKplus_dense_new,Fc_new[d],true);
        if (printMat)
            Fc_new[d].printToFile("Fc_new",folder,d,printCooOrDense);
        Fc_new[d].l2g_i_coo = Bc_new[d].l2g_i_coo;



        /* Gc - constraints matrix */
        if (d == 0)
            cout << "Gc[i] (for each subdom) is being created ... \n" ;
        int n_rowGc = Bc_new[d].n_row_cmprs;
        int n_colGc = R_new[d].n_col;
        Gc_new[d].zero_dense(n_rowGc, n_colGc );
        Bc_new[d].mult(R_new[d],Gc_new[d],true);
        if (printMat)
            Gc_new[d].printToFile("Gc_new",folder,d,printCooOrDense);
        Gc_new[d].l2g_i_coo = Bc_new[d].l2g_i_coo;



//        Gc[i].printToFile("Gc",folder,i,printCooOrDense);
//        Gc[i].l2g_i_coo = Bc[i].l2g_i_coo;


//        for (int i = 0; i < R_new[d].n_row_cmprs; i++){
//            for (int j = 0; j < R_new[d].n_col; j++){
//                cout << R_new[d].dense[i + R_new[d].n_row_cmprs*j] << " ";
//            }
//            cout << endl;
//        }
        if (printMat)
            R_new[d].printToFile("R_new",folder,d,printCooOrDense);

//
    }



#if 1
    for (int i_sub = 0; i_sub < nSubClst; i_sub++){
        Matrix Y;
        Y.zero_dense(K_new[i_sub].n_row_cmprs , R_new[i_sub].n_col);
        K_new[i_sub].mult(R_new[i_sub],Y,true);

        double normK = K_new[i_sub].norm2();
        double normR = R_new[i_sub].norm2();
        double normKR = Y.norm2();

        printf("|K| = %.6e, ", normK);
        printf("|R| = %.6e, ", normR);
        printf("|KR|/(|K|*|R|) = %.6e \n", normKR / (normK*normR));
    }
#endif



/*_NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW__NEW_*/


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                         READING DATA (TXT)                              */
/*                 matrices in coo format are one-based                    */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
 
//    for (int i = 0 ; i < nS ; i++){
//        /* K - stiffness matrix */
//        string path2matrix = options.path2data+"/K"+to_string(i)+".txt";
//        symmetric = 2; format = 1; offset = 1;
//        K[i].readCooFromFile(path2matrix,symmetric,format,offset);
//
//        /* R - kernel matrix */
//        path2matrix = options.path2data+"/R1"+to_string(i)+".txt";
//        symmetric = 0; format = 2; offset = 1;
//        R[i].readCooFromFile(path2matrix,symmetric,format,offset);
//
//        /* rhs - right-hand-side vector */
//        path2matrix = options.path2data+"/f"+to_string(i)+".txt";
//        rhs[i].readCooFromFile(path2matrix,K[i].n_row);
////
//        /* Bf - constraints matrix */
//        path2matrix = options.path2data+"/B1"+to_string(i)+".txt";
//        symmetric = 0; format = 1; offset = 1;
//        Bf[i].readCooFromFile(path2matrix,symmetric,format,offset);
////
//        /* Bc - constraints matrix */
//        path2matrix = options.path2data+"/B0"+to_string(i)+".txt";
//        symmetric = 0; format = 1; offset = 1;
//        Bc[i].readCooFromFile(path2matrix,symmetric,format,offset);
//
////
//        /* Bc_dense - constraints matrix */
//        Bc_dense[i] = Bc[i];
//        reduceZeroRows = true;
//        transpose = true;
//        Bc_dense[i].CSRorCOO2DNS(reduceZeroRows,transpose);
//
//  }








//#if 0
//    for (int i_sub = 0; i_sub < nS; i_sub++){
//        Matrix Y;
//        Y.zero_dense(K[i_sub].n_row_cmprs , R[i_sub].n_col);
//        K[i_sub].mult(R[i_sub],Y,true);
//
//        double normK = K[i_sub].norm2();
//        double normR = R[i_sub].norm2();
//        double normKR = Y.norm2();
//
//        printf("|K| = %.6e, ", normK);
//        printf("|R| = %.6e, ", normR);
//        printf("|KR|/(|K|*|R|) = %.6e \n", normKR / (normK*normR));
//    }
//#endif
//

//    for (int i = 0; i < nS ; i++ ){
//
//        vector <int > nullPivots;
//        R[i].getNullPivots(nullPivots);
//
//        K_reg[i] = K[i];
//        K_reg[i].factorization(nullPivots);
//
//
//
//        Matrix BcK_dense;
//        BcK_dense = Bc_dense[i];
//        BcK_dense.setZero();
//
//
//        K[i].mult(Bc_dense[i],BcK_dense,true);
//        BcK_dense.printToFile("BcK_dense",folder,i,printCooOrDense);
//        Lumped[i].zero_dense(Bc[i].n_row_cmprs,  Bc[i].n_row_cmprs);
//        Bc[i].mult(BcK_dense,Lumped[i],true);
//        Lumped[i].printToFile("Lumped",folder,i,printCooOrDense);
//
//        Matrix BcKplus_dense;
//        BcKplus_dense = Bc_dense[i];
//        BcKplus_dense.setZero();
//
//        K_reg[i].solve(Bc_dense[i],BcKplus_dense);
//        BcKplus_dense.printToFile("BcKplus",folder,i,printCooOrDense);
//        Fc[i].zero_dense(Bc[i].n_row_cmprs,  Bc[i].n_row_cmprs);
//        Bc[i].mult(BcKplus_dense,Fc[i],true);
//        Fc[i].printToFile("Fc",folder,i,printCooOrDense);
//        Fc[i].l2g_i_coo = Bc[i].l2g_i_coo;
//
//        /* Gf - constraints matrix */
//        int n_rowGf = Bf[i].n_row_cmprs;
//        int n_colGf = R[i].n_col;
//        Gf[i].zero_dense(n_rowGf, n_colGf );
//        Bf[i].mult(R[i],Gf[i],true);
//
//        /* Gc - constraints matrix */
//        int n_rowGc = Bc[i].n_row_cmprs;
//        int n_colGc = R[i].n_col;
//        Gc[i].zero_dense(n_rowGc, n_colGc );
//        Bc[i].mult(R[i],Gc[i],true);
//
//
//        /* print */
//        /* K, R, etc.*/
//#if VERBOSE_LEVER>3
//        cout << "   ... printing matrices start  ...";
//#endif
//        K[i].printToFile("K",folder,i,printCooOrDense);
//        K_reg[i].printToFile("K_reg",folder,i,printCooOrDense);
//        R[i].printToFile("R",folder,i,printCooOrDense);
//        Bc[i].printToFile("Bc",folder,i,printCooOrDense);
//        Bc_dense[i].printToFile("Bc_dense",folder,i,printCooOrDense);
//        Bf[i].printToFile("Bf",folder,i,printCooOrDense);
//        Gc[i].printToFile("Gc",folder,i,printCooOrDense);
//        Gc[i].l2g_i_coo = Bc[i].l2g_i_coo;
//        Gf[i].printToFile("Gf",folder,i,printCooOrDense);
//        Gf[i].l2g_i_coo = Bf[i].l2g_i_coo;
//#if VERBOSE_LEVER>3
//        cout << "  " << (i + 1) <<"/" << nS <<"  ...\n";
//#endif
//    }

//    Matrix::testPardiso();

//    create_Gf_clust_new();
//    Gf_clust_new.printToFile("Gf_clust",folder,0,printCooOrDense);

    /* Fc_clust  */

    cout << "Fc_clust is being created ... \n" ;
    create_Fc_clust_new();
    if (printMat)
        Fc_clust_new.printToFile("Fc_clust_new",folder,0,printCooOrDense);

    cout << "Gc_clust is being created ... \n" ;
    create_Gc_clust_new();
    cout << "====================================================\n";
    if (printMat)
        Gc_clust_new.printToFile("Gc_clust_new",folder,0,printCooOrDense);

    cout << "Ac_clust is being created ... \n" ;
    create_Ac_clust_new();
    Ac_clust_new.getBasicMatrixInfo();

#if VERBOSE_LEVER>3
    cout << "0 Ac_clust_new.format = " << Ac_clust_new.format << endl;
#endif
    Ac_clust_new.printToFile("Ac_clust_new",folder,0,printCooOrDense);

#ifdef DISSECTION
    {
      uint64_t *_dslv = new uint64_t;
      int num_threads = 1;
      diss_init(*_dslv, 0, 1, 0, num_threads, 1);
      int sym = 1; // isSym = 1 + upper_flag = 0
      int decomposer = 0;
      diss_s_fact(*_dslv,
          Ac_clust_new.n_row, &Ac_clust_new.i_ptr[0], &Ac_clust_new.j_col[0],
		  sym, decomposer);
      
      int indefinite_flag = 1;
      int scaling = 2;
      double eps_pivot = 1.0e-2;
      diss_n_fact(*_dslv, &Ac_clust_new.val[0], scaling, eps_pivot,
		indefinite_flag);
      int n0;
      diss_get_kern_dim(*_dslv, &n0);
      fprintf(stderr, "%s %d : ## kernel dimension = %d\n", __FILE__, __LINE__, n0);
    }
#endif
    for (int i = 0 ; i < mesh.nSubClst; i++){
       K_reg_new[i].FinalizeSolve(i);
    }

}


//void Cluster::create_Gf_clust(){
//
//    Matrix & A_clust = Gf_clust;
//    vector <Matrix> & A_i = Gf;
//    A_clust.symmetric = 0;
//    bool remapCols = false;
//    create_clust_object(A_clust, A_i, remapCols);
//}

void Cluster::create_Fc_clust_new(){

    Matrix & A_clust = Fc_clust_new;
    vector <Matrix> & A_i = Fc_new;
    A_clust.symmetric = 2;
    bool remapCols = true;
    create_clust_object(A_clust, A_i, remapCols);
}

void Cluster::create_Gc_clust_new(){

    Matrix & A_clust = Gc_clust_new;
    vector <Matrix> & A_i = Gc_new;
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

void Cluster::create_Ac_clust_new(){

    Ac_clust_new.symmetric = 2;
    Ac_clust_new.format= 0;


    Ac_clust_new.nnz = Fc_clust_new.nnz + Gc_clust_new.nnz;

    vector < int_int_dbl > tmpVec;
    tmpVec.resize(Ac_clust_new.nnz);

    int cnt = 0;
    /* Fc_clust_new */
    for (int i = 0; i < Fc_clust_new.n_row; i++) {
        for (int j = Fc_clust_new.i_ptr[i]; j < Fc_clust_new.i_ptr[i + 1]; j++) {
            tmpVec[cnt].I = i;
            tmpVec[cnt].J = Fc_clust_new.j_col[j];
            tmpVec[cnt].V = Fc_clust_new.val[j];
            cnt++;
        }
    }
    for (int i = 0; i < Gc_clust_new.n_row; i++) {
        for (int j = Gc_clust_new.i_ptr[i]; j < Gc_clust_new.i_ptr[i + 1]; j++) {
            tmpVec[cnt].I = i;
            tmpVec[cnt].J = Gc_clust_new.j_col[j] + Fc_clust_new.n_col;
            tmpVec[cnt].V = Gc_clust_new.val[j];
            cnt++;
        }
    }

    Ac_clust_new.sortAndUniqueCOO(tmpVec);
    tmpVec.clear();
    tmpVec.shrink_to_fit();

    Ac_clust_new.n_row = Ac_clust_new.n_row_cmprs;
    Ac_clust_new.n_col = Fc_clust_new.n_col + Gc_clust_new.n_col;

    Ac_clust_new.COO2CSR();

    /* only upper triangular part of symmetric matrix
     *      Ac = [Fc     Gc]    =   [ Ac[0,0] Ac[0,1] ]
     *           [Gc^T   O ]    =   [ Ac[1,0] Ac[1,1] ]
     * is kept in the memeory.
     * Due to dissection, zero diagonal matrix is placed
     * in the block (1,1)                                   */
    int init_nnz = Ac_clust_new.nnz;
    int new_nnz = init_nnz + Gc_clust_new.n_col;
    /* n_row = n_row_cmprs = n_col */
    Ac_clust_new.nnz = new_nnz;
    Ac_clust_new.n_row = Ac_clust_new.n_col;
    Ac_clust_new.n_row_cmprs = Ac_clust_new.n_col;
    /* update sizes of sparse structures*/
    Ac_clust_new.i_ptr.resize(Ac_clust_new.n_col + 1);
    Ac_clust_new.j_col.resize(new_nnz);
    Ac_clust_new.val.resize(new_nnz);
    Ac_clust_new.l2g_i_coo.resize(Ac_clust_new.n_col);


    for (int i = 0; i < Gc_clust_new.n_col; i++){
        Ac_clust_new.i_ptr[Fc_clust_new.n_col + 1 + i] = init_nnz + 1 + i;
        Ac_clust_new.j_col[init_nnz + i] = Fc_clust_new.n_row + i;
        Ac_clust_new.val[init_nnz + i] = 0;
        Ac_clust_new.l2g_i_coo[Fc_clust_new.n_row + i] = Fc_clust_new.n_row + i;
    }
    Ac_clust_new.i_ptr[Ac_clust_new.n_row] = new_nnz;

}



void Cluster::create_cluster_constraints(){

    int nSubClst = mesh.nSubClst;
    vector < vector < int > > subDOFset;
    vector<int>::iterator it;
    subDOFset.resize(nSubClst);
    data.interface.resize(nSubClst);

    for (int d = 0 ; d <  nSubClst; d++){
        subDOFset[d].insert(subDOFset[d].begin(), data.l2g[d].begin(), data.l2g[d].end());
    }


    /* filling subDOFset be DOF set over each subdomain */
    int maxSize=0;
    for (int i = 0 ; i < nSubClst ; i++){
        sort(subDOFset[i].begin(), subDOFset[i].end());
        it = unique (subDOFset[i].begin(), subDOFset[i].end());
        subDOFset[i].resize( distance(subDOFset[i].begin(),it));
        if (subDOFset[i].size() > maxSize)
            maxSize = subDOFset[i].size();
    }

    vector<int> v(2 * maxSize);


    bool print2cmd = false;

    for (int i = 0; i < nSubClst - 1 ; i++){
        for (int j = i + 1; j < nSubClst ; j++){
            v.resize( 2 * maxSize );
            it=set_intersection (subDOFset[i].begin(), subDOFset[i].end(),
                                 subDOFset[j].begin(), subDOFset[j].end(), v.begin());
            v.resize(it-v.begin());

            if (v.size() > 0){
                for (int k = 0 ; k < v.size(); k++){
                    data.interface[i][v[k]].push_back(j);
                    data.interface[j][v[k]].push_back(i);
                }
            }

            if (print2cmd){
                cout << "(" << i << ":"  << j << ")[" << v.size()<< "]\t";
                for (int k = 0 ; k < v.size(); k++)
                    cout << v[k] << " ";
                cout << endl;
            }
        }
        if (print2cmd){
            cout << "\n";
        }
    }

    int global_DOF;
    int ind_neigh_sub;
    int cntLam = 0;
    int j_col_Bc_curr;
    int j_col_Bc_neigh;
    for (int d = 0; d < nSubClst; d++){
        for ( auto it1 = data.interface[d].begin(); it1 != data.interface[d].end(); ++it1  ){
            global_DOF = it1->first;
            if (print2cmd){
                cout << "["  << d << "]: ";
                cout << global_DOF << " ";
            }
            for (int k = 0; k < it1->second.size(); k++){
                ind_neigh_sub = it1->second[k];
                if (print2cmd){
                    cout << ind_neigh_sub << " ";
                }
                if (ind_neigh_sub > d){
//                if ( it1->second.size() == 1 || (d + 1) == ind_neigh_sub){
                    j_col_Bc_curr = data.g2l[d][global_DOF];
                    Bc_new[d].l2g_i_coo.push_back(cntLam);
                    Bc_new[d].j_col.push_back(j_col_Bc_curr);
                    Bc_new[d].val.push_back(1);

                    j_col_Bc_neigh = data.g2l[ind_neigh_sub][global_DOF];
                    Bc_new[ind_neigh_sub].l2g_i_coo.push_back(cntLam);
                    Bc_new[ind_neigh_sub].j_col.push_back(j_col_Bc_neigh);
                    Bc_new[ind_neigh_sub].val.push_back(-1);

                    //cout << j_col_Bc_curr <<" -------- " << j_col_Bc_neigh << endl;

                    cntLam++;
                    break;
                }
            }
            if (print2cmd){
                cout << endl;
            }
        }
    }

    int _nnz;
    for (int d = 0 ; d < nSubClst; d++){
        _nnz = Bc_new[d].val.size();
        Bc_new[d].nnz = _nnz;
        Bc_new[d].n_row = cntLam;
        Bc_new[d].n_row_cmprs = _nnz;
        Bc_new[d].n_col = data.l2g[d].size();
        Bc_new[d].symmetric = 0;
        Bc_new[d].format = 0;

        Bc_new[d].i_coo_cmpr.resize(Bc_new[d].nnz);
        for (int i = 0; i < Bc_new[d].nnz; i++){
            Bc_new[d].i_coo_cmpr[i] = i;
        }
        Bc_new[d].COO2CSR();

//        for (int i = 0 ; i < Bc_new[d].l2g_i_coo.size();i++)
//            cout << Bc_new[d].l2g_i_coo[i] << " ";
//        cout << endl;
    }
}
