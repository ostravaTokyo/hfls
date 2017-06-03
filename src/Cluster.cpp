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
bool reduceZeroRows, transpose, printCooOrDense, checkOrthogonality;


int printMat = options.print_matrices;

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

    Bc.resize(nSubClst);
    BcT_dense.resize(nSubClst);
    K.resize(nSubClst);
    K_reg.resize(nSubClst);
    Fc.resize(nSubClst);
    Gc.resize(nSubClst);
    R.resize(nSubClst);
    Lumped.resize(nSubClst);

    cout << "assembling of K, f ... \n" ;
    data.fe_assemb_local_K_f(mesh);
    cout << "symbolic factorization etc. ... \n" ;
    data.feti_symbolic(mesh,K);
    data.feti_numeric(mesh,K);
    cout << "ker(K) is being created ... \n" ;

    if (options.solver_opt.solver == 0){
        data.create_analytic_ker_K(mesh,R);
    }

    string folder = options.path2data;
    printCooOrDense = true;
    checkOrthogonality = false;




    for (int d = 0; d < K.size(); d++){
        if (options.solver_opt.solver == 0){
            vector < int > nullPivots;
            R[d].getNullPivots(nullPivots);
            K_reg[d] = K[d];
            K_reg[d].factorization(nullPivots);
            if (printMat > 1){
                K_reg[d].printToFile("K_reg",folder,d,printCooOrDense);
                K_reg[d].getBasicMatrixInfo();
            }
        }
        else if (options.solver_opt.solver == 1){
            K[d].symbolic_factorization();
            R[d].label = "kerK";
            K[d].numeric_factorization(R[d],checkOrthogonality);
        }
    }

    // constraints must wait for factorization (to use R)
    cout << "boolean matrix is being created. ... \n" ;
    create_cluster_constraints(options);



    cout << "subdomain:  ";
    for (int d = 0; d < K.size(); d++){
    cout << d <<" ";
        if (printMat > 1){
            K[d].printToFile("K",folder,d,printCooOrDense);
            K[d].getBasicMatrixInfo();
        }
        if (printMat > 1){
            Bc[d].printToFile("Bc",folder,d,printCooOrDense);
            Bc[d].getBasicMatrixInfo();
        }
//
//
        Matrix Bc_tmp;
        Bc_tmp = Bc[d];

        reduceZeroRows = true;
        transpose = true;
        Bc_tmp.CSRorCOO2DNS(reduceZeroRows,transpose);


        BcT_dense[d].zero_dense(Bc_tmp.n_col,Bc_tmp.n_row_cmprs);
        BcT_dense[d].dense = Bc_tmp.dense;
        BcT_dense[d].label = "BcT_dense";



        if (printMat > 2){
            BcT_dense[d].printToFile("BcT_dense",folder,d,printCooOrDense);
            BcT_dense[d].getBasicMatrixInfo();
        }
//

        // generalized inverse test || K * Kplus * K - K || / || K ||
        bool testGenInv = false;

        if (testGenInv){
            Matrix K_dense,Kplus_K,K_Kplus_K;
            K_dense = K[d]; K_dense.label = "K_dense";



// if dissection is activated, code finishes with CORE DUMP
// with paridios it finishes without error
    cout <<"############################################################################\n";
    fprintf(stderr, "%s %d : code crashes if dissection and following block ('#ifdef 1) are active .\n",
                                                                      __FILE__, __LINE__);
    cout <<"############################################################################\n";
#if 1
            K_dense.CSRorCOO2DNS(false,false);
            K_dense.printToFile("K_dense",folder,d,true);
//
            Kplus_K = K_dense;
            Kplus_K.label = "Kplus_K";
            Kplus_K.setZero();


            if (options.solver_opt.solver == 0){
                K_reg[d].solve(K_dense,Kplus_K);
            }
            else if (options.solver_opt.solver == 1){
                K[d].diss_solve(K_dense,Kplus_K);
            }
//
//
            Kplus_K.printToFile("Kplus_K",folder,d,true);
//
            K[d].mult(Kplus_K,K_Kplus_K,true);
            K_Kplus_K.label = "K_Kplus_K";
//
            K_Kplus_K.printToFile("K_Kplus_K",folder,d,true);

            double norm_Kplus_K = Kplus_K.norm2();
            double norm_K = K[d].norm2();
            double norm_K_Kplus_K = K_Kplus_K.norm2();


            for (int i = 0; i < K_dense.numel; i++){
                K_Kplus_K.dense[i] -= K_dense.dense[i];
            }
            double norm_K_minus_K_Kplus_K = K_Kplus_K.norm2();


            cout << " =============================================  \n";
            printf("       || K ||                   = %3.15e\n", norm_K);
            printf("   || Kplus *_K ||               = %3.15e\n", norm_Kplus_K);
            printf("  || K * Kplus * K ||            = %3.15e\n", norm_K_Kplus_K);
            printf("|| K - K * Kplus * K || / || K|| = %3.15e\n", norm_K_minus_K_Kplus_K / norm_K);
            cout << " =============================================  \n";
#endif
        }





        Matrix KplusBcT_dense;
        KplusBcT_dense.zero_dense(BcT_dense[d].n_row_cmprs, BcT_dense[d].n_col);

        if (d == 0)
            cout << "Fc[i] (for each subdom) is being created ... \n" ;

        if (options.solver_opt.solver == 0){
            K_reg[d].solve(BcT_dense[d],KplusBcT_dense);
        }
        else if (options.solver_opt.solver == 1 ){
            K[d].diss_solve(BcT_dense[d],KplusBcT_dense);
        }
        if (printMat > 2){
            KplusBcT_dense.printToFile("KplusBcT",folder,d,printCooOrDense);
            KplusBcT_dense.getBasicMatrixInfo();
        }


        Fc[d].zero_dense(Bc[d].n_row_cmprs,  Bc[d].n_row_cmprs);
        Bc[d].mult(KplusBcT_dense,Fc[d],true);
        if (printMat > 1){
            Fc[d].printToFile("Fc",folder,d,printCooOrDense);
            Fc[d].getBasicMatrixInfo();
        }
        Fc[d].l2g_i_coo = Bc[d].l2g_i_coo;

        /* Gc - constraints matrix */
        if (d == 0)
            cout << "Gc[i] (for each subdom) is being created ... \n" ;
        int n_rowGc = Bc[d].n_row_cmprs;
        int n_colGc = R[d].n_col;
        Gc[d].zero_dense(n_rowGc, n_colGc );
        Bc[d].mult(R[d],Gc[d],true);
        if (printMat > 1){
            Gc[d].printToFile("Gc",folder,d,printCooOrDense);
            Gc[d].getBasicMatrixInfo();
        }

        Gc[d].l2g_i_coo = Bc[d].l2g_i_coo;
        Gc[d].g2l_i_coo = Bc[d].g2l_i_coo;

        if (printMat > 1){
            R[d].printToFile("R",folder,d,printCooOrDense);
            R[d].getBasicMatrixInfo();
        }
    }
    cout << endl;

    /* Fc_clust  */
    cout << "Fc_clust is being created ... \n" ;
    create_Fc_clust();
    if (printMat > 0)
        Fc_clust.printToFile("Fc_clust",folder,0,printCooOrDense);

    Matrix S_Fc_clust;
    Matrix::getEigVal_DNS(Fc_clust,S_Fc_clust,8,5);

    cout << "Gc_clust is being created ... \n" ;
    create_Gc_clust();
    cout << "====================================================\n";
    if (printMat > 0)
        Gc_clust.printToFile("Gc_clust",folder,0,printCooOrDense);

    cout << "GcTGc is being created ... \n" ;
    if (options.solver_opt.GcTGc_assembl_block_by_block){
        cout << " ... block version "<< endl;
        create_GcTGc_clust_sparse();
        GcTGc_clust = GcTGc_sparse_clust;
    }
    else{
        create_GcTGc();
        GcTGc_clust.getBasicMatrixInfo();
    }

    Matrix S_GcTGc;
    Matrix::getEigVal_DNS(GcTGc_clust,S_GcTGc,8,4);
    Matrix::getSingularVal_DNS(GcTGc_clust,S_GcTGc,4,8);


    GcTGc_sparse_clust.getBasicMatrixInfo();
    if (printMat > 1)
        GcTGc_clust.printToFile("GcTGc_clust",folder,0,printCooOrDense);


//    Matrix ker_GcT;
    kerGc.label = "kerGc";
    GcTGc_clust.order_number = 1000;
    GcTGc_clust.symbolic_factorization();
    checkOrthogonality = true;
    GcTGc_clust.numeric_factorization(kerGc,checkOrthogonality);
    fprintf(stderr, "%s %d : ## GcTGc_clust: kernel dimension = %d\n",
                                                    __FILE__, __LINE__, kerGc.n_col);
    if (printMat > 1)
        kerGc.printToFile("kerGc",folder,0,printCooOrDense);

    cout << "Ac_clust is being created ... \n" ;
    create_Ac_clust(options.solver_opt.Ac_extended_by_kerGc);
    Matrix S_Ac_clust;
    Matrix::getEigVal_DNS(Ac_clust,S_Ac_clust,10);

    Ac_clust.getBasicMatrixInfo();





    if (printMat > 1)
        Ac_clust.printToFile("Ac_clust",folder,0,printCooOrDense);

    Ac_clust.order_number = 2000;
    Ac_clust.symbolic_factorization();
    Ac_clust.diss_scaling = 1;

    if (options.solver_opt.Ac_extended_by_kerGc){
        if (options.solver_opt.solver == 0){
            Ac_clust.msglvl = 0;
            Ac_clust.factorization();
        }
        else if (options.solver_opt.solver == 0){
            Ac_clust.diss_scaling = 1;
            Ac_clust.numeric_factorization();
        }
    }
    else{
        Matrix ker_Ac;
        ker_Ac.label = "ker_Ac";
        Ac_clust.numeric_factorization(ker_Ac,checkOrthogonality);
        fprintf(stderr, "%s %d : ## Ac_clust: kernel dimension = %d\n",
                                                    __FILE__, __LINE__, ker_Ac.n_col);
    }


    if (options.solver_opt.solver == 0){
        for (int i = 0 ; i < mesh.nSubClst; i++){
           K_reg[i].FinalizeSolve(i);
        }
    }
    else if (options.solver_opt.solver == 1){
        for (int i = 0 ; i < mesh.nSubClst; i++){
           K[i].FinalizeSolve(i);
        }
    }
}


void Cluster::create_GcTGc(){
    Matrix &GcTGc = GcTGc_clust;
    vector<int> GcTGcmask;
    vector<double> GcTGcdense;
    const int nncol = Gc_clust.n_col;
    GcTGcmask.resize(nncol * nncol, 0);
    GcTGcdense.resize(nncol * nncol, 0.0);
    vector<int> &i_ptr = Gc_clust.i_ptr;
    vector<int> &j_col = Gc_clust.j_col;
    for (int k = 0; k < Gc_clust.n_row; k++) {
      for (int l = i_ptr[k]; l < i_ptr[k + 1]; l++) {
    for (int m = l; m < i_ptr[k + 1]; m++) {
      int itmp = j_col[l] + j_col[m] * nncol;
      GcTGcmask[itmp] = 1;
      GcTGcdense[itmp] += Gc_clust.val[l] * Gc_clust.val[m];
    }
      }
    }
    int cnt = 0;
    for (int i = 0; i < nncol * nncol; i++) {
      if (GcTGcmask[i] == 1) {
    cnt++;
      }
    }
    fprintf(stderr, "%s %d : %d / %d = %d * %d\n", __FILE__, __LINE__,
        cnt, (nncol * nncol), nncol, nncol);
    GcTGc.symmetric = 2;
    GcTGc.format= 0;
    GcTGc.nnz = cnt;

    vector < int_int_dbl > tmpVec;

    tmpVec.resize(cnt);
    cnt = 0;
    for (int i = 0; i < nncol ; i++) {
      for (int j = 0; j < nncol ; j++) {
    if (GcTGcmask[i + j * nncol] == 1) {
      tmpVec[cnt].I = i;
      tmpVec[cnt].J = j;
      tmpVec[cnt].V = GcTGcdense[i + j * nncol];
      cnt++;
    }
      }
    }
    GcTGc.sortAndUniqueCOO(tmpVec);
    GcTGc.n_row = nncol;
    GcTGc.n_row_cmprs = nncol;
    GcTGc.n_col = nncol;
    GcTGc.COO2CSR();
}



void Cluster::create_GcTGc_clust_sparse(){
//
    vector<int>::iterator it;
    vector<int> overlap(n_interf_c_max);
    int ind_neigh;
    GcTGc_sparse_clust.label = "GcTGc_sparse";
    GcTGc_sparse_clust.nnz = 0;


    GcTGc_sparse_clust.symmetric = 2;
    GcTGc_sparse_clust.format= 0;


    vector < int_int_dbl > tmpVec;
    vector < int > blockPointer;
    blockPointer.resize(mesh.nSubClst);
    blockPointer[0] = 0;
    for (int i = 1; i < mesh.nSubClst; i ++)
        blockPointer[i] = blockPointer[i-1] + Gc[i-1].n_col;

    for (int i = 0 ; i < mesh.nSubClst ; i++){
        GcTGc_sparse_clust.nnz += (Gc[i].n_col + 1) * Gc[i].n_col;
        Matrix Gii;
        Matrix Gc_i = Gc[i];
        Gc_i.CSRorCOO2DNS(true,false);
        Gii.mat_mult_dense(Gc_i,"T",Gc_i,"N");
        Gii.label = "Gii";
        Matrix::updateCOOstructure( tmpVec,
                                Gii,blockPointer[i],blockPointer[i]);
        for (int j = 0 ; j < neighbours[i].size(); j++){
            ind_neigh = neighbours[i][j];
            overlap.resize( n_interf_c_max );
            it=set_intersection (Gc[i].l2g_i_coo.begin(), Gc[i].l2g_i_coo.end(),
                                 Gc[ind_neigh].l2g_i_coo.begin(), Gc[ind_neigh].l2g_i_coo.end(),
                                 overlap.begin());
            overlap.resize(it-overlap.begin());
            int n_overlap = overlap.size();
            if (n_overlap > 0){
                Matrix G_i, G_j, GiTGj;
                G_i.submatrix_row_selector(Gc[i],overlap);
                G_j.submatrix_row_selector(Gc[ind_neigh],overlap);
                GiTGj.mat_mult_dense(G_i,"T",G_j,"N");
                Matrix::updateCOOstructure( tmpVec,
                                            GiTGj,
                                            blockPointer[i],
                                            blockPointer[ind_neigh]);
            }
            GcTGc_sparse_clust.nnz +=
                Gc[i].n_col * Gc[ind_neigh].n_col;
        }
    }
    GcTGc_sparse_clust.sortAndUniqueCOO(tmpVec);
    GcTGc_sparse_clust.n_col = GcTGc_sparse_clust.n_row;
    GcTGc_sparse_clust.COO2CSR();
}

void Cluster::create_Fc_clust(){

    Fc_clust.label = "Fc_clust";
    Matrix & A_clust = Fc_clust;
    vector <Matrix> & A_i = Fc;
    A_clust.symmetric = 2;
    bool remapCols = true;
    create_clust_object(A_clust, A_i, remapCols);
}

void Cluster::create_Gc_clust(){

    Gc_clust.label = "Gc_clust";
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

void Cluster::create_Ac_clust(bool Ac_nonsingular){


    /* only upper triangular part of symmetric matrix
     *      Ac = [Fc     Gc]    =   [ Ac[0,0] Ac[0,1] ]
     *           [Gc^T   O ]    =   [ Ac[1,0] Ac[1,1] ]
     * is kept in the memeory.
     * Due to dissection, zero diagonal matrix is placed
     * in the block (1,1)

     *  if Ac_nonsingular = true then
     *      Ac = [Fc     Gc   0]    =   [ Ac[0,0] Ac[0,1] Ac[0,2] ]
     *           [Gc^T   O    H]    =   [ Ac[1,0] Ac[1,1] Ac[1,2] ]
     *           [ O     H^T  0]    =   [ Ac[2,0] Ac[2,1] Ac[2,2] ]
     *   where H = ker(GcTGc)
     */

    Ac_clust.label = "Ac_clust";
    Ac_clust.symmetric = 2;
    Ac_clust.format = 0;

//    bool Ac_singular = false;


    Ac_clust.nnz = Fc_clust.nnz + Gc_clust.nnz + Gc_clust.n_col;
    if (Ac_nonsingular)
        Ac_clust.nnz += kerGc.numel + kerGc.n_col;

//    if (Ac_nonsingular){
//        Ac_clust.nnz = Fc_clust.nnz + Gc_clust.nnz + Gc_clust.n_col +
//                                                            kerGc.numel + kerGc.n_col;
//    }
//    else
//    {
//        Ac_clust.nnz = Fc_clust.nnz + Gc_clust.nnz + Gc_clust.n_col;
//    }

    vector < int_int_dbl > tmpVec;
    tmpVec.resize(Ac_clust.nnz);

    int cnt = 0;


    //  Ac[0,1]
    for (int i = 0; i < Fc_clust.n_row; i++) {
        for (int j = Fc_clust.i_ptr[i]; j < Fc_clust.i_ptr[i + 1]; j++) {
            tmpVec[cnt].I = i;
            tmpVec[cnt].J = Fc_clust.j_col[j];
            tmpVec[cnt].V = Fc_clust.val[j];
            cnt++;
        }
    }
    //  Ac[0,1]
    for (int i = 0; i < Gc_clust.n_row; i++) {
        for (int j = Gc_clust.i_ptr[i]; j < Gc_clust.i_ptr[i + 1]; j++) {
            tmpVec[cnt].I = i;
            tmpVec[cnt].J = Gc_clust.j_col[j] + Fc_clust.n_col;
            tmpVec[cnt].V = Gc_clust.val[j];
            cnt++;
        }
    }
    //  Ac[1,1]
    for (int i = 0; i < Gc_clust.n_col; i++){
        tmpVec[cnt].I = i + Fc_clust.n_row;
        tmpVec[cnt].J = i + Fc_clust.n_row;
        tmpVec[cnt].V = 0.0;
        cnt++;
    }

    if (Ac_nonsingular){
    //  Ac[1,2]
        for (int i = 0; i < kerGc.n_row; i++) {
            for (int j = 0; j < kerGc.n_col; j++) {
                tmpVec[cnt].I = i + Fc_clust.n_row;
                tmpVec[cnt].J = j + Fc_clust.n_col + Gc_clust.n_col;
                tmpVec[cnt].V = kerGc.dense[i + j * kerGc.n_row];
                cnt++;
            }
        }
    //  Ac[2,2]
        for (int i = 0; i < kerGc.n_col; i++){
            tmpVec[cnt].I = i + Fc_clust.n_row + Gc_clust.n_col;
            tmpVec[cnt].J = i + Fc_clust.n_row + Gc_clust.n_col;
            tmpVec[cnt].V = 0.0;
            cnt++;
        }
    }
    //else {
    //}

    Ac_clust.sortAndUniqueCOO(tmpVec);
    tmpVec.clear();
    tmpVec.shrink_to_fit();

    if (Ac_nonsingular){
        Ac_clust.n_row = Fc_clust.n_row_cmprs + Gc_clust.n_col + kerGc.n_col;
    }
    else{
        Ac_clust.n_row = Fc_clust.n_row_cmprs + Gc_clust.n_col;
    }

    Ac_clust.n_row_cmprs = Ac_clust.n_row;
    Ac_clust.n_col = Ac_clust.n_row;

    Ac_clust.COO2CSR();
#if 0
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
#endif
}



void Cluster::create_cluster_constraints(const Options &options){

    int nSubClst = mesh.nSubClst;
    vector < vector < int > > subDOFset;
    vector<int>::iterator it;
    subDOFset.resize(nSubClst);
    data.interface.resize(nSubClst);
    int cntLam = 0;
    int j_col_Bc_curr;
    int j_col_Bc_neigh;


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
    neighbours.resize(nSubClst);
    bool print2cmd = false;

    for (int i = 0; i < nSubClst - 1 ; i++){
        for (int j = i + 1; j < nSubClst ; j++){
            v.resize( 2 * maxSize );
            it=set_intersection (subDOFset[i].begin(), subDOFset[i].end(),
                                 subDOFset[j].begin(), subDOFset[j].end(), v.begin());
            v.resize(it-v.begin());

            if (v.size() > 0){

                // TODO !!! temporarly
                Matrix Bc_from_Rt;
                if (options.solver_opt.typeBc == 1){
                   Bc_from_Rt.zero_dense(R[i].n_col,v.size());
                }


                neighbours[i].push_back(j);
                for (int k = 0 ; k < v.size(); k++){
                    data.interface[i][v[k]].push_back(j);
                    data.interface[j][v[k]].push_back(i);

                    if (options.solver_opt.typeBc == 1){
                        for (int l = 0; l < R[i].n_col;l++){
                            int g2ldof = data.g2l[i][v[k]];
                            Bc_from_Rt.dense[l + k * Bc_from_Rt.n_row_cmprs] =
                                    R[i].dense[g2ldof + l * R[i].n_row_cmprs];
                        }
                    }
                }
                if (options.solver_opt.typeBc == 1){

                    bool is_full_column_rank =
                            Matrix::test_of_Bc_constraints(Bc_from_Rt);

                    if (is_full_column_rank){
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        for (int k = 0; k < Bc_from_Rt.n_row_cmprs;k++){
                            for (int l = 0; l < Bc_from_Rt.n_col;l++){

                                double Bc_from_Rt_lk = Bc_from_Rt.dense[k + l * Bc_from_Rt.n_row_cmprs];

                                // @! i_coo_cmpr is firstly filled by cluster global numbering
                                //    and later remaped to local subdomain numbering
                                j_col_Bc_curr = data.g2l[i][v[l]];
                                Bc[i].i_coo_cmpr.push_back(cntLam);
                                Bc[i].j_col.push_back(j_col_Bc_curr);
                                Bc[i].val.push_back(Bc_from_Rt_lk);

                                j_col_Bc_neigh = data.g2l[j][v[l]];
                                Bc[j].i_coo_cmpr.push_back(cntLam);
                                Bc[j].j_col.push_back(j_col_Bc_neigh);
                                Bc[j].val.push_back(-Bc_from_Rt_lk);
                            }
                            cntLam++;
                        }
                    }
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


    if (options.solver_opt.typeBc != 1 ){
        cntLam = 0;
        for (int d = 0; d < nSubClst; d++){
            for ( auto it1 = data.interface[d].begin(); it1 != data.interface[d].end(); ++it1  ){
                global_DOF = it1->first;
                if (print2cmd){
                    cout << "["  << d << "]: ";
                    cout << global_DOF << " ";
                }

                if (options.solver_opt.typeBc == 0){
                    it = find (mesh.cornerDOFs.begin(), mesh.cornerDOFs.end(), global_DOF);
                    if (it == mesh.cornerDOFs.end())
                        continue;
                }
                //TODO: if 'global_DOF' is found in vector 'cornerDOFs', the entry
                //                                                      should be deleted.

                for (int k = 0; k < it1->second.size(); k++){
                    ind_neigh_sub = it1->second[k];
                    if (print2cmd){
                        cout << ind_neigh_sub << " ";
                    }
                    if (ind_neigh_sub > d){
    //                if ( it1->second.size() == 1 || (d + 1) == ind_neigh_sub){
                        j_col_Bc_curr = data.g2l[d][global_DOF];
                        Bc[d].l2g_i_coo.push_back(cntLam);
                        Bc[d].j_col.push_back(j_col_Bc_curr);
                        Bc[d].val.push_back(1);


                        j_col_Bc_neigh = data.g2l[ind_neigh_sub][global_DOF];
                        Bc[ind_neigh_sub].l2g_i_coo.push_back(cntLam);
                        Bc[ind_neigh_sub].j_col.push_back(j_col_Bc_neigh);
                        Bc[ind_neigh_sub].val.push_back(-1);

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
    }

    int _nnz;
    n_interf_c_max = 0;
    for (int d = 0 ; d < nSubClst; d++){
        _nnz = Bc[d].val.size();
        Bc[d].nnz = _nnz;
        Bc[d].n_row = cntLam;
        Bc[d].n_col = data.l2g[d].size();
        Bc[d].symmetric = 0;
        Bc[d].format = 0;


        if (options.solver_opt.typeBc == 1){
            vector<int>::iterator it;
            vector < int > l2g_(Bc[d].i_coo_cmpr);
            sort(l2g_.begin(), l2g_.end());
            it = unique (l2g_.begin(), l2g_.end());
            l2g_.resize( distance(l2g_.begin(),it));
            Bc[d].n_row_cmprs = l2g_.size();


            for (int i = 0; i < l2g_.size(); i++){
                Bc[d].g2l_i_coo.insert ( pair < int, int > ( l2g_[i] , i ));
            }

            for (int i = 0; i < Bc[d].i_coo_cmpr.size(); i++){
                Bc[d].i_coo_cmpr[i] =
                        Bc[d].g2l_i_coo[ Bc[d].i_coo_cmpr[i] ];
            }
            Bc[d].l2g_i_coo = l2g_;

        }
        else{

            Bc[d].n_row_cmprs = Bc[d].l2g_i_coo.size();
            Bc[d].i_coo_cmpr.resize(Bc[d].n_row_cmprs);
            for (int i = 0; i < Bc[d].n_row_cmprs; i++){
                Bc[d].i_coo_cmpr[i] = i;
                Bc[d].g2l_i_coo.insert ( pair < int, int > ( Bc[d].l2g_i_coo[i] , i ));
            }
        }




        Bc[d].COO2CSR();
//        Bc[d].getBasicMatrixInfo();

        // max number of 'c-type' constraints on one subdomain
        if (Bc[d].n_row_cmprs > n_interf_c_max)
            n_interf_c_max = Bc[d].n_row_cmprs;

    }
}
