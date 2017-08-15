#include "Cluster.hpp"
#include <unistd.h>
#include <cstdlib>
#ifdef DISSECTION
#include "Dissection.hpp"
#endif
using namespace std;

Cluster::Cluster(Options options_)
{

    options = options_;

    int symmetric, format, offset;
    bool reduceZeroRows, transpose, printCooOrDense, checkOrthogonality;
    int printMat = options.print_matrices;
    folder = options.path2data;

    /* mesh belonging to i-th cluster (currently, i=0 only)*/
    mesh.createMesh(options);
    nSubClst = mesh.nSubClst;

    Bc.resize(nSubClst);
    KplusBcT.resize(nSubClst);
    Bf.resize(nSubClst);
    BcT_dense.resize(nSubClst);
    K.resize(nSubClst);

    Preconditioner.resize(nSubClst);

    rhs.resize(nSubClst);
    Fc.resize(nSubClst);
    Gc.resize(nSubClst);
    R.resize(nSubClst);
    Rf.resize(nSubClst);
    Gf.resize(nSubClst);


    cout << "assembling of K, f ... \n" ;
    data.fe_assemb_local_K_f(mesh);
    cout << "symbolic factorization etc. ... \n" ;
    data.feti_symbolic(mesh,K);
    data.feti_numeric(mesh,K,rhs);
    cout << "ker(K) is being created ... \n" ;

    if (options.solver_opt.create_analytic_ker_K){
        data.create_analytic_ker_K(mesh,R);
    }

    printCooOrDense = true;
    checkOrthogonality = true;

    neqClst = 0;
    for (int d = 0; d < K.size(); d++){
        neqClst += K[d].n_row_cmprs;
    }

//    Matrix::testSolver(folder, 1);

    for (int d = 0; d < K.size(); d++){
        K[d].order_number = d;
        K[d].options = options;
        K[d].sym_factor(options.solver_opt.solver);
        R[d].label = "kerK";
//        Matrix Ksing = K[d];
        K[d].num_factor(R[d],checkOrthogonality);
//        K[d].test_K_Kp_K_condition(Ksing);
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

            rhs[d].printToFile("rhs",folder,d,printCooOrDense);
            rhs[d].label = "rhs";
            rhs[d].getBasicMatrixInfo();
            Bc[d].printToFile("Bc",folder,d,printCooOrDense);
            Bc[d].getBasicMatrixInfo();
            Bf[d].printToFile("Bf",folder,d,printCooOrDense);
            Bf[d].getBasicMatrixInfo();
        }
//
        Matrix Bc_tmp;
        Bc_tmp = Bc[d];


        Preconditioner[d].createDirichletPreconditioner(Bf[d],K[d]);

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
        Matrix KplusBcT_dense;
        KplusBcT_dense.zero_dense(BcT_dense[d].n_row_cmprs, BcT_dense[d].n_col);

        if (d == 0)
            cout << "Fc[i] (for each subdom) is being created ... \n" ;

        K[d].solve_system(BcT_dense[d],KplusBcT_dense);
        KplusBcT[d] = KplusBcT_dense;
        KplusBcT[d].symmetric = 0;

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


        //
        for (int i = 0; i < Gc[d].numel; i++)
            Gc[d].dense[i] *= -1;

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
    Fc_clust.getBasicMatrixInfo();
    if (printMat > 0)
        Fc_clust.printToFile("Fc_clust",folder,0,printCooOrDense);

    Matrix S_Fc_clust;
    Matrix::getEigVal_DNS(Fc_clust,S_Fc_clust,15,3);


    Fc_clust.sym_factor(options.solver_opt.solver);
//    Matrix Fc_copy = Fc_clust;
    Fc_clust.num_factor();
//    Fc_clust.test_K_Kp_K_condition(Fc_copy);



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
    nRBM_c = GcTGc_clust.n_row_cmprs;

    Matrix S_GcTGc;
    Matrix::getEigVal_DNS(GcTGc_clust,S_GcTGc,10,3);
//    Matrix::getSingularVal_DNS(GcTGc_clust,S_GcTGc,3,10);


    if (printMat > 0)
        GcTGc_clust.printToFile("GcTGc_clust",folder,0,printCooOrDense);


//    Matrix ker_GcT;
    kerGc.label = "kerGc";
    GcTGc_clust.order_number = 1000;
    GcTGc_clust.sym_factor(1);
    checkOrthogonality = true;
    GcTGc_clust.num_factor(kerGc,checkOrthogonality );
    nRBM_f = kerGc.n_col;
    fprintf(stderr, "%s %d : ## GcTGc_clust: kernel dimension = %d\n",
                                                    __FILE__, __LINE__, kerGc.n_col);

    GcTGc_clust.getBasicMatrixInfo();

    if (printMat > 0)
        kerGc.printToFile("kerGc",folder,0,printCooOrDense);

    cout << "Ac_clust is being created ... \n" ;
    create_Ac_clust(options.solver_opt.Ac_extended_by_kerGc);
    Matrix S_Ac_clust;
    Matrix::getEigVal_DNS(Ac_clust,S_Ac_clust,10,3);


    if (printMat > 0)
        Ac_clust.printToFile("Ac_clust",folder,0,printCooOrDense);

    Ac_clust.order_number = 2000;

    Ac_clust.diss_scaling = 1;
    Ac_clust.sym_factor(options.solver_opt.solver);



//    Matrix Ac_clust_copy;
//    Ac_clust_copy = Ac_clust;
    if (options.solver_opt.Ac_extended_by_kerGc){
        Ac_clust.diss_scaling = 1;

        if (options.solver_opt.solver == 0){
            Ac_clust.iparm[9] = 8;
        }
        Ac_clust.num_factor();
    }
    else{
        Matrix ker_Ac;
        ker_Ac.label = "ker_Ac";
        checkOrthogonality = true;
        Ac_clust.num_factor(ker_Ac,checkOrthogonality);
        //Ac_clust.num_factor();
        fprintf(stderr, "%s %d : ## Ac_clust: kernel dimension = %d\n",
                                                    __FILE__, __LINE__, ker_Ac.n_col);

        if (printMat > 0)
            ker_Ac.printToFile("ker_Ac",folder,0,printCooOrDense);
    }


//    Ac_clust.test_K_Kp_K_condition(Ac_clust_copy);



      //kerGc;
      create_Rf_and_Gf();
      create_GfTGf();
      if (printMat>2){
          Gf_clust.printToFile("Gf_clust",folder,0,true);
          GfTGf.printToFile("GfTGf",folder,0,true);
          invGfTGf.printToFile("iGfTGf",folder,0,true);
          for (int d = 0; d < nSubClst; d++)
            Rf[d].printToFile("Rf",folder,d,true);


      }
      compute_invGfTGf();

      // conjugate gradient //
      pcpg2();




/* ################################################################################### */
/* ########################### FINALIZING ############################################ */
/* ################################################################################### */

    for (int i = 0 ; i < nSubClst; i++)
        K[i].FinalizeSolve(i);
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

    // max number of 'c-type' constraints on one subdomain
    int n_interf_c_max = 0;

    for (int d = 0; d < nSubClst;d++){
        if (Bc[d].n_row_cmprs > n_interf_c_max)
            n_interf_c_max = Bc[d].n_row_cmprs;
    }

    vector<int>::iterator it;
    vector<int> overlap(n_interf_c_max);
    int ind_neigh;
    GcTGc_sparse_clust.label = "GcTGc_sparse";
    GcTGc_sparse_clust.nnz = 0;
    GcTGc_sparse_clust.symmetric = 2;
    GcTGc_sparse_clust.format= 0;
    vector < int_int_dbl > tmpVec;
    vector < int > blockPointer;
    blockPointer.resize(nSubClst);
    blockPointer[0] = 0;
    for (int i = 1; i < nSubClst; i ++)
        blockPointer[i] = blockPointer[i-1] + Gc[i-1].n_col;

    for (int i = 0 ; i < nSubClst ; i++){
        GcTGc_sparse_clust.nnz += (Gc[i].n_col + 1) * Gc[i].n_col;
        Matrix Gii("Gii");
        Matrix Gc_i = Gc[i];
        Gc_i.CSRorCOO2DNS(true,false);
        Gii.mat_mult_dense(Gc_i,"T",Gc_i,"N");
	//        Gii.label = "Gii";
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

    vector < vector < int > > subDOFset;
    vector<int>::iterator it;
    subDOFset.resize(nSubClst);
    data.interface.resize(nSubClst);
    data.interfaces.resize(nSubClst);

    /* filling subDOFset be DOF set over each subdomain */
    for (int d = 0 ; d <  nSubClst; d++){
        subDOFset[d].insert(subDOFset[d].begin(), data.l2g[d].begin(), data.l2g[d].end());
    }
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

    for (int i = 0; i < nSubClst - 1 ; i++){
        for (int j = i + 1; j < nSubClst ; j++){
            v.resize( 2 * maxSize );
            it=set_intersection (subDOFset[i].begin(), subDOFset[i].end(),
                                 subDOFset[j].begin(), subDOFset[j].end(), v.begin());
            v.resize(it-v.begin());

            if (v.size() > 0){
                Interfaces interfaces_;
                data.interfaces[i].push_back(interfaces_);
                data.interfaces[i].back().IdNeighSub = j;
                data.interfaces[i].back().dofs.resize(v.size());
                neighbours[i].push_back(j);
                for (int k = 0 ; k < v.size(); k++){
                    data.interface[i][v[k]].push_back(j);
                    data.interface[j][v[k]].push_back(i);
                    data.interfaces[i].back().dofs[k] = v[k];
                }
            }
        }
    }


    if (options.solver_opt.typeBc == 1 ){
        create_Bc_weightedAverages_in_COO(Bc,options.solver_opt.Bc_fullRank);
    }
    else{
        bool cornersOnlyOrAllDof;
        if (options.solver_opt.typeBc == 0)
            cornersOnlyOrAllDof = true;
        else
            cornersOnlyOrAllDof = false;

        create_Bc_or_Bf_in_CSR(Bc,options.solver_opt.Bc_fullRank, cornersOnlyOrAllDof);
    }


    bool Bf_full_rank           = false;
    bool Bf_cornersOnlyOrAllDof = false;
    bool Bf_addDirConstr        = true;
    create_Bc_or_Bf_in_CSR(Bf,Bf_full_rank,Bf_cornersOnlyOrAllDof,Bf_addDirConstr);


    nLam_c = Bc[0].n_row;
    nLam_f = Bf[0].n_row;

}


void Cluster::create_Bc_weightedAverages_in_COO(vector <Matrix> &Bc_, bool Bc_fullRank){


    int cntLam = 0;
    int j_col_Bc_curr;
    int j_col_Bc_neigh;

    for (int i = 0; i < data.interfaces.size(); i++){
        for (int j = 0;j < data.interfaces[i].size();j++){
            int IdNeighSub = data.interfaces[i][j].IdNeighSub;
            int n_com_dof = data.interfaces[i][j].dofs.size();
            Matrix Bc_from_Rt;
            Bc_from_Rt.zero_dense(R[i].n_col,n_com_dof);

            for (int k = 0 ; k < n_com_dof; k++){
                int dofs_k = data.interfaces[i][j].dofs[k];
                for (int l = 0; l < R[i].n_col;l++){
                    int g2ldof = data.g2l[i][dofs_k];
                    Bc_from_Rt.dense[l + k * Bc_from_Rt.n_row_cmprs] =
                            R[i].dense[g2ldof + l * R[i].n_row_cmprs];
                }
            }

            bool is_local_Bc_full_column_rank =
                    Matrix::test_of_Bc_constraints(Bc_from_Rt);

            if (!Bc_fullRank || is_local_Bc_full_column_rank){
                for (int k = 0; k < Bc_from_Rt.n_row_cmprs;k++){
                    for (int l = 0; l < Bc_from_Rt.n_col;l++){
                        int dofs_l = data.interfaces[i][j].dofs[l];
                        double Bc_from_Rt_lk = Bc_from_Rt.dense[k + l * Bc_from_Rt.n_row_cmprs];
                        // @! i_coo_cmpr is firstly filled by cluster global numbering
                        //    and later remaped to local subdomain numbering
                        j_col_Bc_curr = data.g2l[i][dofs_l];
                        Bc_[i].i_coo_cmpr.push_back(cntLam);
                        Bc_[i].j_col.push_back(j_col_Bc_curr);
                        Bc_[i].val.push_back(Bc_from_Rt_lk);
                        //_
                        j_col_Bc_neigh = data.g2l[IdNeighSub][dofs_l];
                        Bc_[IdNeighSub].i_coo_cmpr.push_back(cntLam);
                        Bc_[IdNeighSub].j_col.push_back(j_col_Bc_neigh);
                        Bc_[IdNeighSub].val.push_back(-Bc_from_Rt_lk);
                    }
                    cntLam++;
                }
            }
        }
    }
    matrix_Bx_COO2CSR(Bc_,cntLam);
}

void Cluster::create_Bc_or_Bf_in_CSR(vector < Matrix > &Bc,
                                     bool full_rank_or_redundant,
                                     bool cornersOnlyOrAllDof)
{
    bool addDirConstr = false;
    create_Bc_or_Bf_in_CSR(Bc, full_rank_or_redundant,cornersOnlyOrAllDof, addDirConstr);
}


void Cluster::create_Bc_or_Bf_in_CSR(vector < Matrix > &Bc_,
                                     bool full_rank_or_redundant,
                                     bool cornersOnlyOrAllDof,
                                     bool addDirConstr)
{

    int global_DOF;
    int ind_neigh_sub;
    int cntLam = 0;
    int j_col_Bc_curr;
    int j_col_Bc_neigh;

    vector<int>::iterator it;

    for (int d = 0; d < Bc_.size(); d++){
        for ( auto it1 = data.interface[d].begin(); it1 != data.interface[d].end(); ++it1  ){
            global_DOF = it1->first;

            if (cornersOnlyOrAllDof){
                it = find (mesh.cornerDOFs.begin(), mesh.cornerDOFs.end(), global_DOF);
                if (it == mesh.cornerDOFs.end())
                    continue;
            }

            if (addDirConstr){
                it = find (mesh.DirichletDOFs.begin(), mesh.DirichletDOFs.end(), global_DOF);
                if (it != mesh.DirichletDOFs.end()){

//                    j_col_Bc_curr = data.g2l[d][global_DOF];
//                    Bc_[d].l2g_i_coo.push_back(cntLam);
//                    Bc_[d].i_coo_cmpr.push_back(cntLam);
//                    Bc_[d].j_col.push_back(j_col_Bc_curr);
//                    Bc_[d].val.push_back(1);
//                    cntLam++;
//
//                    *it = -(*it);
                    continue;
                }
            }


            //TODO: if 'global_DOF' is found in vector 'cornerDOFs', the entry
            //                                                      should be deleted.

            for (int k = 0; k < it1->second.size(); k++){
                ind_neigh_sub = it1->second[k];
                if (ind_neigh_sub > d){
                    j_col_Bc_curr = data.g2l[d][global_DOF];
                    Bc_[d].l2g_i_coo.push_back(cntLam);
                    Bc_[d].i_coo_cmpr.push_back(cntLam);
                    Bc_[d].j_col.push_back(j_col_Bc_curr);
                    Bc_[d].val.push_back(1);

                    j_col_Bc_neigh = data.g2l[ind_neigh_sub][global_DOF];
                    Bc_[ind_neigh_sub].l2g_i_coo.push_back(cntLam);
                    Bc_[ind_neigh_sub].i_coo_cmpr.push_back(cntLam);
                    Bc_[ind_neigh_sub].j_col.push_back(j_col_Bc_neigh);
                    Bc_[ind_neigh_sub].val.push_back(-1);

                    if (addDirConstr){
                        weigth.push_back(1. / ( (double) it1->second.size() + 1 ) );
                    }
                    cntLam++;

                    if (full_rank_or_redundant)
                        break;
                }
            }
        }
    }

    if (addDirConstr){
        for (int d = 0; d < Bc_.size(); d++){
            for (int i = 0; i < mesh.DirichletDOFs.size();i++){
                int dirDof = mesh.DirichletDOFs[i];
                if (dirDof >= 0){
                    it = find (data.l2g[d].begin(), data.l2g[d].end(), dirDof);
                    if (it != data.l2g[d].end()){
                        j_col_Bc_curr = data.g2l[d][dirDof];
                        Bc_[d].l2g_i_coo.push_back(cntLam);
                        Bc_[d].i_coo_cmpr.push_back(cntLam);
                        Bc_[d].j_col.push_back(j_col_Bc_curr);
                        Bc_[d].val.push_back(1);
                        weigth.push_back(1);
                        cntLam++;
                    }
                }
            }
        }
    }

    matrix_Bx_COO2CSR(Bc_,cntLam);

}


void Cluster::matrix_Bx_COO2CSR(vector <Matrix> &Bc_, int cntLam){
    int _nnz;
    for (int d = 0 ; d < Bc_.size(); d++){
        _nnz = Bc_[d].val.size();
        Bc_[d].nnz = _nnz;
        Bc_[d].n_row = cntLam;
        Bc_[d].n_col = data.l2g[d].size();
        Bc_[d].symmetric = 0;
        Bc_[d].format = 0;

        vector<int>::iterator it;
        vector < int > l2g_(Bc_[d].i_coo_cmpr);
        sort(l2g_.begin(), l2g_.end());
        it = unique (l2g_.begin(), l2g_.end());
        l2g_.resize( distance(l2g_.begin(),it));
        Bc_[d].n_row_cmprs = l2g_.size();

        for (int i = 0; i < l2g_.size(); i++){
            Bc_[d].g2l_i_coo.insert ( pair < int, int > ( l2g_[i] , i ));
        }

        for (int i = 0; i < Bc_[d].i_coo_cmpr.size(); i++){
            Bc_[d].i_coo_cmpr[i] =
                    Bc_[d].g2l_i_coo[ Bc_[d].i_coo_cmpr[i] ];
        }
        Bc_[d].l2g_i_coo = l2g_;
        Bc_[d].COO2CSR();
    }
}

void Cluster::mult_Kplus_f(vector <Vector> & rhs_in , vector <Vector> & x_out){

    // gc will be allocated onece at the beginning of Cluster.cpp //
    Vector gc;
    int ngc = nLam_c + nRBM_c;

    if (options.solver_opt.Ac_extended_by_kerGc)
        ngc += nRBM_f;

    gc.zero_dense(ngc);

    for (int d = 0; d < nSubClst; d++){
        Vector BcKplusfc;
        BcKplusfc.mat_mult_dense(KplusBcT[d],"T",rhs_in[d],"N");
        for(int i = 0; i < BcKplusfc.n_row_cmprs; i++){
            gc.dense[Bc[d].l2g_i_coo[i]] += BcKplusfc.dense[i];
        }
    }

    int cntR = 0;
    for (int d = 0; d < nSubClst; d++){
        Vector Rtfc;
        Rtfc.mat_mult_dense(R[d],"T",rhs_in[d],"N");
        for (int i = 0 ; i < R[d].n_col; i++){
            gc.dense[nLam_c + cntR + i] = -Rtfc.dense[i];
        }
        cntR += R[d].n_col;
    }

//    gc.printToFile("gc",folder,0,true);

    Vector lam_alpha;
    lam_alpha = gc;
    lam_alpha.setZero();


    Ac_clust.solve_system(gc,lam_alpha);

//    lam_alpha.printToFile("lam_alpha",folder,0,true);



    cntR = 0;
    for (int d = 0; d < nSubClst ; d++){
        Vector lamc;
        lamc.zero_dense(Bc[d].n_row_cmprs);
        for (int i = 0; i < Bc[d].n_row_cmprs;i++){
           lamc.dense[i] = lam_alpha.dense[Bc[d].l2g_i_coo[i]];
        }
        Vector KplusBcTlamc;
        KplusBcTlamc.mat_mult_dense(KplusBcT[d],"N",lamc,"N");


        Vector Kplusfc;
        K[d].solve_system(rhs_in[d],Kplusfc);


        Vector Ralphac;
        Vector alphac_d; alphac_d.zero_dense(R[d].n_col);
        for (int i = 0; i < R[d].n_col; i++){
            alphac_d.dense[i] = lam_alpha.dense[nLam_c + cntR + i];
        }

//        alphac_d.printToFile("alpha",folder,d,true);
        Ralphac.mat_mult_dense(R[d],"N",alphac_d,"N");
        cntR += R[d].n_col;


        if (x_out[d].numel == 0)
            x_out[d].zero_dense(K[d].n_row_cmprs);


        for (int i = 0; i < x_out[d].n_row_cmprs; i++){
            x_out[d].dense[i] =  Kplusfc.dense[i] - KplusBcTlamc.dense[i] + Ralphac.dense[i];
        }
//        x_out[d].printToFile("x_out",folder,d,true);
    }
}


void Cluster::Preconditioning(Vector const & w_in , Vector & w_out){

    // gc will be allocated onece at the beginning of Cluster.cpp //

    if (w_out.numel < 0)
        w_out.zero_dense(w_in.n_row_cmprs);

    w_out = w_in;

    vector < Vector > xx, yy;
    xx.resize(nSubClst); yy.resize(nSubClst);

    scale(w_out);
    mult_BfT(w_in,xx);

    for (int d = 0; d < nSubClst; d++)
        K[d].mult(xx[d],yy[d],true);

    mult_Bf(yy,w_out);
    scale(w_out);
}

void Cluster::mult_Ff(Vector const & w , Vector & Fw){

    // gc will be allocated onece at the beginning of Cluster.cpp //


    vector < Vector > xx, yy;
    xx.resize(nSubClst); yy.resize(nSubClst);

    mult_BfT(w,xx);

    for (int d = 0; d < nSubClst; d++)
        K[d].solve_system(xx[d],yy[d]);
    mult_Bf(yy,Fw);

    Vector gc;
    int ngc = nLam_c + nRBM_c;

    if (options.solver_opt.Ac_extended_by_kerGc)
        ngc += nRBM_f;

    gc.zero_dense(ngc);

    for (int d = 0; d < nSubClst; d++){
        Vector BcKplusfc;
        BcKplusfc.mat_mult_dense(KplusBcT[d],"T",xx[d],"N");
        for(int i = 0; i < BcKplusfc.n_row_cmprs; i++){
            gc.dense[Bc[d].l2g_i_coo[i]] += BcKplusfc.dense[i];
        }
    }

    int cntR = 0;
    for (int d = 0; d < nSubClst; d++){
        Vector Rtfc;
        Rtfc.mat_mult_dense(R[d],"T",xx[d],"N");
        for (int i = 0 ; i < R[d].n_col; i++){
            gc.dense[nLam_c + cntR + i] = -Rtfc.dense[i];
        }
        cntR += R[d].n_col;
    }

    Vector lam_alpha;
    lam_alpha = gc;
    lam_alpha.setZero();

    Ac_clust.solve_system(gc,lam_alpha);


    cntR = 0;
    for (int d = 0; d < nSubClst ; d++){
        Vector lamc;
        lamc.zero_dense(Bc[d].n_row_cmprs);
        for (int i = 0; i < Bc[d].n_row_cmprs;i++){
           lamc.dense[i] = lam_alpha.dense[Bc[d].l2g_i_coo[i]];
        }
        Vector KplusBcTlamc;
        KplusBcTlamc.mat_mult_dense(KplusBcT[d],"N",lamc,"N");

        Vector alphac_d; alphac_d.zero_dense(R[d].n_col);
        for (int i = 0; i < R[d].n_col; i++){
            alphac_d.dense[i] = lam_alpha.dense[nLam_c + cntR + i];
        }
        xx[d].mat_mult_dense(R[d],"N",alphac_d,"N");
        for (int i = 0; i < xx[d].n_row_cmprs; i++){
            xx[d].dense[i] -= KplusBcTlamc.dense[i];
        }
        cntR += R[d].n_col;
    }



    Vector w1;
    mult_Bf(xx,w1);
    Fw.add(w1,1);
}




void Cluster::mult_Bf(vector < Vector > const &x_in, Vector &lambda ){


    if (lambda.numel == 0)
        lambda.zero_dense(nLam_f);
    else
        lambda.setZero();

    for (int d = 0; d < nSubClst; d++ ){
        Vector Bf_x_d;
        Bf[d].mult(x_in[d],Bf_x_d, true);
        for (int i = 0; i < Bf_x_d.n_row_cmprs; i++){
            lambda.dense[Bf[d].l2g_i_coo[i]] +=  Bf_x_d.dense[i];
        }
    }

}

void Cluster::mult_BfT(Vector const &lambda, vector < Vector > &x_out){

    for (int d = 0; d < nSubClst; d++){
        if (x_out[d].numel == 0)
            x_out[d].zero_dense(K[d].n_row_cmprs);
        else
            x_out[d].setZero();
        Vector lam_;
        lam_.zero_dense(Bf[d].n_row_cmprs);
        for (int i = 0; i < Bf[d].n_row_cmprs;i++){
            lam_.dense[i] = lambda.dense[Bf[d].l2g_i_coo[i]];
        }
        Bf[d].mult(lam_,x_out[d],false);
    }
}


void Cluster::create_GfTGf(){

    GfTGf.zero_dense(Gf[0].n_col, Gf[0].n_col);
    GfTGf.label = "GfTGf";

    Gf_clust.zero_dense(nLam_f,Rf[0].n_col);

    for (int d = 0; d < nSubClst; d++){
        for (int j = 0 ; j < Gf[d].n_col; j++){
            for (int i = 0 ; i < Gf[d].n_row_cmprs; i++){
                Gf_clust.dense[Bf[d].l2g_i_coo[i] + j * Gf_clust.n_row_cmprs] +=
                    Gf[d].dense[i + j * Gf[d].n_row_cmprs];
            }
        }
    }
    GfTGf.mat_mult_dense(Gf_clust,"T",Gf_clust,"N");
//
}

void Cluster::compute_invGfTGf(){

    char uplo = 'L';
    int n = GfTGf.n_row_cmprs;
    int lda = GfTGf.n_row_cmprs;
    int info;

    invGfTGf = GfTGf;
    info = LAPACKE_dpotrf (LAPACK_COL_MAJOR , uplo , n ,&(invGfTGf.dense[0]) , lda );

    if (info != 0)
        fprintf(stderr, "factorized matrix GfTGf failed. \n");



    info = LAPACKE_dpotri (LAPACK_COL_MAJOR, uplo , n , &(invGfTGf.dense[0]) , lda );

    if (info != 0)
        fprintf(stderr, "inverse of GfTGf failed. \n");


    for (int j = 0 ; j < invGfTGf.n_col - 1; j++) {
        for (int i = j ; i < invGfTGf.n_row_cmprs; i++) {
            invGfTGf.dense[j + i * invGfTGf.n_row_cmprs] =
                     invGfTGf.dense[i + j * invGfTGf.n_row_cmprs];
        }
    }

}

void Cluster::mult_Gf(Vector const & alpha, Vector & lambda){

    if (lambda.numel == 0){
        lambda.zero_dense(nRBM_f);
    }
    else{
        lambda.setZero();
    }
    lambda.mat_mult_dense(Gf_clust,"N",alpha,"N");
}

void Cluster::mult_GfT(Vector const & lambda, Vector & alpha){
    if (alpha.numel == 0){
        alpha.zero_dense(nRBM_f);
    }
    else{
        alpha.setZero();
    }
    alpha.mat_mult_dense(Gf_clust,"T",lambda,"N");
}

void Cluster::mult_RfT(vector <Vector> const &x_in, Vector & alpha){

    if (alpha.numel == 0){
        alpha.zero_dense(nRBM_f);
    }
    else{
        alpha.setZero();
    }

    for (int d = 0 ; d < nSubClst; d++){
        Matrix alpha_d;
        alpha_d.mat_mult_dense(Rf[d],"T",x_in[d],"N");
        for (int i = 0; i < alpha_d.n_row_cmprs; i++){
            alpha.dense[i] += alpha_d.dense[i];
        }
    }
}

void Cluster::Projection(Vector const & x, Vector & Px, Vector & alpha){

    Px = x;
    Vector GfTx;
    mult_GfT(x,GfTx);
    alpha.mat_mult_dense(invGfTGf,"N",GfTx,"N");
    mult_Gf(alpha,Px);

    for (int i = 0; i < x.n_row_cmprs; i++)
        Px.dense[i] = x.dense[i] - Px.dense[i];

    for (int i = 0 ; i < alpha.n_row_cmprs;i++)
        alpha.dense[i] *= -1;


}


void Cluster::create_Rf_and_Gf(){

    int cntR = 0;
    for (int d = 0; d < nSubClst; d++){
        Matrix Hd;
        Hd.zero_dense(R[d].n_col,kerGc.n_col);
        for (int j = 0; j < kerGc.n_col;j++){
            for (int i = 0; i < R[d].n_col;i++){
                Hd.dense[i + j * R[d].n_col] = kerGc.dense[cntR + i + j * kerGc.n_row_cmprs];
            }
        }
        cntR += R[d].n_col;
        Rf[d].mat_mult_dense(R[d],"N",Hd,"N");
        int n_rowGf = Bf[d].n_row_cmprs;
        int n_colGf = Rf[d].n_col;
        //TODO Gf (as vector) is no need to keep in the memory
        // later replace by Gf_clust (is already in create GfTGf)
        Gf[d].zero_dense(n_rowGf, n_colGf );
        Bf[d].mult(Rf[d],Gf[d],true);
    }
}



void Cluster::scale(Vector & x){

    for (int i = 0 ; i < x.n_row_cmprs; i++)
        x.dense[i] *= weigth[i];
}




void Cluster::pcpg2(){

    double gPz, gPz_prev, wFw, rho, gamma, wFPg, norm_gPz0;
    double eps_iter = 1e-4;
    Vector g0, d_rhs, e, iGTG_e, lambda, z, Pz;
    Vector Fw, Pg, g, w, w_prev;
    vector < Vector > xx, yy;
    Vector beta, alpha;
    xx.resize(nSubClst); yy.resize(nSubClst);

    xx[0].label = "test";
    mult_Kplus_f(rhs,xx);

    mult_Bf(xx,d_rhs);

    // e = Rt * f
    mult_RfT(rhs,e);

    // lambda0 = G * inv(GtG) * e
    iGTG_e.mat_mult_dense(invGfTGf,"N",e,"N");
    mult_Gf(iGTG_e, lambda);


    // F * lambda0
    mult_Ff(lambda,g0);
    // g0 = F * lambda0 - d_rhs
    g0.add(d_rhs,-1);

    g = g0;

    // Pg0
    Projection(g,Pg,alpha);
    Preconditioning(Pg,z);
    Projection(z,Pz,beta);
    gPz = Matrix::dot(g,Pz);

    norm_gPz0 = sqrt(gPz);

    w = Pz;


    int max_iter = 200;

    for (int it = 0; it < max_iter; it++){

        printf("it: %4d,      |Pg| = %3.9e \n",it + 1, sqrt(gPz) / norm_gPz0);
        if (sqrt(gPz) < eps_iter * norm_gPz0)
            break;

        mult_Ff(w,Fw);
        wFw = Matrix::dot(w,Fw);
        rho = -gPz / wFw;

        lambda.add(w,rho);
        g.add(Fw,rho);

        Projection(g,Pg,alpha);
        Preconditioning(Pg,z);
        Projection(z,Pz,beta);

        gPz_prev = gPz;
        gPz = Matrix::dot(g,Pz);

        gamma = gPz / gPz_prev;

        w_prev = w;
        w = Pz;
        w.add(w_prev,gamma);


        printVTK(yy, xx, lambda, alpha, it);

    }

    lambda.printToFile("lambda_",folder,0,true);

    // final solution (parameter 1000 is number > maxIter)
    printVTK(yy, xx, lambda, alpha, 1000);




}



void Cluster::pcpg(){

    double gPg, wFw, rho, gamma, wFPg, norm_gPg0;
    double eps_iter = 1e-4;
    Vector g0, d_rhs, e, iGTG_e, lambda, z, Pz;
    Vector Fw, Pg, g, w, w_prev;
    vector < Vector > xx, yy;
    Vector alpha;
    xx.resize(nSubClst); yy.resize(nSubClst);
    // d_rhs = B * Kplus * f

    xx[0].label = "test";
    mult_Kplus_f(rhs,xx);

//    for (int d = 0; d < nSubClst; d++){
//        xx[d].printToFile("xx_drhs",folder,d,true);
//    }


    mult_Bf(xx,d_rhs);


//    d_rhs.printToFile("d_rhs",folder,0,true);

    bool new_feti_operator = true;

#if 1
    // e = -Rt * f
    mult_RfT(rhs,e);


    // lambda0 = G * inv(GtG) * e
    iGTG_e.mat_mult_dense(invGfTGf,"N",e,"N");
    mult_Gf(iGTG_e, lambda);


    // F * lambda0
    if (new_feti_operator){
        mult_Ff(lambda,g0);
    }
    else{
        mult_BfT(lambda,xx);
        mult_Kplus_f(xx,yy);
        mult_Bf(yy,g0);
    }
    // g0 = F * lambda0 - d_rhs
    g0.add(d_rhs,-1);

    g = g0;

    // Pg0
    Projection(g,Pg,alpha);


    norm_gPg0 = sqrt(Matrix::dot(Pg,Pg));
    // initial conjugate vector
    w = Pg;

    int max_iter = 200;
    for (int it = 0; it < max_iter; it++){


        gPg = Matrix::dot(Pg,Pg);
        printf("it: %4d,      |Pg| = %3.9e \n",it+1, sqrt(gPg) / norm_gPg0);

        if (sqrt(gPg) < eps_iter * norm_gPg0)
            break;
        // F * w
        if (new_feti_operator){
            mult_Ff(w,Fw);
        }
        else{
            mult_BfT(w,xx);
            mult_Kplus_f(xx,yy);
            mult_Bf(yy,Fw);
        }

        wFw = Matrix::dot(w,Fw);

        rho = -gPg / wFw;

        lambda.add(w,rho);

        // new gradient
        g.add(Fw,rho);

        // P * g    where P = I - G * inv(GtG) *Gt

        Projection(g,Pg,alpha);






// //////////////////////////////////////////////////////////////////////////////////
#if 0
        mult_BfT(lambda,yy);
        for (int d = 0; d < nSubClst;d++){
            for (int i = 0; i < yy[d].n_row_cmprs; i++){
                yy[d].dense[i] = rhs[d].dense[i] - yy[d].dense[i];
            }
        }
        mult_Kplus_f(yy,xx);
        for (int d = 0; d < nSubClst;d++){
            Vector Rf_alpha;
            Rf_alpha.mat_mult_dense(Rf[d],"N",alpha,"N");
            for (int i = 0; i < xx[d].n_row_cmprs; i++){
                xx[d].dense[i] -= Rf_alpha.dense[i];
            }
        }



        vector < double > solution;
        solution.resize(mesh.nPoints * 3, 0 );

        for (int d = 0 ; d < nSubClst; d++){
            for (int i = 0; i < xx[d].n_row_cmprs; i++){
                solution[data.l2g[d][i]] = xx[d].dense[i];
            }
        }
        mesh.SaveVTK(solution, folder,it);
#endif
// //////////////////////////////////////////////////////////////////////////////////

        wFPg = Matrix::dot(Fw,Pg);

        gamma = -wFPg / wFw;

        w_prev = w;

        w = Pg;
        w.add(w_prev,gamma);
    }

    lambda.printToFile("lambda_",folder,0,true);

    Vector uDecomp;
    uDecomp.zero_dense(rhs[0].n_row_cmprs * nSubClst);

    int cnt = 0;
    for (int d = 0 ; d < nSubClst; d++){
        for (int i = 0; i < xx[d].n_row_cmprs; i++){
            uDecomp.dense[cnt] = xx[d].dense[i];
            cnt++;
        }
    }
    uDecomp.printToFile("uDecomp",folder,0,true);


// ================================================================
    if (1){
    mult_BfT(lambda,yy);
    for (int d = 0; d < nSubClst;d++){
        for (int i = 0; i < yy[d].n_row_cmprs; i++){
            yy[d].dense[i] = rhs[d].dense[i] - yy[d].dense[i];
        }
    }
    mult_Kplus_f(yy,xx);
    for (int d = 0; d < nSubClst;d++){
        Vector Rf_alpha;
        Rf_alpha.mat_mult_dense(Rf[d],"N",alpha,"N");
        for (int i = 0; i < xx[d].n_row_cmprs; i++){
            xx[d].dense[i] -= Rf_alpha.dense[i];
        }
    }

    vector < double > solution;
    solution.resize(mesh.nPoints * 3, 0 );

    for (int d = 0 ; d < nSubClst; d++){
        for (int i = 0; i < xx[d].n_row_cmprs; i++){
            solution[data.l2g[d][i]] = xx[d].dense[i];
        }
    }
    mesh.SaveVTK(solution, folder,0);
    }
// ================================================================




#endif
}





void Cluster::printVTK(vector < Vector > & yy, vector <Vector > & xx, Vector &lambda, Vector & alpha , int it){


        mult_BfT(lambda,yy);
        for (int d = 0; d < nSubClst;d++){
            for (int i = 0; i < yy[d].n_row_cmprs; i++){
                yy[d].dense[i] = rhs[d].dense[i] - yy[d].dense[i];
            }
        }
        mult_Kplus_f(yy,xx);
        for (int d = 0; d < nSubClst;d++){
            Vector Rf_alpha;
            Rf_alpha.mat_mult_dense(Rf[d],"N",alpha,"N");
            for (int i = 0; i < xx[d].n_row_cmprs; i++){
                xx[d].dense[i] -= Rf_alpha.dense[i];
            }
        }



        vector < double > solution;
        solution.resize(mesh.nPoints * 3, 0 );

        for (int d = 0 ; d < nSubClst; d++){
            for (int i = 0; i < xx[d].n_row_cmprs; i++){
                solution[data.l2g[d][i]] = xx[d].dense[i];
            }
        }
        mesh.SaveVTK(solution, folder,it);
}
