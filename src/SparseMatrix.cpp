#include "SparseMatrix.h"


using namespace std;

SparseMatrix::SparseMatrix(){ 
    this->nnz = 0;
    this->n_row = 0;
    this->n_col = 0;
    this->isCOO = false;
    this->isCSR = false;
    this->isDNS = false;
}

void SparseMatrix::readCooFromFile( std::string path2matrix)
{

    using namespace std;

    std::ifstream input(path2matrix.c_str());
    std::cout << path2matrix.c_str() <<"\n";

    input >> this->n_row;
    input >> this->n_col;
    input >> this->nnz;

    std::cout<< "n_row "<< this->n_row <<  std::endl;
    std::cout<< "n_col "<< this->n_col <<  std::endl;
    std::cout<< "nnz   "<< this->nnz   <<  std::endl;

    this->i_coo.resize(this->nnz );
    this->j_coo.resize(this->nnz );
    this->v_coo.resize(this->nnz );

    int int1;
    double double1;
    for (int i = 0; i < this->nnz; i++) {

        input >> int1;
        this->i_coo[i ] = int1;

        input >> int1;
        this->j_coo[i ] = int1;

        input >> double1;
        this->v_coo[i ] = double1;
    }
}

void SparseMatrix::COO2CSR(){ 



    MKL_INT mkl_nnz = this->i_coo.size();
    MKL_INT job[6] = {1,//If job(1)=1, the matrix in the coordinate format is converted to the CRS format.
                      1,//If job(2)=0, zero-based indexing for the matrix in CRS format is used;
                      1,//If job(3)=0, zero-based indexing for the matrix in coordinate format is used;
                      0,
                      mkl_nnz,//job(5)=nnz-sets number of the non-zero elements of the matrix A if job(1)=1.
                      0 //If job(6)=0, all arrays acsr, ja, ia are filled in for the output storage.
                     };
    m = this->n_row;
    this->a.resize(mkl_nnz);
    this->ia.resize(m+1);
    this->ja.resize(mkl_nnz);
    MKL_INT info;
    mkl_dcsrcoo(job , &m, &a[0] , &ja[0] , &ia[0] , &mkl_nnz ,
            &v_coo[0], &i_coo[0] , &j_coo[0] , &info );
    if (info != 0){
        std::cout << " PROBLEM: CONVERSION WAS NOT SUCCESSFULL! " <<
                     info << " ++++++++++++++++++++ \n ";
    }

    i_coo.clear();
    j_coo.clear();
    v_coo.clear();

    this->nnz = mkl_nnz;
    this->isCOO = false;
    this->isCSR = true;
    this->isDNS = false;
}

void SparseMatrix::CSR2COO(){ 


    MKL_INT job[6] = {0,//If job(1)=1, the matrix in the coordinate format is converted to the CRS format.
                      1,//If job(2)=0, zero-based indexing for the matrix in CRS format is used;
                      1,//If job(3)=0, zero-based indexing for the matrix in coordinate format is used;
                      0,
                      this->nnz,//job(5)=nnz-sets number of the non-zero elements of the matrix A if job(1)=1.
                      3 //If job(6)=0, all arrays acsr, ja, ia are filled in for the output storage.
                     };
    MKL_INT mkl_nnz = this->nnz;
    //    std::cout<< "---------------------------------------------" <<  mkl_nnz << "\n";
    m = this->n_row;
    this->i_coo.resize(mkl_nnz);
    this->j_coo.resize(mkl_nnz);
    this->v_coo.resize(mkl_nnz);
    MKL_INT info;
    mkl_dcsrcoo(job , &m, &a[0] , &ja[0] , &ia[0] , &mkl_nnz ,
            &v_coo[0], &i_coo[0] , &j_coo[0] , &info );

    if (info != 0){
        std::cout << " PROBLEM: CONVERSION WAS NOT SUCCESSFULL! " <<
                     info << " ++++++++++++++++++++ \n ";
    }
    a.clear();
    ja.clear();
    ia.clear();
    
    this->isCOO = true;
    this->isCSR = false;
    this->isDNS = false;

}

void SparseMatrix::CSR2DNS(int isUpperLowerOrBoth){

    MKL_INT job[6] = {1,// job[0] = 0, the rectangular matrix A is converted to the CSR format;.
                      1,// job[1] = 0, zero-based indexing for the rectangular matrix
                      1,// job[2] = 0, zero-based indexing for the matrix in CSR format is used;
                      isUpperLowerOrBoth,// job[3] = 0-lower, 1-upper, 2-whole
                      this->nnz,//job[4]=nnz-sets number of the non-zero elements of the matrix A if job(1)=1.
                      3 // job[5] = 0, relates to DNS -> CSR (not
                     };
    //    MKL_INT mkl_nnz = this->nnz;
    //    std::cout<< "---------------------------------------------" <<  mkl_nnz << "\n";
    m = this->n_row;
    this->adns.resize(this->n_row * this->n_col);
    MKL_INT info;
    MKL_INT lda = this->n_row;
    MKL_INT m = this->n_row;
    MKL_INT n = this->n_col;

    mkl_ddnscsr (job, &m , &n, &adns[0], &lda, &a[0], &ja[0], &ia[0], &info);

    a.clear();
    ja.clear();
    ia.clear();

    this->isCOO = false;
    this->isCSR = false;
    this->isDNS = true;
}

void SparseMatrix::DNS2CSR(int isUpperLowerOrBoth){

    MKL_INT job[6] = {0,// job[0] = 0, the rectangular matrix A is converted to the CSR format;.
                      1,// job[1] = 0, zero-based indexing for the rectangular matrix
                      1,// job[2] = 0, zero-based indexing for the matrix in CSR format is used;
                      isUpperLowerOrBoth,// job[3] = 0-lower, 1-upper, 2-whole
                      this->nnz,//job[4]=nnz-sets number of the non-zero elements of the matrix A if job(1)=1.
                      3 // job[5] = 0, relates to DNS -> CSR (not
                     };
    MKL_INT info;
    MKL_INT lda = this->n_row;
    MKL_INT m = this->n_row;
    MKL_INT n = this->n_col;

    this->a.resize(this->nnz);
    this->ia.resize(this->n_row+1);
    this->ja.resize(this->nnz);



    mkl_ddnscsr (job, &m , &n, &adns[0], &lda, &a[0], &ja[0], &ia[0], &info);


    adns.clear();

    this->isCOO = false;
    this->isCSR = false;
    this->isDNS = true;
}

void SparseMatrix::dcsradd(SparseMatrix& A, SparseMatrix& B, SparseMatrix& C){
    
    char  trans = 'N';
    MKL_INT request = 1;
    MKL_INT sort    = 0;
    MKL_INT mA = A.n_row;
    MKL_INT nA = A.n_col;
    double beta = 1;
    MKL_INT nzmax;
    MKL_INT info;

    C.ia.resize(mA + 1);



    mkl_dcsradd (&trans, &request , &sort , &mA , &nA ,
                 &A.a[0] , &A.ja[0], &A.ia[0],
            &beta,
            &B.a[0] , &B.ja[0], &B.ia[0],
            &C.a[0] , &C.ja[0], &C.ia[0],
            &nzmax , &info );

    MKL_INT dimOf_a_and_ja = C.ia[mA] - 1;
    std::cout << " dimOf_a_and_ja = "<< dimOf_a_and_ja <<" \n ";
    C.a.resize(dimOf_a_and_ja);
    C.ja.resize(dimOf_a_and_ja);
    request = 2;

    mkl_dcsradd (&trans, &request , &sort , &mA , &nA ,
                 &A.a[0] , &A.ja[0], &A.ia[0],
            &beta,
            &B.a[0] , &B.ja[0], &B.ia[0],
            &C.a[0] , &C.ja[0], &C.ia[0],
            &nzmax , &info );
    if (info != 0){
        std::cout << " PROBLEM: CONVERSION WAS NOT SUCCESSFULL! " <<
                     info << " ++++++++++++++++++++ \n";
    }
    C.nnz = dimOf_a_and_ja;
    C.n_row = A.n_row;
    C.n_col = A.n_col;
    std::cout << " A --------------" << A.ja.size() << "\n";
    std::cout << " B --------------" << B.ja.size() << "\n";
    std::cout << " C --------------" << C.ja.size() << "\n";

}

void SparseMatrix::Acsr_mult_Bdns_is_Cdns(SparseMatrix& A, SparseMatrix& B, SparseMatrix& C){
    char transa = 'N';
    MKL_INT m = A.n_row;
    MKL_INT n = B.n_col;
    MKL_INT k = A.n_col;
    double alpha = 1;
    char matdescra[6] = {'G',' ',' ','F'};
    MKL_INT ldb = B.n_row;
    double beta = 0;
    C.adns.resize(m * n);
    MKL_INT ldc = A.n_row;
    mkl_dcsrmm (&transa , &m , &n , &k , &alpha , matdescra , &A.a[0], &A.ja[0] ,
            &A.ia[0] , &A.ia[1] , &B.adns[0] , &ldb , &beta , &C.adns[0], &ldc );
    C.n_row = m;
    C.n_col = n;
    C.isCOO = false;
    C.isCSR = false;
    C.isDNS = true;
}



void SparseMatrix::InitializeSolve(){
    mtype = -2;              /* Real symmetric matrix */
    nrhs = 1;               /* Number of right hand sides. */
    MKL_INT i;
    double ddum;            /* Double dummy */
    MKL_INT idum;           /* Integer dummy. */
    /* -------------------------------------------------------------------- */
    /* .. Setup Pardiso control parameters. */
    /* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Not in use */
    iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off
                             (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    iparm[34] = 0;
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Which factorization to use. */
    msglvl = 3;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */
    /* -------------------------------------------------------------------- */
    /* .. Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
    n = this->n_row;
    MKL_INT n1 = this->n_row;
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n1, &this->a[0], &this->ia[0], &this->ja[0],
            &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during symbolic factorization: %d", error);
        exit (1);
    }

    phase = 22;

    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n1, &this->a[0], &this->ia[0], &this->ja[0],
            &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);


    //	PARDISO (id->pt, &(id->maxfct), &(id->mnum), &(id->mtype), &phase,
    //	&(id->n), id->a, id->ia, id->ja, &(id->idum), &(id->nrhs),
    //	id->iparm, &(id->msglvl),  &(id->ddum),  &(id->ddum), &(id->error));




}


void SparseMatrix::solve(SparseMatrix& B, SparseMatrix& X){ 
    /* -------------------------------------------------------------------- */
    /* .. Back substitution and iterative refinement. */
    /* -------------------------------------------------------------------- */
    phase = 33;
    iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
    MKL_INT idum;         /* Integer dummy. */
    nrhs = B.n_col;
    X.adns.resize(B.adns.size());
    X.n_row = B.n_row;
    X.n_col = B.n_col;
    X.nnz = B.nnz;

    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &this->a[0], &this->ia[0], &this->ja[0], &idum, &nrhs, iparm, &msglvl,
            &B.adns[0], &X.adns[0], &error);

//            &B.adns[0], &X.adns[0], &error);
    if ( error != 0 )
    {
        printf ("\nERROR during solution: %d", error);
        exit (3);
    }
    printf ("\nSolve completed ... ");
}

void SparseMatrix::FinalizeSolve(int _i){

    double ddum;          /* Double dummy */
    MKL_INT idum;         /* Integer dummy. */
    phase = -1;           /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, &this->ia[0], &this->ja[0], &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);

    if (_i == 0 ) {
        for (int i = 0 ; i < 64 ; i++ ){
            if (iparm[i] < 0 ){
                std::cout <<"iparm[" << i <<"] = " << iparm[i] << "\n";
            }
        }
    }

}



void SparseMatrix::RemoveLower(){

    this->CSR2COO();
    int i = this->i_coo.size() - 1;
    for (i = this->i_coo.size() - 1; i>=1 ; i--){
        if (this->i_coo[i] > this->j_coo[i]){
            this->i_coo.erase(this->i_coo.begin() + i);
            this->j_coo.erase(this->j_coo.begin() + i);
            this->v_coo.erase(this->v_coo.begin() + i);
            this->nnz -= 1;
        }
    }
    this->COO2CSR();
}

void SparseMatrix::printToFile(std::string nameOfMat, int indOfMat){
//
    std::string path2matrix = "data/dump_" + nameOfMat + "_" + std::to_string(indOfMat)+".txt";
    this->CSR2COO();
    std::ofstream myfile_K_reg (path2matrix);
    std::cout.precision(10);
    myfile_K_reg << this->n_row << "\t"
                 << this->n_col << "\t"
                 << this->nnz << "\n" ;
    for (int i = 0; i < this->nnz ; i++ ){
        myfile_K_reg << this->i_coo[i] << "\t"
                     << this->j_coo[i] << "\t" << std::setprecision (17)
                     << this->v_coo[i] << "\n" ;
    }
    myfile_K_reg.close();
    this->COO2CSR();
}


bool cmp_int_int(int_int_dbl a,int_int_dbl b) { return (a.I < b.I); }

void SparseMatrix::compresRows(){
    std::vector < int_int_dbl > tmpVec;

    tmpVec.resize(this->i_coo.size());
    l2g_i_coo.resize(this->i_coo.size());

    for (int i = 0; i < this->i_coo.size();i++){
        tmpVec[i].I = this->i_coo[i];
        tmpVec[i].J = this->j_coo[i];
        tmpVec[i].V = this->v_coo[i];
    }
    sort(tmpVec.begin(),tmpVec.end(),cmp_int_int);

    int prevInd = -1;
    int counter = 0;
    for (int i = 0 ; i < this->i_coo.size(); i ++){
        if (prevInd != tmpVec[i].I){
            this->l2g_i_coo[counter] = tmpVec[i].I;
            counter++;
        }
        prevInd = tmpVec[i].I;
        this->i_coo[i] = counter;      // OK if one-based numbering, otherwise = (counter - 1)
        this->j_coo[i] = tmpVec[i].J;
        this->v_coo[i] = tmpVec[i].V;
    }
    tmpVec.clear();
}

void SparseMatrix::CreateCopyFrom(const SparseMatrix&AtoBeCopied, int is_N_or_T){
  // must be in COO format


    if (is_N_or_T == 1){
        this->i_coo = AtoBeCopied.i_coo;
        this->j_coo = AtoBeCopied.j_coo;
        this->n_row = AtoBeCopied.n_row;
        this->n_col = AtoBeCopied.n_col;
    }
    else
    {
        this->i_coo = AtoBeCopied.j_coo;
        this->j_coo = AtoBeCopied.i_coo;
        this->n_row = AtoBeCopied.n_col;
        this->n_col = AtoBeCopied.n_row;
    }
    this->nnz = AtoBeCopied.nnz;


}

void SparseMatrix::testPardiso(){
/*------------------------------ TEST ------------------------------*/
    SparseMatrix A;
    int n = 21;
    int cnt = 0;
    int nrhs = 5;
    A.nnz = 2 * n - 1;
    A.n_row = n;
    A.n_col = n;
    SparseMatrix B;
    B.n_row = n;
    B.n_col = nrhs;
    B.adns.resize(B.n_row * B.n_col);
    B.nnz = B.n_row * B.n_col;
    double h = 1. / (n + 1);

    for (int i = 0; i < n; i ++){
        for (int j = 0; j < nrhs; j ++){
            B.adns[ nrhs * i + j] =  h * h;
        }
        for (int j = 0; j < n; j ++){
            if (i == j){
                A.i_coo.push_back(i + 1);
                A.j_coo.push_back(j + 1);
                A.v_coo.push_back(2);
                cnt++;
            }
            else if ((j-i) == 1){
                A.i_coo.push_back(i + 1);
                A.j_coo.push_back(j + 1);
                A.v_coo.push_back(-1);
                cnt++;
            }
        }
    }
    A.COO2CSR();
    A.printToFile("A",1000);
    A.InitializeSolve();
    SparseMatrix U;
    A.solve(B, U);
    double x_current = h;
    double u_current, delta = 0;
    for (int i = 0; i < n; i++){
        u_current = 0.5 * x_current * (1. - x_current);
        delta += (u_current - U.adns[i]) * (u_current - U.adns[i]);
        if (n < 50){
            std::cout << u_current << " " << U.adns[i] <<"\n";
        }
        x_current += h;

    }
    std::cout << "|| exact - numerical solution || = " << sqrt(delta) << "\n";

    U.DNS2CSR(2);
    U.printToFile("X",1000);
    B.DNS2CSR(2);
    B.printToFile("B",1000);
    A.FinalizeSolve(0);
/*------------------------------ TEST ------------------------------*/
}
