#include "../include/Matrix.hpp"


using namespace std;

bool cmp_int_int_I(int_int_dbl a,int_int_dbl b) { return (a.I < b.I); }
bool cmp_int_int_J(int_int_dbl a,int_int_dbl b) { return (a.J < b.J); }


Matrix::Matrix(){
    this->nnz = 0;
    this->n_row = 0;
    this->n_col = 0;
    this->n_row_cmprs= 0;
    this->isCOO = false;
    this->isCSR = false;
    this->isDNS = false;
    this->DNS_reducedZeroRows = false;
    this->DNS_transposed = false;
}

void Matrix::readCooFromFile(string path2matrix, int symmetric, int _format,
                                                                     int offset)
{
/* MATRIX MARKET FORMAT                                                       */
/* ASSUMPTION! SIZE OF I_COO_CMPR (J_COL, VAL) IS EQUAL TO NNZ                */
/* TODO OTHERWISE USE SORT + UNIQUE !!!                                       */

/* SYMMETRIC = 0: NOT SYMMETRIC (GENERAL)                                     */
/* SYMMETRIC = 1: LOWER TRIANGULAR                                            */
/* SYMMETRIC = 2: UPPER TRIANGULAR                                            */

/* FORMAT = 0: COO                                                            */
/* FORMAT = 1: CSR                                                            */
/* FORMAT = 2: DNS                                                            */

    ifstream input(path2matrix.c_str());
    cout << path2matrix.c_str() <<"\n";
    input >> this->n_row;
    input >> this->n_col;
    input >> this->nnz;
    cout<< "n_row "<< this->n_row << endl;
    cout<< "n_col "<< this->n_col << endl;
    cout<< "nnz   "<< this->nnz   << endl;
    this->symmetric = symmetric;

    if (_format == 2)    /*DNS*/
    {
        this->dense.resize(this->n_row * this->n_col);
        this->n_row_cmprs = this->n_row;
        this->DNS_reducedZeroRows = false;
        this->DNS_transposed = false;
        this->l2g_i_coo.resize(this->n_row_cmprs);
        for (int i = 0 ; i <this->n_row_cmprs; i++){
            this->l2g_i_coo[i] = i;
        }
    }
    else                /*COO or CSR*/
    {
        this->i_coo_cmpr.resize(this->nnz);
        this->j_col.resize(this->nnz);
        this->val.resize(this->nnz);
    }
/* Reading I, J, V ----------------------------------------------------------*/
    int I,J;
    int cnt = 0;
    double V;
    double _nnz = this->nnz; /* new nnz */
    for (int i = 0; i < this->nnz; i++) {
        input >> I;
        input >> J;
        input >> V;

        I-=offset;
        J-=offset;

        if (_format == 2) /* dense */
        {
            this->dense[I + J * this->n_row] = V;
        }
        else             /* COO or CSR */
        {
/* removing lower triangular or upper triangular-----------------------------*/
            if ((I > J && symmetric == 2) || (I < J && symmetric == 1))
            {
/* reduction of nnz (due to removed lower/upper part)------------------------*/
                _nnz--;
            }
            else
            {
                this->i_coo_cmpr[cnt] = I;
                this->j_col[cnt] = J;
                this->val[cnt] = V;
                cnt++;
            }
        }
    }
/* Update of 'nnz' if lower/upper triang. stored only ------------------------*/
/* and update size of vectors (i_coo_cmpr, j_col, ....)-----------------------*/
    this->nnz = _nnz;
/*                                                                            */
    this->format = _format;
    if (_format != 2){
        if (this->i_coo_cmpr.size() > this->nnz){
            this->i_coo_cmpr.resize(this->nnz);
            this->i_coo_cmpr.shrink_to_fit();

            this->j_col.resize(this->nnz);
            this->j_col.shrink_to_fit();

            this->val.resize(this->nnz);
            this->val.shrink_to_fit();
        }
/* Sorting [I, J, V]  --------------------------------------------------------*/
        if (_format != 2){
            vector < int_int_dbl > tmpVec;
            tmpVec.resize(this->i_coo_cmpr.size());
            for (int i = 0; i < this->i_coo_cmpr.size();i++){
                tmpVec[i].I = this->i_coo_cmpr[i];
                tmpVec[i].J = this->j_col[i];
                tmpVec[i].V = this->val[i];
            }
/* sort according to index I -------------------------------------------------*/
            sort(tmpVec.begin(),tmpVec.end(),cmp_int_int_I);
/* partial sorting according to index J (per sets with same I)----------------*/
            int startInd = 0, endInd = 1;
            int tmpVecIprev = tmpVec[0].I;
            for (int i = 1 ; i < this->i_coo_cmpr.size(); i ++){
                if (tmpVec[i].I == tmpVecIprev){
                    endInd++;
                }
                if (tmpVec[i].I != tmpVecIprev || (tmpVec[i].I == tmpVecIprev &&
                         i == this->i_coo_cmpr.size() - 1))
                {
                    sort(tmpVec.begin() + startInd,
                                      tmpVec.begin() + (endInd  ),cmp_int_int_J);
                    startInd = i;
                    endInd = i + 1;
                }
                tmpVecIprev = tmpVec[i].I;
            }
    /* Cumulating duplicated A[I,J] elements--------------------------------------*/
            int prevInd_I = -1;
            int prevInd_J = -1;
            int counter = 0;
            int cnt_j = -1;
            l2g_i_coo.resize(this->i_coo_cmpr.size());
            for (int i = 0 ; i < this->i_coo_cmpr.size(); i ++){
                if (prevInd_I != tmpVec[i].I){
                    this->l2g_i_coo[counter] = tmpVec[i].I;
                    counter++;
                }
                if (prevInd_I == tmpVec[i].I && prevInd_J == tmpVec[i].J){
                    this->val[cnt_j] += tmpVec[i].V;
                    this->nnz--;
                }
                else {
                    cnt_j++;
                    this->i_coo_cmpr[cnt_j] = counter - 1;
                    this->j_col[cnt_j] = tmpVec[i].J;
                    this->val[cnt_j] = tmpVec[i].V;
                }
                prevInd_I = tmpVec[i].I;
                prevInd_J = tmpVec[i].J;
            }
            this->l2g_i_coo.resize(counter );
            this->l2g_i_coo.shrink_to_fit();
            this->n_row_cmprs = counter;
            tmpVec.clear();
            tmpVec.shrink_to_fit();

        }
        if (_format == 1){
            this->COO2CSR();
        }
    }
}

void Matrix::COO2CSR(){
    if (this->format == 0){
        return;
    }
    int pi = -1;
    int cnt = 0;
    this->i_ptr.resize(this->n_row_cmprs + 1);
    for (int i = 0; i < this->nnz; i++){
        if (this->i_coo_cmpr[i] != pi){
            this->i_ptr[cnt] = i;
            cnt++;
        }
        pi = this->i_coo_cmpr[i];
    }

    this->i_ptr[this->n_row_cmprs] = this->nnz;
    this->i_coo_cmpr.clear();
    this->i_coo_cmpr.shrink_to_fit();
    this->format = 1;
#if 0
/* test - print 'i_ptr' to cmd                                                */
    for (int i = 0 ; i < this->n_row_cmprs + 1; i++){
        cout << " i_ptr["<< i<<"]="
             <<  this->i_ptr[i] << ":  "
             << this->n_row_cmprs + 1 <<"\n";
    }
#endif
}
void Matrix::CSR2COO(){
    if (this->format == 0){
        return;
    }
/* matrix does not contain zero rows                                          */
    this->i_coo_cmpr.resize(this->nnz);
//
    for (int i = 0; i < this->n_row_cmprs; i++) {
        for (int j = this->i_ptr[i]; j < this->i_ptr[i + 1]; j++) {
            this->i_coo_cmpr[j] = i;
        }
    }
    this->i_ptr.clear();
    this->i_ptr.shrink_to_fit();
    this->format = 0;
}



void Matrix::dummyFunction(int _n_row, int I, int J, double V)
{
    if (!this->DNS_reducedZeroRows){
        I = this->l2g_i_coo[I];
    }
    if (!this->symmetric){
        if (this->DNS_transposed){
            this->dense[I * this->n_col + J] = V;
        }
        else{
            this->dense[ I + J * _n_row] = V;
        }
    }
    else{
        this->dense[ I + J * _n_row] = V;
        if (I != J){
            this->dense[ J + I * _n_row] = V;
        }
    }
}



void Matrix::CSRorCOO2DNS(bool _reduceZeroRows, bool _transpose){
/*                                                                            */
    if (this->format == 2){
        /* matrix is in required format, nothing to do                        */
        return;
    }
/*----------------------------------------------------------------------------*/
    int _format = this->format;
    int I, J;
    double V;
    this->DNS_reducedZeroRows = _reduceZeroRows;
    this->DNS_transposed = _transpose;
    int _n_row;
    if (this->DNS_reducedZeroRows ){
        _n_row = this->n_row_cmprs;
    }
    else{
        _n_row = this->n_row;
    }

    this->dense.resize(n_row * this->n_col);
    for (int i = 0; i < _n_row * this->n_col;i++){
        this->dense[i] = 0;
    }

    if (this->format == 0){
        for (int i = 0; i < this->nnz; i++) {
            I = this->i_coo_cmpr[i];
            J = this->j_col[i];
            V = this->val[i];
            dummyFunction(_n_row,I,J,V);
        }
    }
    else if (this->format == 1){
        for (int i = 0; i < this->n_row_cmprs; i++) {
            for (int j = this->i_ptr[i]; j < this->i_ptr[i + 1]; j++) {
                I = i;
                J = this->j_col[j];
                V = this->val[j];
                dummyFunction(_n_row,I,J,V);
            }
        }
    }
    else {
        cerr << "Unexpected format of matrix to be converted to dense" <<
                                                           " format!!!" << "\n";
    }

    if(_format == 0){
        this->i_coo_cmpr.clear();
        this->i_coo_cmpr.shrink_to_fit();
    }
    else if (_format == 1){
        this->i_ptr.clear();
        this->i_ptr.shrink_to_fit();
    }

    this->j_col.clear();
    this->j_col.shrink_to_fit();

    this->val.clear();
    this->val.shrink_to_fit();

    this->format = 2;
}

void Matrix::DNS2CSR(){
    if (this->format == 1){
/* matrix is in required format, nothing to do                                */
        return;
    }
    this->DNS2COO();
    this->COO2CSR();
}

void Matrix::DNS2COO(){
    if (this->format == 0){
        return;
    }
    int _n_row;
    if (this->DNS_reducedZeroRows){
        _n_row = this->n_row_cmprs;
    }
    else{
        _n_row = this->n_row;
    }
    this->nnz = 0;
    for (int i = 0; i < _n_row * this->n_col;i++) {
       if (this->dense[i] != 0){
           this->nnz++;
       }
    }
    this->i_coo_cmpr.resize(this->nnz);
    this->j_col.resize(this->nnz);
    this->val.resize(this->nnz);
    double _val;
    int cnt = 0;
    cout << "_n_row " << _n_row << "\n";
    for (int i = 0; i < _n_row; i++) {
        for (int j = 0; j < this->n_col; j++) {
            if (this->DNS_transposed){
                _val = this->dense[i * this->n_col + j];
            }
            else{
                _val = this->dense[i + j * _n_row];
            }
            if (_val != 0){
                this->i_coo_cmpr[cnt] = i;
                this->j_col[cnt] = j;
                this->val[cnt] = _val;
                cnt++;
            }
        }
    }
    this->dense.clear();
    this->dense.shrink_to_fit();
    this->format = 0;
    this->n_row_cmprs = this->nnz;
}

//void Matrix::RemoveLower(){
// //
// //    this->CSR2COO();
// //    int i = this->i_coo_cmpr.size() - 1;
// //    for (i = this->i_coo_cmpr.size() - 1; i>=1 ; i--){
// //        if (this->i_coo_cmpr[i] > this->j_col[i]){
// //            this->i_coo_cmpr.erase(this->i_coo_cmpr.begin() + i);
// //            this->j_col.erase(this->j_col.begin() + i);
// //            this->val.erase(this->val.begin() + i);
// //            this->nnz -= 1;
// //        }
// //    }
// //    this->COO2CSR();
//}

void Matrix::printToFile(string nameOfMat, int indOfMat,
                                                        bool _printCooOrDense)
{
    const int _format = this->format;
    string path2matrix = "../data/dump_" + nameOfMat + "_" +
                                                     to_string(indOfMat)+".txt";
    FILE *fp = NULL;
    fp = fopen(path2matrix.c_str(), "w");
    if (fp == NULL) {
        cerr << "fail to open " << path2matrix.c_str() << "\n";
    }

    double V;
    if (_printCooOrDense){
        if (this->format == 2){
            this->DNS2COO();
        }
        int I, J,Itmp;
        fprintf(fp, "%d\t%d\t%d\n", this->n_row, this->n_col, this->nnz);
        if (this->format == 0){
            for (int i = 0; i < this->nnz; i++) {
                I = this->l2g_i_coo[this->i_coo_cmpr[i]];
                J = this->j_col[i];
                V = this->val[i];
                fprintf(fp, "%d\t%d\t%.16e\n", I,J,V);
            }
        }
        else if (this->format == 1){
            for (int i = 0; i < this->n_row_cmprs; i++) {
                for (int j = this->i_ptr[i]; j < this->i_ptr[i + 1]; j++) {
                    I = this->l2g_i_coo[i];
                    J = this->j_col[j];
                    V = this->val[j];
                    if (this->DNS_transposed){
                        Itmp = I;
                        J = I;
                        I = Itmp;
                    }
                    fprintf(fp, "%d\t%d\t%.16e\n", I,J,V);
                }
            }
        }
        else {
            cerr << "Unexpected format of matrix to be printed !!!" << "\n";
        }
        if (this->format == 0 && _format == 1){
            this->COO2CSR();
            cout << " AAAAAAAAAAAAAAA \n";
        }
        if (this->format == 1 && _format == 0){
            this->CSR2COO();
            cout << " BBBBBBBBBBBBBBB \n";
        }
        if (_format == 2){
            this->CSRorCOO2DNS(this->DNS_reducedZeroRows,this->DNS_transposed);
        }
        cout << " this->format=" << this->format << "_format =" << _format << "\n";
    }
    else{
        cout << " not implemented yet " << endl;
//        if (this->format < 2) {
//            this->CSRorCOO2DNS(true,true);
//        }
//        cout << "=====this->format======" << this->format << "\n";
//
//        int _n_row, _ni, _nj, _nm;
//        if (this->DNS_reducedZeroRows) {
//            _n_row = this->n_row_cmprs;
//        }
//        else{
//            _n_row = this->n_row;
//        }
//        if (this->DNS_transposed){
//            _ni = this->n_col;
//            _nj = _n_row;
//            _nm = this->n_col;
//        }
//        else{
//            _ni = _n_row;
//            _nj = this->n_col;
//            _nm = _n_row;
//        }
//
//        for (int i = 0; i < _ni; i++) {
//            for (int j = 0; j < _nj; j++){
//                V = this->dense[i + j * _nm];
//                fprintf(fp, "%.16f\t", V);
//            }
//            fprintf(fp, "\n");
//        }
//
//        if (_format != 2){
//            if (_format == 0) {
//                this->DNS2COO();
//            }
//            else if (_format == 1){
//                this->DNS2CSR();
//            }
//        }
    }

    fclose(fp);
}



/*
void Matrix::compresRows(){
    vector < int_int_dbl > tmpVec;

    tmpVec.resize(this->i_coo_cmpr.size());
    l2g_i_coo.resize(this->i_coo_cmpr.size());

    for (int i = 0; i < this->i_coo_cmpr.size();i++){
        tmpVec[i].I = this->i_coo_cmpr[i];
        tmpVec[i].J = this->j_col[i];
        tmpVec[i].V = this->val[i];
    }
    sort(tmpVec.begin(),tmpVec.end(),cmp_int_int_I);

    int prevInd = -1;
    int counter = 0;
    for (int i = 0 ; i < this->i_coo_cmpr.size(); i ++){
        if (prevInd != tmpVec[i].I){
            this->l2g_i_coo[counter] = tmpVec[i].I;
            counter++;
        }
        prevInd = tmpVec[i].I;
        this->i_coo_cmpr[i] = counter - 1;      // OK if one-based numbering, otherwise = (counter - 1)
        this->j_col[i] = tmpVec[i].J;
        this->val[i] = tmpVec[i].V;
    }
    tmpVec.clear();
}*/

Matrix Matrix::CreateCopyFrom(const Matrix&AtoBeCopied){
    Matrix Aout = AtoBeCopied;
    return Aout;
}


//void Matrix::mv_csr(Matrix X, Matrix &AX, bool matA_NorT, bool matX_NorT){
////
//    int n_rhs = 1;
//    bool NorT = false;
//    for (int i = 0; i < this->n_row_cmprs; i++) {
//        for (int j = this->i_ptr[i]; j < this->i_ptr[i + 1]; j++) {
//            if (this->symmetric > 0 ){
//                for (int k = 0; k < n_rhs; k++){
//                    AX.dense[i + k * n_row_cmprs] +=
//                            this->val[j] * X.dense[this->j_col[j] + k * n_row_cmprs];
//                    if (this->j_col[j] != i){
//                        AX.dense[this->j_col[j] + k * n_row_cmprs] +=
//                                      this->val[j] * X.dense[i + k * n_row_cmprs];
//                    }
//                }
//            }
//            else {
//                for (int k = 0; k < n_rhs; k++){
//                    if (NorT){
//                        AX.dense[i + k * n_row_cmprs] +=
//                                this->val[j] * X.dense[this->j_col[j] + k * n_row_cmprs];
//                    }
//                    else {
//                        AX.dense[this->j_col[j] + k * n_row_cmprs] +=
//                                      this->val[j] *
//                                X.dense[i + k * n_row_cmprs];
//                   }
//                }
//            }
//        }
//    }
//}

void Matrix::mv_csr(const double x[], double  Ax[], bool NorT, int n_rhs){
//





    for (int i = 0; i < this->n_row_cmprs; i++) {
        for (int j = this->i_ptr[i]; j < this->i_ptr[i + 1]; j++) {
            if (this->symmetric > 0 ){
                for (int k = 0; k < n_rhs; k++){
                    Ax[i + k * n_row_cmprs] +=
                            this->val[j] * x[this->j_col[j] + k * n_row_cmprs];
                    if (this->j_col[j] != i){
                        Ax[this->j_col[j] + k * n_row_cmprs] +=
                                      this->val[j] * x[i + k * n_row_cmprs];
                    }
                }
            }
            else {
                for (int k = 0; k < n_rhs; k++){
                    if (NorT){
                        Ax[i + k * n_row_cmprs] +=
                                this->val[j] * x[this->j_col[j] + k * n_row_cmprs];
                    }
                    else {
                        Ax[this->j_col[j] + k * n_row_cmprs] +=
                                      this->val[j] *
                                x[i + k * n_row_cmprs];
                   }
                }
            }
        }
    }
}








void Matrix::CsrElementByElement(){
//
//    bool symMatrix=true;
//    int *indK;
//    vector < vector < int > > forCSRformat;
//    forCSRformat.resize(fem->domain->neqSub,vector<int>(0));
//    int elemDOFs = 24;
//    for (int i=0;i<fem->domain->n_elementsSub;i++){
//      indK = stif_glob_number[i].ieq;
//      for (int j=0;j<elemDOFs;j++){
//        for (int k=0;k<elemDOFs;k++){
//          if (symMatrix && (indK[k]>=indK[j]) || !symMatrix){
//            forCSRformat[indK[j]].push_back(indK[k]);
//          }
//        }
//      }
//    }
//    nnz_K = 0;
//    for (int i = 0;i<forCSRformat.size();i++){
//      sort(forCSRformat[i].begin(),forCSRformat[i].end());
//      itv = unique (forCSRformat[i].begin(), forCSRformat[i].end());
//      forCSRformat[i].resize( distance(forCSRformat[i].begin(),itv) );
//      nnz_K+=forCSRformat[i].size();
//    }

}



void Matrix::InitializeSolve(){
#ifdef USE_PARDISO
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
             &n1, &this->val[0], &this->i_ptr[0], &this->j_col[0],
            &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during symbolic factorization: %d", error);
        exit (1);
    }

    phase = 22;

    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n1, &this->val[0], &this->i_ptr[0], &this->j_col[0],
            &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

    //	PARDISO (id->pt, &(id->maxfct), &(id->mnum), &(id->mtype), &phase,
    //	&(id->n), id->a, id->i_ptr, id->j_col, &(id->idum), &(id->nrhs),
    //	id->iparm, &(id->msglvl),  &(id->ddum),  &(id->ddum), &(id->error));
#endif
}


void Matrix::solve(Matrix& B, Matrix& X){
#ifdef USE_PARDISO
    /* -------------------------------------------------------------------- */
    /* .. Back substitution and iterative refinement. */
    /* -------------------------------------------------------------------- */
    phase = 33;
    iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
    MKL_INT idum;         /* Integer dummy. */
    nrhs = B.n_col;
    X.dense.resize(B.dense.size());
    X.n_row = B.n_row;
    X.n_col = B.n_col;
    X.nnz = B.nnz;

    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &this->val[0], &this->i_ptr[0], &this->j_col[0], &idum, &nrhs, iparm, &msglvl,
            &B.dense[0], &X.dense[0], &error);

    if ( error != 0 )
    {
        printf ("\nERROR during solution: %d", error);
        exit (3);
    }
    printf ("\nSolve completed ... ");
#endif
}

void Matrix::FinalizeSolve(int _i)
{
#ifdef USE_PARDISO

    double ddum;          /* Double dummy */
    MKL_INT idum;         /* Integer dummy. */
    phase = -1;           /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, &this->i_ptr[0], &this->j_col[0], &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);

    if (_i == 0 ) {
        for (int i = 0 ; i < 64 ; i++ ){
            if (iparm[i] < 0 ){
                cout <<"iparm[" << i <<"] = " << iparm[i] << "\n";
            }
        }
    }
#endif
}

void Matrix::testPardiso(){
#ifdef USE_PARDISO
/*------------------------------ TEST ------------------------------*/
    Matrix A;
    int n = 21;
    int cnt = 0;
    int nrhs = 5;
    A.nnz = 2 * n - 1;
    A.n_row = n;
    A.n_col = n;
    Matrix B;
    B.n_row = n;
    B.n_col = nrhs;
    B.dense.resize(B.n_row * B.n_col);
    B.nnz = B.n_row * B.n_col;
    double h = 1. / (n + 1);

    for (int i = 0; i < n; i ++){
        for (int j = 0; j < nrhs; j ++){
            B.dense[ nrhs * i + j] =  h * h;
        }
        for (int j = 0; j < n; j ++){
            if (i == j){
                A.i_coo_cmpr.push_back(i + 1);
                A.j_col.push_back(j + 1);
                A.val.push_back(2);
                cnt++;
            }
            else if ((j-i) == 1){
                A.i_coo_cmpr.push_back(i + 1);
                A.j_col.push_back(j + 1);
                A.val.push_back(-1);
                cnt++;
            }
        }
    }
    A.COO2CSR();
    A.printToFile("A",1000,false);
    A.InitializeSolve();
    Matrix U;
    A.solve(B, U);
    double x_current = h;
    double u_current, delta = 0;
    for (int i = 0; i < n; i++){
        u_current = 0.5 * x_current * (1. - x_current);
        delta += (u_current - U.dense[i]) * (u_current - U.dense[i]);
        if (n < 50){
            cout << u_current << " " << U.dense[i] <<"\n";
        }
        x_current += h;

    }
    cout << "|| exact - numerical solution || = " << sqrt(delta) << "\n";

    U.DNS2CSR();
    U.printToFile("X",1000,false);
    B.DNS2CSR();
    B.printToFile("B",1000,false);
    A.FinalizeSolve(0);
/*------------------------------ TEST ------------------------------*/
#endif
}

//void Matrix::dcsradd(Matrix& A, Matrix& B, Matrix& C){
//
//    char  trans = 'N';
//    MKL_INT request = 1;
//    MKL_INT sort    = 0;
//    MKL_INT mA = A.n_row;
//    MKL_INT nA = A.n_col;
//    double beta = 1;
//    MKL_INT nzmax;
//    MKL_INT info;
//
//    C.i_ptr.resize(mA + 1);
//
//
//
//    mkl_dcsradd (&trans, &request , &sort , &mA , &nA ,
//                 &A.a[0] , &A.j_col[0], &A.i_ptr[0],
//            &beta,
//            &B.a[0] , &B.j_col[0], &B.i_ptr[0],
//            &C.a[0] , &C.j_col[0], &C.i_ptr[0],
//            &nzmax , &info );
//
//    MKL_INT dimOf_a_and_ja = C.i_ptr[mA] - 1;
//    cout << " dimOf_a_and_ja = "<< dimOf_a_and_ja <<" \n ";
//    C.a.resize(dimOf_a_and_ja);
//    C.j_col.resize(dimOf_a_and_ja);
//    request = 2;
//
//    mkl_dcsradd (&trans, &request , &sort , &mA , &nA ,
//                 &A.a[0] , &A.j_col[0], &A.i_ptr[0],
//            &beta,
//            &B.a[0] , &B.j_col[0], &B.i_ptr[0],
//            &C.a[0] , &C.j_col[0], &C.i_ptr[0],
//            &nzmax , &info );
//    if (info != 0){
//        cout << " PROBLEM: CONVERSION WAS NOT SUCCESSFULL! " <<
//                     info << " ++++++++++++++++++++ \n";
//    }
//    C.nnz = dimOf_a_and_ja;
//    C.n_row = A.n_row;
//    C.n_col = A.n_col;
//    cout << " A --------------" << A.j_col.size() << "\n";
//    cout << " B --------------" << B.j_col.size() << "\n";
//    cout << " C --------------" << C.j_col.size() << "\n";
//
//}

//void Matrix::Acsr_mult_Bdns_is_Cdns(Matrix& A, Matrix& B, Matrix& C){
//    char transa = 'N';
//    MKL_INT m = A.n_row;
//    MKL_INT n = B.n_col;
//    MKL_INT k = A.n_col;
//    double alpha = 1;
//    char matdescra[6] = {'G',' ',' ','C'};
//    MKL_INT ldb = B.n_row;
//    double beta = 0;
//    C.dense.resize(m * n);
//    MKL_INT ldc = A.n_row;
//    mkl_dcsrmm (&transa , &m , &n , &k , &alpha , matdescra , &A.a[0], &A.j_col[0] ,
//            &A.i_ptr[0] , &A.i_ptr[1] , &B.dense[0] , &ldb , &beta , &C.dense[0], &ldc );
//    C.n_row = m;
//    C.n_col = n;
//    C.isCOO = false;
//    C.isCSR = false;
//    C.isDNS = true;
//}
