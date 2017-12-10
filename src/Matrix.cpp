#include "Matrix.hpp"
#include <math.h>
#include <iostream>


using namespace std;

Matrix::Matrix()
{
    Matrix::init();
}


Matrix:: Matrix(string label_)
{
    Matrix::init();
    label = label_;
}

Matrix::Matrix(int n_row_, int n_col_, bool NorT)
{
    Matrix::init();
    n_row = n_row_;
    n_col = n_col_;
    DNS_transposed = !NorT;
}

Matrix::Matrix(int n_row_, int n_col_, int nnz_, bool NorT)
{
    Matrix::init();
    n_row = n_row_;
    n_col = n_col_;
    nnz   = nnz_;
    DNS_transposed = !NorT;

}


void Matrix::init()
{
  label = "";
  tmp_label = "";
  n_row = 0;
  n_col = 0;
  DNS_transposed = false;
  nnz = 0;
  numel = 0;
  n_row_cmprs = 0;
  DNS_reducedZeroRows = false;
  diss_scaling = -1;
  msglvl = -1;
  solver = -1;
  linear_solver = "";
}


Matrix::~Matrix()
{}


bool Matrix::cmp_int_int_I(TRIPLET a,TRIPLET b)
{ return (a.I < b.I); }
bool Matrix::cmp_int_int_J(TRIPLET a,TRIPLET b)
{ return (a.J < b.J); }

bool Matrix::compareDouble(double i, double j)
{ return fabs(i)<=fabs(j); }



void Matrix::zero_dense(int n_row,int n_col){
    bool NorT = true;
    zero_dense(n_row, n_col, NorT);
}


void Matrix::zero_dense(int n_row_,int n_col_, bool NorT){
    nnz  = n_row_ * n_col_;
    numel = n_row_ * n_col_;
    n_row = n_row_;
    n_col = n_col_;
    n_row_cmprs = n_row_;
    DNS_reducedZeroRows = false;
    DNS_transposed = !NorT;
    symmetric = 0;
    format = 2;
    l2g_i_coo.resize(n_row_cmprs);
    for (int i = 0; i < n_row_cmprs; i++){
        l2g_i_coo[i] = i;
    }
    dense.resize(numel);
    for (int i = 0 ; i < numel; i++){
        dense[i] = 0;
    }
}


void Matrix::readCooFromFile(string path2matrix, int symmetric_, int format_,
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
    input >> n_row;
    input >> n_col;
    input >> nnz;
    symmetric = symmetric_;

#if VERBOSE_LEVEL>3
    cout << path2matrix.c_str() <<"\n";
    cout<< "n_row "<< n_row << endl;
    cout<< "n_col "<< n_col << endl;
    cout<< "nnz   "<< nnz   << endl;
#endif

    if (format_ == 2)    /*DNS*/
    {
        dense.resize(n_row * n_col);
        n_row_cmprs = n_row;
        DNS_reducedZeroRows = false;
        DNS_transposed = false;
        l2g_i_coo.resize(n_row_cmprs);
        for (int i = 0 ; i <n_row_cmprs; i++){
            l2g_i_coo[i] = i;
        }
    }
    else                /*COO or CSR*/
    {
        i_coo_cmpr.resize(nnz);
        j_col.resize(nnz);
        val.resize(nnz);
    }
/* Reading I, J, V ----------------------------------------------------------*/
    int I,J;
    int cnt = 0;
    double V;
    double _nnz = nnz; /* new nnz */
    for (int i = 0; i < nnz; i++) {
        input >> I;
        input >> J;
        input >> V;

        I-=offset;
        J-=offset;

        if (format_ == 2) /* dense */
        {
            dense[I + J * n_row] = V;
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
                i_coo_cmpr[cnt] = I;
                j_col[cnt] = J;
                val[cnt] = V;
                cnt++;
            }
        }
    }
/* Update of 'nnz' if lower/upper triang. stored only ------------------------*/
/* and update size of vectors (i_coo_cmpr, j_col, ....)-----------------------*/
    nnz = _nnz;
/*                                                                            */
    format = format_;
    if (format_ != 2){
        if (i_coo_cmpr.size() > nnz){
            i_coo_cmpr.resize(nnz);
//            i_coo_cmpr.shrink_to_fit();

            j_col.resize(nnz);
//            j_col.shrink_to_fit();

            val.resize(nnz);
//            val.shrink_to_fit();
        }
/* Sorting [I, J, V]  --------------------------------------------------------*/
        vector < TRIPLET > triplet;
        triplet.resize(i_coo_cmpr.size());
        for (int i = 0; i < i_coo_cmpr.size();i++){
            triplet[i].I = i_coo_cmpr[i];
            triplet[i].J = j_col[i];
            triplet[i].V = val[i];
        }
        sortAndUniqueCOO(triplet);
        triplet.clear();
//        triplet.shrink_to_fit();

        if (format_ == 1){
            COO2CSR();
        }
    }
    numel = n_row_cmprs * n_col;
}


void Matrix::sortAndUniqueCOO(vector <TRIPLET> &triplet){
/* sort according to index I -------------------------------------------------*/
    sort(triplet.begin(),triplet.end(),Matrix::cmp_int_int_I);
/* partial sorting according to index J (per sets with same I)----------------*/

  nnz = triplet.size();

    int startInd = 0, endInd = 1;
    int tripletIprev = triplet[0].I;
    for (int i = 1 ; i < nnz; i ++){
        if (triplet[i].I == tripletIprev){
            endInd++;
        }
        if (triplet[i].I != tripletIprev || (triplet[i].I == tripletIprev &&
                 i == nnz - 1))
        {
            sort(triplet.begin() + startInd,
                              triplet.begin() + (endInd  ),Matrix::cmp_int_int_J);
            startInd = i;
            endInd = i + 1;
        }
        tripletIprev = triplet[i].I;
    }


/* Cumulating duplicated A[I,J] elements--------------------------------------*/
    int prevInd_I = -1;
    int prevInd_J = -1;
    int counter = 0;
    int cnt_j = -1;

    if (l2g_i_coo.size() != triplet.size()){
        l2g_i_coo.resize(triplet.size());
    }
    if (i_coo_cmpr.size() != triplet.size()){
        i_coo_cmpr.resize(triplet.size());
    }
    if (j_col.size() != triplet.size()){
        j_col.resize(triplet.size());
    }
    if (val.size() != triplet.size()){
        val.resize(triplet.size());
    }


    int init_nnz = nnz;
    for (int i = 0 ; i < init_nnz; i ++){
        if (prevInd_I != triplet[i].I){
            l2g_i_coo[counter] = triplet[i].I;
            counter++;
        }
        if (prevInd_I == triplet[i].I && prevInd_J == triplet[i].J){
            val[cnt_j] += triplet[i].V;
            nnz--;
        }
        else {
            cnt_j++;
            i_coo_cmpr[cnt_j] = triplet[i].I ;
            j_col[cnt_j] = triplet[i].J;
            val[cnt_j] = triplet[i].V;
        }
        prevInd_I = triplet[i].I;
        prevInd_J = triplet[i].J;
    }

    l2g_i_coo.resize(counter );
//    l2g_i_coo.shrink_to_fit();

    i_coo_cmpr.resize(nnz);
//    i_coo_cmpr.shrink_to_fit();

    j_col.resize(nnz);
//    j_col.shrink_to_fit();

    val.resize(nnz);
//    val.shrink_to_fit();

    n_row_cmprs = counter;
    n_row = counter;

}


void Matrix::COO2CSR(){
//    if (format == 1){
//        return;
//    }

// COO is already is sorted without redundant entries
    int pi = -1;
    int cnt = 0;
    i_ptr.resize(n_row_cmprs + 1);
    for (int i = 0; i < nnz; i++){
        if (i_coo_cmpr[i] != pi){
            i_ptr[cnt] = i;
            cnt++;
        }
        pi = i_coo_cmpr[i];
    }

    i_ptr[n_row_cmprs] = nnz;
    i_coo_cmpr.clear();
//    i_coo_cmpr.shrink_to_fit();
    format = 1;
#if VERBOSE_LEVEL>4
/* test - print 'i_ptr' to cmd                                                */
    for (int i = 0 ; i < n_row_cmprs + 1; i++){
        cout << " i_ptr["<< i<<"]="
             <<  i_ptr[i] << ":  "
             << n_row_cmprs + 1 <<"\n";
    }
#endif
}
void Matrix::CSR2COO(){
    if (format == 0){
        return;
    }
/* matrix does not contain zero rows                                          */
    i_coo_cmpr.resize(nnz);
//
    for (int i = 0; i < n_row_cmprs; i++) {
        for (int j = i_ptr[i]; j < i_ptr[i + 1]; j++) {
            i_coo_cmpr[j] = i;
        }
    }
    i_ptr.clear();
//    i_ptr.shrink_to_fit();
    format = 0;
}



void Matrix::dummyFunction(int n_row_, int I, int J, double V)
{
    if (!DNS_reducedZeroRows){
        I = l2g_i_coo[I];
    }
    if (!symmetric){
        if (DNS_transposed){
            dense[I * n_col + J] = V;
        }
        else{
            dense[ I + J * n_row_] = V;
        }
    }
    else{
        dense[ I + J * n_row_] = V;
        if (I != J){
            dense[ J + I * n_row_] = V;
        }
    }
}


void Matrix::CSRorCOO2DNS(bool reduceZeroRows_, bool transpose_){
   int dbg_ = -1;
   CSRorCOO2DNS(reduceZeroRows_, transpose_,dbg_);

}


void Matrix::CSRorCOO2DNS(bool reduceZeroRows_, bool transpose_,int dbg){
/*                                                                            */

    if (dbg > 0){
        cout << "=============================================="<< endl;
        cout << " " << label << " " << endl;
        cout << "=============================================="<< endl;
    }

    if (format == 2){
        /* matrix is in required format, nothing to do                        */
        return;
    }
/*----------------------------------------------------------------------------*/
    int format_ = format;
    int I, J;
    double V;
    DNS_reducedZeroRows = reduceZeroRows_;
    DNS_transposed = transpose_;
    int n_row_;
    if (DNS_reducedZeroRows ){
        n_row_ = n_row_cmprs;
    }
    else{
        n_row_ = n_row;
    }


    if (l2g_i_coo.size() == 0){
        l2g_i_coo.resize(n_row_cmprs);
        for (int i = 0; i < n_row_cmprs; i++){
            l2g_i_coo[i] = i;
        }
    }


    dense.resize(n_row * n_col);
    for (int i = 0; i < n_row_ * n_col;i++){
        dense[i] = 0;
    }

    if (format == 0){
        for (int i = 0; i < nnz; i++) {
            I = i_coo_cmpr[i];
            J = j_col[i];
            V = val[i];
            dummyFunction(n_row_,I,J,V);
        }
    }
    else if (format == 1){
        for (int i = 0; i < n_row_cmprs; i++) {
            for (int j = i_ptr[i]; j < i_ptr[i + 1]; j++) {
                I = i;
                J = j_col[j];
                V = val[j];
                dummyFunction(n_row_,I,J,V);
            }
        }
    }
    else {
        cout<< "############ "<< label << " ############" << endl;
        cout << "Unexpected format of matrix to be converted to dense" <<
                                                           " format!!!" << "\n";
    }

    if(format_ == 0){
        i_coo_cmpr.clear();
//        i_coo_cmpr.shrink_to_fit();
    }
    else if (format_ == 1){
        i_ptr.clear();
//        i_ptr.shrink_to_fit();
    }

    j_col.clear();
//    j_col.shrink_to_fit();

    val.clear();
//    val.shrink_to_fit();
    numel = n_row_cmprs * n_col;

    format = 2;
}

void Matrix::DNS2CSR(){
    if (format == 1){
/* matrix is in required format, nothing to do                                */
        return;
    }
    DNS2COO();
    COO2CSR();
}

void Matrix::DNS2COO(){

    bool KEEP_ZERO_ENTRIES = true;
    if (format == 0){
        return;
    }
    int _n_row;
    if (DNS_reducedZeroRows){
        _n_row = n_row_cmprs;
    }
    else{
        _n_row = n_row;
    }
    nnz = 0;
    for (int i = 0; i < _n_row * n_col;i++) {
       if (KEEP_ZERO_ENTRIES || dense[i] != 0){
           nnz++;
       }
    }

    if (symmetric){
        nnz = (nnz - n_row_cmprs) * 0.5 + n_row_cmprs;
    }

    i_coo_cmpr.resize(nnz);
    j_col.resize(nnz);
    val.resize(nnz);
    double _val;
    int cnt = 0;
    for (int i = 0; i < _n_row; i++) {
        for (int j = 0; j < n_col; j++) {
            if (DNS_transposed){
                _val = dense[i * n_col + j];
            }
            else{
                _val = dense[i + j * _n_row];
            }
            if (KEEP_ZERO_ENTRIES || _val != 0){
                if (symmetric == 1 && j > i)
                    continue;
                if (symmetric == 2 && i > j)
                    continue;

                i_coo_cmpr[cnt] = i;
                j_col[cnt] = j;
                val[cnt] = _val;
                cnt++;

            }
        }
    }
    dense.clear();
    format = 0;
}

void Matrix::print1dArray(double u[], int nu, string nameOfMat,string folder){

    char _str[256];
    sprintf(_str,"%s/dump_%s.txt",folder.c_str(),nameOfMat.c_str());
    string path2matrix = _str;


    FILE *fp = NULL;
    fp = fopen(path2matrix.c_str(), "w");
    if (fp == NULL) {
        cerr << "fail to open " << path2matrix.c_str() << "\n";
    }

    for (int i = 0 ; i < nu ; i++)
        fprintf(fp, "%.16e\n", u[i]);



    fclose(fp);

}




void Matrix::printToFile(string nameOfMat,string folder, int indOfMat,
                                                        bool _printCooOrDense)
{
    if (l2g_i_coo.size() == 0){
        l2g_i_coo.resize(nnz);
        for (int i = 0; i < n_row_cmprs; i++){
            l2g_i_coo[i] = i;
        }
    }

    const int format_ = format;
    char _str[256];
    sprintf(_str,"%s/dump_%s_%d.txt",folder.c_str(),nameOfMat.c_str(),indOfMat);
    string path2matrix = _str;


    FILE *fp = NULL;
    fp = fopen(path2matrix.c_str(), "w");
    if (fp == NULL) {
        cerr << "fail to open " << path2matrix.c_str() << "\n";
    }

    double V;
    if (_printCooOrDense){
        if (format == 2){
            DNS2COO();
        }
        int I, J,Itmp;
        fprintf(fp, "%%%%MatrixMarket matrix coordinate real %s\n",
        symmetric == 0 ? "general" : "symmetric");
        fprintf(fp, "%d\t%d\t%d\n", n_row, n_col, nnz);
        if (format == 0){
            for (int i = 0; i < nnz; i++) {
                I = l2g_i_coo[i_coo_cmpr[i]] + 1;
                J = j_col[i] + 1;
                V = val[i];
                fprintf(fp, "%d\t%d\t%.16e\n", I,J,V);
            }
        }
        else if (format == 1){
            for (int i = 0; i < n_row_cmprs; i++) {
                for (int j = i_ptr[i]; j < i_ptr[i + 1]; j++) {
                    I = l2g_i_coo[i] + 1;
                    J = j_col[j] + 1;
                    V = val[j];
                    if (DNS_transposed){
                        Itmp = I;
                        J = I;
                        I = Itmp;
                    }
                    fprintf(fp, "%d\t%d\t%.16e\n", I,J,V);
                }
            }
        }
        else {
            cout<< "############ "<< label << " ############" << endl;
            cout << "format = " << format << endl;
            cerr << "Unexpected format of matrix to be printed !!!" << "\n";
        }
        if (format == 0 && format_ == 1){
            COO2CSR();
        }
        if (format == 1 && format_ == 0){
            CSR2COO();
        }
        if (format_ == 2){
            CSRorCOO2DNS(DNS_reducedZeroRows,DNS_transposed);
        }
//        cout << " format=" << format << "format_ =" << format_ << "\n";
    }
    else{
        cout << " print matrix directly from dense is " << endl;
        cout << " not implemented yet " << endl;
    }

    fclose(fp);
}

void Matrix::setZero(){
    if (format == 2){
        for (int i = 0; i < numel; i++){
            dense[i] = 0;
        }
    }
    else{
        for (int i = 0; i < nnz; i++){
            val[i] = 0;
        }
    }
}

void Matrix::submatrix_row_selector(Matrix &A_, vector <int> &v){

    Matrix A;
    A = A_;
    A.CSRorCOO2DNS(true,false);
    zero_dense(v.size(), A.n_col);


    for (int j = 0 ; j < n_col; j++){
        for (int i = 0 ; i < n_row_cmprs; i++){
            dense[i + j * n_row_cmprs ] =  A.dense[ A.g2l_i_coo[ v[i] ] + j * A.n_row_cmprs ];
        }
    }

    A.DNS2CSR();

}


void Matrix::mat_mult_dense(Matrix const &A, string A_NorT, Matrix const &B, string B_NorT){

    int n_k;
    if (A_NorT == "N"){
        n_row_cmprs = A.n_row_cmprs;
        n_k = A.n_col;
    }
    else
    {
        n_row_cmprs = A.n_col;
        n_k = A.n_row_cmprs;
    }
    if (B_NorT == "N")
    {
        n_col = B.n_col;
    }
    else
    {
        n_col = B.n_row_cmprs;
    }

    zero_dense(n_row_cmprs,n_col);


    double A_ik, B_kj;

    for (int i = 0; i < n_row_cmprs; i++){
        for (int j = 0 ; j < n_col; j++){
            for (int k = 0; k < n_k; k++){
                if (A_NorT == "N")
                    A_ik = A.dense[i + k * A.n_row_cmprs];
                else
                    A_ik = A.dense[k + i * A.n_row_cmprs];
                if (B_NorT == "N")
                    B_kj = B.dense[k + j * B.n_row_cmprs];
                else
                    B_kj = B.dense[j + k * B.n_row_cmprs];

                dense[i + j * n_row_cmprs] += A_ik * B_kj;
            }
        }
    }

}


void Matrix::mult(const Matrix& X_in,  Matrix& X_out, bool NorT){

    if (format !=1)
        fprintf(stderr, "%s %d : matrix is not in CSR format .\n",__FILE__,__LINE__);



    int dbg_flag = -1;
    int _n_col_rhs;
    int _n_rws;




    if (NorT)
        _n_rws = n_row_cmprs;
    else
        _n_rws = n_col;


    if (X_in.DNS_transposed){
       _n_col_rhs = X_in.n_row_cmprs;
    }
    else{
        _n_col_rhs = X_in.n_col;
    }



    if (X_out.nnz == 0) {
       X_out.zero_dense(_n_rws,_n_col_rhs);
    }

    mult(&(X_in.dense[0]),&(X_out.dense[0]), NorT, _n_col_rhs, dbg_flag);



}


void Matrix::mult(const double x_in[], double  x_out[], bool NorT, int n_rhs, int dbg_flag){

    for (int i = 0; i < n_row_cmprs; i++) {
        for (int j = i_ptr[i]; j < i_ptr[i + 1]; j++) {
            if (symmetric > 0 ){
                for (int k = 0; k < n_rhs; k++){
                    x_out[i + k * n_row_cmprs] += val[j] * x_in[j_col[j] + k * n_col];
                    if (j_col[j] != i){
                        x_out[j_col[j] + k * n_row_cmprs] += val[j] * x_in[i + k * n_col];
                    }
                }
            }
            else {
                for (int k = 0; k < n_rhs; k++){
                    if (NorT){
                        x_out[i + k * n_row_cmprs] += val[j] * x_in[j_col[j] + k * n_col];
                    }
                    else {
                        x_out[j_col[j] + k * n_col] += val[j] * x_in[i + k * n_row_cmprs];
                    }
                }
            }
        }
    }
}


double Matrix::norm2()
{
    double _norm2 = 0;
    if (format == 2){
        for (unsigned int i = 0 ; i < dense.size(); i++){
            _norm2 += dense[i] * dense[i];
        }
    }
    else{
        for (unsigned int i = 0 ; i < nnz; i++){
            _norm2 += val[i] * val[i];
        }
    }
    return sqrt(_norm2);
}

void Matrix::getBasicMatrixInfo(string str0){
    tmp_label = str0;
    getBasicMatrixInfo();
}


void Matrix::getBasicMatrixInfo(){

    if (tmp_label.empty())
        printf("< %s >              \n", label.c_str());
    else
        printf("< %s >              \n", tmp_label.c_str());
    printf("symmetric           %d (0: unsym, 1: sym. lower tr., 2: sym. upper tr.)\n", symmetric);
    printf("format              %d (0: coo, 1: csr, 2: dense)\n", format);
    printf("nnz:                %d \n", nnz);
    printf("numel:              %d \n", numel);
    printf("n_row_cmprs:        %d \n", n_row_cmprs);
    printf("n_row:              %d \n", n_row);
    printf("n_col:              %d \n", n_col);
    printf("i_ptr.size():       %lu \n", i_ptr.size());
    printf("i_coo_cmpr.size():  %lu \n", i_coo_cmpr.size());
    printf("j_col.size():       %lu \n", j_col.size());
    printf("val.size():         %lu \n", val.size());
    printf("dense.size():       %lu \n", dense.size());
    if (solver > -1)
        printf("solver:             %d (-1: None, 0: pardiso, 1: dissection)\n", solver);
    else
        printf("            :       %s \n", options2["linear_solver"].c_str());


    double density = 100 * nnz / (double)( n_row_cmprs * n_col );

    printf("density:            %1.2f %% \n", density);
}

int Matrix::ij(int i_, int j_){
   return i_ + j_ * n_row_cmprs;
}


void Matrix::getNullPivots(vector < int > & null_pivots_){


  int n_col_ = n_col;
  int n_row_ = n_row_cmprs;

  vector <double> N(dense);
  vector <double>::iterator  it;
  int I,colInd,rowInd;
  double *tmpV = new double[n_row_];
  int *_nul_piv = new int [n_row_];
  double pivot;
  int tmp_int;
  for (int i = 0;i < n_row_;i++){
      _nul_piv[i]=i;
  }
  for (int j=0;j<n_col_;j++){
    it = max_element(N.begin(),N.end()-j*n_row_, compareDouble);
    I = it - N.begin();
    colInd = I/n_row_;
    rowInd = I-colInd*n_row_;
    for (int k=0;k<n_col_-j;k++){
      tmpV[k] = N[ij(n_row_-1-j,k)];
      N[ij(n_row_-1-j,k)] = N[ij(rowInd,k)];
      N[ij(rowInd,k)]= tmpV[k];
    }
    tmp_int = _nul_piv[rowInd];
    _nul_piv[rowInd] = _nul_piv[n_row_-1-j];
    _nul_piv[n_row_-1-j] = tmp_int;
    if (n_col_ - 1 - j != colInd) {
        memcpy( tmpV, &(N[ij(0,n_col_-1-j)]) , sizeof( double ) * n_row_);
        memcpy( &(N[ij(0, n_col_- 1 - j)]), &(N[ij(0, colInd)]), sizeof( double ) * n_row_);
        memcpy( &(N[ij(0,colInd)]),tmpV , sizeof( double ) * n_row_);
    }
    pivot = N[ij(n_row_-1-j,n_col_-1-j)];
    for (int J=0;J<n_col_-j-1;J++){
      for (int I=0;I<n_row_-j;I++){
        N[ij(I,J)] -= N[ij(I,n_col_-1-j)]*N[ij(n_row_-1-j,J)]/pivot;
      }
    }
  }
  null_pivots_.resize(n_col_);
  for (int i = 0;i<n_col_;i++){
    null_pivots_[i] = _nul_piv[n_row_-1-i];
  }
  sort(null_pivots_.begin(),null_pivots_.end());
  delete [] _nul_piv;
  delete [] tmpV;
//
}

void Matrix::sym_factor_dissection(){

#ifdef DISSECTION
      //diss_dslv = new uint64_t;
      diss_num_threads = 1;
      diss_dslv = new uint64_t;
      diss_init(*diss_dslv, order_number, 1, 0, diss_num_threads, 1);
      int diss_sym = 1; // isSym = 1 + upper_flag = 0
      int diss_decomposer = 0;

      diss_s_fact(*diss_dslv,
          n_row, &i_ptr[0], &j_col[0],
          diss_sym, diss_decomposer);
//      printCooOrDense = true;
//      ker_Ac.printToFile("ker_Ac",folder,0,printCooOrDense);

#endif


}

void Matrix::num_factor_dissection(){
#ifdef DISSECTION
    diss_indefinite_flag = 1;
    if (diss_scaling < 0){
        diss_scaling = 1;
    }
    diss_eps_pivot = 1.0e-2;
    diss_n_fact(*diss_dslv, &val[0], diss_scaling, diss_eps_pivot,
      diss_indefinite_flag);
#endif
}



void Matrix::num_factor_dissection(Matrix &R){
//


#ifdef DISSECTION
    num_factor_dissection();
    int n0;
    diss_get_kern_dim(*diss_dslv, &n0);
    R.zero_dense(n_row,n0);
    diss_get_kern_vecs(*diss_dslv, &R.dense[0]);

#endif

}

void Matrix::solve_system_dissection(Matrix const &B, Matrix &X){

#ifdef DISSECTION
    int _n_row,_n_col;
    if (B.DNS_transposed){
        _n_row = B.n_col;
        _n_col = B.n_row_cmprs;
    }
    else
    {
        _n_row = B.n_row_cmprs;
        _n_col = B.n_col;
    }

    X.zero_dense(_n_row, _n_col);
    X.dense = B.dense;
//    cout << "B.n_row = " << B.n_row <<", ";
//    cout << "B.n_row_cmprs = " << B.n_row_cmprs <<", ";
//    cout << "B.n_col = " << B.n_col <<"\n";

    int projection = 0;
    int trans = 0;

    if (_n_col == 1){
        diss_solve_1(*diss_dslv,&X.dense[0],
                      projection, trans);
    }
    else {
        diss_solve_n(*diss_dslv,&X.dense[0],
                     _n_col, projection, trans);
    }
#endif
}

void Matrix::sym_factor(string linear_solver_){

    linear_solver = linear_solver_;

    if (linear_solver.compare("pardiso") == 0){
        sym_factor_pardiso();
    }
    else if (linear_solver.compare("dissection") == 0){
        sym_factor_dissection();
    }
    else{
        cout << "ERROR >>>>>>>>>>>>>>>>>>>>>>>>>>> "  << endl;
        return;
    }
    solver = -1;
}
void Matrix::sym_factor(int pardiso_0_dissection_1){
    solver = pardiso_0_dissection_1;
    if (solver == 0){
        sym_factor_pardiso();
    }
    else if (solver == 1){
        sym_factor_dissection();
    }
    else{
        cout << "ERROR >>>>>>>>>>>>>>>>>>>>>>>>>>> "  << endl;
        return;
    }
}


void Matrix::num_factor(){

    if (solver > -1){
        if (solver == 0){
            num_factor_pardiso();
        }
        else if (solver == 1){
            num_factor_dissection();
        }
        else{
            cout << "ERROR >>>>>>>>>>>>>>>>>>>>>>>>>>> "  << endl;
            return;
        }
    }
    else {
        if (linear_solver.compare("pardiso") == 0){
            num_factor_pardiso();
        }
        else if (linear_solver.compare("dissection") == 0){
            num_factor_dissection();
        }
        else{
            cout << "ERROR >>>>>>>>>>>>>>>>>>>>>>>>>>> "  << endl;
            return;
        }
    }
}

void Matrix::num_factor(Matrix &R, bool checkOrthogonality_){

    // use R to regularize matrix then call factorization
    if (R.numel != 0){
        R.getNullPivots(nullPivots);
        int j, cnt = 0;
    //    cout << "null_pivot: ";
        vector <int > tmp;
        tmp.resize(nullPivots.size());
        for (int i = 0; i < n_row_cmprs; i++) {
            j = i_ptr[i];
            if (i == nullPivots[cnt]){
    //            cout << i << ", ";
                val[j] *= 2;
                tmp[cnt] = j;
                cnt++;
            }
        }

        num_factor();

        bool printCooOrDense = true;
        if (!printed && atoi(options2["print_matrices"].c_str()) > 1){
            printToFile("K_reg",options2["path2data"],order_number,printCooOrDense);
            printed = true;
        }

        for (int i = 0 ; i < tmp.size();i++)
            val[tmp[i]] *= 0.5;


    }
    else{
        if (solver > -1){
            if (solver == 0){
                num_factor_pardiso();
            }
            else if (solver == 1){
                num_factor_dissection(R);
            }
            else {
                cout << "ERROR >>>>>>>>>>>>>>>>>>>>>>>>>>> "  << endl;
                return;
            }
        }
        else{
            if (linear_solver.compare("pardiso") == 0){
                num_factor_pardiso();
            }
            else if (linear_solver.compare("dissection") == 0){
                num_factor_dissection(R);
            }
            else {
                cout << "ERROR >>>>>>>>>>>>>>>>>>>>>>>>>>> "  << endl;
                return;
            }
        }
    }

    if (R.numel !=0 && checkOrthogonality_){

        Matrix Y;
        Y.zero_dense(n_row_cmprs , R.n_col);
        mult(R,Y,true);

        double normK = norm2();
        double normR = R.norm2();
        double normKR = Y.norm2();

        printf("=========    %s    =========\n",R.label.c_str());
        printf("|A| = %.3e, ", normK);
        printf("|N| = %.3e, ", normR);
        printf("|A*N|/(|A|*|N|) = %.3e \n", normKR / (normK*normR));

    }

}

void Matrix::num_factor(Matrix &R){
    bool checkOrthogonality_ = false;
    num_factor(R,checkOrthogonality_);
}




void Matrix::sym_factor_pardiso(){


#ifdef USE_PARDISO
    if (msglvl<0)
        msglvl = 0;
    mtype = -2;              /* Real symmetric matrix */
    n1_p = n_row_cmprs;
    MKL_INT i;
    nrhs = 1;               /* Number of right hand sides. */
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
    iparm[9] = 8;        /* Perturb the pivot elements with 1E-13 */
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
    iparm[34] = 1;        /* One- or zero-based indexing of columns and rows.*/
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Which factorization to use. */
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
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n1_p, &val[0], &i_ptr[0], &j_col[0],
            &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during symbolic factorization: %d", error);
        exit (1);
    }


#endif


}


void Matrix::num_factor_pardiso()
{



#ifdef USE_PARDISO



    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n1_p, &val[0], &i_ptr[0], &j_col[0],
            &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
#endif

}


void Matrix::solve_system(Matrix &B, Matrix & X){

    if (solver > -1){
        if (solver == 0){
            solve_system_pardiso(B,X);
        }
        else if(solver == 1){
            solve_system_dissection(B,X);
        }
        else{
            cout << "ERROR >>>>>>>>>>>>>>>>>>>>>>>>>>> "  << endl;
            return;
        }
    }
    else{
        if (linear_solver.compare("pardiso") == 0){
            solve_system_pardiso(B,X);
        }
        else if(linear_solver.compare("dissection") == 0){
            solve_system_dissection(B,X);
        }
        else{
            cout << "ERROR >>>>>>>>>>>>>>>>>>>>>>>>>>> "  << endl;
            return;
        }

    }

}



void Matrix::solve_system_pardiso(Matrix& B, Matrix& X){
#ifdef USE_PARDISO

    /* -------------------------------------------------------------------- */
    /* .. Back substitution and iterative refinement. */
    /* -------------------------------------------------------------------- */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
    MKL_INT idum;         /* Integer dummy. */
    int _n_row, _n_col;
    if (B.DNS_transposed){
        cout << "T" <<endl;
        _n_row = B.n_col;
        _n_col = B.n_row_cmprs;
    }
    else{
        _n_row = B.n_row_cmprs;
        _n_col = B.n_col;
    }
    nrhs = _n_col;

//    if (X.label == "Kplus_K"){
//        X.getBasicMatrixInfo();
//    }
//    else{
//        cout << "\n\n\n========================= \n\n\n";
//
//    }

    if (X.numel == 0 ){
        X.zero_dense(_n_row,_n_col);
    }
    else
    {
        X.setZero();
    }


    phase = 33;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n1_p, &val[0], &i_ptr[0], &j_col[0], &idum, &nrhs, iparm, &msglvl,
            &B.dense[0], &X.dense[0], &error);

    if ( error != 0 )
    {
        printf ("\nERROR during solution: %d", error);
        exit (3);
    }
//    printf ("\nSolve completed ... ");
#endif

}

void Matrix::FinalizeSolve(int i_)
{
    if (solver < -1){
        if (solver == 0){
#ifdef USE_PARDISO
            double ddum;          /* Double dummy */
            MKL_INT idum;         /* Integer dummy. */
            phase = -1;           /* Release internal memory. */
            PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                     &n, &ddum, &i_ptr[0], &j_col[0], &idum, &nrhs,
                    iparm, &msglvl, &ddum, &ddum, &error);
#endif
        }
        else if (solver == 1)
        {
#ifdef DISSECTION
            diss_free(*diss_dslv);
#endif
        }
    }
    else{
        if (linear_solver.compare("pardiso") == 0 ){
#ifdef USE_PARDISO
            double ddum;          /* Double dummy */
            MKL_INT idum;         /* Integer dummy. */
            phase = -1;           /* Release internal memory. */
            PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                     &n, &ddum, &i_ptr[0], &j_col[0], &idum, &nrhs,
                    iparm, &msglvl, &ddum, &ddum, &error);
#endif
        }
        else if (linear_solver.compare("dissection") == 0)
        {
#ifdef DISSECTION
            diss_free(*diss_dslv);
#endif
        }
    }
}

void Matrix::testSolver(string folder, int pardiso_0_dissection_1){
#if 0
/*------------------------------ TEST ------------------------------*/
    Matrix A, B, U;
    int n = 11;
    int cnt = 0;
    int nrhs = 5;
    A.nnz = 2 * n - 1;
    A.n_row = n;
    A.n_row_cmprs = n;
    A.n_col = n;

    B.zero_dense(n , nrhs);
    double h = 1. / (n + 1);

    for (int i = 0; i < n; i ++){
        for (int j = 0; j < nrhs; j ++){
            B.dense[ nrhs * i + j] =  h * h;
        }
        for (int j = 0; j < n; j ++){
            if (i == j){
                A.i_coo_cmpr.push_back(i);
                A.j_col.push_back(j);
                A.val.push_back(2);
                cnt++;
            }
            else if ((j-i) == 1){
                A.i_coo_cmpr.push_back(i);
                A.j_col.push_back(j);
                A.val.push_back(-1);
                cnt++;
            }
        }
    }
    A.format = 0;
    A.symmetric = 2;

//    A.printToFile("A",folder,1000,true);
    A.COO2CSR();

    cout << "A_III  " << A.i_ptr.size() << endl;
    cout << "A_III  " << A.val.size() << endl;
    cout << "A_III  " << A.j_col.size() << endl;



    A.sym_factor(pardiso_0_dissection_1);
    A.num_factor();
    A.solve_system(B, U);
    double x_current = h;
    double u_current, delta = 0;
    for (int i = 0; i < n; i++){
        u_current = 0.5 * x_current * (1. - x_current);
        delta += (u_current - U.dense[i]) * (u_current - U.dense[i]);
        if (n < 50){
//            cout << u_current << " " << U.dense[i] <<"\n";
            printf("%3.3e %3.3e    del = %1.16e \n",
                   u_current, U.dense[i],fabs( u_current - U.dense[i] )  );
        }
        x_current += h;

    }
    cout << "|| exact - numerical solution || = " << sqrt(delta) << "\n";

    U.DNS2CSR();
    U.printToFile("X",folder,1000,true);
    B.DNS2CSR();
    B.printToFile("B",folder,1000,true);
    A.FinalizeSolve(0);
/*------------------------------ TEST ------------------------------*/
#endif
}

void Matrix::updateCOOstructure(vector <TRIPLET >& vec_, Matrix &denseMat,int i_shift, int j_shift){


    TRIPLET i_triplet;
    for (int j = 0; j < denseMat.n_col;j++){
        for (int i = 0; i < denseMat.n_row_cmprs;i++){
            i_triplet.I = i_shift + i;
            i_triplet.J = j_shift + j;
            i_triplet.V = denseMat.dense[i + j * denseMat.n_row_cmprs ];
            if (i_triplet.J >= i_triplet.I && i_triplet.V != 0){
                vec_.push_back(i_triplet);
            }
        }
    }
}

void Matrix::print_int_vector(vector <double> &dS, int print_first_n, int print_last_n){

    if (print_first_n >0 || print_last_n > 0){
        cout <<"Print ";
    }

    if (print_first_n > 0){
        cout << "first " << print_first_n;
    }

    if (print_first_n >0 && print_last_n > 0){
       cout <<" and " ;
    }

    if (print_last_n > 0){
        cout << "last " << print_last_n;
    }


    if (print_first_n >0 || print_last_n > 0){
        cout <<" entries";
    }

    cout << endl;



    int nPrint = print_first_n > dS.size() ? dS.size(): print_first_n;
    for (int i = 0; i < nPrint;i++){
        cout << "d("<< i << ")=  " << dS[i] << endl;
    }
    nPrint = print_last_n > dS.size() ? dS.size() : print_last_n;
    for (int i = dS.size() - nPrint; i < dS.size() ;i++){
        cout << "d("<< i << ")=  " << dS[i] << endl;
    }
    cout << "-- end --" << endl;
}


void Matrix::getEigVal_DNS(Matrix &A, Matrix &S, int print_first_n){
    Matrix::getEigVal_DNS(A, S, print_first_n, 0);
}

void Matrix::getEigVal_DNS(Matrix A_, Matrix &S, int print_first_n, int print_last_n){
    map <string, string> c_options2 = A_.options2;

    if (A_.n_row_cmprs > 3000){
        cout << "Matrix is too big to get SVD" << endl;
        return;
    }
    if (A_.format != 2)
        A_.CSRorCOO2DNS(false,false);

    Matrix A;
    int A_dim = int(1 + A_.n_row_cmprs) * A_.n_row_cmprs / 2;
    A.zero_dense(1, A_dim);

    int cnt = 0;
    for (int j = 0; j < A_.n_col;j++){
        for (int i = 0; i <= j; i++) {
            A.dense[cnt] = A_.dense[i + j * A_.n_row_cmprs];
            cnt++;
        }
    }
    S.zero_dense(A_.n_row_cmprs,1);
    char JOBZ_ = 'V';
    char UPLO_ = 'U';
    double *ZK_modif = new double[A_.n_row_cmprs * A_.n_row_cmprs];
    MKL_INT info;
    MKL_INT ldzA = A_.n_row_cmprs;
    info = LAPACKE_dspev (LAPACK_COL_MAJOR, JOBZ_, UPLO_,
            A_.n_row_cmprs, &(A.dense[0]), &(S.dense[0]), ZK_modif, ldzA);
    S.printToFile(A_.label,c_options2["path2data"],444,true);
    delete [] ZK_modif;

    for (int i = 0; i < S.n_row_cmprs; i++){
        S.dense[i] =  fabs(S.dense[i]);
    }

    sort(S.dense.begin(),S.dense.end());


    cout << "\n\n#######  Sorted Abs val of Eigenvalues of " << A_.label << " #######\n";
    Matrix::print_int_vector(S.dense,print_first_n,print_last_n);
}

void Matrix::getSingularVal_DNS(Matrix &A, Matrix &S, int print_first_n){
    Matrix::getSingularVal_DNS(A, S, print_first_n, 0);
}

void Matrix::getSingularVal_DNS(Matrix A, Matrix &S, int print_first_n, int print_last_n){
    map <string, string> c_options2 = A.options2;

    if (A.n_row_cmprs > 2000){
        cout << "Matrix is too big to get SVD" << endl;
        return;
    }
    if (A.format != 2)
        A.CSRorCOO2DNS(false,false);
    S.zero_dense(A.n_row_cmprs,1);
    double *U     = new double[A.n_row_cmprs * A.n_row_cmprs];
    double *Vt    = new double[A.n_col * A.n_col];
    double *superb  = new double[A.n_col-1];
    MKL_INT info;
    MKL_INT lds = A.n_row_cmprs;
    MKL_INT Srows = A.n_row_cmprs;
    MKL_INT Scols= A.n_col;
    info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'A', Scols, Srows, &(A.dense[0]), lds,
                          &(S.dense)[0], U, lds, Vt, lds, superb );

    if (print_first_n != 0 && print_last_n != 0){
        cout << "\n\n#######  Singular values of " << A.label << " #######\n";
        Matrix::print_int_vector(S.dense,print_first_n,print_last_n);
        S.printToFile(A.label,c_options2["path2data"],555,true);
    }
    delete [] U ; delete [] Vt; delete [] superb;
}


bool Matrix::test_of_Bc_constraints(Matrix &A){


  Matrix ATA,S_S;
  ATA.mat_mult_dense(A,"N",A,"T");

  double jump_in_sing_values_alerting_singularity = 1e-4;
  int  sc_size = ATA.n_row_cmprs;
  int ind_U_V;
  int defect_K_in = 0;

  Matrix::getSingularVal_DNS(ATA, S_S, 0);
  double ratio;

  for (int i = 0; i < sc_size-1;i++){
    ratio = fabs(S_S.dense[i+1]/S_S.dense[i]);
    if (ratio < jump_in_sing_values_alerting_singularity){
      ind_U_V = i+1;
      defect_K_in=sc_size-(i+1);
      break;
    }
  }

//  cout << " defect BcTBc = " << defect_K_in;
//
//  cout << "\t\t";
//  for (int i = 0; i < sc_size; i++){
//      cout << S_S.dense[i] << ", " ;
//  }
//  cout << endl;

  return defect_K_in == 0;

}



double Matrix::dot(Matrix const &x, Matrix const &y){
    double xy = 0;
    for (int i = 0; i < x.n_row_cmprs; i++)
        xy += x.dense[i] * y.dense[i];
    return xy;
}


double Matrix::dot(double *x, double *y, int n){
  double dot_xy = 0.0;
  for (int i = 0; i< n; i++){
    dot_xy+=x[i]*y[i];
  }
  return dot_xy;
}



void Matrix::add(Matrix &v, double a){

    for (int i = 0; i < n_row_cmprs; i++)
        dense[i] += v.dense[i] * a;

}



void Matrix::test_K_Kp_K_condition(Matrix &Ks){

    bool use_Ks = true;

    if (solver == -1){
        fprintf(stderr, "%s %d :  matrix is not factorized... .\n", __FILE__, __LINE__);
        return;
    }

    Matrix K_dense,Kplus_K,K_Kplus_K;

    K_dense.label = "K_dense";
    K_dense.zero_dense(n_row_cmprs,n_col);


    if (use_Ks){
        Ks.CSRorCOO2DNS(false,false);
    }
    else{
        CSRorCOO2DNS(false,false);
    }

    for (int i = 0; i < K_dense.numel;i++){
        if (use_Ks){
            K_dense.dense[i] = Ks.dense[i];
        }
        else{
            K_dense.dense[i] = dense[i];
        }
    }

    if (!use_Ks){
        DNS2CSR();
    }



//    K_dense.printToFile("K_dense",options.path2data,order_number,true);

    Kplus_K.label = "Kplus_K";
    solve_system(K_dense,Kplus_K);


//    if (label == "stiffness matrix"){
//
//        Matrix ONES, KplusONES;
//        int d = order_number;
//        bool printCooOrDense = true;
//        string folder = "../data/";
//        ONES.zero_dense(n_row_cmprs,1);
//        for (int i =0; i < n_row_cmprs; i++)
//            ONES.dense[i] = 1;
//        solve_system(ONES,KplusONES);
//        ONES.printToFile("ONES",folder,d,printCooOrDense);
//        KplusONES.printToFile("KplusONES",folder,d,printCooOrDense);
//
//        if (d == 0){
//            for (int i = 0; i < n_row_cmprs; i++)
//                cout << KplusONES.dense[i] << endl;
//
//        }
//    }






    if (order_number == 0){
        cout << "is " << label  <<" transposed ??? " << K_dense.DNS_transposed << endl;
    }
//
//    Kplus_K.printToFile("Kplus_K",options.path2data,order_number,true);
//
    mult(Kplus_K,K_Kplus_K,true);
//    K_Kplus_K.mat_mult_dense(K_dense,"N",Kplus_K,"N");
//    K_Kplus_K.label = "K_Kplus_K";
//
//    K_Kplus_K.printToFile("K_Kplus_K",options.path2data,order_number,true);

    double norm_K = K_dense.norm2();
    double norm_Kplus_K = Kplus_K.norm2();
    double norm_K_Kplus_K = K_Kplus_K.norm2();

    for (int i = 0; i < K_dense.numel; i++){
        K_Kplus_K.dense[i] -= K_dense.dense[i];
    }



    double norm_K_minus_K_Kplus_K = K_Kplus_K.norm2();


    cout << " =============================================  \n";
    cout << " A = " << label << " \n";
    printf(" || A ||                          = %3.15e\n", norm_K);
    printf(" || Aplus *_A ||                  = %3.15e\n", norm_Kplus_K);
    printf(" || A * Aplus * A ||              = %3.15e\n", norm_K_Kplus_K);
    printf(" || A - A * Aplus * A || / || A|| = %3.15e\n", norm_K_minus_K_Kplus_K / norm_K);
    cout << " =============================================  \n";

}


int compareVecVec(const vector <int >& a, const vector <int>& b)  { return a[0] < b[0]; }

void Matrix::createDirichletPreconditioner(const Matrix &Bf,
                                           const Matrix & K,
                                           Matrix &Krr_,
                                           Matrix &Krs_,
                                           Matrix &Kss_,
                                           Matrix &Precond){


    map <string, string> c_options2 = K.options2;

    vector<int>::iterator it;
    vector < int > perm_vec = Bf.j_col;


    sort(perm_vec.begin(), perm_vec.end());
    it = unique (perm_vec.begin(), perm_vec.end());
    perm_vec.resize( distance(perm_vec.begin(),it));



    vector <int> perm_vec_full (  K.n_row_cmprs);
    vector <int> perm_vec_diff (  K.n_row_cmprs);
    vector <int> I_row_indices_p (K.nnz);
    vector <int> J_col_indices_p (K.nnz);

    vector < TRIPLET > triplet(K.nnz);




    for (int i = 0; i < perm_vec_full.size(); i++) {
        perm_vec_full[i] = i;
    }

    it = std::set_difference( perm_vec_full.begin(), perm_vec_full.end(), perm_vec.begin(), perm_vec.end(), perm_vec_diff.begin() );
    perm_vec_diff.resize(it - perm_vec_diff.begin());

    perm_vec_full = perm_vec_diff;
    perm_vec_full.insert(perm_vec_full.end(), perm_vec.begin(), perm_vec.end());

    Matrix K_modif = K;
    K_modif.label = "K_modif";

    vector <vector <int > > vec_I1_i2(K_modif.n_row_cmprs, vector <int >(2, 1));

    for (int i = 0; i < K_modif.n_row_cmprs;i++){
        vec_I1_i2[i][0] = perm_vec_full[i];
        vec_I1_i2[i][1] = i; // position to create reverse permutation
    }

    sort(vec_I1_i2.begin(), vec_I1_i2.end(), compareVecVec);

    // permutations made on matrix in COO format
    K_modif.CSR2COO();

//     K_modif.printToFile("K_modif_v1",c_options2["path2data"].c_str(),order_number,printCooOrDense_);


    int I_index,J_index;
    bool unsymmetric = false;
    for (int i = 0;i<K_modif.nnz;i++){
       I_index = vec_I1_i2[K_modif.i_coo_cmpr[i]][1];
       J_index = vec_I1_i2[K_modif.j_col[i]][1];
       if (unsymmetric || I_index<=J_index){
         I_row_indices_p[i]=I_index;
         J_col_indices_p[i]=J_index;
         triplet[i].I =  I_index;
         triplet[i].J =  J_index;
       }
       else{
         I_row_indices_p[i]=J_index;
         J_col_indices_p[i]=I_index;
         triplet[i].I =  J_index;
         triplet[i].J =  I_index;
       }
       triplet[i].V =  K_modif.val[i];
     }



    for (int i = 0; i < K_modif.nnz; i++){
       K_modif.i_coo_cmpr[i] = triplet[i].I;
       K_modif.j_col[i]      = triplet[i].J;
    }

//     K_modif.printToFile("K_modif_v2",c_options2["path2data"],order_number,printCooOrDense_);


//
    K_modif.sortAndUniqueCOO(triplet);

//     K_modif.printToFile("K_modif_v3",c_options2["path2data"],order_number,printCooOrDense_);


    K_modif.COO2CSR();


    int sc_size = perm_vec.size();


    Matrix K_rr("K_rr");
    Matrix K_rs("K_rs");
    Matrix K_ss("K_ss");
    Matrix &KsrInvKrrKrs = Precond;
    KsrInvKrrKrs.label = "KsrInvKrrKrs";

//    Matrix InvKrrKrs("InvKrrKrs");

    int i_start = 0;
    int nonsing_size = K_modif.n_row - sc_size - i_start;
    int j_start = nonsing_size;

    K_ss.getSubDiagBlockmatrix(K_modif,K_ss,nonsing_size,sc_size);
    K_rs.getSubBlockmatrix_rs(K_modif,K_rs,i_start, nonsing_size,j_start,sc_size);



    K_rr.getSubDiagBlockmatrix(K_modif,K_rr,i_start, nonsing_size);

    if (atoi(c_options2["print_matrices"].c_str()) > 3){
        bool printCooOrDense = true;
        K_rs.printToFile("Krs",c_options2["path2data"],order_number,printCooOrDense);
        K_ss.printToFile("Kss",c_options2["path2data"],order_number,printCooOrDense);
        K_rr.printToFile("Krr",c_options2["path2data"],order_number,printCooOrDense);
    }



    K_rr.options2 = c_options2;
    K_rr.sym_factor("pardiso");
    Matrix R;
    bool checkOrthogonality_ = false;
    K_rr.num_factor(R,checkOrthogonality_);

    Matrix K_rs_copy = K_rs;
    K_rs_copy.label = "K_rs_copy";
    K_rs.label = "o";


    if (c_options2.at("preconditioner").compare("Dirichlet_implicit") == 0 ){
        Krr_ = K_rr;
        Kss_ = K_ss;
        Krs_ = K_rs;

    }
    else if (c_options2.at("preconditioner").compare("Dirichlet_explicit") == 0 ){
        Matrix InvKrrKrs;
        K_rs.CSRorCOO2DNS(false,false);
        K_rr.solve_system(K_rs,InvKrrKrs);
        K_rs_copy.mult(InvKrrKrs,KsrInvKrrKrs,false);

        for (int i = 0; i < KsrInvKrrKrs.numel;i++)
            KsrInvKrrKrs.dense[i] *= -1;

        for (int i = 0; i < n_row_cmprs; i++) {
            for (int j = K_ss.i_ptr[i]; j < K_ss.i_ptr[i + 1]; j++) {
                Precond.dense[i + K_ss.j_col[j] * K_ss.n_row_cmprs ] += K_ss.val[j];
                if (K_ss.j_col[j] != i){
                    Precond.dense[K_ss.j_col[j] + i * K_ss.n_row_cmprs ] += K_ss.val[j];
                }
            }
        }
    }
}




void Matrix::getSubBlockmatrix_rs( Matrix & A_in, Matrix & A_out,
                                         int  i_start, int i_size,
                                         int  j_start, int j_size){


//
// Original matrix A_in is assembled from 4 submatrices
//
//      A_in = [A_in(r,r)  A_in(r,s)]
//             [A_in(s,r)  A_in(s,s)].
//
// Function 'getSubBlockmatrix_rs' returns square matrix A_in(r,s) in CSR format.
//
// rev. 2015-10-10 (A.M.)
//
// step 1: getting nnz of submatrix

  int nnz_new=0;
  int offset = A_in.i_ptr[0] ? 1 : 0;
  for (int i = 0;i<i_size;i++){
    for (int j = A_in.i_ptr[i+i_start];j<A_in.i_ptr[i+i_start+1];j++){
      if ((A_in.j_col[j-offset]-offset)>=j_start &&
                      (A_in.j_col[j-offset]-offset)<(j_start+j_size)){
        nnz_new++;
      }
    }
  }

// step 2: allocation 1d arrays
  A_out.val.resize(nnz_new);
  A_out.j_col.resize(nnz_new);
  A_out.i_ptr.resize(i_size+1);
  A_out.n_row =i_size;
  A_out.n_row_cmprs =i_size;
  A_out.n_col =j_size;
  A_out.nnz=nnz_new;
  A_out.numel=nnz_new;
  A_out.format = 1;
  A_out.symmetric = 0;
// step 3: filling 1d arrays
  int ijcnt=0;
  A_out.i_ptr[0]=offset;
  for (int i = 0;i<i_size;i++){
    for (int j = A_in.i_ptr[i+i_start];j<A_in.i_ptr[i+i_start+1];j++){
      if ((A_in.j_col[j-offset]-offset)>=j_start && (A_in.j_col[j-offset]-offset)<(j_start+j_size)){
        A_out.j_col[ijcnt] = (A_in.j_col[j-offset]) - j_start;
        A_out.val[ijcnt]=A_in.val[j-offset];
        ijcnt++;
      }
    }
    A_out.i_ptr[i+1]=offset+ijcnt;
  }
}


void Matrix::getSubDiagBlockmatrix(const Matrix & A_in, Matrix & A_out, int i_start, int size_rr){
//
// Function 'getSubDiagBlockmatrix' returns the diagonal block A_in(r,r) from original A_in,
// where r = { i_start , i_start+1 , i_start+2 , ... , istart + size_rr - 1 }
//
// rev. 2015-10-10 (A.M.)
//
// step 1: getting nnz of submatrix
  int nnz_new=0;
  int offset = A_in.i_ptr[0] ? 1 : 0;
  for (int i = 0;i<size_rr;i++){
    for (int j = A_in.i_ptr[i+i_start];j<A_in.i_ptr[i+i_start+1];j++){
      if ((A_in.j_col[j-offset]-offset)>=i_start &&
                      (A_in.j_col[j-offset]-offset)<(i_start+size_rr)){
        nnz_new++;
      }
    }
  }
// step 2: allocation 1d arrays
  A_out.val.resize(nnz_new);
  A_out.j_col.resize(nnz_new);
  A_out.i_ptr.resize(size_rr+1);
  A_out.n_row = size_rr;
  A_out.n_row_cmprs = size_rr;
  A_out.n_col =size_rr;
  A_out.nnz=nnz_new;
  A_out.symmetric = A_in.symmetric;
  A_out.format = A_in.format;
// step 3: filling 1d arrays
  int ijcnt=0;
  A_out.i_ptr[0]=offset;
  for (int i = 0;i<size_rr;i++){
    for (int j = A_in.i_ptr[i+i_start];j<A_in.i_ptr[i+i_start+1];j++){
      if ((A_in.j_col[j-offset]-offset)>=i_start &&
                    (A_in.j_col[j-offset]-offset)<(i_start+size_rr)){
        A_out.j_col[ijcnt] = (A_in.j_col[j-offset]) - i_start;
        A_out.val[ijcnt]=A_in.val[j-offset];
        ijcnt++;
      }
    }
    A_out.i_ptr[i+1]=offset+ijcnt;
  }
}


void Matrix::get_kernel_from_K(Matrix &K, Matrix &Kplus_R){
    int dbg = -1;
    get_kernel_from_K(K, Kplus_R,dbg);
}


void Matrix::get_kernel_from_K(Matrix &K, Matrix &Kplus_R, int dbg){

//
// Routine calculates kernel Kplus_R of K satisfied euqality K * Kplus_R = O,
// where O is zero matrix, and it makes the matrix K non-singular (K_reg)
// utilizing spectral conditions of Schur complement. Then ||K-K*inv(K_reg)*K||=0.0
//
//
// rev. 2016-02-03 (A.M.)
//==============================================================================
//

typedef int  eslocal;
#define SEQ_VECTOR vector
//
//    1) diagonalScaling
//  reducing of big jump coefficient effect (TODO include diagonal scaling into whole ESPRESO)
//BOOL DIAGONALSCALING                                  = true;
  bool diagonalScaling                                  = true;

//    2) permutVectorActive
//  random selection of singular DOFs
// 0 - no permut., 1 - std::vector shuffle, 2 - generating own random sequence -
//ESLOCAL PERMUTVECTORACTIVE                            = 1;
  eslocal permutVectorActive                            = 1;

//    3) use_null_pivots_or_s_set
  // NtN_Mat from null pivots or fixing DOFs
//BOOL USE_NULL_PIVOTS_OR_S_SET                         = TRUE;
  bool use_null_pivots_or_s_set                         = true;

//    4) diagonalRegularization
//  regularization only on diagonal elements (big advantage: patern of K and K_regular is the same !!!)
//  size of set 's' = defect(K)
//  It's is active, only if and only if 'use_null_pivots_or_s_set = true'
//BOOL DIAGONALREGULARIZATION                           = TRUE;
  bool diagonalRegularization                           = true;


//    6) get_n_first_and_n_last_eigenvals_from_dense_S
// get and print 2*n S eigenvalues
//ESLOCAL GET_N_FIRST_AND_N_LAST_EIGENVALS_FROM_DENSE_S = 0;
  eslocal get_n_first_and_n_last_eigenvals_from_dense_S = 0;


//    8) fixing_nodes_or_dof
// non-singular part determined by fixing nodes (FN),
// min(fixing_nodes_or_dof)>=3; if variable is nonzero,
// parameter sc_size is set to fixing_nodes_or_dof*dofPerNode
//ESLOCAL FIXING_NODES_OR_DOF                           = 0;
  eslocal fixing_nodes_or_dof = 0;
//ESLOCAL DOFPERNODE                                    = 3;
  eslocal dofPerNode                                    = 3;
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
//

//    2) check_nonsing
// if check_nonsing>0, checking of K_rr non-singularity is activated and it is repeated
// (check_nonsing) times.
//ESLOCAL CHECK_NONSING                                 = 0;
  eslocal check_nonsing                                 = 0;

//    3) max_size_of_dense_matrix_to_get_eigs
// if size of K is less then CHECK_N..., K is converted to dense format to get eigenvalues.
//ESLOCAL MAX_SIZE_OF_DENSE_MATRIX_TO_GET_EIGS          = 2500;
  eslocal max_size_of_dense_matrix_to_get_eigs          = 2500;

//    4) sc_size
// specification of size of Schur complement used for detection of zero eigenvalues.
//eslocal  sc_size >= expected defect 'd' (e.g. in elasticity d=6).
//ESLOCAL SC_SIZE                                       = 50;
  eslocal sc_size                                       = 50;

//    5) twenty
// testing last twenty eigenvalues of S to distinguish, if d-last ones are zero or not.
//ESLOCAL TWENTY                                        = 20;
  eslocal twenty                                        = 20;
  // twenty eigenvalues are ascendly ordered in d = d[0],d[1], ..., d[n-2],d[n-1]

//    6) jump_in_eigenvalues_alerting_singularity
// if d[i]/d[i+1]< jump_in_eigenvalues_alerting_singularity, d[i] is last nonzero eigenvalue
//DOUBLE JUMP_IN_EIGENVALUES_ALERTING_SINGULARITY       = 1.0E-5;
  double jump_in_eigenvalues_alerting_singularity       = 1.0e-5;

// ///////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////




  map <string, string> c_options2 = K.options2;



  if (!use_null_pivots_or_s_set) diagonalRegularization=false;


  //TODO if K.rows<=sc_size, use directly input K instead of S
  //
  int n_nodsSub = 0;
//  double rho = K.val[0];
  if (fixing_nodes_or_dof>0){
    sc_size = fixing_nodes_or_dof*dofPerNode;
    n_nodsSub = round(K.n_row_cmprs/dofPerNode);
  }
  //
  //##########################################################################################
  //
  Matrix S1;
  Matrix K_rr;
  Matrix K_rs;
  eslocal i_start = 0;

  if (K.n_row_cmprs < sc_size){
      sc_size = K.n_row;
    fixing_nodes_or_dof=0;
  }




  eslocal nonsing_size = K.n_row_cmprs- sc_size - i_start;
  eslocal j_start = nonsing_size;
  SEQ_VECTOR <eslocal > permVec;

  permVec.resize(K.n_row_cmprs);
  SEQ_VECTOR < SEQ_VECTOR < eslocal > > vec_I1_i2(K.n_row_cmprs, SEQ_VECTOR<eslocal >(2, 1));
  eslocal offset = 0;
  //

  eslocal *I_row_indices_p = new eslocal [K.nnz] ;
  eslocal *J_col_indices_p = new eslocal [K.nnz] ;
  SEQ_VECTOR <eslocal > tmp_vec_s;
  tmp_vec_s.resize(sc_size);
  eslocal v1, n_mv, cnt_permut_vec;
  SEQ_VECTOR <eslocal >::iterator it;
  SEQ_VECTOR <eslocal > fix_dofs;
  fix_dofs.resize(sc_size);
  Matrix K_modif;

  double di=1,dj=1;
  eslocal cnt_iter_check_nonsing=0;

  K_modif = K; // TODO not necessary to do



  SEQ_VECTOR <double> tmp_approx_max_eig;
  tmp_approx_max_eig.resize(K.n_row_cmprs);

  // diagonal scaling of K_modif:
  // K_modif[i,j] = K_modif[i,j]/sqrt(K_modif[i,i]*K_modif[j,j]);
  for (eslocal i = 0;i<K_modif.n_row_cmprs;i++){
    if (diagonalScaling) {
      di=K_modif.val[K_modif.i_ptr[i]-offset];
    }
    for (eslocal j = K_modif.i_ptr[i];j<K_modif.i_ptr[i+1];j++){
      if (diagonalScaling) {
        dj=K_modif.val[
          K_modif.i_ptr[K_modif.j_col[j-offset]-offset]-offset];
      }
      K_modif.val[j-offset]/=sqrt(di*dj);
      tmp_approx_max_eig[i]+=fabs(K.val[j-offset]);
    }
  }


  if (permutVectorActive<3){
    // set row permVec = {0,1,2,3,4,...,K.rows};
    if (fixing_nodes_or_dof==0 || (permutVectorActive==0)){
      for (eslocal i=0; i<K.n_row_cmprs; ++i) { permVec[i]=i;} // 0 1 2 K.rows-1
    }
    else
    {
      for (eslocal i=0; i<n_nodsSub; ++i) { permVec[i]=i;} // 0 1 2 n_nodsSub-1
    }
  }
//
  if (permutVectorActive==1){
    srand(time(NULL));
//    srand(0); // random will be constant until next compiling

    if (fixing_nodes_or_dof==0){
      random_shuffle ( permVec.begin(), permVec.end() );
    }
    else
    {
      std::srand(std::time(0));
      std::random_shuffle ( permVec.begin(), permVec.begin()+n_nodsSub);
      for (eslocal i=n_nodsSub;i>0;i--){
        for (eslocal j=0;j<dofPerNode;j++){
          permVec[dofPerNode*i-1-j] = dofPerNode*permVec[i-1]+j;
        }
      }
    }

    sort(permVec.begin(),permVec.begin()+nonsing_size);
    sort(permVec.begin()+nonsing_size,permVec.end());
  }
  else if (permutVectorActive==2){
    // random permutation
    n_mv = 0;                     // n_mv = size(unique(tmp_vec_s)) has to be equal to sc_size
    cnt_permut_vec=0;
    srand(time(NULL));
    // loop controls, if series 'tmp_vec_s' with unique integers has suffisciant dimension.
    // If not, missing numbers are added and checked again.
    do {
      for (eslocal i = 0;i<(sc_size-n_mv);i++){
        v1 = rand() % K_modif.n_row_cmprs;
        tmp_vec_s[n_mv+i]=v1;
      }
      it=tmp_vec_s.begin();
      std::sort(tmp_vec_s.begin(), tmp_vec_s.end());
      it = std::unique(tmp_vec_s.begin(), tmp_vec_s.end());
      n_mv = distance(tmp_vec_s.begin(),it);
      cnt_permut_vec++;
   } while (n_mv != sc_size && cnt_permut_vec < 100);
    //
    eslocal ik=0,cnt_i=0;
    for (eslocal i = 0;i<permVec.size();i++){
      if (i==tmp_vec_s[ik]){
        permVec[ik+nonsing_size]=tmp_vec_s[ik];
        ik++;
      }
      else{
        permVec[cnt_i]=i;
        cnt_i++;
      }
    }
  }


  //      r = permVec[0:nonsing_size-1]     (singular DOFs)
  //      s = permVec[nonsing_size:end-1]   (non-singular DOFs)
  if (permutVectorActive>0){
//
    for (eslocal i = 0; i < K.n_row_cmprs;i++){
      vec_I1_i2[i][0]=permVec[i];
      vec_I1_i2[i][1]=i; // position to create revers permutation
    }
//






    sort(vec_I1_i2.begin(), vec_I1_i2.end(), compareVecVec);

    // permutations made on matrix in COO format
    K_modif.CSR2COO();


    vector < TRIPLET > triplet(K.nnz);

    int I_index,J_index;
    bool unsymmetric = false;
    for (int i = 0;i<K_modif.nnz;i++){
       I_index = vec_I1_i2[K_modif.i_coo_cmpr[i]][1];
       J_index = vec_I1_i2[K_modif.j_col[i]][1];
       if (unsymmetric || I_index<=J_index){
         I_row_indices_p[i]=I_index;
         J_col_indices_p[i]=J_index;
         triplet[i].I =  I_index;
         triplet[i].J =  J_index;
       }
       else{
         I_row_indices_p[i]=J_index;
         J_col_indices_p[i]=I_index;
         triplet[i].I =  J_index;
         triplet[i].J =  I_index;
       }
       triplet[i].V =  K_modif.val[i];
     }


    for (int i = 0; i < K_modif.nnz; i++){
       K_modif.i_coo_cmpr[i] = triplet[i].I;
       K_modif.j_col[i]      = triplet[i].J;
    }

    K_modif.sortAndUniqueCOO(triplet);
    K_modif.COO2CSR();

  }

//
    for (eslocal i = 0;i<sc_size;i++) fix_dofs[i]=permVec[nonsing_size + i] + offset;
    K_rr.getSubDiagBlockmatrix(K_modif,K_rr,i_start, nonsing_size);
    cnt_iter_check_nonsing++;

  delete [] I_row_indices_p;
  delete [] J_col_indices_p;
//
  K_rs.getSubBlockmatrix_rs(K_modif,K_rs,i_start, nonsing_size,j_start,sc_size);
  K_rr.options2 = c_options2;
  Matrix K_rs_copy = K_rs;


  if (dbg > 0){
    cout << "\t\t\toptions2.size() "  << c_options2.size() << endl;
    for (map<string,string>::const_iterator it=c_options2.begin(); it!=c_options2.end(); ++it)
        cout << "\t\t\t" << it->first << " => " << it->second << '\n';

  }



    S1.getSubDiagBlockmatrix(K_modif,S1,nonsing_size,sc_size);
    Matrix KsrInvKrrKrs;

    if (dbg>0){
        fprintf(stderr, "\n\n\n%s %d : \n\n\n\n", __FILE__, __LINE__);
        cout << K_rr.n_col << endl;
    }


  if (K_rr.n_col!=0)
  {
      cout << "   =====  " << c_options2["linear_solver"] << endl;
      K_rr.sym_factor(c_options2["linear_solver"]);
      Matrix _R;
      bool checkOrthogonality_ = true;
      K_rr.num_factor(_R,checkOrthogonality_);



      Matrix InvKrrKrs;
      K_rs_copy.label = "K_rs_copy";
      K_rs.label = "o";

      K_rs.CSRorCOO2DNS(false,false);
      K_rr.solve_system(K_rs,InvKrrKrs);
      K_rs_copy.mult(InvKrrKrs,KsrInvKrrKrs,false);
    }
  S1.CSRorCOO2DNS(false,false);

    if (K_rr.n_col!=0){
        for (int i = 0; i < S1.n_row_cmprs; i++) {
            for (int j = 0; j < S1.n_row_cmprs; j++) {
                S1.dense[i + j * S1.n_row_cmprs ] -=
                    KsrInvKrrKrs.dense[i + j * S1.n_row_cmprs ];
            }
        }
    }

  int n_S2 = (S1.n_row_cmprs + 1) * S1.n_row_cmprs * 0.5;
  double *S2 = new double[n_S2];

  int cnt = 0;
  for (int j = 0; j<S1.n_col; j++){
    for (int i = 0; i<=j; i++){
        S2[cnt] = S1.dense[i + j * S1.n_row_cmprs];
//        printf("%d %d %e \n", i, j,S1.dense[i + j * S1.n_row_cmprs]);
        cnt++;
    }
  }


// EIGENVALUES AND EIGENVECTORS OF SCHUR COMPLEMENT

  char JOBZ = 'V';
  char UPLO = 'U';
  double *W = new double[S1.n_col];
  double *Z = new double[S1.n_col*S1.n_col];
  MKL_INT info;
  MKL_INT ldz = S1.n_col;
  info = LAPACKE_dspev (LAPACK_COL_MAJOR, JOBZ, UPLO, S1.n_col, &(S2[0]), W, Z, ldz);
  if (info){
    cout <<"info = " << info << " something wrong with Schur complement in SparseSolverCPU::generalIinverse" << endl;
  }
  delete [] S2;

// IDENTIFICATIONS OF ZERO EIGENVALUES
  eslocal defect_K_in;// R_s_cols;
  double ratio;
  eslocal itMax = twenty < S1.n_row_cmprs ? twenty : S1.n_row_cmprs;
  for (eslocal i = itMax-1; i > 0;i--){
    ratio = fabs(W[i-1]/W[i]);
    if (ratio < jump_in_eigenvalues_alerting_singularity){
      defect_K_in=i;
      break;
    }
  }
  if (get_n_first_and_n_last_eigenvals_from_dense_S!=0){
    int i1i = get_n_first_and_n_last_eigenvals_from_dense_S;
    if (i1i>S1.n_row_cmprs){i1i=S1.n_row_cmprs;}
    cout<<"eigenvals of S d{1:" << i1i << "} and d{" <<
         S1.n_row_cmprs-get_n_first_and_n_last_eigenvals_from_dense_S+2 << ":"<< S1.n_row_cmprs<< "}\n";

    for (eslocal i = 0 ; i < S1.n_row_cmprs; i++){
      if (i < get_n_first_and_n_last_eigenvals_from_dense_S ||
            i > S1.n_row_cmprs-get_n_first_and_n_last_eigenvals_from_dense_S){
        cout<< i+1 <<":\t"<< W[i] << "\n";
      }
    }
  }
// --------------- CREATING KERNEL R_s FOR SINGULAR PART (SCHUR COMPLEMENT)
  Matrix R_s;
  R_s.zero_dense(S1.n_row_cmprs, defect_K_in);
  eslocal cntR=0;
  for (eslocal j = 0; j < defect_K_in; j++){
    for (eslocal i = 0; i < R_s.n_row_cmprs; i++){
        R_s.dense[j*R_s.n_row_cmprs+ i] = Z[j*R_s.n_row_cmprs+ i];
        cntR++;
    }
  }
// --------------- CREATING KERNEL R_r FOR NON-SINGULAR PART
  Matrix R_r, KrsRs;
  R_r.zero_dense(K_rr.n_row_cmprs,defect_K_in);
  KrsRs.zero_dense(K_rr.n_row_cmprs,defect_K_in);
  if (K_rr.n_col!=0){
    K_rs_copy.mult(R_s,KrsRs,true);
    K_rr.solve_system(KrsRs,R_r);
  }
               //                                               |
//// --------------- CREATING WHOLE KERNEL Kplus_R = [ (R_r)^T (R_s)^T ]^T

  Kplus_R.zero_dense(K.n_row_cmprs,defect_K_in);
  cntR=0;
  int R_r_rows = R_r.n_row_cmprs;
  for (eslocal j = 0; j < Kplus_R.n_col; j++){
    for (eslocal i = 0; i < R_r_rows; i++){
      if (diagonalScaling){
        di=K.val[K.i_ptr[permVec[i]]-offset];
      }
      Kplus_R.dense[j*Kplus_R.n_row_cmprs + permVec[i]] = R_r.dense[j*R_r_rows + i]/sqrt(di);
      cntR++;
    }
    for (eslocal i = 0; i < R_s.n_row_cmprs; i++){
      if (diagonalScaling){
        di=K.val[K.i_ptr[permVec[i+R_r_rows]]-offset];
      }
      Kplus_R.dense[j*Kplus_R.n_row_cmprs + permVec[i+R_r_rows]] =-R_s.dense[j*R_s.n_row_cmprs + i]/sqrt(di);
      cntR++;
    }
  }


  Matrix::GramSchmidtOrtho(Kplus_R);
  SEQ_VECTOR <eslocal > null_pivots1;
  Kplus_R.getNullPivots(null_pivots1);

// norm of product K*R: second matrix has to be in dense format!!!
//11 - max(eig(K))
//
//  std::vector <double>::iterator  it2;
//  it2 = std::max_element(tmp_approx_max_eig.begin(),tmp_approx_max_eig.end(),compareDouble);
//  double lmx_K_approx       = *it2;
//  double tmp_Norm_K_R       = Matrix::getNorm_K_R(K,Kplus_R);
//  double norm_KR_d_pow_2_approx   = (tmp_Norm_K_R*tmp_Norm_K_R)/(lmx_K_approx*lmx_K_approx);
//  int defect_d                 = Kplus_R.n_col;


//  double rho = lmx_K_approx;
//  cout << "max sigma approx:        " << rho << "\n";
//  cout << "defect(K):               " << defect_K_in <<"\n";
//  cout << "norm_KR_approx:          " << sqrt(norm_KR_d_pow_2_approx) <<"\n";


//
  delete [] W;
  delete [] Z;
}

void Matrix::GramSchmidtOrtho(Matrix &R_in){
    int rows = R_in.n_row_cmprs;
    int cols = R_in.n_col;
  double *w = new double [rows];
  double *R = new double [cols*cols];
  memset(R,0,(cols*cols) * sizeof(double));

  for (int j = 0;j<cols;j++){
    memcpy( w, &(R_in.dense[j*rows]) , sizeof( double ) * rows);
    for (int i = 0;i<j;i++){
      R[j*cols+i] = Matrix::dot(w, &(R_in.dense[i*rows]),rows);
      for (int k=0;k<rows;k++){
        w[k]-=R_in.dense[i*rows+k]*R[j*cols+i];
      }
    }
    R[j*cols+j] = sqrt(Matrix::dot(w,w,rows));
    for (int k=0;k<rows;k++){
      R_in.dense[j*rows+k] = w[k]/R[j*cols+j];
    }
  }
  delete [] w;
  delete [] R;
}




double Matrix::getNorm_K_R(Matrix & K, Matrix &R_in_dense_format){

      double norm_AR_row=0,norm_AR = 0;
      double * AR =  new double [K.n_row_cmprs];


      for (int i = 0;i<R_in_dense_format.n_col;i++){
        memset(AR,0,R_in_dense_format.n_row_cmprs * sizeof(double));
        K.mult(R_in_dense_format.dense.data(), AR, true , 1 , -1);
        norm_AR_row=0.0;
        for (int j = 0; j < R_in_dense_format.n_row_cmprs;j++){
          norm_AR_row+=AR[j]*AR[j];
        }
        norm_AR+=norm_AR_row;
      }
      norm_AR=sqrt(norm_AR);
      delete [] AR;
      return norm_AR;
}
