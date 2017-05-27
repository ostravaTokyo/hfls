#include "Matrix.hpp"
#include <math.h>


using namespace std;


Matrix::Matrix(){
    nnz = 0;
    numel = 0;
    n_row = 0;
    n_col = 0;
    n_row_cmprs= 0;
    DNS_reducedZeroRows = false;
    DNS_transposed = false;
    diss_scaling = -1;
    msglvl = -1;
}

Matrix::~Matrix()
{

}


bool Matrix::cmp_int_int_I(int_int_dbl a,int_int_dbl b)
{ return (a.I < b.I); }
bool Matrix::cmp_int_int_J(int_int_dbl a,int_int_dbl b)
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
            i_coo_cmpr.shrink_to_fit();

            j_col.resize(nnz);
            j_col.shrink_to_fit();

            val.resize(nnz);
            val.shrink_to_fit();
        }
/* Sorting [I, J, V]  --------------------------------------------------------*/
//        if (format_ != 2){
        vector < int_int_dbl > tmpVec;
        tmpVec.resize(i_coo_cmpr.size());
        for (int i = 0; i < i_coo_cmpr.size();i++){
            tmpVec[i].I = i_coo_cmpr[i];
            tmpVec[i].J = j_col[i];
            tmpVec[i].V = val[i];
        }
/*AAA*/
        sortAndUniqueCOO(tmpVec);
/*BBB*/
        tmpVec.clear();
        tmpVec.shrink_to_fit();

        //}
        if (format_ == 1){
            COO2CSR();
        }
    }
    numel = n_row_cmprs * n_col;
}


void Matrix::sortAndUniqueCOO(vector <int_int_dbl> &tmpVec){
/* sort according to index I -------------------------------------------------*/
    sort(tmpVec.begin(),tmpVec.end(),Matrix::cmp_int_int_I);
/* partial sorting according to index J (per sets with same I)----------------*/



//    if (label == "Fc_clust"){
//        cout << "A ===========================================\n";
//        for (int i = 0; i < tmpVec.size();i++){
//           cout << i << ":::  "
//                << tmpVec[i].I << " "
//                <<  tmpVec[i].J  << " "
//                <<  tmpVec[i].V  << "\n";
//        }
//    }


    nnz = tmpVec.size();

    int startInd = 0, endInd = 1;
    int tmpVecIprev = tmpVec[0].I;
    for (int i = 1 ; i < nnz; i ++){
        if (tmpVec[i].I == tmpVecIprev){
            endInd++;
        }
        if (tmpVec[i].I != tmpVecIprev || (tmpVec[i].I == tmpVecIprev &&
                 i == nnz - 1))
        {
            sort(tmpVec.begin() + startInd,
                              tmpVec.begin() + (endInd  ),Matrix::cmp_int_int_J);
            startInd = i;
            endInd = i + 1;
        }
        tmpVecIprev = tmpVec[i].I;
    }



//    if (label == "Fc_clust"){
//        cout << "B ===========================================\n";
//        for (int i = 0; i < tmpVec.size();i++){
//           cout << i << ":::  "
//                <<  tmpVec[i].I << " "
//                <<  tmpVec[i].J  << " "
//                <<  tmpVec[i].V  << "\n";
//        }
//    }


/* Cumulating duplicated A[I,J] elements--------------------------------------*/
    int prevInd_I = -1;
    int prevInd_J = -1;
    int counter = 0;
    int cnt_j = -1;

    if (l2g_i_coo.size() != tmpVec.size()){
        l2g_i_coo.resize(tmpVec.size());
    }
    if (i_coo_cmpr.size() != tmpVec.size()){
        i_coo_cmpr.resize(tmpVec.size());
    }
    if (j_col.size() != tmpVec.size()){
        j_col.resize(tmpVec.size());
    }
    if (val.size() != tmpVec.size()){
        val.resize(tmpVec.size());
    }


    int init_nnz = nnz;
    for (int i = 0 ; i < init_nnz; i ++){
        if (prevInd_I != tmpVec[i].I){
            l2g_i_coo[counter] = tmpVec[i].I;
            counter++;
        }
        if (prevInd_I == tmpVec[i].I && prevInd_J == tmpVec[i].J){
            val[cnt_j] += tmpVec[i].V;
            nnz--;
        }
        else {
            cnt_j++;
            i_coo_cmpr[cnt_j] = tmpVec[i].I ;
            j_col[cnt_j] = tmpVec[i].J;
            val[cnt_j] = tmpVec[i].V;
        }
        prevInd_I = tmpVec[i].I;
        prevInd_J = tmpVec[i].J;
    }


//    if (label == "Fc_clust"){
//        cout << "C ===========================================\n";
//        for (int i = 0; i < init_nnz;i++){
//           cout << i << ":::  "
//                <<  i_coo_cmpr[i]  << " "
//                <<  j_col[i]  << " "
//                <<  val[i]  << "\n";
//        }
//    }







    l2g_i_coo.resize(counter );
    l2g_i_coo.shrink_to_fit();

    i_coo_cmpr.resize(nnz);
    i_coo_cmpr.shrink_to_fit();

    j_col.resize(nnz);
    j_col.shrink_to_fit();

    val.resize(nnz);
    val.shrink_to_fit();

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
    i_coo_cmpr.shrink_to_fit();
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
    i_ptr.shrink_to_fit();
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
/*                                                                            */
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
        i_coo_cmpr.shrink_to_fit();
    }
    else if (format_ == 1){
        i_ptr.clear();
        i_ptr.shrink_to_fit();
    }

    j_col.clear();
    j_col.shrink_to_fit();

    val.clear();
    val.shrink_to_fit();
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
       if (dense[i] != 0){
           nnz++;
       }
    }
    i_coo_cmpr.resize(nnz);
    j_col.resize(nnz);
    val.resize(nnz);
    double _val;
    int cnt = 0;
//    cout << "_n_row " << _n_row << "\n";
    for (int i = 0; i < _n_row; i++) {
        for (int j = 0; j < n_col; j++) {
            if (DNS_transposed){
                _val = dense[i * n_col + j];
            }
            else{
                _val = dense[i + j * _n_row];
            }
            if (_val != 0){
                i_coo_cmpr[cnt] = i;
                j_col[cnt] = j;
                val[cnt] = _val;
                cnt++;
            }
        }
    }
    dense.clear();
    dense.shrink_to_fit();
    format = 0;
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
    string path2matrix = folder+"/dump_" + nameOfMat + "_" +
                                                     to_string(indOfMat)+".txt";
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
//        if (format < 2) {
//            CSRorCOO2DNS(true,true);
//        }
//        cout << "=====format======" << format << "\n";
//
//        int _n_row, _ni, _nj, _nm;
//        if (DNS_reducedZeroRows) {
//            _n_row = n_row_cmprs;
//        }
//        else{
//            _n_row = n_row;
//        }
//        if (DNS_transposed){
//            _ni = n_col;
//            _nj = _n_row;
//            _nm = n_col;
//        }
//        else{
//            _ni = _n_row;
//            _nj = n_col;
//            _nm = _n_row;
//        }
//
//        for (int i = 0; i < _ni; i++) {
//            for (int j = 0; j < _nj; j++){
//                V = dense[i + j * _nm];
//                fprintf(fp, "%.16f\t", V);
//            }
//            fprintf(fp, "\n");
//        }
//
//        if (format_ != 2){
//            if (format_ == 0) {
//                DNS2COO();
//            }
//            else if (format_ == 1){
//                DNS2CSR();
//            }
//        }
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


void Matrix::mat_mult_dense(Matrix&A, string A_NorT, Matrix&B, string B_NorT){

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


    mult(&(X_in.dense[0]),&(X_out.dense[0]), NorT, _n_col_rhs);
}


void Matrix::mult(const double x_in[], double  x_out[], bool NorT, int n_rhs){


//    if (label == "K_")
//        cout << "KpKpKpKKpKpKpKpKpKpKpKpKpKKKKKKKKKKKKpK" << endl;

    for (int i = 0; i < n_row_cmprs; i++) {
        for (int j = i_ptr[i]; j < i_ptr[i + 1]; j++) {
            if (symmetric > 0 ){
                for (int k = 0; k < n_rhs; k++){
                    x_out[i + k * n_row_cmprs] +=
                            val[j] * x_in[j_col[j] + k * n_col];
//                    if (label == "K_")
//                        cout << val[j] * x_in[j_col[j] + k * n_col];
                    if (j_col[j] != i){
                        x_out[j_col[j] + k * n_row_cmprs] +=
                                      val[j] * x_in[i + k * n_col];
//                    if (label == "K_")
//                        cout << val[j] * x_in[i + k * n_col];
                    }
                }
            }
            else {
                for (int k = 0; k < n_rhs; k++){
                    if (NorT){
                        x_out[i + k * n_row_cmprs] +=
                                val[j] * x_in[j_col[j] + k * n_col];
                    }
                    else {
                        x_out[j_col[j] + k * n_row_cmprs] +=
                                      val[j] * x_in[i + k * n_col];
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


void Matrix::getBasicMatrixInfo(){
    printf("< %s >              \n", label.c_str());
    printf("symmetric           %d (0: unsym, 1: sym. lower tr., 2: sym. upper tr.)\n", symmetric);
    printf("format              %d (0: coo, 1: csr, 2: dense)\n", format);
    printf("nnz:                %d \n", nnz);
    printf("n_row_cmprs:        %d \n", n_row_cmprs);
    printf("n_row:              %d \n", n_row);
    printf("n_col:              %d \n", n_col);
    printf("i_ptr.size():       %lu \n", i_ptr.size());
    printf("i_coo_cmpr.size():  %lu \n", i_coo_cmpr.size());
    printf("j_col.size():       %lu \n", j_col.size());
    printf("val.size():         %lu \n", val.size());

    double density = 100 * nnz / (double)( n_row_cmprs * n_col );

    printf("sparsity:           %1.2f %% \n", density);
}

int Matrix::ij(int i_, int j_){
   return i_ + j_ * n_row_cmprs;
}


void Matrix::getNullPivots(vector < int > & null_pivots){


  int cols = n_col;
  int rows = n_row_cmprs;

  vector <double> N(dense);
  vector <double>::iterator  it;
  int I,colInd,rowInd;
  double *tmpV = new double[rows];
  int *_nul_piv = new int [rows];
  double pivot;
  int tmp_int;
  for (int i = 0;i<rows;i++){
      _nul_piv[i]=i;
  }

//  auto ij = [&]( int ii, int jj ) -> int {
//      return ii + rows * jj;
//  };

  for (int j=0;j<cols;j++){
    it = max_element(N.begin(),N.end()-j*rows, compareDouble);
    I = it - N.begin();
    colInd = I/rows;
    rowInd = I-colInd*rows;
    for (int k=0;k<cols-j;k++){
      tmpV[k] = N[ij(rows-1-j,k)];
      N[ij(rows-1-j,k)] = N[ij(rowInd,k)];
      N[ij(rowInd,k)]= tmpV[k];
    }
    tmp_int = _nul_piv[rowInd];
    _nul_piv[rowInd] = _nul_piv[rows-1-j];
    _nul_piv[rows-1-j] = tmp_int;
    if (cols - 1 - j != colInd) {
        memcpy( tmpV, &(N[ij(0,cols-1-j)]) , sizeof( double ) * rows);
        memcpy( &(N[ij(0, cols - 1 - j)]), &(N[ij(0, colInd)]), sizeof( double ) * rows);
        memcpy( &(N[ij(0,colInd)]),tmpV , sizeof( double ) * rows);
    }
    pivot = N[ij(rows-1-j,cols-1-j)];
    for (int J=0;J<cols-j-1;J++){
      for (int I=0;I<rows-j;I++){
        N[ij(I,J)] -= N[ij(I,cols-1-j)]*N[ij(rows-1-j,J)]/pivot;
      }
    }
  }
  null_pivots.resize(cols);
  for (int i = 0;i<cols;i++){
    null_pivots[i] = _nul_piv[rows-1-i];
  }
  sort(null_pivots.begin(),null_pivots.end());
//
  delete [] _nul_piv;
  delete [] tmpV;
//
}

void Matrix::symbolic_factorization(){

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

void Matrix::numeric_factorization(){
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

void Matrix::numeric_factorization(Matrix &R,bool checkOrthogonality){
//
#ifdef DISSECTION
    numeric_factorization();
    int n0;
    diss_get_kern_dim(*diss_dslv, &n0);
    R.zero_dense(n_row,n0);
    diss_get_kern_vecs(*diss_dslv, &R.dense[0]);
    if (checkOrthogonality){
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




#endif

}

void Matrix::diss_solve(Matrix &B, Matrix &X){

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
    diss_solve_n(*diss_dslv,&X.dense[0],
                 _n_col, projection, trans);
#endif
}

void Matrix::factorization(){
#ifdef USE_PARDISO
    InitializeSolve();
#endif
}

void Matrix::factorization(vector <int> & _nullPivots){
#ifdef USE_PARDISO
    // regularization:
    int j, cnt = 0;
    for (int i = 0; i < n_row_cmprs; i++) {
        j = i_ptr[i];
        if (i == _nullPivots[cnt]){
            val[j] *= 100;
            cnt++;
        }
    }
    factorization();
#endif
}



void Matrix::CsrElementByElement(){
#ifdef buildCSR_elemByElem

    bool symMatrix=true;
    int *indK;
    vector < vector < int > > forCSRformat;
    forCSRformat.resize(fem->domain->neqSub,vector<int>(0));
    int elemDOFs = 24;
    for (int i=0;i<fem->domain->n_elementsSub;i++){
      indK = stif_glob_number[i].ieq;
      for (int j=0;j<elemDOFs;j++){
        for (int k=0;k<elemDOFs;k++){
          if (symMatrix && (indK[k]>=indK[j]) || !symMatrix){
            forCSRformat[indK[j]].push_back(indK[k]);
          }
        }
      }
    }
    nnz_K = 0;
    for (int i = 0;i<forCSRformat.size();i++){
      sort(forCSRformat[i].begin(),forCSRformat[i].end());
      itv = unique (forCSRformat[i].begin(), forCSRformat[i].end());
      forCSRformat[i].resize( distance(forCSRformat[i].begin(),itv) );
      nnz_K+=forCSRformat[i].size();
    }

#endif
}



void Matrix::InitializeSolve()
{

#ifdef USE_PARDISO
    if (msglvl<0)
        msglvl = 0;
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

//    iparm[0] = 1;         /* No solver default */
//    iparm[1] = 2;         /* Fill-in reordering from METIS */
//    iparm[3] = 0;         /* No iterative-direct algorithm */
//    iparm[4] = 0;         /* No user fill-in reducing permutation */
//    iparm[5] = 0;         /* Write solution into x */
//    iparm[7] = 2;         /* Max numbers of iterative refinement steps */
//    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
//    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
//    iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
//    iparm[13] = 0;        /* Output: Number of perturbed pivots */
//    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
//    iparm[18] = -1;       /* Output: Mflops for LU factorization */
//    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
//    iparm[34] = 1;        /* PARDISO use C-style indexing for ia and ja arrays */
//    maxfct = 1;           /* Maximum number of numerical factorizations. */
//    mnum = 1;         /* Which factorization to use. */
//    msglvl = 1;           /* Print statistical information in file */
//    error = 0;            /* Initialize error flag */



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
    n = n_row_cmprs;
    MKL_INT n1 = n_row_cmprs;
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n1, &val[0], &i_ptr[0], &j_col[0],
            &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during symbolic factorization: %d", error);
        exit (1);
    }

    phase = 22;

    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n1, &val[0], &i_ptr[0], &j_col[0],
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
    int _n_row, _n_col;
    if (B.DNS_transposed){
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

    if (X.nnz == 0 ){
        X.zero_dense(_n_row,_n_col);
    }


    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &val[0], &i_ptr[0], &j_col[0], &idum, &nrhs, iparm, &msglvl,
            &B.dense[0], &X.dense[0], &error);

    if ( error != 0 )
    {
        printf ("\nERROR during solution: %d", error);
        exit (3);
    }
//    printf ("\nSolve completed ... ");
#endif
}

void Matrix::FinalizeSolve(int _i)
{
    if (solver == 0){
#ifdef USE_PARDISO
        double ddum;          /* Double dummy */
        MKL_INT idum;         /* Integer dummy. */
        phase = -1;           /* Release internal memory. */
        PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
                 &n, &ddum, &i_ptr[0], &j_col[0], &idum, &nrhs,
                iparm, &msglvl, &ddum, &ddum, &error);

    //    if (_i == 0 ) {
    //        for (int i = 0 ; i < 64 ; i++ ){
    //            if (iparm[i] < 0 ){
    //                cout <<"iparm[" << i <<"] = " << iparm[i] << "\n";
    //            }
    //        }
    //    }
#endif
    }
    else if (solver == 1)
    {
#ifdef DISSECTION
        diss_free(*diss_dslv);
#endif
    }
}

void Matrix::testPardiso(string folder){
#ifdef USE_PARDISO
/*------------------------------ TEST ------------------------------*/
    Matrix A;
    int n = 11;
    int cnt = 0;
    int nrhs = 5;
    A.nnz = 2 * n - 1;
    A.n_row = n;
    A.n_row_cmprs = n;
    A.n_col = n;

    Matrix B;
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


    A.InitializeSolve();
    Matrix U;
    A.solve(B, U);
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

void Matrix::updateCOOstructure(vector <int_int_dbl >& vec_, Matrix &denseMat,int i_shift, int j_shift){


    int_int_dbl i_tmpvec;
    for (int j = 0; j < denseMat.n_col;j++){
        for (int i = 0; i < denseMat.n_row_cmprs;i++){
            i_tmpvec.I = i_shift + i;
            i_tmpvec.J = j_shift + j;
            i_tmpvec.V = denseMat.dense[i + j * denseMat.n_row_cmprs ];
            if (i_tmpvec.J >= i_tmpvec.I && i_tmpvec.V != 0){
                vec_.push_back(i_tmpvec);
            }
        }
    }

}
