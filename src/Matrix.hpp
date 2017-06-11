#ifndef __SPARSEMATRIX_HPP_INCLUDED
#define __SPARSEMATRIX_HPP_INCLUDED
#include <string>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <fstream> 
#include "mkl.h"
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include <iomanip>
#include "math.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "map"

#include <unistd.h>
#ifdef DISSECTION
#include "Dissection.hpp"
#endif




//    symmetric = 0;  // 0-not, 1-lower tr., 2-upper tr.
//    format = 1;     // 0-coo, 1-csr, 2-dense
//    offset = 0;


using namespace std;

struct int_int_dbl
{
    int I;
    int J;
    double V;
};

bool myfunction (int i,int j);

class Matrix{
public:
    Matrix();
    Matrix(string );
    ~Matrix();
    void zero_dense(int,int);
    void zero_dense(int,int,bool);
    vector < double >  val;
    vector < int >     i_ptr; // new MKL_INT[crs->m+1];
    vector < int >     j_col; // new MKL_INT[crs->nnz];
    vector < int >     i_coo_cmpr;
    vector < int >     l2g_i_coo;
    map <int,int>      g2l_i_coo;
    vector < double >  dense;
    int format;               /* COO, CSR, DNS */
    int symmetric;            /* 0-no, 1-lower triang., 2-upper triang. */
    int nnz;
    int n_row;
    int n_row_cmprs;
    int n_col;
    int numel;
    bool DNS_reducedZeroRows;
    bool DNS_transposed;
    string label;
    int order_number;

    int ij(int, int);

//    bool printCooOrDense;
    MKL_INT m;
    //coo row compressed

    //csr
    //dns


    void mult(const Matrix &, Matrix &, bool);
    void mult(const double[], double [] , bool, int,int);
    void CsrElementByElement();
    void COO2CSR();
    void CSR2COO();
    void CSRorCOO2DNS(bool, bool);
    void dummyFunction(int,int, int, double);
    void DNS2CSR();
    void DNS2COO();
    void compresRows();
    void readCooFromFile(string, int,int,int);
    void printToFile(string,string,int,bool);
    void InitializeSolve();
    void solve(Matrix&, Matrix&);
    void FinalizeSolve(int);
    void setZero();
    void getNullPivots(vector < int > & );
    void sortAndUniqueCOO(vector < int_int_dbl > &);
    void submatrix_row_selector(Matrix &, vector <int> &);

    void getBasicMatrixInfo();

    MKL_INT iparm[64];
    MKL_INT mtype;
    MKL_INT nrhs;
    MKL_INT n;
    void *pt[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    static void testPardiso(string);
    void factorization(vector < int > &);
    void factorization();

    uint64_t *diss_dslv;
    int diss_num_threads;
    int diss_sym;
    int diss_decomposer;
    int diss_indefinite_flag;
    int diss_scaling;
    int diss_int0;
    int solver;
    double diss_eps_pivot;
    void symbolic_factorization();
    void numeric_factorization();
    void numeric_factorization(Matrix &,bool);
    void diss_solve(Matrix &, Matrix &);

    void mat_mult_dense(Matrix&, string, Matrix&, string);


    double norm2();

    static bool cmp_int_int_I(int_int_dbl ,int_int_dbl );
    static bool cmp_int_int_J(int_int_dbl ,int_int_dbl );
    static bool compareDouble(double, double);
    static void updateCOOstructure(vector <int_int_dbl > &, Matrix &, int, int);

    static bool test_of_Bc_constraints(Matrix &);

    static void print_int_vector(vector <double> &, int, int);

    static void getSingularVal_DNS(Matrix &, Matrix &, int);
    static void getSingularVal_DNS(Matrix , Matrix &, int, int);

    static void getEigVal_DNS(Matrix &, Matrix &, int );
    static void getEigVal_DNS(Matrix  , Matrix &, int, int);




};
#endif
