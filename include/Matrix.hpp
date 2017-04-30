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
//    bool printCooOrDense;
    MKL_INT m;
    //coo row compressed

    //csr
    //dns


    void mult(const Matrix &, Matrix &, bool);
    void mult(const double[], double [] , bool, int);
    void CsrElementByElement();
    void COO2CSR();
    void CSR2COO();
    void CSRorCOO2DNS(bool, bool);
    void dummyFunction(int,int, int, double);
    void DNS2CSR();
    void DNS2COO();
    void compresRows();
    static Matrix CreateCopyFrom(const Matrix&);
    void readCooFromFile(string, int,int,int);
    void printToFile(string,int,bool);
    void InitializeSolve();
    void solve(Matrix&, Matrix&);
    void FinalizeSolve(int);
    void setZero();
    MKL_INT iparm[64];
    MKL_INT mtype;
    MKL_INT nrhs;
    MKL_INT n;
    void *pt[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    static void testPardiso();
    void factorization();


    double norm2();

};
#endif





