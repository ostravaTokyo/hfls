#ifndef __SPARSEMATRIX_H_INCLUDED
#define __SPARSEMATRIX_H_INCLUDED
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



struct int_int_dbl
{
    int I;
    int J;
    int V;
};

bool myfunction (int i,int j);



class SparseMatrix{

public:
    SparseMatrix();
    SparseMatrix(std::string path2matrix);
    //coo
    std::vector < int > i_coo;
    std::vector < int > j_coo;
    std::vector < double > v_coo;

    int nnz;
    int n_row;
    int n_col;
    MKL_INT m;
    //coo row compressed
    std::vector < int >   i_coo_orig;
    std::vector < int >   l2g_i_coo;



    //csr
    std::vector < double > a;   // = new double[crs->nnz];
    std::vector < MKL_INT > ia; // new MKL_INT[crs->m+1];
    std::vector < MKL_INT > ja; // new MKL_INT[crs->nnz];
    //dns
    std::vector < double > adns; // ()

    bool isCOO;
    bool isCSR;
    bool isDNS;


    void COO2CSR();
    void CSR2COO();
    void CSR2DNS(int);
    void DNS2CSR(int);
    void RemoveLower();
    void compresRows();
    void CreateCopyFrom(const SparseMatrix&, int);
    static void dcsradd(SparseMatrix&, SparseMatrix&, SparseMatrix&);
    static void Acsr_mult_Bdns_is_Cdns(SparseMatrix&, SparseMatrix&, SparseMatrix&);
    void readCooFromFile(std::string);
    void printToFile(std::string,int);
    void InitializeSolve();
    void solve(SparseMatrix&, SparseMatrix&);
    void FinalizeSolve(int);
    MKL_INT iparm[64];
    MKL_INT mtype;
    MKL_INT nrhs;
    MKL_INT n;
    void *pt[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    //        void printASCII(std::string, int);
    static void testPardiso();


};
#endif





