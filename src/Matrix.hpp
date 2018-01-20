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
#ifndef WIN32 
#include <unistd.h>
#else
#include "stdint.h"
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include "omp.h"
#include <ctime>
#include <stack>
#include <time.h>


#ifdef DISSECTION
#include "Dissection.hpp"
#endif

//#include "Options.hpp"


//    symmetric = 0;  // 0-not, 1-lower tr., 2-upper tr.
//    format = 1;     // 0-coo, 1-csr, 2-dense
//    offset = 0;


using namespace std;

struct TRIPLET
{
    int I;
    int J;
    double V;
};

bool myfunction (int i,int j);

class Matrix{
public:
    Matrix();
    Matrix(int n_row_, int n_col_, bool NorT);
    Matrix(int n_row_, int n_col_, int nnz_, bool NorT);
    Matrix(string );
    ~Matrix();



    int get_n_row_cmprs();
    int get_n_row();
    int get_n_col();
    int get_numel();
    int get_order_number();
    int get_symmetric();
    int get_format();
    int get_nnz();

    void set_n_row_cmprs(int n){n_row_cmprs = n;}
    void set_n_row(int n){n_row = n;}
    void set_n_col(int n){n_col = n;}
    void set_order_number(int n){order_number = n;}
    void set_label(string s){label = s;}
    void set_symmetric(int i){symmetric = i;}
    void set_format(int i){format = i;}
    void set_nnz(int i){nnz = i;}

    void update_n_col(int n){n_col += n;}
    void update_n_row(int n){n_row += n;}
    void update_n_row_cmprs(int n){n_row_cmprs += n;}
    void update_nnz(int n){nnz += n;}


  //    void zero_dense(int);
    void init();
    void zero_dense(int,int);
    void zero_dense(int,int,bool);

    void mult(const Matrix &, Matrix &, bool);
    void mult(const double[], double [] , bool, int,int);
    void COO2CSR();
    void CSR2COO();
    void CSRorCOO2DNS(bool, bool);
    void CSRorCOO2DNS(bool, bool, int);
    void dummyFunction(int,int, int, double);
    void DNS2CSR();
    void DNS2COO();
    void compresRows();
    void readCooFromFile(string, int,int,int);
    void printToFile(string, string,int,bool);
    void FinalizeSolve(int);
    void setZero();
    void getNullPivots(vector < int > & );
    void sortAndUniqueCOO(vector < TRIPLET > &);
    void submatrix_row_selector(Matrix &, vector <int> &);
    void getSubBlockmatrix_rs(Matrix &, Matrix &, int, int, int, int);
    void getSubDiagBlockmatrix(const Matrix &, const Matrix &, Matrix &,int , int);

    void getBasicMatrixInfo();
    void getBasicMatrixInfo(string);


    static void testSolver(string,int);

    uint64_t *diss_dslv;
    int diss_num_threads;
    int diss_sym;
    int diss_decomposer;
    int diss_indefinite_flag;
    int diss_scaling;
    int diss_int0;
    int solver;
    string linear_solver;
    double diss_eps_pivot;


    void sym_factor();


    void sym_factor_pardiso();
    void num_factor_pardiso();

    void sym_factor_dissection();
    void num_factor_dissection();
    void num_factor_dissection(Matrix &);
    void num_factor_dissection(Matrix &, bool);


    void sym_factor(int);
    void sym_factor(string);
    void num_factor();
    void num_factor(Matrix&);
    void num_factor(Matrix&, bool);
    void solve_system_pardiso(Matrix&, Matrix&);


    void solve_system(Matrix &, Matrix &);

    void solve_system_dissection(Matrix const &, Matrix &);
    vector <int> nullPivots;


    void mat_mult_dense(Matrix const &, string, Matrix const &, string);


    double norm2();

    static bool cmp_int_int_I(TRIPLET ,TRIPLET );
    static bool cmp_int_int_J(TRIPLET ,TRIPLET );
    static bool compareDouble(double, double);
    static void updateCOOstructure(vector <TRIPLET > &, Matrix &, int, int);

    static bool test_of_Bc_constraints(Matrix &);

    static void print_int_vector(vector <double> &, int, int);

    static void getSingularVal_DNS(Matrix &, Matrix &, int);
    static void getSingularVal_DNS(Matrix , Matrix &, int, int);

    static void getEigVal_DNS(Matrix &, Matrix &, int );
    static void getEigVal_DNS(Matrix  , Matrix &, int, int);
    static void print1dArray(double [], int, string, string);

    static double dot(Matrix const &, Matrix const &);
    void add(Matrix &,double);
    void test_K_Kp_K_condition(Matrix &);


    void createDirichletPreconditioner(const Matrix & ,
                                       const Matrix & ,
                                       Matrix &,
                                       Matrix &,
                                       Matrix &,
                                       Matrix &
                                       );
    void getSubDiagBlockmatrix(const Matrix & , Matrix & , int , int );
    static void get_kernel_from_K(Matrix &,  Matrix &, int);
    static void get_kernel_from_K(Matrix &,  Matrix & );
    static void GramSchmidtOrtho(Matrix &);
    static double dot(double *, double*, int);

    static double getNorm_K_R(Matrix & , Matrix &);
    map < string, string > options2;




    vector < int >     l2g_i_coo;
    map <int,int>      g2l_i_coo;
    vector < double >  dense;
    vector < int >     j_col_cmpr;
    vector < int >     i_coo_cmpr;
    vector < int >     j_col; // new MKL_INT[crs->nnz];
    vector < double >  val;
    vector < int >     i_ptr; // new MKL_INT[crs->m+1];


    MKL_INT iparm[64];

protected:

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
    string tmp_label;
    int order_number;
    bool printed;


    int ij(int, int);

    MKL_INT m;
    MKL_INT mtype;
    MKL_INT nrhs;
    MKL_INT n;
    void *pt[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    MKL_INT n1_p;

    double  ddum;            /* Double dummy */
    MKL_INT idum;           /* Integer dummy. */


};


class Vector : public Matrix
{
public:
  Vector(int n_row = 0) : Matrix(n_row, 1, true) { }
  void zero_dense(int n_row)
  {
    bool NorT = true;
    int n_col = 1;
    Matrix::zero_dense(n_row, n_col, NorT);
  }
};
#endif


