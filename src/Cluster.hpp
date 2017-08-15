#ifndef __CLUSTER_HPP_INCLUDED
#define __CLUSTER_HPP_INCLUDED

#include "vector"
#include "map"
#include "Options.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include "Mesh.hpp"
#include "Matrix.hpp"
#include "Array1d_double64.hpp"
#include "Data.hpp"

class Cluster 
{

public:
    Cluster(Options options);
    /* 'new' data created inside the app */
    int nSubClst;
    int neqClst;
    int nLam_c;
    int nLam_f;
    int nRBM_f;
    int nRBM_c;
    string folder;
    std::vector< Matrix> K;
    std::vector< Vector > rhs;
    std::vector< Matrix> Preconditioner;
    std::vector< Matrix> R;
    std::vector< Matrix> Fc;
    std::vector< Matrix> Bc;
    std::vector< Matrix> KplusBcT;
    std::vector< Matrix> Bf;
    std::vector< Matrix> BcT_dense;
    std::vector< Matrix> Gc;
    std::vector< Matrix > Rf;
    std::vector< Matrix > Gf;


    Mesh mesh;
    Data data;

    Matrix Fc_clust;
    Matrix Gc_clust;
    Matrix Ac_clust;
    Matrix Gf_clust;
    Matrix GcTGc_clust;
    Matrix kerGc;

    Matrix GcTGc_sparse_clust;
    Matrix GfTGf;

    Matrix invGfTGf;

    Options options;

    /* FETI */
    void create_Gf_clust();

    /* HFETI */
    void create_clust_object(Matrix &, vector <Matrix> & , bool);
    void create_Fc_clust();
    void create_Gc_clust();
    void create_Ac_clust(bool);
    void create_GcTGc();
    void create_GcTGc_clust_sparse();
    void create_GfTGf();
    void compute_invGfTGf();

    vector < double > weigth;

    /* dual */

    /* constraints */
    vector < vector < int > > neighbours;

    void create_cluster_constraints(const Options &);
    void create_Bc_or_Bf_in_CSR(vector <Matrix> &, bool , bool);
    void create_Bc_or_Bf_in_CSR(vector <Matrix> &, bool , bool, bool);
    void create_Bc_weightedAverages_in_COO(vector <Matrix> &, bool);
    void create_Rf_and_Gf();

    void matrix_Bx_COO2CSR(vector <Matrix> &, int);

    void mult_Kplus_f(vector < Vector > & , vector < Vector > &);
    void mult_Bf(vector < Vector > const &, Vector &);
    void mult_BfT(Vector  const &, vector < Vector > &);
    void mult_Gf(Vector const &, Vector &);
    void mult_GfT(Vector const &, Vector &);
    void mult_RfT(vector < Vector > const &, Vector&);
    void Projection(Vector const & , Vector &, Vector &);

    void mult_Ff(Vector const &, Vector &);
    void Preconditioning(Vector const &, Vector &);
    void scale(Vector &);

    void pcpg();
    void pcpg2();
    void printVTK(vector < Vector > &, vector < Vector > &, Vector &, Vector &, int);

//    static double dot(Matrix &, Matrix &);
//    static double apb(Matrix &, Matrix &, double, double);


};

#endif
