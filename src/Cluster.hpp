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
    std::vector< Matrix> rhs;
    std::vector< Matrix> Lumped;
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


    /* dual */

    /* constraints */
    vector < vector < int > > neighbours;

    void create_cluster_constraints(const Options &);
    void create_Bc_or_Bf_in_CSR(vector <Matrix> &, bool , bool);
    void create_Bc_or_Bf_in_CSR(vector <Matrix> &, bool , bool, bool);
    void create_Bc_weightedAverages_in_COO(vector <Matrix> &, bool);
    void create_Rf_and_Gf();

    void matrix_Bx_COO2CSR(vector <Matrix> &, int);

    void mult_Kplus_f(vector < Matrix > & , vector < Matrix > &);
    void mult_Bf(vector < Matrix > const &, Matrix &);
    void mult_BfT(Matrix const &, vector < Matrix > &);
    void mult_Gf(Matrix const &, Matrix &);
    void mult_GfT(Matrix const &, Matrix &);
    void mult_RfT(vector < Matrix > const &, Matrix &);
    void Project(Matrix const & , Matrix &, Matrix &);

    void mult_Ff(Matrix const &, Matrix &);

    void pcpg();

//    static double dot(Matrix &, Matrix &);
//    static double apb(Matrix &, Matrix &, double, double);


};

#endif
