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

    std::vector< Matrix> K;
    std::vector< Matrix> K_reg;
    std::vector< Matrix> Lumped;
    std::vector< Matrix> R;
    std::vector< Matrix> Fc;
    std::vector< Matrix> Bc;
    std::vector< Matrix> BcT_dense;
    std::vector< Matrix> Gc;


    Mesh mesh;
    Data data;

    Matrix Fc_clust;
    Matrix Gc_clust;
    Matrix Ac_clust;
    Matrix GcTGc_clust;
    Matrix kerGc;

    Matrix GcTGc_sparse_clust;

    /* FETI */
    void create_Gf_clust();

    /* HFETI */
    void create_clust_object(Matrix &, vector <Matrix> & , bool);
    void create_Fc_clust();
    void create_Gc_clust();
    void create_Ac_clust(bool);
    void create_GcTGc();
    void create_GcTGc_clust_sparse();


    /* dual */
//    int n_inerf_c;
    int n_interf_c_max;

    /* constraints */
    void create_cluster_constraints(const Options &);

    void create_Bc_or_Bf(vector <Matrix> &, int);
    void create_Bc_weightedAverages(vector <Matrix> &);

    void matrix_Bx_COO2CSR(vector <Matrix> &, int);
    vector < vector < int > > neighbours;
};

#endif
