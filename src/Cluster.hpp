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
    std::vector< Matrix> K_new;
    std::vector< Matrix> K_reg_new;
    std::vector< Matrix> Lumped_new;
    std::vector< Matrix> R_new;
    std::vector< Matrix> Fc_new;
    std::vector< Matrix> Bc_new;
    std::vector< Matrix> Bc_dense_new;
    std::vector< Matrix> Gc_new;


    Mesh mesh;
    Data data;

    Matrix Fc_clust_new;
    Matrix Gc_clust_new;
    Matrix Ac_clust_new;
    Matrix GcTGc_clust;
    Matrix kerGc;

    Matrix GcTGc_sparse_clust;

    /* FETI */
    void create_Gf_clust_new();

    /* HFETI */
    void create_clust_object(Matrix &, vector <Matrix> & , bool);
    void create_Fc_clust_new();
    void create_Gc_clust_new();
    void create_Ac_clust_new(bool);
    void create_GcTGc();
    void create_GcTGc_clust_sparse();


    /* dual */
//    int n_inerf_c;
    int n_inerf_c_max;

    /* constraints */
    void create_cluster_constraints(const Options &);
    vector < vector < int > > neighbours;
};

#endif
