#ifndef __CLUSTER_HPP_INCLUDED
#define __CLUSTER_HPP_INCLUDED

#include "vector"
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
    /* standard data */
//    std::vector< Matrix> K;
//    std::vector< Matrix> K_reg;
//    std::vector< Matrix> Fc;
//    std::vector< Matrix> Lumped;
//    std::vector< Array1d_double64 > rhs;


    /* 'new' data created inside the app */
    std::vector< Matrix> K_new;
    std::vector< Matrix> K_reg_new;
    std::vector< Matrix> Lumped_new;
    std::vector< Matrix> R_new;
    std::vector< Matrix> Fc_new;
    std::vector< Matrix> Bc_new;
    std::vector< Matrix> Bc_dense_new;
    std::vector< Matrix> Gc_new;

    /* FETI (first ...) */
//    std::vector< Matrix> R;
//    std::vector< Matrix> Bc;
//    std::vector< Matrix> Bc_dense;
//    std::vector< Matrix> Bf;
//    /* FETI  (second ...)*/
//    std::vector< Matrix> Gc;
//    std::vector< Matrix> Gf;

    Mesh mesh;
    Data data;






//    Matrix Gf_clust;

    Matrix Fc_clust_new;
    Matrix Gc_clust_new;
    Matrix Ac_clust_new;


    /* FETI */
    void create_Gf_clust_new();

    /* HFETI */
    void create_clust_object(Matrix &, vector <Matrix> & , bool);
    void create_Fc_clust_new();
    void create_Gc_clust_new();
    void create_Ac_clust_new();

    /* constraints */
    void create_cluster_constraints(const Options &);

};

#endif
