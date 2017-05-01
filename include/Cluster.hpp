#ifndef __CLUSTER_HPP_INCLUDED
#define __CLUSTER_HPP_INCLUDED

#include "vector"
#include "Options.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include "Matrix.hpp"
#include "Array1d_double64.hpp"


//struct int_int_dbl
//{
//    int I;
//    int J;
//    double V;
//};


class Cluster 
{
public:
    Cluster(Options options);
    /* standard data */
    std::vector< Matrix> K;
    std::vector< Matrix> K_reg;
    std::vector< Matrix> Fc;
    std::vector< Matrix> Lumped;
    std::vector< Array1d_double64 > rhs;

    /* FETI (first ...) */
    std::vector< Matrix> R;
    std::vector< Matrix> Bc;
    std::vector< Matrix> Bc_dense;
    std::vector< Matrix> Bf;
    /* FETI  (second ...)*/
    std::vector< Matrix> Gc;


    Matrix Fc_clust;

    void createFc_clust();

};

#endif
