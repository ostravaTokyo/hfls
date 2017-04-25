#ifndef __CLUSTER_HPP_INCLUDED
#define __CLUSTER_HPP_INCLUDED

#include "vector"
#include "Options.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include "Matrix.hpp"
#include "Array1d_double64.hpp"


class Cluster 
{
public:
    Cluster(Options options);
    /* standard data */
    std::vector< Matrix> K;
    std::vector< Array1d_double64 > rhs;

    /* FETI (first ...) */
    std::vector< Matrix> R;
    std::vector< Matrix> K_regularized;
    std::vector< Matrix> Bc;
    std::vector< Matrix> Bct;
    std::vector< Matrix> Bf;
    std::vector< Matrix> RegMat;
    /* FETI  (second ...)*/
    std::vector< Matrix> Gc;
    std::vector< Matrix > Fc_sub;


};

#endif
