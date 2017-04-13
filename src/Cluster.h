#ifndef __CLUSTER_H_INCLUDED
#define __CLUSTER_H_INCLUDED

#include "Options.h"
#include "vector"
#include "SparseMatrix.h"
#include "Array1d_double64.h"
#include <iostream>
#include <string>
#include <fstream>


class Cluster 
{
public:
    Cluster(Options options);
    /* standard data */
    std::vector< SparseMatrix> K;
    std::vector< Array1d_double64 > rhs;

    /* FETI (first ...) */
    std::vector< SparseMatrix> R;
    std::vector< SparseMatrix> K_regularized;
    std::vector< SparseMatrix> Bc;
    std::vector< SparseMatrix> Bct;
    std::vector< SparseMatrix> Bf;
    std::vector< SparseMatrix> RegMat;
    /* FETI  (second ...)*/
    std::vector< SparseMatrix> Gc;
    std::vector< SparseMatrix > Fc_sub;


};

#endif
