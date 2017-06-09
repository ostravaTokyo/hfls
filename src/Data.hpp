#ifndef __DATA_HPP_INCLUDED
#define __DATA_HPP_INCLUDED
#include "Mesh.hpp"
#include "Matrix.hpp"
#include <map>


class local_K_f
{
public:
   double val_K[576];
   double val_f[24];
   int    ieq[24];
   int    nDOF;
};


class Interfaces
{
public:
    int IdNeighSub;
    vector <int> dofs;
};


class Data
{
public:
    Data();
    ~Data();


    void fe_assemb_local_K_f(Mesh &);

    void feti_symbolic(Mesh &, vector <Matrix> &);
    void feti_numeric(Mesh &, vector <Matrix> &);
    void feti_numeric_element(Matrix &, local_K_f &);
    void stf_mtrx_solid45(local_K_f &, Point *, double ,double);
    static double inverse_matrix_3x3(double *, double *);

    vector < local_K_f >  local_K_f_clust;
    vector < vector <int> >  l2g;
    vector < map<int,int> > g2l;
    vector < vector <int> > selectorOfElemPartitId;
    vector < map < int,vector < int > > > interface;
    vector < vector < Interfaces > > interfaces;




    void create_analytic_ker_K(Mesh &, vector <Matrix> &);

};

#endif // __DATA_HPP_INCLUDED
