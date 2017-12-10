#ifndef __DATA_HPP_INCLUDED
#define __DATA_HPP_INCLUDED
#include "Mesh.hpp"
#include "Matrix.hpp"
#include <map>


class local_K_f
{
public:
    void set_val_K(int, int,double);
    void set_val_f(int,double);
    double get_val_K(int i,int j){return val_K[i + j * 24];}
    double get_val_f(int i){return val_f[i];}
    void setIeq(int i,int j){ieq[i] = j;}
    int getIeq(int i){return ieq[i];}
    int get_nDOF(){return nDOF;}
    int set_nDOF(int n){nDOF = n;}
    void setZero_val_K();
    void setZero_val_f();
private:
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


    void fe_assemb_local_K_f(Mesh &, map <string, string> &);

    void feti_symbolic(Mesh &, vector <Matrix> &);
    void feti_numeric(Mesh &, vector <Matrix> &, vector <Vector> &);
    void feti_numeric_element(Matrix &, Vector &, local_K_f &);
    void stf_mtrx_solid45(local_K_f &, vector<Point> , double ,double);
    void create_analytic_ker_K(Mesh &, vector <Matrix> &);
    static double inverse_matrix_3x3(double *, double *);
    void resize_local_K_f_clust(int);

    vector < map<int,int> > g2l;
    vector < map < int,vector < int > > > interface;
    vector < vector < Interfaces > > interfaces;
    void buildMappingStruct(Mesh&);
    vector < vector <int> > l2g;


private:
    vector < local_K_f >  local_K_f_clust;
    vector < vector <int> > selectorOfElemPartitId;

};

#endif // __DATA_HPP_INCLUDED
