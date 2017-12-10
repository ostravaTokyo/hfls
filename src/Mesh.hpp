#ifndef __MESH_HPP_INCLUDED
#define __MESH_HPP_INCLUDED
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include "math.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "Options.hpp"

#include "metis.h"


using namespace std;

class Element
{
public:
   Element(){}
   void setInd(int, int);
   int getInd(int);
   void setPartitionId(int);
   void setMaterialId(int);
   int getPartitionId(){return PartitionId;}
   int getMaterialId(){return MaterialId;}
private:
   int ind[8];
   int PartitionId;
   int MaterialId;

};

class Point
{
public:
    Point(double x_,double y_,double z_):x(x_),y(y_),z(z_){}
    void set(double,double,double);
    void set_x(double );
    void set_y(double );
    void set_z(double );
    double get_x(){return x;}
    double get_y(){return y;}
    double get_z(){return z;}
private:
   double x;
   double y;
   double z;
};

class Material
{
public:
    double young_modulus;
    double poissons_ratio;
};

class Mesh{
public:
    Mesh();
    ~Mesh();

    void allocatePoints(int);
    void allocateElements(int);

    int nElementsClst;
    int nSubClst;
    vector < int > nElementsSub;
    int nPoints;
    int nSub;
    vector < Element> elements;
    vector < Point  > points;
    vector < int > cornerNodes;
    vector < int > cornerDOFs;
    vector < int > DirichletDOFs;
    void createMesh(map< string, string> & );
    Material material;
    void SaveVTK(vector < double > ,string, int);
    void ddm_metis(map <string, string> &);
    int nElSubXYZ[3];
    int nSubXYZ[3];
};




#endif // MESH_H
