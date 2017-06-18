#ifndef __MESH_HPP_INCLUDED
#define __MESH_HPP_INCLUDED
#include <string>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include "math.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "Options.hpp"


using namespace std;

class Element
{
public:
   int ind[8];
   int PartitionId;
   int MaterialId;

};

class Point
{
public:
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

    int nElementsClst;
    int nSubClst;
    vector < int > nElementsSub;
    int nPoints;
    int nSub;
    vector < Element > elements;
    vector < Point  > points;
    vector < int > cornerNodes;
    vector < int > cornerDOFs;
    vector < int > DirichletDOFs;
    void createMesh(const Options &);
    Material material;
    void SaveVTK(vector < double > ,string, int);
};




#endif // MESH_H
