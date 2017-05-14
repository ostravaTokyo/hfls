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
    void createMesh();
};




#endif // MESH_H
