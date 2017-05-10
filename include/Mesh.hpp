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

struct Element
{
   int ind[8];
   int PartitionId;
   int MaterialId;
   int nSub;

};

struct Point
{
   double x;
   double y;
   double z;
};


class Mesh{
public:
    Mesh();
    ~Mesh();

    int nElements;
    int nPoints;
    int nSub;
    vector < Element > elements;
    vector < Point  > points;
    void createMesh();
};




#endif // MESH_H
