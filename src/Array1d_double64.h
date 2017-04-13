#ifndef __ARRAY1D_DOUBLE64_H_INCLUDED
#define __ARRAY1D_DOUBLE64_H_INCLUDED
#include "mkl.h"
#include <vector>
#include <string>
#include <stdio.h>
#include <fstream> 
#include <iostream>

class Array1d_double64{

public:
    Array1d_double64();
    std::vector < double > val;
    void readCooFromFile(std::string, int);
};
#endif
