#ifndef __DRIVER_HPP_INCLUDED
#define __DRIVER_HPP_INCLUDED


#include "Mesh.hpp"
#include "map"
#include <string>
#include "Cluster.hpp"
#include "Solver.hpp"

class Driver
{

public:
    Driver();
    Driver(map <string, string>);
    void initialization();
private:
    Mesh mesh;
    Cluster cluster;
    map <string, string> options2;
    string folder;
    Solver solver;
};

#endif
