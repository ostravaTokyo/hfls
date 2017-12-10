#ifndef __SOLVER_HPP_INCLUDED
#define __SOLVER_HPP_INCLUDED

#include <map>
#include "Matrix.hpp"
#include <string>
#include "Cluster.hpp"
#include <ctime>



class Solver
{
public:
    Solver();
    void pcpg(map <string,string>, Cluster&);
private:
    double time_solver;
    double time_total;
};





#endif
