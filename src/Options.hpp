#ifndef __OPTIONS_HPP_INCLUDED
#define __OPTIONS_HPP_INCLUDED
#include <string>
#include <iostream>

class SolverOptions{
public:
    int solver; // 0 - pardiso, 1 - dissection
};


class CubeMesh{
public:
    int N[3];
    int n[3];
};

class Options {
    public:
        int print_matrices;
        std::string  path2data;
        void set_values(std::string, int, char*[], double, double, int,int);
        CubeMesh meshSetting;
        double young_modulus;
        double poissons_ratio;
        SolverOptions solver_opt;
};
#endif
