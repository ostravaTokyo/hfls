#ifndef __OPTIONS_HPP_INCLUDED
#define __OPTIONS_HPP_INCLUDED
#include <string>
#include <iostream>

class CubeMesh{
public:
    int N[3];
    int n[3];
};

class Options {
    public:
        std::string  path2data;
        void set_values(std::string, int, char*[], double, double);
        CubeMesh meshSetting;
        double young_modulus;
        double poissons_ratio;

};
#endif
