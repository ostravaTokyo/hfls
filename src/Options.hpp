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
        int n_subdomOnCluster;
        std::string  path2data;
        void set_values(int, std::string);
        CubeMesh meshSetting;

};
#endif
