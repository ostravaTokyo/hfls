#include "../include/Options.hpp"
//
void Options::set_values(int _n_subdomOnCluster, std::string _path2data)
{
    std::cout << "in Option.cpp: n_subdomOnCluster = " << _n_subdomOnCluster << std::endl;
    this->n_subdomOnCluster = _n_subdomOnCluster;
    this->path2data = _path2data;
}
