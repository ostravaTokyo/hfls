#include "Options.hpp"
//
void Options::set_values(int n_subdomOnCluster_, std::string path2data_)
{
//  std::cout << "in Option.cpp: n_subdomOnCluster = " << n_subdomOnCluster_ << std::endl;
    n_subdomOnCluster = n_subdomOnCluster_;
    path2data = path2data_;

    meshSetting.N[0] = 2;
    meshSetting.N[1] = 2;
    meshSetting.N[2] = 2;
    meshSetting.n[0] = 2;
    meshSetting.n[1] = 2;
    meshSetting.n[2] = 2;
}
