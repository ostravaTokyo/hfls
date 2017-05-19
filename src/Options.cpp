#include "Options.hpp"
//
void Options::set_values(std::string path2data_, int argc, char * argv[],
                           double young_modulus_, double poissons_ratio_)
{
    path2data = path2data_;

    meshSetting.N[0] = 2;
    meshSetting.N[1] = 2;
    meshSetting.N[2] = 2;
    meshSetting.n[0] = 2;
    meshSetting.n[1] = 2;
    meshSetting.n[2] = 2;

    young_modulus = young_modulus_;
    poissons_ratio = poissons_ratio_;



    if (argc > 1) {
       meshSetting.N[0] =  atoi(argv[1]);
    }
    if (argc > 2) {
       meshSetting.N[1] =  atoi(argv[2]);
    }
    if (argc > 3) {
       meshSetting.N[2] =  atoi(argv[3]);
    }
    if (argc > 4) {
       meshSetting.n[0] =  atoi(argv[4]);
    }
    if (argc > 5) {
       meshSetting.n[1] =  atoi(argv[5]);
    }
    if (argc > 6) {
       meshSetting.n[2] =  atoi(argv[6]);
    }





}
