#include "Array1d_double64.hpp"



Array1d_double64::Array1d_double64() 
{
}


void Array1d_double64::readCooFromFile( std::string path2matrix,int n_val)
{
    using namespace std;

    std::ifstream input(path2matrix.c_str());
#if VERBOSE_LEVEL>3
    std::cout << path2matrix.c_str() <<"\n";
    std::cout << "n_val = " << n_val<< "\n";
    std::cout << "f.size() = " << val.size() << "\n";
#endif

    val.resize(n_val);
    double double1;
    for (int i = 0; i < n_val; i++) {
        input >> double1;
        val[i] = double1;
    }
}
