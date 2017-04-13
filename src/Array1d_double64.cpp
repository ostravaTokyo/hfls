#include "Array1d_double64.h"



Array1d_double64::Array1d_double64() 
{
}


void Array1d_double64::readCooFromFile( std::string path2matrix,int n_val)
{
    using namespace std;

    std::ifstream input(path2matrix.c_str());
    std::cout << path2matrix.c_str() <<"\n";
    std::cout << "n_val = " << n_val<< "\n";

    this->val.resize(n_val);
    double double1;
    for (int i = 0; i < n_val; i++) {
        input >> double1;
        this->val[i] = double1;
    }
    std::cout << "f.size() = " << this->val.size() << "\n";

}
