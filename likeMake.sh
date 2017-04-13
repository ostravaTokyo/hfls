#/bin/bash



filecpp="src/Array1d_double64.cpp src/Cluster.cpp src/SparseMatrix.cpp   src/htfeti.cpp" 


#icpc -std=gnu++98 -fPIC -g -O3 -I/data_space/france/dissection/src -DNO_TO_STRING $filecpp -o htfeti.so \
#    -L/data_space/france/dissection/lib/ -lDissection \
#    -L/opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64 \
#    -lmkl_intel_lp64 -lmkl_sequential -lmkl_core  \
#    -L/opt/intel/compilers_and_libraries_2017/linux/compiler/lib/intel64 -lpthread 


#compiler=icpc
compiler=g++
$compiler -std=gnu++11 -fPIC -g -O3 $filecpp -o htfeti.so \
    -L/opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64 \
    -lmkl_intel_lp64 -lmkl_sequential -lmkl_core  \
    -L/opt/intel/compilers_and_libraries_2017/linux/compiler/lib/intel64 -lpthread 
