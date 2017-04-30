#CXX	= icpc -std=c++11
CXX	= g++ -std=c++11
CCFLAGS   = -g0 -DUSE_PARDISO=1
LD      = $(CXX)

LIB_DIR_MKL = /opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64 
SHARED_COMPILER =  -Xlinker -rpath=$(LIB_DIR_MKL) -L$(LIB_DIR_MKL) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core 
