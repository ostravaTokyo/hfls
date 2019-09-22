PROJ_DIR = /home/mar440/WorkSpace/hfls
CXX	= g++ -std=c++11

UNAME := $(shell uname -s) 
CCFLAGS   = -c -g0  -DUSE_PARDISO=1 -DVERBOSE_LEVEL=3 -I$(PROJ_DIR)/include
LD      = $(CXX) 

METIS_DIR = /home/mar440/Metis
METIS_INC =  $(METIS_DIR)/include
LIB_DIR_METIS =  $(METIS_DIR)/build/libmetis
CCFLAGS += -I$(METIS_INC) 

ifeq ($(UNAME), Linux)
    LIB_METIS = -Xlinker -rpath=$(LIB_DIR_METIS) -L$(LIB_DIR_METIS) -lmetis
endif

LIB_DIR_MKL = /home/mar440/intel/mkl/lib/

#SHARED_COMPILER = -L$(LIB_DIR_MKL) -lmkl_intel_lp64 -lmkl_sequential \
#                  -lmkl_core \
#                  -Xlinker -rpath=/home/mar440/WorkSpace/hfls
SHARED_COMPILER = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 \
                    -lmkl_sequential -lmkl_core -lpthread -lm -ldl



