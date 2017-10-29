#PROJ_DIR = $(HOME)/Desktop/hfls_app_
PROJ_DIR = /data_space/WorkSpace/htfeti_app/gitRepo/hfls
#DISSECTION_INC = /data_space/france/dissection/include
DISSECTION_INC = /data_space/france/dissection_20170530/include


#CXX	= icpc -std=c++11
#CXX	= g++ -std=c++11
CXX	= g++ -std=c++98



UNAME := $(shell uname -s)

#CXX = g++ -std=c++0x



CCFLAGS   = -c -g0 -DDISSECTION -DUSE_PARDISO=1 -DVERBOSE_LEVEL=3 -I$(PROJ_DIR)/include -I$(DISSECTION_INC)
LD      = $(CXX)



METIS_DIR = /data_space/france
METIS_INC =  $(METIS_DIR)/metis-5.1.0/include
LIB_DIR_METIS =  $(METIS_DIR)/metis-5.1.0/build/Linux-x86_64/libmetis
CCFLAGS += -I$(METIS_INC)

CCFLAGS += -I$(METIS_INC)
ifeq ($(UNAME), Linux)
    LIB_METIS = -Xlinker -rpath=$(LIB_DIR_METIS) -L$(LIB_DIR_METIS) -lmetis
endif


#LIB_DIR_MKL = /opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64 
LIB_DIR_MKL = /opt/intel/composer_xe_2015/mkl/lib/ 
#SHARED_COMPILER = -L$(LIB_DIR_MKL) -lmkl_intel_lp64 -lmkl_sequential \
#                  -lmkl_core \
#                  -L. -lDissection

#SHARED_COMPILER = -L$(LIB_DIR_MKL) -lmkl_intel_lp64 -lmkl_sequential \
#                  -lmkl_core \
#                  -Xlinker -rpath=/data_space/france/dissection/lib -L. -lDissection

SHARED_COMPILER = -L$(LIB_DIR_MKL) -lmkl_intel_lp64 -lmkl_sequential \
                  -lmkl_core \
                  -Xlinker -rpath=/data_space/WorkSpace/htfeti_app/gitRepo/hfls -L. -lDissection


