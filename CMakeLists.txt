cmake_minimum_required(VERSION 2.8)
PROJECT(htfeti C CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/")
MESSAGE( STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH} )
set(MKLROOT "/opt/intel/compilers_and_libraries_2017/linux" CACHE PATH "MKL Directory" FORCE)


#
file(GLOB htfeti_SRC
    "htfeti.cpp"
    "src/*.cpp"
)


add_definitions(-DNOT_METIS)
add_definitions(-DNOT_MKL)


if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -fPIC")
endif()
# -lmkl_intel_lp64 -lmkl_sequential -lmkl_core


find_package( MKL REQUIRED )
add_definitions(-DHAVE_MKL)
set(HAVE_MKL 1)
include_directories(${MKL_INCLUDE_DIR})
link_directories(${MKL_LIBRARY_DIR})

# REMOVED BY ME
#add_definitions("/Qmkl:parallel")
# ADDED BY ME
set(MKL_COMPILER_LIB_DIR ${MKL_COMPILER_LIB_DIR} m dl pthread)
set(MKL_LIB_DIR ${MKL_LIB_DIR} m dl mkl_sequential mkl_core)




add_executable(htfeti ${htfeti_SRC})
