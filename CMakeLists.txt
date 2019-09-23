cmake_minimum_required(VERSION 2.8)
PROJECT(htfeti C CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/")
MESSAGE( STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH} )
set(MKLROOT "/home/mar440/intel/" CACHE PATH "MKL Directory" FORCE)

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_BUILD_TYPE Debug)
 
file(GLOB htfeti_SRC "htfeti.cpp" "src/*.cpp")

SET (METIS_LIBRARIES
    "/home/mar440/Metis/build/libmetis/libmetis.a" 
    CACHE PATH "mydoc")
SET (METIS_INCLUDEDIRS 
    "/home/mar440/Metis/include/" 
    CACHE PATH "mydoc")

#add_definitions(-DNOT_METIS)
#add_definitions(-DNOT_MKL) 

if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -fPIC")
endif() 

find_package( MKL REQUIRED ) 
find_package(Threads)

add_definitions(-DHAVE_MKL)
add_definitions(-DUSE_PARDISO)
set(HAVE_MKL 1)
include_directories(${MKL_INCLUDE_DIR})
link_directories(${MKL_LIBRARY_DIR})

include_directories(${METIS_INCLUDEDIRS} ${CMAKE_CURRENT_SOURCE_DIR})

set(BLAS_VENDOR Intel10_64lp) 
find_package(BLAS REQUIRED) 
find_package(LAPACK REQUIRED)

set(MKL_COMPILER_LIB_DIR ${MKL_COMPILER_LIB_DIR} m dl pthread mkl_core mkl_sequential)
set(MKL_LIB_DIR ${MKL_LIB_DIR} m dl mkl_sequential mkl_core pthread)
 

add_executable(${PROJECT_NAME} ${htfeti_SRC})

target_link_libraries(${PROJECT_NAME} ${METIS_LIBRARIES} )
target_link_libraries(${PROJECT_NAME} ${BLAS_LIBRARIES})
