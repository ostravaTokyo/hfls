cmake_minimum_required(VERSION 2.8) 
PROJECT(htfeti CXX) 


#SET (METIS_LIBRARIES 
#    "/data_space/france/metis-5.1.0/build/Linux-x86_64/libmetis/libmetis.a" 
#    CACHE PATH "mydoc")
#SET (METIS_INCLUDEDIRS 
#    "/data_space/france/metis-5.1.0/include/" 
#    CACHE PATH "mydoc")


#configure_file (src/Array1d_double64.hpp src/Options.hpp src/Cluster.hpp src/Matrix.hpp)


include_directories(src)

#htfeti.cpp
#Array1d_double64.cpp
#Options.cpp
#Cluster.cpp
#Matrix.cpp

#include_directories(${METIS_INCLUDEDIRS} ${CMAKE_CURRENT_SOURCE_DIR})


#find_package(VTK 7.1.0 REQUIRED NO_MODULE)
#include(${VTK_USE_FILE})
add_executable(htfeti src/htfeti.cpp)


#target_link_libraries(metisVTK ${VTK_LIBRARIES} ${METIS_LIBRARIES})
