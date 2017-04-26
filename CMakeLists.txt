cmake_minimum_required(VERSION 2.8)
PROJECT(htfeti C CXX)
#
#file(GLOB SRC_FILES RELATIVE "src" "*.cpp")
#add_executable(htfeti ${SRC_FILES})



file(GLOB htfeti_SRC
    "htfeti.cpp"
    "include/*.hpp"
    "src/*.cpp"
)

add_executable(htfeti ${htfeti_SRC})
