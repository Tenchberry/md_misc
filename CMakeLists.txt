cmake_minimum_required(VERSION 3.7)
project(mcmc)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES vVerlet_LJ.cpp)
add_executable(mcmc ${SOURCE_FILES})