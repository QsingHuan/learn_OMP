cmake_minimum_required(VERSION 3.20)
project(learn_OMP)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
message(${CMAKE_CXX_FLAGS})
find_package(OpenMP REQUIRED)

add_executable(learn_OMP main.cpp)
target_link_libraries(learn_OMP OpenMP::OpenMP_CXX)