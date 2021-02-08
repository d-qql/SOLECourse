cmake_minimum_required(VERSION 3.17)
project(SOLECourse)

set(CMAKE_CXX_STANDARD 20)
include_directories(/usr/boost/include)
link_directories(/usr/boost/include)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a -O3 -march=native")
add_executable(SOLECourse main.cpp Dense/DenseMatrix.h utility/Triplet.h utility/Consts.h Dense/GaussMethod.h Dense/HeadGaussMethod.h Generators/RandMatrixGenerator.h utility/Overloads.h Dense/ThomasAlgorithm.h Generators/RandVectorGenerator.h Generators/RandTridiagonalMatrixGen.h Dense/BackSubstitution.h Dense/LU.h Sparse/CSR.h Sparse/SimpleIteration.h Tests/Tests.h Sparse/Gauss-Seidel.h Dense/Householder.h Sparse/Yacobi.h)
target_link_libraries(SOLECourse -lboost_iostreams -lboost_system -lboost_filesystem)