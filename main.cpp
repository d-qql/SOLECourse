#include <iostream>
//#include "utility/Triplet.h"
#include "Sparse/CSR.h"
#include "Sparse/SimpleIteration.h"
#include "utility/Overloads.h"
#include <set>
#include "Dense/DenseMatrix.h"
#include "Dense/GaussMethod.h"
#include "Generators/RandMatrixGenerator.h"
#include "Dense/HeadGaussMethod.h"
#include "utility/gnuplot-iostream.h"
#include "Generators/RandVectorGenerator.h"
#include "Generators/RandTridiagonalMatrixGen.h"
#include "Dense/ThomasAlgorithm.h"
#include "Dense/LU.h"
#include "Tests/Tests.h"
#include <ctime>
#include <fstream>

int main() {

    /*std::set<Triplet<int>> in = {{1, 0, 0}, {2, 0, 1}, {3, 0, 3}, {4, 1, 2}};
    CSR<int> A = CSR<int>(3, 4, in);
    A.print();
    std::set<Triplet<int>> in1 = {{1, 0, 0}, {2, 0, 1}, {4, 1, 1}, {2, 2, 1}, {6, 2, 2}};
    CSR<int> K = CSR<int>(3, 3, in1);
    K.print();
    */
    testSimpleIteration_Speed();
    return 0;
}

