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
#include "Dense/Householder.h"
#include "Sparse/Gauss-Seidel.h"
#include <ctime>
#include <fstream>
#include "Sparse/Yacobi.h"
#include "Chebyshev/Chebyshev.h"
#include "Sparse/SOR.h"
int main() {
    CSR<double> A = CSR<double>(300, 300, GenerateMatrixDiagDominant<double>(300));
    std::vector<double> b = GenerateVector<double>(300, -1, 1);
    std::pair<double, double> k = A.localizeEigenVals();
    double M = k.second/k.first;
    double w = 1 + pow(M/(1+sqrt(1-M*M)), 2);
    std::cout<<GaussSeidel(A, b)<<SOR(A, b, w);

    return 0;
}

