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
#include "Sparse/SSOR.h"
#include "Sparse/GradientDescent.h"
#include "Krylov/Krylov.h"
#include "Sparse/CG.h"
int main() {
    std::set<Triplet<double>> in;
    for(size_t i = 0; i < 300; ++i){
            in.insert({1. + double(i)/33.22, i, i});
    }
    CSR<double> A = CSR<double>(300, 300, in);

    std::vector<double> b = GenerateVector<double>(300, -1, 1);

    std::cout<<CGmethod(A, b)<<GradientDescent(A, b)<<FastSimpleIteration(A, b, ChebyshevRoots<double>({1, 10}, 7));
    Gnuplot gp;
    gp<<"set xlabel 'Номер итерации' \n"
        "set ylabel 'Логарифм невязки'\n"
        "set grid\n"
        "set title 'Зависимость логарифма невязки от номера итерации. Обусловленность 10' font 'Helvetica Bold, 10'\n";
    gp << "plot '../PlotData/CG/10.txt' with lines title 'CG' lc rgb 'blue',"
          "     '../PlotData/SimpleIteration/10.txt' with lines title 'Cheb. Simple Iter' lc rgb 'green',"
          "     '../PlotData/GradientDescent/10.txt' with lines title 'GD' lc rgb 'red'\n";

    return 0;
}

