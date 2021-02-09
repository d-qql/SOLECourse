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

int main() {

    DenseMatrix<double> A = DenseMatrix<double>(10, 10, GenerateMatrix<double>(10, -100, 100));
    std::cout<<LUdecomp(A).first<<LUdecomp(A).second;
/*    CSR<double> A = CSR<double>(100, 100, GenerateMatrixDiagDominant<double>(100));
    std::vector<double> b = GenerateVector<double>(100, -1, 1);
    printSystem(A, b);
    std::cout<<GaussSeidel(A, b);
    std::cout<<Yacobi(A, b);
    std::cout<<SimpleIteration(A, b, 1.);
    Gnuplot gp;
    gp<<"set xlabel 'Номер итерации' \n"
        "set ylabel 'Логарифм невязки'\n"
        "set grid\n"
        "set title 'Зависимость логарифма нормы невязки от номера итерации' font 'Helvetica Bold, 10'\n";
    gp << "plot '../PlotData/SimpleIteration/YacobiNorm.txt' with lines title 'Якоби' lc rgb 'blue',"
          "     '../PlotData/SimpleIteration/SimpleIterNorm.txt' with lines title 'ПМИ' lc rgb 'green',"
          "     '../PlotData/SimpleIteration/GaussSeidelNorm.txt' with lines title 'Гаусс-Зейдель' lc rgb 'red'\n";*/
    return 0;
}

