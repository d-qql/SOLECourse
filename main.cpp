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

int main() {
    std::set<Triplet<double>> in;
    for(size_t i = 0; i < 300; ++i){
            in.insert({double (1 + i/598.), i, i});
    }
    CSR<double> A = CSR<double>(300, 300, in);
    std::cout<<A;
    std::vector<double> b = GenerateVector<double>(300, -1, 1);
    std::vector<double> roots = ChebyshevRoots<double>({1, 1.5}, 7);
    SimpleIteration(A, b, 2./2.5);
    SimpleIteration(A, b, 1.);
    FastSimpleIteration(A, b, roots);


    Gnuplot gp;
    gp<<"set xlabel 'Номер итерации' \n"
        "set ylabel 'ln(|r|))'\n"
        "set grid\n"
        "set title 'Зависимость логарифма модуля невязки от номера итерации ПМИ' font 'Helvetica Bold, 10'\n";
    gp << "plot '../PlotData/SimpleIteration/1SimpleIterNorm.txt' u 1:2 with lines title 'tao = 1' lc rgb 'black', "
          "     '../PlotData/SimpleIteration/OptSimpleIterNorm.txt' u 1:2 with lines title 'Optimal' lc rgb 'purple',"
          "     '../PlotData/SimpleIteration/FastSimpleIterNorm.txt' u 1:2 with lines title 'Chebyshev' lc rgb 'green'\n";

    return 0;
}

