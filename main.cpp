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
#include "Givens/Rotations.h"
#include "GMRES/GMRES.h"
int main() {
    std::set<Triplet<double>> in;
   /* size_t p, r;
    double v;
    std::fstream fin("/home/d-qql/CLionProjects/SOLECourse/cmake-build-debug/MT.txt");
    fin>>p>>r>>v;
    for(int i = 0; i < 1069; ++i){
        fin>>p>>r>>v;
        in.insert({v, p, r});
    }
    std::cout<<in.size();
*/
    size_t SZ = 50;

    std::vector<double> b;

        b = GenerateVector<double>(SZ, -1, 1);
        in = GenerateMatrixDiagDominant<double>(SZ);


        CSR<double> A = CSR<double>(SZ, SZ, in);
        std::cout<<GMRES(A, b)<<GaussSeidel(A, b);
    Gnuplot gp;
    gp<<"set xlabel 'Номер итерации' \n"
        "set ylabel 'Логарифм нормы невязки'\n"
        "set grid\n"
        "set title 'Зависимость логарифма нормы невязки от номера итерации' font 'Helvetica Bold, 10'\n";
    gp << "plot '../PlotData/SimpleIteration/GaussSeidelNorm.txt' with lines title 'GaussSeidel' lc rgb 'blue',"
          "     '../PlotData/GMRES/Norm.txt' with lines title 'GMRES' lc rgb 'red'\n";

    //std::vector<double> b = GenerateVector<double>(5, -1, 1);
    //std::vector<double> v = GenerateVector<double>(10, -10000, 10000);


    /*auto [a, c] = KrylovSubSpace(A, b, v, 10);

    std::cout<<c;

    std::cout<<A;
    for(int i = 0; i < 10; ++i){
        for(int j = 0; j < 9; ++j){
            std::cout<<a[i]*(A*a[j])<<" "<<c(i, j)<<std::endl;
        }
    }
    for(auto i : a) std::cout<<i;
     */
    return 0;
}

