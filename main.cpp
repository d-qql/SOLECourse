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
int main() {
    std::set<Triplet<double>> in;
    for(size_t i = 0; i < 200; ++i){
            in.insert({i+ 1., i, i});
    }
    CSR<double> A = CSR<double>(200, 200, in);
    std::vector<double> b = GenerateVector<double>(200, -10, 10);
    std::vector<double> x(200);
    x = b;
    x[0] = 1.;
    auto [Basis, H] = KrylovSubSpace(A, b, b, 200);
    for(const auto& i : Basis) std::cout<<i;
    DenseMatrix<double> C = DenseMatrix<double>(200,200, H);
    std::cout<<C;

    std::cout<<Basis[10] * (A * Basis[10])<< "\n"<<C(10, 10);
  /*  std::vector<std::string> color = {"red", "green", "yellow", "blue", "orange", "black"};
    Gnuplot gp;
    gp << "set xlabel 'Номер итерации' \n"
          "set ylabel 'Логарифм невязки'\n"
          "set grid\n"
          "set title 'Логарифм нормы невязки от номера итерации при разных числах обусловленности' font 'Helvetica Bold, 10'\n";
    std::set<Triplet<double>> in;
    std::vector<double> b(300);
    for(int k = 5; k < 31; k+=5) {
        for (int i = 0; i < 300; ++i){
            in.insert({ 1 + (i * (k-1)) / 300. , size_t(i), size_t(i)});
            b[i] = 300-i;
        }
        CSR<double> A = CSR<double>(300, 300, in);

        in.clear();
        std::cout<<GradientDescent(A,b, k)<<GaussSeidel(A,b)<<std::endl;
    }
    gp << "plot '../PlotData/GradientDescent/5.txt' u 1:2 with lines title 'Обусловленность: 5' lc rgb 'red'";
    for(int k = 10; k < 31; k+=5) gp << ", '../PlotData/GradientDescent/" + std::to_string(k) + ".txt' u 1:2 with lines title 'Обусловленность: " + std::to_string(k) + "' lc rgb '" + color[k/5-1]+"' ";
    gp << "\n";
*/
    /*Gnuplot gp;
    gp << "set xlabel 'Число обусловленности' \n"
          "set ylabel 'Количество итераций'\n"
          "set grid\n"
          "set title 'Зависимость числа итераций от числа обусловленности' font 'Helvetica Bold, 10'\n";
    std::vector<std::pair<double, double>> plotData;
    std::set<Triplet<double>> in;
    std::vector<double> b(300);
    for(int k = 1; k < 300; ++k) {
        for (int i = 0; i < 300; ++i){
            in.insert({ 1 + (i * (k-1)) / 300. , size_t(i), size_t(i)});
            b[i] = 300-i;
        }
        CSR<double> A = CSR<double>(300, 300, in);
        in.clear();
        int its = GradientDescent(A, b);
        std::cout<<GradientDescent(A,b)<<GaussSeidel(A,b)<<std::endl;
        plotData.emplace_back(k, its);
    }
    gp << "plot '-' u 1:2 with lines title 'Обусловленность -> число итераций' lc rgb 'red'\n";
    gp.send1d(plotData);*/


    /* std::set<Triplet<double>> in;
    for(size_t i = 0; i < 300; ++i){
        for(size_t j = 0; j < 300; ++j) {
           if(i == j) in.insert({300 + static_cast<double>(i) / 59.8, i, i});
           else in.insert({(-1 + static_cast <double > (rand()) /( static_cast <double > (RAND_MAX/(2))))/10, i, j});
        }
    }
    CSR<double> A = CSR<double>(300, 300, in);
    std::vector<double> b = GenerateVector<double>(300, -1, 1);
    double w = 1;
    int N = 10;
    std::cout << SOR(A, b, 1.)<<SSOR(A,b, 1.)<<GradientDescent(A, b);
    Gnuplot gp;
    gp << "set xlabel 'Номер итерации' \n"
          "set ylabel 'ln(|r|))'\n"
          "set grid\n"
          "set title 'Зависимость логарифма модуля невязки от номера итерации ПМИ' font 'Helvetica Bold, 10'\n";

    for(size_t i = 0; i < N; ++i) {
        w = 0.1 + i * 18. / 100;
        std::cout << SOR(A, b, w)<<SSOR(A,b, w);
    }
    std::vector<std::string> color = {"red", "green", "purple", "yellow", "blue", "orange", "gray", "pink", "brown", "black"};
    gp << "plot '../PlotData/FasterSimple/SOR0.100000.txt' with lines title 'w' lc rgb 'red',";
    for(size_t i = 1; i < N-1; ++i) {
        w = 0.1 + i * 18. / 100;
        gp << " '../PlotData/FasterSimple/SOR" + std::to_string(w) + ".txt' with lines title 'w' lc rgb '" + color[i] +"',";
    }
    gp << " '../PlotData/FasterSimple/SOR1.720000.txt' with lines title 'w' lc rgb 'black'\n";
*/
    return 0;
}

