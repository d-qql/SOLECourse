#include <iostream>
//#include "utility/Triplet.h"
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
#include <ctime>
#include <fstream>
std::vector<std::pair<double, double>> testHeadGaussMethod_Time();
std::vector<std::pair<double, double>> testHeadGaussMethod_Norm();
std::vector<std::pair<double, double>> testGaussMethod_Time();
std::vector<std::pair<double, double>> testGaussMethod_Norm();
void CompareGaussNorm();
int main() {
    DenseMatrix<double> A = DenseMatrix<double>(7, 7, GenerateTridiagonal<double>(7, -200, 200));
    std::cout<<A;
    std::vector<double> b = GenerateVector<double>(7, -1, 1);

    std::cout<<b<<ThomasAlgorithm(A, b);

    return 0;
}

std::vector<std::pair<double, double>> testGaussMethod_Norm(){
    //ЛОГАРИФМ МОДУЛЯ НЕВЯЗКИ ОТ РАЗМЕРА МАТРИЦЫ headGaussMethod
    std::vector<std::pair<double, double>> plotData;
    for(int i = 3; i < 3000; ++i){
        DenseMatrix<double> A = DenseMatrix<double>(i, i, GenerateMatrix<double>(i, -3000, 3000));
        std::vector<double> b = GenerateVector<double>(i, -1, 1);
        std::cout<<A<<b;
        std::vector<double> res = GaussMethod(A, b);
        plotData.emplace_back(i, log(norm(b-A*res)));
        std::cout<<res;
    }
    Gnuplot gp;
    gp<<"set xlabel 'Matrix Dimension' \n"
        "set ylabel 'ln(|r|))'\n"
        "set grid\n"
        "set title 'Зависимость логарифма модуля невязки от размера матрицы GaussMethod' font 'Helvetica Bold, 10'\n";
    gp << "plot '-' with lines title 'dim->ln(|r|)' lc rgb 'blue'\n";
    gp.send1d(plotData);
    return plotData;
}

std::vector<std::pair<double, double>> testGaussMethod_Time(){
    //ИССЛЕДОВАНИЕ ВРЕМЕНИ ВЫПОЛНЕНИЯ ОТ РАЗМЕРА МАТРИЦЫ headGaussMethod
    std::vector<std::pair<double, double>> plotData;
    for(int i = 3; i < 500; ++i){
        DenseMatrix<double> A = DenseMatrix<double>(i, i, GenerateMatrix<double>(i, -300, 300));
        std::vector<double> b = GenerateVector<double>(i, -300, 300);
        clock_t end, start = clock();
        std::vector<double> res = GaussMethod(A, b);
        end = clock();
        plotData.emplace_back(log(i), log(double(end-start)/CLOCKS_PER_SEC));
        std::cout<<res;
    }
    Gnuplot gp;
    gp<<"set xlabel 'ln(Matrix Dimension)' \n"
        "set ylabel 'ln(Time)'\n"
        "set grid\n"
        "set title 'Зависимость времени выполнения от размера матрицы GaussMethod' font 'Helvetica Bold, 10'\n";
    gp << "plot '-' with lines title 'ln(dim)->ln(time)' lc rgb 'blue'\n";
    gp.send1d(plotData);
    return plotData;
}

std::vector<std::pair<double, double>> testHeadGaussMethod_Norm(){
    //ЛОГАРИФМ МОДУЛЯ НЕВЯЗКИ ОТ РАЗМЕРА МАТРИЦЫ headGaussMethod
    std::vector<std::pair<double, double>> plotData;
    for(int i = 3; i < 3000; ++i){
        DenseMatrix<double> A = DenseMatrix<double>(i, i, GenerateMatrix<double>(i, -3000, 3000));
        std::vector<double> b = GenerateVector<double>(i, -1, 1);
        std::vector<double> res = headGaussMethod(A, b);
        plotData.emplace_back(i, log(norm(b-A*res)));
        std::cout<<res;
    }
    Gnuplot gp;
    gp<<"set xlabel 'Matrix Dimension' \n"
        "set ylabel 'ln(|r|))'\n"
        "set grid\n"
        "set title 'Зависимость логарифма модуля невязки от размера матрицы headGaussMethod' font 'Helvetica Bold, 10'\n";
    gp << "plot '-' with lines title 'dim->ln(|r|)' lc rgb 'blue'\n";
    gp.send1d(plotData);
    return plotData;
}

std::vector<std::pair<double, double>> testHeadGaussMethod_Time(){
    //ИССЛЕДОВАНИЕ ВРЕМЕНИ ВЫПОЛНЕНИЯ ОТ РАЗМЕРА МАТРИЦЫ headGaussMethod
    std::vector<std::pair<double, double>> plotData;
    for(int i = 3; i < 500; ++i){
        DenseMatrix<double> A = DenseMatrix<double>(i, i, GenerateMatrix<double>(i, -300, 300));
        std::vector<double> b = GenerateVector<double>(i, -300, 300);
        clock_t end, start = clock();
        std::vector<double> res = headGaussMethod(A, b);
        end = clock();
        plotData.emplace_back(log(i), log(double(end-start)/CLOCKS_PER_SEC));
        std::cout<<res;
    }
    Gnuplot gp;
    gp<<"set xlabel 'ln(Matrix Dimension)' \n"
        "set ylabel 'ln(Time)'\n"
        "set grid\n"
        "set title 'Зависимость времени выполнения от размера матрицы headGaussMethod' font 'Helvetica Bold, 10'\n";
    gp << "plot '-' with lines title 'ln(dim)->ln(time)' lc rgb 'blue'\n";
    gp.send1d(plotData);
    return plotData;
}
void CompareGaussNorm(){
    std::vector<std::pair<double, double>> pureGauss = testGaussMethod_Norm();
    std::vector<std::pair<double, double>> headGauss = testHeadGaussMethod_Norm();
    std::ofstream fout("pureGaussNormPlotData.txt");
    for(auto i: pureGauss){
        fout<<i.first<<" "<<i.second<<"\n";
    }
    fout.close();
    fout = std::ofstream("headGaussNormPlotData.txt");
    for(auto i: headGauss){
        fout<<i.first<<" "<<i.second<<"\n";
    }
    fout.close();
    Gnuplot gp;
    gp<<"set xlabel 'Matrix Dimension' \n"
        "set ylabel 'ln(|r|))'\n"
        "set grid\n"
        "set title 'Зависимость логарифма модуля невязки от размера матрицы GaussMethod' font 'Helvetica Bold, 10'\n";
    gp << "plot ";
    gp << "     '-' with lines title 'dim->ln(|r|), pureGauss' ls 11 lc rgb 'blue',"
          "     '-' with lines title 'dim->ln(|r|), headGauss' ls 12 lc rgb 'red'\n";
    gp.send1d(pureGauss);
    gp.send1d(headGauss);
}
