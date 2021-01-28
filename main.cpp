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
#include "Dense/LU.h"
#include <ctime>
#include <fstream>
std::vector<std::pair<double, double>> testHeadGaussMethod_Time(size_t dim);
std::vector<std::pair<double, double>> testHeadGaussMethod_Norm(size_t dim);
std::vector<std::pair<double, double>> testGaussMethod_Time(size_t dim);
std::vector<std::pair<double, double>> testGaussMethod_Norm(size_t dim);
std::vector<std::pair<double, double>> testLUMethod_Time(size_t dim);
std::vector<std::pair<double, double>> testLUMethod_Norm(size_t dim);
void compareStraightMethodsTime(size_t dim);
void compareStraightMethodsNorm(size_t dim);
void testThomasAlgo(size_t dim);
int main() {
    compareStraightMethodsTime(1000);
    return 0;
}
void compareStraightMethodsTime(size_t dim){
    Gnuplot gp;
    gp<<"set xlabel 'Matrix Dimension' \n"
        "set ylabel 'Time'\n"
        "set grid\n"
        "set title 'Зависимость времени от размера матрицы для разных методов' font 'Helvetica Bold, 10'\n";
    gp << "plot ";
    gp << "     '-' with lines title 'dim->time, pureGauss' ls 11 lc rgb 'blue',"
          "     '-' with lines title 'dim->time, headGauss' ls 11 lc rgb 'green',"
          "     '-' with lines title 'dim->time, LU' ls 12 lc rgb 'red'\n";
    gp.send1d(testGaussMethod_Time(dim));
    gp.send1d(testHeadGaussMethod_Time(dim));
    gp.send1d(testLUMethod_Time(dim));
}
std::vector<std::pair<double, double>> testLUMethod_Norm(size_t dim){
    //ЛОГАРИФМ МОДУЛЯ НЕВЯЗКИ ОТ РАЗМЕРА МАТРИЦЫ headGaussMethod
    std::vector<std::pair<double, double>> plotData;
    for(int i = 3; i < dim; ++i){
        DenseMatrix<double> A = DenseMatrix<double>(i, i, GenerateMatrix<double>(i, -1000, 1000));
        std::vector<double> b = GenerateVector<double>(i, -1, 1);
        std::cout<<A<<b;
        std::vector<double> res = solveByLU(A, b);
        plotData.emplace_back(i, log(norm(b-A*res)));
        std::cout<<res;
    }
/*    Gnuplot gp;
    gp<<"set xlabel 'Matrix Dimension' \n"
        "set ylabel 'ln(|r|))'\n"
        "set grid\n"
        "set title 'Зависимость логарифма модуля невязки от размера матрицы GaussMethod' font 'Helvetica Bold, 10'\n";
    gp << "plot '-' with lines title 'dim->ln(|r|)' lc rgb 'blue'\n";
    gp.send1d(plotData);*/
    return plotData;
}
std::vector<std::pair<double, double>> testLUMethod_Time(size_t dim){
    std::vector<std::pair<double, double>> plotData;
    for(int i = 3; i <= dim; ++i){
        DenseMatrix<double> A = DenseMatrix<double>(i, i, GenerateMatrix<double>(i, -1000, 1000));
        std::vector<double> b = GenerateVector<double>(i, -1, 1);
        clock_t end, start = clock();
        std::vector<double> res = solveByLU(A, b);
        end = clock();
        plotData.emplace_back(i, double(end-start)/CLOCKS_PER_SEC);
        std::cout<<res;
    }
    return plotData;
}

void testThomasAlgo(size_t dim){
    DenseMatrix<double> A = DenseMatrix<double>(dim, dim, GenerateTridiagonal<double>(dim, -200, 200));
    std::cout<<A;
    std::vector<double> b = GenerateVector<double>(dim, -1, 1);
    std::cout<<b<<ThomasAlgorithm(A, b);
}
std::vector<std::pair<double, double>> testGaussMethod_Norm(size_t dim){
    //ЛОГАРИФМ МОДУЛЯ НЕВЯЗКИ ОТ РАЗМЕРА МАТРИЦЫ headGaussMethod
    std::vector<std::pair<double, double>> plotData;
    for(int i = 3; i < dim; ++i){
        DenseMatrix<double> A = DenseMatrix<double>(i, i, GenerateMatrix<double>(i, -1000, 1000));
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

std::vector<std::pair<double, double>> testGaussMethod_Time(size_t dim){
    //ИССЛЕДОВАНИЕ ВРЕМЕНИ ВЫПОЛНЕНИЯ ОТ РАЗМЕРА МАТРИЦЫ headGaussMethod
    std::vector<std::pair<double, double>> plotData;
    for(int i = 3; i < dim; ++i){
        DenseMatrix<double> A = DenseMatrix<double>(i, i, GenerateMatrix<double>(i, -1000, 1000));
        std::vector<double> b = GenerateVector<double>(i, -1, 1);
        clock_t end, start = clock();
        std::vector<double> res = GaussMethod(A, b);
        end = clock();
        plotData.emplace_back(i,double(end-start)/CLOCKS_PER_SEC);
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

std::vector<std::pair<double, double>> testHeadGaussMethod_Norm(size_t dim){
    //ЛОГАРИФМ МОДУЛЯ НЕВЯЗКИ ОТ РАЗМЕРА МАТРИЦЫ headGaussMethod
    std::vector<std::pair<double, double>> plotData;
    for(int i = 3; i < dim; ++i){
        DenseMatrix<double> A = DenseMatrix<double>(i, i, GenerateMatrix<double>(i, -1000, 1000));
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

std::vector<std::pair<double, double>> testHeadGaussMethod_Time(size_t dim){
    //ИССЛЕДОВАНИЕ ВРЕМЕНИ ВЫПОЛНЕНИЯ ОТ РАЗМЕРА МАТРИЦЫ headGaussMethod
    std::vector<std::pair<double, double>> plotData;
    for(int i = 3; i < dim; ++i){
        DenseMatrix<double> A = DenseMatrix<double>(i, i, GenerateMatrix<double>(i, -1000, 1000));
        std::vector<double> b = GenerateVector<double>(i, -1, 1);
        clock_t end, start = clock();
        std::vector<double> res = headGaussMethod(A, b);
        end = clock();
        plotData.emplace_back(i, double(end-start)/CLOCKS_PER_SEC);
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
void compareStraightMethodsNorm(size_t dim){
    std::vector<std::pair<double, double>> pureGauss = testGaussMethod_Norm(dim);
    std::vector<std::pair<double, double>> headGauss = testHeadGaussMethod_Norm(dim);
    std::vector<std::pair<double, double>> LU = testLUMethod_Norm(dim);
    std::ofstream fout("pureGaussNormPlotData.txt");
    /*for(auto i: pureGauss){
        fout<<i.first<<" "<<i.second<<"\n";
    }
    fout.close();
    fout = std::ofstream("headGaussNormPlotData.txt");
    for(auto i: headGauss){
        fout<<i.first<<" "<<i.second<<"\n";
    }
    fout.close();*/
    Gnuplot gp;
    gp<<"set xlabel 'Matrix Dimension' \n"
        "set ylabel 'ln(|r|))'\n"
        "set grid\n"
        "set title 'Зависимость логарифма модуля невязки от размера матрицы разных методов' font 'Helvetica Bold, 10'\n";
    gp << "plot ";
    gp << "     '-' with lines title 'dim->ln(|r|), pureGauss' ls 11 lc rgb 'blue',"
          "     '-' with lines title 'dim->ln(|r|), headGauss' ls 11 lc rgb 'green',"
          "     '-' with lines title 'dim->ln(|r|), LU' ls 12 lc rgb 'red'\n";
    gp.send1d(pureGauss);
    gp.send1d(headGauss);
    gp.send1d(LU);
}
