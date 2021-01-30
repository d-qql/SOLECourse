//
// Created by d-qql on 31.01.2021.
//

#ifndef SOLECOURSE_TESTS_H
#define SOLECOURSE_TESTS_H

#include <vector>
#include <ctgmath>
#include "../Dense/DenseMatrix.h"
#include "../Generators/RandMatrixGenerator.h"
#include "../Generators/RandVectorGenerator.h"
#include "../Dense/LU.h"
#include "../utility/Consts.h"
#include "../Generators/RandTridiagonalMatrixGen.h"
#include "../Dense/ThomasAlgorithm.h"
#include "../Dense/GaussMethod.h"
#include "../Dense/HeadGaussMethod.h"
#include "../Sparse/CSR.h"
#include "../Sparse/SimpleIteration.h"

void testSimpleIteration_Speed(){
    std::set<Triplet<double>> in;
    std::vector<double> b(300);
    for(size_t i = 1; i <= 300; ++i){
        in.insert({double(i)/30, i-1, i-1});
        b[i-1] = double(301-i)/30;
    }
    CSR<double> A = CSR<double>(300, 300, in);
    //std::vector<double> tao = {0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.13, 0.15, 0.17, 0.18, 0.185, 0.19, 0.20, 0.23, 0.26, 0.31, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    std::cout<<SimpleIteration(A, b, 0.175);
    std::cout<<SimpleIteration(A, b, 0.18);
    std::cout<<SimpleIteration(A, b, 0.1818);
    std::cout<<SimpleIteration(A, b, 0.183);
    std::cout<<SimpleIteration(A, b, 0.188);

    Gnuplot gp;
    gp<<"set xlabel 'Номер итерации' \n"
        "set ylabel 'ln(|r|))'\n"
        "set grid\n"
        "set title 'Зависимость логарифма модуля невязки от номера итерации ПМИ' font 'Helvetica Bold, 10'\n";
    gp << "plot '../PlotData/SimpleIteration/SimpleIteration_Speed0.175000.txt' u 1:2 with lines title 'tao = 0.175' lc rgb 'blue', "
          "     '../PlotData/SimpleIteration/SimpleIteration_Speed0.180000.txt' u 1:2 with lines title 'tao = 0.18' lc rgb 'red', "
          "     '../PlotData/SimpleIteration/SimpleIteration_Speed0.181800.txt' u 1:2 with lines title 'tao = 0.1818' lc rgb 'green', "
          "     '../PlotData/SimpleIteration/SimpleIteration_Speed0.183000.txt' with lines title 'tao = 0.183' lc rgb 'black', "
          "     '../PlotData/SimpleIteration/SimpleIteration_Speed0.188000.txt' with lines title 'tao = 0.188' lc rgb 'purple'\n";

}

std::vector<std::pair<double, double>> testLUMethod_Norm(size_t dim);

std::vector<std::pair<double, double>> testLUMethod_Time(size_t dim);

std::vector<std::pair<double, double>> testGaussMethod_Norm(size_t dim);

std::vector<std::pair<double, double>> testGaussMethod_Time(size_t dim);

std::vector<std::pair<double, double>> testHeadGaussMethod_Norm(size_t dim);

std::vector<std::pair<double, double>> testHeadGaussMethod_Time(size_t dim);

void testThomasAlgo(size_t dim);

void compareStraightMethodsNorm(size_t dim);

void compareStraightMethodsTime(size_t dim);

void compareStraightMethodsTime(size_t dim){
    Gnuplot gp;
    gp<<"set xlabel 'Matrix Dimension' \n"
        "set ylabel 'Time'\n"
        "set grid\n"
        "set title 'Зависимость времени от размера матрицы для разных методов' font 'Helvetica Bold, 10'\n";
    gp << "plot ";
    gp << "     '-' with lines title 'dim->time, pureGauss' ls 11 lc rgb 'blue',"
          "     '-' with lines title 'dim->time, headGauss' ls 11 lc rgb 'green',"
         // "     '-' with lines title 'dim->time, LU' ls 12 lc rgb 'red'"
          "\n";
    gp.send1d(testGaussMethod_Time(dim));
    gp.send1d(testHeadGaussMethod_Time(dim));
    //gp.send1d(testLUMethod_Time(dim));
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
        plotData.emplace_back(log(i), log(double(end-start)/CLOCKS_PER_SEC));
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
        plotData.emplace_back(log(i),log(double(end-start)/CLOCKS_PER_SEC));
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



#endif //SOLECOURSE_TESTS_H
