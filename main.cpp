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

int main() {

    testMultiply_Time_DifferentDim(1000);
    /*std::set<Triplet<int>> in = {{1, 0, 0}, {2, 0, 1}, {3, 0, 3}, {4, 1, 2}};
    CSR<int> A = CSR<int>(3, 4, in);
    A.print();
    std::set<Triplet<int>> in1 = {{1, 0, 0}, {2, 0, 1}, {4, 1, 1}, {2, 2, 1}, {6, 2, 2}};
    CSR<int> K = CSR<int>(3, 3, in1);
    K.print();
    */
/*    Gnuplot gp;
    gp<<"set xlabel 'Число ненулевых элементов' \n"
        "set ylabel 'Время умножения'\n"
        "set xrange [0:150000]\n"
        "set yrange [0:0.003]\n"
        "set grid\n"
        "set title 'Зависимость времени умножения матрицы на плотный вектор от заполненности' font 'Helvetica Bold, 10'\n";
    gp << "plot '../PlotData/MultiplyTime/Dense.txt' u 1:2 with lines title 'Dense' lc rgb 'blue',"
          "     '../PlotData/MultiplyTime/CSR.txt' u 1:2 with lines title 'CSR' lc rgb 'red'\n";*/
/*    Gnuplot gp;
    gp<<"set xlabel 'Тао' \n"
        "set ylabel 'Число итераций'\n"
        "set grid\n"
        "set title 'Зависимость числа итераций от Тао для схождения с точностью 1e-12' font 'Helvetica Bold, 10'\n";
    gp << "plot '../PlotData/SimpleIteration/EpsilonNumbTao.txt' with lines title 'Taо -> IterNumb' lc rgb 'blue'\n";*/
    return 0;
}

