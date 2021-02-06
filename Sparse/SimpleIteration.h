//
// Created by d-qql on 30.01.2021.
//

#ifndef SOLECOURSE_SIMPLEITERATION_H
#define SOLECOURSE_SIMPLEITERATION_H
#include "CSR.h"
#include "../utility/Consts.h"
#include "../utility/gnuplot-iostream.h"
template<typename T>
std::vector<T> SimpleIteration(const CSR<T> &A, const std::vector<T> &b, const T &tao){

    std::vector<T> x(A.sizeH());
    std::vector<T> r = A * x - b;
    int i = 0;
    T N = norm(r);
    std::vector<std::pair<int, double>> plotData;
    while (N > tolerance<T>) {
        x = -tao * (A * x - b) + x;
        r = A * x - b;
        plotData.emplace_back(i, log(N));
        ++i;
        N = norm(r);
        std::cout<<N<<std::endl;
    }


    std::ofstream fout;
    fout.open("../PlotData/SimpleIteration/IterNumbTao.txt", std::ios::app);
    fout<<tao<<" "<<i<<"\n";
    fout.close();

/*
    std::ofstream fout("../PlotData/SimpleIteration/SimpleIteration_Speed" + std::to_string(tao) + ".txt");
    for(auto k: plotData){
        fout<<k.first<<" "<<k.second<<"\n";
    }
    fout.close();

    Gnuplot gp;
    gp<<"set xlabel 'Номер итерации' \n"
        "set ylabel 'ln(|r|))'\n"
        "set grid\n"
        "set title 'Зависимость логарифма модуля невязки от номера итерации ПМИ' font 'Helvetica Bold, 10'\n";
    gp << "plot '-' with lines title 'Iter->ln(|r|)' lc rgb 'blue'\n";
    gp.send1d(plotData);*/
    return x;
}

#endif //SOLECOURSE_SIMPLEITERATION_H
