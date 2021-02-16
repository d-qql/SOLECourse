//
// Created by d-qql on 14.02.2021.
//

#ifndef SOLECOURSE_CHEBYSHEV_H
#define SOLECOURSE_CHEBYSHEV_H
#include "vector"
#include <ctgmath>
template<typename T>
std::vector<T> ChebyshevRoots(std::pair<T, T> section, size_t PowOf2){
    size_t PolyOrder = std::pow(2, PowOf2);
    std::vector<T> roots(PolyOrder);
    for(size_t i = 1; i <= PolyOrder; ++i){
        roots[i-1] = (section.first + section.second) / 2 +
                (section.second - section.first) *
                static_cast<T>(cos(static_cast<double>(2 * i - 1) * M_PI_2 / PolyOrder)) / 2;
    }
    std::vector<size_t> idx;
    idx.push_back(0);
    idx.push_back(1);
    std::vector<size_t> next_idx;
    size_t curOrder;
    for(size_t i = 1; i < PowOf2; ++i){     //по всем порядкам полинома до текущего
        next_idx.resize(pow(2, double(i)+1));
        curOrder = pow(2, double (i) + 1);
        for(size_t j = 0; j < curOrder-1; j+=2){     //по всем элементам нового вектора прыгая через 1
            next_idx[j] = idx[j/2];
            next_idx[j+1] = curOrder - 1 - next_idx[j];
        }
        idx = next_idx;
    }
    std::vector<T> res(PolyOrder);
    for(size_t i = 0; i < PolyOrder; ++i){
        res[i] = roots[idx[i]];
    }
    std::cout<<roots<<res;
    return res;
}

#endif //SOLECOURSE_CHEBYSHEV_H
