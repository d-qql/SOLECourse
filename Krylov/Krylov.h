//
// Created by d-qql on 02.03.2021.
//

#ifndef SOLECOURSE_KRYLOV_H
#define SOLECOURSE_KRYLOV_H
#include "../utility/Consts.h"
#include "../utility/Overloads.h"
#include "../Sparse/CSR.h"
#include <tuple>

template<typename T>
std::tuple<std::vector<std::vector<T>>, std::set<Triplet<T>>> KrylovSubSpace(const CSR<T>& A, const std::vector<T>& b, std::vector<T> x, size_t N){
    std::set<Triplet<T>> H;
    std::vector<std::vector<T>> Basis;



    std::vector<T> v = b - A * x;
    v = 1. / norm(v) * v;
    Basis.push_back(v);
    T h;

    for(size_t j = 0; j < N; ++j){
        v = A * Basis[j];
        for(size_t i = 0; i <= j; ++i){
            h = Basis[i] * v;
            if( h > tolerance<T>) H.insert({h, i, j});
            v = v - h * Basis[i];
        }
        h = norm(v);
        if( h > tolerance<T>) H.insert({h, j+1, j});
        v = 1. / h * v;
        Basis.push_back(v);
    }
    return std::tuple<std::vector<std::vector<T>>, std::set<Triplet<T>>>(Basis, H);
}

#endif //SOLECOURSE_KRYLOV_H
