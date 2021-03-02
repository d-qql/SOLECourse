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



    std::vector<T> r = b - A * x;
    r = 1. / norm(r) * r;
    Basis.push_back(r);
    T h;

    for(size_t j = 0; j < N - 1; ++j){
        r = A * Basis[j];
        for(size_t i = 0; i <= j; ++i){
            h = Basis[i] * r;
            H.insert({h, i, j});
            r = r - h * Basis[i];
        }
        h = norm(r);
        H.insert({h, j+1, j});
        r = 1. / h * r;
        Basis.push_back(r);
    }
    return std::tuple<std::vector<std::vector<T>>, std::set<Triplet<T>>>(Basis, H);
}

#endif //SOLECOURSE_KRYLOV_H
