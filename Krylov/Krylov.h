//
// Created by d-qql on 02.03.2021.
//

#ifndef SOLECOURSE_KRYLOV_H
#define SOLECOURSE_KRYLOV_H
#include "../utility/Consts.h"
#include "../utility/Overloads.h"
#include "../Sparse/CSR.h"
#include <tuple>
#include "../Dense/DenseMatrix.h"

template<typename T>
std::tuple<DenseMatrix<T>, DenseMatrix<T>> KrylovSubSpace(const CSR<T>& A, const std::vector<T>& r, size_t N){
    std::set<Triplet<T>> H;
    std::vector<std::vector<T>> Basis;
    DenseMatrix<T> V = DenseMatrix<T>(r.size(), N-1);

    std::vector<T> v = r;
    v = 1. / norm(v) * v;
    Basis.push_back(v);
    T h;

    for(size_t j = 0; j < N-1; ++j){
        v = A * Basis[j];
        for(size_t i = 0; i <= j; ++i){
            h = Basis[i] * v;
            if( Tabs(h) > tolerance<T>) H.insert({h, i, j});
            v = v - h * Basis[i];
        }
        h = norm(v);
        if( Tabs(h) > tolerance<T>) H.insert({h, j+1, j});
        v = 1. / h * v;
        Basis.push_back(v);
    }
    for(size_t j = 0; j < N-1; ++j){
        for(size_t i = 0; i < r.size(); ++i){
            V(i, j) = Basis[j][i];
        }
    }

    return std::tuple<DenseMatrix<T>, DenseMatrix<T>>(V, DenseMatrix<T>(N, N-1, H));
}

template<typename T>
std::tuple<std::vector<std::vector<T>>, std::set<Triplet<T>>> GoodKrylovSubSpace(const CSR<T>& A, const std::vector<T>& b, std::vector<T> v, size_t N){
    std::set<Triplet<T>> H;
    std::vector<std::vector<T>> Basis;



    v = b - A * v;
    v = 1. / norm(v) * v;
    Basis.push_back(v);
    T h;
    T tao;
    T p;
    std::vector<T> temp_h;
    for(size_t j = 0; j < N-1; ++j){
        v = A * Basis[j];           // v is t here
        tao = norm(v);
        temp_h.resize(j+1);
        for(size_t i = 0; i <= j; ++i){
            h = Basis[i] * v;
            v = v - h * Basis[i];
            temp_h[i] = h;
        }
        if(abs(norm(v) / tao) <= tolerance<T>){
            for(size_t i = 0; i <=j; ++i){
                p = Basis[i]*v;
                v = v - p*Basis[i];
                temp_h[i] = temp_h[i] + p;
            }
        }
        for(size_t i = 0; i <= j; ++i) {
            if (Tabs(temp_h[i]) > tolerance<T>) H.insert({h, i, j});
        }
        temp_h.clear();
        h = norm(v);
        if( Tabs(h) > tolerance<T>) H.insert({h, j+1, j});
        v = 1. / h * v;
        Basis.push_back(v);
    }
    return std::tuple<std::vector<std::vector<T>>, std::set<Triplet<T>>>(Basis, H);
}


#endif //SOLECOURSE_KRYLOV_H
