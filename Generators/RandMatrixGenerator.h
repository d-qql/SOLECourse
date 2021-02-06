//
// Created by d-qql on 27.01.2021.
//

#ifndef SOLECOURSE_RANDMATRIXGENERATOR_H
#define SOLECOURSE_RANDMATRIXGENERATOR_H

#include <set>
#include "../utility/Triplet.h"

template<typename T>
std::set<Triplet<T>> GenerateMatrix(const size_t& size, const int &a, const int &b){
    srand(time(0));
    std::set<Triplet<T>> out;
    for(size_t i = 0; i < size; ++i){
        for(size_t j = 0; j < size; ++j){
            out.insert({static_cast<T>(a + rand()%(b-a+1)), i, j});
        }
    }
    return out;
}

template<typename T>
std::set<Triplet<T>> GenerateMatrix_FilledNumber(const size_t &size, const size_t& filledNumber, const int &a, const int&b){
    srand(time(0));
    std::set<Triplet<T>> out;
    while(out.size() < filledNumber){
        out.insert({static_cast<T>(a + rand()%(b-a+1)), static_cast<size_t>(rand() % size), static_cast<size_t>(rand() % size)});
    }
    return out;
}
#endif //SOLECOURSE_RANDMATRIXGENERATOR_H
