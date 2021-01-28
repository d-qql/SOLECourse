//
// Created by d-qql on 28.01.2021.
//

#ifndef SOLECOURSE_OVERLOADS_H
#define SOLECOURSE_OVERLOADS_H
#include <vector>

template<typename T>
std::vector<T> operator-(const std::vector<T>& fVec, const std::vector<T>& sVec){
    std::vector<T> res;
    if(fVec.size() == sVec.size()){
        res.resize(fVec.size());
        for(size_t i = 0; i < fVec.size(); ++i){
            res[i]=fVec[i]-sVec[i];
        }
    }
    return res;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& b){
    os<<"( ";
    for(auto i: b){
        os<<i<<" ";
    }
    os<<")\n";
    return os;
}

#endif //SOLECOURSE_OVERLOADS_H
