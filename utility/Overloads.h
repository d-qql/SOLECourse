//
// Created by d-qql on 28.01.2021.
//

#ifndef SOLECOURSE_OVERLOADS_H
#define SOLECOURSE_OVERLOADS_H
#include <vector>

template<typename T>
std::vector<T> operator*(const T &k, std::vector<T> const &vec){
    std::vector<T> result;
    result.resize(vec.size());
    for(size_t i = 0; i<result.size(); ++i){
        result[i] = k*vec[i];
    }
    return result;
}

template<typename T>
std::vector<T> operator+(const std::vector<T>& fVec, const std::vector<T>& sVec){
    std::vector<T> res;
    if(fVec.size() == sVec.size()){
        res.resize(fVec.size());
        for(size_t i = 0; i < fVec.size(); ++i){
            res[i]=fVec[i]+sVec[i];
        }
    }
    return res;
}

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
T operator*(const std::vector<T>& a, const std::vector<T>& b){
    auto res = static_cast<T>(0);
    for(size_t i = 0; i < a.size(); ++i){
        res += a[i]*b[i];
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
