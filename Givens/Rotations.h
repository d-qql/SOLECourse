//
// Created by d-qql on 28.03.2021.
//

#ifndef SOLECOURSE_ROTATIONS_H
#define SOLECOURSE_ROTATIONS_H

#include "../Dense/DenseMatrix.h"
#include "../utility/Consts.h"
#include "../utility/Overloads.h"
template<typename T>
void HessenbergRot(DenseMatrix<T>& H, std::vector<T>& b){
    T cos;
    T sin;
    T up;
    T down;
    T Bup;
    T Bdown;
    for(size_t k = 0; k < H.sizeW(); ++k){
        cos = H(k, k)/sqrt(pow(H(k, k), 2) + pow(H(k + 1, k), 2));
        sin = H(k + 1, k)/sqrt(pow(H(k, k), 2) + pow(H(k + 1, k), 2));
        Bup = cos * b[k] + sin * b[k+1];
        Bdown = -sin * b[k] + cos * b[k+1];
        b[k] = Bup;
        b[k+1] = Bdown;
        for(size_t j = k; j < H.sizeW(); ++j){
            up = cos * H(k, j) + sin * H(k + 1, j);
            down = -sin * H(k, j) + cos * H(k + 1, j);
            if(Tabs(down) < tolerance<T>) down = 0;
            H(k, j) = up;
            H(k + 1, j) = down;
        }
    }
   // std::cout<<H;
}

#endif //SOLECOURSE_ROTATIONS_H
