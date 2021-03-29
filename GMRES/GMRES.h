//
// Created by d-qql on 28.03.2021.
//

#ifndef SOLECOURSE_GMRES_H
#define SOLECOURSE_GMRES_H
#include "../Sparse/CSR.h"
#include "../Dense/BackSubstitution.h"
#include "../Givens/Rotations.h"
#include "../Krylov/Krylov.h"
template<typename T>
std::vector<T> GMRES(const CSR<T>& A, const std::vector<T>& b ){
    std::vector<T> x(b.size());
    std::vector<T> r = b - A * x;
    int m = A.sizeW();
    T N = norm(r);
    std::vector<T> B;
    std::vector<T> temp_x;
    int K = 0;
    while ( Tabs(N) > tolerance<T> ){
        //std::cout<<x;
        for(int i = 2; i <= m+1; ++i){
            B.resize(i, 0);
            B[0] = N;
            auto [V, H] = KrylovSubSpace(A, r, i);
           // std::cout<<V<<std::endl;
           // std::cout<<H<<std::endl;
            HessenbergRot(H, B);
           // std::cout<<H<<std::endl;
            std::cout<<"итерация: "<<i<<std::endl;
           // std::cout<<H<<std::endl;
            std::cout<<"НЕВЯЗКА: "<<Tabs(B.back())<<"\n";
            if(Tabs(B.back()) < tolerance<T> || i == m+1){
                H.__deleteLastRow();
                B.pop_back();
                // std::cout<<B;
                temp_x = backSubstTopTriangular(H, B);
                //std::vector<T> dx = V * temp_x;
                //std::cout<<"Обратный ход: "<<B-H*temp_x<<std::endl;
                // std::cout<<"Dx: "<<V*temp_x;
                x = x + V * temp_x;
                r = b - A * x;
                N = norm(r);
                B.clear();
                ++K;
                break;
            }
            B.clear();
        }
    }
    std::cout<<"НЕВЯЗКА: "<<N<<"\n";
    return x;
}
#endif //SOLECOURSE_GMRES_H
