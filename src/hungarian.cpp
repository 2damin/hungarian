#include "hungarian.hpp"

#include <iostream>
#include <algorithm>
#include <vector>

void Hungarian::init(const float **cost, const int N, const int M)
{
    lx.resize(N);
    ly.resize(M);

    int max_NM = N > M ? N : M;
    float max_value = 0;
    for(int i = 0; i < N; ++i)
        for(int j =0; j< M; ++j)
            max_value = std::max(max_value, (*cost)[i * M + j]);

    //resize square matrix
    if(N != M)
    {
        matrix.resize(max_NM*max_NM);
        for(int i = 0; i < N;++i)
            for(int j = 0; j < M;++j)
                matrix[j + i * max_NM] = (*cost)[j + i * M];

        //fill max value in new element of square matrix
        for(int i = N; i < max_NM; ++i)
            for(int j = 0; j < M; ++j)
                matrix[i * max_NM + j] = max_value;

        for(int j = M; j < max_NM; ++j)
            for(int i = 0; i < N; ++i)
                matrix[i * max_NM + j] = max_value;
    }

    //reset mask
    lx_mask.resize(max_NM,false);
    ly_mask.resize(max_NM,false);

    //infinity_check
    auto inf = std::numeric_limits<double>::infinity();
    auto max_elem = matrix[0];
    for(auto m : matrix)
    {
        if(max_elem == inf)
            max_elem = m;
        else
            max_elem = std::max(max_elem, m);
    }

    if(max_elem == inf)
        max_elem = 0;
    else
        max_elem++;
    
    for(auto &m : matrix)
        if(max_elem == inf)
            m = max_elem;

    std::cout << "--PRINT RESIZED MATRIX--" << std::endl;
    for(int i = 0; i < max_NM;++i)
    {
        for(int j = 0; j < max_NM;++j)
            std::cout << matrix[j + i * max_NM] << " ";
        std::cout << std::endl;
    }

    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j < M; ++j)
        {
            lx[i] = std::max(lx[i], (*cost)[i * M + j]);
        }
        std::cout << "lx " << lx[i] << std::endl;
    }

}

bool Hungarian::augment(const float **cost, const int N, const int M)
{
    std::vector<int> A;
    std::vector<int> B;

    for(int i = 0; i < N; ++i)
    {
        for(int j = 0 ; j < M; ++j)
        {
            int score = lx[i] + ly[j] - (*cost)[i * M + j];
            if(score == 0)
            {
                A.push_back((int)i);
                B.push_back((int)j);
                break;
            }
        }
    }
    for(auto a : A)
        std::cout << "A : " << a << std::endl;
    for(auto b : B)
        std::cout << "B : " << b << std::endl;

    return false;
}

float Hungarian::solve(const float **cost, const int N, const int M, const int _MODE, float *assignment_index){
    std::cout << "solve" << std::endl;
    MODE = _MODE;

    std::cout << __LINE__ << std::endl;
    init(cost, N, M);

    std::cout << __LINE__ << std::endl;

    augment(cost, N, M);

    std::cout << __LINE__ << std::endl;

    return 0;
}