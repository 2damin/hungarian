#include <iostream>
#include <vector>
#include <time.h>

#include "hungarian.hpp"

void genmat(int n, int m, std::vector<float>& mat)
{
    srand(time(0));
    mat.resize(n * m);
    for (int i=0; i < mat.size(); i++) mat[i] = rand() % 100;
}

void dumpmat(int n, int m, std::vector<float>& mat)
{
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<m; j++)
            printf("%f ", mat[i * m + j]);
        printf("\n");
    }
}

int main(int argc, char** argv)
{
    std::cout << "Hungarian Algorithm" << std::endl;
    int n = 10;
    int m = 10;

    Hungarian hu;

    std::vector<float> costMat;
    std::vector<float> assignment_idx;

    genmat(n,m,costMat);

    dumpmat(n,m,costMat);

    const float* costMatPtr = &costMat[0];

    hu.solve(&costMatPtr, n, m, 0, &assignment_idx[0]);

    return 0;
}