#include <iostream>
#include <time.h>
#include <vector>

#include "hungarian.hpp"

void genmat(int n, int m, std::vector<float> &mat) {
  srand(time(0));
  mat.resize(n * m);
  for (int i = 0; i < mat.size(); i++)
    mat[i] = rand() % 100;
}

void dumpmat(int n, int m, std::vector<float> &mat) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++)
      printf("%f ", mat[i * m + j]);
    printf("\n");
  }
}

int main(int argc, char **argv) {
  std::cout << "Hungarian Algorithm" << std::endl;
  int n = 4;
  int m = 5;
  int MODE = 1; // 0 is minimize, 1 is maximize

  int nm_min = n < m ? n : m;

  Hungarian hu;

  std::vector<float> costMat;
  std::vector<float> assignment_idx;

  genmat(n, m, costMat);

  dumpmat(n, m, costMat);

  const float *costMatPtr = &costMat[0];

  hu.solve(&costMatPtr, n, m, MODE, assignment_idx);

  // dumpmat(n, m, costMat);

  // dumpmat(1, nm_min, assignment_idx);

  return 0;
}