#ifndef _hungarian_hpp_
#define _hungarian_hpp_

#include <vector>

class Hungarian {
  std::vector<float> lx, ly;
  std::vector<bool> lx_mask, ly_mask;
  std::vector<float> matrix;
  std::vector<bool> mask_matrix;
  int matched_jobs;
  std::vector<int> alter_paths;
  int MODE;
  int N, M, max_NM;

public:
  void init(const float **cost, const int _N, const int _M);

  // MODE0 is maximum, MODE1 is minimum
  float solve(const float **cost, const int N, const int M, const int _MODE,
              std::vector<float> &assignment_index);

private:
  int step1();

  int step2();

  int step3();

  int step4();

  int step5(const float **cost, float *score,
            std::vector<float> *assignment_idx);
};

#endif //_hungarian_hpp_