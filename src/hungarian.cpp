#include "hungarian.hpp"

#include <algorithm>
#include <iostream>
#include <queue>
#include <vector>

void Hungarian::init(const float **cost, const int _N, const int _M) {
  lx.resize(N, 100);
  ly.resize(M, 0);

  N = _N;
  M = _M;

  max_NM = N > M ? N : M;

  // reset mask
  lx_mask.resize(max_NM, false);
  ly_mask.resize(max_NM, false);

  float max_value = 100;
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < M; ++j)
      max_value = std::max(max_value, (*cost)[i * M + j]);

  // resize square matrix
  if (N != M) {
    matrix.resize(max_NM * max_NM);
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
        matrix[j + i * max_NM] = (*cost)[j + i * M];

    // fill max value in new element of square matrix
    for (int i = N; i < max_NM; ++i)
      for (int j = 0; j < M; ++j)
        matrix[i * max_NM + j] = max_value;
    for (int j = M; j < max_NM; ++j)
      for (int i = 0; i < N; ++i)
        matrix[i * max_NM + j] = max_value;
  } else {
    matrix.resize(N * M);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < M; ++j)
        matrix[j + i * M] = (*cost)[j + i * M];
    }
  }

  // check infinity value
  auto inf = std::numeric_limits<double>::infinity();
  auto max_elem = matrix[0];
  for (auto m : matrix) {
    if (max_elem == inf)
      max_elem = m;
    else
      max_elem = std::max(max_elem, m);
  }

  if (max_elem == inf)
    max_elem = 0;
  else
    max_elem++;

  for (auto &m : matrix)
    if (max_elem == inf)
      m = max_elem;

  std::cout << "--PRINT RESIZED MATRIX--" << std::endl;
  for (int i = 0; i < max_NM; ++i) {
    for (int j = 0; j < max_NM; ++j)
      std::cout << matrix[j + i * max_NM] << " ";
    std::cout << std::endl;
  }

  for (int i = 0; i < max_NM; ++i) {
    for (int j = 0; j < max_NM; ++j) {
      lx[i] = std::min(lx[i], matrix[j + i * max_NM]);
    }
  }
}

// columnwise
int Hungarian::step1() {

  for (int i = 0; i < max_NM; ++i) {
    auto minCol = *std::min_element(matrix.begin() + max_NM * i,
                                    matrix.begin() + max_NM * (i + 1));
    for (int j = 0; j < max_NM; ++j)
      matrix[j + i * max_NM] -= minCol;
  }

  std::cout << "--PRINT STEP1--" << std::endl;
  for (int i = 0; i < max_NM; ++i) {
    for (int j = 0; j < max_NM; ++j)
      std::cout << matrix[j + i * max_NM] << " ";
    std::cout << std::endl;
  }

  return 2;
}

// rowwise
int Hungarian::step2() {
  for (int j = 0; j < max_NM; ++j) {
    float minVal = 100;
    for (int i = 0; i < max_NM; ++i)
      minVal = std::min(minVal, matrix[i * max_NM + j]);
    std::cout << minVal << std::endl;
    for (int i = 0; i < max_NM; ++i) {
      matrix[i * max_NM + j] -= minVal;
      if (matrix[i * max_NM + j] == 0)
        mask_matrix[i * max_NM + j] = true;
    }
  }

  std::cout << "--PRINT STEP2--" << std::endl;
  for (int i = 0; i < max_NM; ++i) {
    for (int j = 0; j < max_NM; ++j)
      std::cout << matrix[j + i * max_NM] << " ";
    std::cout << std::endl;
  }

  return 3;
}

struct cmp {
  bool operator()(std::pair<int, int> a, std::pair<int, int> b) {
    return a.second < b.second;
  }
};

bool dfs_zero(int *count, std::vector<bool> mask_matrix,
              std::vector<bool> lx_mask, std::vector<bool> ly_mask, const int N,
              const int M) {
  int max_count_idx = -1;
  int max_count = 0;
  std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>,
                      cmp>
      cnts;
  for (int i = 0; i < N; ++i) {
    auto cnt = std::count(mask_matrix.begin() + M * i,
                          mask_matrix.begin() + M * (i + 1), 1);
    cnts.push(std::make_pair(i, cnt));
  }
  for (int j = 0; j < M; ++j) {
    int cnt = 0;
    for (int i = 0; i < N; ++i)
      cnt += mask_matrix[j + i * M] == 1 ? 1 : 0;
    cnts.push(std::make_pair(j + M, cnt));
  }
  //   for (int i = 0; i < N; ++i) {
  //     for (int j = 0; j < M; ++j)
  //       std::cout << mask_matrix[j + i * N] << " ";
  //     std::cout << std::endl;
  //   }

  // erase zero line.
  // ----------------
  // cnt.first is index of x,y line, cnt.second is the number of zero value in
  // cnt.first's line.

  int before_cnt = 0;
  std::vector<bool> mask_origin;
  mask_origin.assign(mask_matrix.begin(), mask_matrix.end());
  while (!cnts.empty()) {
    auto cnt = cnts.top();
    cnts.pop();
    // std::cout << cnt.first << " " << cnt.second << std::endl;
    if (cnt.second == 0)
      return true;

    // columnwise (cnt.first / N == 1), rowwise (cnt.first / N == 0)
    if (cnt.first / N == 1) {
      int j_idx = cnt.first % N;
      ly_mask[j_idx] = true;
      for (int i = 0; i < M; ++i)
        mask_matrix[j_idx + i * N] = 0;
    } else {
      int i_idx = cnt.first % N;
      lx_mask[i_idx] = true;
      for (int j = 0; j < N; ++j)
        mask_matrix[j + i_idx * N] = 0;
    }

    int count_here = *count;
    ++*count;
    // std::cout << "count : " << *count << std::endl;
    auto status = dfs_zero(count, mask_matrix, lx_mask, ly_mask, N, M);
    if (status) {
      if (*count < N) {
        std::cout << "Minimum lines : " << *count << std::endl;
        return false;
      } else if (*count == N) {
        // std::cout << "Final count : " << *count << std::endl;
      }
    }

    auto next_cnt = cnts.top();
    if (next_cnt.second < cnt.second)
      break;

    // std::cout << "RECOVER" << std::endl;
    mask_matrix.assign(mask_origin.begin(), mask_origin.end());
    *count = count_here;

    // for (int i = 0; i < N; ++i) {
    //   for (int j = 0; j < M; ++j)
    //     std::cout << mask_matrix[j + i * N] << " ";
    //   std::cout << std::endl;
    // }
  }
  return true;
}

// find zero line
int Hungarian::step3() {
  std::cout << "step3" << std::endl;
  int count = 0;
  auto status = dfs_zero(&count, mask_matrix, lx_mask, ly_mask, max_NM, max_NM);
  std::cout << " Count : " << count << std::endl;

  if (!status)
    return 4;
  else
    return 5;
}

bool min_comp(float a, float b) {
  if (a == 0)
    return false;
  return a < b;
}

int Hungarian::step4() {
  std::cout << "step4" << std::endl;

  std::cout << "--PRINT STEP4--" << std::endl;
  for (int i = 0; i < max_NM; ++i) {
    for (int j = 0; j < max_NM; ++j)
      std::cout << matrix[j + i * max_NM] << " ";
    std::cout << std::endl;
  }

  float delta = 100;
  for (int i = 0; i < max_NM; ++i) {
    for (int j = 0; j < max_NM; ++j) {
      auto tmp = std::min(delta, matrix[j + i * max_NM]);
      delta = tmp == 0 ? delta : tmp;
    }
  }

  //   for (int i = 0; i < max_NM; ++i)
  //     std::cout << "xline : " << lx_mask[i] << " ";
  //   std::cout << std::endl;
  //   for (int j = 0; j < max_NM; ++j)
  //     std::cout << "yline : " << ly_mask[j] << " ";
  //   std::cout << std::endl;
  //   std::cout << "--------" << std::endl;
  std::cout << "delta : " << delta << std::endl;

  //   std::transform(matrix.begin(), matrix.end(), matrix.begin(),
  //                  [delta](float c) { return c == 0 ? c : c - delta; });
  std::transform(matrix.begin(), matrix.end(), matrix.begin(),
                 [delta](float c) { return c - delta; });

  for (int i = 0; i < max_NM; ++i)
    for (int j = 0; j < max_NM; ++j)
      matrix[j + i * max_NM] = std::max(
          0.f, matrix[j + i * max_NM] + delta * (lx_mask[i] + ly_mask[j]));

  std::cout << "--PRINT STEP4 AFTER--" << std::endl;
  for (int i = 0; i < max_NM; ++i) {
    for (int j = 0; j < max_NM; ++j)
      std::cout << matrix[j + i * max_NM] << " ";
    std::cout << std::endl;
  }

  return 2;
}

bool DFS_step5_columnwise(float *score, std::vector<float> &mask_matrix,
                          const float **cost, std::vector<bool> m_mask,
                          const int N, const int M, const int N_idx) {
  if (N_idx >= N) {
    std::cout << "Succed : score_" << *score << std::endl;
    return true;
  }

  for (int j = 0; j < M; ++j) {
    if (mask_matrix[j + M * N_idx] == 0 && !m_mask[j]) {
      m_mask[j] = true;
      auto subscore = *score + (*cost)[j + M * N_idx];
      //*score += (*cost)[j + M * N_idx];
      // std::cout << "s : " << *score << std::endl;
      auto status = DFS_step5_columnwise(&subscore, mask_matrix, cost, m_mask,
                                         N, M, N_idx + 1);
      if (!status) {
        m_mask[j] = false;
        continue;
      } else {
        *score = subscore;
        std::cout << "Succed : score_" << *score << std::endl;
        return true;
      }
    }
  }
  return false;
}

bool DFS_step5_rowwise(float *score, std::vector<float> &mask_matrix,
                       const float **cost, std::vector<bool> n_mask,
                       const int N, const int M, const int M_idx) {
  if (M_idx >= M) {
    std::cout << "Succed : score_" << *score << std::endl;
    return true;
  }

  for (int i = 0; i < N; ++i) {
    if (mask_matrix[M_idx + M * i] == 0 && !n_mask[i]) {
      n_mask[i] = true;
      auto subscore = *score + (*cost)[M_idx + M * i];

      std::cout << "s : " << *score << std::endl;
      auto status = DFS_step5_rowwise(&subscore, mask_matrix, cost, n_mask, N,
                                      M, M_idx + 1);
      if (!status) {
        n_mask[i] = false;
        continue;
      } else {
        *score = subscore;
        std::cout << "Succed : score_" << *score << std::endl;
        return true;
      }
    }
  }
  return false;
}

int Hungarian::step5(const float **cost, float *score, float **assignment_idx) {
  std::cout << "step5" << std::endl;
  std::cout << "--PRINT STEP5 --" << std::endl;
  for (int i = 0; i < max_NM; ++i) {
    for (int j = 0; j < max_NM; ++j)
      std::cout << matrix[j + i * max_NM] << " ";
    std::cout << std::endl;
  }

  std::transform(lx_mask.begin(), lx_mask.end(), lx_mask.begin(),
                 [](bool num) { return false; });
  std::transform(ly_mask.begin(), ly_mask.end(), ly_mask.begin(),
                 [](bool num) { return false; });

  // auto m_mask = std::vector<bool>(M, false);
  std::cout << "=====COLUMNWISE====" << std::endl;
  auto status = DFS_step5_columnwise(score, matrix, cost, lx_mask, N, M, 0);
  if (!status) {
    std::cout << "=====ROWISE====" << std::endl;
    DFS_step5_rowwise(score, matrix, cost, ly_mask, N, M, 0);
  }

  std::cout << "score : " << *score << std::endl;

  return 0;
}

float Hungarian::solve(const float **cost, const int N, const int M,
                       const int _MODE, float *assignment_index) {
  std::cout << "solve" << std::endl;
  MODE = _MODE;

  init(cost, N, M);

  mask_matrix.resize(N * M, false);

  float score = 0.f;
  int step = 1;
  do {
    switch (step) {
    case 1:
      step = step1();
      break;
    case 2:
      step = step2();
      break;
    case 3:
      step = step3();
      break;
    case 4:
      step = step4();
      break;
    case 5:
      step = step5(cost, &score, &assignment_index);
      break;
    }
  } while (step != 0);

  return 0;
}