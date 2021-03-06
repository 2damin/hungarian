#include "hungarian.hpp"

#include <algorithm>
#include <iostream>
#include <queue>
#include <vector>

void Hungarian::init(const float **cost, const int _N, const int _M,
                     const int _MODE) {
  N = _N;
  M = _M;
  MODE = _MODE;
  max_NM = N > M ? N : M;

  int mode_type = MODE == 0 ? 1 : -1;

  // reset masks
  lx_mask.resize(max_NM, false);
  ly_mask.resize(max_NM, false);
  mask_matrix.resize(max_NM * max_NM, false);

  float max_value = 0;
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < M; ++j)
      max_value = std::max(max_value, (*cost)[i * M + j]);

  // if MODE is minimizing, max_value is max value in matrix.
  // MODE is maximizing, max_value is min value in matrix.
  max_value *= mode_type;

  // resize square matrix
  if (N != M) {
    matrix.resize(max_NM * max_NM);
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
        matrix[j + i * max_NM] = (*cost)[j + i * M] * mode_type;

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
        matrix[j + i * M] = (*cost)[j + i * M] * mode_type;
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
}

// (columnwise) subtract minimal value in each nth column
int Hungarian::step1() {
  for (int i = 0; i < max_NM; ++i) {
    auto minCol = *std::min_element(matrix.begin() + max_NM * i,
                                    matrix.begin() + max_NM * (i + 1));
    for (int j = 0; j < max_NM; ++j)
      matrix[j + i * max_NM] -= minCol;
  }

  return 2;
}

// (rowwise) subtract minimal value in each nth row
int Hungarian::step2() {
  std::transform(mask_matrix.begin(), mask_matrix.end(), mask_matrix.begin(),
                 [](bool value) { return false; });

  auto maxVal = *std::max_element(matrix.begin(), matrix.end());
  for (int j = 0; j < max_NM; ++j) {
    float minVal = maxVal;
    for (int i = 0; i < max_NM; ++i)
      minVal = std::min(minVal, matrix[i * max_NM + j]);
    for (int i = 0; i < max_NM; ++i) {
      matrix[i * max_NM + j] -= minVal;
      if (matrix[i * max_NM + j] == 0)
        mask_matrix[i * max_NM + j] = true;
    }
  }

  return 3;
}

struct cmp {
  bool operator()(std::pair<int, int> a, std::pair<int, int> b) {
    return a.second < b.second;
  }
};

bool dfs_zero(int *count, int count_tmp, std::vector<bool> &mask_matrix,
              std::vector<bool> lx_mask, std::vector<bool> ly_mask, const int N,
              const int M) {
  std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>,cmp> cnts;
  int before_cnt = 0;
  std::vector<bool> mask_origin;
  mask_origin.assign(mask_matrix.begin(), mask_matrix.end());

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

  // erase zero line.
  // ----------------
  // cnt.first is index of x,y line, cnt.second is the number of zero value in
  // cnt.first's line.
  while (!cnts.empty()) {
    auto cnt = cnts.top();
    cnts.pop();

    if (cnt.second == 0)
      return true;

    // columnwise (cnt.first / N == 1), rowwise (cnt.first / N == 0)
    if (cnt.first / N == 1) {
      int j_idx = cnt.first % M;
      ly_mask[j_idx] = true;
      for (int i = 0; i < N; ++i)
        mask_matrix[j_idx + i * M] = 0;
    } else {
      int i_idx = cnt.first % M;
      lx_mask[i_idx] = true;
      for (int j = 0; j < M; ++j)
        mask_matrix[j + i_idx * M] = 0;
    }

    auto status = dfs_zero(count, count_tmp + 1, mask_matrix, lx_mask, ly_mask, N, M);
    if (status && count_tmp + 1 < N) {
      // minimal counts of lines is smaller than N
      *count = count_tmp + 1;
      break;
    } else if (status && count_tmp + 1 == N) {
      // minimal counts of lines is same with N
      *count = *count > count_tmp + 1 ? count_tmp + 1
                                      : (*count != 0 ? *count : count_tmp + 1);
      return false;
    } else {
      auto next_cnt = cnts.top();
      if (next_cnt.second < cnt.second || next_cnt.second <= 1)
        break;
      // recover before matrix
      mask_matrix.assign(mask_origin.begin(), mask_origin.end());
    }
  }
  return false;
}

// find zero line
int Hungarian::step3() {
  int count = 0;
  auto status =
      dfs_zero(&count, 0, mask_matrix, lx_mask, ly_mask, max_NM, max_NM);

  if (count != max_NM)
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
  float delta = *std::max_element(matrix.begin(), matrix.end());
  for (int i = 0; i < max_NM; ++i) {
    for (int j = 0; j < max_NM; ++j) {
      auto tmp = std::min(delta, matrix[j + i * max_NM]);
      delta = tmp == 0 ? delta : tmp;
    }
  }
  //subtract delta from each element.
  std::transform(matrix.begin(), matrix.end(), matrix.begin(),
                 [delta](float c) { return c - delta; });
  //add delta to covered element.
  for (int i = 0; i < max_NM; ++i)
    for (int j = 0; j < max_NM; ++j)
      matrix[j + i * max_NM] = std::max(
          0.f, matrix[j + i * max_NM] + delta * (lx_mask[i] + ly_mask[j]));

  return 2;
}

bool DFS_step5_columnwise(std::vector<float> &mask_matrix,
                          std::vector<bool> m_mask, const int N, const int M,
                          const int N_idx, std::vector<int> &assignment_idx) {

  if (N_idx >= N || std::count(m_mask.begin(), m_mask.end(), true) == M)
    return true;

  const int NM_MAX = std::max(N, M);

  for (int j = 0; j < NM_MAX; ++j) {
    if (mask_matrix[j + NM_MAX * N_idx] == 0 && !m_mask[j]) {
      m_mask[j] = true;
      assignment_idx.push_back(j);

      auto status = DFS_step5_columnwise(mask_matrix, m_mask, N, M, N_idx + 1,
                                         assignment_idx);
      if (!status) {
        m_mask[j] = false;
        if (!assignment_idx.empty())
          assignment_idx.pop_back();
        continue;
      } else {
        return true;
      }
    }
  }
  return false;
}

bool DFS_step5_rowwise(std::vector<float> &mask_matrix,
                       std::vector<bool> n_mask, const int N, const int M,
                       const int M_idx, std::vector<int> &assignment_idx) {

  if (M_idx >= M || std::count(n_mask.begin(), n_mask.end(), true) == N)
    return true;

  const int NM_MAX = std::max(N, M);

  for (int i = 0; i < NM_MAX; ++i) {
    if (mask_matrix[M_idx + NM_MAX * i] == 0 && !n_mask[i]) {
      n_mask[i] = true;
      assignment_idx.push_back(M_idx);

      auto status = DFS_step5_rowwise(mask_matrix, n_mask, N, M, M_idx + 1,
                                      assignment_idx);
      if (!status) {
        n_mask[i] = false;
        if (!assignment_idx.empty())
          assignment_idx.pop_back();
        continue;
      } else {
        return true;
      }
    }
  }
  return false;
}

int Hungarian::step5(const float **cost, float *score,
                     std::vector<float> *assignment_idx) {

  // Fill lx_mask and ly_mask with false
  std::transform(lx_mask.begin(), lx_mask.end(), lx_mask.begin(),
                 [](bool num) { return false; });
  std::transform(ly_mask.begin(), ly_mask.end(), ly_mask.begin(),
                 [](bool num) { return false; });

  std::vector<int> assign_idx;
  auto status = DFS_step5_columnwise(matrix, lx_mask, N, M, 0, assign_idx);

  if (!status) {
    assign_idx.clear();
    DFS_step5_rowwise(matrix, ly_mask, N, M, 0, assign_idx);
  }

  int i = 0;
  for (auto idx : assign_idx) {
    assignment_idx->push_back(idx);
    *score += (*cost)[i * M + idx];
    ++i;
  }

  return 0;
}

float Hungarian::solve(const float **cost, const int N, const int M,
                       const int _MODE, std::vector<float> &assignment_index) {

  init(cost, N, M, _MODE);

  std::vector<float> assign_idx;

  std::vector<float> *assign_idx_ptr = &assignment_index;

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
      step = step5(cost, &score, assign_idx_ptr);
      break;
    }
  } while (step != 0);

  std::cout << "Total Cost : " << score << ", assignment_index = { ";
  for (auto a : assignment_index)
    std::cout << a << " ";
  std::cout << "}" << std::endl;

  return 0;
}