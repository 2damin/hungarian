#ifndef _hungarian_hpp_
#define _hungarian_hpp_

#include <vector>

class Hungarian{
    std::vector<float> lx, ly;
    std::vector<int> lx_mask, ly_mask;
    std::vector<float> matrix;
    int matched_jobs;
    std::vector<int> alter_paths;
    int MODE;

    public:
    void init(const float **cost, const int N, const int M);

    float solve(const float **cost, const int N, const int M, const int _MODE, float *assignment_index);

    private:
    bool augment(const float **cost, const int N, const int M);
};

#endif //_hungarian_hpp_