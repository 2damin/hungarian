#ifndef _hungarian_hpp_
#define _hungarian_hpp_

class Hungarian{
    public:
    float solve(const float **Cost, const int N, const int M, const int MODE, float *assignment_index);
};

#endif //_hungarian_hpp_