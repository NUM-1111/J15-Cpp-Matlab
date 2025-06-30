#ifndef PRONAME_TOOL_H
#define PRONAME_TOOL_H

#include <arrayfire.h>
#include <vector>
#include <Eigen/Dense>

using namespace af;

// af::array -> std::vector<T> 1D
template <typename T>
std::vector<T> af_to_vector_1d(const af::array& arr);
// af::array -> std::vector<std::vector<T>> 2D
template <typename T>
std::vector<std::vector<T>> af_to_vector_2d(const af::array& arr);
// af::array -> std::vector<std::vector<std::vector<T>>> 3D
template <typename T>
std::vector<std::vector<std::vector<T>>> af_to_vector_3d(const af::array& arr);

// sinc函数
array sinc(const array& x);
// Kaiser窗
array kaiser_window(int len, float beta = 8.6f);
// 低通滤波器设计
array design_filter(int filter_len, float cutoff, int p, float beta = 8.6f);
// 零插值上采样
array zero_insert(const array& x, int p);

// 3D中值滤波器（MATLAB兼容）
array medfilt3_matlab_compatible(const array& input_array, int m, int n, int p);

// 基于Eigen的fir1实现
Eigen::VectorXd fir1(int N, const std::vector<double>& Wn);

#endif // PRONAME_TOOL_H 