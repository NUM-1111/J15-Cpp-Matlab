#include "tool.h"
#include <cmath>
#include <stdexcept>
#include <Eigen/Dense>
#include <vector>

// 1D

template <typename T>
std::vector<T> af_to_vector_1d(const af::array& arr) {
    if (arr.dims(1) > 1 || arr.dims(2) > 1 || arr.dims(3) > 1)
        throw std::invalid_argument("Array is not 1D");
    std::vector<T> result(arr.elements());
    arr.as(f32).host(result.data());
    return result;
}

template std::vector<float> af_to_vector_1d<float>(const af::array&);
template std::vector<double> af_to_vector_1d<double>(const af::array&);

// 2D

template <typename T>
std::vector<std::vector<T>> af_to_vector_2d(const af::array& arr) {
    if (arr.dims(2) > 1 || arr.dims(3) > 1)
        throw std::invalid_argument("Array is not 2D");
    dim4 d = arr.dims();
    std::vector<T> flat(d[0] * d[1]);
    arr.as(f32).host(flat.data());
    std::vector<std::vector<T>> result(d[1], std::vector<T>(d[0]));
    for (dim_t y = 0; y < d[1]; ++y)
        for (dim_t x = 0; x < d[0]; ++x)
            result[y][x] = flat[y * d[0] + x];
    return result;
}

template std::vector<std::vector<float>> af_to_vector_2d<float>(const af::array&);
template std::vector<std::vector<double>> af_to_vector_2d<double>(const af::array&);

// 3D

template <typename T>
std::vector<std::vector<std::vector<T>>> af_to_vector_3d(const af::array& arr) {
    if (arr.dims(3) > 1)
        throw std::invalid_argument("Array is not 3D");
    dim4 d = arr.dims();
    std::vector<T> flat(arr.elements());
    arr.as(f32).host(flat.data());
    std::vector<std::vector<std::vector<T>>> result(
            d[2], std::vector<std::vector<T>>(d[1], std::vector<T>(d[0]))
    );
    for (dim_t z = 0; z < d[2]; ++z)
        for (dim_t y = 0; y < d[1]; ++y)
            for (dim_t x = 0; x < d[0]; ++x) {
                size_t idx = x + y * d[0] + z * d[0] * d[1];
                result[z][y][x] = flat[idx];
            }
    return result;
}

template std::vector<std::vector<std::vector<float>>> af_to_vector_3d<float>(const af::array&);
template std::vector<std::vector<std::vector<double>>> af_to_vector_3d<double>(const af::array&);

// sinc
array sinc(const array& x) {
    const float PI = 3.14159265358979323846f;
    array pix = PI * x;
    array result = sin(pix) / (pix + 1e-12f);
    result(x == 0) = 1.0f;
    return result;
}

// besselI0_approx
static float besselI0_approx(float x) {
    float ax = std::abs(x);
    float y, result;
    if (ax < 3.75f) {
        y = (x / 3.75f) * (x / 3.75f);
        result = 1.0f + y * (3.5156229f + y * (3.0899424f + y * (1.2067492f + y * (0.2659732f + y * (0.0360768f + y * 0.0045813f)))));
    } else {
        y = 3.75f / ax;
        result = (std::exp(ax) / std::sqrt(ax)) * (0.39894228f + y * (0.01328592f + y * (0.00225319f + y * (-0.00157565f + y * (0.00916281f + y * (-0.02057706f + y * (0.02635537f + y * (-0.01647633f + y * 0.00392377f))))))));
    }
    return result;
}

// kaiser_window
array kaiser_window(int len, float beta) {
    array n = iota(dim4(len));
    n -= (len - 1) / 2.0f;
    array ratio = 2.0f * n / (len - 1);
    array w(len);
    float denom = besselI0_approx(beta);
    for (int i = 0; i < len; ++i) {
        float val = beta * std::sqrt(1.0f - ratio(i).scalar<float>() * ratio(i).scalar<float>());
        w(i) = besselI0_approx(val) / denom;
    }
    return w;
}

// design_filter
array design_filter(int filter_len, float cutoff, int p, float beta) {
    array n = iota(dim4(filter_len)) - (filter_len - 1) / 2.0f;
    array h = 2 * cutoff * sinc(2 * cutoff * n);
    array w = kaiser_window(filter_len, beta);
    h *= w;
    return h * p / sum<float>(h);
}

// zero_insert
array zero_insert(const array& x, int p) {
    int len = x.dims(0);
    int new_len = len * p;
    array y = constant(0, new_len);
    array idx = iota(dim4(len)) * p;
    y(idx) = x;
    return y;
}

// 3D中值滤波器（MATLAB兼容）
array medfilt3_matlab_compatible(const array& input_array, int m, int n, int p) {
    if (m % 2 == 0 || n % 2 == 0 || p % 2 == 0) {
        throw std::invalid_argument("medfilt3 kernel size must be odd in all dimensions.");
    }
    dim4 dims = input_array.dims();
    int rows  = dims[0];
    int cols  = dims[1];
    int pages = dims[2];
    std::vector<double> host_data(rows * cols * pages);
    input_array.host(host_data.data());
    std::vector<double> result_data = host_data;
    std::vector<double> neighborhood_values;
    neighborhood_values.reserve(m * n * p);
    for (int z = 0; z < pages; ++z) {
        for (int y = 0; y < rows; ++y) {
            for (int x = 0; x < cols; ++x) {
                neighborhood_values.clear();
                for (int dz = -p / 2; dz <= p / 2; ++dz) {
                    for (int dy = -n / 2; dy <= n / 2; ++dy) {
                        for (int dx = -m / 2; dx <= m / 2; ++dx) {
                            int actual_x = x + dx;
                            int actual_y = y + dy;
                            int actual_z = z + dz;
                            if (actual_x < 0) actual_x = -actual_x - 1;
                            else if (actual_x >= cols) actual_x = 2 * cols - actual_x - 1;
                            if (actual_y < 0) actual_y = -actual_y - 1;
                            else if (actual_y >= rows) actual_y = 2 * rows - actual_y - 1;
                            if (actual_z < 0) actual_z = -actual_z - 1;
                            else if (actual_z >= pages) actual_z = 2 * pages - actual_z - 1;
                            int idx = actual_y + actual_x * rows + actual_z * rows * cols;
                            neighborhood_values.push_back(host_data[idx]);
                        }
                    }
                }
                std::nth_element(
                        neighborhood_values.begin(),
                        neighborhood_values.begin() + neighborhood_values.size() / 2,
                        neighborhood_values.end()
                );
                double median = neighborhood_values[neighborhood_values.size() / 2];
                int out_idx = y + x * rows + z * rows * cols;
                result_data[out_idx] = median;
            }
        }
    }
    array output(rows, cols, pages, result_data.data());
    return output;
}

// 基于Eigen的fir1实现
// N: 滤波器阶数，Wn: 归一化频率（[w1, w2]，范围0~1，带通/低通/高通）
Eigen::VectorXd fir1(int N, const std::vector<double>& Wn) {
    // MATLAB fir1 默认窗函数为汉明窗
    int M = N;
    Eigen::VectorXd h = Eigen::VectorXd::Zero(M + 1);
    double fc1 = 0, fc2 = 0;
    bool bandpass = (Wn.size() == 2);
    bool lowpass = (Wn.size() == 1);
    bool highpass = false;
    if (bandpass) {
        fc1 = Wn[0] / 2.0; // MATLAB归一化到[0,1]，fir1内部再/2
        fc2 = Wn[1] / 2.0;
    } else if (lowpass) {
        fc2 = Wn[0] / 2.0;
        fc1 = 0;
    } else if (Wn.size() == 1 && Wn[0] > 0.5) { // 高通
        highpass = true;
        fc1 = Wn[0] / 2.0;
        fc2 = 0.5;
    } else {
        throw std::invalid_argument("Wn必须为1个或2个元素");
    }
    // 生成理想脉冲响应
    for (int n = 0; n <= M; ++n) {
        int k = n - M / 2;
        double val = 0.0;
        if (bandpass) {
            if (k == 0) {
                val = 2 * (fc2 - fc1);
            } else {
                val = (sin(2 * M_PI * fc2 * k) - sin(2 * M_PI * fc1 * k)) / (M_PI * k);
            }
        } else if (lowpass) {
            if (k == 0) {
                val = 2 * fc2;
            } else {
                val = sin(2 * M_PI * fc2 * k) / (M_PI * k);
            }
        } else if (highpass) {
            if (k == 0) {
                val = 1 - 2 * fc1;
            } else {
                val = -sin(2 * M_PI * fc1 * k) / (M_PI * k);
            }
        }
        h(n) = val;
    }
    // 汉明窗
    for (int n = 0; n <= M; ++n) {
        h(n) *= 0.54 - 0.46 * cos(2 * M_PI * n / M);
    }
    return h;
} 