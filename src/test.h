#ifndef PRONAME_TEST_H
#define PRONAME_TEST_H
#include "tool.h"
#endif // PRONAME_TEST_H

#include <arrayfire.h>
#include <vector>
#include <stdexcept>

using namespace af;
using namespace std;

// 1D 转换：af::array -> std::vector<T>
template <typename T>
vector<T> af_to_vector_1d(const af::array& arr) {
    if (arr.dims(1) > 1 || arr.dims(2) > 1 || arr.dims(3) > 1)
        throw invalid_argument("Array is not 1D");

    vector<T> result(arr.elements());
    arr.as(f32).host(result.data()); // 可改为 as(f64) 或不强制转换
    return result;
}

// 2D 转换：af::array -> vector<vector<T>> （行优先）
template <typename T>
vector<vector<T>> af_to_vector_2d(const af::array& arr) {
    if (arr.dims(2) > 1 || arr.dims(3) > 1)
        throw invalid_argument("Array is not 2D");

    dim4 d = arr.dims();
    vector<T> flat(d[0] * d[1]);
    arr.as(f32).host(flat.data());

    vector<vector<T>> result(d[1], vector<T>(d[0]));
    for (dim_t y = 0; y < d[1]; ++y)
        for (dim_t x = 0; x < d[0]; ++x)
            result[y][x] = flat[y * d[0] + x];

    return result;
}

// 3D 转换：af::array -> vector<vector<vector<T>>> （Z-Y-X顺序）
template <typename T>
vector<vector<vector<T>>> af_to_vector_3d(const af::array& arr) {
    if (arr.dims(3) > 1)
        throw invalid_argument("Array is not 3D");

    dim4 d = arr.dims(); // x, y, z
    vector<T> flat(arr.elements());
    arr.as(f32).host(flat.data());

    vector<vector<vector<T>>> result(
            d[2], vector<vector<T>>(d[1], vector<T>(d[0]))
    );

    for (dim_t z = 0; z < d[2]; ++z)
        for (dim_t y = 0; y < d[1]; ++y)
            for (dim_t x = 0; x < d[0]; ++x) {
                size_t idx = x + y * d[0] + z * d[0] * d[1];
                result[z][y][x] = flat[idx];
            }

    return result;
}
