//
// Created by 29436 on 2025/5/18.
//

#include "resample.h"
#include <cmath>

array resample(const af::array& x, int p, int q) {
    int len = x.dims(0);
    int filter_len = 10 * std::max(p, q) + 1;
    float cutoff = 1.0f / (2.0f * std::max(p, q));
    array h = design_filter(filter_len, cutoff, p);
    array up = zero_insert(x, p);
    array filtered = convolve1(up, h, AF_CONV_EXPAND);
    int delay = (filter_len - 1) / 2;
    filtered = filtered(seq(delay, delay + up.dims(0) - 1));
    int L = static_cast<int>(std::ceil(len * p / static_cast<float>(q)));
    array idx = (iota(dim4(L)) * q).as(u32);
    unsigned int max_idx = max(idx).scalar<unsigned int>();
    if (max_idx >= filtered.dims(0)) {
        idx = idx(idx < filtered.dims(0));
    }
    return filtered(idx);
}
