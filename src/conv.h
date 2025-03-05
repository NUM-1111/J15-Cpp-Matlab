//
// Created by 29436 on 2025/3/5.
//

#ifndef PRONAME_CONV_H
#define PRONAME_CONV_H

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <stdexcept>

using namespace Eigen;
using namespace std;

namespace ConvLibrary {
    VectorXd conv(const VectorXd& firstInput,const VectorXd& secondInput){
        // 获取两信号的长度
        int n = firstInput.size() + secondInput.size() - 1;

        // 扩展两个信号的大小以适应卷积结果
        VectorXd extendedFirstInput = firstInput;
        VectorXd extendedSecondInput = secondInput;

        extendedFirstInput.conservativeResize(n);
        extendedSecondInput.conservativeResize(n);

        // 在频域中进行FFT
        FFT<double> fft;
        VectorXcd F = fft.fwd(extendedFirstInput);
        VectorXcd S = fft.fwd(extendedSecondInput);

        // 在频域中进行点乘
        VectorXcd FResult = F.array() * S.array();

        // 逆FFT得到卷积结果
        VectorXd TResult = fft.inv(FResult).real();

        // 返回时域卷积结果
        return TResult;
    }
}

#endif //PRONAME_CONV_H
