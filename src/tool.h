#ifndef PRONAME_TOOL_H
#define PRONAME_TOOL_H

#include <arrayfire.h>
#include <vector>
#include <Eigen/Dense>
#include <complex>

const double pi = 3.14159265358979323846;
// 基于Eigen的fir1实现
Eigen::VectorXd fir1(int N, const std::vector<double>& Wn);

// hilbert变换（对每列）
Eigen::MatrixXcd hilbert(const Eigen::MatrixXd& x);
// FIR滤波器（对每列）
Eigen::MatrixXcd filter(const Eigen::VectorXd& b, const Eigen::MatrixXcd& x);
// 简单峰值检测，返回最大num_peaks个峰的索引
std::vector<int> findpeaks(const Eigen::VectorXd& x, int num_peaks);
// AoA估计主函数
void AoA_estimate(
    const Eigen::MatrixXd& xl,
    const Eigen::VectorXd& theta,
    double Frequency_of_Center,
    double Bandwidth,
    double c,
    double d,
    int Number_of_Angle,
    double fs,
    Eigen::VectorXd& Angle_of_Arrival,
    Eigen::VectorXd& Angle_of_Spectrum
);
// 复数double向量与实数double向量卷积
Eigen::VectorXcd conv(const Eigen::VectorXcd& A, const Eigen::VectorXd& B);
// 高效FFT卷积（复数double向量与实数double向量）
Eigen::VectorXcd fft_conv(const Eigen::VectorXcd& A, const Eigen::VectorXd& B);

#endif // PRONAME_TOOL_H 