#include <arrayfire.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <Eigen/Dense>
#include <complex>
#include "tool.h"
#include <chrono>

using namespace af;

int main() {
    /*
    // AoA_estimate 测试
    int Number_of_All = 200;
    int Number_of_Sensor = 4;
    int Number_of_Angle = 1;
    double Frequency_of_Center = 3000;
    double Bandwidth = 1600;
    double c = 340.0;
    double d = 0.05;
    double fs = 48000;
    Eigen::VectorXd theta = Eigen::VectorXd::LinSpaced(181, -90, 90); // -90~90度
    // 构造一个简单的信号（单一入射角）
    double true_angle = 30.0;
    Eigen::MatrixXd xl(Number_of_All, Number_of_Sensor);
    for (int n = 0; n < Number_of_All; ++n) {
        for (int m = 0; m < Number_of_Sensor; ++m) {
            double t = n / fs;
            double phase = 2 * M_PI * Frequency_of_Center * t + 2 * M_PI * Frequency_of_Center / c * d * m * std::sin(true_angle * M_PI / 180.0);
            xl(n, m) = std::sin(phase);
        }
    }
    Eigen::VectorXd Angle_of_Arrival, Angle_of_Spectrum;
    AoA_estimate(xl, theta, Frequency_of_Center, Bandwidth, c, d, Number_of_Angle, fs, Angle_of_Arrival, Angle_of_Spectrum);
    std::cout << "估计到的入射角: ";
    for (int i = 0; i < Angle_of_Arrival.size(); ++i) {
        std::cout << Angle_of_Arrival(i) << " ";
    }
    std::cout << std::endl;
    std::cout << "对应谱值: ";
    for (int i = 0; i < Angle_of_Spectrum.size(); ++i) {
        std::cout << Angle_of_Spectrum(i) << " ";
    }
    std::cout << std::endl;
    */

    // conv测试
    // Eigen::VectorXcd A(5);
    // A << std::complex<double>(1,1), std::complex<double>(2,-1), std::complex<double>(0,2), std::complex<double>(-1,0), std::complex<double>(0,-1);
    // Eigen::VectorXd B(3);
    // B << 1, 2, 3;
    // Eigen::VectorXcd C = conv(A, B);
    // std::cout << "conv(A, B) = ";
    // for (int i = 0; i < C.size(); ++i) {
    //     std::cout << C(i) << " ";
    // }
    // std::cout << std::endl;

    // fft_conv大规模测试
    // int N = 100000;
    // Eigen::VectorXcd A_big = Eigen::VectorXcd::Random(N);
    // Eigen::VectorXd B_big = Eigen::VectorXd::Random(N);
    // std::cout << "开始大规模fft_conv..." << std::endl;
    // auto start = std::chrono::high_resolution_clock::now();
    // Eigen::VectorXcd C_big = fft_conv(A_big, B_big);
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed = end - start;
    // std::cout << "fft_conv(A_big, B_big) 完成，输出长度: " << C_big.size() << ", 用时: " << elapsed.count() << " 秒" << std::endl;
    // std::cout << "前3个结果: " << C_big(0) << ", " << C_big(1) << ", " << C_big(2) << std::endl;

    // fft_conv小规模测试
    Eigen::VectorXcd A_small(10);
    for (int i = 0; i < 10; ++i) A_small(i) = std::complex<double>(i + 1, 0.5 * i);
    Eigen::VectorXd B_small(4);
    B_small << 1, -1, 2, 0.5;
    Eigen::VectorXcd C_small = fft_conv(A_small, B_small);
    std::cout << "fft_conv(A_small, B_small) = ";
    for (int i = 0; i < C_small.size(); ++i) {
        std::cout << C_small(i) << " ";
    }
    std::cout << std::endl;

    return 0;
}

