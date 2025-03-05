#include <iostream>
#include <Eigen/Dense>
#include<unsupported/Eigen/FFT>
#include "fft.h"
#include "ifft.h"
#include "conv.h"

using namespace IFFTLibrary;
using namespace FFTLibrary;
using namespace ConvLibrary;

int main() {
    /*
    // 测试一维FFT
    VectorXd vec(4);
    vec << 1, 2, 3, 4;
    VectorXcd vec_fft = FFTLibrary::fft(vec,6);
    std::cout << "一维FFT（n=6）:\n" << vec_fft << std::endl;

    // 测试矩阵按列FFT
    MatrixXd mat(3, 2);
    mat << 1, 4,
            2, 5,
            3, 6;
    MatrixXcd mat_fft_col = FFTLibrary::fft(mat, 4, 1);
    std::cout << "矩阵按列FFT（n=4）:\n" << mat_fft_col << std::endl;

     */

    // ifftxian测试数据
//    VectorXcd input(4);
//    input << std::complex<double>(1, 0), std::complex<double>(2, 0), std::complex<double>(3, 0), std::complex<double>(4, 0);
//
//    // 打印原始数据
//    cout << "Original Input:\n" << input << "\n\n";
//
//    // 计算 IFFT
//    VectorXcd ifft_result = ifft(input);
//    cout << "IFFT Result:\n" << ifft_result << "\n";
//
//    // 指定长度 IFFT
//    VectorXcd ifft_result_n = ifft(input, 6);
//    cout << "IFFT Result with specified length 6:\n" << ifft_result_n << "\n";


// 设计测试数据：3x4 矩阵
//    MatrixXcd input(3, 4);
//    input << complex<double>(1, 0), complex<double>(2, 0), complex<double>(3, 0), complex<double>(4, 0),
//            complex<double>(5, 0), complex<double>(6, 0), complex<double>(7, 0), complex<double>(8, 0),
//            complex<double>(9, 0), complex<double>(10, 0), complex<double>(11, 0), complex<double>(12, 0);
//
//    // 打印原始矩阵
//    cout << "Original Input (3x4 Matrix):\n" << input << "\n\n";
//
//    // 计算 IFFT（按行）
//    MatrixXcd ifft_result_row = ifft(input, 4, 2);  // 4 是指定的长度
//    cout << "IFFT Result (Row-wise):\n" << ifft_result_row << "\n";
//
//    // 计算 IFFT（按列）
//    MatrixXcd ifft_result_col = ifft(input, 3, 1);  // 4 是指定的长度
//    cout << "IFFT Result (Column-wise):\n" << ifft_result_col << "\n";


    // 输入信号和核
    VectorXd x(4);
    x << 1, 2, 3, 4;  // 输入信号 x(t)

    VectorXd h(4);
    h << 1, 1, 1, 1;  // 卷积核 h(t)

    VectorXd result = conv(x,h);

    cout << result << endl;

    return 0;
}