#include <iostream>
#include <Eigen/Dense>
#include<unsupported/Eigen/FFT>
#include "fft.h"

int main() {
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

    return 0;
}