//
// Created by 29436 on 2025/3/5.
//

#ifndef PRONAME_IFFT_H
#define PRONAME_IFFT_H
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <stdexcept>
using namespace Eigen;
using namespace std;

namespace IFFTLibrary {

    // 一维IFFT（无指定长度）
    VectorXcd ifft(const VectorXcd& input) {
        FFT<double> fft;
        VectorXcd output;
        fft.inv(output, input);  // 直接使用 Eigen 提供的 IFFT
        return output;
    }

    // 一维IFFT（指定长度n）
    VectorXcd ifft(const VectorXcd& input, int n) {
        VectorXcd temp = VectorXcd::Zero(n);
        temp.head(min(n, (int)input.size())) = input.head(min(n, (int)input.size()));

        FFT<double> fft;
        VectorXcd output;
        fft.inv(output, temp);
        return output;
    }


    // 多维IFFT（指定长度n和维度dim）
    MatrixXcd ifft(const MatrixXcd& input, int n, int dim) {
        if (dim != 1 && dim != 2) {
            throw invalid_argument("dim must be 1 (column-wise) or 2 (row-wise)");
        }

        FFT<double> fft;
        MatrixXcd output;

        if (dim == 1) {  // 按列处理
            output.resize(n, input.cols());
            for (int col = 0; col < input.cols(); ++col) {
                VectorXcd column = input.col(col);
                VectorXcd temp = VectorXcd::Zero(n);
                temp.head(min(n, (int)column.size())) = column.head(min(n, (int)column.size()));

                VectorXcd ifft_col;
                fft.inv(ifft_col, temp);
                output.col(col) = ifft_col;
            }
        } else {  // 按行处理
            output.resize(input.rows(), n);
            for (int row = 0; row < input.rows(); ++row) {
                VectorXcd row_data = input.row(row);
                VectorXcd temp = VectorXcd::Zero(n);
                temp.head(min(n, (int)row_data.size())) = row_data.head(min(n, (int)row_data.size()));

                VectorXcd ifft_row;
                fft.inv(ifft_row, temp);
                output.row(row) = ifft_row;
            }
        }

        return output;
    }

}

#endif //PRONAME_IFFT_H
