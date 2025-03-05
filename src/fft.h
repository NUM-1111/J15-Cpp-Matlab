#ifndef PRONAME_FFT_H
#define PRONAME_FFT_H

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <stdexcept>

using namespace Eigen;

namespace FFTLibrary {

    // 一维FFT（无指定长度）
    VectorXcd fft(const VectorXd& input) {
        FFT<double> fft;
        VectorXcd output;
        fft.fwd(output, input);
        return output;
    }

    // 一维FFT（指定长度n）
    VectorXcd fft(const VectorXd& input, int n) {
        VectorXd temp = VectorXd::Zero(n);
        temp.head(std::min(n, (int)input.size())) = input.head(std::min(n, (int)input.size())); //

        FFT<double> fft;
        VectorXcd output;
        fft.fwd(output, temp);
        return output;
    }


    // 多维FFT（指定长度n和维度dim）
    MatrixXcd fft(const MatrixXd& input, int n, int dim) {
        if (dim != 1 && dim != 2) {
            throw std::invalid_argument("dim must be 1 (column-wise) or 2 (row-wise)");
        }

        FFT<double> fft;
        MatrixXcd output;

        if (dim == 1) {  // 按列处理
            output.resize(n, input.cols());
            for (int col = 0; col < input.cols(); ++col) {
                VectorXd column = input.col(col);
                VectorXd temp = VectorXd::Zero(n);
                temp.head(std::min(n, (int)column.size())) = column.head(std::min(n, (int)column.size()));

                VectorXcd fft_col;
                fft.fwd(fft_col, temp);
                output.col(col) = fft_col;
            }
        } else {  // 按行处理
            output.resize(input.rows(), n);
            for (int row = 0; row < input.rows(); ++row) {
                VectorXd row_data = input.row(row);
                VectorXd temp = VectorXd::Zero(n);
                temp.head(std::min(n, (int)row_data.size())) = row_data.head(std::min(n, (int)row_data.size()));

                VectorXcd fft_row;
                fft.fwd(fft_row, temp);
                output.row(row) = fft_row;
            }
        }

        return output;
    }
}

#endif // PRONAME_FFT_H