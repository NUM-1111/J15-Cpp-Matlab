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
        VectorXd processed = input;
        int m = n - input.size();

        if (m > 0) {
            processed.conservativeResize(n);
            processed.tail(m).setZero();
        } else if (m < 0) {
            processed.conservativeResize(n);
        }

        FFT<double> fft;
        VectorXcd output;
        fft.fwd(output, processed);
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
                VectorXcd fft_col;

                // 调整列长度到n
                if (n > column.size()) {
                    column.conservativeResize(n);
                    column.tail(n - column.size()).setZero();
                } else if (n < column.size()) {
                    column.conservativeResize(n);
                }

                fft.fwd(fft_col, column);
                output.col(col) = fft_col;
            }
        } else {  // 按行处理
            output.resize(input.rows(), n);
            for (int row = 0; row < input.rows(); ++row) {
                VectorXd row_data = input.row(row);
                VectorXcd fft_row;

                // 调整行长度到n
                if (n > row_data.size()) {
                    row_data.conservativeResize(n);
                    row_data.tail(n - row_data.size()).setZero();
                } else if (n < row_data.size()) {
                    row_data.conservativeResize(n);
                }

                fft.fwd(fft_row, row_data);
                output.row(row) = fft_row;
            }
        }

        return output;
    }
}

#endif // PRONAME_FFT_H