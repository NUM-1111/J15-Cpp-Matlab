#include "tool.h"
#include <cmath>
#include <stdexcept>
#include <Eigen/Dense>
#include <vector>
#include <unsupported/Eigen/FFT>
#include <algorithm>



// 基于Eigen的fir1实现
// N: 滤波器阶数，Wn: 归一化频率（[w1, w2]，范围0~1，带通/低通/高通）
Eigen::VectorXd fir1(int N, const std::vector<double>& Wn) {
    // MATLAB fir1 默认窗函数为汉明窗
    int M = N;
    Eigen::VectorXd h = Eigen::VectorXd::Zero(M + 1);
    double fc1 = 0, fc2 = 0;
    bool bandpass = (Wn.size() == 2);
    bool lowpass = (Wn.size() == 1);
    bool highpass = false;
    if (bandpass) {
        fc1 = Wn[0];
        fc2 = Wn[1];
    } else if (lowpass) {
        fc2 = Wn[0];
        fc1 = 0;
    } else if (Wn.size() == 1 && Wn[0] > 0.5) {
        highpass = true;
        fc1 = Wn[0];
        fc2 = 0.5;
    } else {
        throw std::invalid_argument("Wn必须为1个或2个元素");
    }
    // 生成理想脉冲响应
    for (int n = 0; n <= M; ++n) {
        double k = (double)n - (double)M / 2.0;
        double val = 0.0;
        if (bandpass) {
            if (k == 0.0) {
                val = fc2 - fc1;
            } else {
                val = (sin(pi * fc2 * k) - sin(pi * fc1 * k)) / (pi * k);
            }
        } else if (lowpass) {
            if (k == 0.0) {
                val = fc2;
            } else {
                val = sin(pi * fc2 * k) / (pi * k);
            }
        } else if (highpass) {
            if (k == 0.0) {
                val = 1 - fc1;
            } else {
                val = -sin(pi * fc1 * k) / (pi * k);
            }
        }
        h(n) = val;
    }
    // 汉明窗
    for (int n = 0; n <= M; ++n) {
        h(n) *= 0.54 - 0.46 * cos(2 * pi * n / M);
    }
    return h;
}

// hilbert变换（对每列）
Eigen::MatrixXcd hilbert(const Eigen::MatrixXd& x) {
    // 对每列做Hilbert变换
    int N = x.rows();
    int M = x.cols();
    Eigen::MatrixXcd out(N, M);
    for (int col = 0; col < M; ++col) {
        Eigen::VectorXd v = x.col(col);
        Eigen::FFT<double> fft;
        Eigen::VectorXcd V;
        fft.fwd(V, v);
        Eigen::VectorXcd H = Eigen::VectorXcd::Zero(N);
        if (N % 2 == 0) {
            H(0) = 1;
            H(N/2) = 1;
            for (int i = 1; i < N/2; ++i) H(i) = 2;
        } else {
            H(0) = 1;
            for (int i = 1; i <= (N-1)/2; ++i) H(i) = 2;
        }
        V = V.cwiseProduct(H);
        Eigen::VectorXcd v_analytic;
        fft.inv(v_analytic, V);
        out.col(col) = v_analytic;
    }
    return out;
}

// FIR滤波器（对每列）
Eigen::MatrixXcd filter(const Eigen::VectorXd& b, const Eigen::MatrixXcd& x) {
    int N = x.rows();
    int M = x.cols();
    int L = b.size();
    Eigen::MatrixXcd y = Eigen::MatrixXcd::Zero(N, M);
    for (int m = 0; m < M; ++m) {
        for (int n = 0; n < N; ++n) {
            std::complex<double> acc = 0.0;
            for (int k = 0; k < L; ++k) {
                if (n - k >= 0)
                    acc += b(k) * x(n - k, m);
            }
            y(n, m) = acc;
        }
    }
    return y;
}

// 简单峰值检测，返回最大num_peaks个峰的索引
std::vector<int> findpeaks(const Eigen::VectorXd& x, int num_peaks) {
    std::vector<std::pair<double, int>> peaks;
    for (int i = 1; i < x.size() - 1; ++i) {
        if (x(i) > x(i-1) && x(i) > x(i+1)) {
            peaks.emplace_back(x(i), i);
        }
    }
    std::sort(peaks.begin(), peaks.end(), std::greater<>());
    std::vector<int> idx;
    for (int i = 0; i < min(num_peaks, (int)peaks.size()); ++i) {
        idx.push_back(peaks[i].second);
    }
    return idx;
}

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
) {
    double db = Bandwidth / 5.0;
    int Number_of_All = xl.rows();
    int Number_of_Sensor = xl.cols();
    int N_theta = theta.size();
    Eigen::MatrixXcd xlh = hilbert(xl); // shape: (Number_of_All, Number_of_Sensor)
    Eigen::VectorXd P_music = Eigen::VectorXd::Zero(N_theta);
    for (double f = Frequency_of_Center - Bandwidth / 2.0; f <= Frequency_of_Center + Bandwidth / 2.0; f += db) {
        std::vector<double> Wn = { (f - db / 2.0) / (fs / 2.0), (f + db / 2.0) / (fs / 2.0) };
        Eigen::VectorXd b = fir1(128, Wn);
        Eigen::MatrixXcd yl = filter(b, xlh); // shape: (Number_of_All, Number_of_Sensor)
        Eigen::MatrixXcd R = Eigen::MatrixXcd::Zero(Number_of_Sensor, Number_of_Sensor);
        for (int i = 0; i < Number_of_All; ++i) {
            R += yl.row(i).adjoint().transpose() * yl.row(i);
        }
        R /= Number_of_All;
        Eigen::MatrixXcd atheta(Number_of_Sensor, N_theta);
        for (int t = 0; t < N_theta; ++t) {
            for (int m = 0; m < Number_of_Sensor; ++m) {
                atheta(m, t) = std::exp(std::complex<double>(0, 2 * pi * f / c * d * m * std::sin(theta(t) * pi / 180.0)));
            }
        }
        Eigen::JacobiSVD<Eigen::MatrixXcd> svd(R, Eigen::ComputeFullU);
        Eigen::MatrixXcd U = svd.matrixU();
        Eigen::MatrixXcd Un = U.rightCols(Number_of_Sensor - Number_of_Angle);
        for (int t = 0; t < N_theta; ++t) {
            Eigen::VectorXcd a = atheta.col(t);
            std::complex<double> denom = (a.adjoint() * Un * Un.adjoint() * a)(0, 0);
            P_music(t) += 1.0 / std::abs(denom);
        }
    }
    std::vector<int> peak_indices = findpeaks(P_music, Number_of_Angle);
    Angle_of_Arrival.resize(Number_of_Angle);
    Angle_of_Spectrum.resize(Number_of_Angle);
    for (int i = 0; i < Number_of_Angle; ++i) {
        Angle_of_Arrival(i) = theta(peak_indices[i]);
        Angle_of_Spectrum(i) = P_music(peak_indices[i]);
    }
}

// 复数double向量与实数double向量卷积
Eigen::VectorXcd conv(const Eigen::VectorXcd& A, const Eigen::VectorXd& B) {
    int nA = A.size();
    int nB = B.size();
    int nC = nA + nB - 1;
    Eigen::VectorXcd C = Eigen::VectorXcd::Zero(nC);
    for (int i = 0; i < nC; ++i) {
        for (int j = max(0, i - nB + 1); j <= min(i, nA - 1); ++j) {
            C(i) += A(j) * B(i - j);
        }
    }
    return C;
}

// 高效FFT卷积（复数double向量与实数double向量）
Eigen::VectorXcd fft_conv(const Eigen::VectorXcd& A, const Eigen::VectorXd& B) {
    int nA = A.size();
    int nB = B.size();
    int nC = nA + nB - 1;
    int nFFT = 1;
    while (nFFT < nC) nFFT <<= 1;
    Eigen::VectorXcd A_pad = Eigen::VectorXcd::Zero(nFFT);
    Eigen::VectorXcd B_pad = Eigen::VectorXcd::Zero(nFFT);
    A_pad.head(nA) = A;
    for (int i = 0; i < nB; ++i) B_pad(i) = B(i);
    Eigen::FFT<double> fft;
    Eigen::VectorXcd FA, FB;
    fft.fwd(FA, A_pad);
    fft.fwd(FB, B_pad);
    Eigen::VectorXcd FC = FA.cwiseProduct(FB);
    Eigen::VectorXcd C_pad;
    fft.inv(C_pad, FC);
    return C_pad.head(nC);
} 