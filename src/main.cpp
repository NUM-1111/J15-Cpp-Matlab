#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <Eigen/Dense>
#include "tool.h"

using namespace af;

int main() {
        // fir1 测试
        int fc = 3000;
        int fsr = 48000;
        int B = 1600;
        std::vector<double> Wn = {
            2.0 * (fc - B / 2.0) / fsr,
            2.0 * (fc + B / 2.0) / fsr
        };
        Eigen::VectorXd Blo = fir1(1023, Wn);
        std::cout << "fir1 测试: Blo 前10个系数:" << std::endl;
        for (int i = 0; i < 10; ++i) {
            std::cout << Blo(i) << " ";
        }
        std::cout << std::endl;
    return 0;
}

