#include <iostream>
#include "Eigen/Dense"

using namespace Eigen;
using namespace std;
// TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or
// click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.
int main() {
    cout << "Hello,world" << endl;
    VectorXd aa = VectorXd ::LinSpaced(5,1,5);
    cout << aa << endl;
    return 0;
}