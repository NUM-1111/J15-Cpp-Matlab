//
// Created by 29436 on 2025/3/4.
//
#include <iostream>
#include <Eigen/Dense>  // 引入 Eigen 头文件
using namespace Eigen;
using namespace std;

int main() {
    //矩阵与向量的操作
    MatrixXd A(3,3);
    A << 1,2,3,
        4,5,6,
        7,8,9;
    MatrixXd B = A.transpose();  //转置
    MatrixXd C = A.inverse();  //求逆
    double detA = A.determinant(); //行列式--?
    MatrixXd D = A*B;  //矩阵乘法
    VectorXd  v = A.row(0); //获取第一行

    cout << A << endl;
    cout << B << endl;
    cout << C << endl;
    cout << D << endl;
    cout << v << endl;
    cout << detA << endl;

    cout <<"---------------------------------------------"<< endl;

    //特征值的分解
    EigenSolver<MatrixXd> solver(A);
    VectorXd eigenvalues = solver.eigenvalues().real();//特征值--?
    MatrixXd eigenvectors = solver.eigenvectors().real();//特征向量--?
    cout << eigenvalues << endl << eigenvectors<<endl;

    cout << "--------------------------------------------"<< endl;

    //SVD分解--?
    JacobiSVD<MatrixXd> svd(A,ComputeThinU | ComputeThinV);
    MatrixXd U = svd.matrixU();
    VectorXd S = svd.singularValues();
    MatrixXd V = svd.matrixV();
    cout << U << endl << S << endl << V << endl;
    return 0;
}
