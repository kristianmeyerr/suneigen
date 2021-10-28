#include <iostream>
#include <Eigen/Dense>

int main() {
    double* leak = new double[10];
    std::cout << leak << std::endl;
    std::cout << "Hello, World!" << std::endl;

    Eigen::MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    std::cout << m << std::endl;

    return 0;
}
