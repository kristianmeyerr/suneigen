#include <iostream>

int main() {
    double* leak = new double[10];
    std::cout << leak << std::endl;
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
