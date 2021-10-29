#include <iostream>
#include "NVectorEigen.hpp"


int main() {

    auto ret = suneigen::constructEmptyVectorXd(5);
    std::cout << "Hello, World!" << ret << std::endl;
    return 0;
}
