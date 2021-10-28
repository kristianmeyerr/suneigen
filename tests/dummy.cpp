#include "doctest/doctest.h"
#include "example.h"
#include <iostream>

// Tests that don't naturally fit in the headers/.cpp files directly
// can be placed in a tests/*.cpp file. Integration tests are a good example.

TEST_CASE("complicated integration tests could be here")
{
    double* leak = new double[10];
    std::cout << leak << std::endl;
    Dummy d;
    CHECK(d.doSomething() == true);
}
