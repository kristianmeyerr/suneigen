find_package(Eigen3 3.4 REQUIRED)
find_package(SUNDIALS 5.8 REQUIRED COMPONENTS sundials_cvodes sundials_nvecserial)
set(GSL_LITE_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/thirdparty/gsl/")