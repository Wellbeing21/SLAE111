cmake_minimum_required(VERSION 3.24)
project(SLAE)

set(CMAKE_CXX_STANDARD 20)

add_executable(SLAE Tridiagonal_matrix/tridiagonal.cpp Tridiagonal_matrix/tridiagonal.h solvers/solve_tridiagonal.hpp)
set(SRC_FILE ../Tridiagonal_matrix/tridiagonal.cpp ../Tridiagonal_matrix/tridiagonal.h ../solvers/solve_tridiagonal.hpp)

enable_testing()
add_subdirectory(tests)
#target_link_libraries(SLAE gtest gmock pthread)
#target_link_libraries(SLAE gtest_main gmock_main)

#add_subdirectory(gtest)