# 'Google_test' is the subproject name

# 'lib' is the folder with Google Test sources
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(test_test test_test.cpp)
add_test(NAME test_test COMMAND test_test.cpp)
target_link_libraries(test_test SLAE gtest gtest_main pthread)