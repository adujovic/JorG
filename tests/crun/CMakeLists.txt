cmake_minimum_required(VERSION 3.13)
 
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
 
add_executable(runTests.x gtest.cpp 
                          hello/hello.h hello/hello.cpp 
                          loop/loop.h loop/loop.cpp 
                          pi/pi.h pi/pi.cpp 
                          sum/sum.h sum/sum.cpp)
target_compile_options(runTests.x PRIVATE -std=c++17 -Wall -pedantic -Wextra -fopenmp)
target_link_libraries( runTests.x ${GTEST_LIBRARIES} pthread)
target_link_options(   runTests.x PRIVATE -fopenmp)
