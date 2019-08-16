#include "hello/hello.h"
#include "pi/pi.h"
#include "loop/loop.h"
#include "sum/sum.h"
#include <gtest/gtest.h>

TEST(MandelbrotTest, JustPrint) {
    mandelbrot(78,44);
    ASSERT_EQ(1,1);
}

TEST(HelloTest, JustPrint) {
    hello();
    ASSERT_EQ(1,1);
}

TEST(LoopTest, Addition) {
    for(int i=0; i<15; ++i) ASSERT_EQ(loop(i),(i*(i-1))/2);
}

TEST(SumTest, Addition) {
    for(int i=0; i<15; ++i)
        for(int j=i+1; j<15; ++j){
            auto n = j-i+1;
            ASSERT_EQ(sum(i,j),n*i + (n*(n-1))/2);
        }
}

TEST(SumTest, WrongOrder) {
    for(int i=1; i<15; ++i)
        for(int j=i-1; j>=0; --j){
            ASSERT_EQ(sum(i,j),0);
        }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
