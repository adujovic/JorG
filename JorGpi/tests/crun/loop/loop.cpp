#include "loop.h"

int loop(int N){
    int sum = 0;
    for(int i=0; i<N; ++i){
        sum += i;
        std::cout<<i<<" ";
    }
    std::cout<<std::endl;
    return sum;
}
