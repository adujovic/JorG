#include "sum.h"

int sum(int N,int M){
    int sum = 0;
    for(int i=N; i<=M; ++i)
        sum += i;
    return sum;
}
