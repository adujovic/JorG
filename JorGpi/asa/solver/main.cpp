#include <iostream>
#include <cstdlib>
#include "solver.h"

int main(int argc, char** argv){
    if (argc < 6) {
        std::cerr<<"main() missing some of the required arguments: basis, supercell, flippable, reference, number of unique spins."<<std::endl;
        exit(-1);
    }
    if (argc == 6) {
        return solver(argv[1],argv[2],argv[3],std::atoi(argv[4]),std::atoi(argv[5]));
    }
    return solver(argv[1],argv[2],argv[3],std::atoi(argv[4]),std::atoi(argv[5]),std::atoi(argv[6]));
}
