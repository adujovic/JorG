#include "../../../asa/solver/solver.h"

int main(int argc, char** argv){
    if(argc < 6){
        std::cerr<<"No files given!"<<std::endl;
        std::cerr<<"One should provide:"<<std::endl;
        std::cerr<<"    (  i) basis"<<std::endl;
        std::cerr<<"    ( ii) supercell"<<std::endl;
        std::cerr<<"    (iii) flippable"<<std::endl;
        std::cerr<<"    ( iv) reference (i.e. 48)"<<std::endl;
        std::cerr<<"    (  v) number of unique flips (i.e. 2)"<<std::endl;
        exit(-1);
    }

    return solver(argv[1],argv[2],argv[3],std::atoi(argv[4]),std::atoi(argv[5]));
}
