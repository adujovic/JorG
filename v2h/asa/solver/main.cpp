#include <iostream>
#include <iomanip>
#include <fstream>

#include <vector>
#include <array>
#include <numeric>
#include <random>
#include <regex>

#include <cmath>
#include <cstring>

#include "../asa.h"
#include "../ising.h"
#include "../arithmeticvector.h"

int main(int argc, char** argv){
    if(argc < 5){
        std::cerr<<"No files given!"<<std::endl;
        std::cerr<<"One should provide:"<<std::endl;
        std::cerr<<"    (  i) basis"<<std::endl;
        std::cerr<<"    ( ii) supercell"<<std::endl;
        std::cerr<<"    (iii) flippable"<<std::endl;
        std::cerr<<"    ( iv) reference (i.e. 48)"<<std::endl;
        exit(-1);
    }

#ifdef _SITESNUMBER
    constexpr size_t SITESNUMBER = _SITESNUMBER;
#else
    constexpr size_t SITESNUMBER = 64;
#endif

    std::size_t pos;
    std::size_t buff;
    std::string line;
    
    /* Reading data from files:
     * crystal direction (basis) in argv[1]
     * supercell                 in argv[2]
     * unique flippable spins    in argv[3]
     * new reference point       in argv[4]
     */ 
    std::array<typename ising::IsingModel<SITESNUMBER>::VectorType,3> basis;
    std::ifstream directions(argv[1]);
    if(directions.is_open()) {
        size_t id = 0;
        while (getline(directions,line)){
            buff = 0;
            std::array<double,3> position;
            for(int i=0; i<3; ++i){
              position[i] = std::stod(line.substr(buff),&pos);
              buff += pos;
            }
            basis[id++] = position;
        }
        directions.close();
    }

    std::vector<std::tuple<size_t,typename ising::IsingModel<SITESNUMBER>::VectorType,double>> supercell; 
    std::ifstream inSupercell(argv[2]);
    if(inSupercell.is_open()) {
        while (getline(inSupercell,line)){
            size_t idx = std::stoi(line,&pos);
            buff = pos;
            std::array<double,3> position;
            for(int i=0; i<3; ++i){
              position[i] = std::stod(line.substr(buff),&pos);
              buff += pos;
            }
            auto moment = std::stod(line.substr(buff));
            supercell.push_back(std::make_tuple(idx,position,moment));
        }
        inSupercell.close();
    }

    std::vector<std::tuple<size_t,typename ising::IsingModel<SITESNUMBER>::VectorType,double>> flippable; 
    std::ifstream inFlippable(argv[3]);
    if(inFlippable.is_open()) {
        while (getline(inFlippable,line)){
            size_t idx = std::stoi(line,&pos);
            buff = pos;
            std::array<double,3> position;
            for(int i=0; i<3; ++i){
              position[i] = std::stod(line.substr(buff),&pos);
              buff += pos;
            }
            auto moment = std::stod(line.substr(buff));
            flippable.push_back(std::make_tuple(idx,position,moment));
        }
        inFlippable.close();
    }

    size_t reference = std::stoi(argv[4]);

    std::cout<<"########################################"<<std::endl;
    std::cout<<"##                                    ##"<<std::endl;
    std::cout<<"##                                    ##"<<std::endl;
    std::cout<<"##         Ising Model Solver         ##"<<std::endl;
    std::cout<<"##                                    ##"<<std::endl;
    std::cout<<"##                                    ##"<<std::endl;
    std::cout<<"########################################"<<std::endl;

    std::random_device randomDevice{};
    std::mt19937 generator{randomDevice()};
    std::normal_distribution<> gauss{0.0,1.0};

    ising::IsingModel<SITESNUMBER> model;

    model.set_basis(basis);
    model.set_supercell(flippable);
    model.set_reference(reference);

    double minDistance = (basis[0]^basis[1])*basis[2];
    double maxDistance = 0.0;
    for(const auto& atom1 : flippable){
        for(const auto& atom2 : flippable){
             auto d = ising::IsingModel<SITESNUMBER>::VectorType::norm(
                               std::get<1>(atom1) - std::get<1>(atom2) );
            if (d > 1e-7 && d < minDistance) minDistance = d;
            if (d > maxDistance)             maxDistance = d;
        }
    }
    auto decayCoeff = 24*log10(2)*log(10)/(maxDistance-minDistance);


    std::vector<std::tuple<unsigned,unsigned,double>> d;

    std::bitset<SITESNUMBER> mask(std::string(SITESNUMBER,'0'));

    for (const auto& atom : flippable) mask.set(std::get<0>(atom));

    mask.reset(reference);

    for(auto n = 0.5; n<128; n*=2){
        model.reset();
        d.clear();
        for (const auto& atom1 : flippable) {
            for (const auto& atom2 : supercell) {
                if(std::get<2>(atom2) > 0.0 && std::get<0>(atom1) != std::get<0>(atom2)) 
                    d.push_back(std::make_tuple(std::get<0>(atom1),
                                                std::get<0>(atom2),
                                                -exp(-n*decayCoeff*ising::IsingModel<SITESNUMBER>::VectorType::norm(
                                                    std::get<1>(atom1) - std::get<1>(atom2) ))));
            }
        }
        model.add_interaction(d);
        
        for(unsigned m = 0; m<5; ++m){
            model.randomize_state();
#ifdef _VEBOSE    
            std::cout<<"               ";
            for(int i=1; i<(SITESNUMBER-reference); ++i) std::cout<<" ";
            std::cout<<"| reference"<<std::endl;
            std::cout<<"Mask:          "<<mask<<std::endl;
#endif
            std::cout<<n*decayCoeff<<"  \t("<<m<<")\t";
            model.run(&mask);
        }
    }

    return 0;
}


