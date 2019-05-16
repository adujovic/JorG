#include <iostream>
#include <iomanip>
#include <fstream>

#include <vector>
#include <array>
#include <unordered_set>
#include <numeric>
#include <random>
#include <regex>

#include <cmath>
#include <cstring>

#include "../asa.h"
#include "../ising.h"
#include "../arithmeticvector.h"

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
     * number of unique flips    in argv[5]
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

    size_t unique_flips = std::stoi(argv[5]);

    std::cout<<"******************************************************************************************"<<std::endl;
    std::cout<<"******************************************************************************************"<<std::endl;
    std::cout<<"***                                                                                    ***"<<std::endl; 
    std::cout<<"***       ▄▄▄▄         ▀▀█             ▀                                               ***"<<std::endl; 
    std::cout<<"***      █▀   ▀  ▄▄▄     █    ▄   ▄  ▄▄▄    ▄ ▄▄    ▄▄▄▄                               ***"<<std::endl; 
    std::cout<<"***      ▀█▄▄▄  █▀ ▀█    █    ▀▄ ▄▀    █    █▀  █  █▀ ▀█                               ***"<<std::endl; 
    std::cout<<"***          ▀█ █   █    █     █▄█     █    █   █  █   █                               ***"<<std::endl; 
    std::cout<<"***      ▀▄▄▄█▀ ▀█▄█▀    ▀▄▄    █    ▄▄█▄▄  █   █  ▀█▄▀█                               ***"<<std::endl; 
    std::cout<<"***                                                 ▄  █                               ***"<<std::endl; 
    std::cout<<"***                                                  ▀▀                                ***"<<std::endl; 
    std::cout<<"***                                                                                    ***"<<std::endl; 
    std::cout<<"***      ▄▄▄▄▄           ▀                       ▄    ▄            █         ▀▀█       ***"<<std::endl; 
    std::cout<<"***        █     ▄▄▄   ▄▄▄    ▄ ▄▄    ▄▄▄▄       ██  ██  ▄▄▄    ▄▄▄█   ▄▄▄     █       ***"<<std::endl; 
    std::cout<<"***        █    █   ▀    █    █▀  █  █▀ ▀█       █ ██ █ █▀ ▀█  █▀ ▀█  █▀  █    █       ***"<<std::endl; 
    std::cout<<"***        █     ▀▀▀▄    █    █   █  █   █       █ ▀▀ █ █   █  █   █  █▀▀▀▀    █       ***"<<std::endl; 
    std::cout<<"***      ▄▄█▄▄  ▀▄▄▄▀  ▄▄█▄▄  █   █  ▀█▄▀█       █    █ ▀█▄█▀  ▀█▄██  ▀█▄▄▀    ▀▄▄     ***"<<std::endl; 
    std::cout<<"***                                   ▄  █                                             ***"<<std::endl; 
    std::cout<<"***                                    ▀▀                                              ***"<<std::endl; 
    std::cout<<"******************************************************************************************"<<std::endl;
    std::cout<<"******************************************************************************************"<<std::endl;
                                      
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

#ifdef _VERBOSE    
            std::cout<<"               ";
            for(int i=1; i<(SITESNUMBER-reference); ++i) std::cout<<" ";
            std::cout<<"| reference"<<std::endl;
            std::cout<<"Mask:          "<<mask<<std::endl;
#endif

#ifdef _LOWERLIMIT
    size_t lowerlimit = _LOWERLIMIT;
#else
    size_t lowerlimit = 0.25*mask.count(); // doesn't allow less than 25% flipps
#endif
#ifdef _UPPERLIMIT
    size_t upperlimit = _UPPERLIMIT;
#else
    size_t upperlimit = 0.75*mask.count(); // doesn't allow more than 75% flipps
#endif


    std::unordered_set<std::bitset<SITESNUMBER>> solutions;
    size_t iteration = 0U;
    for(auto n = 4.0; n<1024.0; n*=2){
        model.reset();
        d.clear();
        for (const auto& atom1 : flippable) {
            for (const auto& atom2 : supercell) {
                if(std::fabs(std::get<2>(atom2)) > 1e-9 && std::get<0>(atom1) != std::get<0>(atom2)) 
                    d.push_back(std::make_tuple(std::get<0>(atom1),
                                                std::get<0>(atom2),
                                                -exp(-n*decayCoeff*ising::IsingModel<SITESNUMBER>::VectorType::norm(
                                                    std::get<1>(atom1) - std::get<1>(atom2) ))));
            }
        }
        model.add_interaction(d);

#ifdef _VERBOSE    
        std::cout<<"Limits are set to be: ["<<lowerlimit<<","<<upperlimit<<']'<<std::endl;
#endif

        for(unsigned m = 0; m<unique_flips; ++m){
            model.randomize_state();

            auto x = model.run(&mask);
	    if(x.count() < lowerlimit) continue;
	    if(x.count() > upperlimit) continue;

            std::cout<<"("<<iteration<<")\t"<<n*decayCoeff<<"\t";
	    for (unsigned i=0; i<reference; ++i){
		if(mask[i])  std::cout<<"\033[32m"<<x[i]<<"\033[39m";
		else         std::cout<<x[i];
	    }
            std::cout<<"\033[1m"<<x[reference]<<"\033[0m";
	    for (unsigned i=reference+1; i<SITESNUMBER; ++i){
		if(mask[i])  std::cout<<"\033[91m"<<x[i]<<"\033[39m";
		else         std::cout<<x[i];
	    }
	    std::cout<<std::endl;
	    solutions.insert(x);
            ++iteration;
	    if(solutions.size() >= unique_flips) break;
        }
    if(solutions.size() >= unique_flips) break;
    }

    std::ofstream ostrm("best.flips", std::ios::out);
    for(const auto& x : solutions){
        for(unsigned b=0U; b<SITESNUMBER; ++b){
          ostrm<<x[b]<<" ";
        }
        ostrm<<std::endl;
    }

    return 0;
}

// BF    = '\033[1m'
// IT    = '\033[3m'
// UN    = '\033[4m'
// BLINK = '\033[5m'
// END   = '\033[0m'
// INV   = '\033[7m'
// HID   = '\033[8m'
// DEFAULT       = "\033[39m"
// BLACK         = "\033[30m"
// DARKRED       = "\033[31m"
// DARKGREEN     = "\033[32m"
// DARKYELLOW    = "\033[33m"
// DARKBLUE      = "\033[34m"
// DARKMAGENTA   = "\033[35m"
// DARKCYAN      = "\033[36m"
// GRAY          = "\033[37m"
// DARKGRAY      = "\033[90m"
// RED           = "\033[91m"
// GREEN         = "\033[92m"
// YELLOW        = "\033[93m"
// BLUE          = "\033[94m"
// MAGENTA       = "\033[95m"
// CYAN          = "\033[96m"
// WHITE         = "\033[97m"


