#ifndef _ISING_H
#define _ISING_H

#include <iostream>

#include <tuple>
#include <array>
#include <vector>
#include <bitset>

#include <algorithm>
#include <random>
#include <limits>
#include <string>
#include <type_traits>

#include "asa.h"
#include "arithmeticvector.h"

namespace auxiliary{

template<size_t N>
class IsingModel{
public:
	IsingModel(bool QUICK_START=true):randomDevice(),randomEngine(randomDevice()),flip(0U,1U),
	                                  uniform(
 					    std::numeric_limits<unsigned long long int>::min(),
                                            std::numeric_limits<unsigned long long int>::max()),
					  lattice(randomEngine,uniform){
          if(QUICK_START) hamiltonian.reserve(N*N/2); // maximum number of interactions ( each site with each )
          else            hamiltonian.reserve(N);     // minimum number of interactions ( once  per each site )
	  randomize(12);

	}
	virtual ~IsingModel(){}

	typedef std::tuple<unsigned,unsigned,double>      TwoSiteInteraction;
	typedef std::vector<TwoSiteInteraction>       HamiltonianType;
	typedef celerium::ArithmeticVector                VectorType;
private:
	std::random_device randomDevice;
        std::mt19937 randomEngine;
        std::uniform_int_distribution<unsigned> flip;
        std::uniform_int_distribution<unsigned long long int>
		uniform;

	struct LatticeType{
	  std::bitset<N>           nodes;
	  std::array<VectorType,N> positions;

	  LatticeType(std::mt19937& engine,
		      std::uniform_int_distribution<unsigned long long int>& distribution)
	              :nodes(distribution(engine)){
	    constexpr auto seedSize    = 8*sizeof(unsigned long long int);
	    auto currentSize = seedSize;
	    while (currentSize < N){
	      nodes <<= seedSize;
	      nodes  |= std::bitset<N>(distribution(engine));
	      currentSize += seedSize;
	    }
	    std::cout<<nodes<<std::endl;
	  }
	} lattice;

	gsl::SimulatedAnnealing solver;
	HamiltonianType         hamiltonian;

public:
	unsigned randomize(size_t maxNumberOfFlips=1){
	  unsigned ones  = std::min(maxNumberOfFlips,N);
	  unsigned zeros = N - ones;
	  std::string mask(ones,'1');
	  mask += std::string(zeros,'0');
	  std::shuffle(mask.begin(),mask.end(),randomEngine);
	  std::cout<<mask<<std::endl;

	  lattice.nodes ^= std::bitset<N>(mask);
          std::cout<<lattice.nodes<<std::endl;
	  return 0;
  	}	  

public:
	void add_interaction(...){
	  std::cerr<<"Wrong input for auxiliary::IsingModel::add_interaction:"<<std::endl;
	  std::cerr<<"           Either non-iterable or iterable of non <int,int,float> tuples."<<std::endl;
	  std::cerr<<"           Nothing was added!"<<std::endl;
	}

	// For SFINAE compiler-time evaulation
	template<class T>
	T tester(T t){
	  if(std::is_integral<T>::value) return static_cast<unsigned>(t);
	  return t;
	};

	template<class intlike, class floatlike>
	auto add_interaction(intlike i, intlike j, floatlike J) -> decltype((unsigned)(tester<intlike>)(i),void()){
	  hamiltonian.push_back(std::make_tuple(i,j,J));
	}

	template<class T>
	auto add_interaction(T interaction) -> decltype((TwoSiteInteraction&)(tester<T>)(interaction),void()){
	  hamiltonian.push_back(interaction);
	}

	template<class T>
	auto add_interaction(T interaction) -> decltype((TwoSiteInteraction&&)(tester<T>)(interaction),void()){
	  hamiltonian.push_back(interaction);
	}

	template<class Iterable>
	auto add_interaction(const Iterable& interactions) -> decltype((Iterable::iterator)(interactions.begin())(),void()){
	    for(auto& interaction : interactions)
	      hamiltonian.push_back(interaction);
	}

};	

}//end of namespace auxiliary

#endif
