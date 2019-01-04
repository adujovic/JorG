#include <iostream>
#include <iomanip>

#include <vector>
#include <array>
#include <numeric>
#include <random>

#include <cmath>
#include <cstring>

#include "asa.h"
#include "ising.h"

template<size_t N>
double parabolic_plus_cos(void *input) {
  double* x = static_cast<double *>(input);
  double result = 0.0;
  double alpha = 0.5;

  for(unsigned i = 0U; i<N; ++i)
    result += 1*pow(x[i],2) - alpha*cos(10.0*x[i]) + alpha;
  return result;
}

template<size_t N>
double measure(void *inX, void *inY) {
  double* x = static_cast<double *>(inX);
  double* y = static_cast<double *>(inY);

  double err = 0.0;
  for(unsigned i = 0U; i<N; ++i)
    err += pow(x[i]-y[i],2);

  return sqrt(err)/N;
}

template<size_t N>
void step(const gsl_rng* r, void* inX, double step_size) {
  double* x = static_cast<double*>(inX);

  for(unsigned i = 0U; i<N; ++i){
  double randomOne = gsl_rng_uniform(r);
  double randomTwo = gsl_rng_uniform(r);
     x[i] += sqrt(-2.0*log(randomOne))*cos(M_2_PI*randomTwo)*step_size;
  }
}

void P1(void *xp) {
	std::cout<<"\t"<<*((double *) xp);
}

//int main (int argc, char **argv){
int main (){
	auxiliary::IsingModel<100> model;
	model.add_interaction(1,1,0.1);
	model.add_interaction(std::make_tuple(1U,1U,1.0));
	std::vector<std::tuple<unsigned,unsigned,double>> d;
	d.push_back(std::make_tuple(1U,1U,1.0));
	d.push_back(std::make_tuple(0U,0U,1.0));
	model.add_interaction(d);
	exit(0);
	constexpr size_t M = 3;
	std::cout << std::setprecision(9);
	gsl::SimulatedAnnealing asa;
	asa.set_energy(parabolic_plus_cos<M>);
	asa.set_measure(measure<M>);
	asa.set_step(step<M>);
	//asa.set_print(P1);

	std::random_device randomDevice{};
	std::mt19937 generator{randomDevice()};
    	std::normal_distribution<> gauss{0.2,1};

	std::vector<double>  stdvec(M);
	double              carr[M];
	std::array<double,M> stdarr;
	for(size_t i=0U; i<M; ++i){
		carr[i]   = gauss(generator);
		stdarr[i] = gauss(generator);
		stdvec[i] = gauss(generator);
	}
	asa.run(&carr,M*sizeof(double),1000);
	asa.run(stdarr,1000);
	asa.run(stdvec,1000);
	std::cout<<"Cstyle array solution: ";
	for(size_t i=0U; i<M; ++i) std::cout<<" "<<carr[i];
	std::cout<<std::endl;
	std::cout<<"  std::array solution: ";
	for(size_t i=0U; i<M; ++i) std::cout<<" "<<stdarr[i];
	std::cout<<std::endl;
	std::cout<<" std::vector solution: ";
	for(size_t i=0U; i<M; ++i) std::cout<<" "<<stdvec[i];
	std::cout<<std::endl;

	return 0;
}
