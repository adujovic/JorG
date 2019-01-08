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

template<size_t N> double parabolic_plus_cos (       void* input      ); 
template<size_t N> double measure            (       void* inX, 
                                                     void* inY        );
template<size_t N> void step                 ( const gsl_rng* r,
                                                     void* inX,
                                                     double step_size );
template<size_t N> void print_state          (       void *xp         );

int main (){
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

    constexpr size_t SITESNUMBER = 64;
    ising::IsingModel<SITESNUMBER> model;

    std::vector<std::tuple<unsigned,unsigned,double>> d;
    for(size_t i = 0U; i<SITESNUMBER-1; ++i){
        d.push_back(std::make_tuple(i,i+1,-fabs(gauss(generator))));
    }
    model.add_interaction(d);

    model.run();

    constexpr size_t M = 3;
    std::cout << std::setprecision(9);
    gsl::SimulatedAnnealing asa;
    asa.set_energy  ( parabolic_plus_cos<M> );
    asa.set_measure ( measure<M>            );
    asa.set_step    ( step<M>               );
//    asa.set_print   ( print_state<M>        );

    gauss = std::normal_distribution<>{0.2,1};

    std::vector<double>  stdvec(M);
    double               carr[M];
    std::array<double,M> stdarr;

    std::cout<<"########################################"<<std::endl;
    std::cout<<"##                                    ##"<<std::endl;
    std::cout<<"##                                    ##"<<std::endl;
    std::cout<<"##        Simulated Annealing         ##"<<std::endl;
    std::cout<<"##                                    ##"<<std::endl;
    std::cout<<"##                                    ##"<<std::endl;
    std::cout<<"########################################"<<std::endl;

    std::cout<<"Starting from:         \t";
    for(size_t i=0U; i<M; ++i){
        auto tmp  = gauss(generator);
        std::cout<<tmp<<"\t";
        carr[i]   = tmp; 
        stdarr[i] = tmp;
        stdvec[i] = tmp;
    }
    std::cout<<std::endl;

    asa.run(&carr,  M*sizeof(double), 1000);
    asa.run(stdarr,                   1000);
    asa.run(stdvec,                   1000);

    std::cout<<"Cstyle array solution:";
    for(size_t i=0U; i<M; ++i) std::cout<<"\t"<<carr[i];
    std::cout<<std::endl;

    std::cout<<"  std::array solution:";
    for(size_t i=0U; i<M; ++i) std::cout<<"\t"<<stdarr[i];
    std::cout<<std::endl;

    std::cout<<" std::vector solution:";
    for(size_t i=0U; i<M; ++i) std::cout<<"\t"<<stdvec[i];
    std::cout<<std::endl;

    return 0;
}

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
double measure(void* inX, void* inY) {
    double* x = static_cast<double*>(inX);
    double* y = static_cast<double*>(inY);

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

template<size_t N>
void print_state(void* inX) {
    double* x = static_cast<double*>(inX);
    for(unsigned i = 0U; i<N; ++i)
        std::cout<<" "<<x[i];
}
