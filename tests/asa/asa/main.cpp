#include <iostream>
#include <iomanip>

#include <vector>
#include <array>
#include <numeric>
#include <random>

#include <cmath>
#include <cstring>

#include "../../../asa/asa.h"
#include "../../../asa/ising.h"

template<size_t N> double parabolic_plus_cos(void* input);
template<size_t N> double measure(void* inX, void* inY);
template<size_t N> void   step(const gsl_rng* r, void* inX, double step_size);
template<size_t N> void   print_state(void *xp);
template<size_t M, typename ArrayLike> 
void print_result(const ArrayLike& arr);

int main (){
    constexpr size_t M = 3;
    std::random_device randomDevice{};
    std::mt19937 generator{randomDevice()};
    std::normal_distribution<> gauss(0.1,0.2);

    std::cout << std::setprecision(9);
    gsl::SimulatedAnnealing asa;
    asa.set_energy  ( parabolic_plus_cos<M> );
    asa.set_measure ( measure<M>            );
    asa.set_step    ( step<M>               );
#ifdef _VERBOSE
    asa.set_print   ( print_state<M>        );
#endif

    std::vector<double>  stdvec(M);
    double               carr[M];
    std::array<double,M> stdarr;

    std::cout<<"########################################"<<std::endl;
    std::cout<<"##        Simulated Annealing         ##"<<std::endl;
    std::cout<<"########################################"<<std::endl;

    typename gsl::SimulatedAnnealing::Parameters params;

    params.n_tries       = 100;      /* how many points
                                      * to try for each step
                                      */
    params.iters_fixed_T = 1000;     // how many iterations at each temperature?
    params.step_size     = 1e-3;     /* max step size in the random walk;
                                      * !!! unused for ising::IsingModel !!!
                                      */
    params.k             = 1.0;      /* Boltzman constant;
                                      * scaling factor for "energy";
                                      * should be of the order of "energy"
                                      */
    params.t_initial     = 1e4;      // initial "temperature"
    params.mu_t          = 1.0+1e-1; // "temperature" step
    params.t_min         = 1e-5;     // final "temperature"
    asa.set_parameters(params);

    std::cout<<"Starting from:\t";
    for(size_t i=0U; i<M; ++i){
        auto tmp  = gauss(generator);
        carr[i]   = tmp;
        stdarr[i] = tmp;
        stdvec[i] = tmp;
    }
    print_result<M>(carr);

    asa.run   (&carr,  M*sizeof(double), 100000);
    asa.run<M>(stdarr,                   100000);
    asa.run   (stdvec,                   100000);

    std::cout<<"Cstyle array solution:";
    print_result<M>(carr);

    std::cout<<"  std::array solution:";
    print_result<M>(stdarr);

    std::cout<<" std::vector solution:";
    print_result<M>(stdvec);

    return 0;
}

template<size_t M, typename ArrayLike>
void print_result(const ArrayLike& arr){
    for(size_t i=0U; i<M; ++i) std::cout<<"\t"<<arr[i];
    std::cout<<std::endl;
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
        x[i] += step_size*2.0*(randomOne-0.5);
    }
}

template<size_t N>
void print_state(void* inX) {
    double* x = static_cast<double*>(inX);
    for(unsigned i = 0U; i<N; ++i)
        std::cout<<" "<<x[i];
}
