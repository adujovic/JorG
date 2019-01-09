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

template<size_t N> double parabolic_plus_cos (       void* input      ); 
template<size_t N> double measure            (       void* inX, 
                                                     void* inY        );
template<size_t N> void step                 ( const gsl_rng* r,
                                                     void* inX,
                                                     double step_size );
template<size_t N> void print_state          (       void *xp         );

int main(int argc, char** argv){
    if(argc < 4){
        std::cerr<<"No files given!"<<std::endl;
        exit(-1);
    }

    constexpr size_t SITESNUMBER = 64;

    std::size_t pos;
    std::size_t buff;
    std::string line;
    
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

    std::vector<std::pair<typename ising::IsingModel<SITESNUMBER>::VectorType,double>> supercell; 
    std::ifstream crystal(argv[2]);
    if(crystal.is_open()) {
        while (getline(crystal,line)){
            std::stoi(line,&pos);
            buff = pos;
            std::array<double,3> position;
            for(int i=0; i<3; ++i){
              position[i] = std::stod(line.substr(buff),&pos);
              buff += pos;
            }
            auto moment = std::stod(line.substr(buff));
            supercell.push_back(std::make_pair(position,moment));
        }
        crystal.close();
    }

    typename ising::IsingModel<SITESNUMBER>::VectorType reference;
    std::ifstream ref(argv[1]);
    if(ref.is_open()) {
        getline(ref,line);
        buff = 0;
        for(int i=0; i<3; ++i){
          reference[i] = std::stod(line.substr(buff),&pos);
          buff += pos;
        }
        ref.close();
    }

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
    model.set_supercell(supercell);
    model.set_reference(reference);

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
#ifdef _VERBOSE
    asa.set_print   ( print_state<M>        );
#endif

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

    gauss = std::normal_distribution<>{0.1,0.2};

    std::cout<<"Starting from:         \t";
    for(size_t i=0U; i<M; ++i){
        auto tmp  = gauss(generator);
        std::cout<<tmp<<"\t";
        carr[i]   = tmp; 
        stdarr[i] = tmp;
        stdvec[i] = tmp;
    }
    std::cout<<std::endl;

    asa.run   (&carr,  M*sizeof(double), 100000);
    asa.run<M>(stdarr,                   100000);
    asa.run   (stdvec,                   100000);

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
        x[i] += step_size*2.0*(randomOne-0.5);
    }
}

template<size_t N>
void print_state(void* inX) {
    double* x = static_cast<double*>(inX);
    for(unsigned i = 0U; i<N; ++i)
        std::cout<<" "<<x[i];
}
