#ifndef _ASA_H
#define _ASA_H

#include <iostream>
#include <memory>
#include <type_traits>
#include <vector>
#include <array>

#include <gsl/gsl_siman.h>

namespace gsl{

class SimulatedAnnealing{
public:
    SimulatedAnnealing();
    virtual ~SimulatedAnnealing();

protected:
    const gsl_rng_type * workspace;
    gsl_rng * randomNumberGenerator;

    double (*energy) (void*);
    double (*measure)(void*,
                      void*);
    void   (*step)   (const gsl_rng*,
                            void*,
                            double);
    void   (*print)  (void*);

    gsl_siman_params_t asaParameters;

public:
    void run(void* init, 
             size_t byteSize, 
             size_t _n_tries = 1000);
    void run(std::vector<double>& init,
             size_t _n_tries = 1000);

    template<size_t N>
    void run(std::array<double,N>& init,
             size_t _n_tries = 1000){

        asaParameters.n_tries = _n_tries;
        void* data = static_cast<void*>(init.data());
        gsl_siman_solve(randomNumberGenerator, data,
                        energy, step,
                        measure, print,
                        NULL, NULL, NULL,
                        sizeof(double)*N, asaParameters);
    }

    template<class State>
    void run(State& init, unsigned _size){
        void* data = static_cast<void*>(&init);
        gsl_siman_solve(randomNumberGenerator, data,
                        energy, step,
                        measure, print,
                        NULL, NULL, NULL,
                        _size, asaParameters);
    }

    void set_energy (double (*)(void*));
    void set_measure(double (*)(void*,
                                void*));
    void set_step   (void   (*)(const gsl_rng*,
                                      void*,
                                      double));
    void set_print  (void   (*)(void*));

}; // end of class SimulatedAnnealing

template<typename T>
T myself(T t){
    return t;
}

template<typename R,typename ...T>
R pass(T... args __attribute__((unused))){
	if constexpr(std::is_void<R>::value) return;	
	if constexpr(std::is_arithmetic<R>::value){
		if constexpr(sizeof...(args)==2) return R(1);
		else return R(0);
	}
	return R();
}

} // end of namespace gsl

#endif
