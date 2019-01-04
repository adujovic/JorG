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
	double (*measure)(void*,void*);
	void   (*step)   (const gsl_rng*,void*,double);
	void   (*print)  (void*);

	gsl_siman_params_t asaParameters;

public:
	void run(void* init, size_t byteSize, size_t _n_tries = 1000);
	void run(std::vector<double>& init,   size_t _n_tries = 1000);

	template<size_t N>
	void run(std::array<double,N>& init,  size_t _n_tries = 1000){
		asaParameters.n_tries = _n_tries;
		void* data = static_cast<void*>(init.data());
		gsl_siman_solve(randomNumberGenerator, data,
				energy, step,
				measure, print,
				NULL, NULL, NULL,
				sizeof(double)*N, asaParameters);
	}
	
	void set_energy (double (*)(void*));
	void set_measure(double (*)(void*,void*));
	void set_step   (void   (*)(const gsl_rng*,void*,double));
	void set_print  (void   (*)(void*));







};

template<typename R,typename ...T>
R pass(T...);


} // end of namespace gsl


#endif
