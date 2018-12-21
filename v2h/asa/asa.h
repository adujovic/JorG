#ifndef _ASA_H
#define _ASA_H

#include <iostream>
#include <memory>
#include <type_traits>

#include <gsl/gsl_siman.h>

namespace gsl{

class SimulatedAnnealing{
public:
	SimulatedAnnealing();
	~SimulatedAnnealing();

protected:
	size_t systemDimension;

	const gsl_rng_type * T;
	gsl_rng * r;

	double (*energy) (void*);
	double (*measure)(void*,void*);
	void   (*step)   (const gsl_rng* r,void*,double);
	void   (*print)  (void*);

	gsl_siman_params_t asaParameters;

public:
	
	void set_energy (double (*)(void*));
	void set_measure(double (*)(void*,void*));
	void set_step   (void   (*)(const gsl_rng* r,void*,double));
	void set_print  (void   (*)(void*));







};

template<typename R,typename ...T>
R pass(T...);
} // end of namespace gsl


#endif
