#include "asa.h"

gsl::SimulatedAnnealing::SimulatedAnnealing(){
	gsl_rng_env_setup();
	
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	energy  = gsl::pass<double,void*>;
	measure = gsl::pass<double,void*,void*>;
	step    = gsl::pass<void,const gsl_rng*,void*,double>;
	print   = gsl::pass<void,void*>;

	asaParameters.n_tries       = 100;
	asaParameters.iters_fixed_T = 5;
	asaParameters.step_size     = 1e-2;
	asaParameters.k             = 1.0;
	asaParameters.t_initial     = 1e-2;
	asaParameters.mu_t          = 1.0 + 3e-4;
	asaParameters.t_min         = 1e-6;

	std::cout<<"energy"<<"  = "<<energy(NULL)<<std::endl;
	std::cout<<"measure"<<" = "<<measure(NULL,NULL)<<std::endl;
}

gsl::SimulatedAnnealing::~SimulatedAnnealing(){
	gsl_rng_free(r);
}

void gsl::SimulatedAnnealing::set_energy (double (*_new)(void*)){
	energy = _new;
}

void gsl::SimulatedAnnealing::set_measure(double (*_new)(void*,void*)){
	measure = _new;
}

void gsl::SimulatedAnnealing::set_step   (void   (*_new)(const gsl_rng* r,void*,double)){
	step = _new;
}

void gsl::SimulatedAnnealing::set_print  (void   (*_new)(void*)){
	print = _new;
}

template<typename R,typename ...T>
R gsl::pass(T... args __attribute__((unused))){
	if constexpr(std::is_void<R>::value) return;	
	if constexpr(std::is_arithmetic<R>::value){
		if constexpr(sizeof...(args)==2) return R(1);
		else return R(0);
	}
	return R();
}

