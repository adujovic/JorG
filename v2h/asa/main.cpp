#include <iostream>
#include <iomanip>

#include <cmath>
#include <cstring>

#include "asa.h"

double E1(void *xp) {
  double x = * ((double *) xp);

  return exp(-pow((x-1.0),2.0))*sin(8*x);
}

double M1(void *xp, void *yp) {
  double x = *((double *) xp);
  double y = *((double *) yp);

  return fabs(x - y);
}

void S1(const gsl_rng * r, void *xp, double step_size) {
  double old_x = *((double *) xp);
  double new_x;

  double u = gsl_rng_uniform(r);
  new_x = u * 2 * step_size - step_size + old_x;

  std::memcpy(xp, &new_x, sizeof(new_x));
}

void P1(void *xp) {
	std::cout<<"\t"<<*((double *) xp);
}

int main (int argc, char **argv){
	std::cout << std::setprecision(9);
	gsl::SimulatedAnnealing asa;
	asa.set_energy(E1);
	asa.set_measure(M1);
	asa.set_step(S1);
	asa.set_print(P1);

	double init[1] = {1.0};
	asa.run(&init,sizeof(double));

	return 0;
}
