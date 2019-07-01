#ifndef _PI_MANDELBROT
#define _PI_MANDELBROT
#include <complex>
#include <cstdio>
 
typedef std::complex<double> complex;
int MandelbrotCalculate(const complex& c, int maxiter);

extern "C"
void mandelbrot(int width = 78, int height = 44);
#endif
