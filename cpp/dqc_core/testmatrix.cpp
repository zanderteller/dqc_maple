#include "matrix.h"
#include <complex>  // we were using our own homegrown "complex.h"
#include <cstdio>

using std::complex;

int main(int argc, char **argv)
{
	Matrix<complex<double>> foo(3,3);
	Matrix<double> bar(3,3);
	foo.set(0,0, complex<double>(1, 0));
	foo.set(1,0, complex<double>(2, 1));
	foo.set(2,0, complex<double>(5, 2));
	foo.set(0,1, complex<double>(3, 0.5));
	foo.set(1,1, complex<double>(4, -1));
	foo.set(2,1, complex<double>(3.97854, .00111234));
	foo.set(0,2, complex<double>(6, 1));
	foo.set(1,2, complex<double>(9.238478234, 0.23491230));
	foo.set(2,2, complex<double>(-19.34832, 1.32345345));
	bar.set(0,0, 5);
	bar.set(1,0, 6);
	bar.set(2,0, 9);
	bar.set(0,1, 7);
	bar.set(1,1, 8);
	bar.set(2,1, 10.2341234);
	bar.set(0,2, 1);
	bar.set(1,2, 3.14159);
	bar.set(2,2, 8.12349);
	auto baz(foo * bar);
	for (int y = 0; y < baz.rows(); y++)
	{
		printf("[");
		for (int x = 0; x < baz.cols(); x++)
		{
			printf(x ? ", (%f + %fi)" : "(%f + %fi)", baz.get(x,y).real(), baz.get(x,y).imag());
		}
		printf("]\n");
	}
	return 0;
}