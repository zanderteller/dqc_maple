#ifndef UTILS_H_
#define UTILS_H_

#ifdef _WIN32
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
#endif

#include "matrix.h"

Matrix<double> L2NormRows(Matrix<double> mat);

void EigValsVecs(Matrix<double>& mat, Matrix<double>& eig_vals, Matrix<double>& eig_vecs);

#endif //UTILS_H_
