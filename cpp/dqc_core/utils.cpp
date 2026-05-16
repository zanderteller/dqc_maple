/********************

General utility functions.

********************/

#include <assert.h>
#include "utils.h"
#include "matrix.h"
#include <lapacke.h>
#include <cmath>

// given a matrix mat, return matrix mat_norm with each row of mat L2-normalized
Matrix<double> L2NormRows(Matrix<double> mat)
{
	int nrows = mat.rows();
	int ncols = mat.cols();
	Matrix<double> mat_norm(ncols, nrows);

	#pragma omp parallel for shared(mat, nrows, ncols, mat_norm)
	for (int row_idx = 0; row_idx < nrows; row_idx++) {
		double val;
		double norm = 0.0;
		for (int col_idx = 0; col_idx < ncols; col_idx++) {
			val = mat.get(col_idx, row_idx);
			norm += val * val;
		}
		norm = sqrt(norm);
		for (int col_idx = 0; col_idx < ncols; col_idx++) {
			mat_norm.set(col_idx, row_idx, mat.get(col_idx, row_idx) / norm);
		}
	} // end for each row
	return mat_norm;
} // end function L2NormRows

// eig_vals will contain the eigenvalues, in ascending order
// eig_vecs will contain the associated eigenvectors as columns of the matrix
// note: if calculations fail to converge, eig_vals and eig_vecs will contain all zeros
// after this function returns.
void EigValsVecs(Matrix<double>& mat, Matrix<double>& eig_vals, Matrix<double>& eig_vecs)
{
	// verify correct size of matrices
	lapack_int n = mat.rows();
	assert(n == mat.cols() && n == eig_vals.rows() && n == eig_vecs.rows() && n == eig_vecs.cols());

	// we ignore these arguments
	lapack_int* m = (lapack_int*)malloc(sizeof(lapack_int));
	lapack_int* isuppz = (lapack_int*)malloc(2 * n * sizeof(lapack_int));

	double tolerance = 1e-12;
	lapack_int info = 1;
	while (info != 0) {
		// LAPACKE_dsyevr overwrites parts of the input matrix, so make a copy
		Matrix<double> mat_copy(mat);

		// following advice in the MKL documentation, we use dsyevr instead of dsyev -- supposedly it's faster and uses less memory.
		// dsyevr has some extra arguments that seem to be ignorable for our purposes (see above).
		info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'A', 'U', n, mat_copy.raw(), n, 0, 0, 0, 0, tolerance, m, eig_vals.raw(), eig_vecs.raw(), n, isuppz);

		if (info != 0) {
			// increase tolerance (up to a limit) until calculations succeed
			tolerance *= 2;
			if (tolerance > 1e-10)
				break;
		} // end if increasing tolerance
	} // end while increasing tolerance until calculations succeed

	free(m);
	free(isuppz);

	// OLD:
	//assert(info == 0); // assert success of the eigenvector calculations
	// NEW:
	if (info != 0) {
		// eigenvector calculations failed, probably due to too many zero entries in mat.
		// return all zeros in eig_vals and eig_vecs.
		// (note: this is better than raising an assertion error, which Python/Maple code calling C++ code
		// can't catch. calling code needs to check for and deal with the all-zero case.)
		
		// put zeros into eig_vals
		for (int row_idx = 0; row_idx < n; row_idx++)
			eig_vals.set(0, row_idx, 0);

		// put zeros into eig_vecs
		for (int col_idx = 0; col_idx < n; col_idx++)
			for (int row_idx = 0; row_idx < n; row_idx++)
				eig_vecs.set(col_idx, row_idx, 0);
	} // end if calculcations failed to converge

} // end function EigValsVecs

extern "C"
{
	DLLEXPORT void SVDC(double* mat, int num_rows, int num_cols, double* u, double* s, double* vt)
	{
		/*
		NOTES
		- use DGESDD rather than DGESVD (divide-and-conquer algorithm is faster in most cases.  there may
		  be corner cases where one version or the other yields more accurate results, but we're punting
		  on there here.)
		- from the documentation: "for an <m x n> input matrix and job type 'S', the first min(m, n) columns
		  of U and the first min(m, n) rows of VT are returned in the arrays u and vt."
		- so, u should be <m x min(m, n)>, s should be <min(m, n)> (1-D vector), and vt should be <min(m, n) x n>.
		- singular values in s are returned in descending order.
		- the contents of the input matrix are destroyed, so pass in a copy.
		*/
		double* mat_copy = new double[num_rows * num_cols];
		memcpy(mat_copy, mat, num_rows * num_cols * sizeof(double));
		lapack_int ldu = lapack_int(num_rows < num_cols ? num_rows : num_cols);
		lapack_int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', lapack_int(num_rows), lapack_int(num_cols), mat_copy,
			lapack_int(num_cols), s, u, ldu, vt, lapack_int(num_cols));
		delete[] mat_copy;
		assert(info == 0); // assert success of the SVD calculations
	} // end function SVDC

	// this version takes a simple array instead of a Matrix object -- for external use (from Maple, etc.)
	DLLEXPORT void L2NormRowsC(double* mat, int num_rows, int num_cols, double* mat_out)
	{
		#pragma omp parallel for shared(mat, num_rows, num_cols, mat_out)
		for (int row_idx = 0; row_idx < num_rows; row_idx++) {
			double val;
			double norm = 0.0;
			int col_idx;
			int base_idx = row_idx * num_cols;
			int arr_idx = base_idx;
			for (col_idx = 0; col_idx < num_cols; col_idx++) {
				val = mat[arr_idx++];
				norm += val * val;
			}
			norm = sqrt(norm);
			arr_idx = base_idx;
			for (col_idx = 0; col_idx < num_cols; col_idx++) {
				mat_out[arr_idx] = mat[arr_idx] / norm;
				arr_idx++;
			}
		} // end for (each row)
	} // end function L2NormRowsC

} // end extern C
