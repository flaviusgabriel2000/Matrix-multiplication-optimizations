/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"
#include "cblas.h"

/* 
 * Add your BLAS implementation here
 */
double* my_solver(int N, double *A, double *B) {
	double alpha = 1.0, beta = 0.0;
	double *C = (double *)malloc(N * N * sizeof(double));

	// C matrix looks like:
	// C = A * B * A_t + B_t * B_t;


	// We firstly overwrite C with the value of B_t * B_t
	// C = B_t * B_t
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, N, N, N, alpha, B, N, B, N, beta, C, N);

	// Now we compute A * B (needed for A * B * A_t) and since we don't need the original
	// value of B anymore (we stored it in C), we can overwrite it with the result
    // B = A * B
    cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, N, N, alpha, A, N, B, N);

	// We compute B * A_t (which expands to A * B * A_t) and add the result to the existing
	// values of C (we specify this using a value of 1.0 for the beta parameter)
	// C = B * A_t + C
	beta = alpha;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, alpha, B, N, A, N, beta, C, N);

	return C;
}
