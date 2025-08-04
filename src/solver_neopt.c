/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"
#include <string.h>

/*
 * Returns the transpose of a matrix
 */
double *transpose(double *mat, int N, int isTriangular) {
	if (!mat) {
		return NULL;
	}

	double *mat_t;

	if (isTriangular) {
		mat_t = (double *)calloc(N * N, sizeof(double));
		for (int i = 0; i < N; i++) {
			for (int j = i; j < N; j++) {
				mat_t[j * N + i] = mat[i * N + j];
			}
		}
	} else {
		mat_t = (double *)malloc(N * N * sizeof(double));
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				mat_t[j * N + i] = mat[i * N + j];
			}
		}
	}

	return mat_t;
}


/*
 * Add your unoptimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	double *C = (double *)calloc(N * N, sizeof(double));
	double *A_t = transpose(A, N, 1);
	double *B_t = transpose(B, N, 0);

	// Firstly, we compute B_t * B_t and store it in C
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				C[i * N + j] += B_t[i * N + k] * B_t[k * N + j];
			}
		}
	}


	// Now, we compute A * B and store it in B_t since we don't need it anymore
	// Keep in mind that A is an upper triangular matrix
	memset(B_t, 0, N * N * sizeof(*B_t));
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = i; k < N; k++) {
				B_t[i * N + j] += A[i * N + k] * B[k * N + j];
			}
		}
	}

	// We now have to compute A * B * A_t which is the overwritten B_t * A_t.
	// We can store the result in either A or B since we don't need them anymore
	// Let's assume we store it in A
	memset(A, 0, N * N * sizeof(*A));
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				A[i * N + j] += B_t[i * N + k] * A_t[k * N + j];
			}
		}
	}

	// Finally, we have to compute C = A * B * A_t + B_t * B_t, which
	// is equivalent to C += A in our case
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			C[i * N + j] += A[i * N + j];
		}
	}

	free(A_t);
	A_t = NULL;
	free(B_t);
	B_t = NULL;

	return C;
}
