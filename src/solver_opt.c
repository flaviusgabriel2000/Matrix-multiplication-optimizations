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
 * Add your optimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	double *C = (double *)calloc(N * N, sizeof(double));
	register double *A_t = transpose(A, N, 1);
	register double *B_t = transpose(B, N, 0);
	register int i, j, k;

	// Firstly, we compute B_t * B_t and store it in C
	for (i = 0; i < N; i++) {
		register double *orig_pb_t = &B_t[i * N]; // i * N + 0
        for (j = 0; j < N; j++) {
			register double *pb_t_left = orig_pb_t;
			register double *pb_t_right = &B_t[j]; // 0 * N + j
			register double sum = 0;
            for (k = 0; k < N; k++) {
				sum += *pb_t_left * *pb_t_right;
				pb_t_left++;
				pb_t_right += N;
            }
			C[i * N + j] = sum;
        }
    }

	// Now, we compute A * B and store it in B_t since we don't need it anymore
	// Keep in mind that A is an upper triangular matrix
	memset(B_t, 0, N * N * sizeof(*B_t));
	for (i = 0; i < N; i++) {
		register double *orig_pa = &A[i * N + i]; // ignore zeros, start at the diagonal
		for (j = 0; j < N; j++) {
			register double *pa = orig_pa;
			register double *pb = &B[i * N + j];
			register double sum = 0;
			for (k = i; k < N; k++) {
				sum += *pa * *pb;
				pa++;
				pb += N;
			}
			B_t[i * N + j] = sum;
		}
	}

	// We now have to compute A * B * A_t which is the overwritten B_t * A_t.
	// We can store the result in either A or B since we don't need them anymore
	// Let's assume we store it in A
	memset(A, 0, N * N * sizeof(*A));
	for (i = 0; i < N; i++) {
		register double *orig_pb_t = &B_t[i * N]; // i * N + 0
		for (j = 0; j < N; j++) {
			register double *pb_t = orig_pb_t;
			register double *pa_t = &A_t[j]; // 0 * N + j
			register double sum = 0;
			for (k = 0; k < N; k++) {
				sum += *pb_t * *pa_t;
				pb_t++;
				pa_t += N;
			}
			A[i * N + j] = sum;
		}
	}

	// Finally, we have to compute C = A * B * A_t + B_t * B_t, which
	// is equivalent to C += A in our case
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			C[i * N + j] += A[i * N + j];
		}
	}

	free(A_t);
	A_t = NULL;
	free(B_t);
	B_t = NULL;

	return C;
}
