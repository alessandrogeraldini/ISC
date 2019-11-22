#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>

/* calculates detereminant, trace, eigenvalues and eigenvectors (the latter only if the eigenvalues are real)
   if eigenvalues are complex (or imaginary) the output eigenvectors are zero 
   if eigenvalues are complex they are output in two columns, the first is the mod and the second is the angle 
   the two eigenvalues are mod*exp(plum/minus i*angle) */
void linalg2x2(double **input, double evecs[2][2], double evals[2], double *determinant, double *trace) {
	double detval, trval;
	double angle, absq, unnorm_ycomp_evect, normevect;
	//evals = malloc(2*sizeof(double)); evals[0] = malloc(2*sizeof(double)); evals[1] = malloc(2*sizeof(double));
	//evecs = calloc(2,sizeof(double)); evecs[0] = calloc(2,sizeof(double)); evecs[1] = calloc(2,sizeof(double));
	detval = input[0][0]*input[1][1] - input[1][0]*input[0][1];
	trval 	     = input[0][0] + input[1][1];
	if (detval > 0.25*trval*trval) {
		absq = sqrt(detval - 0.25*trval*trval);
		angle = atan2(absq, 0.5*trval);
		evals[0] = detval;
		evals[1]= angle;
	} 
	else {
		absq = sqrt(0.25*trval*trval - detval);
		evals[0] = 0.5*trval + absq;
		evals[1] = 0.5*trval - absq;
		unnorm_ycomp_evect = (evals[0] - input[0][0])/input[0][1];
		normevect = sqrt(1.0 + pow(unnorm_ycomp_evect, 2.0));
		evecs[0][0] = 1.0/normevect;
		evecs[1][0] = unnorm_ycomp_evect/normevect;
		unnorm_ycomp_evect = (evals[1] - input[0][0])/input[0][1];
		normevect = sqrt(1.0 + pow(unnorm_ycomp_evect, 2.0));
		evecs[0][1] = 1.0/normevect;
		evecs[1][1] = unnorm_ycomp_evect/normevect;
	}
	*determinant = detval;
	*trace = trval;
}
 
int test_linalg() {
	double **matrix;
	double evec[2][2], eval[2], det, tr;
	matrix = malloc(2*sizeof(double*));
	matrix[0] = malloc(2*sizeof(double)); matrix[1] = malloc(2*sizeof(double));
	matrix[0][0] = 1.0; matrix[1][0] = 1.0; matrix[0][1] = 1.0; matrix[1][1] = 1.0;
	linalg2x2(matrix, evec, eval, &det, &tr);
	printf("determinant = %f\n", det);
	printf("trace = %f\n", tr);
	printf("evecs = (%f,%f) and (%f, %f)\n", evec[0][0], evec[0][1], evec[1][0], evec[1][1]);
	printf("lambdas = (%f,%f)\n", eval[0], eval[1]);
	return 0;
}

double **set_identity() {
	double **matrix_identity;
	matrix_identity = malloc(2*sizeof(double));
	matrix_identity[0] = malloc(2*sizeof(double));
	matrix_identity[1] = malloc(2*sizeof(double));
	matrix_identity[0][0] = 1.0; matrix_identity[1][1] = 1.0;
	matrix_identity[1][0] = 0.0; matrix_identity[0][1] = 0.0;
	return matrix_identity;
}
