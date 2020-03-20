// Author: Alessandro Geraldini
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>

/* calculates detereminant, trace, eigenvalues and eigenvectors (the latter only if the eigenvalues are real)
   if eigenvalues are complex (or imaginary) the output eigenvectors are zero 
   if eigenvalues are complex they are output in two columns, the first is the mod and the second is the angle 
   the two eigenvalues are mod*exp(plum/minus i*angle) */
void linalg2x2(double **input, double **evecs, double *evals, double *determinant, double *trace) {
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
 
//int test_linalg() {
//	double **matrix;
//	double evec[2][2], eval[2], det, tr;
//	matrix = malloc(2*sizeof(double*));
//	matrix[0] = malloc(2*sizeof(double)); matrix[1] = malloc(2*sizeof(double));
//	matrix[0][0] = 1.0; matrix[1][0] = 1.0; matrix[0][1] = 1.0; matrix[1][1] = 1.0;
//	linalg2x2(matrix, evec, eval, &det, &tr);
//	printf("determinant = %f\n", det);
//	printf("trace = %f\n", tr);
//	printf("evecs = (%f,%f) and (%f, %f)\n", evec[0][0], evec[0][1], evec[1][0], evec[1][1]);
//	printf("lambdas = (%f,%f)\n", eval[0], eval[1]);
//	return 0;
//}

double **set_identity() {
	double **matrix_identity;
	matrix_identity = malloc(2*sizeof(double));
	matrix_identity[0] = malloc(2*sizeof(double));
	matrix_identity[1] = malloc(2*sizeof(double));
	matrix_identity[0][0] = 1.0; matrix_identity[1][1] = 1.0;
	matrix_identity[1][0] = 0.0; matrix_identity[0][1] = 0.0;
	return matrix_identity;
}

double **set_zeros() {
	double **matrix_zeros;
	matrix_zeros = calloc(2,sizeof(double));
	matrix_zeros[0] = calloc(2,sizeof(double));
	matrix_zeros[1] = calloc(2,sizeof(double));
	//matrix_zeros[0][0] = 0.0; matrix_zeros[1][1] = 0.0;
	//matrix_zeros[1][0] = 0.0; matrix_zeros[0][1] = 0.0;
	//printf("%f %f %f %f\n", matrix_zeros[0][0], matrix_zeros[0][1], matrix_zeros[1][0], matrix_zeros[1][1]);
	return matrix_zeros;
}

double **invert2x2(double **input2x2, double* pdet_tangent) {
	double **result=malloc(2*sizeof(double*));
	int i=0;
	for (i=0;i<2;i++) {
		result[i] = malloc(2*sizeof(double));
	}
	//printf("input[0][0] =%f, input[0][1] =%f, input[1][0] =%f, input[1][1] =%f\n", input2x2[0][0], input2x2[0][1], input2x2[1][0], input2x2[1][1]);
	*pdet_tangent = input2x2[0][0]*input2x2[1][1] - input2x2[0][1]*input2x2[1][0];
	result[0][0] =  input2x2[1][1]/(*pdet_tangent);
	result[1][1] =  input2x2[0][0]/(*pdet_tangent);
	result[0][1] = -input2x2[0][1]/(*pdet_tangent);
	result[1][0] = -input2x2[1][0]/(*pdet_tangent);
	//printf("det =%f\n", *pdet_tangent);
	//printf("result[0][0] =%f, result[0][1] =%f, result[1][0] =%f, result[1][1] =%f\n", result[0][0], result[0][1], result[1][0], result[1][1]);
	return result;
}

void printmat(char *name, double **input, int dimrows, int dimcols) {
	int row, col;
	printf("matrix %s =\n", name);
	for (row=0;row<dimrows;row++) {
		printf(" |");
		for (col=0;col<dimcols;col++) {
			printf(" %9.5f", input[row][col]);
		}
		printf(" |\n");
	}
}

double **multiply2x2(double **input2x2, double **input2xdims, int dims) {
	/* has been checked to work for 2x2 matrix multiplying a 2x1 column vector (and I think for matrix multiplication too) */
	double **result; //=malloc(dims*sizeof(double*)); // this was wrong?
	int i=0;
	result = malloc(2*sizeof(double*));
	for (i=0;i<2;i++) {
		result[i] = malloc(dims*sizeof(double*));
	}
	//printmat("input2x2", input2x2, 2, 2);
	//printmat("input2xdims", input2xdims, 2, dims);
	for (i=0;i<dims;i++) {
		result[0][i] = input2x2[0][0]*input2xdims[0][i] + input2x2[0][1]*input2xdims[1][i]; //input2xdims[1][i];
		result[1][i] = input2x2[1][0]*input2xdims[0][i] + input2x2[1][1]*input2xdims[1][i]; //input2xdims[1][i];
	}
	//printmat("result = input2x2*input2xdims", result, 2, dims);
	return result;
}

void multiply2x2reassign(double **input1, double **input2, int argument_number) {
	/* has been checked to work for 2x2 matrix multiplying a 2x1 column vector (and I think for matrix multiplication too) */
	double a, b, c, d;
	//printmat("input2x2", input2x2, 2, 2);
	//printmat("input2xdims", input2xdims, 2, dims);
	a = input1[0][0]*input2[0][0] + input1[0][1]*input2[1][0]; 
	b = input1[0][0]*input2[0][1] + input1[0][1]*input2[1][1]; 
	c = input1[1][0]*input2[0][0] + input1[1][1]*input2[1][0]; 
	d = input1[1][0]*input2[0][1] + input1[1][1]*input2[1][1]; 
	if (argument_number == 1) {
		input1[0][0] = a;
		input1[0][1] = b;
		input1[1][0] = c;
		input1[1][1] = d;
	}
	else {
		input2[0][0] = a;
		input2[0][1] = b;
		input2[1][0] = c;
		input2[1][1] = d;
	}
}

double *multiply(double **input1, double *input2) {
	/* has been checked to work for 2x2 matrix multiplying a 2x1 column vector (and I think for matrix multiplication too) */
	double a, b, *result=malloc(2*sizeof(double));
	//printmat("input2x2", input2x2, 2, 2);
	//printmat("input2xdims", input2xdims, 2, dims);
	a = input1[0][0]*input2[0] + input1[0][1]*input2[1]; 
	b = input1[1][0]*input2[0] + input1[1][1]*input2[1]; 
	result[0] = a;
	result[1] = b;
	return result;
}

double **add2x2(double c1, double **input2x2, double c2, double **input2xdims, int dims) {
	/* has been checked to work for 2x2 matrix multiplying a 2x1 column vector (and I think for matrix multiplication too) */
	double **result=malloc(dims*sizeof(double*));
	int i=0;
	for (i=0;i<dims;i++) {
		result[i] = malloc(dims*sizeof(double*));
	}
	printmat("input2x2", input2x2, 2, 2);
	printmat("input2xdims", input2xdims, 2, dims);
	for (i=0;i<dims;i++) {
		result[0][i] = c1*input2x2[0][i] + c2*input2xdims[0][i];
		result[1][i] = c1*input2x2[1][i] + c2*input2xdims[1][i];
	}
	printmat("sum =", result, 2, dims);
	return result;
}

void add2x2reassign(double c1, double **input1, double c2, double **input2, int argument_number) {
	/* has been checked to work for 2x2 matrix multiplying a 2x1 column vector (and I think for matrix multiplication too) */
	double a, b, c, d;
	//printmat("input2x2", input2x2, 2, 2);
	//printmat("input2xdims", input2xdims, 2, dims);
	a = c1*input1[0][0] + c2*input2[0][0];
	b = c1*input1[0][1] + c2*input2[0][1];
	c = c1*input1[1][0] + c2*input2[1][0];
	d = c1*input1[1][1] + c2*input2[1][1];
	if (argument_number == 1) {
		input1[0][0] = a;
		input1[0][1] = b;
		input1[1][0] = c;
		input1[1][1] = d;
	}
	else {
		input2[0][0] = a;
		input2[0][1] = b;
		input2[1][0] = c;
		input2[1][1] = d;
	}
}

void symmeigs(double **input, double *evec_largeeval, double *evec_smalleval, double *evals) {
	double **evecs=malloc(2*sizeof(double*)), **symm=malloc(2*sizeof(double*));
	double det, tr, crossprod;
	evecs[0] = malloc(2*sizeof(double)); evecs[1] = malloc(2*sizeof(double));
	symm[0] = malloc(2*sizeof(double)); symm[1] = malloc(2*sizeof(double));
	symm[0][0] = input[1][0];                     symm[1][1] = -input[0][1];
	symm[0][1] = 0.5*(input[1][1] - input[0][0]); symm[1][0] = symm[0][1];
	linalg2x2(symm, evecs, evals, &det, &tr);
	//printmat("result = input2x2*input2xdims", result, 2, dims);
	if (evals[0] < evals[1]) {
		evec_largeeval[0] = evecs[0][1]; evec_largeeval[1] = evecs[1][1];	
		evec_smalleval[0] = evecs[0][0]; evec_smalleval[1] = evecs[1][0];	
	}
	else {
		evec_smalleval[0] = evecs[0][1]; evec_smalleval[1] = evecs[1][1];	
		evec_largeeval[0] = evecs[0][0]; evec_largeeval[1] = evecs[1][0];	
	}
	crossprod = evec_largeeval[0]*evec_smalleval[1] - evec_largeeval[1]*evec_smalleval[0];
	if (crossprod < 0.0) {
		printf("\n\n\n\n\n\n\n\nYOLO\n\n\n\n\n\n\n");
		evec_smalleval[0] *= (-1);
		evec_smalleval[1] *= (-1);
	}
	free(symm[0]);
	free(symm[1]);
	//free(evecs[0]);
	//free(evecs[1]);
	free(symm);
	//free(evecs);
}

//void gradsymmeigs(double **input, double **gradinput, double *evec_largeeval, double *evec_smalleval, double *gradevec_largeeval, double *gradevec_smalleval) {
//	double *evals=malloc(2*sizeof(double)), **evecs=malloc(2*sizeof(double*)), **symm=malloc(2*sizeof(double*));
//	double det, tr;
//	evecs[0] = malloc(2*sizeof(double)); evecs[1] = malloc(2*sizeof(double));
//	symm[0] = malloc(2*sizeof(double)); symm[1] = malloc(2*sizeof(double));
//	symm[0][0] = input[1][0];                     symm[1][1] = -input[0][1];
//	symm[0][1] = 0.5*(input[1][1] - input[0][0]); symm[1][0] = symm[0][1];
//	linalg2x2(symm, evecs, evals, &det, &tr);
//	//printmat("result = input2x2*input2xdims", result, 2, dims);
//	if (evals[0] < evals[1]) {
//		evec_largeeval[0] = evecs[0][1]; evec_largeeval[1] = evecs[1][1];	
//		evec_smalleval[0] = evecs[0][0]; evec_smalleval[1] = evecs[1][0];	
//	}
//	else {
//		evec_smalleval[0] = evecs[0][1]; evec_smalleval[1] = evecs[1][1];	
//		evec_largeeval[0] = evecs[0][0]; evec_largeeval[1] = evecs[1][0];	
//	}
//	//free(symm[0]);
//	//free(symm[1]);
//	//free(evecs[0]);
//	//free(evecs[1]);
//	//free(symm);
//	//free(evecs);
//	//free(evals);
//}

double inner(double *left_vect, double **matrix, double *right_vect) {
	double answer;
	answer = left_vect[0]*matrix[0][0]*right_vect[0] + left_vect[0]*matrix[0][1]*right_vect[1] 
	       + left_vect[1]*matrix[1][0]*right_vect[0] + left_vect[1]*matrix[1][1]*right_vect[1];
	return answer;
}
