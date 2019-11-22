#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_sf_bessel.h>
#include "isc.h"
#include "input.h"

double **invert2x2(double input2x2[2][2], double* pdet_tangent) {
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
	/* has been checked to work for 2x2 matrix multiplying a 2x1 column vector */
	double **result=malloc(dims*sizeof(double*));
	int i=0;
	result = malloc(2*sizeof(double*));
	for (i=0;i<2;i++) {
		result[i] = malloc(dims*sizeof(double*));
	}
	printmat("input2x2", input2x2, 2, 2);
	printmat("input2xdims", input2xdims, 2, dims);
	for (i=0;i<dims;i++) {
		result[0][i] = input2x2[0][0]*input2xdims[0][i] + input2x2[0][1]*input2xdims[1][i]; //input2xdims[1][i];
		result[1][i] = input2x2[1][0]*input2xdims[0][i] + input2x2[1][1]*input2xdims[1][i]; //input2xdims[1][i];
	}
	printmat("result = input2x2*input2xdims", result, 2, dims);
	return result;
}

void printstruct(char *name, struct position *input) {
	printf("structure position %s:\nposition = (%f, %f)\ntangent = |%10.5f %10.5f|\n          |%10.5f %10.5f|\n", name, input->loc[0], input->loc[1], input->tangent[0][0], input->tangent[0][1], input->tangent[1][0], input->tangent[1][1]);
	//printf("for structure position %s:\nposition = (%f, %f)\n", name, input->loc[0], input->loc[1]);
	//printmat("tangent", input->tangent, 2, 2);
}


struct position *findcentre(double ***coils, int *n_coils, int **n_segs, struct position *fieldline) {
	// declarations
	//clock_t start = clock();
	int i=0;
	struct position fieldline_start, deltafieldline;
	double varphi=0.0, **jumptocentre=calloc(1,sizeof(double*)); //, **testunity;
	double dvarphi = 2.0*M_PI/(N_gridphi_field_period*field_periods);
	//double **matrix;
	double **pdeltafieldline, **pTminusI, det_tangent;
	double **inverseTminusI, **inverseT;
	double error=1.0, errorlimit=0.0000001, factor=1.0;
	int N_line= (int) 3*N_gridphi_field_period;
	struct position *centre=malloc(N_line*sizeof(struct position));
	do {
		fieldline->tangent[0][0]=1.0; fieldline->tangent[0][1]=0.0; fieldline->tangent[1][0]=0.0; fieldline->tangent[1][1]=1.0;
		fieldline_start = *fieldline;
		varphi = 0.0;
		printf("RR=%f\n", fieldline->loc[0]);
		for (i=0; i<N_line; i++)
		{
			centre[i] = *fieldline;
			//if (i%N_gridphi_field_period==0) {
			//	printf("varphi = %f\n", varphi);
			//	printstruct("fieldline[i]\n", fieldline);
			//}
			//printstruct("fieldline[i]\n", fieldline);
			RK4(fieldline, varphi, dvarphi, coils, n_coils, n_segs);
			varphi += dvarphi;
		}
		deltafieldline = addstructs(1.0, fieldline, -1.0, &fieldline_start); 
		//printf("fieldline->tangent[1][1] = %f\n", fieldline->tangent[1][0]);
		inverseTminusI = invert2x2(deltafieldline.tangent, &det_tangent);
		inverseT = invert2x2(fieldline->tangent, &det_tangent);
		//printstruct("fieldline", fieldline); //printf("det(tangent)=%f\n", det_tangent);
		pdeltafieldline = calloc(2,sizeof(double*)); 
		pTminusI = calloc(2,sizeof(double*)); 
		*pdeltafieldline     = &deltafieldline.loc[0];
		*(pdeltafieldline+1) = &deltafieldline.loc[1];
		*pTminusI     = &deltafieldline.tangent[0][0];
		*(pTminusI+1) = &deltafieldline.tangent[1][0];
		//matrix = multiply2x2(inverseTminusI, pTminusI, 2); //printmat("matrix", matrix, 2, 2);
		jumptocentre = multiply2x2(inverseTminusI, pdeltafieldline, 1);
		jumptocentre[0][0] *= factor;
		jumptocentre[1][0] *= factor;
		//printstruct("deltafieldline", &deltafieldline);
		//printmat("jumptocentre", jumptocentre, 2, 1);
		//printstruct("fieldline_start", &fieldline_start);
		printstruct("fieldline", fieldline);
		fieldline->loc[0] = fieldline_start.loc[0] - jumptocentre[0][0]; 
		fieldline->loc[1] = fieldline_start.loc[1] - jumptocentre[1][0];
		error = sqrt(pow(deltafieldline.loc[0], 2.0) + pow(deltafieldline.loc[1], 2.0)); 
		free(jumptocentre);
	} while(error>errorlimit);
	//clock_t int3 = clock();
	return centre;
}

struct position *findisland(double ***coils, int *n_coils, int **n_segs, struct position *fieldline, int tor_mode, int pol_mode) {
	// declarations
	//clock_t start = clock();
	int i=0, N_field_periods;
	struct position fieldline_start, deltafieldline;
	double varphi=0.0, **jumptocentre=calloc(1,sizeof(double*)); //, **testunity;
	double dvarphi = 2.0*M_PI/(N_gridphi_field_period*field_periods);
	double **matrix;
	double **pdeltafieldline, **pTminusI, det_tangent;
	double **inverseTminusI, **inverseT;
	double error=1.0, errorlimit = 0.00000001;
	double factor = 1.0;
	int N_line=0;
	struct position *centre=malloc(N_line*sizeof(struct position));
	if (pol_mode % field_periods == 0) N_field_periods = tor_mode;	
	else 				   N_field_periods = field_periods*tor_mode;
	N_line = N_field_periods*N_gridphi_field_period;
	do {
		fieldline_start = *fieldline;
		varphi = 0.0;
		printf("RR=%f\n", fieldline->loc[0]);
		for (i=0; i<N_line; i++)
		{
			centre[i] = *fieldline;
			if (i%N_gridphi_field_period==0)
			{
				printf("varphi = %f\n", varphi);
				printstruct("fieldline[i]\n", fieldline);
			}
			//printstruct("fieldline[i]\n", fieldline);
			RK4(fieldline, varphi, dvarphi, coils, n_coils, n_segs);
			varphi += dvarphi;
		}
		deltafieldline = addstructs(1.0, fieldline, -1.0, &fieldline_start); 
		//printf("fieldline->tangent[1][1] = %f\n", fieldline->tangent[1][0]);
		inverseTminusI = invert2x2(deltafieldline.tangent, &det_tangent);
		inverseT = invert2x2(fieldline->tangent, &det_tangent);
		printstruct("fieldline", fieldline);
		printf("det(tangent)=%f\n", det_tangent);
//		printf("tangent[0][0] address is %p\n", &&deltafieldline
		
		pdeltafieldline = calloc(2,sizeof(double*)); 
		pTminusI = calloc(2,sizeof(double*)); 
		//p11deltafieldline = &; p12deltafieldline = &deltafieldline->loc[1];
		//*pdeltafieldline = p11deltafieldline;
		//*(pdeltafieldline+1) = p12deltafieldline;
		*pdeltafieldline     = &deltafieldline.loc[0];
		*(pdeltafieldline+1) = &deltafieldline.loc[1];
		*pTminusI     = &deltafieldline.tangent[0][0];
		*(pTminusI+1) = &deltafieldline.tangent[1][0];
		matrix = multiply2x2(inverseTminusI, pTminusI, 2);
		printmat("matrix", matrix, 2, 2);
		jumptocentre = multiply2x2(inverseTminusI, pdeltafieldline, 1);
		jumptocentre[0][0] *= (factor);
		jumptocentre[1][0] *= (factor);
		printstruct("deltafieldline", &deltafieldline);
		printmat("jumptocentre", jumptocentre, 2, 1);
		printstruct("fieldline_start", &fieldline_start);
		printstruct("fieldline", fieldline);
		fieldline->loc[0] = fieldline_start.loc[0] - jumptocentre[0][0]; 
		fieldline->loc[1] = fieldline_start.loc[1] - jumptocentre[1][0];
		fieldline->tangent[0][0]=1.0; fieldline->tangent[0][1]=0.0; fieldline->tangent[1][0]=0.0; fieldline->tangent[1][1]=1.0;
		error = sqrt(pow(deltafieldline.loc[0], 2.0) + pow(deltafieldline.loc[1], 2.0)); 
		printf("error=%f\n", error);
		//free(jumptocentre);
	} while (error > errorlimit);
	//clock_t int3 = clock();
	return centre;
}

//struct ext_position *alongcentre(double RR, double ZZ, int tor_mode, int pol_mode, double ***coils, int *n_coils, int **n_segs) {
//	int i=0, N_field_periods;
//	struct position fieldline_start, fieldline, *centre;
//	double varphi=0.0;
//	double dvarphi = 2.0*M_PI/(N_gridphi_field_period*field_periods);
//	double **inverted, detcentre, evecs[2][2], evals[2], det, trace;
//	double **ptarray;
//	double **full_orbit;
//	int N_line=0, centre_ind, sec_ind;
//	struct ext_position *ext_centre;
//
//	ptarray = malloc(2*sizeof(double*));
//	if (pol_mode % field_periods == 0) N_field_periods = tor_mode;	
//	else 				   N_field_periods = field_periods*tor_mode;
//	N_line = N_field_periods*N_gridphi_field_period;
//	centre = malloc(N_field_periods*sizeof(struct position));
//	ext_centre = malloc(N_field_periods*sizeof(struct ext_position));
//	fieldline_start.loc[0] = RR; fieldline_start.loc[1] = ZZ;
//	fieldline_start.tangent[0][0] = 1.0; fieldline_start.tangent[0][1] = 0.0;
//	fieldline_start.tangent[1][0] = 0.0; fieldline_start.tangent[1][1] = 1.0;
//	fieldline = fieldline_start;
//	centre[0] = fieldline;
//	varphi = 0.0;
//	for (i=1; i<N_line+1; i++) {
//		RK4(&fieldline, varphi, dvarphi, coils, n_coils, n_segs);
//		varphi += dvarphi;
//		if (i%N_gridphi_field_period==0) {
//			centre_ind = (i / N_gridphi_field_period);
//			centre[centre_ind] = fieldline;
//			printf("varphi = %f\n", varphi);
//			printstruct("fieldline[i]\n", &fieldline);
//			inverted = invert2x2(centre[(centre_ind-1) % N_field_periods].tangent, &detcentre);
//			ptarray[0] = &centre[centre_ind % N_field_periods].tangent[0][0];
//			ptarray[1] = &centre[centre_ind % N_field_periods].tangent[1][0];
//			ext_centre[centre_ind % N_field_periods].part_tangent = multiply2x2(ptarray, inverted, 2);
//		}
//	}
//	for (centre_ind=0;centre_ind<N_field_periods;centre_ind++) {
//		ext_centre[centre_ind % N_field_periods].full_tangent = multiply2x2(ext_centre[(centre_ind +1)%N_field_periods].part_tangent, ext_centre[(centre_ind+2)%N_field_periods].part_tangent, 2);
//		for (sec_ind=3;sec_ind<N_field_periods+1;sec_ind++) {
//			ext_centre[centre_ind % N_field_periods].full_tangent = multiply2x2(ext_centre[centre_ind % N_field_periods].full_tangent, ext_centre[(centre_ind + sec_ind)%N_field_periods].part_tangent, 2);
//		}
//	}
//	for (centre_ind=0;centre_ind<N_field_periods;centre_ind++) {
//		linalg2x2(ext_centre[centre_ind].full_tangent, evecs, evals, &det, &trace);
//		if ((fabs(evecs[0][0]) < small) || (fabs(evecs[0][1]) < small)) {
//			ext_centre[centre_ind].angle = evals[1];
//		}
//	}
//	//ext_centre[centre_index % N_field_periods].part_tangent = multiply2x2(centre[centre_index % N_field_periods].tangent, inverted, 2);
//	
//	//clock_t int3 = clock();
//	return ext_centre;
//} 
