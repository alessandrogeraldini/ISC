#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_sf_bessel.h>
#include "isc.h"
#include "input.h"

void printstruct(char *name, struct position *input) {
	printf("structure position %s:\nposition = (%f, %f)\ntangent = |%10.5f %10.5f|\n          |%10.5f %10.5f|\n", 
	 name, input->loc[0], input->loc[1], input->tangent[0][0], input->tangent[0][1], input->tangent[1][0], input->tangent[1][1]);
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
	double **pdeltafieldline, det_tangent;
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
		*pdeltafieldline     = &deltafieldline.loc[0];
		*(pdeltafieldline+1) = &deltafieldline.loc[1];
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
	double **pdeltafieldline, det_tangent;
	double **inverseTminusI, **inverseT;
	double error=1.0, errorlimit = 0.00000001;
	double factor = 1.0;
	int N_line=0;
	struct position *centre=malloc(N_line*sizeof(struct position));
	if (pol_mode % field_periods == 0) N_field_periods = pol_mode;	
	else 				   N_field_periods = field_periods*pol_mode;
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
		//pTminusI = calloc(2,sizeof(double*)); 
		//p11deltafieldline = &; p12deltafieldline = &deltafieldline->loc[1];
		//*pdeltafieldline = p11deltafieldline;
		//*(pdeltafieldline+1) = p12deltafieldline;
		*pdeltafieldline     = &deltafieldline.loc[0];
		*(pdeltafieldline+1) = &deltafieldline.loc[1];
		//*pTminusI     = &deltafieldline.tangent[0][0];
		//*(pTminusI+1) = &deltafieldline.tangent[1][0];
		matrix = multiply2x2(inverseTminusI, deltafieldline.tangent, 2);
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

struct ext_position *alongcentre(double RR, double ZZ, int tor_mode, int pol_mode, double ***coils, int *n_coils, int **n_segs) {
	int i=0, N_field_periods, q0;
	struct position fieldline_start, fieldline, *centre;
	double varphi=0.0, omega;
	double dvarphi = 2.0*M_PI/(N_gridphi_field_period*field_periods);
	double **inverted, detcentre, **evecs1, *evals1, det=0, trace=0;
	double circumference;
	int N_line=0, main_ind, centre_ind, sec_ind;
	struct ext_position *ext_centre;

	evecs1 = malloc(2*sizeof(double*));
	evecs1[0] = malloc(2*sizeof(double)); evecs1[1] = malloc(2*sizeof(double));
	evals1 = malloc(2*sizeof(double));

	if (tor_mode % field_periods == 0) N_field_periods = pol_mode;	
	else 				   N_field_periods = field_periods*pol_mode;
	N_line = N_field_periods*N_gridphi_field_period;
	centre = malloc(N_field_periods*sizeof(struct position));
	ext_centre = malloc(N_field_periods*sizeof(struct ext_position));
	fieldline_start.loc[0] = RR; fieldline_start.loc[1] = ZZ;
	fieldline_start.tangent = set_identity();
	fieldline = fieldline_start;
	centre[0] = fieldline;
	varphi = 0.0;
	printf("varphi = %f\n", varphi);
	printstruct("fieldline[i]\n", &fieldline);
	for (i=1; i<N_line+1; i++) {
		RK4(&fieldline, varphi, dvarphi, coils, n_coils, n_segs);
		varphi += dvarphi;
		if (i%N_gridphi_field_period==0) {
			centre_ind = (i / N_gridphi_field_period);
			centre[centre_ind % N_field_periods] = fieldline;
			printf("varphi = %f\n", varphi);
			printstruct("fieldline[i]\n", &fieldline);
			inverted = invert2x2(centre[(centre_ind-1) % N_field_periods].tangent, &detcentre);
			printmat("inverted", inverted, 2, 2);
			//ptarray[0] = &centre[centre_ind % N_field_periods].tangent[0][0];
			//ptarray[1] = &centre[centre_ind % N_field_periods].tangent[1][0];
			ext_centre[centre_ind % N_field_periods].loc[0] = fieldline.loc[0];
			ext_centre[centre_ind % N_field_periods].loc[1] = fieldline.loc[1];
			ext_centre[centre_ind % N_field_periods].part_tangent = multiply2x2(centre[centre_ind % N_field_periods].tangent, inverted, 2);
			printmat("ext_centre.part_tangent", ext_centre[centre_ind % N_field_periods].part_tangent, 2, 2);
		}
	}
	circumference = 0.0;
	for (centre_ind=0;centre_ind<N_field_periods;centre_ind++) {
		circumference += sqrt(pow(ext_centre[(centre_ind+1)%N_field_periods].loc[0] 
					        - ext_centre[centre_ind].loc[0], 2.0)
				    	    + pow(ext_centre[(centre_ind+1)%N_field_periods].loc[1] 
						- ext_centre[centre_ind].loc[1], 2.0));
		ext_centre[centre_ind].epar = malloc(2*sizeof(double)); ext_centre[centre_ind].eperp = malloc(2*sizeof(double));
		ext_centre[centre_ind].full_tangent = multiply2x2(ext_centre[(centre_ind +2)%N_field_periods].part_tangent, ext_centre[(centre_ind+1)%N_field_periods].part_tangent, 2);
		for (sec_ind=3;sec_ind<N_field_periods+1;sec_ind++) {
			multiply2x2reassign(ext_centre[(centre_ind + sec_ind)%N_field_periods].part_tangent, ext_centre[centre_ind].full_tangent, 2);
		}
		printmat("ext_centre.fulltangent", ext_centre[centre_ind].full_tangent, 2, 2);
		symmeigs(ext_centre[centre_ind].full_tangent, ext_centre[centre_ind].eperp, ext_centre[centre_ind].epar);
		printf("eigenvectors are (%f, %f) and (%f, %f)\n", ext_centre[centre_ind].eperp[0], ext_centre[centre_ind].eperp[1], ext_centre[centre_ind].epar[0], ext_centre[centre_ind].epar[1]);
	}
	for (main_ind=0;main_ind<N_field_periods; main_ind++) {
		linalg2x2(ext_centre[main_ind].full_tangent, evecs1, evals1, &det, &trace);
		ext_centre[main_ind].circumference = circumference;
		if ((fabs(evecs1[0][0]) < small) || (fabs(evecs1[0][1]) < small)) {
			ext_centre[main_ind].angle = evals1[1];
			omega = field_periods*evals1[1]/(2.0*M_PI*N_field_periods); // I have changed it compared to Cary and Hanson's paper
			q0 = (int) (field_periods/(4.0*omega) - N_field_periods/2.0 + 0.5);
			ext_centre[main_ind].q0_index = q0;
		}
		printf("angle=%f\n", ext_centre[main_ind].angle);
		printf("q0=%d\n", ext_centre[main_ind].q0_index);
		ext_centre[main_ind].long_tangent = malloc(N_field_periods*sizeof(double));
		for (centre_ind=0;centre_ind<N_field_periods;centre_ind++) {
			//ext_centre[centre_ind].long_tangent = set_identity();
			//for (sec_ind=0; sec_ind<q0/N_field_periods; sec_ind++) {
			//	multiply2x2reassign(ext_centre[0].full_tangent, ext_centre[centre_ind].long_tangent, 2); 
			//}
			//for (sec_ind=0; sec_ind<q0%N_field_periods; sec_ind++) {
			//	multiply2x2reassign(ext_centre[sec_ind+1].part_tangent, ext_centre[centre_ind].long_tangent, 2); 
			//}
			ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%N_field_periods] = set_identity();
			for (sec_ind=0; sec_ind<(q0+centre_ind)/N_field_periods; sec_ind++) {
				multiply2x2reassign(ext_centre[main_ind].full_tangent, ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%N_field_periods], 2); 
			}
			for (sec_ind=0; sec_ind<(q0+centre_ind)%N_field_periods; sec_ind++) {
				multiply2x2reassign(ext_centre[(sec_ind+1+main_ind)%N_field_periods].part_tangent, ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%N_field_periods], 2); 
			}
		}
	}
	//ext_centre[centre_index % N_field_periods].part_tangent = multiply2x2(centre[centre_index % N_field_periods].tangent, inverted, 2);
	//clock_t int3 = clock();
	return ext_centre;
} 

double *islandwidth(struct ext_position *ext_fieldline, int tor_mode, int pol_mode) {
	int main_index, centre_index, N_field_periods;
	double *wperp, sum_matrix_elements, matrix_element;
	double circumference=ext_fieldline[0].circumference;
	if (tor_mode % field_periods == 0) N_field_periods = pol_mode;	
	else 				   N_field_periods = field_periods*pol_mode;
	wperp = malloc(N_field_periods*sizeof(double));
	for (main_index=0;main_index<N_field_periods; main_index++) {
		//circumference = 0.0;		
		sum_matrix_elements = 0.0;		
		for (centre_index=0;centre_index<N_field_periods; centre_index++) {
			matrix_element = inner(ext_fieldline[centre_index].epar, 
					       ext_fieldline[main_index].long_tangent[centre_index], 
					       ext_fieldline[main_index].eperp);	
			sum_matrix_elements += fabs(matrix_element);
			//circumference += sqrt(pow(ext_fieldline[(centre_index+1)%N_field_periods].loc[0] 
			//		        - ext_fieldline[centre_index].loc[0], 2.0)
			//	    	    + pow(ext_fieldline[(centre_index+1)%N_field_periods].loc[1] 
			//			- ext_fieldline[centre_index].loc[1], 2.0));
			//printf("circumference = %f\nmatrix_element = %f\n", circumference, matrix_element);
			//printf("other circumference = %f\n", ext_fieldline[centre_index].circumference);
		}
		wperp[main_index] = 2.0*N_field_periods*circumference/(M_PI*pol_mode*sum_matrix_elements);
		printf("width = %f for index= %d\n", wperp[main_index], main_index);
	}
	return wperp;
}
