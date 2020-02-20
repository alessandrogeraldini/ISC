/* Author: Alessandro Geraldini */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_sf_bessel.h>
#include "isc.h"

void printstruct(char *name, struct position *input) {
	printf("structure position %s:\nposition = (%10.8f, %10.8f)\ntangent = |%10.8f %10.8f|\n          |%10.8f %10.8f|\n", 
	 name, input->loc[0], input->loc[1], input->tangent[0][0], input->tangent[0][1], input->tangent[1][0], input->tangent[1][1]);
	//printf("for structure position %s:\nposition = (%f, %f)\n", name, input->loc[0], input->loc[1]);
	//printmat("tangent", input->tangent, 2, 2);
}


struct position *findcentre(double ***coils, int *n_coils, int **n_segs, struct position *fieldline, int N_gridphi_toroidal) {
	// declarations
	//clock_t start = clock();
	int i=0;
	struct position fieldline_start, deltafieldline;
	double varphi=0.0, **jumptocentre=calloc(1,sizeof(double*)); //, **testunity;
	double dvarphi = 2.0*M_PI/(N_gridphi_toroidal);
	//double **matrix;
	double **pdeltafieldline, det_tangent;
	double **inverseTminusI, **inverseT;
	double error=1.0, errorlimit=0.00000001, factor=1.0;
	int ll=1;
	struct position *centre=malloc(N_gridphi_toroidal*sizeof(struct position));
	do {
		fieldline->tangent[0][0]=1.0; fieldline->tangent[0][1]=0.0; fieldline->tangent[1][0]=0.0; fieldline->tangent[1][1]=1.0;
		fieldline_start = *fieldline;
		varphi = 0.0;
		printf("RR=%f\n", fieldline->loc[0]);
		for (i=0; i<N_gridphi_toroidal; i++)
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
		*pdeltafieldline      = &deltafieldline.loc[0];
		*(pdeltafieldline+1)  = &deltafieldline.loc[1];
		jumptocentre          = multiply2x2(inverseTminusI, pdeltafieldline, 1);
		jumptocentre[0][0]   *= factor;
		jumptocentre[1][0]   *= factor;
		//printstruct("deltafieldline", &deltafieldline);
		//printmat("jumptocentre", jumptocentre, 2, 1);
		//printstruct("fieldline_start", &fieldline_start);
		printstruct("fieldline", fieldline);
		fieldline->loc[0] = fieldline_start.loc[0] - jumptocentre[0][0]; 
		fieldline->loc[1] = fieldline_start.loc[1] - jumptocentre[1][0];
		error = sqrt(pow(deltafieldline.loc[0], 2.0) + pow(deltafieldline.loc[1], 2.0)); 
		//free(jumptocentre[0]);
		//free(jumptocentre[1]);
		//free(jumptocentre);
		i = 1;
	} while(error>errorlimit && ll==0);
	//clock_t int3 = clock();
	return centre;
}

struct position *findisland(double ***coils, int *n_coils, int **n_segs, struct position *fieldline, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode) {
	// declarations
	//clock_t start = clock();
	int i=0;
	struct position fieldline_start, deltafieldline;
	double varphi=0.0, **jumptocentre=calloc(1,sizeof(double*)); //, **testunity;
	double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
	double **matrix;
	double **pdeltafieldline, det_tangent;
	double **inverseTminusI, **inverseT;
	double error = 1.0, errorlimit = 0.000000000001;
	double factor = 1.0;
	int N_line=0, L_field_periods;
	struct position *centre=malloc(N_line*sizeof(struct position));
	if (tor_mode % m0_symmetry == 0) 
		L_field_periods = pol_mode;	
	else 				 
		L_field_periods = m0_symmetry*pol_mode;
	N_line = L_field_periods*N_gridphi_per_field_period;
	printf("N_line = %d\n", N_line);
	printf("field_periods = %d\n", m0_symmetry);
	do {
		fieldline_start = *fieldline;
		varphi = 0.0;
		printf("RR=%f\n", fieldline->loc[0]);
		for (i=0; i<N_line; i++)
		{
			centre[i] = *fieldline;
			//printf("where do I seg fault?\n");
			if (i%N_gridphi_per_field_period==0)
			{
				printf("varphi = %f\n", varphi);
				printstruct("fieldline\n", fieldline);
				//printf("where do I seg fault?\n");
			}
			//printstruct("fieldline[i]\n", fieldline);
			RK4(fieldline, varphi, dvarphi, coils, n_coils, n_segs);
			varphi += dvarphi;
			//printf("where do I seg fault?\n");
		}
		//printf("where do I seg fault?\n");
		deltafieldline = addstructs(1.0, fieldline, -1.0, &fieldline_start); 
		//printf("fieldline->tangent[1][1] = %f\n", fieldline->tangent[1][0]);
		inverseTminusI = invert2x2(deltafieldline.tangent, &det_tangent);
		inverseT = invert2x2(fieldline->tangent, &det_tangent);
		printf("varphi = %f\n", varphi);
		printstruct("fieldline", fieldline);
		printf("det(tangent)=%f\n", det_tangent);
//		printf("tangent[0][0] address is %p\n", &&deltafieldline
		printmat("inverseTminusI", inverseTminusI, 2, 2);
		
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
	} while ((deltafieldline.loc[0] > errorlimit) || (deltafieldline.loc[1] > errorlimit));
	//clock_t int3 = clock();
	return centre;
}

struct ext_position *alongcentre(double RR, double ZZ, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode, double ***coils, int *n_coils, int **n_segs) {
	int i=0, q0, L_field_periods;
	struct position fieldline_start, fieldline, *centre;
	double varphi=0.0, omega;
	double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
	double **inverted, detcentre, **evecs1, *evals1, det=0, trace=0;
	double circumference;
	double *evals=malloc(2*sizeof(double));
	int N_line=0, main_ind, centre_ind, sec_ind;
	struct ext_position *ext_centre;

	evecs1 = malloc(2*sizeof(double*));
	evecs1[0] = malloc(2*sizeof(double)); evecs1[1] = malloc(2*sizeof(double));
	evals1 = malloc(2*sizeof(double));

	if (tor_mode % m0_symmetry == 0) 
		L_field_periods = pol_mode;	
	else 				   
		L_field_periods = m0_symmetry*pol_mode;
	N_line = L_field_periods*N_gridphi_per_field_period;
	centre = malloc(L_field_periods*sizeof(struct position));
	ext_centre = malloc(L_field_periods*sizeof(struct ext_position));
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
		if (i%N_gridphi_per_field_period==0) {
			centre_ind = (i / N_gridphi_per_field_period);
			centre[centre_ind % L_field_periods] = fieldline;
			printf("varphi = %f\n", varphi);
			printstruct("fieldline[i]\n", &fieldline);
			inverted = invert2x2(centre[(centre_ind-1) % L_field_periods].tangent, &detcentre);
			printmat("inverted", inverted, 2, 2);
			//ptarray[0] = &centre[centre_ind % L_field_periods].tangent[0][0];
			//ptarray[1] = &centre[centre_ind % L_field_periods].tangent[1][0];
			ext_centre[centre_ind % L_field_periods].loc[0] = fieldline.loc[0];
			ext_centre[centre_ind % L_field_periods].loc[1] = fieldline.loc[1];
			ext_centre[centre_ind % L_field_periods].part_tangent = multiply2x2(centre[centre_ind % L_field_periods].tangent, inverted, 2);
			printmat("ext_centre.part_tangent", ext_centre[centre_ind % L_field_periods].part_tangent, 2, 2);
		}
	}
	circumference = 0.0;
	for (centre_ind=0;centre_ind<L_field_periods;centre_ind++) {
		circumference += sqrt(pow(ext_centre[(centre_ind+1)%L_field_periods].loc[0] 
					        - ext_centre[centre_ind].loc[0], 2.0)
				    	    + pow(ext_centre[(centre_ind+1)%L_field_periods].loc[1] 
						- ext_centre[centre_ind].loc[1], 2.0));
		ext_centre[centre_ind].epar = malloc(2*sizeof(double)); ext_centre[centre_ind].eperp = malloc(2*sizeof(double));
		ext_centre[centre_ind].full_tangent = multiply2x2(ext_centre[(centre_ind +2)%L_field_periods].part_tangent, ext_centre[(centre_ind+1)%L_field_periods].part_tangent, 2);
		for (sec_ind=3;sec_ind<L_field_periods+1;sec_ind++) {
			multiply2x2reassign(ext_centre[(centre_ind + sec_ind)%L_field_periods].part_tangent, ext_centre[centre_ind].full_tangent, 2);
		}
		printmat("ext_centre.fulltangent", ext_centre[centre_ind].full_tangent, 2, 2);
		symmeigs(ext_centre[centre_ind].full_tangent, ext_centre[centre_ind].eperp, ext_centre[centre_ind].epar, evals);
		printf("eigenvectors are (%f, %f) and (%f, %f)\n", ext_centre[centre_ind].eperp[0], ext_centre[centre_ind].eperp[1], ext_centre[centre_ind].epar[0], ext_centre[centre_ind].epar[1]);
	}
	for (main_ind=0;main_ind<L_field_periods; main_ind++) {
		linalg2x2(ext_centre[main_ind].full_tangent, evecs1, evals1, &det, &trace);
		ext_centre[main_ind].circumference = circumference;
		if ((fabs(evecs1[0][0]) < small) || (fabs(evecs1[0][1]) < small)) {
			ext_centre[main_ind].angle = evals1[1];
			omega = m0_symmetry*evals1[1]/(2.0*M_PI*L_field_periods); // I have changed it compared to Cary and Hanson's paper
			q0 = (int) (m0_symmetry/(4.0*omega) - L_field_periods/2.0 + 0.5);
			ext_centre[main_ind].q0_index = q0;
		}
		printf("angle=%f\n", ext_centre[main_ind].angle);
		printf("q0=%d\n", ext_centre[main_ind].q0_index);
		ext_centre[main_ind].long_tangent = malloc(L_field_periods*sizeof(double));
		for (centre_ind=0;centre_ind<L_field_periods;centre_ind++) {
			//ext_centre[centre_ind].long_tangent = set_identity();
			//for (sec_ind=0; sec_ind<q0/L_field_periods; sec_ind++) {
			//	multiply2x2reassign(ext_centre[0].full_tangent, ext_centre[centre_ind].long_tangent, 2); 
			//}
			//for (sec_ind=0; sec_ind<q0%L_field_periods; sec_ind++) {
			//	multiply2x2reassign(ext_centre[sec_ind+1].part_tangent, ext_centre[centre_ind].long_tangent, 2); 
			//}
			ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_field_periods] = set_identity();
			for (sec_ind=0; sec_ind<(q0+centre_ind)/L_field_periods; sec_ind++) {
				multiply2x2reassign(ext_centre[main_ind].full_tangent, ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_field_periods], 2); 
			}
			for (sec_ind=0; sec_ind<(q0+centre_ind)%L_field_periods; sec_ind++) {
				multiply2x2reassign(ext_centre[(sec_ind+1+main_ind)%L_field_periods].part_tangent, ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_field_periods], 2); 
			}
		}
	}
	//ext_centre[centre_index % L_field_periods].part_tangent = multiply2x2(centre[centre_index % L_field_periods].tangent, inverted, 2);
	//clock_t int3 = clock();
	return ext_centre;
} 

double *islandwidth(struct ext_position *ext_fieldline, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode) {
	int main_index, centre_index, L_field_periods;
	double *wperp, sum_matrix_elements, matrix_element;
	double circumference=ext_fieldline[0].circumference;
	printf("tor_mode=%d\n", tor_mode);
	if (tor_mode % m0_symmetry == 0) L_field_periods = pol_mode;	
	else 				   L_field_periods = m0_symmetry*pol_mode;
	printf("L_field_periods=%d\n", L_field_periods);
	wperp = malloc(L_field_periods*sizeof(double));
	for (main_index=0;main_index<L_field_periods; main_index++) {
		//circumference = 0.0;		
		sum_matrix_elements = 0.0;		
		for (centre_index=0;centre_index<L_field_periods; centre_index++) {
			matrix_element = inner(ext_fieldline[centre_index].epar, 
					       ext_fieldline[main_index].long_tangent[centre_index], 
					       ext_fieldline[main_index].eperp);	
			sum_matrix_elements += fabs(matrix_element);
			//sum_matrix_elements += fabs(matrix_element)/sin(ext_fieldline[main_index].angle*ext_fieldline[main_index].q0_index/L_field_periods); //This has the sin piece in the denominator, which can be safely set to unity
			//circumference += sqrt(pow(ext_fieldline[(centre_index+1)%L_field_periods].loc[0] 
			//		        - ext_fieldline[centre_index].loc[0], 2.0)
			//	    	    + pow(ext_fieldline[(centre_index+1)%L_field_periods].loc[1] 
			//			- ext_fieldline[centre_index].loc[1], 2.0));
			//printf("circumference = %f\nmatrix_element = %f\n", circumference, matrix_element);
			//printf("other circumference = %f\n", ext_fieldline[centre_index].circumference);
		}
		//circumference = 2.0*M_PI*0.2102;
		//circumference = 2.0*M_PI*0.2094;
		printf("circumference = %f\n", circumference);
		wperp[main_index] = 2.0*L_field_periods*circumference/(M_PI*pol_mode*sum_matrix_elements);
		printf("width = %f for index= %d\n", wperp[main_index], main_index);
	}
	return wperp;
}

struct position *adjfindisland(double ***coils, int *n_coils, int **n_segs, struct position *fieldline, struct position *lambda, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode) {
	// declarations
	//clock_t start = clock();
	int i=0;
	struct position lambda_start, deltalambda, fieldline_start;
	double varphi=0.0, **jumptocentre=calloc(1,sizeof(double*)); //, **testunity;
	double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
	double **matrix;
	double **pdeltalambda, det_tangent;
	double **inverseTminusI, **inverseT;
	double error = 1.0, errorlimit = 0.000000000001;
	double factor = 1.0;
	int N_line=0, L_field_periods;
	struct position *centrelambda=malloc(N_line*sizeof(struct position));
	if (tor_mode % m0_symmetry == 0) 
		L_field_periods = pol_mode;	
	else 				 
		L_field_periods = m0_symmetry*pol_mode;
	N_line = L_field_periods*N_gridphi_per_field_period;
	printf("N_line = %d\n", N_line);
	printf("field_periods = %d\n", m0_symmetry);
	fieldline_start = *fieldline;
	do {
		fieldline->loc[0] = fieldline_start.loc[0];
		fieldline->loc[1] = fieldline_start.loc[1];
		lambda_start = *lambda;
		varphi = 0.0;
		printf("RR=%f\n", fieldline->loc[0]);
		for (i=0; i<N_line; i++)
		{
			centrelambda[i] = *lambda;
			if (i%N_gridphi_per_field_period==0)
			{
				printf("varphi = %f\n", varphi);
				//printstruct("lambda\n", lambda);
				if (i/N_gridphi_per_field_period==0) lambda->loc[0] += 1.0;
			}
			RK4_wadj(fieldline, lambda, varphi, dvarphi, coils, n_coils, n_segs);
			varphi += dvarphi;
		}
		deltalambda = addstructs(1.0, lambda, -1.0, &lambda_start); 
		//printf("fieldline->tangent[1][1] = %f\n", fieldline->tangent[1][0]);
		inverseTminusI = invert2x2(deltalambda.tangent, &det_tangent);
		inverseT = invert2x2(lambda->tangent, &det_tangent);
		//printf("varphi = %f\n", varphi);
		//printstruct("fieldline", lambda);
		//printf("det(tangent)=%f\n", det_tangent);
//		printf("tangent[0][0] address is %p\n", &&deltalambda
		//printmat("inverseTminusI", inverseTminusI, 2, 2);
		
		pdeltalambda = calloc(2,sizeof(double*)); 
		*pdeltalambda     = &deltalambda.loc[0];
		*(pdeltalambda+1) = &deltalambda.loc[1];
		matrix = multiply2x2(inverseTminusI, deltalambda.tangent, 2);
		printmat("matrix", matrix, 2, 2);
		jumptocentre = multiply2x2(inverseTminusI, pdeltalambda, 1);
		jumptocentre[0][0] *= (factor);
		jumptocentre[1][0] *= (factor);
		printstruct("deltalambda", &deltalambda);
		printmat("jumptocentre", jumptocentre, 2, 1);
		printstruct("lambda_start", &lambda_start);
		printstruct("lambda", lambda);
		lambda->loc[0] = lambda_start.loc[0] - jumptocentre[0][0]; 
		lambda->loc[1] = lambda_start.loc[1] - jumptocentre[1][0];
		lambda->tangent[0][0] = 1.0; lambda->tangent[0][1]=0.0; lambda->tangent[1][0]=0.0; lambda->tangent[1][1]=1.0;
		error = sqrt(pow(deltalambda.loc[0], 2.0) + pow(deltalambda.loc[1], 2.0)); 
		printf("%f\n", deltalambda.loc[0]);
		printf("%f\n", deltalambda.loc[1]);
		//free(jumptocentre);
	} while ((fabs(deltalambda.loc[0]) > errorlimit) || (fabs(deltalambda.loc[1]) > errorlimit));
	//clock_t int3 = clock();
	return centrelambda;
}

double *adjevaluate(double ***coils, int *n_coils, int **n_segs, struct position *centre, struct position *lambda_centre, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode) {
	// declarations
	//clock_t start = clock();
	int i=0;
	double varphi=0.0;
	double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
	double *number, initial=0.0;
	int N_line=0, L_field_periods;
	if (tor_mode % m0_symmetry == 0) 
		L_field_periods = pol_mode;	
	else 				 
		L_field_periods = m0_symmetry*pol_mode;

	number = &initial;
	N_line = L_field_periods*N_gridphi_per_field_period;
	printf("N_line = %d\n", N_line);
	printf("field_periods = %d\n", m0_symmetry);
	varphi = 0.0;
	for (i=0; i<N_line; i++)
	{
		if (i%N_gridphi_per_field_period==0)
		{
			printf("varphi = %f\n", varphi);
			printstruct("lambda",lambda_centre);
			printstruct("Xp",centre);
			printf("number = %f\n", *number);
			if (i/N_gridphi_per_field_period==0) lambda_centre->loc[0] += 1.0;
		}
		RK4_adjeval(number, centre, lambda_centre, varphi, dvarphi, coils, n_coils, n_segs);
		varphi += dvarphi;
		//printf("number = %f\n", *number);
	}
	printf("varphi = %f\n", varphi);
	printstruct("lambda",lambda_centre);
	printstruct("Xp",centre);
	printf("number = %f\n", *number);
	return number;
}

struct position *gradcentre(double RR, double ZZ, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode, double ***coils, int *n_coils, int **n_segs) {
	int i=0, L_field_periods;
	struct position fieldline_start, fieldline, gradfieldline, *gradfieldlineret, deltafieldline;
	double varphi=0.0, gradX[2], **pdeltafieldline, **inverseT, **inverseTminusI, det_tangent=0.0;
	double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
	int N_line=0;

	//evecs1 = malloc(2*sizeof(double*));
	//evecs1[0] = malloc(2*sizeof(double)); evecs1[1] = malloc(2*sizeof(double));
	//evals1 = malloc(2*sizeof(double));

	if (tor_mode % m0_symmetry == 0) 
		L_field_periods = pol_mode;	
	else 				   
		L_field_periods = m0_symmetry*pol_mode;

	N_line = L_field_periods*N_gridphi_per_field_period;
	fieldline_start.loc[0] = RR; fieldline_start.loc[1] = ZZ;
	fieldline_start.tangent = set_identity();
	printf("fieldline_start = (%f, %f)\n", RR, ZZ);
	printf("RR = %f, ZZ = %f\n", RR, ZZ);
	fieldline = fieldline_start;
	gradfieldline = fieldline_start;
	gradfieldline.loc[0] = 0.0; gradfieldline.loc[1] = 0.0;
	gradfieldline.tangent = set_zeros();
	varphi = 0.0;
	printf("varphi = %f\n", varphi);
	//for (i=0; i<8*N_line + N_line/2; i++) 
	for (i=0; i<8*N_line; i++) {
		RK4_wgrad(&fieldline, &gradfieldline, varphi, dvarphi, coils, n_coils, n_segs);
		varphi += dvarphi;
		printstruct("gradfieldline\n", &gradfieldline);
		//if (i%N_gridphi_per_field_period==0) {
		//	centre_ind = (i / N_gridphi_per_field_period);
		//	centre[centre_ind % L_field_periods] = fieldline;
		//	printf("varphi = %f\n", varphi);
		//	printstruct("fieldline[i]\n", &fieldline);
		//	inverted = invert2x2(centre[(centre_ind-1) % L_field_periods].tangent, &detcentre);
		//	printmat("inverted", inverted, 2, 2);
		//	//ptarray[0] = &centre[centre_ind % L_field_periods].tangent[0][0];
		//	//ptarray[1] = &centre[centre_ind % L_field_periods].tangent[1][0];
		//}
	}
	//gradfieldline.loc[0]  = 0.00133;
	//gradfieldline.loc[1]  = -0.09218;
	gradfieldlineret = &gradfieldline;
	printf("(%f, %f)\n", gradfieldlineret->loc[0], gradfieldlineret->loc[1]);
	deltafieldline = addstructs(-1.0, &fieldline, 1.0, &fieldline_start); 
	//printmat("deltafieldline", deltafieldline.tangent, 2, 2);
	//printf("fieldline->tangent[1][1] = %f\n", fieldline->tangent[1][0]);
	inverseT = invert2x2(fieldline.tangent, &det_tangent);
	printf("det = %f\n", det_tangent);
	inverseTminusI = invert2x2(deltafieldline.tangent, &det_tangent);
	printf("det = %f\n", det_tangent);
	
	pdeltafieldline = calloc(2,sizeof(double*)); 
	*pdeltafieldline      = &gradfieldline.loc[0];
	*(pdeltafieldline+1)  = &gradfieldline.loc[1];

	//gradX = multiply2x2(inverseTminusI, pdeltafieldline, 1);
	gradX[0] = inverseTminusI[0][0]*gradfieldline.loc[0] + inverseTminusI[0][1]*gradfieldline.loc[1]; 
	gradX[1] = inverseTminusI[1][0]*gradfieldline.loc[0] + inverseTminusI[1][1]*gradfieldline.loc[1]; 
	printf("%f \n", inverseTminusI[0][0]*gradfieldline.loc[0]); 
	printf("%f \n", inverseTminusI[0][1]*gradfieldline.loc[1]); 
	printf("%f \n", inverseTminusI[1][0]*gradfieldline.loc[0]); 
	printf("%f \n", inverseTminusI[1][1]*gradfieldline.loc[1]); 
	printf("%f \n", inverseTminusI[0][0]); 
	printf("%f \n", inverseTminusI[0][1]); 
	printf("%f \n", inverseTminusI[1][0]); 
	printf("%f \n", inverseTminusI[1][1]); 

	printf("(%f, %f)\n", gradX[0], gradX[1]);
	gradfieldline.loc[0] = gradX[0];
	gradfieldline.loc[1] = gradX[1];
	gradfieldlineret = &gradfieldline;
	printf("(%f, %f)\n", gradfieldlineret->loc[0], gradfieldlineret->loc[1]);
	return gradfieldlineret;
}

//struct ext_position *gradalongcentre(double RR, double ZZ, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode, double ***coils, int *n_coils, int **n_segs) {
//	int i=0, q0, L_field_periods;
//	struct position fieldline_start, fieldline, gradfieldline, *centre, *gradcentre;
//	double varphi=0.0, omega, deltadet;
//	double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
//	double **inverted, detcentre, **evecs1, *evals1, det=0, trace=0, **spare_matrix;
//	double circumference;
//	double *evals=malloc(2*sizeof(double));
//	int N_line=0, main_ind, centre_ind, sec_ind;
//	struct ext_position *ext_centre, *grad_ext_centre;
//
//	evecs1 = malloc(2*sizeof(double*));
//	evecs1[0] = malloc(2*sizeof(double)); evecs1[1] = malloc(2*sizeof(double));
//	evals1 = malloc(2*sizeof(double));
//
//	if (tor_mode % m0_symmetry == 0) 
//		L_field_periods = pol_mode;	
//	else 				   
//		L_field_periods = m0_symmetry*pol_mode;
//	N_line = L_field_periods*N_gridphi_per_field_period;
//	centre = malloc(L_field_periods*sizeof(struct position));
//	gradcentre = malloc(L_field_periods*sizeof(struct position));
//	ext_centre = malloc(L_field_periods*sizeof(struct ext_position));
//	grad_ext_centre = malloc(L_field_periods*sizeof(struct ext_position));
//	fieldline_start.loc[0] = RR; fieldline_start.loc[1] = ZZ;
//	fieldline_start.tangent = set_identity();
//	fieldline = fieldline_start;
//	gradfieldline.loc[0] = 0.0; gradfieldline.loc[1] = 0.0;
//	gradfieldline.tangent = set_zeros();
//	centre[0] = fieldline;
//	gradcentre[0] = gradfieldline;
//	varphi = 0.0;
//	printf("varphi = %f\n", varphi);
//	printstruct("fieldline[i]\n", &fieldline);
//	for (i=1; i<N_line+1; i++) {
//		//RK4(&fieldline, varphi, dvarphi, coils, n_coils, n_segs);
//		RK4_wgrad(&fieldline, &gradfieldline, varphi, dvarphi, coils, n_coils, n_segs);
//		varphi += dvarphi;
//		if (i%N_gridphi_per_field_period==0) {
//			centre_ind = (i / N_gridphi_per_field_period);
//			centre[centre_ind % L_field_periods] = fieldline;
//			gradcentre[centre_ind % L_field_periods] = gradfieldline;
//			printf("varphi = %f\n", varphi);
//			printstruct("fieldline[i]\n", &fieldline);
//			inverted = invert2x2(centre[(centre_ind-1) % L_field_periods].tangent, &detcentre);
//			printmat("inverted", inverted, 2, 2);
//			//ptarray[0] = &centre[centre_ind % L_field_periods].tangent[0][0];
//			//ptarray[1] = &centre[centre_ind % L_field_periods].tangent[1][0];
//			ext_centre[centre_ind % L_field_periods].loc[0] = fieldline.loc[0];
//			ext_centre[centre_ind % L_field_periods].loc[1] = fieldline.loc[1];
//			ext_centre[centre_ind % L_field_periods].part_tangent = multiply2x2(centre[centre_ind % L_field_periods].tangent, inverted, 2);
//			//printmat("ext_centre.part_tangent", ext_centre[centre_ind % L_field_periods].part_tangent, 2, 2);
//			grad_ext_centre[centre_ind % L_field_periods].loc[0] = gradfieldline.loc[0];
//			grad_ext_centre[centre_ind % L_field_periods].loc[1] = gradfieldline.loc[1];
//			grad_ext_centre[centre_ind % L_field_periods].part_tangent = add2x2( 
//			1.0 , multiply2x2(gradcentre[centre_ind % L_field_periods].tangent, inverted, 2), 
//			-1.0, multiply2x2(ext_centre[centre_ind % L_field_periods].part_tangent, multiply2x2(gradcentre[(centre_ind-1) % L_field_periods].tangent, inverted, 2), 2), 2 );
//
//			//printmat("ext_centre.part_tangent", ext_centre[centre_ind % L_field_periods].part_tangent, 2, 2);
//		}
//	}
//	// The part below calculates the circumference and the full orbit tangent maps
//	circumference = 0.0;
//	for (centre_ind=0;centre_ind<L_field_periods;centre_ind++) {
//		circumference += sqrt(pow(ext_centre[(centre_ind+1)%L_field_periods].loc[0] 
//					        - ext_centre[centre_ind].loc[0], 2.0)
//				    	    + pow(ext_centre[(centre_ind+1)%L_field_periods].loc[1] 
//						- ext_centre[centre_ind].loc[1], 2.0));
//		ext_centre[centre_ind].epar = malloc(2*sizeof(double)); ext_centre[centre_ind].eperp = malloc(2*sizeof(double));
//		ext_centre[centre_ind].full_tangent = multiply2x2(ext_centre[(centre_ind +2)%L_field_periods].part_tangent, ext_centre[(centre_ind+1)%L_field_periods].part_tangent, 2);
//		grad_ext_centre[centre_ind].full_tangent = add2x2( 1.0, multiply2x2(grad_ext_centre[(centre_ind +2)%L_field_periods].part_tangent, ext_centre[(centre_ind+1)%L_field_periods].part_tangent, 2), 1.0, multiply2x2(ext_centre[(centre_ind +2)%L_field_periods].part_tangent, grad_ext_centre[(centre_ind+1)%L_field_periods].part_tangent, 2), 2);
//
//		for (sec_ind=3;sec_ind<L_field_periods+1;sec_ind++) {
//			multiply2x2reassign(ext_centre[(centre_ind + sec_ind)%L_field_periods].part_tangent, ext_centre[centre_ind].full_tangent, 2);
//			spare_matrix = multiply2x2(grad_ext_centre[(centre_ind + sec_ind)%L_field_periods].part_tangent, ext_centre[centre_ind].full_tangent, 2);
//			multiply2x2reassign(ext_centre[(centre_ind + sec_ind)%L_field_periods].part_tangent, grad_ext_centre[centre_ind].full_tangent, 2);
//			add2x2reassign(1.0, spare_matrix, 1.0, grad_ext_centre[centre_ind].full_tangent, 2);
//			free(spare_matrix[0]); free(spare_matrix[1]);
//			free(spare_matrix);
//		}
//		printmat("ext_centre.fulltangent", ext_centre[centre_ind].full_tangent, 2, 2);
//		symmeigs(ext_centre[centre_ind].full_tangent, ext_centre[centre_ind].eperp, ext_centre[centre_ind].epar, evals);
//		printf("eigenvectors are (%f, %f) and (%f, %f)\n", ext_centre[centre_ind].eperp[0], ext_centre[centre_ind].eperp[1], ext_centre[centre_ind].epar[0], ext_centre[centre_ind].epar[1]);
//	}
//	for (main_ind=0;main_ind<L_field_periods; main_ind++) {
//		deltadet = ext_centre[main_ind].full_tangent[0][0]*grad_ext_centre[main_ind].full_tangent[1][1] + grad_ext_centre[main_ind].full_tangent[1][1]*ext_centre[main_ind].full_tangent[1][1] - ext_centre[main_ind].full_tangent[0][1]*grad_ext_centre[main_ind].full_tangent[1][0] - grad_ext_centre[main_ind].full_tangent[0][1]*ext_centre[main_ind].full_tangent[1][0]; 
//		printf("deltadet=%f\n\n", deltadet);
//		linalg2x2(ext_centre[main_ind].full_tangent, evecs1, evals1, &det, &trace);
//		ext_centre[main_ind].circumference = circumference;
//		if ((fabs(evecs1[0][0]) < small) || (fabs(evecs1[0][1]) < small)) {
//			ext_centre[main_ind].angle = evals1[1];
//			omega = m0_symmetry*evals1[1]/(2.0*M_PI*L_field_periods); // I have changed it compared to Cary and Hanson's paper
//			q0 = (int) (m0_symmetry/(4.0*omega) - L_field_periods/2.0 + 0.5);
//			ext_centre[main_ind].q0_index = q0;
//		}
//		printf("angle=%f\n", ext_centre[main_ind].angle);
//		printf("q0=%d\n", ext_centre[main_ind].q0_index);
//		ext_centre[main_ind].long_tangent = malloc(L_field_periods*sizeof(double));
//		grad_ext_centre[main_ind].long_tangent = malloc(L_field_periods*sizeof(double));
//		for (centre_ind=0;centre_ind<L_field_periods;centre_ind++) {
//			//ext_centre[centre_ind].long_tangent = set_identity();
//			//for (sec_ind=0; sec_ind<q0/L_field_periods; sec_ind++) {
//			//	multiply2x2reassign(ext_centre[0].full_tangent, ext_centre[centre_ind].long_tangent, 2); 
//			//}
//			//for (sec_ind=0; sec_ind<q0%L_field_periods; sec_ind++) {
//			//	multiply2x2reassign(ext_centre[sec_ind+1].part_tangent, ext_centre[centre_ind].long_tangent, 2); 
//			//}
//			ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_field_periods] = set_identity();
//			grad_ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_field_periods] = set_zeros();
//			for (sec_ind=0; sec_ind<(q0+centre_ind)/L_field_periods; sec_ind++) {
//				multiply2x2reassign(ext_centre[main_ind].full_tangent, ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_field_periods], 2); 
//				spare_matrix = multiply2x2(grad_ext_centre[main_ind].full_tangent, ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_field_periods], 2);
//				multiply2x2reassign(ext_centre[main_ind].full_tangent, grad_ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_field_periods], 2);
//				add2x2reassign(1.0, spare_matrix, 1.0, grad_ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_field_periods], 2);
//				free(spare_matrix[0]); free(spare_matrix[1]);
//				free(spare_matrix);
//			}
//			for (sec_ind=0; sec_ind<(q0+centre_ind)%L_field_periods; sec_ind++) {
//				multiply2x2reassign(ext_centre[(sec_ind+1+main_ind)%L_field_periods].part_tangent, ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_field_periods], 2); 
//				spare_matrix = multiply2x2(grad_ext_centre[(sec_ind+1+main_ind)%L_field_periods].part_tangent, ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_field_periods], 2);
//				multiply2x2reassign(ext_centre[(sec_ind+1+main_ind)%L_field_periods].part_tangent, grad_ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_field_periods], 2);
//				add2x2reassign(1.0, spare_matrix, 1.0, grad_ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_field_periods], 2);
//				free(spare_matrix[0]); free(spare_matrix[1]);
//				free(spare_matrix);
//			}
//		}
//	}
//	//ext_centre[centre_index % L_field_periods].part_tangent = multiply2x2(centre[centre_index % L_field_periods].tangent, inverted, 2);
//	//clock_t int3 = clock();
//	return grad_ext_centre;
//} 

// NEW FUNCTION because previous one is incorrect
struct ext_position *gradalongcentrealt(double RR, double ZZ, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode, double ***coils, int *n_coils, int **n_segs) {
	int i=0, q0, L_field_periods;
	struct position fieldline_start, fieldline, gradfieldline, *centre, *gradcentre;
	double evals[2], **gradsigmaM, **gradsymm, **sigma;
	double varphi=0.0, omega, actualdet, graddet;
	double deltaepar, deltaeperp;
	double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
	double **inverted, detcentre, **evecs1, *evals1, det=0, trace=0, **spare_matrix;
	double circumference;
	int N_line=0, main_ind, centre_ind, sec_ind;
	struct ext_position *ext_centre, *grad_ext_centre;

	evecs1 = malloc(2*sizeof(double*));
	evecs1[0] = malloc(2*sizeof(double)); evecs1[1] = malloc(2*sizeof(double));
	evals1 = malloc(2*sizeof(double));

	if (tor_mode % m0_symmetry == 0) 
		L_field_periods = pol_mode;	
	else 				   
		L_field_periods = m0_symmetry*pol_mode;
	N_line = L_field_periods*N_gridphi_per_field_period;
	centre = malloc(L_field_periods*sizeof(struct position));
	gradcentre = malloc(L_field_periods*sizeof(struct position));
	ext_centre = malloc(L_field_periods*sizeof(struct ext_position));
	grad_ext_centre = malloc(L_field_periods*sizeof(struct ext_position));
	fieldline_start.loc[0] = RR; fieldline_start.loc[1] = ZZ;
	fieldline_start.tangent = set_identity();
	fieldline = fieldline_start;
	gradfieldline.loc[0] = 0.0; gradfieldline.loc[1] = 0.0;
	gradfieldline.tangent = set_zeros();
	centre[0] = fieldline;
	gradcentre[0] = gradfieldline;
	varphi = 0.0;
	printf("varphi = %f\n", varphi);
	printstruct("fieldline[i]\n", &fieldline);
	for (i=1; i<N_line+1; i++) {
		//RK4(&fieldline, varphi, dvarphi, coils, n_coils, n_segs);
		RK4_wgrad(&fieldline, &gradfieldline, varphi, dvarphi, coils, n_coils, n_segs);
		varphi += dvarphi;
		if (i%N_gridphi_per_field_period==0) {
			centre_ind = (i / N_gridphi_per_field_period);
			centre[centre_ind % L_field_periods] = fieldline;
			gradcentre[centre_ind % L_field_periods] = gradfieldline;
			printf("varphi = %f\n", varphi);
			printstruct("fieldline[i]\n", &fieldline);
			inverted = invert2x2(centre[(centre_ind-1) % L_field_periods].tangent, &detcentre);
			printmat("inverted", inverted, 2, 2);
			//ptarray[0] = &centre[centre_ind % L_field_periods].tangent[0][0];
			//ptarray[1] = &centre[centre_ind % L_field_periods].tangent[1][0];
			ext_centre[centre_ind % L_field_periods].loc[0] = fieldline.loc[0];
			ext_centre[centre_ind % L_field_periods].loc[1] = fieldline.loc[1];
			ext_centre[centre_ind % L_field_periods].part_tangent = multiply2x2(centre[centre_ind % L_field_periods].tangent, inverted, 2);
			//printmat("ext_centre.part_tangent", ext_centre[centre_ind % L_field_periods].part_tangent, 2, 2);
			//printmat("ext_centre.part_tangent", ext_centre[centre_ind % L_field_periods].part_tangent, 2, 2);
		}
	}
	// The part below calculates the circumference and the full orbit tangent maps
	circumference = 0.0;
	sigma = malloc(2*sizeof(double)); sigma[0] = malloc(2*sizeof(double)); sigma[1] = malloc(2*sizeof(double));
	sigma[0][0] = 0.0;
	sigma[0][1] = -1.0;
	sigma[1][0] = 1.0;
	sigma[1][1] = 0.0;
	//gradsymm = malloc(2*sizeof(double)); gradsymm[0] = malloc(2*sizeof(double)); gradsymm[1] = malloc(2*sizeof(double));
	for (main_ind=0;main_ind<L_field_periods;main_ind++) {
		circumference += sqrt(pow(ext_centre[(main_ind+1)%L_field_periods].loc[0] 
					        - ext_centre[main_ind].loc[0], 2.0)
				    	    + pow(ext_centre[(main_ind+1)%L_field_periods].loc[1] 
						- ext_centre[main_ind].loc[1], 2.0));
		ext_centre[main_ind].epar = malloc(2*sizeof(double)); ext_centre[main_ind].eperp = malloc(2*sizeof(double));
		//grad_ext_centre[main_ind].epar = malloc(2*sizeof(double)); grad_ext_centre[main_ind].eperp = malloc(2*sizeof(double));
		ext_centre[main_ind].full_tangent = multiply2x2(ext_centre[(main_ind +2)%L_field_periods].part_tangent, ext_centre[(main_ind+1)%L_field_periods].part_tangent, 2);

		fieldline = centre[main_ind];
		fieldline.tangent[0][0] = 1.0;
		fieldline.tangent[0][1] = 0.0;
		fieldline.tangent[1][0] = 0.0;
		fieldline.tangent[1][1] = 1.0;
		gradfieldline = gradcentre[main_ind];
		gradfieldline.loc[0] = 0.0; gradfieldline.loc[1] = 0.0;
		gradfieldline.tangent[0][0] = 0.0;
		gradfieldline.tangent[0][1] = 0.0;
		gradfieldline.tangent[1][0] = 0.0;
		gradfieldline.tangent[1][1] = 0.0;
		varphi = main_ind*2.0*M_PI/m0_symmetry;
		for (i=0; i<N_gridphi_per_field_period*L_field_periods; i++) {
			centre_ind = (i / N_gridphi_per_field_period);
			RK4_wgrad(&fieldline, &gradfieldline, varphi, dvarphi, coils, n_coils, n_segs);
			varphi += dvarphi;
		}
		ext_centre[main_ind].full_tangent = fieldline.tangent;
		grad_ext_centre[main_ind].full_tangent = gradfieldline.tangent;
		printmat("ext_centre.fulltangent", ext_centre[main_ind].full_tangent, 2, 2);
		printmat("grad_ext_centre.fulltangent1", grad_ext_centre[main_ind].full_tangent, 2, 2);
		symmeigs(ext_centre[main_ind].full_tangent, ext_centre[main_ind].eperp, ext_centre[main_ind].epar, evals);
		printmat("grad_ext_centre.fulltangent2", grad_ext_centre[main_ind].full_tangent, 2, 2);
		graddet = ext_centre[main_ind].full_tangent[0][0]*grad_ext_centre[main_ind].full_tangent[1][1] + grad_ext_centre[main_ind].full_tangent[1][1]*ext_centre[main_ind].full_tangent[1][1] - ext_centre[main_ind].full_tangent[0][1]*grad_ext_centre[main_ind].full_tangent[1][0] - grad_ext_centre[main_ind].full_tangent[0][1]*ext_centre[main_ind].full_tangent[1][0]; 
		printf("graddet=%f\n\n", graddet);
		actualdet = ext_centre[main_ind].full_tangent[0][0]*ext_centre[main_ind].full_tangent[1][1] -  ext_centre[main_ind].full_tangent[0][1]*ext_centre[main_ind].full_tangent[1][0]; 
		printf("det=%f\n\n", actualdet);
		//gradsigmaM = multiply2x2(sigma, grad_ext_centre[main_ind].full_tangent, 2);
		//gradsymm[0][0] = gradsigmaM[0][0];
		//gradsymm[1][1] = gradsigmaM[1][1];
		//gradsymm[0][1] = 0.5*gradsigmaM[1][0] + 0.5*gradsigmaM[0][1];
		//gradsymm[1][0] = gradsymm[0][1];
		
		//deltaepar = inner(ext_centre[main_ind].eperp,  gradsymm, ext_centre[main_ind].epar)/(evals[0]-evals[1]);
		//deltaeperp = - inner(ext_centre[main_ind].epar,  gradsymm, ext_centre[main_ind].eperp)/(evals[0]-evals[1]);
		//grad_ext_centre[main_ind].epar[0] = deltaepar*ext_centre[main_ind].eperp[0]; grad_ext_centre[main_ind].epar[1] = deltaepar*ext_centre[main_ind].eperp[1];
		//grad_ext_centre[main_ind].eperp[0] = deltaeperp*ext_centre[main_ind].epar[0]; grad_ext_centre[main_ind].eperp[1] = deltaeperp*ext_centre[main_ind].epar[1];
		//free(gradsigmaM);
		//grad_ext_centre[main_ind].eperp[0] 
		printmat("grad_ext_centre.fulltangent3", grad_ext_centre[main_ind].full_tangent, 2, 2);
	}
	for (main_ind=0;main_ind<L_field_periods; main_ind++) {
		printmat("grad_ext_centre.fulltangent4", grad_ext_centre[main_ind].full_tangent, 2, 2);
		linalg2x2(ext_centre[main_ind].full_tangent, evecs1, evals1, &det, &trace);
		ext_centre[main_ind].circumference = circumference;
		if ((fabs(evecs1[0][0]) < small) || (fabs(evecs1[0][1]) < small)) {
			ext_centre[main_ind].angle = evals1[1];
			omega = m0_symmetry*evals1[1]/(2.0*M_PI*L_field_periods); // I have changed it compared to Cary and Hanson's paper
			q0 = (int) (m0_symmetry/(4.0*omega) - L_field_periods/2.0 + 0.5);
			ext_centre[main_ind].q0_index = q0;
			grad_ext_centre[main_ind].angle = -0.5*(grad_ext_centre[main_ind].full_tangent[0][0]+grad_ext_centre[main_ind].full_tangent[1][1])/sin(ext_centre[main_ind].angle);
		}
		printmat("grad_ext_centre.fulltangent", grad_ext_centre[main_ind].full_tangent, 2, 2);
		//printf("gradangle pieces are %f %f\n", grad_ext_centre[main_ind].full_tangent[0][0]+grad_ext_centre[main_ind].full_tangent[1][1]);
		printf("gradangle=%f\n", grad_ext_centre[main_ind].angle);
		printf("angle=%f\n", ext_centre[main_ind].angle);
		printf("q0=%d\n", ext_centre[main_ind].q0_index);
		varphi = 0.0; //modify this 
		ext_centre[main_ind].long_tangent = malloc(L_field_periods*sizeof(double));
		grad_ext_centre[main_ind].long_tangent = malloc(L_field_periods*sizeof(double));
		fieldline = centre[main_ind];
		fieldline.tangent[0][0] = 1.0;
		fieldline.tangent[0][1] = 0.0;
		fieldline.tangent[1][0] = 0.0;
		fieldline.tangent[1][1] = 1.0;
		gradfieldline.loc[0] = 0.0; gradfieldline.loc[1] = 0.0;
		gradfieldline.tangent[0][0] = 0.0;
		gradfieldline.tangent[0][1] = 0.0;
		gradfieldline.tangent[1][0] = 0.0;
		gradfieldline.tangent[1][1] = 0.0;
		for (i=0; i<N_gridphi_per_field_period*(q0+L_field_periods); i++) {
			centre_ind = (i / N_gridphi_per_field_period);
			if ( (i%N_gridphi_per_field_period==0) && (centre_ind/q0 == 1) ) {
				ext_centre[main_ind].long_tangent[(main_ind+centre_ind)%L_field_periods] = fieldline.tangent;
				grad_ext_centre[main_ind].long_tangent[(main_ind+centre_ind)%L_field_periods] = gradfieldline.tangent;
				printf("index = %d\n", centre_ind%q0);
			}
			RK4_wgrad(&fieldline, &gradfieldline, varphi, dvarphi, coils, n_coils, n_segs);
			varphi += dvarphi;
		}
	}
	//ext_centre[centre_index % L_field_periods].part_tangent = multiply2x2(centre[centre_index % L_field_periods].tangent, inverted, 2);
	//clock_t int3 = clock();
	return grad_ext_centre;
} 

double *gradislandwidth(struct ext_position *ext_fieldline, struct ext_position *grad_ext_fieldline, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode) {
	int main_index, centre_index, L_field_periods;
	double *wperp, *gradwperp, sum_matrix_elements, grad_sum_matrix_elements, matrix_element, grad_matrix_element;
	double circumference=ext_fieldline[0].circumference;
	if (tor_mode % m0_symmetry == 0) L_field_periods = pol_mode;	
	else 				   L_field_periods = m0_symmetry*pol_mode;
	wperp = malloc(L_field_periods*sizeof(double));
	gradwperp = malloc(L_field_periods*sizeof(double));
	for (main_index=0;main_index<L_field_periods; main_index++) {
		//circumference = 0.0;		
		sum_matrix_elements = 0.0;		
		grad_sum_matrix_elements = 0.0;		
		for (centre_index=0;centre_index<L_field_periods; centre_index++) {
			matrix_element = inner(ext_fieldline[centre_index].epar, 
					       ext_fieldline[main_index].long_tangent[centre_index], 
					       ext_fieldline[main_index].eperp);	
			grad_matrix_element = ( inner(ext_fieldline[centre_index].epar, grad_ext_fieldline[main_index].long_tangent[centre_index], ext_fieldline[main_index].eperp)/sin( ext_fieldline[main_index].angle*(ext_fieldline[main_index].q0_index+centre_index)/L_field_periods) 
			- cos( ext_fieldline[main_index].angle*(ext_fieldline[main_index].q0_index+centre_index)/L_field_periods ) * ( grad_ext_fieldline[main_index].angle*(ext_fieldline[main_index].q0_index+centre_index)/L_field_periods) * inner( ext_fieldline[centre_index].epar, ext_fieldline[main_index].long_tangent[centre_index], ext_fieldline[main_index].eperp ) / pow(sin( ext_fieldline[main_index].angle*(ext_fieldline[main_index].q0_index+centre_index)/L_field_periods), 2.0 ) ) * matrix_element/fabs(matrix_element);
			printf("%f %f\n", cos( ext_fieldline[main_index].angle*(ext_fieldline[main_index].q0_index+centre_index)/L_field_periods) ,  (grad_ext_fieldline[main_index].angle*(ext_fieldline[main_index].q0_index+centre_index)/L_field_periods)) ;
			printf("%f\n", grad_ext_fieldline[main_index].angle);
				//		+inner(ext_fieldline[centre_index].epar, 
				//	       ext_fieldline[main_index].long_tangent[centre_index], 
				//	       grad_ext_fieldline[main_index].eperp)
				//		+inner(grad_ext_fieldline[centre_index].epar, 
				//	       ext_fieldline[main_index].long_tangent[centre_index], 
				//	       ext_fieldline[main_index].eperp))*matrix_element/fabs(matrix_element);	
			sum_matrix_elements += fabs(matrix_element);
			grad_sum_matrix_elements += grad_matrix_element;
			//printmat("symmetrized invariant matrix", ext_fieldline[main_index].long_tangent[centre_index], 2, 2); 
			printmat("symmetrized invariant matrix", grad_ext_fieldline[main_index].long_tangent[centre_index], 2, 2); 
			printf("%f\n", inner(ext_fieldline[centre_index].epar, 
					       grad_ext_fieldline[main_index].long_tangent[centre_index], 
					       ext_fieldline[main_index].eperp));	
			//printf("%f\n", inner(grad_ext_fieldline[centre_index].epar, 
			//		       ext_fieldline[main_index].long_tangent[centre_index], 
			//		       ext_fieldline[main_index].eperp));	
			//printf("%f\n", inner(ext_fieldline[centre_index].epar, 
			//		       ext_fieldline[main_index].long_tangent[centre_index], 
			//		       grad_ext_fieldline[main_index].eperp));	
			//sum_matrix_elements += fabs(matrix_element)/sin(ext_fieldline[main_index].angle*ext_fieldline[main_index].q0_index/L_field_periods); //This has the sin piece in the denominator, which can be safely set to unity
			//printf("cosine = %f\n", cos(ext_fieldline[main_index].angle*ext_fieldline[main_index].q0_index/L_field_periods));
			//circumference += sqrt(pow(ext_fieldline[(centre_index+1)%L_field_periods].loc[0] 
			//		        - ext_fieldline[centre_index].loc[0], 2.0)
			//	    	    + pow(ext_fieldline[(centre_index+1)%L_field_periods].loc[1] 
			//			- ext_fieldline[centre_index].loc[1], 2.0));
			//printf("circumference = %f\nmatrix_element = %f\n", circumference, matrix_element);
			//printf("other circumference = %f\n", ext_fieldline[centre_index].circumference);
		}
		//circumference = 2.0*M_PI*0.2097;
		printf("circumference = %f\n", circumference);
		wperp[main_index] = 2.0*L_field_periods*circumference/(M_PI*pol_mode*sum_matrix_elements);
		gradwperp[main_index] = - ( grad_sum_matrix_elements / sum_matrix_elements ) * wperp[main_index];
		printf("width = %.10f for index= %d\n", wperp[main_index], main_index);
		printf("gradwidth = %f for index= %d\n", gradwperp[main_index], main_index);
	}
	return gradwperp;
}
