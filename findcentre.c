/* Author: Alessandro Geraldini */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_fit.h>
#include "isc.h"
#define smallerror 1.0e-13
#define largeerror 1.0e-13
#define FACTOR 1.0
#define MAX_COUNT 20

/*Jargon
Fixed point: intersection of the periodic field line with a poloidal symmetry plane
*/

// This function returns the poloidally reordered index of the fixed point  
int rho(int fieldaligned_index, int n_turns, int L_fixedpoints) {
	int reordered_index;
	reordered_index = (n_turns*fieldaligned_index)%L_fixedpoints;
	return reordered_index;
}

// This function returns the index of the fixed point from the poloidally reordered one 
int inverserho(int reordered_index, int n_turns, int L_fixedpoints) {
	int fieldaligned_index;
	//printf("reordered_index=%d, n_turns=%d, L_fixedpoints=%d\n", reordered_index, n_turns, L_fixedpoints);
	fieldaligned_index = ( ( (reordered_index%n_turns)*L_fixedpoints + reordered_index ) / n_turns ) % L_fixedpoints;
	return fieldaligned_index;
}

// Print structures with a column 2-vector and a 2x2 tangent map
void printstructposition(char *name, struct position *input) {
	printf("structure position %s:\nposition = (%10.8f, %10.8f)\ntangent = |%10.8f %10.8f|\n          |%10.8f %10.8f|\n", 
	 name, input->loc[0], input->loc[1], input->tangent[0][0], input->tangent[0][1], input->tangent[1][0], input->tangent[1][1]);
	//printf("for structure position %s:\nposition = (%f, %f)\n", name, input->loc[0], input->loc[1]);
	//printmat("tangent", input->tangent, 2, 2);
}

// Print structures with a column 3-vector, a 3x2 matrix and a 3x2x2 tensor (B-field value, first and second derivatives) 
void printstructfield(char *name, struct field *input) {
	printf("structure field %s:\nlocation = (%10.10f, %10.10f, %10.10f)\nderivative = |%10.10f %10.10f |\n          |%10.10f %10.10f |\n          |%10.10f %10.10f |", 
	 name, input->value[0], input->value[1], input->value[2], input->derivative[0][0], input->derivative[0][1], input->derivative[1][0], input->derivative[1][1], input->derivative[2][0], input->derivative[2][1]);
}

// Find the magnetic axis 
struct position *solve_magneticaxis(struct fieldparams allparams, struct field **Bfieldsaved, struct position *fieldline, int N_gridphi_fieldperiod) {
	int count = 0, i=0, N_lines = N_gridphi_fieldperiod, m0_fieldperiods, alternate=1;
	struct position fieldline_start, deltafieldline;
	double varphi=0.0, *jumptocentre;
	double **evecs1=malloc(2*sizeof(double*)), *evals1=malloc(2*sizeof(double)), det, trace;
	double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*allparams.m0_fieldperiods);
	//double **matrix;
	double det_tangent;
	double **inverseTminusI; // **inverseT;
	double error=1.0, factor=FACTOR, olderror = 99999999999.9;
	struct position *centre=malloc(N_lines*sizeof(struct position));
	evecs1[0] = malloc(2*sizeof(double)); evecs1[1] = malloc(2*sizeof(double));
	m0_fieldperiods = allparams.m0_fieldperiods;
	fieldline_start.tangent = set_identity();
	fieldline_start.loc[0] = fieldline->loc[0];
	fieldline_start.loc[1] = fieldline->loc[1];
	//printf("Finding the magnetic axis position at varphi = 0 and symmetry planes 2*pi*k/m_0 for integer k:\n");
	do {
		fieldline->tangent[0][0]=1.0; fieldline->tangent[0][1]=0.0; 
		fieldline->tangent[1][0]=0.0; fieldline->tangent[1][1]=1.0;
		varphi = 0.0;
		printf("iteration number %d = (%14.10f, %14.10f)\n", count, fieldline->loc[0], fieldline->loc[1]);
		for (i=0; i<N_lines; i++)
		{
			//printf("i=%d/%d\n", i, N_lines);
			centre[i] = *fieldline;
			//printstructposition("fieldline", fieldline);
			RK4step(fieldline, varphi, dvarphi, allparams, Bfieldsaved[i]);
			varphi += dvarphi;
			//printstructposition("fieldline", fieldline);
			//printstructfield("field", Bfieldsaved[i]);
			//printstructfield("field", Bfieldsaved[i]+1);
			//printstructfield("field", Bfieldsaved[i]+2);
			//printstructfield("field", Bfieldsaved[i]+3);
		}
		deltafieldline = addstructs(-1.0, fieldline, 1.0, &fieldline_start); 
		inverseTminusI = invert2x2(deltafieldline.tangent, &det_tangent);
		jumptocentre      = multiply(inverseTminusI, &deltafieldline.loc[0]);
		jumptocentre[0]  *= factor;
		jumptocentre[1]  *= factor;
		if (count == 0) {
			linalg2x2(fieldline->tangent, evecs1, evals1, &det, &trace);
			//q0plushalfL = (M_PI/2.0)/evecs1;
			alternate = (int) (M_PI/2.0)/evals1[1];
			alternate = 1;
		}
		if (count % alternate == 0 || error < largeerror) {
			fieldline->loc[0] = fieldline_start.loc[0] - jumptocentre[0]; 
			fieldline->loc[1] = fieldline_start.loc[1] - jumptocentre[1];
			//fieldline->tangent[0][0]=1.0; fieldline->tangent[0][1]=0.0; 
			//fieldline->tangent[1][0]=0.0; fieldline->tangent[1][1]=1.0;
			fieldline_start.loc[0] = fieldline->loc[0];
			fieldline_start.loc[1] = fieldline->loc[1];
			fieldline_start.tangent[0][0] = 1.0; fieldline_start.tangent[0][1] = 0.0;
			fieldline_start.tangent[1][0] = 0.0; fieldline_start.tangent[1][1] = 1.0;
		}
		error = sqrt( (deltafieldline.loc[0]*deltafieldline.loc[0] + deltafieldline.loc[1]*deltafieldline.loc[1]) / (fieldline->loc[0]*fieldline->loc[0] + fieldline->loc[1]*fieldline->loc[1]) ); 
		free(jumptocentre);
		free(inverseTminusI[0]); free(inverseTminusI[1]);
		free(inverseTminusI);
		i = 1;
		//printf("error = %15.15f\n", error);
		//exit(0);
		//factor*=0.99;
		count+=1;
		if (error > olderror) factor *= (olderror/error);
		else factor = FACTOR;
		olderror = error;
	} while(error>smallerror && count < MAX_COUNT);
	if (count == MAX_COUNT) printf("Warning: axis not found\n");
	return centre;
}

void solve_magneticaxis_save(double *axisguess, struct position **centre, struct fieldparams allparams, struct field **Bfieldsaved, int N_gridphi_fieldperiod) {
	// declarations
	//clock_t start = clock();
	int count = 0, i=0, N_lines = N_gridphi_fieldperiod, m0_fieldperiods, alternate=1;
	struct position fieldline_start, deltafieldline, fieldline;
	double varphi=0.0, *jumptocentre;
	double **evecs1=malloc(2*sizeof(double*)), *evals1=malloc(2*sizeof(double)), det, trace;
	double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*allparams.m0_fieldperiods);
	//double **matrix;
	double det_tangent;
	double **inverseTminusI; // **inverseT;
	double error=1.0, factor=FACTOR, olderror=9999999999.9;
	evecs1[0] = malloc(2*sizeof(double)); evecs1[1] = malloc(2*sizeof(double));
	m0_fieldperiods = allparams.m0_fieldperiods;
	//type, m0_fieldperiods, coils, n_coils, n_segs, constparams
	//printf("-->Entered\n");
	//printf("m0_fieldperiods=%d\n", m0_fieldperiods);
	fieldline_start.tangent = set_identity();
	fieldline.tangent = set_identity();
	fieldline.loc[0] = fieldline_start.loc[0] = axisguess[0];
	fieldline.loc[1] = fieldline_start.loc[1] = axisguess[1];
	//printf("Finding the magnetic axis position at varphi = 0 and symmetry planes 2*pi*k/m_0 for integer k:\n");
	do {
		fieldline.tangent[0][0]=1.0; fieldline.tangent[0][1]=0.0; 
		fieldline.tangent[1][0]=0.0; fieldline.tangent[1][1]=1.0;
		varphi = 0.0;
		printf("iteration number %d = (%14.10f, %14.10f)\n", count, fieldline.loc[0], fieldline.loc[1]);
		for (i=0; i<N_lines; i++)
		{
			//printf("i=%d/%d\n", i, N_lines);
			centre[i][0] = fieldline;
			RK4stepsave(&fieldline, centre[i], varphi, dvarphi, allparams, Bfieldsaved[i]);
			varphi += dvarphi;
			//printstructposition("fieldline", fieldline);
		}
		deltafieldline = addstructs(-1.0, &fieldline, 1.0, &fieldline_start); 
		inverseTminusI = invert2x2(deltafieldline.tangent, &det_tangent);
		jumptocentre      = multiply(inverseTminusI, &deltafieldline.loc[0]);
		jumptocentre[0]  *= factor;
		jumptocentre[1]  *= factor;
		if (count == 0) {
			linalg2x2(fieldline.tangent, evecs1, evals1, &det, &trace);
			//q0plushalfL = (M_PI/2.0)/evecs1;
			alternate = (int) (M_PI/2.0)/evals1[1];
			alternate = 1;
		}
		if (count % alternate == 0 || error < largeerror) {
			fieldline.loc[0] = fieldline_start.loc[0] - jumptocentre[0]; 
			fieldline.loc[1] = fieldline_start.loc[1] - jumptocentre[1];
			//fieldline->tangent[0][0]=1.0; fieldline->tangent[0][1]=0.0; 
			//fieldline->tangent[1][0]=0.0; fieldline->tangent[1][1]=1.0;
			fieldline_start.loc[0] = fieldline.loc[0];
			fieldline_start.loc[1] = fieldline.loc[1];
			fieldline_start.tangent[0][0] = 1.0; fieldline_start.tangent[0][1] = 0.0;
			fieldline_start.tangent[1][0] = 0.0; fieldline_start.tangent[1][1] = 1.0;
		}
		error = sqrt( (deltafieldline.loc[0]*deltafieldline.loc[0] + deltafieldline.loc[1]*deltafieldline.loc[1]) / (fieldline.loc[0]*fieldline.loc[0] + fieldline.loc[1]*fieldline.loc[1]) ); 
		free(jumptocentre);
		free(inverseTminusI[0]); free(inverseTminusI[1]);
		free(inverseTminusI);
		i = 1;
		//printf("error = %15.15f\n", error);
		//exit(0);
		//factor*=0.99;
		count+=1;
		if (error > olderror) factor*= (olderror/error);
		else factor = FACTOR;
		olderror = error;
	} while(error>smallerror && count < MAX_COUNT);
	if (count == MAX_COUNT) printf("Warning: axis not found\n");
	//} while(error>smallerror);
	//printf("Leaving<--\n");
	//clock_t int3 = clock();
	return;
}

//struct position *solve_islandcenter(struct fieldparams allparams, struct field **Bfieldsaved, struct position *fieldline, int L_fixedpoints, int N_gridphi_fieldperiod) {
//	// declarations
//	//clock_t start = clock();
//	int i=0, N_line, count = 0, alternate = 1;
//	double **evecs1=malloc(2*sizeof(double*)), *evals1=malloc(2*sizeof(double)), det, trace;
//	struct position fieldline_start, deltafieldline;
//	double varphi=0.0, *jumptocentre;
//	double dvarphi = 2.0*M_PI/(allparams.m0_fieldperiods*N_gridphi_fieldperiod);
//	//double **matrix;
//	double **pdeltafieldline, det_tangent;
//	double **inverseTminusI, **inverseT;
//	double error = 1.0, olderror= 9999999999.9;
//	double factor = FACTOR;
//	struct position *centre;
//	evecs1[0] = malloc(2*sizeof(double)); evecs1[1] = malloc(2*sizeof(double));
//	pdeltafieldline = calloc(2,sizeof(double*)); 
//	N_line = L_fixedpoints*N_gridphi_fieldperiod;
//	centre = malloc(N_line*sizeof(struct position));
//	fieldline_start = *fieldline;
//	printf("Finding the island centre position at varphi = 0:\n");
//	do {
//		varphi = 0.0;
//		printf("iteration #%d = (%14.10f, %14.10f)\n", count, fieldline->loc[0], fieldline->loc[1]);
//		for (i=0; i<N_line; i++)
//		{
//			centre[i] = *fieldline;
//			RK4step(fieldline, varphi, dvarphi, allparams, Bfieldsaved[i]);
//			//printf("B field saved value (BR) = %f\n", Bfieldsaved[i][0].value[0]);
//			varphi += dvarphi;
//		}
//		deltafieldline = addstructs(1.0, fieldline, -1.0, &fieldline_start); 
//		inverseTminusI = invert2x2(deltafieldline.tangent, &det_tangent);
//		inverseT = invert2x2(fieldline->tangent, &det_tangent);
//		jumptocentre = multiply(inverseTminusI, &deltafieldline.loc[0]);
//		jumptocentre[0] *= (factor);
//		jumptocentre[1] *= (factor);
//		if (count == 0) {
//			linalg2x2(fieldline->tangent, evecs1, evals1, &det, &trace);
//			//q0plushalfL = (M_PI/2.0)/evecs1;
//			alternate = (int) (M_PI/2.0)/evals1[1];
//			alternate = 1;
//		}
//		if (count%alternate == 0 || error < largeerror) {
//			fieldline->loc[0] = fieldline_start.loc[0] - jumptocentre[0]; 
//			fieldline->loc[1] = fieldline_start.loc[1] - jumptocentre[1];
//			fieldline->tangent[0][0]=1.0; fieldline->tangent[0][1]=0.0; 
//			fieldline->tangent[1][0]=0.0; fieldline->tangent[1][1]=1.0;
//			fieldline_start = *fieldline;
//		}
//		error = sqrt( (deltafieldline.loc[0]*deltafieldline.loc[0] + deltafieldline.loc[1]*deltafieldline.loc[1] ) / (fieldline->loc[0]*fieldline->loc[0] + fieldline->loc[1]*fieldline->loc[1]) ); 
//		//printf("and error = %15.15f\n", error);
//		free(jumptocentre);
//		//factor*=0.99;
//		if (error > olderror) factor*= (olderror/error);
//		else factor = FACTOR;
//		olderror = error;
//		count++;
//	//} while (error > smallerror);
//	} while(error>smallerror && count < MAX_COUNT);
//	if (count == MAX_COUNT) printf("Warning: periodic field line not found\n");
//	return centre;
//}

void solve_islandcenter_save(struct position **centre, double *Res, struct fieldparams allparams, struct field **Bfieldsaved, int L_fixedpoints, int N_gridphi_fieldperiod) {
	// declarations
	//clock_t start = clock();
	int i=0, N_line, count = 0, alternate = 1;
	double **evecs1=malloc(2*sizeof(double*)), *evals1=malloc(2*sizeof(double)), det, trace;
	struct position fieldline, fieldline_start, deltafieldline;
	double varphi=0.0, *jumptocentre;
	double dvarphi = 2.0*M_PI/(allparams.m0_fieldperiods*N_gridphi_fieldperiod);
	//double **matrix;
	double **pdeltafieldline, det_tangent;
	double **inverseTminusI, **inverseT;
	double error = 1.0, absolerror=1.0, olderror = 9999999999.9;
	double factor = FACTOR;
	evecs1[0] = malloc(2*sizeof(double)); evecs1[1] = malloc(2*sizeof(double));
	pdeltafieldline = calloc(2,sizeof(double*)); 
	N_line = L_fixedpoints*N_gridphi_fieldperiod;
	//centre = malloc(N_line*sizeof(struct position));
	//for (i=0;i<N_line;i++) centre[i] = malloc(4*sizeof(struct position));
	//centre = malloc(N_line*sizeof(struct position));
	//for (i=0; i<N_line; i++) centre[i] = malloc(4*sizeof(struct position));
	fieldline.tangent = set_identity();
	fieldline_start.tangent = set_identity();
	fieldline_start.loc[0] = fieldline.loc[0] = centre[0][0].loc[0];
	fieldline_start.loc[1] = fieldline.loc[1] = centre[0][0].loc[1];
	//printf("Finding the island centre position at varphi = 0:\n");
	//if (axis == NULL) {
		do {
			varphi = 0.0;
			printf("iteration #%d = (%14.10f, %14.10f)\n", count, fieldline.loc[0], fieldline.loc[1]);
			for (i=0; i<N_line; i++)
			{
				//centre[i] = *fieldline;
				centre[i][0] = fieldline;
				RK4stepsave(&fieldline, centre[i], varphi, dvarphi, allparams, Bfieldsaved[i]);
				varphi += dvarphi;
			}
			centre[0][0].tangent[0][0] = fieldline.tangent[0][0];
			centre[0][0].tangent[1][0] = fieldline.tangent[1][0];
			centre[0][0].tangent[0][1] = fieldline.tangent[0][1];
			centre[0][0].tangent[1][1] = fieldline.tangent[1][1];
			*Res = 0.5 - 0.25*(fieldline.tangent[0][0] + fieldline.tangent[1][1]); 
			deltafieldline = addstructs(1.0, &fieldline, -1.0, &fieldline_start); 
			//printstructposition("fieldline", &fieldline); printstructposition("fieldline_start", &fieldline_start);
			//printstructposition("deltafieldline", &deltafieldline); 
			inverseTminusI = invert2x2(deltafieldline.tangent, &det_tangent);
			//printmat("inverseTminusI", inverseTminusI, 2, 2);
			//inverseT = invert2x2(fieldline.tangent, &det_tangent);
			jumptocentre = multiply(inverseTminusI, &deltafieldline.loc[0]);
			//printf("jumptocentre = (%f, %f)\n", jumptocentre[0], jumptocentre[1]);
			jumptocentre[0] *= (factor);
			jumptocentre[1] *= (factor);
			if (count == 0) {
				linalg2x2(fieldline.tangent, evecs1, evals1, &det, &trace);
				//q0plushalfL = (M_PI/2.0)/evecs1;
				alternate = (int) (M_PI/2.0)/evals1[1];
				alternate = 1;
			}
			if (count%alternate == 0 || error < largeerror) {
				fieldline.loc[0] = fieldline_start.loc[0] = fieldline_start.loc[0] - jumptocentre[0]; 
				fieldline.loc[1] = fieldline_start.loc[1] = fieldline_start.loc[1] - jumptocentre[1];
				fieldline.tangent[0][0] = fieldline_start.tangent[0][0] = 1.0; 
				fieldline.tangent[0][1] = fieldline_start.tangent[0][1] = 0.0; 
				fieldline.tangent[1][0] = fieldline_start.tangent[1][0] = 0.0; 
				fieldline.tangent[1][1] = fieldline_start.tangent[1][1] = 1.0;
			}
			absolerror = (deltafieldline.loc[0]*deltafieldline.loc[0] + deltafieldline.loc[1]*deltafieldline.loc[1] );
			error = sqrt( absolerror/ (fieldline.loc[0]*fieldline.loc[0] + fieldline.loc[1]*fieldline.loc[1]) ); 
			//printf("error = %15.15f\n", error);
			//printf("absolerror = %15.15f\n", absolerror);
			printf("%f\n", fieldline.loc[0]*fieldline.loc[0] + fieldline.loc[1]*fieldline.loc[1]);
			free(jumptocentre);
			free(inverseTminusI[0]); free(inverseTminusI[1]);
			free(inverseTminusI);
			if (error > olderror) factor*= (olderror/error);
			else factor = FACTOR;
			olderror = error;
			count++;
		//} while (error > smallerror);//[i%N_gridphi_fieldperiods]
		} while(error>smallerror && count < MAX_COUNT);
		if (count == MAX_COUNT) printf("Warning: periodic field line not found\n");
	//}
	//else {

	//}
	return;
}

void extsolve_periodicfieldline(struct ext_position *ext_centre, struct position *fieldline, struct position *axis, int m0_symmetry, int L_fixedpoints, int pol_mode, int *q0_index, int N_gridphi_fieldperiod) {
	int i=0, q0, q0ver, kminus, kplus, rhominus, rhoplus, number_turns;
	struct position fieldline_start, *adjfieldline;
	double varphi=0.0, omega, rawangle, rawangle_principal, *angle_axis;
	//double c0, c1, cov00, cov01, cov11, 
	double rotation_direction;
	//double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*m0_symmetry);
	double **inverted, **adjinverted, detcentre, **evecs1, *evals1, det=0, trace=0;
	double circumference, ref_angle;// longchord[2], longchordmag, storeeperp[2];
	double *evals=malloc(2*sizeof(double));// temppar[2];
	double sum_matrix_elements, matrix_element;
	double matrix_element_oldway, sum_matrix_elements_oldway=0.0;
	int N_line=0, main_ind, centre_ind, sec_ind;

	adjfieldline = malloc(L_fixedpoints*sizeof(struct position));
	evecs1 = malloc(2*sizeof(double*));
	evecs1[0] = malloc(2*sizeof(double)); evecs1[1] = malloc(2*sizeof(double));
	evals1 = malloc(2*sizeof(double));
	N_line = L_fixedpoints*N_gridphi_fieldperiod;
	angle_axis = malloc(L_fixedpoints*sizeof(double));
	//centre = malloc(L_fixedpoints*sizeof(struct position));
	fieldline_start.tangent = set_identity();
	adjfieldline->loc[0] = 0.0; adjfieldline->loc[1] = 0.0;
	varphi = 0.0;
	ref_angle = atan2(fieldline[0].loc[1] - axis[0].loc[1], fieldline[0].loc[0] - axis[0].loc[0]);
	printf("ref_angle=%f\n", ref_angle);
	angle_axis[0] = 0.0; 
	//printstructposition("fieldline[i]\n", &fieldline);
	number_turns=0;
	for (centre_ind = 0; centre_ind< L_fixedpoints; centre_ind++) 
		adjfieldline[centre_ind].tangent = adj2x2(fieldline[centre_ind].tangent, &detcentre);
	//for (i=1; i<N_line+1; i++) {
	//	RK4step_lambdacirc_mutangent(&fieldline, &adjfieldline, varphi, dvarphi, Bfield_saved[i-1]);
	//	varphi += dvarphi;
	//	if (i%N_gridphi_fieldperiod==0) {
	//		centre_ind = (i / N_gridphi_fieldperiod);
	//		//printf("DEBUG: centre_ind=%d/%d\n", centre_ind%L_fixedpoints, L_fixedpoints);
	//		centre[centre_ind % L_fixedpoints] = fieldline;
	//		adjcentre[centre_ind % L_fixedpoints] = adjfieldline;
	//		//printf("varphi = %f\n", varphi);
	//		//printstructposition("fieldline[i]\n", &fieldline);
	//		rawangle_principal = atan2(fieldline.loc[1] - axis[i%N_gridphi_fieldperiod].loc[1], fieldline.loc[0] - axis[i%N_gridphi_fieldperiod].loc[0]);
	//		//angle_axis[centre_ind%L_fixedpoints] = remainder(rawangle, 2.0*M_PI);
	//		rawangle = rawangle_principal + number_turns*M_PI*2.0;
	//		printf("rawangle_principal = %f\n", rawangle_principal);
	//		if (rawangle - angle_axis[(centre_ind-1)%L_fixedpoints] - M_PI > 0.0) {
	//			number_turns -= 1;
	//		}
	//		else if (rawangle - angle_axis[(centre_ind-1)%L_fixedpoints] + M_PI < 0.0) {
	//			number_turns += 1;
	//		}
	//		if (i/N_gridphi_fieldperiod == 1) {
	//			number_turns=0;
	//		}
	//		//printf("number_turns = %d\n", number_turns);
	//		rawangle = rawangle_principal + number_turns*M_PI*2.0;
	//		angle_axis[centre_ind%L_fixedpoints] = rawangle;
	//		//if (angle_axis[centre_ind%L_fixedpoints] < - M_PI) {
	//		//	angle_axis[centre_ind%L_fixedpoints] += M_PI;	
	//		//}
	//		//else if (angle_axis[centre_ind%L_fixedpoints] > M_PI) {
	//		//	angle_axis[centre_ind%L_fixedpoints] -= M_PI;	
	//		//}
	//		//angle_axis[centre_ind%L_fixedpoints] = atan2(fieldline.loc[1] - axis[1], fieldline.loc[0] - axis[0]) - ref_angle;
	//		//printstructposition("centre", centre + ( (centre_ind-1) % L_fixedpoints ) );
	//		inverted = invert2x2(centre[(centre_ind-1) % L_fixedpoints].tangent, &detcentre);
	//		adjinverted = invert2x2(adjcentre[(centre_ind-1) % L_fixedpoints].tangent, &adjdetcentre);
	//		//printmat("inverted", inverted, 2, 2);
	//		//ptarray[0] = &centre[centre_ind % L_fixedpoints].tangent[0][0];
	//		//ptarray[1] = &centre[centre_ind % L_fixedpoints].tangent[1][0];
	//		ext_centre[centre_ind % L_fixedpoints].loc[0] = fieldline.loc[0];
	//		ext_centre[centre_ind % L_fixedpoints].loc[1] = fieldline.loc[1];
	//		ext_centre[centre_ind % L_fixedpoints].part_tangent = multiply2x2(centre[centre_ind % L_fixedpoints].tangent, inverted, 2);
	//		ext_centre[centre_ind % L_fixedpoints].adj_part_tangent = multiply2x2(adjcentre[centre_ind % L_fixedpoints].tangent, adjinverted, 2);
	//		phidata[centre_ind % L_fixedpoints] = varphi;
	//		//printmat("ext_centre.part_tangent", ext_centre[centre_ind % L_fixedpoints].part_tangent, 2, 2);
	//		free(inverted[0]); free(inverted[1]);
	//		free(inverted);
	//		free(adjinverted[0]); free(adjinverted[1]);
	//		free(adjinverted);
	//	}
	//}
	for (centre_ind=0; centre_ind<L_fixedpoints; centre_ind++) {
		//printf("DEBUG: centre_ind=%d/%d\n", centre_ind%L_fixedpoints, L_fixedpoints);
		//centre[centre_ind] = ;
		//adjcentre[centre_ind % L_fixedpoints] = adjfieldline;
		//printf("varphi = %f\n", varphi);
		//printstructposition("fieldline[i]\n", &fieldline);
		rawangle_principal = atan2(fieldline[centre_ind].loc[1] - axis[centre_ind].loc[1], fieldline[centre_ind].loc[0] - axis[centre_ind].loc[0]);
		//angle_axis[centre_ind%L_fixedpoints] = remainder(rawangle, 2.0*M_PI);
		rawangle = rawangle_principal + number_turns*M_PI*2.0;
		//printf("rawangle_principal = %f\n", rawangle_principal);
		if (rawangle - angle_axis[(centre_ind+L_fixedpoints-1)%L_fixedpoints] - M_PI > 0.0) {
			number_turns -= 1;
		}
		else if (rawangle - angle_axis[(centre_ind+L_fixedpoints-1)%L_fixedpoints] + M_PI < 0.0) {
			number_turns += 1;
		}
		if (i/N_gridphi_fieldperiod == 1) {
			number_turns=0;
		}
		//printf("number_turns = %d\n", number_turns);
		rawangle = rawangle_principal + number_turns*M_PI*2.0;
		angle_axis[centre_ind] = rawangle;
		//printstructposition("fieldline", fieldline + (centre_ind+L_fixedpoints-1) % L_fixedpoints );
		if (centre_ind == 1) { 
		// when evaluating the partial tangent map for 0-->1 (here indexed 1, perhaps 0 would have been better)
		// one must take care to set the tangent map evaluated at 0 to be the identity 
		// instead of the full orbit tangent map which is stored there
		// this explains this seemingly arbitrary (but NECESSARY) exception
			inverted = set_identity();
			adjinverted = set_identity();
		}
		else {
			inverted = invert2x2(fieldline[(centre_ind+L_fixedpoints-1) % L_fixedpoints].tangent, &detcentre);
			adjinverted = invert2x2(adjfieldline[(centre_ind+L_fixedpoints-1) % L_fixedpoints].tangent, &detcentre);
		}
		//adjinverted = invert2x2(adjcentre[(centre_ind+L_fixedpoints-1) % L_fixedpoints].tangent, &adjdetcentre);
		//transpose2x2reassign(adjfieldline[centre_ind].tangent);
		//adjfieldline[centre_ind].tangent[0][1] = 
		//printmat("inverted", inverted, 2, 2);
		ext_centre[centre_ind].loc[0] = fieldline[centre_ind].loc[0];
		ext_centre[centre_ind].loc[1] = fieldline[centre_ind].loc[1];
		ext_centre[centre_ind].part_tangent = multiply2x2(fieldline[centre_ind].tangent, inverted, 2);
		//printmat("part", ext_centre[centre_ind].part_tangent, 2, 2);
		ext_centre[centre_ind].adj_part_tangent = multiply2x2(adjfieldline[centre_ind].tangent, adjinverted, 2);
		//printmat("adj_part", ext_centre[centre_ind].adj_part_tangent, 2, 2);
		//phidata[centre_ind] = varphi;
		//printmat("ext_centre.part_tangent", ext_centre[centre_ind % L_fixedpoints].part_tangent, 2, 2);
		free(inverted[0]); free(inverted[1]);
		free(inverted);
		free(adjinverted[0]); free(adjinverted[1]);
		free(adjinverted);
	}
	//gsl_fit_linear(phidata, 1, angle_axis, 1, L_fixedpoints, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
	//for (i=0; i<L_fixedpoints; i++) printf("islandcentre[%d] = (%f, %f)\n", i, ext_centre[i].loc[0], ext_centre[i].loc[1]);
	//c0 = remainder(13.9,3.0);
	//printf("c0=%f, c1=%f\n\n\n", c0, c1);
	//if (angle_axis[1] > 0.0) {
	//	if ( (angle_axis[2] > angle_axis[1]) || (angle_axis[2] < 0.0) ) clockwise = 0;
	//	else  clockwise = 1;
	//}
	//else {
	//	if ( (angle_axis[2] < angle_axis[1]) || (angle_axis[2] > 0.0) ) clockwise = 1;
	//	else clockwise = 0;
	//}
	//inverted = invert2x2(centre[(centre_ind-1) % L_fixedpoints].tangent, &detcentre);
	circumference = 0.0;
	number_turns = abs(number_turns);
	printf("number_turns = %d\n", number_turns);
	for (centre_ind=0;centre_ind<L_fixedpoints;centre_ind++) {
		//printf("DEBUG: centre_ind=%d/%d\n", centre_ind, L_fixedpoints);
		//if (centre_ind != 0) {
		//	printf("angle_axis[%d]=%f, angle_axis[%d]=%f\n", centre_ind, angle_axis[centre_ind], (centre_ind-1)%L_fixedpoints, angle_axis[(centre_ind-1)%L_fixedpoints]);
		//	//if ( (clockwise == 0) && ( ( (angle_axis[centre_ind] > small) && (angle_axis[centre_ind-1] < - small) ) ) ) (*n_turns) +=1;	
		//	//else if ( (clockwise == 1) && ( ( (angle_axis[centre_ind] < -small) && (angle_axis[centre_ind-1] > small) ) ) ) (*n_turns) +=1;	
		//	*n_turns = number_turns;
		//}
		//angle_axis[centre_ind] = atan2(ext_centre[(centre_ind+1)%L_fixedpoints].loc[1]-axis[1], ext_centre[(centre_ind+1)%L_fixedpoints].loc[0]-axis[0]);
		ext_centre[centre_ind].epar = malloc(2*sizeof(double)); ext_centre[centre_ind].eperp = malloc(2*sizeof(double));
		ext_centre[centre_ind].full_tangent = multiply2x2(ext_centre[(centre_ind +2)%L_fixedpoints].part_tangent, ext_centre[(centre_ind+1)%L_fixedpoints].part_tangent, 2);
		ext_centre[centre_ind].adj_full_tangent = multiply2x2(ext_centre[(centre_ind +2)%L_fixedpoints].adj_part_tangent, ext_centre[(centre_ind+1)%L_fixedpoints].adj_part_tangent, 2);
		for (sec_ind=3;sec_ind<L_fixedpoints+1;sec_ind++) {
			//printmat("on the way ext_centre.fulltangent", ext_centre[centre_ind].full_tangent, 2, 2);
			multiply2x2reassign(ext_centre[(centre_ind + sec_ind)%L_fixedpoints].part_tangent, ext_centre[centre_ind].full_tangent, 2);
			multiply2x2reassign(ext_centre[(centre_ind + sec_ind)%L_fixedpoints].adj_part_tangent, ext_centre[centre_ind].adj_full_tangent, 2);
		}
		//printmat("fulltangent", fieldline[L_fixedpoints-1].tangent, 2, 2);
		//printmat("ext_centre.fulltangent", ext_centre[centre_ind].full_tangent, 2, 2);
		//printmat("ext_centre._adj_fulltangent", ext_centre[centre_ind].adj_full_tangent, 2, 2);
		symmeigs(ext_centre[centre_ind].full_tangent, ext_centre[centre_ind].eperp, ext_centre[centre_ind].epar, evals);

		rhominus = ( L_fixedpoints + (number_turns*centre_ind)%L_fixedpoints - 1) %L_fixedpoints;
		rhoplus =  ( L_fixedpoints + (number_turns*centre_ind)%L_fixedpoints + 1) %L_fixedpoints;
		//kminus = ( ( (rhominus%number_turns)*L_fixedpoints + rhominus ) / number_turns ) %L_fixedpoints;
		//printf("DEBUG: centre_ind = %d/%d\n", centre_ind, L_fixedpoints);
		kminus = inverserho(rhominus, number_turns, L_fixedpoints);
		kplus = inverserho(rhoplus, number_turns, L_fixedpoints);
		ext_centre[centre_ind].chord[0] = ext_centre[centre_ind].loc[0] - ext_centre[kminus].loc[0];
		ext_centre[centre_ind].chord[1] = ext_centre[centre_ind].loc[1] - ext_centre[kminus].loc[1];
		ext_centre[centre_ind].chordplus[0] = ext_centre[kplus].loc[0] - ext_centre[centre_ind].loc[0];
		ext_centre[centre_ind].chordplus[1] = ext_centre[kplus].loc[1] - ext_centre[centre_ind].loc[1];
	}
	linalg2x2(ext_centre[0].full_tangent, evecs1, evals1, &det, &trace);
	for (main_ind=0; main_ind<L_fixedpoints; main_ind++) 
	ext_centre[main_ind].Res = 0.5 - 0.25*(ext_centre[main_ind].full_tangent[0][0] + ext_centre[main_ind].full_tangent[1][1]);
	if (1.0/ext_centre[0].Res > 1.0) { // condition for periodic field line to be island centre
		printf("Periodic field line is O point\n");
		if (*q0_index ==0) {
			// evals1 is the angle of rotation alpha in the paper
			omega = m0_symmetry*evals1[1]/(2.0*M_PI*L_fixedpoints); 
			// omega is different compared to Cary and Hanson's paper, including a factor of m0 
			q0 = (int) (m0_symmetry/(4.0*omega) - L_fixedpoints/2.0 + 0.5);
			*q0_index = q0;
	if ((fabs(evecs1[0][0]) > small) || (fabs(evecs1[0][1]) > small)) printf("WARNING: q0 is probably garbage\n");
			//printf("q0=%d\n", q0); //printf("omega=%f\n", omega); //printf("angle=%f\n", evals[1]);
		}
		else q0 = *q0_index;	
		if (number_turns == 0) {
			printf("ERROR: the code will fail because number_turns = 0\n");
			printf("number_turns= %d\n", (number_turns));
		}
		for (main_ind=0;main_ind<L_fixedpoints; main_ind++) {
			rotation_direction = ( ext_centre[main_ind].eperp[0]*ext_centre[main_ind].epar[1] - ext_centre[main_ind].eperp[1]*ext_centre[main_ind].epar[0] ) * (ext_centre[main_ind].full_tangent[1][0]/ fabs(ext_centre[main_ind].full_tangent[1][0]) );
			ext_centre[main_ind].epar[0] *= rotation_direction;
			ext_centre[main_ind].epar[1] *= rotation_direction;
			linalg2x2(ext_centre[0].full_tangent, evecs1, evals1, &det, &trace);
			omega = m0_symmetry*evals1[1]/(2.0*M_PI*L_fixedpoints); 
			// omega is different compared to Cary and Hanson's paper, including a factor of m0 
			q0ver = (int) (m0_symmetry/(4.0*omega) - L_fixedpoints/2.0 + 0.5);
	//if (q0 != q0ver) printf("WARNING: indices q0 are not the same when calculated from all island centres\n");
			//printf("evecs[0][0] = %f, evecs[0][1] = %f\n", evecs1[0][0], evecs1[0][1]);
			//printf("evecs[1][0] = %f, evecs[1][1] = %f\n", evecs1[0][0], evecs1[0][1]);
			//printf("evals[0] = %f, evals[1] = %f\n", evals1[0], evals1[1]);
			//printf("angle=%f\n", ext_centre[main_ind].angle);
			//printf("q0=%d\n", ext_centre[main_ind].q0_index);
			ext_centre[main_ind].long_tangent = malloc(L_fixedpoints*sizeof(double));
			ext_centre[main_ind].sperp = malloc(L_fixedpoints*sizeof(double));
			for (centre_ind=0;centre_ind<L_fixedpoints;centre_ind++) {
				//ext_centre[centre_ind].long_tangent = set_identity();
				//for (sec_ind=0; sec_ind<q0/L_fixedpoints; sec_ind++) {
				//	multiply2x2reassign(ext_centre[0].full_tangent, ext_centre[centre_ind].long_tangent, 2); 
				//}
				//for (sec_ind=0; sec_ind<q0%L_fixedpoints; sec_ind++) {
				//	multiply2x2reassign(ext_centre[sec_ind+1].part_tangent, ext_centre[centre_ind].long_tangent, 2); 
				//}
				ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints] = set_identity();
				ext_centre[main_ind].sperp[(main_ind+q0+centre_ind)%L_fixedpoints] = malloc(2*sizeof(double));
				for (sec_ind=0; sec_ind<(q0+centre_ind)/L_fixedpoints; sec_ind++) {
					multiply2x2reassign(ext_centre[main_ind].full_tangent, ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints], 2); 
				}
				for (sec_ind=0; sec_ind<(q0+centre_ind)%L_fixedpoints; sec_ind++) {
					multiply2x2reassign(ext_centre[(sec_ind+1+main_ind)%L_fixedpoints].part_tangent, ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints], 2); 
				}
				ext_centre[main_ind].sperp[(main_ind+q0+centre_ind)%L_fixedpoints][0] = ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints][0][0]*ext_centre[main_ind].eperp[0] + ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints][0][1]*ext_centre[main_ind].eperp[1];
				ext_centre[main_ind].sperp[(main_ind+q0+centre_ind)%L_fixedpoints][1] = ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints][1][0]*ext_centre[main_ind].eperp[0] + ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints][1][1]*ext_centre[main_ind].eperp[1];
			}
		}
		for (main_ind=0;main_ind<L_fixedpoints; main_ind++) {
			//printf("main_ind=%d/%d\n", main_ind, L_fixedpoints);
			sum_matrix_elements = sum_matrix_elements_oldway = circumference = 0.0;		
			for (centre_ind=0;centre_ind<L_fixedpoints; centre_ind++) {
				circumference += pow( pow(ext_centre[centre_ind].chord[0], 2.0) + pow(ext_centre[centre_ind].chord[1], 2.0) , 0.5 );
				matrix_element_oldway = inner(ext_centre[centre_ind].epar, 
						       ext_centre[main_ind].long_tangent[centre_ind], 
						       ext_centre[main_ind].eperp);	
				matrix_element = matrix_element_oldway;
				//matrix_element = ext_centre[centre_ind].epar[0]*ext_centre[main_ind].sperp[centre_ind][0] + ext_centre[centre_ind].epar[1]*ext_centre[main_ind].sperp[centre_ind][1];
				sum_matrix_elements_oldway += fabs(matrix_element_oldway);
				//ext_centre[main_ind].sign[centre_ind] = (int) ( matrix_element/fabs(matrix_element) );
				// TEMPORARY
				if ( ( matrix_element/fabs(matrix_element) ) < 0.0 ) { 
					if (main_ind != 0) printf("Warning: Not supposed to be here\n");
					ext_centre[centre_ind].epar[0] *= (-1);
					ext_centre[centre_ind].epar[1] *= (-1);
					ext_centre[centre_ind].eperp[0] *= (-1);
					ext_centre[centre_ind].eperp[1] *= (-1);
				}

				//printf("epar = (%f, %f)\n", ext_centre[centre_ind].epar[0], ext_centre[centre_ind].epar[1]);
				//printf("sperp = (%f, %f)\n", ext_centre[main_ind].sperp[centre_ind][0], ext_centre[main_ind].sperp[centre_ind][1]);
				//printf("sign[%d][%d] = %d\n", main_ind, centre_ind, sign[main_ind][centre_ind]);
				//sum_matrix_elements += fabs(matrix_element);
				//printf("matrix_element_%d,%d = %f\n", main_ind, centre_ind, matrix_element);
				sum_matrix_elements += fabs(matrix_element);
				//new_centre_ind = ( centre_ind - ext_centre[main_ind].q0_ind%L_fixedpoints + L_fixedpoints ) % L_fixedpoints;
				//sum_matrix_elements += fabs(matrix_element)/sin(ext_fieldline[main_ind].angle*(ext_fieldline[main_ind].q0_ind+new_centre_ind)/L_fixedpoints); //This has the sin piece in the denominator, which can be safely set to unity
			}
			//circumference = 2.0*M_PI*0.2094;
			//printf("displacement around = (%f, %f)\n", displacementaround[0], displacementaround[1]);
			//printf("circumference = %f\n", circumference);
			//printf("sum_matrix_elements = %f\n", sum_matrix_elements);
			//printf("sum_matrix_elements_oldway = %f\n", sum_matrix_elements_oldway);
			//wperp[main_ind] = 2.0*L_fixedpoints*circumference/(M_PI*pol_mode*sum_matrix_elements);
			ext_centre[main_ind].width = 2.0*L_fixedpoints*circumference/(M_PI*pol_mode*sum_matrix_elements);
			ext_centre[main_ind].Sigma = sum_matrix_elements;
			ext_centre[main_ind].circ = circumference;
			printf("Residue = %f for ind= %d\n", ext_centre[main_ind].Res, main_ind);
			printf("width = %f for ind= %d\n", ext_centre[main_ind].width, main_ind);
			printf("Sigma = %f for ind= %d\n", ext_centre[main_ind].Sigma, main_ind);
		}
	}
	else {
		printf("Periodic field line is X point\n");
		circumference = 0.0;		
		for (centre_ind=0;centre_ind<L_fixedpoints; centre_ind++) {
			circumference += pow( pow(ext_centre[centre_ind].chord[0], 2.0) + pow(ext_centre[centre_ind].chord[1], 2.0) , 0.5 );
		}
		for (main_ind=0; main_ind< L_fixedpoints; main_ind++)
		ext_centre[main_ind].circ = circumference;
		printf("circumference = %f\n", ext_centre[0].circ);
		printf("Residue = %f for ind= %d\n", ext_centre[0].Res, main_ind);
	}
	//for (centre_ind = 0; centre_ind<L_fixedpoints; centre_ind++) {
	//	longchord[0] = ext_centre[centre_ind].chord[0] + ext_centre[centre_ind].chordplus[0];
	//	longchord[1] = ext_centre[centre_ind].chord[1] + ext_centre[centre_ind].chordplus[1];
	//	longchordmag = pow(longchord[0]*longchord[0] + longchord[1]*longchord[1], 0.5);
	//	longchord[0] /= longchordmag;
	//	longchord[1] /= longchordmag;

	//	//printf("chord = (%f, %f)\n", longchord[0], longchord[1]);
	//	if ( fabs(longchord[0]*ext_centre[centre_ind].epar[0] + longchord[1]*ext_centre[centre_ind].epar[1]) < fabs(longchord[0]*ext_centre[centre_ind].eperp[0] + longchord[1]*ext_centre[centre_ind].eperp[1]) ) {
	//		//printf("\n\n\n\nYOLO\n\n\n\n\n");
	//		storeeperp[0] = ext_centre[centre_ind].eperp[0] ;
	//		storeeperp[1] = ext_centre[centre_ind].eperp[1] ;
	//		ext_centre[centre_ind].eperp[0] = ext_centre[centre_ind].epar[0];
	//		ext_centre[centre_ind].eperp[1] = ext_centre[centre_ind].epar[1];
	//		ext_centre[centre_ind].epar[0]  = storeeperp[0] ;
	//		ext_centre[centre_ind].epar[1]  = storeeperp[1] ;
	//	}
	//	if ( (longchord[0]*ext_centre[centre_ind].epar[0] + longchord[1]*ext_centre[centre_ind].epar[1]) < 0.0 ) {
	//		//printf("\n\n\n\nBROOLO\n\n\n\n\n");
	//		ext_centre[centre_ind].eperp[0] *= (-1);
	//		ext_centre[centre_ind].eperp[1] *= (-1);
	//		ext_centre[centre_ind].epar[0]  *= (-1);
	//		ext_centre[centre_ind].epar[1]  *= (-1);
	//	}
	//      printf("eperp = (%f, %f) and epar = (%f, %f)\n", ext_centre[centre_ind].eperp[0], ext_centre[centre_ind].eperp[1], ext_centre[centre_ind].epar[0], ext_centre[centre_ind].epar[1]);
	//      printf("position is (%f, %f)\n", ext_centre[centre_ind].loc[0], ext_centre[centre_ind].loc[1]);
	//	//if ( (longchord[0]*ext_centre[centre_ind].epar[1] + longchord[1]*ext_centre[centre_ind].epar[0]) > 0.0 ) {
	//	//	printf("\n\n\n\nMAGOOLO\n\n\n\n\n");
	//	//	ext_centre[centre_ind].eperp[0] *= (-1);
	//	//	ext_centre[centre_ind].eperp[1] *= (-1);
	//	//	//ext_centre[centre_ind].epar[0]  *= (-1);
	//	//	//ext_centre[centre_ind].epar[1]  *= (-1);
	//	//}
	//}
	//clock_t int3 = clock();
	//printf("q0 = %d\n", *q0_index);
} 

double *calc_islandwidth(double *circ, double *Sigma, struct ext_position *ext_fieldline, int m0_symmetry, int L_fixedpoints, int pol_mode) {
	int main_index, centre_index;// new_centre_index;
	double *wperp, sum_matrix_elements, matrix_element;
	double circumference;
	double displacementaround[2];
	double matrix_element_oldway, sum_matrix_elements_oldway=0.0;
	//printf("L_fixedpoints=%d\n", L_fixedpoints);
	wperp = malloc(L_fixedpoints*sizeof(double));
	displacementaround[0] = displacementaround[1] = 0.0;
	for (main_index=0;main_index<L_fixedpoints; main_index++) {
		//printf("main_index=%d/%d\n", main_index, L_fixedpoints);
		sum_matrix_elements = sum_matrix_elements_oldway = circumference = 0.0;		
		for (centre_index=0;centre_index<L_fixedpoints; centre_index++) {
			//printf("centre_index=%d/%d\n", centre_index, L_fixedpoints);
			//centre_index = ( ext_fieldline[main_index].q0_index + new_centre_index ) %L_fixedpoints;
			displacementaround[0] += ext_fieldline[centre_index].chord[0];
			displacementaround[1] += ext_fieldline[centre_index].chord[1];
			circumference += pow( pow(ext_fieldline[centre_index].chord[0], 2.0) + pow(ext_fieldline[centre_index].chord[1], 2.0) , 0.5 );
			matrix_element_oldway = inner(ext_fieldline[centre_index].epar, 
					       ext_fieldline[main_index].long_tangent[centre_index], 
					       ext_fieldline[main_index].eperp);	
			matrix_element = matrix_element_oldway;
			//matrix_element = ext_fieldline[centre_index].epar[0]*ext_fieldline[main_index].sperp[centre_index][0] + ext_fieldline[centre_index].epar[1]*ext_fieldline[main_index].sperp[centre_index][1];
			sum_matrix_elements_oldway += fabs(matrix_element_oldway);
			//ext_fieldline[main_index].sign[centre_index] = (int) ( matrix_element/fabs(matrix_element) );
			// TEMPORARY
			if ( ( matrix_element/fabs(matrix_element) ) < 0.0 ) { 
				if (main_index != 0) printf("Warning: Not supposed to be here\n");
				ext_fieldline[centre_index].epar[0] *= (-1);
				ext_fieldline[centre_index].epar[1] *= (-1);
				ext_fieldline[centre_index].eperp[0] *= (-1);
				ext_fieldline[centre_index].eperp[1] *= (-1);
			}

			//printf("epar = (%f, %f)\n", ext_fieldline[centre_index].epar[0], ext_fieldline[centre_index].epar[1]);
			//printf("sperp = (%f, %f)\n", ext_fieldline[main_index].sperp[centre_index][0], ext_fieldline[main_index].sperp[centre_index][1]);
			//printf("sign[%d][%d] = %d\n", main_index, centre_index, sign[main_index][centre_index]);
			//sum_matrix_elements += fabs(matrix_element);
			//printf("matrix_element_%d,%d = %f\n", main_index, centre_index, matrix_element);
			sum_matrix_elements += fabs(matrix_element);
			//new_centre_index = ( centre_index - ext_fieldline[main_index].q0_index%L_fixedpoints + L_fixedpoints ) % L_fixedpoints;
			//sum_matrix_elements += fabs(matrix_element)/sin(ext_fieldline[main_index].angle*(ext_fieldline[main_index].q0_index+new_centre_index)/L_fixedpoints); //This has the sin piece in the denominator, which can be safely set to unity
		}
		//circumference = 2.0*M_PI*0.2094;
		//printf("displacement around = (%f, %f)\n", displacementaround[0], displacementaround[1]);
		//printf("circumference = %f\n", circumference);
		//printf("sum_matrix_elements = %f\n", sum_matrix_elements);
		//printf("sum_matrix_elements_oldway = %f\n", sum_matrix_elements_oldway);
		//wperp[main_index] = 2.0*L_fixedpoints*circumference/(M_PI*pol_mode*sum_matrix_elements);
		wperp[main_index] = 2.0*L_fixedpoints*circumference/(M_PI*pol_mode*sum_matrix_elements);
		//printf("width = %f for index= %d\n", wperp[main_index], main_index);
		Sigma[main_index] = sum_matrix_elements;
		*circ = circumference;
	}
	return wperp;
}

//double ***shapecircumference(double ***coils, int n_coils, int *n_segs, struct ext_position *ext_centre, struct position *lambda_circ, int m0_symmetry, int N_gridphi_fieldperiod, int tor_mode, int pol_mode) {
//	int i=0, q0_ind, main_ind, sec_ind, coil_ind, seg_ind, N_line=0, L_fixedpoints;
//	struct position centre, lambda, mu;
//	double varphi=0.0, chordlength, chordpluslength;
//	double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*m0_symmetry);
//	double *gradcircumference, *gradtangent, initialgradcirc=0.0, *gradwidth;
//	double ***shapecirc, sum_matrix_elements, matrix_element;
//	//double circumference2=ext_centre[0].circumference, circumference;
//	shapecirc = malloc(n_coils*sizeof(double));
//	for (coil_ind=0;coil_ind<n_coils;coil_ind++) {
//		shapecirc[coil_ind] = malloc(n_segs[coil_ind]*sizeof(double));
//		for (seg_ind=0; seg_ind<n_segs[coil_ind]; seg_ind++) {
//			shapecirc[coil_ind][seg_ind] = calloc(3,sizeof(double));
//		}
//	}
//	printf("tor_mode=%d\n", tor_mode);
//	if (tor_mode % m0_symmetry == 0) 
//		L_fixedpoints = pol_mode;	
//	else 				 
//		L_fixedpoints = m0_symmetry*pol_mode;
//	printf("L_fixedpoints=%d\n", L_fixedpoints);
//	centre.loc[0] = ext_centre[0].loc[0]; centre.loc[1] = ext_centre[0].loc[1];
//	//lambda_circ.loc[0] = lambda[0]; lambda_circ.loc[1] = lambda[1];
//	centre.tangent = set_identity();
//	//lambda_circ.tangent = set_identity();
//	gradcircumference = &initialgradcirc;
//	N_line = L_fixedpoints*N_gridphi_fieldperiod;
//	varphi = 0.0;
//	for (i=0; i<N_line; i++)
//	{
//		if (i%N_gridphi_fieldperiod==0)
//		{
//			main_ind = i/N_gridphi_fieldperiod;
//			chordlength= pow( pow(ext_centre[main_ind].chord[0], 2.0) + pow(ext_centre[main_ind].chord[1], 2.0), 0.5);
//			chordpluslength= pow( pow(ext_centre[main_ind].chordplus[0], 2.0) + pow(ext_centre[main_ind].chordplus[1], 2.0), 0.5);
//			lambda_circ->loc[0] += (ext_centre[main_ind].chord[0]/chordlength - ext_centre[main_ind].chordplus[0]/chordpluslength);
//			lambda_circ->loc[1] += (ext_centre[main_ind].chord[1]/chordlength - ext_centre[main_ind].chordplus[1]/chordpluslength);
//		}
//		RK4_adjshapecirc(shapecirc, &centre, lambda_circ, varphi, dvarphi, coils, n_coils, n_segs);
//		varphi += dvarphi;
//		//printf("gradcircumference = %f\n", *gradcircumference);
//	}
//	//printf("varphi = %f\n", varphi);
//	//printstructposition("lambda",lambda_circ);
//	//printstructposition("Xp",&centre);
//	//printf("gradcircumference = %f\n", *gradcircumference);
//	return shapecirc;
//}

struct position solve_lambda_circ(struct ext_position *ext_centre, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod) {
	int centre_ind, count = 0; 
	struct position lambda, lambda_start, deltalambda; //, fieldline_start, fieldline;
	double varphi=0.0, *jumptocentre;
	//double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*m0_symmetry);
	double *vecresult; // **result;
	double det_tangent;
	double **inverseTminusI; // **inverseT;
	double error = 1.0;
	double factor = FACTOR, chordlength, chordpluslength;
	int N_line=0;
	N_line = L_fixedpoints*N_gridphi_fieldperiod;
	printf("N_line = %d\n", N_line);
	printf("field_periods = %d\n", m0_symmetry);
	lambda.loc[0] = 1.0; lambda.loc[1] = 1.0;
	lambda.tangent = set_identity();
	lambda_start.tangent = set_identity();
	do {
		lambda_start.loc[0] = lambda.loc[0]; lambda_start.loc[1] = lambda.loc[1];
		lambda_start.tangent[0][0] = lambda.tangent[0][0]; lambda_start.tangent[0][1] = lambda.tangent[0][1];
		lambda_start.tangent[1][0] = lambda.tangent[1][0]; lambda_start.tangent[1][1] = lambda.tangent[1][1];
		varphi = 0.0;
		for (centre_ind=0; centre_ind<L_fixedpoints; centre_ind++) {
				chordlength= pow( pow(ext_centre[centre_ind].chord[0], 2.0) + pow(ext_centre[centre_ind].chord[1], 2.0), 0.5);
				chordpluslength= pow( pow(ext_centre[centre_ind].chordplus[0], 2.0) + pow(ext_centre[centre_ind].chordplus[1], 2.0), 0.5);
				lambda.loc[0] += (ext_centre[centre_ind].chord[0]/chordlength - ext_centre[centre_ind].chordplus[0]/chordpluslength);
				lambda.loc[1] += (ext_centre[centre_ind].chord[1]/chordlength - ext_centre[centre_ind].chordplus[1]/chordpluslength);
				vecresult = multiply(ext_centre[(centre_ind+1)%L_fixedpoints].adj_part_tangent, lambda.loc);
				multiply2x2reassign(ext_centre[(centre_ind+1)%L_fixedpoints].adj_part_tangent, lambda.tangent, 2);
				lambda.loc[0] = vecresult[0]; lambda.loc[1] = vecresult[1];
		}
		deltalambda = addstructs(1.0, &lambda, -1.0, &lambda_start); 
		inverseTminusI = invert2x2(deltalambda.tangent, &det_tangent);
		jumptocentre = multiply(inverseTminusI, &deltalambda.loc[0]);
		jumptocentre[0] *= (factor);
		jumptocentre[1] *= (factor);
		lambda.loc[0] = lambda_start.loc[0] - jumptocentre[0]; 
		lambda.loc[1] = lambda_start.loc[1] - jumptocentre[1];
		lambda.tangent[0][0] = 1.0; lambda.tangent[0][1]=0.0; lambda.tangent[1][0]=0.0; lambda.tangent[1][1]=1.0;
		error = sqrt( ( deltalambda.loc[0]*deltalambda.loc[0] + deltalambda.loc[1]*deltalambda.loc[1] ) / ( lambda.loc[0]*lambda.loc[0] + lambda.loc[1]*lambda.loc[1] ) ) ; 
		//printstructposition("lambda", &lambda);
		//printf("error = %15.15f\n", error);
		free(jumptocentre);
		free(vecresult);
		free(inverseTminusI[0]); free(inverseTminusI[1]);
		free(inverseTminusI);
		//free(result[0]); free(result[1]); free(result);
	//} while (error > smallerror );
		count+=1;
	} while(error>smallerror && count < MAX_COUNT);
	if (count == MAX_COUNT) printf("Warning: periodic lambda not found\n");
	return lambda;
}

struct position **solve_mu_tangent(struct ext_position *ext_centre, int m0_symmetry, int L_fixedpoints, int q0_ind, int N_gridphi_fieldperiod) {
	int centre_ind, main_ind, Q0_ind, Q_ind, q_ind;
	struct position mu_start, fieldline_start, fieldline, **mu;
	double *vecresult, **inverted, det;// varphi=0.0;
	int N_line;
	N_line = L_fixedpoints*N_gridphi_fieldperiod;
	mu = malloc(L_fixedpoints*sizeof(struct position*));
	for (main_ind=0; main_ind<L_fixedpoints; main_ind++)  {
		Q0_ind = q0_ind/L_fixedpoints + 2;
		//printf("q0_ind = %d\n", q0_ind);
		//printf("L_fixedpoints  = %d\n", L_fixedpoints);
		//printf("Q0_ind = %d\n", Q0_ind);
		mu[main_ind] = malloc((Q0_ind)*sizeof(struct position));
		mu_start.loc[0] = 0.0;  mu_start.loc[1] = 0.0; 
		mu_start.tangent = set_identity();
		//printf("N_line = %d\n", N_line);
		//printf("field_periods = %d\n", m0_symmetry);
		fieldline_start.loc[0] = ext_centre[main_ind].loc[0]; fieldline_start.loc[1] = ext_centre[main_ind].loc[1];
		fieldline.tangent = set_identity();
		fieldline.loc[0] = fieldline_start.loc[0]; fieldline.loc[1] = fieldline_start.loc[1];
		//printf("RR=%f\n", fieldline.loc[0]);
		for (Q_ind=Q0_ind-1; Q_ind>=0; Q_ind--) {
			mu[main_ind][Q_ind].tangent = set_identity(); // = mu_start; 
			mu[main_ind][Q_ind].loc[0] = mu_start.loc[0]; 
			mu[main_ind][Q_ind].loc[1] = mu_start.loc[1]; 
			for (centre_ind=L_fixedpoints-1; centre_ind>=0; centre_ind--) {
				q_ind = Q_ind*L_fixedpoints + centre_ind;
				inverted = invert2x2(ext_centre[(main_ind+centre_ind+1)%L_fixedpoints].adj_part_tangent, &det);
				vecresult = multiply(inverted, mu[main_ind][Q_ind].loc);
				////result = multiply2x2(inverted, mu[main_ind][Q_ind].tangent, 2);
				free(inverted[0]); free(inverted[1]);
				free(inverted);
				mu[main_ind][Q_ind].loc[0] = vecresult[0]; mu[main_ind][Q_ind].loc[1] = vecresult[1];
				if ( ( q_ind >= q0_ind ) && ( q_ind < q0_ind + L_fixedpoints ) ) {
					//mu[main_ind][Q_ind].loc[0] -= (ext_centre[main_ind].sign[(main_ind+centre_ind)%L_fixedpoints]*ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[0]);
					//mu[main_ind][Q_ind].loc[1] -= (ext_centre[main_ind].sign[(main_ind+centre_ind)%L_fixedpoints]*ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[1]);
					mu[main_ind][Q_ind].loc[0] -= ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[0];
					mu[main_ind][Q_ind].loc[1] -= ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[1];
					//printf("epar[0] = %f\n", ext_centre[centre_ind].epar[0]);
					//printf("epar[1] = %f\n", ext_centre[centre_ind].epar[1]);
				}
				//mu[main_ind][Q_ind].tangent[0][0] = result[0][0]; mu[main_ind][Q_ind].tangent[0][1] = result[0][1];
				//mu[main_ind][Q_ind].tangent[1][0] = result[1][0]; mu[main_ind][Q_ind].tangent[1][1] = result[1][1];
				//printf("main_ind=%d, Q_ind=%d, centre_ind=%d\n", main_ind, Q_ind, centre_ind);
				//printstructposition("mu", mu[main_ind]+Q_ind);
			}
			mu_start = mu[main_ind][Q_ind];
			//printf("Q_ind=%d\n", Q_ind);
			//printstructposition("mu\n", mu[main_ind]+Q_ind);
			//printstructposition("mu_start\n", &mu_start);

			//for (i=N_line-1; i>=0; i--)
			//{
			//	RK4_adjsimple(&fieldline, mu[main_ind]+Q_ind, varphi, -dvarphi, coils, n_coils, n_segs);
			//	if (i%N_gridphi_fieldperiod==0)
			//	{
			//		centre_ind = i/N_gridphi_fieldperiod;
			//		//printstructposition("mu\n", mu);
			//		//if (centre_ind==0) mu->loc[0] += 1.0;
			//		q_ind = Q_ind*L_fixedpoints + centre_ind;
			//		if ( ( q_ind >= q0_ind ) && ( q_ind < q0_ind + L_fixedpoints ) ) {
			//			//mu[main_ind][Q_ind].loc[0] = 0.0;
			//			//mu[main_ind][Q_ind].loc[1] = 0.0;
			//			mu[main_ind][Q_ind].loc[0] -= (ext_centre[(main_ind+centre_ind+1)%L_fixedpoints].epar[0]);
			//			mu[main_ind][Q_ind].loc[1] -= (ext_centre[(main_ind+centre_ind+1)%L_fixedpoints].epar[1]);
			//			//printf("epar[0] = %f\n", ext_centre[centre_ind].epar[0]);
			//			//printf("epar[1] = %f\n", ext_centre[centre_ind].epar[1]);
			//		}
			//		if ( centre_ind == 1) {
			//			printf("main_ind=%d, Q_ind=%d, centre_ind=%d\n", main_ind, Q_ind, centre_ind);
			//			printstructposition("mu", mu[main_ind]+Q_ind);
			//			//printstructposition("mu +", mu[main_ind]+Q_ind+1);
			//		}
			//	}
			//	varphi -= dvarphi;
			//}
			//mu_start = mu[main_ind][Q_ind];
		}
		//clock_t int3 = clock();
	}
	printf("exit solve_mu_tangent module \n");
	return mu;
}

struct position **solve_lambda_tangent(struct field **Bfield_saved, struct ext_position *ext_centre, struct position **mu, int m0_symmetry, int L_fixedpoints, int q0_ind, int N_gridphi_fieldperiod) {
	// declarations
	//clock_t start = clock();
	int N_line, i=0, centre_ind, main_ind, Q0_ind, Q_ind, q_ind, ind, count = 0;
	struct position fieldline_start, fieldline, sperp, muvar, lambdavar, deltalambda, **lambda, lambdain; //, fieldline_start, fieldline;
	double varphi=0.0, *jumptocentre; 
	double det_tangent;
	double **inverseTminusI, **inverseT;
	double error = 1.0;
	double **matrix2, *vector;
	double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*m0_symmetry), factor = 1.0;
	N_line = L_fixedpoints*N_gridphi_fieldperiod;
	lambda = malloc(L_fixedpoints*sizeof(struct position)); 
	matrix2 = set_identity();
	sperp.tangent = set_identity();
	lambdavar.tangent = set_identity();
	Q0_ind = q0_ind/L_fixedpoints + 2;
	fieldline.tangent = set_identity();
	//printf("Finding adjoint variable lambda_Q,k:\n");
	for (main_ind=0;main_ind<L_fixedpoints;main_ind++) {
		//printf("main_ind = %d/%d\n", main_ind, L_fixedpoints);
		lambda[main_ind] = malloc(Q0_ind*sizeof(struct position));
		fieldline_start.loc[0] = ext_centre[main_ind].loc[0]; fieldline_start.loc[1] = ext_centre[main_ind].loc[1];
		fieldline.tangent[0][0] = fieldline.tangent[1][1] = 1.0; 
		fieldline.tangent[1][0] = fieldline.tangent[0][1] = 0.0; 
		lambdavar.loc[0] = 0.0; lambdavar.loc[1] = 0.0;
		for (Q_ind=0;Q_ind<Q0_ind;Q_ind++) {
			factor = 1.0;
			lambda[main_ind][Q_ind].loc[0] = 0.0; lambda[main_ind][Q_ind].loc[1] = 0.0;
			lambda[main_ind][Q_ind].tangent = set_identity();
			lambdain.loc[0] = 0.0; lambdain.loc[1] = 0.0;
			lambdain.tangent = set_identity();
			count = 0;
			do {
				lambdavar.tangent[0][0] = 1.0; lambdavar.tangent[0][1] = 0.0; 
				lambdavar.tangent[1][0] = 0.0; lambdavar.tangent[1][1] = 1.0; 
				lambdavar.loc[0] = lambda[main_ind][Q_ind].loc[0]; 
				lambdavar.loc[1] = lambda[main_ind][Q_ind].loc[1]; 
				matrix2[0][0] = matrix2[1][1] = 1.0;
				matrix2[1][0] = matrix2[0][1] = 0.0;
				for (ind=0;ind<Q_ind; ind++) {
					multiply2x2reassign(ext_centre[main_ind].full_tangent, matrix2, 2);
				}
				vector = multiply(matrix2, ext_centre[main_ind].eperp);
				sperp.loc[0] = vector[0];  sperp.loc[1] = vector[1];
				sperp.tangent[0][0] = sperp.tangent[1][1] = 1.0; 
				sperp.tangent[1][0] = sperp.tangent[0][1] = 0.0; 
				fieldline.loc[0] = fieldline_start.loc[0];
				fieldline.loc[1] = fieldline_start.loc[1];
				fieldline.tangent[0][0] = fieldline.tangent[1][1] = 1.0; 
				fieldline.tangent[1][0] = fieldline.tangent[0][1] = 0.0; 
				muvar = mu[main_ind][Q_ind];
				varphi = main_ind*2.0*M_PI/m0_symmetry;
				for (i=0; i<N_line; i++) {
					if (i%N_gridphi_fieldperiod==0) {
						centre_ind = i/N_gridphi_fieldperiod;
						q_ind = Q_ind*L_fixedpoints + centre_ind;
						if ( q_ind >= q0_ind && q_ind < q0_ind + L_fixedpoints) {
							muvar.loc[0] += (ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[0]);
							muvar.loc[1] += (ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[1]);
						}
					}
					//printf("Q_ind = %d/%d, i=%d/%d\n", Q_ind, Q0_ind, i, N_line);
					RK4step_lambdatangent(&fieldline, &lambdain, &sperp, &muvar, varphi, dvarphi, Bfield_saved[(main_ind*N_gridphi_fieldperiod + i) % N_line]);
					varphi += dvarphi;
				}
				deltalambda = addstructs(1.0, &lambdavar, -1.0, &lambdain); 
				inverseTminusI = invert2x2(deltalambda.tangent, &det_tangent);
				inverseT = invert2x2(lambdain.tangent, &det_tangent);
				jumptocentre = multiply(inverseTminusI, &deltalambda.loc[0]);
				jumptocentre[0] *= (factor);
				jumptocentre[1] *= (factor);
				lambdain.loc[0] = lambdavar.loc[0] - jumptocentre[0]; 
				lambdain.loc[1] = lambdavar.loc[1] - jumptocentre[1];
				lambdain.tangent[0][0] = 1.0; lambdain.tangent[0][1]=0.0; lambdain.tangent[1][0]=0.0; lambdain.tangent[1][1]=1.0;
				error = sqrt((pow(deltalambda.loc[0], 2.0) + pow(deltalambda.loc[1], 2.0))/(lambdain.loc[0]*lambdain.loc[0] + lambdain.loc[1]*lambdain.loc[1])); 
				//printstructposition("lambdatangent", &lambdain);
				//printf("error = %15.15f\n", error);
				count++;
				lambda[main_ind][Q_ind].loc[0] = lambdain.loc[0];
				lambda[main_ind][Q_ind].loc[1] = lambdain.loc[1];
				//factor*=0.99;
			//} while (error > smallerror );
			} while(error>smallerror && count < MAX_COUNT);
			if (count == MAX_COUNT) printf("Warning: periodic lambda not found\n");
		}
	}
	return lambda;
}

struct position **solve_lambdaRes(struct field **Bfield_saved, struct position **Xp_saved, double **adjfulltangent, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod) {
	// declarations
	//clock_t start = clock();
	int N_line=L_fixedpoints*N_gridphi_fieldperiod, i=0, count = 0;
	struct position lambdavar, deltalambda, **lambda=malloc(N_line*sizeof(struct position)); 
	double varphi=0.0, *jumptocentre; 
	double det;
	double **lambdatan0, **inverseTminusI, **TminusI=set_identity();
	double error = 1.0;
	double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*m0_symmetry), factor = FACTOR;
	for (i=0; i<N_line; i++)  lambda[i] = malloc(4*sizeof(struct position)); 
	lambdavar.loc[0] = lambdavar.loc[1] = 0.0; 
	//lambdavar.tangent = set_identity();
	lambdatan0 = invert2x2(adjfulltangent, &det);
	lambdavar.tangent = invert2x2(adjfulltangent, &det);
	count = 0;
	//printf("Finding adjoint variable lambdaRes:\n");
	do {
		lambdavar.tangent[0][0] = lambdatan0[0][0];
		lambdavar.tangent[0][1] = lambdatan0[0][1];
		lambdavar.tangent[1][0] = lambdatan0[1][0];
		lambdavar.tangent[1][1] = lambdatan0[1][1];
		//printstructposition("lambdavar (start)", &lambdavar);
		varphi = 0.0;
		for (i=0; i<N_line; i++) {
			//printf("i=%d/%d\n", i, N_line); printstructposition("lambdavar (during)", &lambdavar);
			lambda[i][0] = lambdavar;
			RK4step_lambdaRes(&lambdavar, lambda[i], Xp_saved[i], varphi, dvarphi, Bfield_saved[i]);
			varphi += dvarphi;
		}
		//printstructposition("lambdavar (end)", &lambdavar); //printstructposition("lambda", lambda[0]);
		deltalambda = addstructs(1.0, &lambdavar, -1.0, lambda[0]); 
		TminusI[0][0] = lambdatan0[0][0] - 1.0;
		TminusI[0][1] = lambdatan0[0][1] ;
		TminusI[1][0] = lambdatan0[1][0] ;
		TminusI[1][1] = lambdatan0[1][1] - 1.0;
		inverseTminusI = invert2x2(TminusI, &det);
		jumptocentre = multiply(inverseTminusI, &deltalambda.loc[0]);
		jumptocentre[0] *= (factor);
		jumptocentre[1] *= (factor);
		lambdavar.loc[0] += jumptocentre[0]; 
		lambdavar.loc[1] += jumptocentre[1];
		error = sqrt((pow(deltalambda.loc[0], 2.0) + pow(deltalambda.loc[1], 2.0))/(lambdavar.loc[0]*lambdavar.loc[0] + lambdavar.loc[1]*lambdavar.loc[1])); 
		//printf("error = %15.15f\n", error);
		count++;
		//factor*=0.99;
	//} while (error > smallerror );
	} while(error>smallerror && count < MAX_COUNT);
	if (count == MAX_COUNT) printf("Warning: periodic lambda not found\n");
	//printstructposition("lambdaRes", &lambda);
	return lambda;
}

void solve_gradcirc(double *gradcircumference, struct field **Bfield_island, struct ext_position *ext_centre, struct position *lambda_circ, int L_fixedpoints, int N_gridphi_fieldperiod, struct fieldparams allparams, int diffparams_ind1, int diffparams_ind2) {
	// declarations
	//clock_t start = clock();
	int i=0, main_ind, N_line=0;
	int num_params, m0_fieldperiods = allparams.m0_fieldperiods;
	struct position centre;
	double varphi=0.0, chordlength, chordpluslength, **matrix2;
	double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*m0_fieldperiods);
	//double *gradtangent;

	num_params = allparams.n_diff;
	printf("num_params = %d\n", num_params);
	for (i=0;i<num_params; i++) {
		gradcircumference[i] = 0.0;
	}
	matrix2=set_identity();
	//gradtangent = malloc(L_fixedpoints*sizeof(double));
	centre.loc[0] = ext_centre[0].loc[0]; centre.loc[1] = ext_centre[0].loc[1];
	//lambda_circ.loc[0] = lambda[0]; lambda_circ.loc[1] = lambda[1];
	centre.tangent = set_identity();
	//lambda_circ.tangent = set_identity();
	//for (row=0; row<num_params; row++) {
	//	gradcircumference[row] = 0.0;
	//}
	N_line = L_fixedpoints*N_gridphi_fieldperiod;
	//printf("N_line = %d\n", N_line);
	//printf("field_periods = %d\n", m0_fieldperids);
	varphi = 0.0;
	for (i=0; i<N_line; i++)
	{
		//printf("i=%d/%d\n", i, N_line);
		if (i%N_gridphi_fieldperiod==0)
		{
			main_ind = i/N_gridphi_fieldperiod;
			chordlength= pow( pow(ext_centre[main_ind].chord[0], 2.0) + pow(ext_centre[main_ind].chord[1], 2.0), 0.5);
			chordpluslength= pow( pow(ext_centre[main_ind].chordplus[0], 2.0) + pow(ext_centre[main_ind].chordplus[1], 2.0), 0.5);
			lambda_circ->loc[0] += (ext_centre[main_ind].chord[0]/chordlength - ext_centre[main_ind].chordplus[0]/chordpluslength);
			lambda_circ->loc[1] += (ext_centre[main_ind].chord[1]/chordlength - ext_centre[main_ind].chordplus[1]/chordpluslength);
		}
		RK4step_gradcirc(gradcircumference, &centre, lambda_circ, varphi, dvarphi, Bfield_island[i], allparams, diffparams_ind1, diffparams_ind2);
		varphi += dvarphi;
		//printf("gradcircumference = %f\n", gradcircumference[0]);
	}
	return;
}

void solve_gradtangent(double **number, struct field **Bfield_island, struct ext_position *ext_centre, struct position **lambdaQ, struct position **muQ, int L_fixedpoints, int q0_ind, int N_gridphi_fieldperiod, struct fieldparams allparams, int diffparams_ind1, int diffparams_ind2) {
	int i=0, ind, centre_ind, Q0_ind, Q_ind, q_ind, main_ind;
	int num_params, m0_fieldperiods = allparams.m0_fieldperiods;
	struct position fieldline, lambda, mu, sperp;
	double varphi=0.0, *vecresult, **matrix2;
	double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*m0_fieldperiods);
	double checksign, width;
	//double **number;
	int N_line;
	num_params = allparams.n_diff; 

	matrix2 = set_identity();
	N_line = L_fixedpoints*N_gridphi_fieldperiod;
	fieldline.loc[0] = ext_centre[0].loc[0]; fieldline.loc[1] = ext_centre[0].loc[1];
	fieldline.tangent = set_identity();
	//centre.loc[0] = ext_centre[0].loc[0]; centre.loc[1] = ext_centre[0].loc[1];
	////lambda_centre.loc[0] = lambda[0]; lambda_centre.loc[1] = lambda[1];
	//centre.tangent = set_identity();
	//lambda_centre.tangent = set_identity();
	Q0_ind = q0_ind/L_fixedpoints + 2;
	//number = malloc(L_fixedpoints*sizeof(double));
	sperp.tangent = set_identity();
	//printf("N_line = %d\n", N_line);
	//printf("field_periods = %d\n", m0_fieldperiods);
	for (main_ind=0; main_ind<L_fixedpoints; main_ind++) {
		width = 0.0;
		varphi = main_ind*2.0*M_PI/m0_fieldperiods;
		//number[main_ind] = calloc(num_params, sizeof(double)); 
		for (i=0;i<num_params;i++) {
			number[main_ind][i] = 0.0;
		}
		fieldline.tangent[0][0] = fieldline.tangent[1][1] = 1.0;
		fieldline.tangent[0][1] = fieldline.tangent[1][0] = 0.0;
		//sperp.loc[0] = ext_centre[main_ind].eperp[0]; sperp.loc[1] = ext_centre[main_ind].eperp[1];
		for (Q_ind=0;Q_ind<Q0_ind;Q_ind++) {
			lambda = lambdaQ[main_ind][Q_ind];
			mu = muQ[main_ind][Q_ind];
			//printstructposition("lambda", &lambda);
			//printstructposition("mu", &mu);
			fieldline.loc[0] = ext_centre[main_ind].loc[0]; fieldline.loc[1] = ext_centre[main_ind].loc[1];
			matrix2[0][0] = matrix2[1][1] = 1.0;
			matrix2[0][1] = matrix2[1][0] = 0.0;
			for (ind=0;ind<Q_ind; ind++) {
				multiply2x2reassign(ext_centre[main_ind].full_tangent, matrix2, 2);
			}
			//printstructposition("fieldline", &fieldline);
			//printf("Q_ind=%d\n", Q_ind);
			//printmat("matrix2", matrix2, 2,2 );
			vecresult = multiply(matrix2, ext_centre[main_ind].eperp);
			sperp.tangent[0][0] = sperp.tangent[1][1] = 1.0; 
			sperp.tangent[1][0] = sperp.tangent[0][1] = 0.0; 
			//vecresult = multiply(fieldline.tangent, ext_centre[main_ind].eperp);
			sperp.loc[0] = vecresult[0]; sperp.loc[1] = vecresult[1];
			free(vecresult);
			//printf("sperp=%f\n", sperp.tangent[1][1]);
			//fieldline->tangent[0][0] = 1.0 ; fieldline->tangent[0][1] = 0.0;
			//fieldline->tangent[1][0] = 0.0 ; fieldline->tangent[0][1] = 1.0;
			for (i=0; i<N_line; i++)
			{
				if (i%N_gridphi_fieldperiod==0)
				{
					centre_ind = i/N_gridphi_fieldperiod;
					//printf("varphi = %f\n", varphi);
					//printstructposition("Xp", &fieldline);
					//printf("number = %f\n", number[centre_ind]);
					//if (i/N_gridphi_fieldperiod==0) lambda_centre->loc[0] += 1.0;
					//printstructposition("lambda\n", lambda);
					//if (centre_ind==1) 
					//{
					//	printf("main_ind=%d, Q_ind=%d, centre_ind=%d\n", main_ind, Q_ind, centre_ind); 
					//	printf("number = %f\n", number[main_ind]);
					//	printstructposition("mu", &mu);
					//	printstructposition("lambda", &lambda);
					//}
					q_ind = Q_ind*L_fixedpoints + centre_ind;
					if ( q_ind >= q0_ind && q_ind < q0_ind + L_fixedpoints) {
						//mu.loc[0] += (sign[main_ind][centre_ind]*ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[0]);
						//mu.loc[1] += (sign[main_ind][centre_ind]*ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[1]);
						//checksign = sperp.loc[0]*ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[0] + sperp.loc[1]*ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[1];
						checksign = sperp.loc[0]*ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[0] + sperp.loc[1]*ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[1];
						width += fabs(checksign);
						checksign /= fabs(checksign);
						
						//printf("epar = (%f, %f)\n", ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[0], ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[1]);
//						printf("sperp = (%f, %f)\n", sperp.loc[0], sperp.loc[1]);
						//printf("sign_%d,%d = %d\n", main_ind, centre_ind, sign[main_ind][(main_ind+centre_ind)%L_fixedpoints]);
						//printf("sign_%d,%d = %d\n", main_ind, centre_ind, ext_centre[main_ind].sign[(main_ind+centre_ind)%L_fixedpoints]);
						//printf("checksign_%d,%d = %f\n", main_ind, centre_ind, checksign);
//
						//printf("sperp (stored) = (%f, %f)\n", ext_centre[main_ind].sperp[(centre_ind+main_ind)%L_fixedpoints][0], ext_centre[main_ind].sperp[(centre_ind+main_ind)%L_fixedpoints][1]);
						//printf("sperp (calculated) = (%f, %f)\n", sperp.loc[0], sperp.loc[1]);
//
						//mu.loc[0] += (ext_centre[main_ind].sign[(main_ind+centre_ind)%L_fixedpoints]*ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[0]);
						//mu.loc[1] += (ext_centre[main_ind].sign[(main_ind+centre_ind)%L_fixedpoints]*ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[1]);
						mu.loc[0] += (ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[0]);
						mu.loc[1] += (ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[1]);
					}
				}
				//RK4step_gradtangent(double *number, struct position *Xp, struct position *lambda, struct position *sperp, struct position *mu, double varphi, double dvarphi, struct field *Bfield_saved, struct fieldparams allparams, int diffparam_ind1, int diffparam_ind2) 
				RK4step_gradtangent(number[main_ind], &fieldline, &lambda, &sperp, &mu, varphi, dvarphi, Bfield_island[(main_ind*N_gridphi_fieldperiod + i) % N_line], allparams, diffparams_ind1, diffparams_ind2);
				//RK4step_lambdatangent(&fieldline, &lambda, &sperp, &mu, varphi, dvarphi, Bfield_island[i]);
				//RK4step_lambdacirc_mutangent(&fieldline, &mu, varphi, dvarphi, Bfield_island[i]);
				varphi += dvarphi;
			}
		}
		//printf("varphi = %f\n", varphi);
		//printstructposition("lambda", &lambda);
		//printstructposition("Xp",&fieldline);
	}
	//for (centre_ind=0; centre_ind<L_fixedpoints; centre_ind++) {
	//	//printf("number = %f\n", number[centre_ind]);
	//}
	return;
}

void solve_gradRes(double *number, struct field **Bfield_island, struct position **fieldline, struct position **lambdamu, int L_fixedpoints, int N_gridphi_fieldperiod, struct fieldparams allparams, int diffparams_ind1, int diffparams_ind2) {
	int N_line, i=0;
	int num_params ;
	double varphi=0.0;
	double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*allparams.m0_fieldperiods);
	num_params = allparams.n_diff; 
	N_line = L_fixedpoints*N_gridphi_fieldperiod;
	for (i=0; i<N_line; i++)
	{
		//printf("i=%d/%d\n", i, N_line); printstructposition("lambdamu ", lambdamu[i]);
		//void RK4step_gradRes(double *number, struct position *Xp, struct position *lambdamu, double varphi, double dvarphi, struct field *Bfield_saved, struct fieldparams allparams, int diffparam_ind1, int diffparam_ind2) {
		RK4step_gradRes(number, fieldline[i], lambdamu[i], varphi, dvarphi, Bfield_island[i], allparams, diffparams_ind1, diffparams_ind2);
		//printf("number=%f\n", *number);
		varphi += dvarphi;
	}
	return;
}


