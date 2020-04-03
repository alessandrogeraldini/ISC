/* Author: Alessandro Geraldini */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_sf_bessel.h>
#include "isc.h"

int rho(int fieldaligned_index, int n_turns, int L_fixedpoints) {
	int reordered_index;
	reordered_index = (n_turns*fieldaligned_index)%L_fixedpoints;
	return reordered_index;
	}

int inverserho(int reordered_index, int n_turns, int L_fixedpoints) {
	int fieldaligned_index;
	fieldaligned_index = ( ( (reordered_index%n_turns)*L_fixedpoints + reordered_index ) / n_turns ) %L_fixedpoints;
	return fieldaligned_index;
}

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
	double error=1.0, errorlimit=0.0000000000001, factor=1.0;
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
		//printstruct("fieldline", fieldline);
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
	int N_line=0, L_fixedpoints;
	struct position *centre=malloc(N_line*sizeof(struct position));
	if (tor_mode % m0_symmetry == 0) 
		L_fixedpoints = pol_mode;	
	else 				 
		L_fixedpoints = m0_symmetry*pol_mode;
	N_line = L_fixedpoints*N_gridphi_per_field_period;
	printf("N_line = %d\n", N_line);
	printf("field_periods = %d\n", m0_symmetry);
	do {
		fieldline_start = *fieldline;
		varphi = 0.0;
		//printf("RR=%f\n", fieldline->loc[0]);
		for (i=0; i<N_line; i++)
		{
			centre[i] = *fieldline;
			//printf("where do I seg fault?\n");
			//if (i%N_gridphi_per_field_period==0)
			//{
			//	//printf("varphi = %f\n", varphi);
			//	//printstruct("fieldline\n", fieldline);
			//	//printf("where do I seg fault?\n");
			//}
			//printstruct("fieldline[i]\n", fieldline);
			RK4(fieldline, varphi, dvarphi, coils, n_coils, n_segs);
			varphi += dvarphi;
			//printf("where do I seg fault?\n");
		}
		deltafieldline = addstructs(1.0, fieldline, -1.0, &fieldline_start); 
		//printf("fieldline->tangent[1][1] = %f\n", fieldline->tangent[1][0]);
		inverseTminusI = invert2x2(deltafieldline.tangent, &det_tangent);
		inverseT = invert2x2(fieldline->tangent, &det_tangent);
		//printf("varphi = %f\n", varphi);
		//printstruct("fieldline", fieldline);
		//printf("det(tangent)=%f\n", det_tangent);
//		//printf("tangent[0][0] address is %p\n", &&deltafieldline
		//printmat("inverseTminusI", inverseTminusI, 2, 2);
		
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
		//printmat("matrix", matrix, 2, 2);
		jumptocentre = multiply2x2(inverseTminusI, pdeltafieldline, 1);
		jumptocentre[0][0] *= (factor);
		jumptocentre[1][0] *= (factor);
		//printstruct("deltafieldline", &deltafieldline);
		//printmat("jumptocentre", jumptocentre, 2, 1);
		//printstruct("fieldline_start", &fieldline_start);
		//printstruct("fieldline", fieldline);
		fieldline->loc[0] = fieldline_start.loc[0] - jumptocentre[0][0]; 
		fieldline->loc[1] = fieldline_start.loc[1] - jumptocentre[1][0];
		fieldline->tangent[0][0]=1.0; fieldline->tangent[0][1]=0.0; fieldline->tangent[1][0]=0.0; fieldline->tangent[1][1]=1.0;
		error = sqrt(pow(deltafieldline.loc[0], 2.0) + pow(deltafieldline.loc[1], 2.0)); 
		//printf("error=%f\n", error);
		//free(jumptocentre);
	} while ((deltafieldline.loc[0] > errorlimit) || (deltafieldline.loc[1] > errorlimit));
	//clock_t int3 = clock();
	return centre;
}

struct ext_position *alongcentre(double RR, double ZZ, double *axis, int *n_turns, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode, double ***coils, int *n_coils, int **n_segs) {
	int i=0, q0, L_fixedpoints, clockwise=2, kminus, kplus, rhominus, rhoplus;
	struct position fieldline_start, fieldline, adjfieldline, *centre, *adjcentre;
	double varphi=0.0, omega, *angle_axis;
	double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
	double **inverted, **adjinverted, detcentre, adjdetcentre, **evecs1, *evals1, det=0, trace=0;
	double circumference, ref_angle, **result, **result2;
	double *evals=malloc(2*sizeof(double));
	int N_line=0, main_ind, centre_ind, sec_ind;
	struct ext_position *ext_centre;

	evecs1 = malloc(2*sizeof(double*));
	evecs1[0] = malloc(2*sizeof(double)); evecs1[1] = malloc(2*sizeof(double));
	evals1 = malloc(2*sizeof(double));

	if (tor_mode % m0_symmetry == 0) 
		L_fixedpoints = pol_mode;	
	else 				   
		L_fixedpoints = m0_symmetry*pol_mode;
	N_line = L_fixedpoints*N_gridphi_per_field_period;
	angle_axis = malloc(L_fixedpoints*sizeof(double));
	centre = malloc(L_fixedpoints*sizeof(struct position));
	adjcentre = malloc(L_fixedpoints*sizeof(struct position));
	ext_centre = malloc(L_fixedpoints*sizeof(struct ext_position));
	fieldline_start.loc[0] = RR; fieldline_start.loc[1] = ZZ;
	fieldline_start.tangent = set_identity();
	fieldline = fieldline_start;
	adjfieldline.tangent = set_identity();
	adjfieldline.loc[0] = 0.0; adjfieldline.loc[1] = 0.0;
	centre[0] = fieldline;
	adjcentre[0] = adjfieldline;
	varphi = 0.0;
	ref_angle = atan2(fieldline.loc[1] - axis[1], fieldline.loc[0] - axis[0]);
	angle_axis[0] = 0.0; 
	//printf("varphi = %f\n", varphi);
	//printstruct("fieldline[i]\n", &fieldline);
	for (i=1; i<N_line+1; i++) {
		//RK4(&fieldline, varphi, dvarphi, coils, n_coils, n_segs);
		RK4_adjsimple(&fieldline, &adjfieldline, varphi, dvarphi, coils, n_coils, n_segs);
		varphi += dvarphi;
		if (i%N_gridphi_per_field_period==0) {
			centre_ind = (i / N_gridphi_per_field_period);
			centre[centre_ind % L_fixedpoints] = fieldline;
			adjcentre[centre_ind % L_fixedpoints] = adjfieldline;
			//printf("varphi = %f\n", varphi);
			//printstruct("fieldline[i]\n", &fieldline);
			angle_axis[centre_ind%L_fixedpoints] = atan2(fieldline.loc[1] - axis[1], fieldline.loc[0] - axis[0]) - ref_angle;
			if (angle_axis[centre_ind%L_fixedpoints] < - M_PI) {
				angle_axis[centre_ind%L_fixedpoints] += M_PI;	
			}
			else if (angle_axis[centre_ind%L_fixedpoints] > M_PI) {
				angle_axis[centre_ind%L_fixedpoints] -= M_PI;	
			}
			angle_axis[centre_ind%L_fixedpoints] = atan2(fieldline.loc[1] - axis[1], fieldline.loc[0] - axis[0]) - ref_angle;
			inverted = invert2x2(centre[(centre_ind-1) % L_fixedpoints].tangent, &detcentre);
			adjinverted = invert2x2(adjcentre[(centre_ind-1) % L_fixedpoints].tangent, &adjdetcentre);
			//printmat("inverted", inverted, 2, 2);
			//ptarray[0] = &centre[centre_ind % L_fixedpoints].tangent[0][0];
			//ptarray[1] = &centre[centre_ind % L_fixedpoints].tangent[1][0];
			ext_centre[centre_ind % L_fixedpoints].loc[0] = fieldline.loc[0];
			ext_centre[centre_ind % L_fixedpoints].loc[1] = fieldline.loc[1];
			ext_centre[centre_ind % L_fixedpoints].part_tangent = multiply2x2(centre[centre_ind % L_fixedpoints].tangent, inverted, 2);
			ext_centre[centre_ind % L_fixedpoints].adj_part_tangent = multiply2x2(adjcentre[centre_ind % L_fixedpoints].tangent, adjinverted, 2);
			//printmat("ext_centre.part_tangent", ext_centre[centre_ind % L_fixedpoints].part_tangent, 2, 2);
		}
	}
	if (angle_axis[1] > 0.0) {
		if ( (angle_axis[2] > angle_axis[1]) || (angle_axis[2] < 0.0) ) clockwise = 0;
		else  clockwise = 1;
	}
	else {
		if ( (angle_axis[2] < angle_axis[1]) || (angle_axis[2] > 0.0) ) clockwise = 1;
		else clockwise = 0;
	}
	//inverted = invert2x2(centre[(centre_ind-1) % L_fixedpoints].tangent, &detcentre);
	circumference = 0.0;
	//rearr_index[0] = 0;
	for (centre_ind=0;centre_ind<L_fixedpoints;centre_ind++) {
		if (centre_ind != 0) {
			if ( (clockwise == 0) && ( ( (angle_axis[centre_ind] > small) && (angle_axis[centre_ind-1] < - small) ) ) ) (*n_turns) +=1;	
			else if ( (clockwise == 1) && ( ( (angle_axis[centre_ind] < -small) && (angle_axis[centre_ind-1] > small) ) ) ) (*n_turns) +=1;	
		}
		angle_axis[centre_ind] = atan2(ext_centre[(centre_ind+1)%L_fixedpoints].loc[1]-axis[1], ext_centre[(centre_ind+1)%L_fixedpoints].loc[0]-axis[0]);
		circumference += sqrt(pow(ext_centre[(centre_ind+1)%L_fixedpoints].loc[0] 
					        - ext_centre[centre_ind].loc[0], 2.0)
				    	    + pow(ext_centre[(centre_ind+1)%L_fixedpoints].loc[1] 
						- ext_centre[centre_ind].loc[1], 2.0));
		ext_centre[centre_ind].epar = malloc(2*sizeof(double)); ext_centre[centre_ind].eperp = malloc(2*sizeof(double));
		ext_centre[centre_ind].full_tangent = multiply2x2(ext_centre[(centre_ind +2)%L_fixedpoints].part_tangent, ext_centre[(centre_ind+1)%L_fixedpoints].part_tangent, 2);
		ext_centre[centre_ind].adj_full_tangent = multiply2x2(ext_centre[(centre_ind +2)%L_fixedpoints].adj_part_tangent, ext_centre[(centre_ind+1)%L_fixedpoints].adj_part_tangent, 2);
		for (sec_ind=3;sec_ind<L_fixedpoints+1;sec_ind++) {
			multiply2x2reassign(ext_centre[(centre_ind + sec_ind)%L_fixedpoints].part_tangent, ext_centre[centre_ind].full_tangent, 2);
			multiply2x2reassign(ext_centre[(centre_ind + sec_ind)%L_fixedpoints].adj_part_tangent, ext_centre[centre_ind].adj_full_tangent, 2);
		}
		//printmat("ext_centre.fulltangent", ext_centre[centre_ind].full_tangent, 2, 2);
		//printmat("ext_centre._adj_fulltangent", ext_centre[centre_ind].adj_full_tangent, 2, 2);
		symmeigs(ext_centre[centre_ind].full_tangent, ext_centre[centre_ind].eperp, ext_centre[centre_ind].epar, evals);
		if (ext_centre[centre_ind].eperp[0]*(ext_centre[centre_ind].loc[0] - axis[0]) + ext_centre[centre_ind].eperp[1]*(ext_centre[centre_ind].loc[1] - axis[1]) < 0.0) {
			ext_centre[centre_ind].eperp[0] *= (-1);
			ext_centre[centre_ind].eperp[1] *= (-1);
			ext_centre[centre_ind].epar[0] *= (-1);
			ext_centre[centre_ind].epar[1] *= (-1);
			// this fix will not work in general. You can imagine weird shapes where it wouldn't
		}
		//printf("eigenvectors are (%f, %f) and (%f, %f)\n", ext_centre[centre_ind].eperp[0], ext_centre[centre_ind].eperp[1], ext_centre[centre_ind].epar[0], ext_centre[centre_ind].epar[1]);
		printf("position is (%f, %f)\n", ext_centre[centre_ind].loc[0], ext_centre[centre_ind].loc[1]);
	}
	//printf("n_turns = %d\n", *n_turns);
	for (main_ind=0;main_ind<L_fixedpoints; main_ind++) {
		linalg2x2(ext_centre[main_ind].full_tangent, evecs1, evals1, &det, &trace);
		ext_centre[main_ind].circumference = circumference;
		//rhominus = ( L_fixedpoints + (n_turns*main_ind)%L_fixedpoints - 1) %L_fixedpoints;
		rhominus = ( L_fixedpoints + ((*n_turns)*main_ind)%L_fixedpoints - 1) %L_fixedpoints;
		rhoplus =  ((*n_turns)*main_ind)%L_fixedpoints + 1;
		//kminus = ( ( (rhominus%n_turns)*L_fixedpoints + rhominus ) / n_turns ) %L_fixedpoints;
		kminus = inverserho(rhominus, (*n_turns), L_fixedpoints);
		kplus = inverserho(rhoplus, (*n_turns), L_fixedpoints);
		ext_centre[main_ind].chord[0] = ext_centre[main_ind].loc[0] - ext_centre[kminus].loc[0];
		ext_centre[main_ind].chord[1] = ext_centre[main_ind].loc[1] - ext_centre[kminus].loc[1];
		ext_centre[main_ind].chordplus[0] = ext_centre[kplus].loc[0] - ext_centre[main_ind].loc[0];
		ext_centre[main_ind].chordplus[1] = ext_centre[kplus].loc[1] - ext_centre[main_ind].loc[1];
		if ((fabs(evecs1[0][0]) < small) || (fabs(evecs1[0][1]) < small)) {
			ext_centre[main_ind].angle = evals1[1];
			omega = m0_symmetry*evals1[1]/(2.0*M_PI*L_fixedpoints); // I have changed it compared to Cary and Hanson's paper
			q0 = (int) (m0_symmetry/(4.0*omega) - L_fixedpoints/2.0 + 0.5);
			ext_centre[main_ind].q0_index = q0;
		}
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
			result = multiply2x2(ext_centre[main_ind].full_tangent, ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints], 2); 
			result2 = multiply2x2(ext_centre[main_ind].full_tangent, result, 2); 
			for (sec_ind=0; sec_ind<(q0+centre_ind)%L_fixedpoints; sec_ind++) {
				multiply2x2reassign(ext_centre[(sec_ind+1+main_ind)%L_fixedpoints].part_tangent, ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints], 2); 
			}
			//result = multiply2x2(ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints], &(ext_centre[main_ind].eperp), 2);
			ext_centre[main_ind].sperp[(main_ind+q0+centre_ind)%L_fixedpoints][0] = ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints][0][0]*ext_centre[main_ind].eperp[0] + ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints][0][1]*ext_centre[main_ind].eperp[1];
			ext_centre[main_ind].sperp[(main_ind+q0+centre_ind)%L_fixedpoints][1] = ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints][1][0]*ext_centre[main_ind].eperp[0] + ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints][1][1]*ext_centre[main_ind].eperp[1];
		}
		//ext_centre[main_ind].sperp_final[0] = result[0][0]*ext_centre[main_ind].eperp[0] + result[0][1]*ext_centre[main_ind].eperp[1];
		//ext_centre[main_ind].sperp_final[1] = result[1][0]*ext_centre[main_ind].eperp[0] + result[1][1]*ext_centre[main_ind].eperp[1];
		//sperp_final is probably unnecessary
		//ext_centre[main_ind].sperp_final[0] = result2[0][0]*ext_centre[main_ind].eperp[0] + result2[0][1]*ext_centre[main_ind].eperp[1];
		//ext_centre[main_ind].sperp_final[1] = result2[1][0]*ext_centre[main_ind].eperp[0] + result2[1][1]*ext_centre[main_ind].eperp[1];
	}
	//ext_centre[centre_index % L_fixedpoints].part_tangent = multiply2x2(centre[centre_index % L_fixedpoints].tangent, inverted, 2);
	//clock_t int3 = clock();
	return ext_centre;
} 

//double *islandwidthnew(struct ext_position *ext_fieldline, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode) {
double *islandwidthnew(double ***coils, int *n_coils, int **n_segs, struct ext_position *ext_centre, struct position *lambda_circ, struct position **lambda_tangent, struct position **mu_tangent, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode) {
	// declarations
	//clock_t start = clock();
	int i=0, q0_ind, Q0_ind, Q_ind, q_ind, main_ind, sec_ind, N_line=0, L_fixedpoints;
	struct position centre, lambda, mu, sperp;
	double varphi=0.0, chordlength, chordpluslength, *vecresult;
	double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
	double *gradcircumference, *gradtangent, initialgradcirc=0.0, *gradwidth;
	double *width, sum_matrix_elements, matrix_element;
	double circumference2=ext_centre[0].circumference, circumference;
	double matrix_element_oldway, sum_matrix_elements_oldway=0.0;
	printf("tor_mode=%d\n", tor_mode);
	if (tor_mode % m0_symmetry == 0) 
		L_fixedpoints = pol_mode;	
	else 				 
		L_fixedpoints = m0_symmetry*pol_mode;
	printf("L_fixedpoints=%d\n", L_fixedpoints);
	width = malloc(L_fixedpoints*sizeof(double));
	for (main_ind=0;main_ind<L_fixedpoints; main_ind++) {
		sum_matrix_elements = sum_matrix_elements_oldway = circumference = 0.0;		
		for (sec_ind=0;sec_ind<L_fixedpoints; sec_ind++) {
			circumference += pow( pow(ext_centre[main_ind].chord[0], 2.0) + pow(ext_centre[main_ind].chord[1], 2.0) , 0.5 );
			matrix_element_oldway = inner(ext_centre[sec_ind].epar, 
					       ext_centre[main_ind].long_tangent[sec_ind], 
					       ext_centre[main_ind].eperp);	
			matrix_element = ext_centre[sec_ind].epar[0]*ext_centre[main_ind].sperp[sec_ind][0] + ext_centre[sec_ind].epar[1]*ext_centre[main_ind].sperp[sec_ind][1];
			sum_matrix_elements_oldway += fabs(matrix_element_oldway);
			//sum_matrix_elements += fabs(matrix_element);
			//printf("matrix_element=%f\n", matrix_element);
			sum_matrix_elements += matrix_element;
			//printf("sum_matrix_element=%f\n", sum_matrix_elements);
			//sum_matrix_elements += fabs(matrix_element)/sin(ext_centre[main_ind].angle*ext_centre[main_ind].q0_index/L_fixedpoints); //This has the sin piece in the denominator, which can be safely set to unity
			//circumference += sqrt(pow(ext_centre[(sec_ind+1)%L_fixedpoints].loc[0] 
			//		        - ext_centre[sec_ind].loc[0], 2.0)
			//	    	    + pow(ext_centre[(sec_ind+1)%L_fixedpoints].loc[1] 
			//			- ext_centre[sec_ind].loc[1], 2.0));
			//printf("circumference = %f\nmatrix_element = %f\n", circumference, matrix_element);
			//printf("other circumference = %f\n", ext_centre[sec_ind].circumference);
		}
		//circumference = 2.0*M_PI*0.2102;
		//circumference = 2.0*M_PI*0.2094;
		printf("circumference = %f\n", circumference);
		//printf("circumference_oldway = %f\n", circumference2);
		printf("sum_matrix_elements = %f\n", sum_matrix_elements);
		//printf("sum_matrix_elements_oldway = %f\n", sum_matrix_elements_oldway);
		width[main_ind] = 2.0*L_fixedpoints*circumference/(M_PI*pol_mode*sum_matrix_elements);
		printf("width = %f for index= %d\n", width[main_ind], main_ind);
	}

	gradtangent = malloc(L_fixedpoints*sizeof(double));
	gradwidth = malloc(L_fixedpoints*sizeof(double));
	centre.loc[0] = ext_centre[0].loc[0]; centre.loc[1] = ext_centre[0].loc[1];
	//lambda_circ.loc[0] = lambda[0]; lambda_circ.loc[1] = lambda[1];
	centre.tangent = set_identity();
	//lambda_circ.tangent = set_identity();
	gradcircumference = &initialgradcirc;
	N_line = L_fixedpoints*N_gridphi_per_field_period;
	varphi = 0.0;
	for (i=0; i<N_line; i++)
	{
		if (i%N_gridphi_per_field_period==0)
		{
			main_ind = i/N_gridphi_per_field_period;
			chordlength= pow( pow(ext_centre[main_ind].chord[0], 2.0) + pow(ext_centre[main_ind].chord[1], 2.0), 0.5);
			chordpluslength= pow( pow(ext_centre[main_ind].chordplus[0], 2.0) + pow(ext_centre[main_ind].chordplus[1], 2.0), 0.5);
			lambda_circ->loc[0] += (ext_centre[main_ind].chord[0]/chordlength - ext_centre[main_ind].chordplus[0]/chordpluslength);
			lambda_circ->loc[1] += (ext_centre[main_ind].chord[1]/chordlength - ext_centre[main_ind].chordplus[1]/chordpluslength);
		}
		RK4_adjgrad(gradcircumference, &centre, lambda_circ, varphi, dvarphi, coils, n_coils, n_segs);
		varphi += dvarphi;
		//printf("gradcircumference = %f\n", *gradcircumference);
	}
	//printf("varphi = %f\n", varphi);
	//printstruct("lambda",lambda_circ);
	//printstruct("Xp",&centre);
	printf("gradcircumference = %f\n", *gradcircumference);

	centre.loc[0] = ext_centre[0].loc[0]; centre.loc[1] = ext_centre[0].loc[1];
	centre.tangent[0][0] = centre.tangent[1][1] = 1.0;
	centre.tangent[0][1] = centre.tangent[1][0] = 0.0;
	//centre.loc[0] = ext_centre[0].loc[0]; centre.loc[1] = ext_centre[0].loc[1];
	////lambda_centre.loc[0] = lambda[0]; lambda_centre.loc[1] = lambda[1];
	//centre.tangent = set_identity();
	//lambda_centre.tangent = set_identity();
	q0_ind = ext_centre[0].q0_index; 
	Q0_ind = q0_ind/L_fixedpoints + 2;
	sperp.tangent = set_identity();
	for (main_ind=0; main_ind<L_fixedpoints; main_ind++) {
		varphi = main_ind*2.0*M_PI/m0_symmetry;
		gradtangent[main_ind] = 0.0; 
		centre.loc[0] = ext_centre[main_ind].loc[0]; centre.loc[1] = ext_centre[main_ind].loc[1];
		centre.tangent[0][0] = centre.tangent[1][1] = 1.0;
		centre.tangent[0][1] = centre.tangent[1][0] = 0.0;
		vecresult = multiply(centre.tangent, ext_centre[main_ind].eperp);
		sperp.loc[0] = vecresult[0]; sperp.loc[1] = vecresult[1];
		free(vecresult);
		for (Q_ind=0;Q_ind<Q0_ind;Q_ind++) {
			lambda = lambda_tangent[main_ind][Q_ind];
			mu = mu_tangent[main_ind][Q_ind];
			//printstruct("lambda", &lambda);
			//printstruct("mu", &mu);

			centre.loc[0] = ext_centre[main_ind].loc[0]; centre.loc[1] = ext_centre[main_ind].loc[1];
			vecresult = multiply(centre.tangent, ext_centre[main_ind].eperp);
			sperp.loc[0] = vecresult[0]; sperp.loc[1] = vecresult[1];
			free(vecresult);

			//printf("sperp=%f\n", sperp.tangent[1][1]);
			//centre->tangent[0][0] = 1.0 ; centre->tangent[0][1] = 0.0;
			//centre->tangent[1][0] = 0.0 ; centre->tangent[0][1] = 1.0;
			for (i=0; i<N_line; i++)
			{
				if (i%N_gridphi_per_field_period==0)
				{
					sec_ind = i/N_gridphi_per_field_period;
					//printf("varphi = %f\n", varphi);
					//printstruct("Xp", &centre);
					//printf("gradtangent = %f\n", gradtangent[sec_ind]);
					//if (i/N_gridphi_per_field_period==0) lambda_centre->loc[0] += 1.0;
					//printstruct("lambda\n", lambda);
					//if (sec_ind==0) 
					//{
					//	printf("main_ind=%d, Q_ind=%d, sec_ind=%d\n", main_ind, Q_ind, sec_ind); 
					//	printf("gradtangent = %f\n", gradtangent[main_ind]);
					//	printstruct("mu", &mu);
					//	printstruct("lambda", &lambda);
					//}
					q_ind = Q_ind*L_fixedpoints + sec_ind;
					//check if makes sense
					if ( q_ind >= q0_ind && q_ind < q0_ind + L_fixedpoints) {
						mu.loc[0] += (ext_centre[(main_ind+sec_ind)%L_fixedpoints].epar[0]);
						mu.loc[1] += (ext_centre[(main_ind+sec_ind)%L_fixedpoints].epar[1]);
					}
				}
				//RK4_adjtangent(&centre, &lambda, &sperp, &mu, varphi, dvarphi, coils, n_coils, n_segs);
				//RK4_adjgrad(gradtangent, &centre, &lambda, varphi, dvarphi, coils, n_coils, n_segs);
				RK4_adjgradtangent(gradtangent+main_ind, &centre, &lambda, &sperp, &mu, varphi, dvarphi, coils, n_coils, n_segs);
				varphi += dvarphi;
			}
		}
		//printf("varphi = %f\n", varphi);
		//printstruct("lambda", &lambda);
		//printstruct("Xp",&centre);
	}
	for (main_ind=0; main_ind<L_fixedpoints; main_ind++) {
		gradwidth[main_ind] = width[main_ind]*((*gradcircumference)/circumference - gradtangent[main_ind]/sum_matrix_elements);
		printf("gradcirc = %f\n", *gradcircumference);
		printf("gradtangent = %f\n", gradtangent[main_ind]);
		printf("gradwidth = %f\n", gradwidth[main_ind]);
	}

	return width;
}

//double *islandwidth(struct ext_position *ext_fieldline, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode) {
//	int main_index, centre_index, L_fixedpoints;
//	double *wperp, sum_matrix_elements, matrix_element;
//	double circumference2=ext_fieldline[0].circumference, circumference;
//	double matrix_element_oldway, sum_matrix_elements_oldway=0.0;
//	printf("tor_mode=%d\n", tor_mode);
//	if (tor_mode % m0_symmetry == 0) L_fixedpoints = pol_mode;	
//	else 				   L_fixedpoints = m0_symmetry*pol_mode;
//	printf("L_fixedpoints=%d\n", L_fixedpoints);
//	wperp = malloc(L_fixedpoints*sizeof(double));
//	for (main_index=0;main_index<L_fixedpoints; main_index++) {
//		sum_matrix_elements = sum_matrix_elements_oldway = circumference = 0.0;		
//		for (centre_index=0;centre_index<L_fixedpoints; centre_index++) {
//			circumference += pow( pow(ext_fieldline[main_index].chord[0], 2.0) + pow(ext_fieldline[main_index].chord[1], 2.0) , 0.5 );
//			matrix_element_oldway = inner(ext_fieldline[centre_index].epar, 
//					       ext_fieldline[main_index].long_tangent[centre_index], 
//					       ext_fieldline[main_index].eperp);	
//			matrix_element = ext_fieldline[centre_index].epar[0]*ext_fieldline[main_index].sperp[centre_index][0] + ext_fieldline[centre_index].epar[1]*ext_fieldline[main_index].sperp[centre_index][1];
//			sum_matrix_elements_oldway += fabs(matrix_element_oldway);
//			//sum_matrix_elements += fabs(matrix_element);
//			//printf("matrix_element=%f\n", matrix_element);
//			sum_matrix_elements += matrix_element;
//			//printf("sum_matrix_element=%f\n", sum_matrix_elements);
//			//sum_matrix_elements += fabs(matrix_element)/sin(ext_fieldline[main_index].angle*ext_fieldline[main_index].q0_index/L_fixedpoints); //This has the sin piece in the denominator, which can be safely set to unity
//			//circumference += sqrt(pow(ext_fieldline[(centre_index+1)%L_fixedpoints].loc[0] 
//			//		        - ext_fieldline[centre_index].loc[0], 2.0)
//			//	    	    + pow(ext_fieldline[(centre_index+1)%L_fixedpoints].loc[1] 
//			//			- ext_fieldline[centre_index].loc[1], 2.0));
//			//printf("circumference = %f\nmatrix_element = %f\n", circumference, matrix_element);
//			//printf("other circumference = %f\n", ext_fieldline[centre_index].circumference);
//		}
//		//circumference = 2.0*M_PI*0.2102;
//		//circumference = 2.0*M_PI*0.2094;
//		printf("circumference = %f\n", circumference);
//		//printf("circumference_oldway = %f\n", circumference2);
//		printf("sum_matrix_elements = %f\n", sum_matrix_elements);
//		//printf("sum_matrix_elements_oldway = %f\n", sum_matrix_elements_oldway);
//		wperp[main_index] = 2.0*L_fixedpoints*circumference/(M_PI*pol_mode*sum_matrix_elements);
//		printf("width = %f for index= %d\n", wperp[main_index], main_index);
//	}
//	return wperp;
//}

struct position findadjsimple(double ***coils, int *n_coils, int **n_segs, struct ext_position *ext_centre, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode) {
	// declarations
	//clock_t start = clock();
	int centre_ind; //i=0
	struct position lambda, lambda_start, deltalambda; //, fieldline_start, fieldline;
	double varphi=0.0, **jumptocentre=calloc(1,sizeof(double*)); //, **testunity;
	//double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
	double **matrix, *vecresult, **result;
	double **pdeltalambda, det_tangent;
	double **inverseTminusI, **inverseT;
	double error = 1.0, errorlimit = 0.000000000001;
	double factor = 1.0, chordlength, chordpluslength;
	int N_line=0, L_fixedpoints;
	if (tor_mode % m0_symmetry == 0) 
		L_fixedpoints = pol_mode;	
	else 				 
		L_fixedpoints = m0_symmetry*pol_mode;
	N_line = L_fixedpoints*N_gridphi_per_field_period;
	printf("N_line = %d\n", N_line);
	printf("field_periods = %d\n", m0_symmetry);
	lambda.loc[0] = 1.0; lambda.loc[1] = 1.0;
	lambda.tangent = set_identity();
	//fieldline_start.loc[0] = ext_centre[0].loc[0]; fieldline_start.loc[1] = ext_centre[0].loc[1];
	//fieldline.tangent = set_identity();
	do {
		//fieldline.loc[0] = fieldline_start.loc[0];
		//fieldline.loc[1] = fieldline_start.loc[1];
		lambda_start = lambda;
		varphi = 0.0;
		for (centre_ind=0; centre_ind<L_fixedpoints; centre_ind++) {
				chordlength= pow( pow(ext_centre[centre_ind].chord[0], 2.0) + pow(ext_centre[centre_ind].chord[1], 2.0), 0.5);
				chordpluslength= pow( pow(ext_centre[centre_ind].chordplus[0], 2.0) + pow(ext_centre[centre_ind].chordplus[1], 2.0), 0.5);
				lambda.loc[0] += (ext_centre[centre_ind].chord[0]/chordlength - ext_centre[centre_ind].chordplus[0]/chordpluslength);
				lambda.loc[1] += (ext_centre[centre_ind].chord[1]/chordlength - ext_centre[centre_ind].chordplus[1]/chordpluslength);
				vecresult = multiply(ext_centre[(centre_ind+1)%L_fixedpoints].adj_part_tangent, lambda.loc);
				result = multiply2x2(ext_centre[(centre_ind+1)%L_fixedpoints].adj_part_tangent, lambda.tangent, 2);
				lambda.loc[0] = vecresult[0]; lambda.loc[1] = vecresult[1];
				lambda.tangent = result;
				//printf("lambda = (%f, %f)\n", lambda->loc[0], lambda->loc[1]);
				//printmat("lambda start\n", lambda_start.tangent, 2, 2);
		}
		// temporary (the above is much faster and now works, so I can permanently delete what's below eventually)
		//for (i=0; i<N_line; i++)
		//{
		//	if (i%N_gridphi_per_field_period==0)
		//	{
		//		centre_ind = i/N_gridphi_per_field_period;
		//		printf("varphi = %f\n", varphi);
		//		//printstruct("lambda\n", lambda);
		//		//if (centre_ind==0) lambda->loc[0] += 1.0;
		//		chordlength= pow( pow(ext_centre[centre_ind].chord[0], 2.0) + pow(ext_centre[centre_ind].chord[1], 2.0), 0.5);
		//		chordpluslength= pow( pow(ext_centre[centre_ind].chordplus[0], 2.0) + pow(ext_centre[centre_ind].chordplus[1], 2.0), 0.5);
		//		lambda->loc[0] += (ext_centre[centre_ind].chord[0]/chordlength - ext_centre[centre_ind].chordplus[0]/chordpluslength);
		//		lambda->loc[1] += (ext_centre[centre_ind].chord[1]/chordlength - ext_centre[centre_ind].chordplus[1]/chordpluslength);
		//	}
		//	RK4_adjsimple(&fieldline, lambda, varphi, dvarphi, coils, n_coils, n_segs);
		//	varphi += dvarphi;
		//}
		deltalambda = addstructs(1.0, &lambda, -1.0, &lambda_start); 
		//printstruct("lambda\n", &lambda);
		//printstruct("deltalambda", &deltalambda);
		//printf("fieldline->tangent[1][1] = %f\n", fieldline->tangent[1][0]);
		inverseTminusI = invert2x2(deltalambda.tangent, &det_tangent);
		inverseT = invert2x2(lambda.tangent, &det_tangent);
		//printf("varphi = %f\n", varphi);
		//printstruct("fieldline", lambda);
		//printf("det(tangent)=%f\n", det_tangent);
//		printf("tangent[0][0] address is %p\n", &&deltalambda
		//printmat("inverseTminusI", inverseTminusI, 2, 2);
		
		pdeltalambda = calloc(2,sizeof(double*)); 
		*pdeltalambda     = &deltalambda.loc[0];
		*(pdeltalambda+1) = &deltalambda.loc[1];
		matrix = multiply2x2(inverseTminusI, deltalambda.tangent, 2);
		//printmat("matrix", matrix, 2, 2);
		jumptocentre = multiply2x2(inverseTminusI, pdeltalambda, 1);
		jumptocentre[0][0] *= (factor);
		jumptocentre[1][0] *= (factor);
		//printstruct("deltalambda", &deltalambda);
		//printmat("jumptocentre", jumptocentre, 2, 1);
		//printstruct("lambda_start", &lambda_start);
		//printstruct("lambda", &lambda);
		lambda.loc[0] = lambda_start.loc[0] - jumptocentre[0][0]; 
		lambda.loc[1] = lambda_start.loc[1] - jumptocentre[1][0];
		lambda.tangent[0][0] = 1.0; lambda.tangent[0][1]=0.0; lambda.tangent[1][0]=0.0; lambda.tangent[1][1]=1.0;
		error = sqrt(pow(deltalambda.loc[0], 2.0) + pow(deltalambda.loc[1], 2.0)); 
		//printf("%f\n", deltalambda.loc[0]);
		//printf("%f\n", deltalambda.loc[1]);
		//free(jumptocentre);
	} while ((fabs(deltalambda.loc[0]) > errorlimit) || (fabs(deltalambda.loc[1]) > errorlimit));
	//clock_t int3 = clock();
	return lambda;
}

struct position **findadjmu(double ***coils, int *n_coils, int **n_segs, struct ext_position *ext_centre, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode) {
	// declarations
	//clock_t start = clock();
	int centre_ind, main_ind, q0_ind, Q0_ind, Q_ind, q_ind;
	struct position mu_start, fieldline_start, fieldline, **mu;
	//double **jumptocentre=calloc(1,sizeof(double*)); //, **testunity;
	//double varphi = 0.0, dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
	double **result, *vecresult, **inverted, det;// varphi=0.0;
	int N_line=0, L_fixedpoints;
	if (tor_mode % m0_symmetry == 0) 
		L_fixedpoints = pol_mode;	
	else 				 
		L_fixedpoints = m0_symmetry*pol_mode;
	mu = malloc(L_fixedpoints*sizeof(struct position*));
	mu_start.tangent = set_identity();
	for (main_ind=0; main_ind<L_fixedpoints; main_ind++)  {
		q0_ind = ext_centre[main_ind].q0_index; 
		Q0_ind = q0_ind/L_fixedpoints + 2;
		mu[main_ind] = malloc((Q0_ind)*sizeof(struct position));
		mu_start.loc[0] = 0.0;  mu_start.loc[1] = 0.0; 
		N_line = L_fixedpoints*N_gridphi_per_field_period;
		//printf("N_line = %d\n", N_line);
		//printf("field_periods = %d\n", m0_symmetry);
		fieldline_start.loc[0] = ext_centre[main_ind].loc[0]; fieldline_start.loc[1] = ext_centre[main_ind].loc[1];
		fieldline.tangent = set_identity();
		fieldline.loc[0] = fieldline_start.loc[0]; fieldline.loc[1] = fieldline_start.loc[1];
		//printf("RR=%f\n", fieldline.loc[0]);
		for (Q_ind=Q0_ind-1; Q_ind>=0; Q_ind--) {
			mu[main_ind][Q_ind] = mu_start; 
			for (centre_ind=L_fixedpoints-1; centre_ind>=0; centre_ind--) {
				inverted = invert2x2(ext_centre[(main_ind+centre_ind+1)%L_fixedpoints].adj_part_tangent, &det);
				vecresult = multiply(inverted, mu[main_ind][Q_ind].loc);
				result = multiply2x2(inverted, mu[main_ind][Q_ind].tangent, 2);
				mu[main_ind][Q_ind].loc[0] = vecresult[0]; mu[main_ind][Q_ind].loc[1] = vecresult[1];
				mu[main_ind][Q_ind].tangent = result;
				q_ind = Q_ind*L_fixedpoints + centre_ind;
				if ( ( q_ind >= q0_ind ) && ( q_ind < q0_ind + L_fixedpoints ) ) {
					//mu[main_ind][Q_ind].loc[0] = 0.0;
					//mu[main_ind][Q_ind].loc[1] = 0.0;
					mu[main_ind][Q_ind].loc[0] -= (ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[0]);
					mu[main_ind][Q_ind].loc[1] -= (ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[1]);
					//printf("epar[0] = %f\n", ext_centre[centre_ind].epar[0]);
					//printf("epar[1] = %f\n", ext_centre[centre_ind].epar[1]);
				}
				//printf("main_ind=%d, Q_ind=%d, centre_ind=%d\n", main_ind, Q_ind, centre_ind);
				//printstruct("mu", mu[main_ind]+Q_ind);
			}
			mu_start = mu[main_ind][Q_ind];
			//printf("Q_ind=%d\n", Q_ind);
			//printstruct("mu\n", mu[main_ind]+Q_ind);
			//printstruct("mu_start\n", &mu_start);

			//for (i=N_line-1; i>=0; i--)
			//{
			//	RK4_adjsimple(&fieldline, mu[main_ind]+Q_ind, varphi, -dvarphi, coils, n_coils, n_segs);
			//	if (i%N_gridphi_per_field_period==0)
			//	{
			//		centre_ind = i/N_gridphi_per_field_period;
			//		//printstruct("mu\n", mu);
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
			//			printstruct("mu", mu[main_ind]+Q_ind);
			//			//printstruct("mu +", mu[main_ind]+Q_ind+1);
			//		}
			//	}
			//	varphi -= dvarphi;
			//}
			//mu_start = mu[main_ind][Q_ind];
		}
		//clock_t int3 = clock();
	}
	printf("exit findadjmu module \n");
	return mu;
}

struct position **findadjtangent(double ***coils, int *n_coils, int **n_segs, struct ext_position *ext_centre, struct position **mu, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode) {
	// declarations
	//clock_t start = clock();
	int i=0, centre_ind, main_ind, q0_ind, Q0_ind, Q_ind, q_ind, ind, count = 0;
	struct position fieldline_start, fieldline, sperp, muvar, lambdavar, deltalambda, **lambda, lambdain; //, fieldline_start, fieldline;
	double varphi=0.0, **jumptocentre=calloc(1,sizeof(double*)); //, **testunity;
	//double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
	double **pdeltalambda, det_tangent;
	double **inverseTminusI, **inverseT;
	double error = 1.0, errorlimit = 0.000000001;
	double **matrix, **matrix2, *vector;
	double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry), factor = 1.0;
	int N_line=0, L_fixedpoints;
	if (tor_mode % m0_symmetry == 0) 
		L_fixedpoints = pol_mode;	
	else 				 
		L_fixedpoints = m0_symmetry*pol_mode;
	N_line = L_fixedpoints*N_gridphi_per_field_period;
	printf("N_line = %d\n", N_line);
	printf("field_periods = %d\n", m0_symmetry);
	lambda = malloc(L_fixedpoints*sizeof(struct position*)); 
	matrix2 = set_identity();
	sperp.tangent = set_identity();
	lambdavar.tangent = set_identity();
	//for (main_ind=0;main_ind<1;main_ind++) {
	q0_ind = ext_centre[0].q0_index; 
	Q0_ind = q0_ind/L_fixedpoints + 2;
	printf("Q0_ind=%d\n", Q0_ind);
	for (main_ind=0;main_ind<L_fixedpoints;main_ind++) {
		lambda[main_ind] = malloc(Q0_ind*sizeof(struct position));
		fieldline_start.loc[0] = ext_centre[main_ind].loc[0]; fieldline_start.loc[1] = ext_centre[main_ind].loc[1];
		fieldline.tangent = set_identity();
		lambdavar.loc[0] = 0.0; lambdavar.loc[1] = 0.0;
		//sperp.loc[0] = ext_centre[main_ind].eperp[0]; sperp.loc[1] = ext_centre[main_ind].eperp[1];
		for (Q_ind=0;Q_ind<Q0_ind;Q_ind++) {
			lambda[main_ind][Q_ind].loc[0] = 0.0; lambda[main_ind][Q_ind].loc[1] = 0.0;
			lambda[main_ind][Q_ind].tangent = set_identity();
			lambdain.loc[0] = 0.0; lambdain.loc[1] = 0.0;
			lambdain.tangent = set_identity();
			//fieldline_start.loc[0] = ext_centre[0].loc[0]; fieldline_start.loc[1] = ext_centre[0].loc[1];
			//fieldline.tangent = set_identity();
			//printf("RR=%f\n", fieldline.loc[0]);
			count = 0;
			do {
				//lambdavar.loc[0] = lambda[main_ind][Q_ind].loc[0]; 
				//lambdavar.loc[1] = lambda[main_ind][Q_ind].loc[1]; 
				lambdavar.loc[0] = lambda[main_ind][Q_ind].loc[0]; 
				lambdavar.loc[1] = lambda[main_ind][Q_ind].loc[1]; 
				matrix2[0][0] = 1.0;
				matrix2[0][1] = 0.0;
				matrix2[1][0] = 0.0;
				matrix2[1][1] = 1.0;
				for (ind=0;ind<Q_ind; ind++) {
					multiply2x2reassign(ext_centre[main_ind].full_tangent, matrix2, 2);
				}
				vector = multiply(matrix2, ext_centre[main_ind].eperp);
				sperp.loc[0] = vector[0];  sperp.loc[1] = vector[1];
				sperp.tangent[0][0] = sperp.tangent[1][1] = 1.0; 
				sperp.tangent[1][0] = sperp.tangent[0][1] = 0.0; 
				//printstruct("lambda", lambda[main_ind]+Q_ind);
				//lambda[main_ind][Q_ind].tangent[0][0] = 1.0; lambda[main_ind][Q_ind].tangent[1][0] = 0.0; 
				//lambda[main_ind][Q_ind].tangent[0][1] = 0.0; lambda[main_ind][Q_ind].tangent[1][1] = 1.0; 
				fieldline.loc[0] = fieldline_start.loc[0];
				fieldline.loc[1] = fieldline_start.loc[1];
				fieldline.tangent[0][0] = fieldline.tangent[1][1] = 1.0; 
				fieldline.tangent[1][0] = fieldline.tangent[0][1] = 0.0; 
				//printstruct("fieldline", &fieldline);
				//printstruct("lambdain", &lambdain);
				muvar = mu[main_ind][Q_ind];
				//varphi = 2.0*M_PI*(Q0_ind+1)*L_fixedpoints/m0_symmetry;
				varphi = main_ind*2.0*M_PI/m0_symmetry;
				for (i=0; i<N_line; i++) {
					if (i%N_gridphi_per_field_period==0) {
						//printf("N_line = %d\n", N_line);
						//printf("count = %d\n", count);
						//printf("varphi = %f\n", varphi);
						//printstruct("fieldline", &fieldline);
						//printstruct("lambdain", &lambdain);
						//printstruct("sperp", &sperp);
						centre_ind = i/N_gridphi_per_field_period;
						q_ind = Q_ind*L_fixedpoints + centre_ind;
						
						//if (Q_ind == Q0_ind-1) {
						//	printf("main_ind=%d, Q_ind=%d, centre_ind=%d\n", main_ind, Q_ind, centre_ind);
						//	printf("count = %d\n", count);
						//	printstruct("muvar", &muvar);
						//	printstruct("lambdain", &lambdain);
						//}
						if ( q_ind >= q0_ind && q_ind < q0_ind + L_fixedpoints) {
							muvar.loc[0] += (ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[0]);
							muvar.loc[1] += (ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[1]);
						}
					}
					RK4_adjtangent(&fieldline, &lambdain, &sperp, &muvar, varphi, dvarphi, coils, n_coils, n_segs);
						//printstruct("lambda\n", lambda);
						//if (centre_ind==0) lambda->loc[0] += 1.0;
				//	printstruct("lambda", &lambdain);
					//RK4_adjtangent(&fieldline, lambda[main_ind]+Q_ind, &sperp, &muvar, varphi, dvarphi, coils, n_coils, n_segs);
					//RK4_adjsimple(&fieldline, &lambdain, varphi, dvarphi, coils, n_coils, n_segs);
					//RK4_adjsimple(&fieldline, lambda[main_ind]+Q_ind, varphi, dvarphi, coils, n_coils, n_segs);
					varphi += dvarphi;
				}
				//deltalambda = addstructs(-1.0, &lambdavar, 1.0, lambda[main_ind]+Q_ind); 
				//printf("varphi = %f\n", varphi);
				//printstruct("fieldline", &fieldline);
				//printstruct("lambdain", &lambdain);
				//printf("\n\n");
				deltalambda = addstructs(1.0, &lambdavar, -1.0, &lambdain); 
				//printstruct("deltalambda", &deltalambda);
				// lambda tangent does not seem to be correct
				//printstruct("lambdaQ", lambda[main_ind]+Q_ind);
				//printstruct("lambdavar", &lambdavar);
				//printstruct("lambdaQ", lambda[main_ind]+Q_ind);
				//printf("\n\n\n\n\n\n\n\n");
				//printf("fieldline->tangent[1][1] = %f\n", fieldline->tangent[1][0]);
				inverseTminusI = invert2x2(deltalambda.tangent, &det_tangent);
				//inverseT = invert2x2(lambda[main_ind][Q_ind].tangent, &det_tangent);
				inverseT = invert2x2(lambdain.tangent, &det_tangent);
				//printf("varphi = %f\n", varphi);
				//printstruct("fieldline", lambda);
				//printf("det(tangent)=%f\n", det_tangent);
		//		printf("tangent[0][0] address is %p\n", &&deltalambda
				//printmat("inverseTminusI", inverseTminusI, 2, 2);
				
				pdeltalambda = calloc(2,sizeof(double*)); 
				*pdeltalambda     = &deltalambda.loc[0];
				*(pdeltalambda+1) = &deltalambda.loc[1];
				matrix = multiply2x2(inverseTminusI, deltalambda.tangent, 2);
				//printmat("matrix", matrix, 2, 2);
				jumptocentre = multiply2x2(inverseTminusI, pdeltalambda, 1);
				jumptocentre[0][0] *= (factor);
				jumptocentre[1][0] *= (factor);
				//printstruct("deltalambda", &deltalambda);
				//printmat("jumptocentre", jumptocentre, 2, 1);
				//printstruct("lambda", lambda[main_ind]+Q_ind);
				//lambda[main_ind][Q_ind].loc[0] = lambdavar.loc[0] - jumptocentre[0][0]; 
				//lambda[main_ind][Q_ind].loc[1] = lambdavar.loc[1] - jumptocentre[1][0];
				//lambda[main_ind][Q_ind].tangent[0][0] = 1.0; lambda[main_ind][Q_ind].tangent[0][1]=0.0; lambda[main_ind][Q_ind].tangent[1][0]=0.0; lambda[main_ind][Q_ind].tangent[1][1]=1.0;
				lambdain.loc[0] = lambdavar.loc[0] - jumptocentre[0][0]; 
				lambdain.loc[1] = lambdavar.loc[1] - jumptocentre[1][0];
				lambdain.tangent[0][0] = 1.0; lambdain.tangent[0][1]=0.0; lambdain.tangent[1][0]=0.0; lambdain.tangent[1][1]=1.0;
				error = sqrt(pow(deltalambda.loc[0], 2.0) + pow(deltalambda.loc[1], 2.0)); 
				//printf("%f\n", deltalambda.loc[0]);
				//printf("%f\n", deltalambda.loc[1]);
				//free(jumptocentre);
				//printf("count=%d\n", count);
				count++;
				lambda[main_ind][Q_ind].loc[0] = lambdain.loc[0];
				lambda[main_ind][Q_ind].loc[1] = lambdain.loc[1];
			} while ( ((fabs(deltalambda.loc[0]/lambdain.loc[0]) > errorlimit) || (fabs(deltalambda.loc[1]/lambdain.loc[1]) > errorlimit)) );
		//clock_t int3 = clock();
		}
	}
	return lambda;
}

//double *adjgrad(double ***coils, int *n_coils, int **n_segs, struct ext_position *ext_centre, struct position *lambda_circ, struct position **lambda_tangent, struct position **mu_tangent, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode) {
//	// declarations
//	//clock_t start = clock();
//	int i=0, q0_ind, Q0_ind, Q_ind, q_ind, main_ind, sec_ind, ind, N_line=0, L_fixedpoints;
//	struct position centre, lambda, mu, sperp;
//	double varphi=0.0, chordlength, chordpluslength, *vecresult, **matrix2;
//	double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
//	double *gradcircumference, *gradtangent, initial=0.0, initialgradcirc=0.0, *adjgradwidth;
//	if (tor_mode % m0_symmetry == 0) 
//		L_fixedpoints = pol_mode;	
//	else 				 
//		L_fixedpoints = m0_symmetry*pol_mode;
//
//	matrix2=set_identity();
//	gradtangent = malloc(L_fixedpoints*sizeof(double));
//	adjgradwidth = malloc(L_fixedpoints*sizeof(double));
//	centre.loc[0] = ext_centre[0].loc[0]; centre.loc[1] = ext_centre[0].loc[1];
//	//lambda_circ.loc[0] = lambda[0]; lambda_circ.loc[1] = lambda[1];
//	centre.tangent = set_identity();
//	//lambda_circ.tangent = set_identity();
//	gradcircumference = &initialgradcirc;
//	N_line = L_fixedpoints*N_gridphi_per_field_period;
//	printf("N_line = %d\n", N_line);
//	printf("field_periods = %d\n", m0_symmetry);
//	printf("field_periods = %d\n", m0_symmetry);
//	printf("field_periods = %d\n", m0_symmetry);
//	varphi = 0.0;
//	for (i=0; i<N_line; i++)
//	{
//		if (i%N_gridphi_per_field_period==0)
//		{
//			main_ind = i/N_gridphi_per_field_period;
//			chordlength= pow( pow(ext_centre[main_ind].chord[0], 2.0) + pow(ext_centre[main_ind].chord[1], 2.0), 0.5);
//			chordpluslength= pow( pow(ext_centre[main_ind].chordplus[0], 2.0) + pow(ext_centre[main_ind].chordplus[1], 2.0), 0.5);
//			lambda_circ->loc[0] += (ext_centre[main_ind].chord[0]/chordlength - ext_centre[main_ind].chordplus[0]/chordpluslength);
//			lambda_circ->loc[1] += (ext_centre[main_ind].chord[1]/chordlength - ext_centre[main_ind].chordplus[1]/chordpluslength);
//		}
//		RK4_adjgrad(gradcircumference, &centre, lambda_circ, varphi, dvarphi, coils, n_coils, n_segs);
//		varphi += dvarphi;
//		//printf("gradcircumference = %f\n", *gradcircumference);
//	}
//	//printf("varphi = %f\n", varphi);
//	printstruct("lambda",lambda_circ);
//	printstruct("lambda",lambda_circ);
//	printstruct("Xp",&centre);
//	printf("gradcircumference = %f\n", *gradcircumference);
//	printf("gradcircumference = %f\n", *gradcircumference);
//
//	centre.loc[0] = ext_centre[0].loc[0]; centre.loc[1] = ext_centre[0].loc[1];
//	centre.tangent[0][0] = centre.tangent[1][1] = 1.0;
//	centre.tangent[0][1] = centre.tangent[1][0] = 0.0;
//	//centre.loc[0] = ext_centre[0].loc[0]; centre.loc[1] = ext_centre[0].loc[1];
//	////lambda_centre.loc[0] = lambda[0]; lambda_centre.loc[1] = lambda[1];
//	//centre.tangent = set_identity();
//	//lambda_centre.tangent = set_identity();
//	q0_ind = ext_centre[0].q0_index; 
//	Q0_ind = q0_ind/L_fixedpoints + 2;
//	sperp.tangent = set_identity();
//	for (main_ind=0; main_ind<L_fixedpoints; main_ind++) {
//		varphi = main_ind*2.0*M_PI/m0_symmetry;
//		initial = 0.0;
//		gradtangent[main_ind] = initial; 
//		centre.tangent[0][0] = centre.tangent[1][1] = 1.0;
//		centre.tangent[0][1] = centre.tangent[1][0] = 0.0;
//		for (Q_ind=0;Q_ind<Q0_ind;Q_ind++) {
//			lambda = lambda_tangent[main_ind][Q_ind];
//			mu = mu_tangent[main_ind][Q_ind];
//			//printstruct("lambda", &lambda);
//			//printstruct("mu", &mu);
//			matrix2[0][0] = matrix2[1][1] = 1.0;
//			matrix2[0][1] = matrix2[1][0] = 0.0;
//			for (ind=0;ind<Q_ind; ind++) {
//				multiply2x2reassign(ext_centre[main_ind].full_tangent, matrix2, 2);
//			}
//			vecresult = multiply(matrix2, ext_centre[main_ind].eperp);
//			sperp.loc[0] = vecresult[0];  sperp.loc[1] = vecresult[1];
//			free(vecresult);
//			sperp.tangent[0][0] = sperp.tangent[1][1] = 1.0; 
//			sperp.tangent[1][0] = sperp.tangent[0][1] = 0.0; 
//			centre.loc[0] = ext_centre[main_ind].loc[0]; centre.loc[1] = ext_centre[main_ind].loc[1];
//			//vecresult = multiply(centre.tangent, ext_centre[main_ind].eperp);
//			//sperp.loc[0] = vecresult[0]; sperp.loc[1] = vecresult[1];
//			//free(vecresult);
//			//printf("sperp=%f\n", sperp.tangent[1][1]);
//			//centre->tangent[0][0] = 1.0 ; centre->tangent[0][1] = 0.0;
//			//centre->tangent[1][0] = 0.0 ; centre->tangent[0][1] = 1.0;
//			for (i=0; i<N_line; i++)
//			{
//				if (i%N_gridphi_per_field_period==0)
//				{
//					sec_ind = i/N_gridphi_per_field_period;
//					//printf("varphi = %f\n", varphi);
//					//printstruct("Xp", &centre);
//					//printf("gradtangent = %f\n", gradtangent[sec_ind]);
//					//if (i/N_gridphi_per_field_period==0) lambda_centre->loc[0] += 1.0;
//					//printstruct("lambda\n", lambda);
//					//if (sec_ind==0) 
//					//{
//					//	printf("main_ind=%d, Q_ind=%d, sec_ind=%d\n", main_ind, Q_ind, sec_ind); 
//					//	printf("gradtangent = %f\n", gradtangent[main_ind]);
//					//	printstruct("mu", &mu);
//					//	printstruct("lambda", &lambda);
//					//}
//					q_ind = Q_ind*L_fixedpoints + sec_ind;
//					if ( q_ind >= q0_ind && q_ind < q0_ind + L_fixedpoints) {
//						mu.loc[0] += (ext_centre[(main_ind+sec_ind)%L_fixedpoints].epar[0]);
//						mu.loc[1] += (ext_centre[(main_ind+sec_ind)%L_fixedpoints].epar[1]);
//					}
//				}
//				//RK4_adjtangent(&centre, &lambda, &sperp, &mu, varphi, dvarphi, coils, n_coils, n_segs);
//				//RK4_adjgrad(gradtangent, &centre, &lambda, varphi, dvarphi, coils, n_coils, n_segs);
//				RK4_adjgradtangent(gradtangent+main_ind, &centre, &lambda, &sperp, &mu, varphi, dvarphi, coils, n_coils, n_segs);
//				varphi += dvarphi;
//			}
//		}
//		//printf("varphi = %f\n", varphi);
//		//printstruct("lambda", &lambda);
//		//printstruct("Xp",&centre);
//	}
//	for (sec_ind=0; sec_ind<L_fixedpoints; sec_ind++) {
//		printf("gradtangent = %f\n", gradtangent[sec_ind]);
//	}
//
//	return gradtangent;
//}

//double *adjgradtangent(double ***coils, int *n_coils, int **n_segs, struct ext_position *ext_centre, struct position **lambdaQ, struct position **muQ, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode) {
//	// declarations
//	//clock_t start = clock();
//	int i=0, centre_ind, q0_ind, Q0_ind, Q_ind, q_ind, main_ind;
//	struct position fieldline, lambda, mu, sperp;
//	double varphi=0.0, *vecresult;
//	double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
//	double *number, initial=0.0;
//	int N_line=0, L_fixedpoints;
//	if (tor_mode % m0_symmetry == 0) 
//		L_fixedpoints = pol_mode;	
//	else 				 
//		L_fixedpoints = m0_symmetry*pol_mode;
//
//	fieldline.loc[0] = ext_centre[0].loc[0]; fieldline.loc[1] = ext_centre[0].loc[1];
//	fieldline.tangent = set_identity();
//	//centre.loc[0] = ext_centre[0].loc[0]; centre.loc[1] = ext_centre[0].loc[1];
//	////lambda_centre.loc[0] = lambda[0]; lambda_centre.loc[1] = lambda[1];
//	//centre.tangent = set_identity();
//	//lambda_centre.tangent = set_identity();
//	q0_ind = ext_centre[0].q0_index; 
//	Q0_ind = q0_ind/L_fixedpoints + 2;
//	number = malloc(L_fixedpoints*sizeof(double));
//	sperp.tangent = set_identity();
//	N_line = L_fixedpoints*N_gridphi_per_field_period;
//	printf("N_line = %d\n", N_line);
//	printf("field_periods = %d\n", m0_symmetry);
//	for (main_ind=0; main_ind<L_fixedpoints; main_ind++) {
//		varphi = main_ind*2.0*M_PI/m0_symmetry;
//		initial = 0.0;
//		number[main_ind] = initial; 
//		fieldline.tangent[0][0] = fieldline.tangent[1][1] = 1.0;
//		fieldline.tangent[0][1] = fieldline.tangent[1][0] = 0.0;
//		for (Q_ind=0;Q_ind<Q0_ind;Q_ind++) {
//			lambda = lambdaQ[main_ind][Q_ind];
//			mu = muQ[main_ind][Q_ind];
//			//printstruct("lambda", &lambda);
//			//printstruct("mu", &mu);
//			fieldline.loc[0] = ext_centre[main_ind].loc[0]; fieldline.loc[1] = ext_centre[main_ind].loc[1];
//			vecresult = multiply(fieldline.tangent, ext_centre[main_ind].eperp);
//			sperp.loc[0] = vecresult[0]; sperp.loc[1] = vecresult[1];
//			free(vecresult);
//			//printf("sperp=%f\n", sperp.tangent[1][1]);
//			//fieldline->tangent[0][0] = 1.0 ; fieldline->tangent[0][1] = 0.0;
//			//fieldline->tangent[1][0] = 0.0 ; fieldline->tangent[0][1] = 1.0;
//			for (i=0; i<N_line; i++)
//			{
//				if (i%N_gridphi_per_field_period==0)
//				{
//					centre_ind = i/N_gridphi_per_field_period;
//					//printf("varphi = %f\n", varphi);
//					//printstruct("Xp", &fieldline);
//					//printf("number = %f\n", number[centre_ind]);
//					//if (i/N_gridphi_per_field_period==0) lambda_centre->loc[0] += 1.0;
//					//printstruct("lambda\n", lambda);
//					//if (centre_ind==0) 
//					//{
//					//	printf("main_ind=%d, Q_ind=%d, centre_ind=%d\n", main_ind, Q_ind, centre_ind); 
//					//	printf("number = %f\n", number[main_ind]);
//					//	printstruct("mu", &mu);
//					//	printstruct("lambda", &lambda);
//					//}
//					q_ind = Q_ind*L_fixedpoints + centre_ind;
//					if ( q_ind >= q0_ind && q_ind < q0_ind + L_fixedpoints) {
//						mu.loc[0] += (ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[0]);
//						mu.loc[1] += (ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[1]);
//					}
//				}
//				//RK4_adjtangent(&fieldline, &lambda, &sperp, &mu, varphi, dvarphi, coils, n_coils, n_segs);
//				//RK4_adjgrad(number, &fieldline, &lambda, varphi, dvarphi, coils, n_coils, n_segs);
//				RK4_adjgradtangent(number+main_ind, &fieldline, &lambda, &sperp, &mu, varphi, dvarphi, coils, n_coils, n_segs);
//				varphi += dvarphi;
//			}
//		}
//		//printf("varphi = %f\n", varphi);
//		//printstruct("lambda", &lambda);
//		//printstruct("Xp",&fieldline);
//	}
//	for (centre_ind=0; centre_ind<L_fixedpoints; centre_ind++) {
//		printf("number = %f\n", number[centre_ind]);
//	}
//	return number;
//}
//can delete this eventually
//struct position *gradcentre(double RR, double ZZ, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode, double ***coils, int *n_coils, int **n_segs) {
//	int i=0, L_fixedpoints;
//	struct position fieldline_start, fieldline, gradfieldline, *gradfieldlineret, deltafieldline;
//	double varphi=0.0, gradX[2], **pdeltafieldline, **inverseT, **inverseTminusI, det_tangent=0.0;
//	double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
//	int N_line=0;
//
//	//evecs1 = malloc(2*sizeof(double*));
//	//evecs1[0] = malloc(2*sizeof(double)); evecs1[1] = malloc(2*sizeof(double));
//	//evals1 = malloc(2*sizeof(double));
//
//	if (tor_mode % m0_symmetry == 0) 
//		L_fixedpoints = pol_mode;	
//	else 				   
//		L_fixedpoints = m0_symmetry*pol_mode;
//
//	N_line = L_fixedpoints*N_gridphi_per_field_period;
//	fieldline_start.loc[0] = RR; fieldline_start.loc[1] = ZZ;
//	fieldline_start.tangent = set_identity();
//	printf("fieldline_start = (%f, %f)\n", RR, ZZ);
//	printf("RR = %f, ZZ = %f\n", RR, ZZ);
//	fieldline = fieldline_start;
//	gradfieldline = fieldline_start;
//	gradfieldline.loc[0] = 0.0; gradfieldline.loc[1] = 0.0;
//	gradfieldline.tangent = set_zeros();
//	varphi = 0.0;
//	printf("varphi = %f\n", varphi);
//	//for (i=0; i<8*N_line + N_line/2; i++) 
//	for (i=0; i<8*N_line; i++) {
//		RK4_wgrad(&fieldline, &gradfieldline, varphi, dvarphi, coils, n_coils, n_segs);
//		varphi += dvarphi;
//		printstruct("gradfieldline\n", &gradfieldline);
//		//if (i%N_gridphi_per_field_period==0) {
//		//	centre_ind = (i / N_gridphi_per_field_period);
//		//	centre[centre_ind % L_fixedpoints] = fieldline;
//		//	printf("varphi = %f\n", varphi);
//		//	printstruct("fieldline[i]\n", &fieldline);
//		//	inverted = invert2x2(centre[(centre_ind-1) % L_fixedpoints].tangent, &detcentre);
//		//	printmat("inverted", inverted, 2, 2);
//		//	//ptarray[0] = &centre[centre_ind % L_fixedpoints].tangent[0][0];
//		//	//ptarray[1] = &centre[centre_ind % L_fixedpoints].tangent[1][0];
//		//}
//	}
//	//gradfieldline.loc[0]  = 0.00133;
//	//gradfieldline.loc[1]  = -0.09218;
//	gradfieldlineret = &gradfieldline;
//	printf("(%f, %f)\n", gradfieldlineret->loc[0], gradfieldlineret->loc[1]);
//	deltafieldline = addstructs(-1.0, &fieldline, 1.0, &fieldline_start); 
//	//printmat("deltafieldline", deltafieldline.tangent, 2, 2);
//	//printf("fieldline->tangent[1][1] = %f\n", fieldline->tangent[1][0]);
//	inverseT = invert2x2(fieldline.tangent, &det_tangent);
//	printf("det = %f\n", det_tangent);
//	inverseTminusI = invert2x2(deltafieldline.tangent, &det_tangent);
//	printf("det = %f\n", det_tangent);
//	
//	pdeltafieldline = calloc(2,sizeof(double*)); 
//	*pdeltafieldline      = &gradfieldline.loc[0];
//	*(pdeltafieldline+1)  = &gradfieldline.loc[1];
//
//	//gradX = multiply2x2(inverseTminusI, pdeltafieldline, 1);
//	gradX[0] = inverseTminusI[0][0]*gradfieldline.loc[0] + inverseTminusI[0][1]*gradfieldline.loc[1]; 
//	gradX[1] = inverseTminusI[1][0]*gradfieldline.loc[0] + inverseTminusI[1][1]*gradfieldline.loc[1]; 
//	printf("%f \n", inverseTminusI[0][0]*gradfieldline.loc[0]); 
//	printf("%f \n", inverseTminusI[0][1]*gradfieldline.loc[1]); 
//	printf("%f \n", inverseTminusI[1][0]*gradfieldline.loc[0]); 
//	printf("%f \n", inverseTminusI[1][1]*gradfieldline.loc[1]); 
//	printf("%f \n", inverseTminusI[0][0]); 
//	printf("%f \n", inverseTminusI[0][1]); 
//	printf("%f \n", inverseTminusI[1][0]); 
//	printf("%f \n", inverseTminusI[1][1]); 
//
//	printf("(%f, %f)\n", gradX[0], gradX[1]);
//	gradfieldline.loc[0] = gradX[0];
//	gradfieldline.loc[1] = gradX[1];
//	gradfieldlineret = &gradfieldline;
//	printf("(%f, %f)\n", gradfieldlineret->loc[0], gradfieldlineret->loc[1]);
//	return gradfieldlineret;
//}

//struct ext_position *gradalongcentre(double RR, double ZZ, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode, double ***coils, int *n_coils, int **n_segs) {
//	int i=0, q0, L_fixedpoints;
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
//		L_fixedpoints = pol_mode;	
//	else 				   
//		L_fixedpoints = m0_symmetry*pol_mode;
//	N_line = L_fixedpoints*N_gridphi_per_field_period;
//	centre = malloc(L_fixedpoints*sizeof(struct position));
//	gradcentre = malloc(L_fixedpoints*sizeof(struct position));
//	ext_centre = malloc(L_fixedpoints*sizeof(struct ext_position));
//	grad_ext_centre = malloc(L_fixedpoints*sizeof(struct ext_position));
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
//			centre[centre_ind % L_fixedpoints] = fieldline;
//			gradcentre[centre_ind % L_fixedpoints] = gradfieldline;
//			printf("varphi = %f\n", varphi);
//			printstruct("fieldline[i]\n", &fieldline);
//			inverted = invert2x2(centre[(centre_ind-1) % L_fixedpoints].tangent, &detcentre);
//			printmat("inverted", inverted, 2, 2);
//			//ptarray[0] = &centre[centre_ind % L_fixedpoints].tangent[0][0];
//			//ptarray[1] = &centre[centre_ind % L_fixedpoints].tangent[1][0];
//			ext_centre[centre_ind % L_fixedpoints].loc[0] = fieldline.loc[0];
//			ext_centre[centre_ind % L_fixedpoints].loc[1] = fieldline.loc[1];
//			ext_centre[centre_ind % L_fixedpoints].part_tangent = multiply2x2(centre[centre_ind % L_fixedpoints].tangent, inverted, 2);
//			//printmat("ext_centre.part_tangent", ext_centre[centre_ind % L_fixedpoints].part_tangent, 2, 2);
//			grad_ext_centre[centre_ind % L_fixedpoints].loc[0] = gradfieldline.loc[0];
//			grad_ext_centre[centre_ind % L_fixedpoints].loc[1] = gradfieldline.loc[1];
//			grad_ext_centre[centre_ind % L_fixedpoints].part_tangent = add2x2( 
//			1.0 , multiply2x2(gradcentre[centre_ind % L_fixedpoints].tangent, inverted, 2), 
//			-1.0, multiply2x2(ext_centre[centre_ind % L_fixedpoints].part_tangent, multiply2x2(gradcentre[(centre_ind-1) % L_fixedpoints].tangent, inverted, 2), 2), 2 );
//
//			//printmat("ext_centre.part_tangent", ext_centre[centre_ind % L_fixedpoints].part_tangent, 2, 2);
//		}
//	}
//	// The part below calculates the circumference and the full orbit tangent maps
//	circumference = 0.0;
//	for (centre_ind=0;centre_ind<L_fixedpoints;centre_ind++) {
//		circumference += sqrt(pow(ext_centre[(centre_ind+1)%L_fixedpoints].loc[0] 
//					        - ext_centre[centre_ind].loc[0], 2.0)
//				    	    + pow(ext_centre[(centre_ind+1)%L_fixedpoints].loc[1] 
//						- ext_centre[centre_ind].loc[1], 2.0));
//		ext_centre[centre_ind].epar = malloc(2*sizeof(double)); ext_centre[centre_ind].eperp = malloc(2*sizeof(double));
//		ext_centre[centre_ind].full_tangent = multiply2x2(ext_centre[(centre_ind +2)%L_fixedpoints].part_tangent, ext_centre[(centre_ind+1)%L_fixedpoints].part_tangent, 2);
//		grad_ext_centre[centre_ind].full_tangent = add2x2( 1.0, multiply2x2(grad_ext_centre[(centre_ind +2)%L_fixedpoints].part_tangent, ext_centre[(centre_ind+1)%L_fixedpoints].part_tangent, 2), 1.0, multiply2x2(ext_centre[(centre_ind +2)%L_fixedpoints].part_tangent, grad_ext_centre[(centre_ind+1)%L_fixedpoints].part_tangent, 2), 2);
//
//		for (sec_ind=3;sec_ind<L_fixedpoints+1;sec_ind++) {
//			multiply2x2reassign(ext_centre[(centre_ind + sec_ind)%L_fixedpoints].part_tangent, ext_centre[centre_ind].full_tangent, 2);
//			spare_matrix = multiply2x2(grad_ext_centre[(centre_ind + sec_ind)%L_fixedpoints].part_tangent, ext_centre[centre_ind].full_tangent, 2);
//			multiply2x2reassign(ext_centre[(centre_ind + sec_ind)%L_fixedpoints].part_tangent, grad_ext_centre[centre_ind].full_tangent, 2);
//			add2x2reassign(1.0, spare_matrix, 1.0, grad_ext_centre[centre_ind].full_tangent, 2);
//			free(spare_matrix[0]); free(spare_matrix[1]);
//			free(spare_matrix);
//		}
//		printmat("ext_centre.fulltangent", ext_centre[centre_ind].full_tangent, 2, 2);
//		symmeigs(ext_centre[centre_ind].full_tangent, ext_centre[centre_ind].eperp, ext_centre[centre_ind].epar, evals);
//		printf("eigenvectors are (%f, %f) and (%f, %f)\n", ext_centre[centre_ind].eperp[0], ext_centre[centre_ind].eperp[1], ext_centre[centre_ind].epar[0], ext_centre[centre_ind].epar[1]);
//	}
//	for (main_ind=0;main_ind<L_fixedpoints; main_ind++) {
//		deltadet = ext_centre[main_ind].full_tangent[0][0]*grad_ext_centre[main_ind].full_tangent[1][1] + grad_ext_centre[main_ind].full_tangent[1][1]*ext_centre[main_ind].full_tangent[1][1] - ext_centre[main_ind].full_tangent[0][1]*grad_ext_centre[main_ind].full_tangent[1][0] - grad_ext_centre[main_ind].full_tangent[0][1]*ext_centre[main_ind].full_tangent[1][0]; 
//		printf("deltadet=%f\n\n", deltadet);
//		linalg2x2(ext_centre[main_ind].full_tangent, evecs1, evals1, &det, &trace);
//		ext_centre[main_ind].circumference = circumference;
//		if ((fabs(evecs1[0][0]) < small) || (fabs(evecs1[0][1]) < small)) {
//			ext_centre[main_ind].angle = evals1[1];
//			omega = m0_symmetry*evals1[1]/(2.0*M_PI*L_fixedpoints); // I have changed it compared to Cary and Hanson's paper
//			q0 = (int) (m0_symmetry/(4.0*omega) - L_fixedpoints/2.0 + 0.5);
//			ext_centre[main_ind].q0_index = q0;
//		}
//		printf("angle=%f\n", ext_centre[main_ind].angle);
//		printf("q0=%d\n", ext_centre[main_ind].q0_index);
//		ext_centre[main_ind].long_tangent = malloc(L_fixedpoints*sizeof(double));
//		grad_ext_centre[main_ind].long_tangent = malloc(L_fixedpoints*sizeof(double));
//		for (centre_ind=0;centre_ind<L_fixedpoints;centre_ind++) {
//			//ext_centre[centre_ind].long_tangent = set_identity();
//			//for (sec_ind=0; sec_ind<q0/L_fixedpoints; sec_ind++) {
//			//	multiply2x2reassign(ext_centre[0].full_tangent, ext_centre[centre_ind].long_tangent, 2); 
//			//}
//			//for (sec_ind=0; sec_ind<q0%L_fixedpoints; sec_ind++) {
//			//	multiply2x2reassign(ext_centre[sec_ind+1].part_tangent, ext_centre[centre_ind].long_tangent, 2); 
//			//}
//			ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints] = set_identity();
//			grad_ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints] = set_zeros();
//			for (sec_ind=0; sec_ind<(q0+centre_ind)/L_fixedpoints; sec_ind++) {
//				multiply2x2reassign(ext_centre[main_ind].full_tangent, ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints], 2); 
//				spare_matrix = multiply2x2(grad_ext_centre[main_ind].full_tangent, ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints], 2);
//				multiply2x2reassign(ext_centre[main_ind].full_tangent, grad_ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints], 2);
//				add2x2reassign(1.0, spare_matrix, 1.0, grad_ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints], 2);
//				free(spare_matrix[0]); free(spare_matrix[1]);
//				free(spare_matrix);
//			}
//			for (sec_ind=0; sec_ind<(q0+centre_ind)%L_fixedpoints; sec_ind++) {
//				multiply2x2reassign(ext_centre[(sec_ind+1+main_ind)%L_fixedpoints].part_tangent, ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints], 2); 
//				spare_matrix = multiply2x2(grad_ext_centre[(sec_ind+1+main_ind)%L_fixedpoints].part_tangent, ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints], 2);
//				multiply2x2reassign(ext_centre[(sec_ind+1+main_ind)%L_fixedpoints].part_tangent, grad_ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints], 2);
//				add2x2reassign(1.0, spare_matrix, 1.0, grad_ext_centre[main_ind].long_tangent[(main_ind+q0+centre_ind)%L_fixedpoints], 2);
//				free(spare_matrix[0]); free(spare_matrix[1]);
//				free(spare_matrix);
//			}
//		}
//	}
//	//ext_centre[centre_index % L_fixedpoints].part_tangent = multiply2x2(centre[centre_index % L_fixedpoints].tangent, inverted, 2);
//	//clock_t int3 = clock();
//	return grad_ext_centre;
//} 

// can delete this eventually
//struct ext_position *gradalongcentrealt(double RR, double ZZ, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode, double ***coils, int *n_coils, int **n_segs) {
//	int i=0, q0, L_fixedpoints;
//	struct position fieldline_start, fieldline, gradfieldline, *centre, *gradcentre;
//	double evals[2], **sigma; //**gradsymm, **gradsigmaM, 
//	double varphi=0.0, omega, actualdet, graddet;
//	double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
//	double **inverted, detcentre, **evecs1, *evals1, det=0, trace=0; // **spare_matrix;
//	double circumference;
//	int N_line=0, main_ind, centre_ind;
//	struct ext_position *ext_centre, *grad_ext_centre;
//
//	evecs1 = malloc(2*sizeof(double*));
//	evecs1[0] = malloc(2*sizeof(double)); evecs1[1] = malloc(2*sizeof(double));
//	evals1 = malloc(2*sizeof(double));
//
//	if (tor_mode % m0_symmetry == 0) 
//		L_fixedpoints = pol_mode;	
//	else 				   
//		L_fixedpoints = m0_symmetry*pol_mode;
//	N_line = L_fixedpoints*N_gridphi_per_field_period;
//	centre = malloc(L_fixedpoints*sizeof(struct position));
//	gradcentre = malloc(L_fixedpoints*sizeof(struct position));
//	ext_centre = malloc(L_fixedpoints*sizeof(struct ext_position));
//	grad_ext_centre = malloc(L_fixedpoints*sizeof(struct ext_position));
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
//			centre[centre_ind % L_fixedpoints] = fieldline;
//			gradcentre[centre_ind % L_fixedpoints] = gradfieldline;
//			printf("varphi = %f\n", varphi);
//			printstruct("fieldline[i]\n", &fieldline);
//			inverted = invert2x2(centre[(centre_ind-1) % L_fixedpoints].tangent, &detcentre);
//			printmat("inverted", inverted, 2, 2);
//			//ptarray[0] = &centre[centre_ind % L_fixedpoints].tangent[0][0];
//			//ptarray[1] = &centre[centre_ind % L_fixedpoints].tangent[1][0];
//			ext_centre[centre_ind % L_fixedpoints].loc[0] = fieldline.loc[0];
//			ext_centre[centre_ind % L_fixedpoints].loc[1] = fieldline.loc[1];
//			ext_centre[centre_ind % L_fixedpoints].part_tangent = multiply2x2(centre[centre_ind % L_fixedpoints].tangent, inverted, 2);
//			//printmat("ext_centre.part_tangent", ext_centre[centre_ind % L_fixedpoints].part_tangent, 2, 2);
//			//printmat("ext_centre.part_tangent", ext_centre[centre_ind % L_fixedpoints].part_tangent, 2, 2);
//		}
//	}
//	// The part below calculates the circumference and the full orbit tangent maps
//	circumference = 0.0;
//	sigma = malloc(2*sizeof(double)); sigma[0] = malloc(2*sizeof(double)); sigma[1] = malloc(2*sizeof(double));
//	sigma[0][0] = 0.0;
//	sigma[0][1] = -1.0;
//	sigma[1][0] = 1.0;
//	sigma[1][1] = 0.0;
//	//gradsymm = malloc(2*sizeof(double)); gradsymm[0] = malloc(2*sizeof(double)); gradsymm[1] = malloc(2*sizeof(double));
//	for (main_ind=0;main_ind<L_fixedpoints;main_ind++) {
//		circumference += sqrt(pow(ext_centre[(main_ind+1)%L_fixedpoints].loc[0] 
//					        - ext_centre[main_ind].loc[0], 2.0)
//				    	    + pow(ext_centre[(main_ind+1)%L_fixedpoints].loc[1] 
//						- ext_centre[main_ind].loc[1], 2.0));
//		ext_centre[main_ind].epar = malloc(2*sizeof(double)); ext_centre[main_ind].eperp = malloc(2*sizeof(double));
//		//grad_ext_centre[main_ind].epar = malloc(2*sizeof(double)); grad_ext_centre[main_ind].eperp = malloc(2*sizeof(double));
//		ext_centre[main_ind].full_tangent = multiply2x2(ext_centre[(main_ind +2)%L_fixedpoints].part_tangent, ext_centre[(main_ind+1)%L_fixedpoints].part_tangent, 2);
//
//		fieldline = centre[main_ind];
//		fieldline.tangent[0][0] = 1.0;
//		fieldline.tangent[0][1] = 0.0;
//		fieldline.tangent[1][0] = 0.0;
//		fieldline.tangent[1][1] = 1.0;
//		gradfieldline = gradcentre[main_ind];
//		gradfieldline.loc[0] = 0.0; gradfieldline.loc[1] = 0.0;
//		gradfieldline.tangent[0][0] = 0.0;
//		gradfieldline.tangent[0][1] = 0.0;
//		gradfieldline.tangent[1][0] = 0.0;
//		gradfieldline.tangent[1][1] = 0.0;
//		varphi = main_ind*2.0*M_PI/m0_symmetry;
//		for (i=0; i<N_gridphi_per_field_period*L_fixedpoints; i++) {
//			centre_ind = (i / N_gridphi_per_field_period);
//			RK4_wgrad(&fieldline, &gradfieldline, varphi, dvarphi, coils, n_coils, n_segs);
//			varphi += dvarphi;
//		}
//		ext_centre[main_ind].full_tangent = fieldline.tangent;
//		grad_ext_centre[main_ind].full_tangent = gradfieldline.tangent;
//		printmat("ext_centre.fulltangent", ext_centre[main_ind].full_tangent, 2, 2);
//		printmat("grad_ext_centre.fulltangent1", grad_ext_centre[main_ind].full_tangent, 2, 2);
//		symmeigs(ext_centre[main_ind].full_tangent, ext_centre[main_ind].eperp, ext_centre[main_ind].epar, evals);
//		printmat("grad_ext_centre.fulltangent2", grad_ext_centre[main_ind].full_tangent, 2, 2);
//		graddet = ext_centre[main_ind].full_tangent[0][0]*grad_ext_centre[main_ind].full_tangent[1][1] + grad_ext_centre[main_ind].full_tangent[1][1]*ext_centre[main_ind].full_tangent[1][1] - ext_centre[main_ind].full_tangent[0][1]*grad_ext_centre[main_ind].full_tangent[1][0] - grad_ext_centre[main_ind].full_tangent[0][1]*ext_centre[main_ind].full_tangent[1][0]; 
//		printf("graddet=%f\n\n", graddet);
//		actualdet = ext_centre[main_ind].full_tangent[0][0]*ext_centre[main_ind].full_tangent[1][1] -  ext_centre[main_ind].full_tangent[0][1]*ext_centre[main_ind].full_tangent[1][0]; 
//		printf("det=%f\n\n", actualdet);
//		//gradsigmaM = multiply2x2(sigma, grad_ext_centre[main_ind].full_tangent, 2);
//		//gradsymm[0][0] = gradsigmaM[0][0];
//		//gradsymm[1][1] = gradsigmaM[1][1];
//		//gradsymm[0][1] = 0.5*gradsigmaM[1][0] + 0.5*gradsigmaM[0][1];
//		//gradsymm[1][0] = gradsymm[0][1];
//		
//		//deltaepar = inner(ext_centre[main_ind].eperp,  gradsymm, ext_centre[main_ind].epar)/(evals[0]-evals[1]);
//		//deltaeperp = - inner(ext_centre[main_ind].epar,  gradsymm, ext_centre[main_ind].eperp)/(evals[0]-evals[1]);
//		//grad_ext_centre[main_ind].epar[0] = deltaepar*ext_centre[main_ind].eperp[0]; grad_ext_centre[main_ind].epar[1] = deltaepar*ext_centre[main_ind].eperp[1];
//		//grad_ext_centre[main_ind].eperp[0] = deltaeperp*ext_centre[main_ind].epar[0]; grad_ext_centre[main_ind].eperp[1] = deltaeperp*ext_centre[main_ind].epar[1];
//		//free(gradsigmaM);
//		//grad_ext_centre[main_ind].eperp[0] 
//		printmat("grad_ext_centre.fulltangent3", grad_ext_centre[main_ind].full_tangent, 2, 2);
//	}
//	for (main_ind=0;main_ind<L_fixedpoints; main_ind++) {
//		printmat("grad_ext_centre.fulltangent4", grad_ext_centre[main_ind].full_tangent, 2, 2);
//		linalg2x2(ext_centre[main_ind].full_tangent, evecs1, evals1, &det, &trace);
//		ext_centre[main_ind].circumference = circumference;
//		if ((fabs(evecs1[0][0]) < small) || (fabs(evecs1[0][1]) < small)) {
//			ext_centre[main_ind].angle = evals1[1];
//			omega = m0_symmetry*evals1[1]/(2.0*M_PI*L_fixedpoints); // I have changed it compared to Cary and Hanson's paper
//			q0 = (int) (m0_symmetry/(4.0*omega) - L_fixedpoints/2.0 + 0.5);
//			ext_centre[main_ind].q0_index = q0;
//			grad_ext_centre[main_ind].angle = -0.5*(grad_ext_centre[main_ind].full_tangent[0][0]+grad_ext_centre[main_ind].full_tangent[1][1])/sin(ext_centre[main_ind].angle);
//		}
//		printmat("grad_ext_centre.fulltangent", grad_ext_centre[main_ind].full_tangent, 2, 2);
//		//printf("gradangle pieces are %f %f\n", grad_ext_centre[main_ind].full_tangent[0][0]+grad_ext_centre[main_ind].full_tangent[1][1]);
//		printf("gradangle=%f\n", grad_ext_centre[main_ind].angle);
//		printf("angle=%f\n", ext_centre[main_ind].angle);
//		printf("q0=%d\n", ext_centre[main_ind].q0_index);
//		varphi = 0.0; //modify this 
//		ext_centre[main_ind].long_tangent = malloc(L_fixedpoints*sizeof(double));
//		grad_ext_centre[main_ind].long_tangent = malloc(L_fixedpoints*sizeof(double));
//		fieldline = centre[main_ind];
//		fieldline.tangent[0][0] = 1.0;
//		fieldline.tangent[0][1] = 0.0;
//		fieldline.tangent[1][0] = 0.0;
//		fieldline.tangent[1][1] = 1.0;
//		gradfieldline.loc[0] = 0.0; gradfieldline.loc[1] = 0.0;
//		gradfieldline.tangent[0][0] = 0.0;
//		gradfieldline.tangent[0][1] = 0.0;
//		gradfieldline.tangent[1][0] = 0.0;
//		gradfieldline.tangent[1][1] = 0.0;
//		for (i=0; i<N_gridphi_per_field_period*(q0+L_fixedpoints); i++) {
//			centre_ind = (i / N_gridphi_per_field_period);
//			if ( (i%N_gridphi_per_field_period==0) && (centre_ind/q0 == 1) ) {
//				ext_centre[main_ind].long_tangent[(main_ind+centre_ind)%L_fixedpoints] = fieldline.tangent;
//				grad_ext_centre[main_ind].long_tangent[(main_ind+centre_ind)%L_fixedpoints] = gradfieldline.tangent;
//				printf("index = %d\n", centre_ind%q0);
//			}
//			RK4_wgrad(&fieldline, &gradfieldline, varphi, dvarphi, coils, n_coils, n_segs);
//			varphi += dvarphi;
//		}
//	}
//	//ext_centre[centre_index % L_fixedpoints].part_tangent = multiply2x2(centre[centre_index % L_fixedpoints].tangent, inverted, 2);
//	//clock_t int3 = clock();
//	return grad_ext_centre;
//} 

//can soon delete this
//double *gradislandwidth(struct ext_position *ext_fieldline, struct ext_position *grad_ext_fieldline, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode) {
//	int main_index, centre_index, L_fixedpoints;
//	double *wperp, *gradwperp, sum_matrix_elements, grad_sum_matrix_elements, matrix_element, grad_matrix_element;
//	double circumference=ext_fieldline[0].circumference;
//	if (tor_mode % m0_symmetry == 0) L_fixedpoints = pol_mode;	
//	else 				   L_fixedpoints = m0_symmetry*pol_mode;
//	wperp = malloc(L_fixedpoints*sizeof(double));
//	gradwperp = malloc(L_fixedpoints*sizeof(double));
//	for (main_index=0;main_index<L_fixedpoints; main_index++) {
//		//circumference = 0.0;		
//		sum_matrix_elements = 0.0;		
//		grad_sum_matrix_elements = 0.0;		
//		for (centre_index=0;centre_index<L_fixedpoints; centre_index++) {
//			matrix_element = inner(ext_fieldline[centre_index].epar, 
//					       ext_fieldline[main_index].long_tangent[centre_index], 
//					       ext_fieldline[main_index].eperp);	
//			grad_matrix_element = ( inner(ext_fieldline[centre_index].epar, grad_ext_fieldline[main_index].long_tangent[centre_index], ext_fieldline[main_index].eperp)/sin( ext_fieldline[main_index].angle*(ext_fieldline[main_index].q0_index+centre_index)/L_fixedpoints) 
//			- cos( ext_fieldline[main_index].angle*(ext_fieldline[main_index].q0_index+centre_index)/L_fixedpoints ) * ( grad_ext_fieldline[main_index].angle*(ext_fieldline[main_index].q0_index+centre_index)/L_fixedpoints) * inner( ext_fieldline[centre_index].epar, ext_fieldline[main_index].long_tangent[centre_index], ext_fieldline[main_index].eperp ) / pow(sin( ext_fieldline[main_index].angle*(ext_fieldline[main_index].q0_index+centre_index)/L_fixedpoints), 2.0 ) ) * matrix_element/fabs(matrix_element);
//			printf("%f %f\n", cos( ext_fieldline[main_index].angle*(ext_fieldline[main_index].q0_index+centre_index)/L_fixedpoints) ,  (grad_ext_fieldline[main_index].angle*(ext_fieldline[main_index].q0_index+centre_index)/L_fixedpoints)) ;
//			printf("%f\n", grad_ext_fieldline[main_index].angle);
//				//		+inner(ext_fieldline[centre_index].epar, 
//				//	       ext_fieldline[main_index].long_tangent[centre_index], 
//				//	       grad_ext_fieldline[main_index].eperp)
//				//		+inner(grad_ext_fieldline[centre_index].epar, 
//				//	       ext_fieldline[main_index].long_tangent[centre_index], 
//				//	       ext_fieldline[main_index].eperp))*matrix_element/fabs(matrix_element);	
//			sum_matrix_elements += fabs(matrix_element);
//			grad_sum_matrix_elements += grad_matrix_element;
//			//printmat("symmetrized invariant matrix", ext_fieldline[main_index].long_tangent[centre_index], 2, 2); 
//			printmat("symmetrized invariant matrix", grad_ext_fieldline[main_index].long_tangent[centre_index], 2, 2); 
//			printf("%f\n", inner(ext_fieldline[centre_index].epar, 
//					       grad_ext_fieldline[main_index].long_tangent[centre_index], 
//					       ext_fieldline[main_index].eperp));	
//			//printf("%f\n", inner(grad_ext_fieldline[centre_index].epar, 
//			//		       ext_fieldline[main_index].long_tangent[centre_index], 
//			//		       ext_fieldline[main_index].eperp));	
//			//printf("%f\n", inner(ext_fieldline[centre_index].epar, 
//			//		       ext_fieldline[main_index].long_tangent[centre_index], 
//			//		       grad_ext_fieldline[main_index].eperp));	
//			//sum_matrix_elements += fabs(matrix_element)/sin(ext_fieldline[main_index].angle*ext_fieldline[main_index].q0_index/L_fixedpoints); //This has the sin piece in the denominator, which can be safely set to unity
//			//printf("cosine = %f\n", cos(ext_fieldline[main_index].angle*ext_fieldline[main_index].q0_index/L_fixedpoints));
//			//circumference += sqrt(pow(ext_fieldline[(centre_index+1)%L_fixedpoints].loc[0] 
//			//		        - ext_fieldline[centre_index].loc[0], 2.0)
//			//	    	    + pow(ext_fieldline[(centre_index+1)%L_fixedpoints].loc[1] 
//			//			- ext_fieldline[centre_index].loc[1], 2.0));
//			//printf("circumference = %f\nmatrix_element = %f\n", circumference, matrix_element);
//			//printf("other circumference = %f\n", ext_fieldline[centre_index].circumference);
//		}
//		//circumference = 2.0*M_PI*0.2097;
//		printf("circumference = %f\n", circumference);
//		wperp[main_index] = 2.0*L_fixedpoints*circumference/(M_PI*pol_mode*sum_matrix_elements);
//		gradwperp[main_index] = - ( grad_sum_matrix_elements / sum_matrix_elements ) * wperp[main_index];
//		printf("width = %.10f for index= %d\n", wperp[main_index], main_index);
//		printf("gradwidth = %f for index= %d\n", gradwperp[main_index], main_index);
//	}
//	return gradwperp;
//}

