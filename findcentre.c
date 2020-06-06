/* Author: Alessandro Geraldini */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_fit.h>
#include "isc.h"

int rho(int fieldaligned_index, int n_turns, int L_fixedpoints) {
	int reordered_index;
	reordered_index = (n_turns*fieldaligned_index)%L_fixedpoints;
	return reordered_index;
	}

int inverserho(int reordered_index, int n_turns, int L_fixedpoints) {
	int fieldaligned_index;
	fieldaligned_index = ( ( (reordered_index%n_turns)*L_fixedpoints + reordered_index ) / n_turns ) % L_fixedpoints;
	return fieldaligned_index;
}

void printstructposition(char *name, struct position *input) {
	printf("structure position %s:\nposition = (%10.8f, %10.8f)\ntangent = |%10.8f %10.8f|\n          |%10.8f %10.8f|\n", 
	 name, input->loc[0], input->loc[1], input->tangent[0][0], input->tangent[0][1], input->tangent[1][0], input->tangent[1][1]);
	//printf("for structure position %s:\nposition = (%f, %f)\n", name, input->loc[0], input->loc[1]);
	//printmat("tangent", input->tangent, 2, 2);
}

struct position *solve_magneticaxis(double ***coils, int n_coils, int *n_segs, struct field **Bfieldsaved, struct position *fieldline, int N_gridphi_toroidal) {
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
			//	printstructposition("fieldline[i]\n", fieldline);
			//}
			//printstructposition("fieldline[i]\n", fieldline);
			RK4step(fieldline, varphi, dvarphi, coils, n_coils, n_segs, Bfieldsaved[i]);
			//printf("B field saved value (BR) = %f\n", Bfieldsaved[i][0].value[0]);
			varphi += dvarphi;
		}
		deltafieldline = addstructs(1.0, fieldline, -1.0, &fieldline_start); 
		//printf("fieldline->tangent[1][1] = %f\n", fieldline->tangent[1][0]);
		inverseTminusI = invert2x2(deltafieldline.tangent, &det_tangent);
		inverseT = invert2x2(fieldline->tangent, &det_tangent);
		//printstructposition("fieldline", fieldline); //printf("det(tangent)=%f\n", det_tangent);
		pdeltafieldline = calloc(2,sizeof(double*)); 
		*pdeltafieldline      = &deltafieldline.loc[0];
		*(pdeltafieldline+1)  = &deltafieldline.loc[1];
		jumptocentre          = multiply2x2(inverseTminusI, pdeltafieldline, 1);
		jumptocentre[0][0]   *= factor;
		jumptocentre[1][0]   *= factor;
		//printstructposition("deltafieldline", &deltafieldline);
		//printmat("jumptocentre", jumptocentre, 2, 1);
		//printstructposition("fieldline_start", &fieldline_start);
		//printstructposition("fieldline", fieldline);
		fieldline->loc[0] = fieldline_start.loc[0] - jumptocentre[0][0]; 
		fieldline->loc[1] = fieldline_start.loc[1] - jumptocentre[1][0];
		error = sqrt(pow(deltafieldline.loc[0], 2.0) + pow(deltafieldline.loc[1], 2.0)); 
		//free(jumptocentre[0]);
		//free(jumptocentre[1]);
		//free(jumptocentre);
		i = 1;
	} while(error>errorlimit);
	//clock_t int3 = clock();
	return centre;
}

struct position *solve_islandcenter(double ***coils, int n_coils, int *n_segs, struct field **Bfieldsaved, struct position *fieldline, int m0_fieldperiods, int L_fixedpoints, int N_gridphi_fieldperiod) {
	// declarations
	//clock_t start = clock();
	int i=0, N_line;
	struct position fieldline_start, deltafieldline;
	double varphi=0.0, **jumptocentre=calloc(1,sizeof(double*)); //, **testunity;
	double dvarphi = 2.0*M_PI/(m0_fieldperiods*N_gridphi_fieldperiod);
	double **matrix;
	double **pdeltafieldline, det_tangent;
	double **inverseTminusI, **inverseT;
	double error = 1.0, errorlimit = 0.000000000001;
	double factor = 1.0;
	struct position *centre;
	N_line = L_fixedpoints*N_gridphi_fieldperiod;
	centre = malloc(N_line*sizeof(struct position));
	printf("N_line = %d\n", N_line);
	do {
		fieldline_start = *fieldline;
		varphi = 0.0;
		//printf("RR=%f\n", fieldline->loc[0]);
		for (i=0; i<N_line; i++)
		{
			centre[i] = *fieldline;
			//printstructposition("fieldline[i]\n", fieldline);
			RK4step(fieldline, varphi, dvarphi, coils, n_coils, n_segs, Bfieldsaved[i]);
			//printf("B field saved value (BR) = %f\n", Bfieldsaved[i][0].value[0]);
			varphi += dvarphi;
			//printf("where do I seg fault?\n");
		}
		deltafieldline = addstructs(1.0, fieldline, -1.0, &fieldline_start); 
		//printf("fieldline->tangent[1][1] = %f\n", fieldline->tangent[1][0]);
		inverseTminusI = invert2x2(deltafieldline.tangent, &det_tangent);
		inverseT = invert2x2(fieldline->tangent, &det_tangent);
		//printf("varphi = %f\n", varphi);
		//printstructposition("fieldline", fieldline);
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
		//printstructposition("deltafieldline", &deltafieldline);
		//printmat("jumptocentre", jumptocentre, 2, 1);
		//printstructposition("fieldline_start", &fieldline_start);
		//printstructposition("fieldline", fieldline);
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

struct ext_position *solve_islandcenter_full(double *islandcenter, double *axis, int *n_turns, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod, struct field **Bfield_saved) {
	int i=0, q0, clockwise=2, kminus, kplus, rhominus, rhoplus, number_turns;
	struct position fieldline_start, fieldline, adjfieldline, *centre, *adjcentre;
	double varphi=0.0, omega, rawangle, rawangle_principal, *angle_axis, *phidata=malloc(L_fixedpoints*sizeof(double));
	double c0, c1, cov00, cov01, cov11, sumsq;
	double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*m0_symmetry);
	double **inverted, **adjinverted, detcentre, adjdetcentre, **evecs1, *evals1, det=0, trace=0;
	double circumference, ref_angle, **result, **result2;
	double *evals=malloc(2*sizeof(double));
	int N_line=0, main_ind, centre_ind, sec_ind;
	struct ext_position *ext_centre;

	evecs1 = malloc(2*sizeof(double*));
	evecs1[0] = malloc(2*sizeof(double)); evecs1[1] = malloc(2*sizeof(double));
	evals1 = malloc(2*sizeof(double));
	N_line = L_fixedpoints*N_gridphi_fieldperiod;
	angle_axis = malloc(L_fixedpoints*sizeof(double));
	centre = malloc(L_fixedpoints*sizeof(struct position));
	adjcentre = malloc(L_fixedpoints*sizeof(struct position));
	ext_centre = malloc(L_fixedpoints*sizeof(struct ext_position));
	fieldline_start.loc[0] = islandcenter[0]; fieldline_start.loc[1] = islandcenter[1];
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
	//printstructposition("fieldline[i]\n", &fieldline);
	number_turns=0;
	for (i=1; i<N_line+1; i++) {
		//RK4(&fieldline, varphi, dvarphi, coils, n_coils, n_segs);
		//printf("DEBUG: i=%d/%d\n", i, N_line);
		RK4step_lambdacirc_mutangent(&fieldline, &adjfieldline, varphi, dvarphi, Bfield_saved[i-1]);
		varphi += dvarphi;
		if (i%N_gridphi_fieldperiod==0) {
			centre_ind = (i / N_gridphi_fieldperiod);
			centre[centre_ind % L_fixedpoints] = fieldline;
			adjcentre[centre_ind % L_fixedpoints] = adjfieldline;
			//printf("varphi = %f\n", varphi);
			//printstructposition("fieldline[i]\n", &fieldline);
			rawangle_principal = atan2(fieldline.loc[1] - axis[1], fieldline.loc[0] - axis[0]);
			//angle_axis[centre_ind%L_fixedpoints] = remainder(rawangle, 2.0*M_PI);
			rawangle = rawangle_principal + number_turns*M_PI*2.0;
			if (rawangle - angle_axis[(centre_ind-1)%L_fixedpoints] - M_PI > 0.0) {
				number_turns -= 1;
			}
			else if (rawangle - angle_axis[(centre_ind-1)%L_fixedpoints] + M_PI < 0.0) {
				number_turns += 1;
			}
			if (i/N_gridphi_fieldperiod == 1) {
				number_turns=0;
			}
			rawangle = rawangle_principal + number_turns*M_PI*2.0;
			angle_axis[centre_ind%L_fixedpoints] = rawangle;
			//if (angle_axis[centre_ind%L_fixedpoints] < - M_PI) {
			//	angle_axis[centre_ind%L_fixedpoints] += M_PI;	
			//}
			//else if (angle_axis[centre_ind%L_fixedpoints] > M_PI) {
			//	angle_axis[centre_ind%L_fixedpoints] -= M_PI;	
			//}
			//angle_axis[centre_ind%L_fixedpoints] = atan2(fieldline.loc[1] - axis[1], fieldline.loc[0] - axis[0]) - ref_angle;
			printf("i=%d/%d\n", i, N_line);
			printf("i=%d/%d\n", i, N_line);
			printstructposition("centre", centre + ( (centre_ind-1) % L_fixedpoints ) );
			inverted = invert2x2(centre[(centre_ind-1) % L_fixedpoints].tangent, &detcentre);
			printf("i=%d/%d\n", i, N_line);
			adjinverted = invert2x2(adjcentre[(centre_ind-1) % L_fixedpoints].tangent, &adjdetcentre);
			//printmat("inverted", inverted, 2, 2);
			//ptarray[0] = &centre[centre_ind % L_fixedpoints].tangent[0][0];
			//ptarray[1] = &centre[centre_ind % L_fixedpoints].tangent[1][0];
			ext_centre[centre_ind % L_fixedpoints].loc[0] = fieldline.loc[0];
			ext_centre[centre_ind % L_fixedpoints].loc[1] = fieldline.loc[1];
			ext_centre[centre_ind % L_fixedpoints].part_tangent = multiply2x2(centre[centre_ind % L_fixedpoints].tangent, inverted, 2);
			ext_centre[centre_ind % L_fixedpoints].adj_part_tangent = multiply2x2(adjcentre[centre_ind % L_fixedpoints].tangent, adjinverted, 2);
			phidata[centre_ind % L_fixedpoints] = varphi;
			//printmat("ext_centre.part_tangent", ext_centre[centre_ind % L_fixedpoints].part_tangent, 2, 2);
			free(inverted[0]); free(inverted[1]);
			free(inverted);
		}
	}
	gsl_fit_linear(phidata, 1, angle_axis, 1, L_fixedpoints, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
	//c0 = remainder(13.9,3.0);
	printf("c0=%f, c1=%f\n\n\n", c0, c1);
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
			printf("angle_axis[%d]=%f, angle_axis[%d]=%f\n", centre_ind, angle_axis[centre_ind], (centre_ind-1)%L_fixedpoints, angle_axis[(centre_ind-1)%L_fixedpoints]);
			//if ( (clockwise == 0) && ( ( (angle_axis[centre_ind] > small) && (angle_axis[centre_ind-1] < - small) ) ) ) (*n_turns) +=1;	
			//else if ( (clockwise == 1) && ( ( (angle_axis[centre_ind] < -small) && (angle_axis[centre_ind-1] > small) ) ) ) (*n_turns) +=1;	
			*n_turns = number_turns;
		}
		//angle_axis[centre_ind] = atan2(ext_centre[(centre_ind+1)%L_fixedpoints].loc[1]-axis[1], ext_centre[(centre_ind+1)%L_fixedpoints].loc[0]-axis[0]);
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
	printf("n_turns = %d\n", *n_turns);
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
	//width = malloc(L_fixedpoints*sizeof(double));
			ext_centre[main_ind].angle = evals1[1];
			omega = m0_symmetry*evals1[1]/(2.0*M_PI*L_fixedpoints); // I have changed it compared to Cary and Hanson's paper
			q0 = (int) (m0_symmetry/(4.0*omega) - L_fixedpoints/2.0 + 0.5);
			ext_centre[main_ind].q0_index = q0;
			//printf("q0=%d\n", q0);
			//printf("m0=%d\n", m0_symmetry);
			//printf("omega=%f\n", omega);
			//printf("angle=%f\n", evals[1]);
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

double *islandwidthnew(struct field **Bfield_saved, struct ext_position *ext_centre, struct position *lambda_circ, struct position **lambda_tangent, struct position **mu_tangent, int m0_symmetry, int N_polmode, int L_fixedpoints, int N_gridphi_fieldperiod, double *param, int num_params) {
	// declarations
	//clock_t start = clock();
	int i=0, q0_ind, Q0_ind, Q_ind, q_ind, main_ind, sec_ind, N_line=0;
	struct position centre, lambda, mu, sperp;
	double varphi=0.0, chordlength, chordpluslength, *vecresult;
	double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*m0_symmetry);
	double *gradcircumference, *gradtangent, initialgradcirc=0.0, *gradwidth;
	double ***shape, *width, sum_matrix_elements, matrix_element;
	double circumference;
	shape = malloc(3*sizeof(double));
	//for (main_ind=0;main_ind<3;main_ind++) {
	//	shape[main_ind] = malloc((*n_coils)*sizeof(double));
	//	for (sec_ind=0; sec_ind<(*n_coils);sec_ind++) {
	//		shape[main_ind][sec_ind] = malloc((*(n_segs[main_ind]))*sizeof(double));
	//	}
	//}
	printf("L_fixedpoints=%d\n", L_fixedpoints);
	width = malloc(L_fixedpoints*sizeof(double));
	for (main_ind=0;main_ind<L_fixedpoints; main_ind++) {
		sum_matrix_elements = circumference = 0.0;		
		for (sec_ind=0;sec_ind<L_fixedpoints; sec_ind++) {
			circumference += pow( pow(ext_centre[main_ind].chord[0], 2.0) + pow(ext_centre[main_ind].chord[1], 2.0) , 0.5 );
			matrix_element = ext_centre[sec_ind].epar[0]*ext_centre[main_ind].sperp[sec_ind][0] + ext_centre[sec_ind].epar[1]*ext_centre[main_ind].sperp[sec_ind][1];
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
		printf("circumference = %10.10f\n", circumference);
		printf("sum_matrix_elements = %f\n", sum_matrix_elements);
		width[main_ind] = 2.0*L_fixedpoints*circumference/(M_PI*N_polmode*sum_matrix_elements);
		printf("width = %13.13f for index= %d\n", width[main_ind], main_ind);
	}

	gradtangent = malloc(L_fixedpoints*sizeof(double));
	gradwidth = malloc(L_fixedpoints*sizeof(double));
	centre.loc[0] = ext_centre[0].loc[0]; centre.loc[1] = ext_centre[0].loc[1];
	//lambda_circ.loc[0] = lambda[0]; lambda_circ.loc[1] = lambda[1];
	centre.tangent = set_identity();
	//lambda_circ.tangent = set_identity();
	gradcircumference = &initialgradcirc;
	N_line = L_fixedpoints*N_gridphi_fieldperiod;
	varphi = 0.0;
	for (i=0; i<N_line; i++)
	{
		if (i%N_gridphi_fieldperiod==0)
		{
			main_ind = i/N_gridphi_fieldperiod;
			chordlength= pow( pow(ext_centre[main_ind].chord[0], 2.0) + pow(ext_centre[main_ind].chord[1], 2.0), 0.5);
			chordpluslength= pow( pow(ext_centre[main_ind].chordplus[0], 2.0) + pow(ext_centre[main_ind].chordplus[1], 2.0), 0.5);
			lambda_circ->loc[0] += (ext_centre[main_ind].chord[0]/chordlength - ext_centre[main_ind].chordplus[0]/chordpluslength);
			lambda_circ->loc[1] += (ext_centre[main_ind].chord[1]/chordlength - ext_centre[main_ind].chordplus[1]/chordpluslength);
		}
		RK4step_gradcirc(gradcircumference, &centre, lambda_circ, varphi, dvarphi, Bfield_saved[i], param, 3); /// currently wrong
		varphi += dvarphi;
		//printf("gradcircumference = %f\n", *gradcircumference);
	}
	//printf("varphi = %f\n", varphi);
	//printstructposition("lambda",lambda_circ);
	//printstructposition("Xp",&centre);
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
		//printf("main_ind=%d/%d\n", main_ind, L_fixedpoints);
		varphi = main_ind*2.0*M_PI/m0_symmetry;
		gradtangent[main_ind] = 0.0; 
		centre.loc[0] = ext_centre[main_ind].loc[0]; centre.loc[1] = ext_centre[main_ind].loc[1];
		centre.tangent[0][0] = centre.tangent[1][1] = 1.0;
		centre.tangent[0][1] = centre.tangent[1][0] = 0.0;
		//vecresult = multiply(centre.tangent, ext_centre[main_ind].eperp);
		//free(vecresult);
		sperp.loc[0] = ext_centre[main_ind].eperp[0]; sperp.loc[1] = ext_centre[main_ind].eperp[1];
		for (Q_ind=0;Q_ind<Q0_ind;Q_ind++) {
			lambda = lambda_tangent[main_ind][Q_ind];
			mu = mu_tangent[main_ind][Q_ind];
			//printstructposition("lambda", &lambda);
			//printstructposition("mu", &mu);
			centre.loc[0] = ext_centre[main_ind].loc[0]; centre.loc[1] = ext_centre[main_ind].loc[1];
			vecresult = multiply(centre.tangent, ext_centre[main_ind].eperp);
			sperp.loc[0] = vecresult[0]; sperp.loc[1] = vecresult[1];
			free(vecresult);
			//printf("sperp=%f\n", sperp.tangent[1][1]);
			for (i=0; i<N_line; i++)
			{
				if (i%N_gridphi_fieldperiod==0)
				{
					sec_ind = i/N_gridphi_fieldperiod;
					//printf("varphi = %f\n", varphi);
					//printstructposition("Xp", &centre);
					//printf("gradtangent = %f\n", gradtangent[sec_ind]);
					//if (i/N_gridphi_fieldperiod==0) lambda_centre->loc[0] += 1.0;
					//printstructposition("lambda\n", lambda);
					//if (sec_ind==0) 
					//{
					//	printf("main_ind=%d, Q_ind=%d, sec_ind=%d\n", main_ind, Q_ind, sec_ind); 
					//	printf("gradtangent = %f\n", gradtangent[main_ind]);
					//	printstructposition("mu", &mu);
					//	printstructposition("lambda", &lambda);
					//}
					q_ind = Q_ind*L_fixedpoints + sec_ind;
					//check if makes sense
					if ( q_ind >= q0_ind && q_ind < q0_ind + L_fixedpoints) {
						mu.loc[0] += (ext_centre[(main_ind+sec_ind)%L_fixedpoints].epar[0]);
						mu.loc[1] += (ext_centre[(main_ind+sec_ind)%L_fixedpoints].epar[1]);
					}
				}
				//RK4_adjtangent(&centre, &lambda, &sperp, &mu, varphi, dvarphi, coils, n_coils, n_segs);
				RK4step_gradtangent(gradtangent+main_ind, &centre, &lambda, &sperp, &mu, varphi, dvarphi, Bfield_saved[i], param, num_params);
				varphi += dvarphi;
			}
		}
		//printf("varphi = %f\n", varphi);
		//printstructposition("lambda", &lambda);
		//printstructposition("Xp",&centre);
	}
	for (main_ind=0; main_ind<L_fixedpoints; main_ind++) {
		gradwidth[main_ind] = width[main_ind]*((*gradcircumference)/circumference - gradtangent[main_ind]/sum_matrix_elements);
		printf("gradcirc = %f\n", *gradcircumference);
		printf("gradtangent = %f\n", gradtangent[main_ind]);
		printf("gradwidth = %f\n", gradwidth[main_ind]);
	}

	return width;
}

double *calc_islandwidth(double *circ, double *Sigma, struct ext_position *ext_fieldline, int m0_symmetry, int L_fixedpoints, int pol_mode, int N_gridphi_fieldperiod) {
	int main_index, centre_index;
	double *wperp, sum_matrix_elements, matrix_element;
	double circumference2=ext_fieldline[0].circumference, circumference;
	double matrix_element_oldway, sum_matrix_elements_oldway=0.0;
	printf("L_fixedpoints=%d\n", L_fixedpoints);
	wperp = malloc(L_fixedpoints*sizeof(double));
	for (main_index=0;main_index<L_fixedpoints; main_index++) {
		sum_matrix_elements = sum_matrix_elements_oldway = circumference = 0.0;		
		for (centre_index=0;centre_index<L_fixedpoints; centre_index++) {
			circumference += pow( pow(ext_fieldline[main_index].chord[0], 2.0) + pow(ext_fieldline[main_index].chord[1], 2.0) , 0.5 );
			matrix_element_oldway = inner(ext_fieldline[centre_index].epar, 
					       ext_fieldline[main_index].long_tangent[centre_index], 
					       ext_fieldline[main_index].eperp);	
			matrix_element = ext_fieldline[centre_index].epar[0]*ext_fieldline[main_index].sperp[centre_index][0] + ext_fieldline[centre_index].epar[1]*ext_fieldline[main_index].sperp[centre_index][1];
			sum_matrix_elements_oldway += fabs(matrix_element_oldway);
			//sum_matrix_elements += fabs(matrix_element);
			//printf("matrix_element=%f\n", matrix_element);
			sum_matrix_elements += matrix_element;
			//printf("sum_matrix_element=%f\n", sum_matrix_elements);
			//sum_matrix_elements += fabs(matrix_element)/sin(ext_fieldline[main_index].angle*ext_fieldline[main_index].q0_index/L_fixedpoints); //This has the sin piece in the denominator, which can be safely set to unity
			//circumference += sqrt(pow(ext_fieldline[(centre_index+1)%L_fixedpoints].loc[0] 
			//		        - ext_fieldline[centre_index].loc[0], 2.0)
			//	    	    + pow(ext_fieldline[(centre_index+1)%L_fixedpoints].loc[1] 
			//			- ext_fieldline[centre_index].loc[1], 2.0));
			//printf("circumference = %f\nmatrix_element = %f\n", circumference, matrix_element);
			//printf("other circumference = %f\n", ext_fieldline[centre_index].circumference);
		}
		//circumference = 2.0*M_PI*0.2102;
		//circumference = 2.0*M_PI*0.2094;
		printf("circumference = %f\n", circumference);
		printf("circumference_oldway = %f\n", circumference2);
		printf("sum_matrix_elements = %f\n", sum_matrix_elements);
		printf("sum_matrix_elements_oldway = %f\n", sum_matrix_elements_oldway);
		wperp[main_index] = 2.0*L_fixedpoints*circumference/(M_PI*pol_mode*sum_matrix_elements);
		printf("width = %f for index= %d\n", wperp[main_index], main_index);
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
	// declarations
	//clock_t start = clock();
	int centre_ind; //i=0
	struct position lambda, lambda_start, deltalambda; //, fieldline_start, fieldline;
	double varphi=0.0, **jumptocentre=calloc(1,sizeof(double*)); //, **testunity;
	//double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*m0_symmetry);
	double **matrix, *vecresult, **result;
	double **pdeltalambda, det_tangent;
	double **inverseTminusI, **inverseT;
	double error = 1.0, errorlimit = 0.000000000001;
	double factor = 1.0, chordlength, chordpluslength;
	int N_line=0;
	N_line = L_fixedpoints*N_gridphi_fieldperiod;
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
				//chordpluslength = chordlength; // 
				printf("ext_centre.chord = (%f, %f)\n", ext_centre[centre_ind].chord[0], ext_centre[centre_ind].chord[1]);
				printf("ext_centre.chordplus = (%f, %f)\n", ext_centre[centre_ind].chordplus[0], ext_centre[centre_ind].chordplus[1]);
				printf("chordlength = %f, chordpluslength = %f\n", chordlength, chordpluslength);
				lambda.loc[0] += (ext_centre[centre_ind].chord[0]/chordlength - ext_centre[centre_ind].chordplus[0]/chordpluslength);
				lambda.loc[1] += (ext_centre[centre_ind].chord[1]/chordlength - ext_centre[centre_ind].chordplus[1]/chordpluslength);
				vecresult = multiply(ext_centre[(centre_ind+1)%L_fixedpoints].adj_part_tangent, lambda.loc);
				result = multiply2x2(ext_centre[(centre_ind+1)%L_fixedpoints].adj_part_tangent, lambda.tangent, 2);
				lambda.loc[0] = vecresult[0]; lambda.loc[1] = vecresult[1];
				lambda.tangent = result;
				printf("lambda = (%f, %f)\n", lambda.loc[0], lambda.loc[1]);
				//printmat("lambda start\n", lambda_start.tangent, 2, 2);
		}
		///// temporary (the above is much faster and now works, so I can permanently delete what's below eventually)
		/////for (i=0; i<N_line; i++)
		/////{
		/////	if (i%N_gridphi_fieldperiod==0)
		/////	{
		/////		centre_ind = i/N_gridphi_fieldperiod;
		/////		printf("varphi = %f\n", varphi);
		/////		//printstructposition("lambda\n", lambda);
		/////		//if (centre_ind==0) lambda->loc[0] += 1.0;
		/////		chordlength= pow( pow(ext_centre[centre_ind].chord[0], 2.0) + pow(ext_centre[centre_ind].chord[1], 2.0), 0.5);
		/////		chordpluslength= pow( pow(ext_centre[centre_ind].chordplus[0], 2.0) + pow(ext_centre[centre_ind].chordplus[1], 2.0), 0.5);
		/////		lambda->loc[0] += (ext_centre[centre_ind].chord[0]/chordlength - ext_centre[centre_ind].chordplus[0]/chordpluslength);
		/////		lambda->loc[1] += (ext_centre[centre_ind].chord[1]/chordlength - ext_centre[centre_ind].chordplus[1]/chordpluslength);
		/////	}
		/////	RK4_adjsimple(&fieldline, lambda, varphi, dvarphi, coils, n_coils, n_segs);
		/////	varphi += dvarphi;
		/////}
		deltalambda = addstructs(1.0, &lambda, -1.0, &lambda_start); 
		//printstructposition("lambda\n", &lambda);
		//printstructposition("deltalambda", &deltalambda);
		//printf("fieldline->tangent[1][1] = %f\n", fieldline->tangent[1][0]);
		inverseTminusI = invert2x2(deltalambda.tangent, &det_tangent);
		inverseT = invert2x2(lambda.tangent, &det_tangent);
		//printf("varphi = %f\n", varphi);
		//printstructposition("fieldline", lambda);
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
		//printstructposition("deltalambda", &deltalambda);
		//printmat("jumptocentre", jumptocentre, 2, 1);
		//printstructposition("lambda_start", &lambda_start);
		//printstructposition("lambda", &lambda);
		lambda.loc[0] = lambda_start.loc[0] - jumptocentre[0][0]; 
		lambda.loc[1] = lambda_start.loc[1] - jumptocentre[1][0];
		lambda.tangent[0][0] = 1.0; lambda.tangent[0][1]=0.0; lambda.tangent[1][0]=0.0; lambda.tangent[1][1]=1.0;
		error = sqrt(pow(deltalambda.loc[0], 2.0) + pow(deltalambda.loc[1], 2.0)); 
		//printf("%f\n", deltalambda.loc[0]);
		//printf("%f\n", deltalambda.loc[1]);
		//free(jumptocentre);
	} while ((fabs(deltalambda.loc[0]) > errorlimit) || (fabs(deltalambda.loc[1]) > errorlimit));
	//clock_t int3 = clock();
	printstructposition("lambda", &lambda);
	return lambda;
}

struct position **solve_mu_tangent(struct ext_position *ext_centre, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod) {
	// declarations
	//clock_t start = clock();
	int centre_ind, main_ind, q0_ind, Q0_ind, Q_ind, q_ind;
	struct position mu_start, fieldline_start, fieldline, **mu;
	//double **jumptocentre=calloc(1,sizeof(double*)); //, **testunity;
	//double varphi = 0.0, dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*m0_symmetry);
	double **result, *vecresult, **inverted, det;// varphi=0.0;
	int N_line;
	N_line = L_fixedpoints*N_gridphi_fieldperiod;
	mu = malloc(L_fixedpoints*sizeof(struct position*));
	mu_start.tangent = set_identity();
	for (main_ind=0; main_ind<L_fixedpoints; main_ind++)  {
		q0_ind = ext_centre[main_ind].q0_index; 
		Q0_ind = q0_ind/L_fixedpoints + 2;
		mu[main_ind] = malloc((Q0_ind)*sizeof(struct position));
		mu_start.loc[0] = 0.0;  mu_start.loc[1] = 0.0; 
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

struct position **solve_lambda_tangent(struct field **Bfield_saved, struct ext_position *ext_centre, struct position **mu, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod) {
	// declarations
	//clock_t start = clock();
	int N_line, i=0, centre_ind, main_ind, q0_ind, Q0_ind, Q_ind, q_ind, ind, count = 0;
	struct position fieldline_start, fieldline, sperp, muvar, lambdavar, deltalambda, **lambda, lambdain; //, fieldline_start, fieldline;
	double varphi=0.0, **jumptocentre=calloc(1,sizeof(double*)); //, **testunity;
	//double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*m0_symmetry);
	double **pdeltalambda, det_tangent;
	double **inverseTminusI, **inverseT;
	double error = 1.0, errorlimit = 0.000000001;
	double **matrix, **matrix2, *vector;
	double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*m0_symmetry), factor = 1.0;
	N_line = L_fixedpoints*N_gridphi_fieldperiod;
	printf("N_line = %d\n", N_line);
	printf("field_periods = %d\n", m0_symmetry);
	lambda = malloc(L_fixedpoints*sizeof(struct position)); 
	matrix2 = set_identity();
	sperp.tangent = set_identity();
	lambdavar.tangent = set_identity();
	q0_ind = ext_centre[0].q0_index; 
	Q0_ind = q0_ind/L_fixedpoints + 2;
	printf("Q0_ind=%d\n", Q0_ind);
	fieldline.tangent = set_identity();
	for (main_ind=0;main_ind<L_fixedpoints;main_ind++) {
		lambda[main_ind] = malloc(Q0_ind*sizeof(struct position));
		fieldline_start.loc[0] = ext_centre[main_ind].loc[0]; fieldline_start.loc[1] = ext_centre[main_ind].loc[1];
		fieldline.tangent[0][0] = fieldline.tangent[1][1] = 1.0; 
		fieldline.tangent[1][0] = fieldline.tangent[0][1] = 0.0; 
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
				//printstructposition("lambda", lambda[main_ind]+Q_ind);
				//lambda[main_ind][Q_ind].tangent[0][0] = 1.0; lambda[main_ind][Q_ind].tangent[1][0] = 0.0; 
				//lambda[main_ind][Q_ind].tangent[0][1] = 0.0; lambda[main_ind][Q_ind].tangent[1][1] = 1.0; 
				fieldline.loc[0] = fieldline_start.loc[0];
				fieldline.loc[1] = fieldline_start.loc[1];
				fieldline.tangent[0][0] = fieldline.tangent[1][1] = 1.0; 
				fieldline.tangent[1][0] = fieldline.tangent[0][1] = 0.0; 
				//printstructposition("fieldline", &fieldline);
				//printstructposition("lambdain", &lambdain);
				muvar = mu[main_ind][Q_ind];
				//varphi = 2.0*M_PI*(Q0_ind+1)*L_fixedpoints/m0_symmetry;
				varphi = main_ind*2.0*M_PI/m0_symmetry;
				for (i=0; i<N_line; i++) {
					if (i%N_gridphi_fieldperiod==0) {
						//printf("N_line = %d\n", N_line);
						//printf("count = %d\n", count);
						//printf("varphi = %f\n", varphi);
						//printstructposition("fieldline", &fieldline);
						//printstructposition("lambdain", &lambdain);
						//printstructposition("sperp", &sperp);
						centre_ind = i/N_gridphi_fieldperiod;
						q_ind = Q_ind*L_fixedpoints + centre_ind;
						
						//if ( (count == 1) && (centre_ind == 1000) ) {
						//	printf("main_ind=%d, Q_ind=%d, centre_ind=%d\n", main_ind, Q_ind, centre_ind);
						//	printf("count = %d\n", count);
						//	printstructposition("muvar", &muvar);
						//	printstructposition("lambdain", &lambdain);
						//}
						if ( q_ind >= q0_ind && q_ind < q0_ind + L_fixedpoints) {
							muvar.loc[0] += (ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[0]);
							muvar.loc[1] += (ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[1]);
						}
					}
					RK4step_lambdatangent(&fieldline, &lambdain, &sperp, &muvar, varphi, dvarphi, Bfield_saved[(main_ind*N_gridphi_fieldperiod + i) % N_line]);
					varphi += dvarphi;
				}
				//deltalambda = addstructs(-1.0, &lambdavar, 1.0, lambda[main_ind]+Q_ind); 
				//printf("varphi = %f\n", varphi);
				//printstructposition("fieldline", &fieldline);
				//printstructposition("lambdain", &lambdain);
				deltalambda = addstructs(1.0, &lambdavar, -1.0, &lambdain); 
				//printstructposition("deltalambda", &deltalambda);
				//printstructposition("lambdaQ", lambda[main_ind]+Q_ind);
				//printstructposition("lambdavar", &lambdavar);
				//printf("this above should always be the identity\n");
				inverseTminusI = invert2x2(deltalambda.tangent, &det_tangent);
				inverseT = invert2x2(lambdain.tangent, &det_tangent);
				//printf("det(tangent)=%f\n", det_tangent);
				pdeltalambda = calloc(2,sizeof(double*)); 
				*pdeltalambda     = &deltalambda.loc[0];
				*(pdeltalambda+1) = &deltalambda.loc[1];
				matrix = multiply2x2(inverseTminusI, deltalambda.tangent, 2);
				//printmat("matrix", matrix, 2, 2);
				jumptocentre = multiply2x2(inverseTminusI, pdeltalambda, 1);
				jumptocentre[0][0] *= (factor);
				jumptocentre[1][0] *= (factor);
				//printmat("jumptocentre", jumptocentre, 2, 1);
				lambdain.loc[0] = lambdavar.loc[0] - jumptocentre[0][0]; 
				lambdain.loc[1] = lambdavar.loc[1] - jumptocentre[1][0];
				lambdain.tangent[0][0] = 1.0; lambdain.tangent[0][1]=0.0; lambdain.tangent[1][0]=0.0; lambdain.tangent[1][1]=1.0;
				error = sqrt(pow(deltalambda.loc[0], 2.0) + pow(deltalambda.loc[1], 2.0)); 
				//free(jumptocentre);
				//printf("count=%d\n", count);
				count++;
				lambda[main_ind][Q_ind].loc[0] = lambdain.loc[0];
				lambda[main_ind][Q_ind].loc[1] = lambdain.loc[1];
			} while ( ((fabs(deltalambda.loc[0]/lambdain.loc[0]) > errorlimit) || (fabs(deltalambda.loc[1]/lambdain.loc[1]) > errorlimit)) );
		}
	}
	return lambda;
}

void solve_gradcirc(double *gradcircumference, struct field **Bfield_island, struct ext_position *ext_centre, struct position *lambda_circ, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod, double *param, int num_params) {
	// declarations
	//clock_t start = clock();
	int i=0, main_ind, N_line=0;
	struct position centre;
	double varphi=0.0, chordlength, chordpluslength, **matrix2;
	double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*m0_symmetry);
	//double *gradtangent;

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
	//printf("field_periods = %d\n", m0_symmetry);
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
		RK4step_gradcirc(gradcircumference, &centre, lambda_circ, varphi, dvarphi, Bfield_island[i], param, num_params);
		varphi += dvarphi;
		//printf("gradcircumference = %f\n", gradcircumference[0]);
	}
	//printf("varphi = %f\n", varphi);
	//printstructposition("lambda",lambda_circ);
	//printstructposition("Xp",&centre);
	//printf("gradcircumference = %f\n", *gradcircumference);

	//centre.loc[0] = ext_centre[0].loc[0]; centre.loc[1] = ext_centre[0].loc[1];
	//centre.tangent[0][0] = centre.tangent[1][1] = 1.0;
	//centre.tangent[0][1] = centre.tangent[1][0] = 0.0;
	////centre.loc[0] = ext_centre[0].loc[0]; centre.loc[1] = ext_centre[0].loc[1];
	//////lambda_centre.loc[0] = lambda[0]; lambda_centre.loc[1] = lambda[1];
	////centre.tangent = set_identity();
	////lambda_centre.tangent = set_identity();
	//q0_ind = ext_centre[0].q0_index; 
	//Q0_ind = q0_ind/L_fixedpoints + 2;
	//sperp.tangent = set_identity();
	//for (main_ind=0; main_ind<L_fixedpoints; main_ind++) {
	//	varphi = main_ind*2.0*M_PI/m0_symmetry;
	//	initial = 0.0;
	//	gradtangent[main_ind] = initial; 
	//	centre.tangent[0][0] = centre.tangent[1][1] = 1.0;
	//	centre.tangent[0][1] = centre.tangent[1][0] = 0.0;
	//	for (Q_ind=0;Q_ind<Q0_ind;Q_ind++) {
	//		lambda = lambda_tangent[main_ind][Q_ind];
	//		mu = mu_tangent[main_ind][Q_ind];
	//		//printstructposition("lambda", &lambda);
	//		//printstructposition("mu", &mu);
	//		matrix2[0][0] = matrix2[1][1] = 1.0;
	//		matrix2[0][1] = matrix2[1][0] = 0.0;
	//		for (ind=0;ind<Q_ind; ind++) {
	//			multiply2x2reassign(ext_centre[main_ind].full_tangent, matrix2, 2);
	//		}
	//		vecresult = multiply(matrix2, ext_centre[main_ind].eperp);
	//		sperp.loc[0] = vecresult[0];  sperp.loc[1] = vecresult[1];
	//		free(vecresult);
	//		sperp.tangent[0][0] = sperp.tangent[1][1] = 1.0; 
	//		sperp.tangent[1][0] = sperp.tangent[0][1] = 0.0; 
	//		centre.loc[0] = ext_centre[main_ind].loc[0]; centre.loc[1] = ext_centre[main_ind].loc[1];
	//		//vecresult = multiply(centre.tangent, ext_centre[main_ind].eperp);
	//		//sperp.loc[0] = vecresult[0]; sperp.loc[1] = vecresult[1];
	//		//free(vecresult);
	//		//printf("sperp=%f\n", sperp.tangent[1][1]);
	//		//centre->tangent[0][0] = 1.0 ; centre->tangent[0][1] = 0.0;
	//		//centre->tangent[1][0] = 0.0 ; centre->tangent[0][1] = 1.0;
	//		for (i=0; i<N_line; i++)
	//		{
	//			if (i%N_gridphi_fieldperiod==0)
	//			{
	//				sec_ind = i/N_gridphi_fieldperiod;
	//				//printf("varphi = %f\n", varphi);
	//				//printstructposition("Xp", &centre);
	//				//printf("gradtangent = %f\n", gradtangent[sec_ind]);
	//				//if (i/N_gridphi_fieldperiod==0) lambda_centre->loc[0] += 1.0;
	//				//printstructposition("lambda\n", lambda);
	//				//if (sec_ind==0) 
	//				//{
	//				//	printf("main_ind=%d, Q_ind=%d, sec_ind=%d\n", main_ind, Q_ind, sec_ind); 
	//				//	printf("gradtangent = %f\n", gradtangent[main_ind]);
	//				//	printstructposition("mu", &mu);
	//				//	printstructposition("lambda", &lambda);
	//				//}
	//				q_ind = Q_ind*L_fixedpoints + sec_ind;
	//				if ( q_ind >= q0_ind && q_ind < q0_ind + L_fixedpoints) {
	//					mu.loc[0] += (ext_centre[(main_ind+sec_ind)%L_fixedpoints].epar[0]);
	//					mu.loc[1] += (ext_centre[(main_ind+sec_ind)%L_fixedpoints].epar[1]);
	//				}
	//			}
	//			//RK4_adjtangent(&centre, &lambda, &sperp, &mu, varphi, dvarphi, coils, n_coils, n_segs);
	//			RK4step_gradtangent(gradtangent+main_ind, &centre, &lambda, &sperp, &mu, varphi, dvarphi, coils, n_coils, n_segs);
	//			varphi += dvarphi;
	//		}
	//	}
	//	//printf("varphi = %f\n", varphi);
	//	//printstructposition("lambda", &lambda);
	//	//printstructposition("Xp",&centre);
	//}
	//for (sec_ind=0; sec_ind<L_fixedpoints; sec_ind++) {
	//	printf("gradtangent = %f\n", gradtangent[sec_ind]);
	//}

	return;
}

void solve_gradtangent(double **number, struct field **Bfield_island, struct ext_position *ext_centre, struct position **lambdaQ, struct position **muQ, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod, double *param, int num_params) {

	int i=0, ind, centre_ind, q0_ind, Q0_ind, Q_ind, q_ind, main_ind;
	struct position fieldline, lambda, mu, sperp;
	double varphi=0.0, *vecresult, **matrix2;
	double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*m0_symmetry);
	//double **number;
	int N_line;

	matrix2 = set_identity();
	N_line = L_fixedpoints*N_gridphi_fieldperiod;
	fieldline.loc[0] = ext_centre[0].loc[0]; fieldline.loc[1] = ext_centre[0].loc[1];
	fieldline.tangent = set_identity();
	//centre.loc[0] = ext_centre[0].loc[0]; centre.loc[1] = ext_centre[0].loc[1];
	////lambda_centre.loc[0] = lambda[0]; lambda_centre.loc[1] = lambda[1];
	//centre.tangent = set_identity();
	//lambda_centre.tangent = set_identity();
	q0_ind = ext_centre[0].q0_index; 
	Q0_ind = q0_ind/L_fixedpoints + 2;
	//number = malloc(L_fixedpoints*sizeof(double));
	sperp.tangent = set_identity();
	//printf("N_line = %d\n", N_line);
	//printf("field_periods = %d\n", m0_symmetry);
	for (main_ind=0; main_ind<L_fixedpoints; main_ind++) {
		varphi = main_ind*2.0*M_PI/m0_symmetry;
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
						mu.loc[0] += (ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[0]);
						mu.loc[1] += (ext_centre[(main_ind+centre_ind)%L_fixedpoints].epar[1]);
					}
				}
				RK4step_gradtangent(number[main_ind], &fieldline, &lambda, &sperp, &mu, varphi, dvarphi, Bfield_island[(main_ind*N_gridphi_fieldperiod + i) % N_line], param, num_params);
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


