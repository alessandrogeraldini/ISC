//Author: Alessandro Geraldini
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "isc.h"
#include <gsl/gsl_fit.h>

void iotaprofile(double RRin, double r_interval, int n_points, int m0_symmetry, int N_gridphi_per_field_period, double *minor_radius, double *iota, double ***coils, int n_coils, int *n_segs) {
	/* declarations */
	clock_t initial_t = clock();
	//int N_gridphi_per_field_period=20, field_periods=3; // now in separate file as defined constants
	int i=0, i_periods=0, ind, include_centre=0, number_field_periods = 10000;
	int number_turns = 0, N_gridphi_tor = N_gridphi_per_field_period*m0_symmetry;
	struct position *Xp=calloc(1,sizeof(struct position)), *fieldcentre, startcentre;
	double varphi=0.0, ZZin=0.0, deltar = r_interval/n_points;
	double dvarphi = 2.0*M_PI/(N_gridphi_per_field_period*m0_symmetry);
	double angle_round_axis=0.0, angle_round_axis_principal, angle_old, z_position, r_position;
	double angle[number_field_periods*N_gridphi_per_field_period], phidata[number_field_periods*N_gridphi_per_field_period];
	double c0=0.0, c1=0.0, cov00=0.0, cov01=0.0, cov11=0.0, sumsq=0.0;
	double **evec=malloc(2*sizeof(double)), eval[2], det, trace, iota_axis;
	struct field **Bfield_axis;
	FILE *Poincare=NULL, *IOTA=NULL;
	printf("field periods = %d\n", m0_symmetry);
	printf("N points per field period = %d\n", N_gridphi_per_field_period);
	startcentre.loc[0] = RRin; startcentre.loc[1] = ZZin;
	startcentre.tangent = malloc(2*sizeof(double));
	startcentre.tangent[0] = malloc(2*sizeof(double)); startcentre.tangent[1] = malloc(2*sizeof(double));
	startcentre.tangent[0][0] = 1.0; startcentre.tangent[0][1] = 0.0; startcentre.tangent[1][0] = 0.0; startcentre.tangent[1][1] = 1.0;
	Xp->tangent = malloc(2*sizeof(double));
	Xp->tangent[0] = malloc(2*sizeof(double)); Xp->tangent[1] = malloc(2*sizeof(double));
	evec[0]=malloc(2*sizeof(double)); evec[1]=malloc(2*sizeof(double));
	if ((Poincare = fopen("Poincare_file.txt", "w"))==NULL) printf("Could not open Poincare.txt\n");
	if ((IOTA = fopen("iota_file.txt", "w"))==NULL) printf("Could not open iota.txt\n");
	/* find the magnetic axis */
	if (include_centre == 1) {
		Bfield_axis = malloc(N_gridphi_tor*sizeof(struct field));
		for (ind=0; ind<N_gridphi_tor; ind++) {
			Bfield_axis[ind] = malloc(4*sizeof(struct field));
		}
		fieldcentre = solve_magneticaxis(coils, n_coils, n_segs, Bfield_axis, &startcentre, N_gridphi_tor);
		linalg2x2(startcentre.tangent, evec, eval, &det, &trace);
		printf("evec=%f\n", evec[0][0]);
		iota_axis = eval[1];
		fprintf(Poincare, "%f %f %f\n", 0.0, fieldcentre->loc[0], fieldcentre->loc[1]);
	}
	/* loop over different initial values of (R,Z) */
	for (ind=0;ind<n_points;ind+=1) {
		minor_radius[ind] = ind*deltar;
		printf("here?\n");
		Xp->loc[0] = startcentre.loc[0] + ind*deltar; Xp->loc[1] = startcentre.loc[1];
		Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		varphi = 0.0;
		printf("initial (R,Z)=(%f,%f)\n", Xp->loc[0], Xp->loc[1]);
		angle_old = 0.0;
		number_turns = 0;
		for (i_periods=0; i_periods<number_field_periods; i_periods++)  {
			fprintf(Poincare, "%f %f %f\n", varphi, Xp->loc[0], Xp->loc[1]);
			//printf("at varphi=%f, (R,Z)=(%f,%f)\n", varphi, Xp->loc[0], Xp->loc[1]);
		for (i=0; i<N_gridphi_per_field_period; i++)  {
			if (include_centre == 1) {
				r_position = Xp->loc[0] - (fieldcentre+i)->loc[0];
				z_position = Xp->loc[1] - (fieldcentre+i)->loc[1];
			}
			else {
				r_position = Xp->loc[0] - 1.0;
				z_position = Xp->loc[1];
			}
			//printf("r_position=%f\tz_position=%f\n", r_position, z_position);
			angle_round_axis_principal = atan2(z_position, r_position);
			angle_round_axis = angle_round_axis_principal + number_turns*M_PI*2.0;
			if (angle_round_axis - angle_old - M_PI > 0.0) {
				number_turns -= 1;
			}
			else if (angle_round_axis - angle_old + M_PI < 0.0) {
				number_turns += 1;
			}
			angle_round_axis = angle_round_axis_principal + number_turns*M_PI*2.0;
			angle[i_periods*N_gridphi_per_field_period + i] = angle_round_axis;
			phidata[i_periods*N_gridphi_per_field_period + i] = varphi;
			//printf("%f\t%f\n", angle[i], phidata[i]);
			angle_old = angle_round_axis;
			RK4step(Xp, varphi, dvarphi, coils, n_coils, n_segs, Bfield_axis[i]);
			//printstructposition("fieldcentre", fieldcentre+i);
			varphi += dvarphi;
		}
		}
		gsl_fit_linear(phidata, 1, angle, 1, number_field_periods*N_gridphi_per_field_period, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
		//gsl_fit_linear(phidata, 1, angle, 1, 32768, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
		iota[ind] = c1;
		printf("\niota=%f\tintercept=%f\n", iota[ind], c0);
		fprintf(IOTA, "%f %f\n", minor_radius[ind], iota[ind]);
	}
	fclose(Poincare);
	fclose(IOTA);
	free(evec[0]); free(evec[1]);
	free(Xp->tangent[0]); free(Xp->tangent[1]);
	free(startcentre.tangent[0]); free(startcentre.tangent[1]);
	free(evec);
	free(Xp->tangent);
	free(Xp);
	clock_t final_t = clock();
	printf("Time after obtaining iota profile and Poincare plot is: %f\n", (double) (final_t-initial_t)/CLOCKS_PER_SEC);
}
