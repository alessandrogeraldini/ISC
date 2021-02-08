//Author: Alessandro Geraldini
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "isc.h"
#include <gsl/gsl_fit.h>
#define error 1.0e-10

//void iotaprofile(double RRin, double r_interval, int n_points, int m0_fieldperiods, int N_gridphi_fieldperiod, double *minor_radius, double *iota, double ***coils, int n_coils, int *n_segs);
void iotaprofile(struct position axis, int m0_fieldperiods, int N_gridphi_fieldperiod, double *rmin, double *zmin, double *iota, struct fieldparams allparams) {
	/* declarations */
	clock_t initial_t = clock();
	//int N_gridphi_fieldperiod=20, field_periods=3; // now in separate file as defined constants
	int i=0, i_periods=0, ind, include_centre=1, number_field_periods = 1000;
	int number_turns = 0, N_gridphi_tor = N_gridphi_fieldperiod*m0_fieldperiods;
	int countlines=0, indices = number_field_periods*N_gridphi_fieldperiod, sizestring=100, nn;
	struct position Xp, *fieldcentre;
	double varphi=0.0;
	double dvarphi = 2.0*M_PI/(N_gridphi_fieldperiod*m0_fieldperiods);
	double angle_round_axis=0.0, angle_round_axis_principal, angle_old, z_position, r_position;
	double *angle=malloc(indices*sizeof(double)), *phidata=malloc(indices*sizeof(double));
	double c0=0.0, c1=0.0, cov00=0.0, cov01=0.0, cov11=0.0, sumsq=0.0;
	double **evec=malloc(2*sizeof(double)), eval[2], det, trace, iota_axis;
	double *line_double;
	struct field **Bfield_axis;
	char string[sizestring];
	FILE *gridPoincare=NULL, *Poincare=NULL, *IOTA=NULL;
	printf("field periods = %d\n", m0_fieldperiods);
	printf("N points per field period = %d\n", N_gridphi_fieldperiod);
	Xp.tangent = malloc(2*sizeof(double));
	Xp.tangent[0] = malloc(2*sizeof(double)); Xp.tangent[1] = malloc(2*sizeof(double));
	evec[0]=malloc(2*sizeof(double)); evec[1]=malloc(2*sizeof(double));
	if ((gridPoincare = fopen("gridPoincare.txt", "r"))==NULL) printf("Could not open gridPoincare.txt\n");
	if ((Poincare = fopen("Poincare.txt", "w"))==NULL) printf("Could not open Poincare.txt\n");
	if ((IOTA = fopen("iota.txt", "w")) == NULL) printf("Could not open iota.txt\n");

	while (fgets(string, sizestring, gridPoincare) != NULL) {
		countlines+=1;
	}
	printf("countlines=%d\n", countlines);
	rmin = malloc(countlines*sizeof(double));
	zmin = malloc(countlines*sizeof(double));
	iota = malloc(countlines*sizeof(double));
	rewind(gridPoincare);
	countlines = 0;
	while (fgets(string, sizestring, gridPoincare) != NULL) {
		line_double = linetodata(string, &nn);
		rmin[countlines] = line_double[0];
		zmin[countlines] = line_double[1];
		printf("rmin[%d] = %f, zmin[%d] = %f\n", countlines, rmin[countlines], countlines, zmin[countlines]);
		countlines+=1;
	}
	/* find the magnetic axis */
	if (include_centre == 1) {
		Bfield_axis = malloc(N_gridphi_fieldperiod*sizeof(struct field));
		for (ind=0; ind<N_gridphi_fieldperiod; ind++) {
			Bfield_axis[ind] = malloc(4*sizeof(struct field));
		}
		printf("YOLO\n");
		fieldcentre = solve_magneticaxis(allparams, Bfield_axis, &axis, N_gridphi_fieldperiod, error);
		linalg2x2(axis.tangent, evec, eval, &det, &trace);
		printf("evec=%f\n", evec[0][0]);
		iota_axis = -eval[1]*m0_fieldperiods/(2.0*M_PI);
		printf("iota_axis=%f\n", iota_axis);
		fprintf(Poincare, "%f %f %f\n", 0.0, fieldcentre->loc[0], fieldcentre->loc[1]);
		fprintf(IOTA, "%f %f %f\n", axis.loc[0], axis.loc[1], iota_axis);
	}
	/* loop over different initial values of (R,Z) */
	//axis.loc[0] = 1.0;
	//axis.loc[1] = 0.0;
	for (ind=0;ind<countlines;ind+=1) {
		Xp.loc[0] = axis.loc[0] + rmin[ind]; Xp.loc[1] = axis.loc[1] + zmin[ind];
		Xp.tangent[0][0]=1.0; Xp.tangent[0][1]=0.0; Xp.tangent[1][0]=0.0; Xp.tangent[1][1]=1.0;
		varphi = 0.0;
		printf("initial (R,Z)=(%f,%f)\n", Xp.loc[0], Xp.loc[1]);
		angle_old = 0.0;
		number_turns = 0;
		r_position = rmin[ind]; z_position = zmin[ind];
		for (i_periods=0; i_periods<number_field_periods; i_periods++)  {
			fprintf(Poincare, "%f %f\n", Xp.loc[0], Xp.loc[1]);
			//printf("at varphi=%f, (R,Z)=(%f,%f)\n", varphi, Xp.loc[0], Xp.loc[1]);
			for (i=0; i<N_gridphi_fieldperiod; i++)  {
				if (include_centre == 1) {
					r_position = Xp.loc[0] - fieldcentre[i].loc[0];
					z_position = Xp.loc[1] - fieldcentre[i].loc[1];
				}
				else {
					r_position = Xp.loc[0] - 1.0;
					z_position = Xp.loc[1];
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
				angle[i_periods*N_gridphi_fieldperiod + i] = angle_round_axis;
				phidata[i_periods*N_gridphi_fieldperiod + i] = varphi;
				//printf("%f\t%f\n", angle[i], phidata[i]);
				angle_old = angle_round_axis;
				RK4step(&Xp, varphi, dvarphi, allparams, Bfield_axis[i]);
//void RK4step(struct position *Xp, double varphi, double dvarphi, struct fieldparams allparams, struct field *Bfield_saved) {
				//printstructposition("fieldcentre", fieldcentre+i);
				varphi += dvarphi;
			}
		}
		gsl_fit_linear(phidata, 1, angle, 1, number_field_periods*N_gridphi_fieldperiod, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
		//gsl_fit_linear(phidata, 1, angle, 1, 32768, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
		iota[ind] = c1;
		printf("\niota=%f\tintercept=%f\n", iota[ind], c0);
		fprintf(IOTA, "%f %f %f\n", axis.loc[0] + rmin[ind], axis.loc[1] + zmin[ind], iota[ind]);
	}
	fclose(gridPoincare);
	fclose(Poincare);
	fclose(IOTA);
	free(evec[0]); free(evec[1]);
	free(Xp.tangent[0]); free(Xp.tangent[1]);
	free(evec);
	free(Xp.tangent);
	clock_t final_t = clock();
	printf("Time after obtaining iota profile and Poincare plot is: %f\n", (double) (final_t-initial_t)/CLOCKS_PER_SEC);
}
