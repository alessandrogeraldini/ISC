
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "isc.h"
#include "input.h"
#include <gsl/gsl_sf_bessel.h>

void iotaprofile(double *rmin, double *iota, double ***coils, int *n_coils, int **n_segs) {
	/* declarations */
	clock_t initial_t = clock();
	//int N_gridphi_field_period=20, field_periods=3; // now in separate file as defined constants
	int i=0, ind, index_centre;
	int N_field_periods=200, n_points = 40, number_turns = 0;
	struct position *Xp=calloc(1,sizeof(struct position)), *fieldcentre, startcentre;
	double varphi=0.0, RRin=1.55, ZZin=0.0, deltar=0.002;
	double dvarphi = 2.0*M_PI/(N_gridphi_field_period*field_periods);
	double angle_round_axis=0.0, angle_round_axis_principal, angle_old, z_position, r_position;
	double angle[N_field_periods*N_gridphi_field_period], phidata[N_field_periods*N_gridphi_field_period];
	double c0=0.0, c1=0.0, cov00=0.0, cov01=0.0, cov11=0.0, sumsq=0.0;
	double evec[2][2], eval[2], det, trace, iota_axis;
	FILE *Poincare=NULL, *IOTA;
	startcentre.loc[0] = RRin; startcentre.loc[1] = ZZin;
	startcentre.tangent[0][0] = 1.0; startcentre.tangent[0][1] = 0.0; startcentre.tangent[1][0] = 0.0; startcentre.tangent[1][1] = 1.0;
	/* prepare for loop over different initial values of (R,Z) */
	fieldcentre = findcentre(coils, n_coils, n_segs, &startcentre);
	rmin = malloc(n_points*sizeof(double));
	iota = malloc(n_points*sizeof(double));
	linalg2x2(startcentre.tangent, evec, eval, &det, &trace);
	printf("evec=%f\n", evec[0][0]);
	iota_axis = eval[1];
	
	if ((Poincare = fopen("Poincare.txt", "w"))==NULL) printf("Could not open Poincare.txt\n");
	if ((IOTA = fopen("iota.txt", "w"))==NULL) printf("Could not open iota.txt\n");
	/* loop over different initial values of (R,Z) */
	for (ind=0;ind<n_points;ind+=1) {
		*(rmin+ind) = ind*deltar;
		Xp->loc[0] = startcentre.loc[0] + ind*deltar; Xp->loc[1] = startcentre.loc[1];
		Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		varphi = 0.0;
		printf("initial (R,Z)=(%f,%f)\n", Xp->loc[0], Xp->loc[1]);
		angle_old = 0.0;
		number_turns = 0;
		for (i=0; i<N_field_periods*N_gridphi_field_period; i++)  {
			index_centre = i%N_gridphi_field_period;
			r_position = Xp->loc[0] - (fieldcentre+index_centre)->loc[0];
			z_position = Xp->loc[1] - (fieldcentre+index_centre)->loc[1];
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
			angle[i] = angle_round_axis;
			phidata[i] = varphi;
			//printf("%f\t%f\n", angle[i], phidata[i]);
			angle_old = angle_round_axis;
			if (index_centre==0)  {
				printf("at varphi=%f, (R,Z)=(%f,%f)\n", varphi, Xp->loc[0], Xp->loc[1]);
				fprintf(Poincare, "%f %f %f\n", varphi, Xp->loc[0], Xp->loc[1]);
			}
			RK4(Xp, varphi, dvarphi, coils, n_coils, n_segs);
			//printstruct("fieldcentre", fieldcentre+i);
			varphi += dvarphi;
		}
		gsl_fit_linear(phidata, 1, angle, 1, N_field_periods*N_gridphi_field_period, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
		iota[ind] = c1;
		printf("\niota=%f\tintercept=%f\n", iota[ind], c0);
		fprintf(IOTA, "%f %f\n", *(rmin+ind), *(iota+ind));
	}
	fclose(Poincare);
	fclose(IOTA);
	clock_t final_t = clock();
	printf("Time after obtaining iota profile and Poincare plot is: %f\n", (double) (final_t-initial_t)/CLOCKS_PER_SEC);
}
