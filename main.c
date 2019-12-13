//Author: Alessandro Geraldini
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "isc.h"
#include <gsl/gsl_sf_bessel.h>

int main()
{
	clock_t start = clock();
	int N_gridphi_per_field_period=30, m0_symmetry=5, N_gridphi_tor=N_gridphi_per_field_period*m0_symmetry;
	int *n_coils=NULL, **n_segs=NULL, qq=1, *qq_segs=NULL, qq_segsq=1, pol_mode = 5, tor_mode = 5;
	int make_Poincare_iota=2, find_axis = 2, find_islands=2, n_points=20, N_field_periods_map = 2000;
	//struct field *BB;
	struct position *Xp=calloc(1,sizeof(struct position)), *axis, *island_centre;
	struct ext_position *ext_centre;
	double varphi=0.0; 
	struct field *Bpoint;
	double ***coils=NULL, *width;
	double r_interval = 0.002, *iota=NULL, *minor_radius=NULL;
	double **evec=malloc(2*sizeof(double*)), *eval=malloc(2*sizeof(double)), trace, det, iota_axis;
	evec[0]=malloc(2*sizeof(double)); evec[1]=malloc(2*sizeof(double));
	/* make arrays of coils */
	printf("field periods = %d\n", m0_symmetry);
	printf("N points per field period = %d\n", N_gridphi_per_field_period);
	n_coils = &qq; qq_segs = &qq_segsq; n_segs = &qq_segs; // Initializes all pointers
	coils = coil_grid(n_coils, n_segs);
	clock_t int1 = clock();
	printf("Time after creating arrays of coils: %f\n", (double) (int1-start)/CLOCKS_PER_SEC);
	Xp->loc[0] = 1.002; Xp->loc[1]= 0.0;
	Xp->tangent = malloc(2*sizeof(double*));
	Xp->tangent[0] = malloc(2*sizeof(double)); Xp->tangent[1] = malloc(2*sizeof(double));
	/* find magnetic axis */
	//Bpoint = Bfield(Xp->loc, varphi, coils, n_coils, n_segs);
	//printf("B=(%f, %f, %f)\n", Bpoint->value[0], Bpoint->value[1], Bpoint->value[2]);
	printf("Do you want to find the axis? (1=yes; 0=no) ");
	scanf("%d", &find_axis);
	printf("Do you want to produce a Poincare plot and an iota profile? (1=yes; 0=no) ");
	scanf("%d", &make_Poincare_iota);
	printf("Do you want to find magnetic islands? (1=yes; 0=no) ");
	scanf("%d", &find_islands);
	if (find_axis == 1) {
		//Xp->loc[0] = 1.55; Xp->loc[1]= 0.0;
		Xp->loc[0] = 1.1; Xp->loc[1]= 0.0;
		Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		axis = findcentre(coils, n_coils, n_segs, Xp, N_gridphi_tor);
		linalg2x2(Xp->tangent, evec, eval, &det, &trace);
		printf("evec=%f\n", evec[0][0]);
		iota_axis = eval[1]/(2*M_PI);
		printf("iota_axis=%f\n", iota_axis);
		printf("unity?%f\n", eval[0]);
	}
	/* make Poincare plots (optional) */
	if (make_Poincare_iota==1) {
		//Xp->loc[0] = 1.09; Xp->loc[1]=0.0;
		Xp->loc[0] = 1.009; Xp->loc[1]=0.0;
		Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		//printf("centre is (%f, %f)\n", Xp->loc[0], Xp->loc[1]);
		printf("yolo\n");
		iotaprofile(Xp->loc[0], r_interval, n_points, m0_symmetry, N_gridphi_per_field_period, minor_radius, iota, coils, n_coils, n_segs);
		clock_t int3 = clock();
		printf("Time after filling Poincare plot file: %f\n", (double) (int3-start)/CLOCKS_PER_SEC);
	}

	if (find_islands == 1) {
		tor_mode = 6;
		pol_mode = 3;
		//Xp->loc[0] = 1.3; Xp->loc[1]= -0.515; // island is at(1.299008, -0.515194)
		//Xp->loc[0] = 1.299009; Xp->loc[1]= -0.515194; // NCSX: island is at(1.299008, -0.515194)
		//Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		//island_centre = findisland(coils, n_coils, n_segs, Xp,  3, 12);
		//ext_centre = alongcentre(island_centre[0].loc[0], island_centre[0].loc[1], 3, 12, coils, n_coils, n_segs);
		//width = islandwidth(ext_centre, 3, 12);

		//Xp->loc[0] = 0.945; Xp->loc[1]= 0.0; // Dommaschk (5,2) amp 1.73: island is at(1.299008, -0.515194)
		//Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		//island_centre = findisland(coils, n_coils, n_segs, Xp,  5, 7);
		//ext_centre = alongcentre(island_centre[0].loc[0], island_centre[0].loc[1], 5, 7, coils, n_coils, n_segs);
		//width = islandwidth(ext_centre, 5, 7);

		//Xp->loc[0] = 1.091; Xp->loc[1]= 0.0; // Dommaschk (5,2) amp 0.00001 
		//Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		//island_centre = findisland(coils, n_coils, n_segs, Xp, m0_symmetry, N_gridphi_per_field_period, 5, 2);
		//ext_centre = alongcentre(island_centre[0].loc[0], island_centre[0].loc[1], m0_symmetry, N_gridphi_per_field_period, 5, 2, coils, n_coils, n_segs);
		//width = islandwidth(ext_centre, m0_symmetry, N_gridphi_per_field_period, 5, 2);
		Xp->loc[0] = 1.0097; Xp->loc[1]= 0.0; // Dommaschk (5,2) amp 0.00001 
		Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		island_centre = findisland(coils, n_coils, n_segs, Xp, m0_symmetry, N_gridphi_per_field_period, tor_mode, pol_mode);
		ext_centre = alongcentre(island_centre[0].loc[0], island_centre[0].loc[1], m0_symmetry, N_gridphi_per_field_period, tor_mode, pol_mode, coils, n_coils, n_segs);
		width = islandwidth(ext_centre, m0_symmetry, N_gridphi_per_field_period, tor_mode, pol_mode);
	}

	//BB = Bfield(Xp->loc, varphi, coils, n_coils, n_segs);  //printf("BR(R=%f,Z=%f)=%f\n", Xp->loc[0], Xp->loc[1], BB->value[0]); //clock_t int2 = clock(); //printf("Time after evaluating B field: %f\n", (double) (int2-start)/CLOCKS_PER_SEC);
	
	return 0;
}
