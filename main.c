#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "isc.h"
#include "input.h"
#include <gsl/gsl_sf_bessel.h>

int main()
{
	clock_t start = clock();
	//int N_gridphi_field_period=20, field_periods=3;
	int *n_coils=NULL, **n_segs=NULL, qq=1, *qq_segs=NULL, qq_segsq=1;
	int make_Poincare_iota=2, find_axis = 2;
	//struct field *BB;
	struct position *Xp=calloc(1,sizeof(struct position)), *axis;
	struct ext_position *ext_centre;
	//double varphi=0.0; 
	double ***coils=NULL, *width;
	double *r_min=NULL, *iota=NULL;
	double **evec=malloc(2*sizeof(double*)), *eval=malloc(2*sizeof(double)), trace, det, iota_axis;
	evec[0]=malloc(2*sizeof(double)); evec[1]=malloc(2*sizeof(double));
	/* make arrays of coils */
	n_coils = &qq; qq_segs = &qq_segsq; n_segs = &qq_segs; // Initializes all pointers
	coils = coil_grid(n_coils, n_segs);
	clock_t int1 = clock();
	printf("Time after creating arrays of coils: %f\n", (double) (int1-start)/CLOCKS_PER_SEC);
	Xp->tangent = malloc(2*sizeof(double*));
	Xp->tangent[0] = malloc(2*sizeof(double)); Xp->tangent[1] = malloc(2*sizeof(double));
	/* find magnetic axis */
	printf("Do you want to find the axis? (1=yes; 0=no) ");
	scanf("%d", &find_axis);
	if (find_axis == 1) {
		Xp->loc[0] = 1.55; Xp->loc[1]= 0.0;
		Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		axis = findcentre(coils, n_coils, n_segs, Xp);
		linalg2x2(Xp->tangent, evec, eval, &det, &trace);
		printf("evec=%f\n", evec[0][0]);
		iota_axis = eval[1]/(2*M_PI);
		printf("iota_axis=%f\n", iota_axis);
		printf("unity?%f\n", eval[0]);
	}
	/* make Poincare plots (optional) */
	printf("Do you want to produce a Poincare plot and an iota profile? (1=yes; 0=no) ");
	scanf("%d", &make_Poincare_iota);
	if (make_Poincare_iota==1) {
		Xp->loc[0] = 1.52; Xp->loc[1]=0.0;
		Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		//printf("centre is (%f, %f)\n", Xp->loc[0], Xp->loc[1]);
		iotaprofile(r_min, iota, coils, n_coils, n_segs);
		clock_t int3 = clock();
		printf("Time after filling Poincare plot file: %f\n", (double) (int3-start)/CLOCKS_PER_SEC);
		}
	//Xp->loc[0] = 1.3; Xp->loc[1]= -0.515; // island is at(1.299008, -0.515194)
	Xp->loc[0] = 1.299009; Xp->loc[1]= -0.515194; // island is at(1.299008, -0.515194)
	Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
	//island_centre = findisland(coils, n_coils, n_segs, Xp,  12, 3);
	ext_centre = alongcentre(Xp->loc[0], Xp->loc[1], 1, 4, coils, n_coils, n_segs);
	width = islandwidth(ext_centre, 1, 4);
	//BB = Bfield(Xp->loc, varphi, coils, n_coils, n_segs);  //printf("BR(R=%f,Z=%f)=%f\n", Xp->loc[0], Xp->loc[1], BB->value[0]); //clock_t int2 = clock(); //printf("Time after evaluating B field: %f\n", (double) (int2-start)/CLOCKS_PER_SEC);
	
	return 0;
}
