#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "magfield.h"
#include "extract_coils.h"
#include "linetodata.h"
#include "RungeKutta.h"

int main()
{
	clock_t start = clock();
	int N_gridphi_field_period=20, field_periods=3, i=0, *n_coils=NULL, **n_segs=NULL, qq = 1, *qq_segs=NULL, qq_segsq = 1, ind;
	double *BB=NULL;
	double varphi=0.0, RRin=1.45, ZZin=0.0;
	double *Xp=NULL, *Xpnew = NULL;
	double dvarphi = 2.0*M_PI/(N_gridphi_field_period*field_periods);
	double ***coils=NULL;
	FILE *Poincare=NULL;

	Xp = malloc(2*sizeof(double));
	n_coils = &qq; qq_segs = &qq_segsq; n_segs = &qq_segs; // Initializes all pointers
	coils = coil_grid(n_coils, n_segs);
	clock_t int1 = clock();
	printf("Time after creating arrays of coils: %f\n", (double) (int1-start)/CLOCKS_PER_SEC);
	Xp[0]=1.55;Xp[1]=0.0;
	BB = Bfield(Xp, varphi, coils, n_coils, n_segs); //free(BB);  //printf("BR = %f, BZ= %f, Bvarphi=%f\n", BB[0], BB[1], BB[2]);
	printf("B(R=%f,Z=%f)=%f\n", Xp[0], Xp[1], *BB);
	clock_t int2 = clock(); //printf("Time after evaluating B field: %f\n", (double) (int2-start)/CLOCKS_PER_SEC);

	if ((Poincare = fopen("Poincare.txt", "w"))==NULL)
	{
		printf("Could not open Poincare.txt\n");
	}
	for (ind=0;ind<10;ind+=1)
	{
		Xp[0]=RRin+ind*0.02; Xp[1]=ZZin;
		varphi = 0.0;
		printf("RR=%f\n", Xp[0]);
		for (i=0; i<100*N_gridphi_field_period; i++)
		{
			if (i%N_gridphi_field_period==0)
			{
				printf("varphi = %f\n", varphi);
				fprintf(Poincare, "%f %f %f\n", varphi, Xp[0], Xp[1]);
			}
			Xpnew = RK4(Xp, varphi, dvarphi, coils, n_coils, n_segs);
			free(Xp);
			Xp = Xpnew;
			//free(Xpnew);
			varphi += dvarphi;
		}
	}
	fclose(Poincare);
	clock_t int3 = clock();
	printf("Time after filling Poincare plot file: %f\n", (double) (int3-start)/CLOCKS_PER_SEC);
	return 0;
}
