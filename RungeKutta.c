//Author: Alessandro Geraldini
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "magfield.h"

double *Ffunc(double Bparg[2], double Bvarphiarg, double Rarg)
{
	double *Freturn;
	Freturn = calloc(2,sizeof(double)); 
	Freturn[0] = Rarg*Bparg[0]/Bvarphiarg;
	Freturn[1] = Rarg*Bparg[1]/Bvarphiarg;
	return Freturn;
}

double *RK4(double *Xp, double varphi, double dvarphi, double ***coils, int *num_coils, int **num_segs)
{
	double *FF, dXp1[2], dXp2[2], dXp3[2], dXp4[2], *Xp1, *Xp2, *Xp3;
	double Bp[2], Bvarphi;	
	double *Bpoint;
	double *Xp4;
	int i=0;
	Xp1 = calloc(2,sizeof(double));
	Xp2 = calloc(2,sizeof(double));
	Xp3 = calloc(2,sizeof(double));
	Xp4 = calloc(2,sizeof(double));
	
	Bpoint = Bfield(Xp, varphi, coils, num_coils, num_segs);
	//printf("Bpoint[0] = %f\n", Bpoint[0]);
	Bp[0] = Bpoint[0]; Bp[1] = Bpoint[1]; Bvarphi = Bpoint[2];
	free(Bpoint);
	//printf("Bp[0] = %f\n", Bp[0]);
	for (i=0;i<2;i++)
	{
		FF = Ffunc(Bp, Bvarphi, Xp[0]);
		dXp1[i] = 0.5*dvarphi*FF[i];
		Xp1[i] = Xp[i] + dXp1[i];
	}
	free(FF);
	//printf("Xp1[0] = %f\n", Xp1[0]);
	Bpoint = Bfield(Xp1, varphi+(dvarphi/2.0), coils, num_coils, num_segs);
	//printf("Bpoint[0] = %f\n", Bpoint[0]);
	Bp[0] = Bpoint[0]; Bp[1] = Bpoint[1]; Bvarphi = Bpoint[2];
	//printf("Bp[0] = %f\n", Bp[0]);
	free(Bpoint);
	for (i=0;i<2;i++)
	{
		FF = Ffunc(Bp, Bvarphi, Xp1[0]);
		dXp2[i] = 0.5*dvarphi*FF[i];
		Xp2[i] = Xp[i] + dXp2[i];
	}
	free(FF);
	//printf("Xp2[0] = %f\n", Xp2[0]);
	Bpoint = Bfield(Xp2, varphi+(dvarphi/2.0), coils, num_coils, num_segs);
	//printf("Bpoint[0] = %f\n", Bpoint[0]);
	Bp[0] = Bpoint[0]; Bp[1] = Bpoint[1]; Bvarphi = Bpoint[2];
	//printf("Bp[0] = %f\n", Bp[0]);
	free(Bpoint);
	for (i=0;i<2;i++)
	{
		FF = Ffunc(Bp, Bvarphi, Xp2[0]);
		dXp3[i] = 0.5*dvarphi*FF[i];
		Xp3[i] = Xp[i] + 2.0*dXp3[i];
	}
	free(FF);
	//printf("Xp3[0] = %f\n", Xp3[0]);
	Bpoint = Bfield(Xp3, varphi+dvarphi, coils, num_coils, num_segs);
	//printf("Bpoint[0] = %f\n", Bpoint[0]);
	Bp[0] = Bpoint[0]; Bp[1] = Bpoint[1]; Bvarphi = Bpoint[2];
	//printf("Bp[0] = %f\n", Bp[0]);
	free(Bpoint);
	for (i=0;i<2;i++)
	{
		FF = Ffunc(Bp, Bvarphi, Xp3[0]);
		dXp4[i] = 0.5*dvarphi*FF[i];
		Xp4[i] = Xp[i] + (1.0/3.0)*(dXp1[i] + 2.0*dXp2[i] + 2.0*dXp3[i] + dXp4[i]); 
	}
	free(FF);
	//printf("Xp4[0] = %f\n", Xp4[0]);
	free(Xp1);
	free(Xp2);
	free(Xp3);
	return Xp4;
}
