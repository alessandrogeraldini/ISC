#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "magfield.h"
#include "linetodata.h"
#include "extract_coils.h"

double *Bfield(double *XXp, double varphi, double ***coils, int *num_coils, int **num_segs)
{
	//clock_t start = clock();
	double BZ=0.0, BX=0.0, BY=0.0, BR = 0.0, Bvarphi = 0.0;
	double Xminxc=0.0, Yminyc=0.0, Zminzc=0.0, Rminrc=0.0, dl=0.0, dxc=0.0, dyc=0.0, dzc=0.0;
	double XX=XXp[0]*cos(varphi), YY=XXp[0]*sin(varphi), ZZ = XXp[1];
	double *magfield=NULL;
	int coil_index=0, coilseg_index=0; 
	//int j, k;

	//n_coils = *num_coils;
	//n_segs = calloc(n_coils,sizeof(int));
	//n_segs = *num_segs;

	//printf("n_segs[1]=%d\n", n_segs[1]);

	//j=3; k = 198;
	//printf("Sanity check: print line of coil file\n");
	//printf("coils[0][%d][%d]=%f\n", j, k, coils[0][j][k]);
	//printf("coils[1][%d][%d]=%f\n", j, k, coils[1][j][k]);
	//printf("coils[2][%d][%d]=%f\n", j, k, coils[2][j][k]);
	//printf("coils[3][%d][%d]=%f\n", j, k, coils[3][j][k]);
	//printf("Sanity check: print another line of coil file\n");
	//printf("coils[0][0][1]=%Le\n", coils[0][0][1]);
	//printf("coils[1][0][1]=%Le\n", coils[1][0][1]);
	//printf("coils[2][0][1]=%Le\n", coils[2][0][1]);
	//printf("coils[3][0][1]=%Le\n", coils[3][0][1]);
	for (coil_index=0;coil_index<*num_coils;coil_index++)
	{
		//printf("num_coils = %d\n", *num_coils);
		//printf("coil_index = %d\n", coil_index);
		for (coilseg_index=0;coilseg_index<*(*num_segs+coil_index)-1;coilseg_index++)
		{
			//printf("num_segs[%d] = %d\n", coil_index, *(*num_segs+coil_index));
			//printf("coilseg_index = %d\n", coilseg_index);
			Xminxc = XX - 0.5*coils[0][coil_index][coilseg_index] - 0.5*coils[0][coil_index][coilseg_index+1];
			Yminyc = YY - 0.5*coils[1][coil_index][coilseg_index] - 0.5*coils[1][coil_index][coilseg_index+1]; 
			Zminzc = ZZ - 0.5*coils[2][coil_index][coilseg_index] - 0.5*coils[2][coil_index][coilseg_index+1]; 
			Rminrc = sqrt(pow(Xminxc,  2.0) + pow(Yminyc, 2.0)  + pow(Zminzc, 2.0)); 
			dxc = coils[0][coil_index][coilseg_index+1] - coils[0][coil_index][coilseg_index];
			dyc = coils[1][coil_index][coilseg_index+1] - coils[1][coil_index][coilseg_index];
			dzc = coils[2][coil_index][coilseg_index+1] - coils[2][coil_index][coilseg_index];
			//printf("dyc = %Le\n", dyc);
			dl = sqrt(pow(dxc, 2.0) + pow(dyc, 2.0) + pow(dzc, 2.0));
			BX += (coils[3][coil_index][coilseg_index]*(dyc*Zminzc - dzc*Yminyc)/pow(Rminrc, 3.0));
			BY += (coils[3][coil_index][coilseg_index]*(dzc*Xminxc - dxc*Zminzc)/pow(Rminrc, 3.0));
			BZ += (coils[3][coil_index][coilseg_index]*(dxc*Yminyc - dyc*Xminxc)/pow(Rminrc, 3.0));
		}
	}
	BX*=pow(10.0,-7.0); BY*=pow(10.0,-7.0); BZ*=pow(10.0,-7.0); // multiply by mu_0/4*PI
	//if (BX> 10.0)
	//{
	////printf("BX = %f, BY = %f, BZ = %f\n", BX, BY, BZ);
	//}
	BR = cos(varphi)*BX + sin(varphi)*BY;
	Bvarphi = -sin(varphi)*BX + cos(varphi)*BY;

	magfield = calloc(3,sizeof(double));
	magfield[0] = BR;
	magfield[1] = BZ;
	magfield[2] = Bvarphi;
	//clock_t int3 = clock();
	//printf("Time after evaluating B field at a point: %f\n", (double) (int3-start)/CLOCKS_PER_SEC);
	return magfield;
}

