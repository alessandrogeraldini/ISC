#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "isc.h"

double ***coil_grid(int *num_coils, int **num_segs)
{
clock_t start = clock();
int lenline = 1000, line_num, ii=0;
double **X_coil=NULL, **Y_coil=NULL, **Z_coil=NULL, **I_coil=NULL;
double *storevals;
char line[lenline];
int nn, max_num_segs, num_coils_temp=0, coil_index=0, coilseg_index=0, *num_segs_temp;
double ***XXvect_coil= NULL;
XXvect_coil = calloc(4,sizeof(double));
FILE *file_coils;

clock_t int1 = clock();
printf("Time after declarations: %f\n", (double) (int1-start)/CLOCKS_PER_SEC);
file_coils = fopen("coils.c09r00", "r");
if (file_coils == NULL)
{	
	printf("Cannot open coils.c09r00\n");
	exit(1);
}
line_num = 0;
while (fgets(line, lenline, file_coils) != NULL)
{	
	line_num += 1;
	storevals = linetodata(line, &nn);
	if (nn > 4)
	{
		coil_index += 1;
	}
}
rewind(file_coils);
num_coils_temp = coil_index;
num_segs_temp = calloc(coil_index,sizeof(int));
*num_segs = calloc(coil_index,sizeof(int));
max_num_segs = 0;
coil_index = 0;
while (fgets(line, lenline, file_coils) != NULL)
{	
	storevals = linetodata(line, &nn);
	if (nn == 4)
	{
		coilseg_index += 1; 
	}
	else if (nn > 4)
	{
		num_segs_temp[coil_index] = coilseg_index;
		coil_index += 1;
		if (coilseg_index > max_num_segs)
		{
			max_num_segs = coilseg_index;
		}
		coilseg_index = 0;
	}
}
printf("num_coils_temp = %d\n", num_coils_temp);
coil_index = 0;
coilseg_index = 0;
XXvect_coil = calloc(num_coils_temp,sizeof(double));
for (ii=0;ii<4;ii+=1)
{
	XXvect_coil[ii] = calloc(num_coils_temp,sizeof(double));
	for (coil_index = 0; coil_index < num_coils_temp; coil_index ++)
	{
		XXvect_coil[ii][coil_index] = calloc(max_num_segs,sizeof(double));
	}
}
coil_index = 0;
rewind(file_coils);
while (fgets(line, lenline, file_coils) != NULL)
{	
	storevals = linetodata(line, &nn);
	if (nn > 3)
	{
		//printf("hello world, number of coils = %d, coil_index = %d, coilseg_index = %d\n", num_coils, coil_index, coilseg_index);
		XXvect_coil[0][coil_index][coilseg_index] = *storevals;
		XXvect_coil[1][coil_index][coilseg_index] = *(storevals+1);
		XXvect_coil[2][coil_index][coilseg_index] = *(storevals+2);
		XXvect_coil[3][coil_index][coilseg_index] = *(storevals+3);
		//printf("coil_index=%d, coilseg_index = %d\n", coil_index, coilseg_index);
		//printf("storevals[0] = %Le\n", *storevals);
		//printf("storevals[1] = %Le\n", *(storevals+1));
		//printf("storevals[2] = %Le\n", *(storevals+2));
		//printf("storevals[3] = %Le\n", *(storevals+3));
		coilseg_index += 1; 
	}
	if (nn > 4)
	{
		coil_index += 1;
		coilseg_index = 0;
	}
	free(storevals);
}
fclose(file_coils);
*num_coils = num_coils_temp;
*num_segs = num_segs_temp;
return XXvect_coil;
}

struct field *Bfield(double *Xp, double varphi, double ***coils, int *num_coils, int **num_segs)
{
	//clock_t start = clock();
	double BZ=0.0, BX=0.0, BY=0.0, BR=0.0, Bvarphi=0.0, dBXdR=0.0, dBYdR=0.0, dBZdR=0.0, dBXdZ=0.0, dBYdZ=0.0, dBZdZ=0.0, dBRdR=0.0, dBvarphidR=0.0, dBRdZ=0.0, dBvarphidZ=0.0;
	double Xminxc=0.0, Yminyc=0.0, Zminzc=0.0, Rminrc=0.0, dl=0.0, dxc=0.0, dyc=0.0, dzc=0.0;
	double XX=Xp[0]*cos(varphi), YY=Xp[0]*sin(varphi), ZZ = Xp[1];
	struct field *magfield=calloc(1,sizeof(struct field));
	int coil_index=0, coilseg_index=0; 
	int derivatives[2];
	int check = 0;
	double XpdR, YpdR, RpdRminrc, XpdRminxc, YpdRminyc, BXpdR, BYpdR, BZpdR, dR = 0.0001;
	double RpdZminrc, ZpdZminzc, BXpdZ, BYpdZ, BZpdZ, dZ = 0.0001;
	derivatives[0]=1;
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
			Rminrc = sqrt(pow(Xminxc,  2.0) + pow(Yminyc, 2.0) + pow(Zminzc, 2.0)); 

			dxc = coils[0][coil_index][coilseg_index+1] - coils[0][coil_index][coilseg_index];
			dyc = coils[1][coil_index][coilseg_index+1] - coils[1][coil_index][coilseg_index];
			dzc = coils[2][coil_index][coilseg_index+1] - coils[2][coil_index][coilseg_index];
			dl = sqrt(pow(dxc, 2.0) + pow(dyc, 2.0) + pow(dzc, 2.0));

			BX += (coils[3][coil_index][coilseg_index]*(dyc*Zminzc - dzc*Yminyc)/pow(Rminrc, 3.0));
			BY += (coils[3][coil_index][coilseg_index]*(dzc*Xminxc - dxc*Zminzc)/pow(Rminrc, 3.0));
			BZ += (coils[3][coil_index][coilseg_index]*(dxc*Yminyc - dyc*Xminxc)/pow(Rminrc, 3.0));
			
			if (derivatives[0] == 1)
			{
				dBXdR += (coils[3][coil_index][coilseg_index]*(dyc*0.0         - dzc*sin(varphi))/pow(Rminrc, 3.0));
				dBYdR += (coils[3][coil_index][coilseg_index]*(dzc*cos(varphi) - dxc*0.0        )/pow(Rminrc, 3.0));
				dBZdR += (coils[3][coil_index][coilseg_index]*(dxc*sin(varphi) - dyc*cos(varphi))/pow(Rminrc, 3.0));

				dBXdR -= 3*(coils[3][coil_index][coilseg_index]*(Xminxc*cos(varphi) + Yminyc*sin(varphi))*(dyc*Zminzc - dzc*Yminyc)/pow(Rminrc, 5.0));
				dBYdR -= 3*(coils[3][coil_index][coilseg_index]*(Xminxc*cos(varphi) + Yminyc*sin(varphi))*(dzc*Xminxc - dxc*Zminzc)/pow(Rminrc, 5.0));
				dBZdR -= 3*(coils[3][coil_index][coilseg_index]*(Xminxc*cos(varphi) + Yminyc*sin(varphi))*(dxc*Yminyc - dyc*Xminxc)/pow(Rminrc, 5.0));

				dBXdZ += (coils[3][coil_index][coilseg_index]*( dyc)/pow(Rminrc, 3.0));
				dBYdZ += (coils[3][coil_index][coilseg_index]*(-dxc)/pow(Rminrc, 3.0));
				//dBZdZ += (coils[3][coil_index][coilseg_index]*( 0.0)/pow(Rminrc, 3.0));

				dBXdZ -= 3*(coils[3][coil_index][coilseg_index]*(Zminzc)*(dyc*Zminzc - dzc*Yminyc)/pow(Rminrc, 5.0));
				dBYdZ -= 3*(coils[3][coil_index][coilseg_index]*(Zminzc)*(dzc*Xminxc - dxc*Zminzc)/pow(Rminrc, 5.0));
				dBZdZ -= 3*(coils[3][coil_index][coilseg_index]*(Zminzc)*(dxc*Yminyc - dyc*Xminxc)/pow(Rminrc, 5.0));
			}
		}
	}

	//magfield = calloc(9,sizeof(double)); // first triplet B; if derivatives[0] ==1: second triplet dB/dR, third triplet dB/dZ 
	BX*=pow(10.0,-7.0); BY*=pow(10.0,-7.0); BZ*=pow(10.0,-7.0); // multiply by mu_0/4*PI
	BR = cos(varphi)*BX + sin(varphi)*BY; Bvarphi = -sin(varphi)*BX + cos(varphi)*BY; // BZ = BZ;
	magfield->value[0] = BR; magfield->value[1] = BZ; magfield->value[2] = Bvarphi;
	if (derivatives[0]==1)
	{
		//printf("dBXdR=%f\n", dBXdR);
		dBXdR*=pow(10.0,-7.0); dBYdR*=pow(10.0,-7.0); dBZdR*=pow(10.0,-7.0); // multiply by mu_0/4*PI
		dBXdZ*=pow(10.0,-7.0); dBYdZ*=pow(10.0,-7.0); dBZdZ*=pow(10.0,-7.0); // multiply by mu_0/4*PI
		dBRdR = dBXdR*cos(varphi) + dBYdR*sin(varphi); dBvarphidR = -dBXdR*sin(varphi) + dBYdR*cos(varphi);
		dBRdZ = dBXdZ*cos(varphi) + dBYdZ*sin(varphi); dBvarphidZ = -dBXdZ*sin(varphi) + dBYdZ*cos(varphi);
		//printf("dBXdR=%f attempt 2\n", dBXdR);
		magfield->derivative[0][0] = dBRdR; 
		magfield->derivative[1][0] = dBZdR; 
		magfield->derivative[2][0] = dBvarphidR;
		magfield->derivative[0][1] = dBRdZ; 
		magfield->derivative[1][1] = dBZdZ; 
		magfield->derivative[2][1] = dBvarphidZ;
	}

	if (check==1) {
		BX = 0.0; BY = 0.0; BZ = 0.0;
		BXpdR = 0.0; BYpdR = 0.0, BZpdR = 0.0;
		BXpdZ = 0.0; BYpdZ = 0.0, BZpdZ = 0.0;
		dBXdR = 0.0; dBYdR = 0.0; dBZdR = 0.0;
		dBXdZ = 0.0; dBYdZ = 0.0; dBZdZ = 0.0;
		BXpdR = 0.0; BYpdR = 0.0; BZpdR = 0.0;
		for (coil_index=0;coil_index<*num_coils;coil_index++) {
			for (coilseg_index=0;coilseg_index<*(*num_segs+coil_index)-1;coilseg_index++) {
				Xminxc = XX - 0.5*coils[0][coil_index][coilseg_index] - 0.5*coils[0][coil_index][coilseg_index+1];
				Yminyc = YY - 0.5*coils[1][coil_index][coilseg_index] - 0.5*coils[1][coil_index][coilseg_index+1]; 
				Zminzc = ZZ - 0.5*coils[2][coil_index][coilseg_index] - 0.5*coils[2][coil_index][coilseg_index+1]; 
				Rminrc = sqrt(pow(Xminxc,  2.0) + pow(Yminyc, 2.0) + pow(Zminzc, 2.0)); 
				
				dxc = coils[0][coil_index][coilseg_index+1] - coils[0][coil_index][coilseg_index];
				dyc = coils[1][coil_index][coilseg_index+1] - coils[1][coil_index][coilseg_index];
				dzc = coils[2][coil_index][coilseg_index+1] - coils[2][coil_index][coilseg_index];
				dl = sqrt(pow(dxc, 2.0) + pow(dyc, 2.0) + pow(dzc, 2.0));

				BX += (coils[3][coil_index][coilseg_index]*(dyc*Zminzc - dzc*Yminyc)/pow(Rminrc, 3.0));
				BY += (coils[3][coil_index][coilseg_index]*(dzc*Xminxc - dxc*Zminzc)/pow(Rminrc, 3.0));
				BZ += (coils[3][coil_index][coilseg_index]*(dxc*Yminyc - dyc*Xminxc)/pow(Rminrc, 3.0));

				XpdR = (Xp[0]+dR)*cos(varphi); YpdR = (Xp[0]+dR)*sin(varphi); 
				XpdRminxc = XpdR - 0.5*coils[0][coil_index][coilseg_index] - 0.5*coils[0][coil_index][coilseg_index+1];
				YpdRminyc = YpdR - 0.5*coils[1][coil_index][coilseg_index] - 0.5*coils[1][coil_index][coilseg_index+1]; 
				ZpdZminzc = Xp[1] + dZ - 0.5*coils[2][coil_index][coilseg_index] - 0.5*coils[2][coil_index][coilseg_index+1]; 
				RpdRminrc = sqrt(pow(XpdRminxc,  2.0) + pow(YpdRminyc, 2.0) + pow(Zminzc, 2.0)); 
				RpdZminrc = sqrt(pow(Xminxc,  2.0) + pow(Yminyc, 2.0) + pow(ZpdZminzc, 2.0)); 

				BXpdR += (coils[3][coil_index][coilseg_index]*(dyc*Zminzc - dzc*YpdRminyc)/pow(RpdRminrc, 3.0));
				BYpdR += (coils[3][coil_index][coilseg_index]*(dzc*XpdRminxc - dxc*Zminzc)/pow(RpdRminrc, 3.0));
				BZpdR += (coils[3][coil_index][coilseg_index]*(dxc*YpdRminyc - dyc*XpdRminxc)/pow(RpdRminrc, 3.0));
				BXpdZ += (coils[3][coil_index][coilseg_index]*(dyc*ZpdZminzc - dzc*Yminyc)/pow(RpdZminrc, 3.0));
				BYpdZ += (coils[3][coil_index][coilseg_index]*(dzc*Xminxc - dxc*ZpdZminzc)/pow(RpdZminrc, 3.0));
				BZpdZ += (coils[3][coil_index][coilseg_index]*(dxc*Yminyc - dyc*Xminxc)/pow(RpdZminrc, 3.0));
				}
			}
		dBXdR = (BXpdR - BX)/dR;
		dBYdR = (BYpdR - BY)/dR;
		dBZdR = (BZpdR - BZ)/dR;
		dBXdZ = (BXpdZ - BX)/dZ;
		dBYdZ = (BYpdZ - BY)/dZ;
		dBZdZ = (BZpdZ - BZ)/dZ;
		dBXdR*=pow(10.0,-7.0); dBYdR*=pow(10.0,-7.0); dBZdR*=pow(10.0,-7.0); // multiply by mu_0/4*PI
		dBXdZ*=pow(10.0,-7.0); dBYdZ*=pow(10.0,-7.0); dBZdZ*=pow(10.0,-7.0); // multiply by mu_0/4*PI
		dBRdR = dBXdR*cos(varphi) + dBYdR*sin(varphi); dBvarphidR = -dBXdR*sin(varphi) + dBYdR*cos(varphi);
		dBRdZ = dBXdZ*cos(varphi) + dBYdZ*sin(varphi); dBvarphidZ = -dBXdZ*sin(varphi) + dBYdZ*cos(varphi);

		printf("%f %f\n", magfield->derivative[0][0], dBRdR);
		printf("%f %f\n", magfield->derivative[1][0], dBZdR);
		printf("%f %f\n", magfield->derivative[2][0], dBvarphidR);
		printf("%f %f\n", magfield->derivative[0][1], dBRdZ);
		printf("%f %f\n", magfield->derivative[1][1], dBZdZ);
		printf("%f %f\n", magfield->derivative[2][1], dBvarphidZ);
	}
	//clock_t int3 = clock();
	//printf("Time after evaluating B field at a point: %f\n", (double) (int3-start)/CLOCKS_PER_SEC);
	return magfield;
}

