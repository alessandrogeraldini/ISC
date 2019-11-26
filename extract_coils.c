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
