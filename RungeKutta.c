/*Author: Alessandro Geraldini */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "isc.h"

struct position *func(struct field *Bparg, struct position *Xparg)
{
	struct position *funcreturn=malloc(sizeof(struct position));
	double *Jacobian[2], Dirac;
	double *p1tangent = &(Xparg->tangent[0][0]), **p2tangent = &p1tangent;
	int row, col;
	for (row=0;row<2;row++)
	{	
		Jacobian[row] = malloc(2*sizeof(double));
		for (col=0;col<2;col++)
		{
			if (col==0)	Dirac=1.0;
			else		Dirac=0.0;
			Jacobian[row][col] = Dirac*Bparg->value[row]/Bparg->value[2] 
			+ Xparg->loc[0]*Bparg->derivative[row][col]/Bparg->value[2] 
			- Xparg->loc[0]*Bparg->value[row]*Bparg->derivative[2][col]/pow(Bparg->value[2], 2.0); 
			//printf("Bparg->derivative[0][0] = %f\n", Bparg->derivative[0][0]);
		}
	}
	//funcreturn->tangent = multiply2x2(Jacobian, p2tangent, 2);
	for (row=0;row<2;row++)
	{
		funcreturn->loc[row] = Xparg->loc[0]*Bparg->value[row]/Bparg->value[2];
		for (col=0;col<2;col++)
		{	
			funcreturn->tangent[row][col] = Jacobian[row][0]*Xparg->tangent[0][col] + Jacobian[row][1]*Xparg->tangent[1][col];
		}
	}
	return funcreturn;
}

//struct position *multiplystruct(double num, struct position *structin)
//{
//	struct position *structout;
//	int index1, index2;
//	for (index1=0;index1<2;index1++)
//	{
//		structout->loc[index1] = num*structin->loc[index1];
//		for (index2=0;index2<2;index2++)	
//		structout->tangent[index1][index2] = num*structin->tangent[index1][index2];
//	}
//	return *structout;
//}

struct position addstructs(double num1, struct position *struct1, double num2, struct position *struct2) {
	struct position structsum;
	int index1, index2;
	for (index1=0;index1<2;index1++) {
		structsum.loc[index1] = num1*struct1->loc[index1] + num2*struct2->loc[index1];
		for (index2=0;index2<2;index2++)	
		structsum.tangent[index1][index2] = num1*struct1->tangent[index1][index2] + num2*struct2->tangent[index1][index2];
	}
	return structsum;
}

//struct ext_position addextstructs(double num1, struct ext_position *struct1, double num2, struct ext_position *struct2)
//{
//	struct ext_position structsum;
//	int index1, index2;
//	for (index1=0;index1<2;index1++) {
//		structsum.loc[index1] = num1*struct1->loc[index1] + num2*struct2->loc[index1];
//		for (index2=0;index2<2;index2++) {
//		structsum.comptangent[index1][index2]     = num1*struct1->comptangent[index1][index2] 
//							  + num2*struct2->comptangent[index1][index2];
//		structsum.parttangent[index1][index2]     = num1*struct1->parttangent[index1][index2] 
//							  + num2*struct2->parttangent[index1][index2];
//		structsum.fulltangent[index1][index2]     = num1*struct1->fulltangent[index1][index2] 
//							  + num2*struct2->fulltangent[index1][index2];
//		structsum.symmfulltangent[index1][index2] = num1*struct1->symmfulltangent[index1][index2] 
//							  + num2*struct2->symmfulltangent[index1][index2];
//		}
//	}
//	return structsum;
//}

void RK4(struct position *Xp, double varphi, double dvarphi, double ***coils, int *num_coils, int **num_segs)
{
	int index1, index2;
	double Bp[2], Bvarphi, dBpdX[2][2], dBvarphidX[2];	
	struct field *Bpoint;
	struct position *dXp1, *dXp2, *dXp3, *dXp4, Xpold;
	int i=0;
	Xpold = *Xp;
	Bpoint = Bfield(Xp->loc, varphi, coils, num_coils, num_segs);
	dXp1 = func(Bpoint, Xp);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp1);
	//printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpoint = Bfield(Xp->loc, varphi+(dvarphi/2.0), coils, num_coils, num_segs);
	dXp2 = func(Bpoint, Xp);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp2);
	//printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpoint = Bfield(Xp->loc, varphi+(dvarphi/2.0), coils, num_coils, num_segs);
	dXp3 = func(Bpoint, Xp);
	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp3);
	//printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpoint = Bfield(Xp->loc, varphi+dvarphi, coils, num_coils, num_segs);
	dXp4 = func(Bpoint, Xp);
	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp1);
	*Xp = addstructs(1.0, Xp,     (1.0/3.0)*dvarphi, dXp2);
	*Xp = addstructs(1.0, Xp,     (1.0/3.0)*dvarphi, dXp3);
	*Xp = addstructs(1.0, Xp,     (1.0/6.0)*dvarphi, dXp4);
	//printf("after 4th segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	//free(dXp1); free(dXp2); free(dXp3); free(dXp4);
	return;
}
