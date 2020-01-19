
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
	funcreturn->tangent = malloc(2*sizeof(double*));
	funcreturn->tangent[0] = malloc(2*sizeof(double)); funcreturn->tangent[1] = malloc(2*sizeof(double));
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

struct position *dfunc(struct field *Bparg, struct field *dBparg, struct position *Xparg)
{
	struct position *funcreturn=malloc(sizeof(struct position));
	double *Jacobian[2], Dirac;
	double *p1tangent = &(Xparg->tangent[0][0]), **p2tangent = &p1tangent;
	int row, col;
	funcreturn->tangent = malloc(2*sizeof(double*));
	funcreturn->tangent[0] = malloc(2*sizeof(double)); funcreturn->tangent[1] = malloc(2*sizeof(double));
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
		funcreturn->loc[row] = Xparg->loc[0]*dBparg->value[row]/Bparg->value[2];
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
	structsum.tangent = malloc(2*sizeof(double*));
	structsum.tangent[0] = malloc(2*sizeof(double)); structsum.tangent[1] = malloc(2*sizeof(double));
	for (index1=0;index1<2;index1++) {
		structsum.loc[index1] = num1*struct1->loc[index1] + num2*struct2->loc[index1];
		for (index2=0;index2<2;index2++)	
		structsum.tangent[index1][index2] = num1*struct1->tangent[index1][index2] + num2*struct2->tangent[index1][index2];
	}
	return structsum;
}

struct field addstructsfield(double num1, struct field *struct1, double num2, struct field *struct2) {
	struct field structsum;
	int index1, index2;
	for (index1=0;index1<3;index1++) {
		structsum.value[index1] = num1*struct1->value[index1] + num2*struct2->value[index1];
		for (index2=0;index2<2;index2++)	
		structsum.derivative[index1][index2] = num1*struct1->derivative[index1][index2] + num2*struct2->derivative[index1][index2];
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
	//printf("magnetic field before 1st segment of RK iteratiion: Bfield=%f\n", Bpoint->value[0]);
	dXp1 = func(Bpoint, Xp);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp1);
	//printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
	dXp2 = func(Bpoint, Xp);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp2);
	//printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
	dXp3 = func(Bpoint, Xp);
	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp3);
	//printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpoint = Bfield(Xp->loc, varphi+dvarphi, coils, num_coils, num_segs);
	dXp4 = func(Bpoint, Xp);
	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp1);
	*Xp = addstructs(1.0, Xp,     (1.0/3.0)*dvarphi, dXp2);
	*Xp = addstructs(1.0, Xp,     (1.0/3.0)*dvarphi, dXp3);
	*Xp = addstructs(1.0, Xp,     (1.0/6.0)*dvarphi, dXp4);
	//free(dXp1->tangent[0]); free(dXp1->tangent[1]); free(dXp1->tangent);
	//free(dXp2->tangent[0]); free(dXp2->tangent[1]); free(dXp2->tangent);
	//free(dXp3->tangent[0]); free(dXp3->tangent[1]); free(dXp3->tangent);
	//free(dXp4->tangent[0]); free(dXp4->tangent[1]); free(dXp4->tangent);
	//free(dXp1); free(dXp2); free(dXp3); free(dXp4);
	//printf("after 4th segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	//free(dXp1); free(dXp2); free(dXp3); free(dXp4);
	return;
}


void RK4_wgrad(struct position *Xp, struct position *Xpgrad, double varphi, double dvarphi, double ***coils, int *num_coils, int **num_segs)
{
	int index1, index2;
	double Bp[2], Bvarphi, dBpdX[2][2], dBvarphidX[2];	
	struct field *Bpoint, *Bpointgrad;
	struct position *dXp1, *dXp2, *dXp3, *dXp4, Xpold;
	struct position *dXpgrad1, *dXpgrad2, *dXpgrad3, *dXpgrad4, Xpgradold;
	int i=0;
	Xpold = *Xp;
	Xpgradold = *Xpgrad;
	Bpoint = Bfield(Xp->loc, varphi, coils, num_coils, num_segs);
	Bpointgrad = gradBfield(Xp->loc, varphi, coils, num_coils, num_segs);
	//printf("magnetic field before 1st segment of RK iteratiion: Bfield=%f\n", Bpoint->value[0]);
	printf("YOOOO\n");
	dXp1 = func(Bpoint, Xp);
	dXpgrad1 = dfunc(Bpoint, Bpointgrad, Xp);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp1);
	*Xpgrad = addstructs(1.0, &Xpgradold, 0.5*dvarphi, dXpgrad1);
	printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
	Bpointgrad = gradBfield(Xp->loc, varphi, coils, num_coils, num_segs);
	dXp2 = func(Bpoint, Xp);
	dXpgrad2 = dfunc(Bpoint, Bpointgrad, Xp);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp2);
	*Xpgrad = addstructs(1.0, &Xpgradold, 0.5*dvarphi, dXpgrad2);
	printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
	Bpointgrad = gradBfield(Xp->loc, varphi, coils, num_coils, num_segs);
	dXp3 = func(Bpoint, Xp);
	dXpgrad3 = dfunc(Bpoint, Bpointgrad, Xp);
	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp3);
	*Xpgrad = addstructs(1.0, &Xpgradold, 0.5*dvarphi, dXpgrad3);
	printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpoint = Bfield(Xp->loc, varphi+dvarphi, coils, num_coils, num_segs);
	Bpointgrad = gradBfield(Xp->loc, varphi, coils, num_coils, num_segs);
	dXp4 = func(Bpoint, Xp);
	dXpgrad3 = dfunc(Bpoint, Bpointgrad, Xp);
	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp1);
	*Xp = addstructs(1.0, Xp,     (1.0/3.0)*dvarphi, dXp2);
	*Xp = addstructs(1.0, Xp,     (1.0/3.0)*dvarphi, dXp3);
	*Xp = addstructs(1.0, Xp,     (1.0/6.0)*dvarphi, dXp4);
	*Xpgrad = addstructs(1.0, &Xpgradold, (1.0/6.0)*dvarphi, dXpgrad1);
	*Xpgrad = addstructs(1.0, Xpgrad,     (1.0/3.0)*dvarphi, dXpgrad2);
	*Xpgrad = addstructs(1.0, Xpgrad,     (1.0/3.0)*dvarphi, dXpgrad3);
	*Xpgrad = addstructs(1.0, Xpgrad,     (1.0/6.0)*dvarphi, dXpgrad4);
	printf("after 4th segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	//free(dXp1->tangent[0]); free(dXp1->tangent[1]); free(dXp1->tangent);
	//free(dXp2->tangent[0]); free(dXp2->tangent[1]); free(dXp2->tangent);
	//free(dXp3->tangent[0]); free(dXp3->tangent[1]); free(dXp3->tangent);
	//free(dXp4->tangent[0]); free(dXp4->tangent[1]); free(dXp4->tangent);
	//free(dXp1); free(dXp2); free(dXp3); free(dXp4);
	//free(dXp1); free(dXp2); free(dXp3); free(dXp4);
	return;
}
