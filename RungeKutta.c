/*Author: Alessandro Geraldini 
  Description: this file contains functions that step forward in toroidal angle
               the functions named derivative contain the derivatives
               the functions named RK4 perform the step forward */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "isc.h"

struct position *derivative_center(struct field *Bparg, struct position *Xparg)
{
	struct position *funcreturn=malloc(sizeof(struct position));
	double Jacobian[2][2], Dirac;
	int row, col;
	//printf("Bp[0]=%f\n", Bparg->value[0]);
	//printf("Bp[1]=%f\n", Bparg->value[1]);
	//printf("Bp[2]=%f\n", Bparg->value[2]);
	funcreturn->tangent = malloc(2*sizeof(double*));
	funcreturn->tangent[0] = malloc(2*sizeof(double)); funcreturn->tangent[1] = malloc(2*sizeof(double));
	for (row=0;row<2;row++)
	{	
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
	//printmat("Jacobian", Jacobian, 2, 2);
	//funcreturn->tangent = multiply2x2(Jacobian, p2tangent, 2);
	for (row=0;row<2;row++)
	{
		funcreturn->loc[row] = Xparg->loc[0]*Bparg->value[row]/Bparg->value[2];
		for (col=0;col<2;col++)
		{	
			funcreturn->tangent[row][col] = Jacobian[row][0]*Xparg->tangent[0][col] + Jacobian[row][1]*Xparg->tangent[1][col];
		}
	}
	//printf("funcereturn=%f\n", funcreturn->loc[0]);
	return funcreturn;
}

struct position *derivative_sperp(struct field *Bparg, struct position *Xparg, struct position *sperp)
{
	struct position *funcreturn=malloc(sizeof(struct position));
	double Jacobian[2][2], Dirac;
	int row, col;
	funcreturn->tangent = malloc(2*sizeof(double*));
	funcreturn->tangent[0] = malloc(2*sizeof(double)); funcreturn->tangent[1] = malloc(2*sizeof(double));
	for (row=0;row<2;row++)
	{	
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
	//printmat("Jacobian", Jacobian, 2, 2);
	//funcreturn->tangent = multiply2x2(Jacobian, p2tangent, 2);
	for (row=0;row<2;row++)
	{
		funcreturn->loc[row] = Jacobian[row][0]*sperp->loc[0] + Jacobian[row][1]*sperp->loc[1];
		for (col=0;col<2;col++)
		{	
			//funcreturn->tangent[row][col] = Jacobian[row][0]*sperp->tangent[0][col] + Jacobian[row][1]*sperp->tangent[1][col];
		}
	}
	return funcreturn;
}

struct position *derivative_lambdacirc_mutangent(struct field *Bparg, struct position *Xparg, struct position *adjvariable) {
	struct position *adjfuncreturn=malloc(sizeof(struct position));
	double Rhat[2];
	int row, col;
	Rhat[0] =1.0; Rhat[1] =0.0;
	adjfuncreturn->tangent = malloc(2*sizeof(double*));
	adjfuncreturn->tangent[0] = malloc(2*sizeof(double)); adjfuncreturn->tangent[1] = malloc(2*sizeof(double));
	for (row=0;row<2;row++)
	{
		adjfuncreturn->loc[row] = 0.0;
		for (col=0;col<2;col++)
		{	
			adjfuncreturn->loc[row] -=  ( adjvariable->loc[col] / Bparg->value[2] ) * ( Rhat[row]*Bparg->value[col] + Xparg->loc[0] * ( Bparg->derivative[col][row] - Bparg->derivative[2][row]*Bparg->value[col]/Bparg->value[2] ) );
			adjfuncreturn->tangent[row][col] = - ( Rhat[row]*Bparg->value[0]*adjvariable->tangent[0][col]/Bparg->value[2] + Rhat[row]*Bparg->value[1]*adjvariable->tangent[1][col]/Bparg->value[2] 
			+ (Xparg->loc[0]/Bparg->value[2]) * (Bparg->derivative[0][row]*adjvariable->tangent[0][col] + Bparg->derivative[1][row]*adjvariable->tangent[1][col] 
			 - Bparg->derivative[2][row]*Bparg->value[0]*adjvariable->tangent[0][col]/Bparg->value[2] - Bparg->derivative[2][row]*Bparg->value[1]*adjvariable->tangent[1][col]/Bparg->value[2] ) ) ;
		}
	}
	return adjfuncreturn;
}

struct position *derivative_lambdatangent(struct field *Bp, struct position *Xp, struct position *lambdaQ, struct position *sperp, struct position *mu) {
	struct position *adjfuncreturn=malloc(sizeof(struct position));
	double Rhat[2];
	double object[2][2];// Hessian[2][2][2];
	int row, col, third;
	Rhat[0] =1.0; Rhat[1] =0.0;
	//for (row=0;row<2;row++) {
	//	for (col=0;col<2;col++) {
	//		for (third=0;third<2;third++) {
	//			Hessian[row][col][third] = Rhat[col]*(Bp->derivative[row][third]/Bp->value[2] - Rhat[col]*Bp->value[row] * Bp->derivative[2][third]/ pow(Bp->value[2], 2.0) ) + Xp->loc[0]*(Bp->twoderivative[row][col][third]/(2.0*Bp->value[2]) - Bp->value[row]*Bp->twoderivative[2][third][col]/(2.0*pow(Bp->value[2], 2.0)) + Bp->value[row]*Bp->derivative[2][col]*Bp->derivative[2][third]/pow(Bp->value[2], 3.0) - Bp->derivative[row][col]*Bp->derivative[2][third]/pow(Bp->value[2], 2.0) )
	//		 	+ Rhat[third]*(Bp->derivative[row][col]/Bp->value[2] - Rhat[third]*Bp->value[row] * Bp->derivative[2][col]/ pow(Bp->value[2], 2.0) ) + Xp->loc[0]*(Bp->twoderivative[row][third][col]/(2.0*Bp->value[2]) - Bp->value[row]*Bp->twoderivative[2][col][third]/(2.0*pow(Bp->value[2], 2.0)) + Bp->value[row]*Bp->derivative[2][third]*Bp->derivative[2][col]/pow(Bp->value[2], 3.0) - Bp->derivative[row][third]*Bp->derivative[2][col]/pow(Bp->value[2], 2.0) );
	//		}
	//	}
	//}
	adjfuncreturn->tangent = malloc(2*sizeof(double*));
	adjfuncreturn->tangent[0] = malloc(2*sizeof(double)); adjfuncreturn->tangent[1] = malloc(2*sizeof(double));
	for (row=0;row<2;row++)
	{
		adjfuncreturn->loc[row] = 0.0;
		for (col=0;col<2;col++) {	
			object[col][row] = 0.0;
			for (third=0;third<2;third++) {
				//object[col][row] +=  Hessian[third][col][row]*mu->loc[third];			
				// I suspect the remaining error is in expression below containing the Hessian
				// Hessian seems fine
				object[col][row] 
				+= ( Rhat[col]*( Bp->derivative[third][row]/Bp->value[2] - Bp->value[third]*Bp->derivative[2][row]/ pow(Bp->value[2], 2.0) ) + Xp->loc[0]*( Bp->twoderivative[third][col][row]/(2.0*Bp->value[2]) - Bp->value[third]*Bp->twoderivative[2][row][col]/(2.0*pow(Bp->value[2], 2.0)) + Bp->value[third]*Bp->derivative[2][col]*Bp->derivative[2][row]/pow(Bp->value[2], 3.0) - Bp->derivative[third][col]*Bp->derivative[2][row]/pow(Bp->value[2], 2.0) )
			 	+    Rhat[row]*( Bp->derivative[third][col]/Bp->value[2] - Bp->value[third]*Bp->derivative[2][col]/ pow(Bp->value[2], 2.0) ) + Xp->loc[0]*( Bp->twoderivative[third][row][col]/(2.0*Bp->value[2]) - Bp->value[third]*Bp->twoderivative[2][col][row]/(2.0*pow(Bp->value[2], 2.0)) + Bp->value[third]*Bp->derivative[2][row]*Bp->derivative[2][col]/pow(Bp->value[2], 3.0) - Bp->derivative[third][row]*Bp->derivative[2][col]/pow(Bp->value[2], 2.0) ) ) * mu->loc[third];
				// I suspect the remaining error is in expression above containing the Hessian
				// Hessian seems fine
			}
			adjfuncreturn->loc[row] -=  ( ( lambdaQ->loc[col] / Bp->value[2] ) * ( Rhat[row]*Bp->value[col] + Xp->loc[0] * ( Bp->derivative[col][row] - Bp->derivative[2][row]*Bp->value[col]/Bp->value[2] ) ) + object[col][row]*sperp->loc[col] );
			adjfuncreturn->tangent[row][col] = - ( Rhat[row]*Bp->value[0]*lambdaQ->tangent[0][col]/Bp->value[2] + Rhat[row]*Bp->value[1]*lambdaQ->tangent[1][col]/Bp->value[2] 
			+ (Xp->loc[0]/Bp->value[2]) * (Bp->derivative[0][row]*lambdaQ->tangent[0][col] + Bp->derivative[1][row]*lambdaQ->tangent[1][col] 
			 - Bp->derivative[2][row]*Bp->value[0]*lambdaQ->tangent[0][col]/Bp->value[2] - Bp->derivative[2][row]*Bp->value[1]*lambdaQ->tangent[1][col]/Bp->value[2] ) ) ;
		}
	}
	return adjfuncreturn;
}

double *derivative_gradcirc(struct field *Bparg, struct field *gradBparg, struct position *Xparg, struct position *lambdaarg, int num_params)
{
	double *adjfunc = NULL;
	int row, col;
	adjfunc = calloc(num_params,sizeof(double)); 
	for (row=0;row<num_params;row++) {
		for (col=0;col<2;col++) {
			adjfunc[row] += ( lambdaarg->loc[col] * Xparg->loc[0] * ( - gradBparg[row].value[col] + gradBparg[row].value[2]*Bparg->value[col] / Bparg->value[2] ) / Bparg->value[2] );
		}
	}
	return adjfunc;
}

//double ***derivative_adjshape(struct field *Bparg, struct field ***shapeBparg, struct position *Xparg, struct position *lambdaarg, int num_coils, int *num_segs) {
//	double ***adjfunc = malloc(num_coils*sizeof(double));
//	int coil_ind, seg_ind, dir_ind;
//	for (coil_ind=0; coil_ind<num_coils; coil_ind++) {
//		adjfunc[coil_ind] = malloc(num_segs[coil_ind]*sizeof(double));
//		for (seg_ind=0; seg_ind<num_segs[coil_ind]; seg_ind++) {
//			adjfunc[coil_ind][seg_ind] = malloc(3*sizeof(double));
//			for (dir_ind=0; dir_ind<3; dir_ind++) {
//				adjfunc[coil_ind][seg_ind][dir_ind] = lambdaarg->loc[0] * Xparg->loc[0] * ( - shapeBparg[coil_ind][seg_ind][dir_ind].value[0] + shapeBparg[coil_ind][seg_ind][dir_ind].value[2]*Bparg->value[0] / Bparg->value[2] ) / Bparg->value[2] 
//				+ lambdaarg->loc[1] * Xparg->loc[0] * ( - shapeBparg[coil_ind][seg_ind][dir_ind].value[1] + shapeBparg[coil_ind][seg_ind][dir_ind].value[2]*Bparg->value[1] / Bparg->value[2] ) / Bparg->value[2] ;
//			}
//		}
//	}
//	return adjfunc;
//}

double *derivative_gradtangent(struct field *Bparg, struct field *gradBparg, struct position *Xparg, struct position *lambdaarg, struct position *sperp, struct position *mu, int num_params)
{
	double *adjfunc=calloc(num_params, sizeof(double)), extraterm, Rhat[2];
	int pro, row, col;
	//adjfunc = calloc(num_params,sizeof(double)); 
	Rhat[0] = 1.0; Rhat[1] = 0.0;
	for (pro=0;pro<num_params;pro++) {
		for (row=0;row<2;row++) {
			extraterm = 0.0;
			for (col=0;col<2;col++) {
				extraterm += ( sperp->loc[col] / Bparg->value[2] ) * ( Rhat[col]*gradBparg[pro].value[row] - Rhat[col]*gradBparg[pro].value[2] * Bparg->value[row] / Bparg->value[2] + Xparg->loc[0] * ( gradBparg->derivative[row][col] - gradBparg->derivative[2][col]*Bparg->value[row]/Bparg->value[2] - Bparg->derivative[2][col]*gradBparg->value[row]/Bparg->value[2] ) ) ;
			}
			adjfunc[pro] += ( lambdaarg->loc[row] * Xparg->loc[0] 
			* ( - gradBparg[pro].value[row] + gradBparg[pro].value[2]* Bparg->value[row] / Bparg->value[2] ) 
		/ Bparg->value[2] - mu->loc[row] * extraterm );
		//printf("lambda=%f\n", lambdaarg->loc[row]);
		//printf("mu=%f\n", mu->loc[row]);
		//printf("Xparg->loc=%f\n", Xparg->loc[row]);
		}
	}
	return adjfunc;
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
	//structsum.tangent = malloc(2*sizeof(double));
	//structsum.tangent[0] = malloc(2*sizeof(double)); structsum.tangent[1] = malloc(2*sizeof(double));
	structsum.tangent = set_identity();
	for (index1=0;index1<2;index1++) {
		structsum.loc[index1] = num1*struct1->loc[index1] + num2*struct2->loc[index1];
		for (index2=0;index2<2;index2++) {
			structsum.tangent[index1][index2] = num1*struct1->tangent[index1][index2] + num2*struct2->tangent[index1][index2];
		}
	}
	return structsum;
}

struct position addstructsreassign(double num1, struct position *struct1, double num2, struct position *struct2) {
	struct position structsum;
	int index1, index2;
	structsum.tangent = malloc(2*sizeof(double*));
	structsum.tangent[0] = malloc(2*sizeof(double)); structsum.tangent[1] = malloc(2*sizeof(double));
	for (index1=0;index1<2;index1++) {
		structsum.loc[index1] = num1*struct1->loc[index1] + num2*struct2->loc[index1];
		for (index2=0;index2<2;index2++)	
		structsum.tangent[index1][index2] = num1*struct1->tangent[index1][index2] + num2*struct2->tangent[index1][index2];
	}
	for (index1=0;index1<2;index1++) {
		struct1->loc[index1] = structsum.loc[index1];
		for (index2=0;index2<2;index2++)	
			struct1->tangent[index1][index2] = structsum.tangent[index1][index2];
	}
	free(structsum.tangent[0]); free(structsum.tangent[1]);
	free(structsum.tangent);
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

void RK4step(struct position *Xp, double varphi, double dvarphi, double ***coils, int num_coils, int *num_segs, struct field *Bfield_saved) {
	int row;
	struct field *Bpoint;
	struct position *dXp[4], Xpold;
	Xpold = *Xp;
	Bpoint = Bfield(Xp->loc, varphi, coils, num_coils, num_segs);
	Bfield_saved[0] = *Bpoint;
	//printf("magnetic field before 1st segment of RK iteratiion: Bfield=%f\n", Bpoint->value[0]);
	dXp[0] = derivative_center(Bpoint, Xp);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp[0]);
	//printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
	Bfield_saved[1] = *Bpoint;
	dXp[1] = derivative_center(Bpoint, Xp);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp[1]);
	//printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
	Bfield_saved[2] = *Bpoint;
	dXp[2] = derivative_center(Bpoint, Xp);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp[2]);
	//printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpoint = Bfield(Xp->loc, varphi+dvarphi, coils, num_coils, num_segs);
	Bfield_saved[3] = *Bpoint;
	dXp[3] = derivative_center(Bpoint, Xp);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp[0]);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp[1]);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp[2]);
	addstructsreassign(1.0, Xp,     (1.0/6.0)*dvarphi, dXp[3]);
	//printf("after 4th segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	for (row=0; row<4; row++) {
		free(dXp[row]->tangent[0]); free(dXp[row]->tangent[1]); 
		free(dXp[row]->tangent);
		free(dXp[row]);
	}
	return;
}

void RK4step_withsavedB(struct position *Xp, double varphi, double dvarphi, struct field *Bfield_saved)
{
	int row;
	struct position *dXp[4], Xpold;
	Xpold = *Xp;
	//printf("magnetic field before 1st segment of RK iteration: Bfield=%f\n", Bfield_saved[0].value[0]);
	dXp[0] = derivative_center(Bfield_saved, Xp);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp[0]);
	//printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	dXp[1] = derivative_center(Bfield_saved+1, Xp);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp[1]);
	//printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	dXp[2] = derivative_center(Bfield_saved+2, Xp);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp[2]);
	//printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	dXp[3] = derivative_center(Bfield_saved+3, Xp);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp[0]);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp[1]);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp[2]);
	addstructsreassign(1.0, Xp,     (1.0/6.0)*dvarphi, dXp[3]);
	//printf("after 4th segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	for (row=0; row<4; row++) {
		free(dXp[row]->tangent[0]); free(dXp[row]->tangent[1]); 
		free(dXp[row]->tangent);
		free(dXp[row]);
	}
	return;
}

void RK4step_lambdacirc_mutangent(struct position *Xp, struct position *adjvariable, double varphi, double dvarphi, struct field *Bfield_saved) {
	int row;
	struct field *Bpoint;
	struct position *dXp[4], *dadjvariable[4], Xpold, adjvariableold;
	Xpold.tangent = set_identity();
	adjvariableold.tangent = set_identity();
	for (row = 0; row < 2; row++) {
		Xpold.loc[row] = Xp->loc[row];
		Xpold.tangent[row][0] = Xp->tangent[row][0];
		Xpold.tangent[row][1] = Xp->tangent[row][1];
		adjvariableold.loc[row] = adjvariable->loc[row];
		adjvariableold.tangent[row][0] = adjvariable->tangent[row][0];
		adjvariableold.tangent[row][1] = adjvariable->tangent[row][1];
	}
//	printstructposition("Xpold", &Xpold);
	//printf("magnetic field before 1st segment of RK iteratiion: Bfield=%f\n", Bpoint->value[0]);
	//printf("Xp->loc[0] = %f\n", Xp->loc[0]);
	//printf("Xp->loc[1] = %f\n", Xp->loc[1]);
	//printf("Xp->tangent[0][0] = %f\n", Xp->tangent[0][0]);
	//printf("Xp->tangent[0][1] = %f\n", Xp->tangent[0][1]);
	//printf("Xp->tangent[1][0] = %f\n", Xp->tangent[1][0]);
	//printf("Xp->tangent[1][1] = %f\n", Xp->tangent[1][1]);
	//printf("Bpoint->value[0] = %f\n", Bpoint->value[0]);
	//printf("Bpoint->value[1] = %f\n", Bpoint->value[1]);
	//printf("Bpoint->value[2] = %f\n", Bpoint->value[2]);
	//printf("Bpoint->derivative[0][0] = %f\n", Bpoint->derivative[0][0]);
	//printf("Bpoint->derivative[1][0] = %f\n", Bpoint->derivative[1][0]);
	//printf("Bpoint->derivative[2][0] = %f\n", Bpoint->derivative[2][0]);
	//printf("Bpoint->derivative[0][1] = %f\n", Bpoint->derivative[0][1]);
	//printf("Bpoint->derivative[1][1] = %f\n", Bpoint->derivative[1][1]);
	//printf("Bpoint->derivative[2][1] = %f\n", Bpoint->derivative[2][1]);

	//Bpoint = Bfield(Xp->loc, varphi, coils, num_coils, num_segs);
	Bpoint = Bfield_saved;
	//printf("magnetic field before 1st segment of RK iteratiion: Bfield=%f\n", Bpoint->value[0]);
	dXp[0] = derivative_center(Bpoint, Xp);
	dadjvariable[0] = derivative_lambdacirc_mutangent(Bpoint, Xp, adjvariable);
	// I think this memory needs to be kept (because it corresponds to Xp
	// before the next iteration, which is preserved in the function that
	// calls this function. So the 2 lines below must remain commented out
	//free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	//free(adjvariable->tangent[0]); free(adjvariable->tangent[1]); free(adjvariable->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp[0]);
	*adjvariable = addstructs(1.0, &adjvariableold, 0.5*dvarphi, dadjvariable[0]);
	//printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	//Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
	Bpoint = Bfield_saved+1;
	dXp[1] = derivative_center(Bpoint, Xp);
	dadjvariable[1] = derivative_lambdacirc_mutangent(Bpoint, Xp, adjvariable);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(adjvariable->tangent[0]); free(adjvariable->tangent[1]); free(adjvariable->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp[1]);
	*adjvariable = addstructs(1.0, &adjvariableold, 0.5*dvarphi, dadjvariable[1]);
	//printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	//Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
	Bpoint = Bfield_saved+2;
	dXp[2] = derivative_center(Bpoint, Xp);
	dadjvariable[2] = derivative_lambdacirc_mutangent(Bpoint, Xp, adjvariable);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(adjvariable->tangent[0]); free(adjvariable->tangent[1]); free(adjvariable->tangent);
	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp[2]);
	*adjvariable = addstructs(1.0, &adjvariableold, dvarphi, dadjvariable[2]);
	//printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	//Bpoint = Bfield(Xp->loc, varphi+dvarphi, coils, num_coils, num_segs);
	Bpoint = Bfield_saved+3;
	dXp[3] = derivative_center(Bpoint, Xp);
	dadjvariable[3] = derivative_lambdacirc_mutangent(Bpoint, Xp, adjvariable);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(adjvariable->tangent[0]); free(adjvariable->tangent[1]); free(adjvariable->tangent);
	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp[0]);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp[1]);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp[2]);
	addstructsreassign(1.0, Xp,     (1.0/6.0)*dvarphi, dXp[3]);
	*adjvariable = addstructs(1.0, &adjvariableold, (1.0/6.0)*dvarphi, dadjvariable[0]);
	addstructsreassign(1.0, adjvariable,     (1.0/3.0)*dvarphi, dadjvariable[1]);
	addstructsreassign(1.0, adjvariable,     (1.0/3.0)*dvarphi, dadjvariable[2]);
	addstructsreassign(1.0, adjvariable,     (1.0/6.0)*dvarphi, dadjvariable[3]);
	for (row=0; row<4; row++) {
		free(dXp[row]->tangent[0]); free(dXp[row]->tangent[1]); 
		free(dadjvariable[row]->tangent[0]); free(dadjvariable[row]->tangent[1]); 
		free(dXp[row]->tangent);
		free(dadjvariable[row]->tangent);
		free(dXp[row]);
		free(dadjvariable[row]);
	}
	return;
}

void RK4step_lambdatangent(struct position *Xp, struct position *lambda, struct position *sperp, struct position *mu, double varphi, double dvarphi, struct field *Bfield_saved)
{
	int row;
	struct field *Bpoint;
	struct position *dXp[4], *dlambda[4], *dsperp[4], *dmu[4], Xpold, lambdaold, sperpold, muold;
	Xpold.tangent = set_identity();
	lambdaold.tangent = set_identity();
	sperpold.tangent = set_identity();
	muold.tangent = set_identity();
	for (row = 0; row < 2; row++) {
		Xpold.loc[row] = Xp->loc[row];
		Xpold.tangent[row][0] = Xp->tangent[row][0];
		Xpold.tangent[row][1] = Xp->tangent[row][1];
		lambdaold.loc[row] = lambda->loc[row];
		lambdaold.tangent[row][0] = lambda->tangent[row][0];
		lambdaold.tangent[row][1] = lambda->tangent[row][1];
		muold.loc[row] = mu->loc[row];
		muold.tangent[row][0] = mu->tangent[row][0];
		muold.tangent[row][1] = mu->tangent[row][1];
		sperpold.loc[row] = sperp->loc[row];
		sperpold.tangent[row][0] = sperp->tangent[row][0];
		sperpold.tangent[row][1] = sperp->tangent[row][1];
	}
	//Xpold = *Xp;
	//lambdaold = *lambda;
	//muold = *mu;
	//sperpold = *sperp;
	//printf("magnetic field before 1st segment of RK iteratiion: Bfield=%f\n", Bpoint->value[0]);

	//Bpoint = Bfield(Xp->loc, varphi, coils, num_coils, num_segs);
	Bpoint = Bfield_saved;
	//printf("magnetic field before 1st segment of RK iteratiion: Bfield=%f\n", Bpoint->value[0]);
	dXp[0] = derivative_center(Bpoint, Xp);
	dlambda[0] = derivative_lambdatangent(Bpoint, Xp, lambda, sperp, mu);
	//dlambda[0] = derivative_lambdacirc_mutangent(Bpoint, Xp, lambda);
	dsperp[0] = derivative_sperp(Bpoint, Xp, sperp);
	dmu[0] = derivative_lambdacirc_mutangent(Bpoint, Xp, mu);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	//free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	free(sperp->tangent[0]); free(sperp->tangent[1]); free(sperp->tangent);
	//free(mu->tangent[0]); free(mu->tangent[1]); free(mu->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp[0]);
	*lambda = addstructs(1.0, &lambdaold, 0.5*dvarphi, dlambda[0]);
	*sperp = addstructs(1.0, &sperpold, 0.5*dvarphi, dsperp[0]);
	*mu = addstructs(1.0, &muold, 0.5*dvarphi, dmu[0]);
	//printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	//Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
	Bpoint = Bfield_saved+1;
	dXp[1] = derivative_center(Bpoint, Xp);
	dlambda[1] = derivative_lambdatangent(Bpoint, Xp, lambda, sperp, mu);
	//dlambda[1] = derivative_lambdacirc_mutangent(Bpoint, Xp, lambda);
	dsperp[1] = derivative_sperp(Bpoint, Xp, sperp);
	dmu[1] = derivative_lambdacirc_mutangent(Bpoint, Xp, mu);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	free(sperp->tangent[0]); free(sperp->tangent[1]); free(sperp->tangent);
	free(mu->tangent[0]); free(mu->tangent[1]); free(mu->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp[1]);
	*lambda = addstructs(1.0, &lambdaold, 0.5*dvarphi, dlambda[1]);
	*sperp = addstructs(1.0, &sperpold, 0.5*dvarphi, dsperp[1]);
	*mu = addstructs(1.0, &muold, 0.5*dvarphi, dmu[1]);
	//printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	//Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
	Bpoint = Bfield_saved+2;
	dXp[2] = derivative_center(Bpoint, Xp);
	dlambda[2] = derivative_lambdatangent(Bpoint, Xp, lambda, sperp, mu);
	//dlambda[2] = derivative_lambdacirc_mutangent(Bpoint, Xp, lambda);
	dsperp[2] = derivative_sperp(Bpoint, Xp, sperp);
	dmu[2] = derivative_lambdacirc_mutangent(Bpoint, Xp, mu);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	free(sperp->tangent[0]); free(sperp->tangent[1]); free(sperp->tangent);
	free(mu->tangent[0]); free(mu->tangent[1]); free(mu->tangent);
	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp[2]);
	*lambda = addstructs(1.0, &lambdaold, dvarphi, dlambda[2]);
	*sperp = addstructs(1.0, &sperpold, dvarphi, dsperp[2]);
	*mu = addstructs(1.0, &muold, dvarphi, dmu[2]);
	//printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	//Bpoint = Bfield(Xp->loc, varphi+dvarphi, coils, num_coils, num_segs);
	Bpoint = Bfield_saved+3;
	dXp[3] = derivative_center(Bpoint, Xp);
	dlambda[3] = derivative_lambdatangent(Bpoint, Xp, lambda, sperp, mu);
	//dlambda[3] = derivative_lambdacirc_mutangent(Bpoint, Xp, lambda);
	dsperp[3] = derivative_sperp(Bpoint, Xp, sperp);
	dmu[3] = derivative_lambdacirc_mutangent(Bpoint, Xp, mu);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	free(sperp->tangent[0]); free(sperp->tangent[1]); free(sperp->tangent);
	free(mu->tangent[0]); free(mu->tangent[1]); free(mu->tangent);
	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp[0]);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp[1]);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp[2]);
	addstructsreassign(1.0, Xp,     (1.0/6.0)*dvarphi, dXp[3]);
	*lambda = addstructs(1.0, &lambdaold, (1.0/6.0)*dvarphi, dlambda[0]);
	addstructsreassign(1.0, lambda,     (1.0/3.0)*dvarphi, dlambda[1]);
	addstructsreassign(1.0, lambda,     (1.0/3.0)*dvarphi, dlambda[2]);
	addstructsreassign(1.0, lambda,     (1.0/6.0)*dvarphi, dlambda[3]);
	*sperp = addstructs(1.0, &sperpold, (1.0/6.0)*dvarphi, dsperp[0]);
	addstructsreassign(1.0, sperp,     (1.0/3.0)*dvarphi, dsperp[1]);
	addstructsreassign(1.0, sperp,     (1.0/3.0)*dvarphi, dsperp[2]);
	addstructsreassign(1.0, sperp,     (1.0/6.0)*dvarphi, dsperp[3]);
	*mu = addstructs(1.0, &muold, (1.0/6.0)*dvarphi, dmu[0]);
	addstructsreassign(1.0, mu,     (1.0/3.0)*dvarphi, dmu[1]);
	addstructsreassign(1.0, mu,     (1.0/3.0)*dvarphi, dmu[2]);
	addstructsreassign(1.0, mu,     (1.0/6.0)*dvarphi, dmu[3]);
	for (row=0; row<4; row++) {
		free(dXp[row]->tangent[0]); free(dXp[row]->tangent[1]); 
		free(dlambda[row]->tangent[0]); free(dlambda[row]->tangent[1]); 
		free(dsperp[row]->tangent[0]); free(dsperp[row]->tangent[1]); 
		free(dmu[row]->tangent[0]); free(dmu[row]->tangent[1]); 
		free(dXp[row]->tangent);
		free(dlambda[row]->tangent);
		free(dsperp[row]->tangent);
		free(dmu[row]->tangent);
		free(dXp[row]);
		free(dlambda[row]);
		free(dsperp[row]);
		free(dmu[row]);
	}
	return;
}

void RK4step_gradcirc(double *number, struct position *Xp, struct position *lambda, double varphi, double dvarphi, struct field *Bfield_saved, double *param, int num_params)
{
	int row;
	struct field *Bpoint, *Bpointgrad;
	struct position *dXp[4], Xpold;
	struct position *dlambda[4], lambdaold;
	double *dnumber[4];
	double numberold[num_params]; 
	Xpold.tangent = set_identity();
	lambdaold.tangent = set_identity();
	for (row = 0; row < 2; row++) {
		Xpold.loc[row] = Xp->loc[row];
		Xpold.tangent[row][0] = Xp->tangent[row][0];
		Xpold.tangent[row][1] = Xp->tangent[row][1];
		lambdaold.loc[row] = lambda->loc[row];
		lambdaold.tangent[row][0] = lambda->tangent[row][0];
		lambdaold.tangent[row][1] = lambda->tangent[row][1];
	}
	for (row=0;row<num_params;row++) 
		numberold[row] = number[row];
	Bpoint = Bfield_saved;
	Bpointgrad = gradBfield(Xp->loc, varphi, param);
	//printf("magnetic field before 1st segment of RK iteratiion: Bfield=%f\n", Bpoint->value[0]);
	dXp[0] = derivative_center(Bpoint, Xp);
	//printf("dXp[0]=%f\n", dXp[0]->loc[0]);
	dlambda[0] = derivative_lambdacirc_mutangent(Bpoint, Xp, lambda);
	dnumber[0] = derivative_gradcirc(Bpoint, Bpointgrad, Xp, lambda, num_params);
	free(Bpointgrad);
	//free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	//free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	//printf("1st segment of RK iteratiion: BR=%f, GXBR=%f, dR=%f, dlambdaR=%f, dnumberX=%f\n", Bpoint->value[0], Bpointgrad[0].value[0], dXp[0]->loc[0], dlambda[0]->loc[0], dnumber[0][0]);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp[0]);
	*lambda = addstructs(1.0, &lambdaold, 0.5*dvarphi, dlambda[0]);
	//printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	//Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
	Bpoint = Bfield_saved+1;
	Bpointgrad = gradBfield(Xp->loc, varphi, param);
	dXp[1] = derivative_center(Bpoint, Xp);
	dlambda[1] = derivative_lambdacirc_mutangent(Bpoint, Xp, lambda);
	dnumber[1] = derivative_gradcirc(Bpoint, Bpointgrad, Xp, lambda, num_params);
	free(Bpointgrad);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp[1]);
	*lambda = addstructs(1.0, &lambdaold, 0.5*dvarphi, dlambda[1]);
	//printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	//Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
	Bpoint = Bfield_saved+2;
	Bpointgrad = gradBfield(Xp->loc, varphi, param);
	dXp[2] = derivative_center(Bpoint, Xp);
	dlambda[2] = derivative_lambdacirc_mutangent(Bpoint, Xp, lambda);
	dnumber[2] = derivative_gradcirc(Bpoint, Bpointgrad, Xp, lambda, num_params);
	free(Bpointgrad);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp[2]);
	*lambda = addstructs(1.0, &lambdaold, dvarphi, dlambda[2]);
	//printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	//Bpoint = Bfield(Xp->loc, varphi+dvarphi, coils, num_coils, num_segs);
	Bpoint = Bfield_saved+3;
	Bpointgrad = gradBfield(Xp->loc, varphi, param);
	dXp[3] = derivative_center(Bpoint, Xp);
	dlambda[3] = derivative_lambdacirc_mutangent(Bpoint, Xp, lambda);
	dnumber[3] = derivative_gradcirc(Bpoint, Bpointgrad, Xp, lambda, num_params);
	free(Bpointgrad);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp[0]);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp[1]);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp[2]);
	addstructsreassign(1.0, Xp,     (1.0/6.0)*dvarphi, dXp[3]);
	*lambda = addstructs(1.0, &lambdaold, (1.0/6.0)*dvarphi, dlambda[0]);
	addstructsreassign(1.0, lambda,     (1.0/3.0)*dvarphi, dlambda[1]);
	addstructsreassign(1.0, lambda,     (1.0/3.0)*dvarphi, dlambda[2]);
	addstructsreassign(1.0, lambda,     (1.0/6.0)*dvarphi, dlambda[3]);
	for (row=0;row<num_params;row++) {
		number[row] = numberold[row] + dvarphi*(dnumber[0][row]/6.0 + dnumber[1][row]/3.0 + dnumber[2][row]/3.0 + dnumber[3][row]/6.0);
	}

	for (row=0; row<4; row++) {
		free(dXp[row]->tangent[0]); free(dXp[row]->tangent[1]); 
		free(dlambda[row]->tangent[0]); free(dlambda[row]->tangent[1]); 
		free(dXp[row]->tangent);
		free(dlambda[row]->tangent);
		free(dXp[row]);
		free(dlambda[row]);
		free(dnumber[row]);
	}
	//printf("number=%f\n", *number);
	//printf("numberold=%f\n", numberold);
	return;
}

//void RK4_adjshapecirc(double ***shapecirc, struct position *Xp, struct position *lambda, double varphi, double dvarphi, double ***coils, int num_coils, int *num_segs)
//{
//	struct field *Bpoint, ***Bpointshape;
//	struct position *dXp[4], Xpold;
//	struct position *dlambda[4], lambdaold;
//	double ***dnumber[4];
//	int coil_ind, seg_ind, dir_ind;
//	//printf("START: Runge Kutta \n");
//	Xpold = *Xp;
//	lambdaold = *lambda;
//	Bpoint = Bfield(Xp->loc, varphi, coils, num_coils, num_segs);
//	Bpointshape = shapeBfield(Xp->loc, varphi, coils, num_coils, num_segs);
//	//printf("magnetic field before 1st segment of RK iteratiion: Bfield=%f\n", Bpoint->value[0]);
//	dXp[0] = derivative_center(Bpoint, Xp);
//	dlambda[0] = derivative_lambdacirc_mutangent(Bpoint, Xp, lambda);
//	dnumber[0] = derivative_adjshape(Bpoint, Bpointshape, Xp, lambda, num_coils, num_segs);
//	//printf("dnumber1=%f\n", dnumber1);
//	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp[0]);
//	*lambda = addstructs(1.0, &lambdaold, 0.5*dvarphi, dlambda[0]);
//	//printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
//	Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
//	Bpointshape = shapeBfield(Xp->loc, varphi, coils, num_coils, num_segs);
//	dXp[1] = derivative_center(Bpoint, Xp);
//	dlambda[1] = derivative_lambdacirc_mutangent(Bpoint, Xp, lambda);
//	dnumber[1] = derivative_adjshape(Bpoint, Bpointshape, Xp, lambda, num_coils, num_segs);
//	//printf("dnumber2=%f\n", dnumber2);
//	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp[1]);
//	*lambda = addstructs(1.0, &lambdaold, 0.5*dvarphi, dlambda[1]);
//	//printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
//	Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
//	Bpointshape = shapeBfield(Xp->loc, varphi, coils, num_coils, num_segs);
//	dXp[2] = derivative_center(Bpoint, Xp);
//	dlambda[2] = derivative_lambdacirc_mutangent(Bpoint, Xp, lambda);
//	dnumber[2] = derivative_adjshape(Bpoint, Bpointshape, Xp, lambda, num_coils, num_segs);
//	//printf("dnumber3=%f\n", dnumber3);
//	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp[2]);
//	*lambda = addstructs(1.0, &lambdaold, dvarphi, dlambda[2]);
//	//printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
//	Bpoint = Bfield(Xp->loc, varphi+dvarphi, coils, num_coils, num_segs);
//	Bpointshape = shapeBfield(Xp->loc, varphi, coils, num_coils, num_segs);
//	dXp[3] = derivative_center(Bpoint, Xp);
//	dlambda[3] = derivative_lambdacirc_mutangent(Bpoint, Xp, lambda);
//	dnumber[3] = derivative_adjshape(Bpoint, Bpointshape, Xp, lambda, num_coils, num_segs);
//	//printf("dnumber4=%f\n", dnumber4);
//	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp[0]);
//	*Xp = addstructs(1.0, Xp,     (1.0/3.0)*dvarphi, dXp[1]);
//	*Xp = addstructs(1.0, Xp,     (1.0/3.0)*dvarphi, dXp[2]);
//	*Xp = addstructs(1.0, Xp,     (1.0/6.0)*dvarphi, dXp[3]);
//	*lambda = addstructs(1.0, &lambdaold, (1.0/6.0)*dvarphi, dlambda[0]);
//	*lambda = addstructs(1.0, lambda,     (1.0/3.0)*dvarphi, dlambda[1]);
//	*lambda = addstructs(1.0, lambda,     (1.0/3.0)*dvarphi, dlambda[2]);
//	*lambda = addstructs(1.0, lambda,     (1.0/6.0)*dvarphi, dlambda[3]);
//	for (coil_ind=0; coil_ind<num_coils; coil_ind++) {
//		for (seg_ind=0; seg_ind<num_segs[coil_ind]; seg_ind++) {
//			for (dir_ind=0; dir_ind<3; dir_ind++) {
//				shapecirc[coil_ind][seg_ind][dir_ind] += ( dvarphi*(dnumber[0][coil_ind][seg_ind][dir_ind]/6.0 + dnumber[1][coil_ind][seg_ind][dir_ind]/3.0 + dnumber[2][coil_ind][seg_ind][dir_ind]/3.0 + dnumber[3][coil_ind][seg_ind][dir_ind]/6.0) );
//			}
//			//printf("DONE: seg_ind=%d/%d\n", seg_ind, num_segs[coil_ind]);
//		}
//		//printf("DONE: coil_ind=%d/%d\n", coil_ind, num_coils);
//	}
//	//printf("DONE: Runge Kutta \n");
//	return;
//}

void RK4step_gradtangent(double *number, struct position *Xp, struct position *lambda, struct position *sperp, struct position *mu, double varphi, double dvarphi, struct field *Bfield_saved, double *param, int num_params) {
	int row;
	struct field *Bpoint, *Bpointgrad;
	struct position *dXp[4], *dlambda[4], *dsperp[4], *dmu[4], Xpold, lambdaold, sperpold, muold;
	double *dnumber[4];
	double numberold[num_params];
	Xpold.tangent = set_identity();
	lambdaold.tangent = set_identity();
	muold.tangent = set_identity();
	sperpold.tangent = set_identity();
	for (row = 0; row < 2; row++) {
		Xpold.loc[row] = Xp->loc[row];
		Xpold.tangent[row][0] = Xp->tangent[row][0];
		Xpold.tangent[row][1] = Xp->tangent[row][1];
		lambdaold.loc[row] = lambda->loc[row];
		lambdaold.tangent[row][0] = lambda->tangent[row][0];
		lambdaold.tangent[row][1] = lambda->tangent[row][1];
		muold.loc[row] = mu->loc[row];
		muold.tangent[row][0] = mu->tangent[row][0];
		muold.tangent[row][1] = mu->tangent[row][1];
		sperpold.loc[row] = sperp->loc[row];
		sperpold.tangent[row][0] = sperp->tangent[row][0];
		sperpold.tangent[row][1] = sperp->tangent[row][1];
	}
	for (row = 0; row < num_params; row++) 
		numberold[row] = number[row];
	//printf("magnetic field before 1st segment of RK iteratiion: Bfield=%f\n", Bpoint->value[0]);
	//Bpoint = Bfield(Xp->loc, varphi, coils, num_coils, num_segs);
	Bpoint = Bfield_saved;
	Bpointgrad = gradBfield(Xp->loc, varphi, param);
	//printf("magnetic field before 1st segment of RK iteratiion: Bfield=%f\n", Bpoint->value[0]);
	dXp[0] = derivative_center(Bpoint, Xp);
	dlambda[0] = derivative_lambdatangent(Bpoint, Xp, lambda, sperp, mu);
	dsperp[0] = derivative_sperp(Bpoint, Xp, sperp);
	dmu[0] = derivative_lambdacirc_mutangent(Bpoint, Xp, mu);
	dnumber[0] = derivative_gradtangent(Bpoint, Bpointgrad, Xp, lambda, sperp, mu, num_params);
	free(Bpointgrad);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	//free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	free(sperp->tangent[0]); free(sperp->tangent[1]); free(sperp->tangent);
	//free(mu->tangent[0]); free(mu->tangent[1]); free(mu->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp[0]);
	*lambda = addstructs(1.0, &lambdaold, 0.5*dvarphi, dlambda[0]);
	*sperp = addstructs(1.0, &sperpold, 0.5*dvarphi, dsperp[0]);
	*mu = addstructs(1.0, &muold, 0.5*dvarphi, dmu[0]);
	//printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	//Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
	Bpoint = Bfield_saved+1;
	Bpointgrad = gradBfield(Xp->loc, varphi + 0.5*dvarphi, param);
	dXp[1] = derivative_center(Bpoint, Xp);
	dlambda[1] = derivative_lambdatangent(Bpoint, Xp, lambda, sperp, mu);
	dsperp[1] = derivative_sperp(Bpoint, Xp, sperp);
	dmu[1] = derivative_lambdacirc_mutangent(Bpoint, Xp, mu);
	dnumber[1] = derivative_gradtangent(Bpoint, Bpointgrad, Xp, lambda, sperp, mu, num_params);
	free(Bpointgrad);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	free(sperp->tangent[0]); free(sperp->tangent[1]); free(sperp->tangent);
	free(mu->tangent[0]); free(mu->tangent[1]); free(mu->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp[1]);
	*lambda = addstructs(1.0, &lambdaold, 0.5*dvarphi, dlambda[1]);
	*sperp = addstructs(1.0, &sperpold, 0.5*dvarphi, dsperp[1]);
	*mu = addstructs(1.0, &muold, 0.5*dvarphi, dmu[1]);
	//printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	//Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
	Bpoint = Bfield_saved+2;
	Bpointgrad = gradBfield(Xp->loc, varphi + 0.5*dvarphi, param);
	dXp[2] = derivative_center(Bpoint, Xp);
	dlambda[2] = derivative_lambdatangent(Bpoint, Xp, lambda, sperp, mu);
	dsperp[2] = derivative_sperp(Bpoint, Xp, sperp);
	dmu[2] = derivative_lambdacirc_mutangent(Bpoint, Xp, mu);
	dnumber[2] = derivative_gradtangent(Bpoint, Bpointgrad, Xp, lambda, sperp, mu, num_params);
	free(Bpointgrad);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	free(sperp->tangent[0]); free(sperp->tangent[1]); free(sperp->tangent);
	free(mu->tangent[0]); free(mu->tangent[1]); free(mu->tangent);
	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp[2]);
	*lambda = addstructs(1.0, &lambdaold, dvarphi, dlambda[2]);
	*sperp = addstructs(1.0, &sperpold, dvarphi, dsperp[2]);
	*mu = addstructs(1.0, &muold, dvarphi, dmu[2]);
	//printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	//Bpoint = Bfield(Xp->loc, varphi+dvarphi, coils, num_coils, num_segs);
	Bpoint = Bfield_saved+3;
	Bpointgrad = gradBfield(Xp->loc, varphi, param);
	dXp[3] = derivative_center(Bpoint, Xp);
	dlambda[3] = derivative_lambdatangent(Bpoint, Xp, lambda, sperp, mu);
	dsperp[3] = derivative_sperp(Bpoint, Xp, sperp);
	dmu[3] = derivative_lambdacirc_mutangent(Bpoint, Xp, mu);
	dnumber[3] = derivative_gradtangent(Bpoint, Bpointgrad, Xp, lambda, sperp, mu, num_params);
	free(Bpointgrad);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	free(sperp->tangent[0]); free(sperp->tangent[1]); free(sperp->tangent);
	free(mu->tangent[0]); free(mu->tangent[1]); free(mu->tangent);
	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp[0]);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp[1]);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp[2]);
	addstructsreassign(1.0, Xp,     (1.0/6.0)*dvarphi, dXp[3]);
	*lambda = addstructs(1.0, &lambdaold, (1.0/6.0)*dvarphi, dlambda[0]);
	addstructsreassign(1.0, lambda,     (1.0/3.0)*dvarphi, dlambda[1]);
	addstructsreassign(1.0, lambda,     (1.0/3.0)*dvarphi, dlambda[2]);
	addstructsreassign(1.0, lambda,     (1.0/6.0)*dvarphi, dlambda[3]);
	*sperp = addstructs(1.0, &sperpold, (1.0/6.0)*dvarphi, dsperp[0]);
	addstructsreassign(1.0, sperp,     (1.0/3.0)*dvarphi, dsperp[1]);
	addstructsreassign(1.0, sperp,     (1.0/3.0)*dvarphi, dsperp[2]);
	addstructsreassign(1.0, sperp,     (1.0/6.0)*dvarphi, dsperp[3]);
	*mu = addstructs(1.0, &muold, (1.0/6.0)*dvarphi, dmu[0]);
	addstructsreassign(1.0, mu,     (1.0/3.0)*dvarphi, dmu[1]);
	addstructsreassign(1.0, mu,     (1.0/3.0)*dvarphi, dmu[2]);
	addstructsreassign(1.0, mu,     (1.0/6.0)*dvarphi, dmu[3]);
	//*number = numberold + dvarphi*(dnumber[0]/6.0 + dnumber[1]/3.0 + dnumber[2]/3.0 + dnumber[3]/6.0);
	for (row=0;row<num_params;row++) {
		number[row] = numberold[row] + dvarphi*(dnumber[0][row]/6.0 + dnumber[1][row]/3.0 + dnumber[2][row]/3.0 + dnumber[3][row]/6.0);
	}
	//printf("dnumber = %f\n", dvarphi*(dnumber[0]/6.0 + dnumber[1]/3.0 + dnumber[2]/3.0 + dnumber[3]/6.0));
	for (row=0; row<4; row++) {
		free(dXp[row]->tangent[0]); free(dXp[row]->tangent[1]); 
		free(dlambda[row]->tangent[0]); free(dlambda[row]->tangent[1]); 
		free(dsperp[row]->tangent[0]); free(dsperp[row]->tangent[1]); 
		free(dmu[row]->tangent[0]); free(dmu[row]->tangent[1]); 
		free(dXp[row]->tangent);
		free(dlambda[row]->tangent);
		free(dsperp[row]->tangent);
		free(dmu[row]->tangent);
		free(dXp[row]);
		free(dlambda[row]);
		free(dsperp[row]);
		free(dmu[row]);
		free(dnumber[row]);
	}
	return;
}

//void RK4_adjgradtangent(double *number, struct position *Xp, struct position *lambda, struct position *sperp, struct position *mu, double varphi, double dvarphi, double ***coils, int num_coils, int *num_segs) {
//	struct field *Bpoint, *Bpointgrad;
//	struct position *dXp[4], *dlambda[4], *dsperp[4], *dmu[4], Xpold, lambdaold, sperpold, muold;
//	double numberold, dnumber[4];
//	numberold = *number;
//	Xpold = *Xp;
//	lambdaold = *lambda;
//	muold = *mu;
//	sperpold = *sperp;
//	//printf("magnetic field before 1st segment of RK iteratiion: Bfield=%f\n", Bpoint->value[0]);
//	Bpoint = Bfield(Xp->loc, varphi, coils, num_coils, num_segs);
//	Bpointgrad = gradBfield(Xp->loc, varphi);
//	//printf("magnetic field before 1st segment of RK iteratiion: Bfield=%f\n", Bpoint->value[0]);
//	dXp[0] = forward(Bpoint, Xp);
//	dlambda[0] = forward_adjtangent(Bpoint, Xp, lambda, sperp, mu);
//	dsperp[0] = forward_sperp(Bpoint, Xp, sperp);
//	dmu[0] = forward_adjsimple(Bpoint, Xp, mu);
//	dnumber[0] = forward_adjgradtangent(Bpoint, Bpointgrad, Xp, lambda, sperp, mu);
//	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp[0]);
//	*lambda = addstructs(1.0, &lambdaold, 0.5*dvarphi, dlambda[0]);
//	*sperp = addstructs(1.0, &sperpold, 0.5*dvarphi, dsperp[0]);
//	*mu = addstructs(1.0, &muold, 0.5*dvarphi, dmu[0]);
//	//printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
//	Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
//	Bpointgrad = gradBfield(Xp->loc, varphi + 0.5*dvarphi);
//	dXp[1] = forward(Bpoint, Xp);
//	dlambda[1] = forward_adjtangent(Bpoint, Xp, lambda, sperp, mu);
//	dsperp[1] = forward_sperp(Bpoint, Xp, sperp);
//	dmu[1] = forward_adjsimple(Bpoint, Xp, mu);
//	dnumber[1] = forward_adjgradtangent(Bpoint, Bpointgrad, Xp, lambda, sperp, mu);
//	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp[1]);
//	*lambda = addstructs(1.0, &lambdaold, 0.5*dvarphi, dlambda[1]);
//	*sperp = addstructs(1.0, &sperpold, 0.5*dvarphi, dsperp[1]);
//	*lambda = addstructs(1.0, &muold, 0.5*dvarphi, dmu[1]);
//	//printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
//	Bpoint = Bfield(Xp->loc, varphi + 0.5*dvarphi, coils, num_coils, num_segs);
//	Bpointgrad = gradBfield(Xp->loc, varphi + 0.5*dvarphi);
//	dXp[2] = forward(Bpoint, Xp);
//	dlambda[2] = forward_adjtangent(Bpoint, Xp, lambda, sperp, mu);
//	dsperp[2] = forward_sperp(Bpoint, Xp, sperp);
//	dmu[2] = forward_adjsimple(Bpoint, Xp, mu);
//	dnumber[2] = forward_adjgradtangent(Bpoint, Bpointgrad, Xp, lambda, sperp, mu);
//	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp[2]);
//	*lambda = addstructs(1.0, &lambdaold, dvarphi, dlambda[2]);
//	*sperp = addstructs(1.0, &sperpold, dvarphi, dsperp[2]);
//	*mu = addstructs(1.0, &muold, dvarphi, dmu[2]);
//	//printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
//	Bpoint = Bfield(Xp->loc, varphi+dvarphi, coils, num_coils, num_segs);
//	Bpointgrad = gradBfield(Xp->loc, varphi);
//	dXp[3] = forward(Bpoint, Xp);
//	dlambda[3] = forward_adjtangent(Bpoint, Xp, lambda, sperp, mu);
//	dsperp[3] = forward_sperp(Bpoint, Xp, sperp);
//	dmu[3] = forward_adjsimple(Bpoint, Xp, mu);
//	dnumber[3] = forward_adjgradtangent(Bpoint, Bpointgrad, Xp, lambda, sperp, mu);
//	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp[0]);
//	*Xp = addstructs(1.0, Xp,     (1.0/3.0)*dvarphi, dXp[1]);
//	*Xp = addstructs(1.0, Xp,     (1.0/3.0)*dvarphi, dXp[2]);
//	*Xp = addstructs(1.0, Xp,     (1.0/6.0)*dvarphi, dXp[3]);
//	*lambda = addstructs(1.0, &lambdaold, (1.0/6.0)*dvarphi, dlambda[0]);
//	*lambda = addstructs(1.0, lambda,     (1.0/3.0)*dvarphi, dlambda[1]);
//	*lambda = addstructs(1.0, lambda,     (1.0/3.0)*dvarphi, dlambda[2]);
//	*lambda = addstructs(1.0, lambda,     (1.0/6.0)*dvarphi, dlambda[3]);
//	*sperp = addstructs(1.0, &sperpold, (1.0/6.0)*dvarphi, dsperp[0]);
//	*sperp = addstructs(1.0, sperp,     (1.0/3.0)*dvarphi, dsperp[1]);
//	*sperp = addstructs(1.0, sperp,     (1.0/3.0)*dvarphi, dsperp[2]);
//	*sperp = addstructs(1.0, sperp,     (1.0/6.0)*dvarphi, dsperp[3]);
//	*mu = addstructs(1.0, &muold, (1.0/6.0)*dvarphi, dmu[0]);
//	*mu = addstructs(1.0, mu,     (1.0/3.0)*dvarphi, dmu[1]);
//	*mu = addstructs(1.0, mu,     (1.0/3.0)*dvarphi, dmu[2]);
//	*mu = addstructs(1.0, mu,     (1.0/6.0)*dvarphi, dmu[3]);
//	*number = numberold + dvarphi*(dnumber[0]/6.0 + dnumber[1]/3.0 + dnumber[2]/3.0 + dnumber[3]/6.0);
//	//printf("dnumber = %f\n", dvarphi*(dnumber[0]/6.0 + dnumber[1]/3.0 + dnumber[2]/3.0 + dnumber[3]/6.0));
//	return;
//}
