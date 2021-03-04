/*Author: Alessandro Geraldini 
  Description: this file contains functions that step forward in toroidal angle
               the functions named derivative contain the derivatives
               the functions named RK4 perform the step forward */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "isc.h"

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

struct position derivative_center(struct field *Bparg, struct position *Xparg)
{
	struct position funcreturn;
	double Jacobian[2][2], Dirac;
	int row, col;
	funcreturn.tangent = malloc(2*sizeof(double*));
	funcreturn.tangent[0] = malloc(2*sizeof(double)); funcreturn.tangent[1] = malloc(2*sizeof(double));
	for (row=0;row<2;row++)
	{	
		for (col=0;col<2;col++)
		{
			if (col==0)	Dirac=1.0;
			else		Dirac=0.0;
			Jacobian[row][col] = Dirac*Bparg->value[row]/Bparg->value[2] 
			+ Xparg->loc[0]*Bparg->derivative[row][col]/Bparg->value[2] 
			- Xparg->loc[0]*Bparg->value[row]*Bparg->derivative[2][col]/pow(Bparg->value[2], 2.0); 
		}
	}
	for (row=0;row<2;row++)
	{
		funcreturn.loc[row] = Xparg->loc[0]*Bparg->value[row]/Bparg->value[2];
		for (col=0;col<2;col++)
		{	
			funcreturn.tangent[row][col] = Jacobian[row][0]*Xparg->tangent[0][col] + Jacobian[row][1]*Xparg->tangent[1][col];
		}
	}
	return funcreturn;
}

struct position derivative_sperp(struct field *Bparg, struct position *Xparg, struct position *sperp)
{
	struct position funcreturn;
	double Jacobian[2][2], Dirac;
	int row, col;
	funcreturn.tangent = malloc(2*sizeof(double*));
	funcreturn.tangent[0] = malloc(2*sizeof(double)); funcreturn.tangent[1] = malloc(2*sizeof(double));
	for (row=0;row<2;row++)
	{	
		for (col=0;col<2;col++)
		{
			if (col==0)	Dirac=1.0;
			else		Dirac=0.0;
			Jacobian[row][col] = Dirac*Bparg->value[row]/Bparg->value[2] 
			+ Xparg->loc[0]*Bparg->derivative[row][col]/Bparg->value[2] 
			- Xparg->loc[0]*Bparg->value[row]*Bparg->derivative[2][col]/pow(Bparg->value[2], 2.0); 
		}
	}
	for (row=0;row<2;row++)
	{
		funcreturn.loc[row] = Jacobian[row][0]*sperp->loc[0] + Jacobian[row][1]*sperp->loc[1];
		for (col=0;col<2;col++)
		{	
			// tangent matrix unnecessary here
			//funcreturn->tangent[row][col] = Jacobian[row][0]*sperp->tangent[0][col] + Jacobian[row][1]*sperp->tangent[1][col];
		}
	}
	return funcreturn;
}

struct position derivative_lambdacirc_mutangent(struct field *Bparg, struct position *Xparg, struct position *adjvariable) {
	struct position adjfuncreturn;
	double Rhat[2];
	int row, col;
	Rhat[0] =1.0; Rhat[1] =0.0;
	adjfuncreturn.tangent = malloc(2*sizeof(double*));
	adjfuncreturn.tangent[0] = malloc(2*sizeof(double)); adjfuncreturn.tangent[1] = malloc(2*sizeof(double));
	for (row=0;row<2;row++)
	{
		adjfuncreturn.loc[row] = 0.0;
		for (col=0;col<2;col++)
		{	
			adjfuncreturn.loc[row] -=  ( adjvariable->loc[col] / Bparg->value[2] ) * ( Rhat[row]*Bparg->value[col] + Xparg->loc[0] * ( Bparg->derivative[col][row] - Bparg->derivative[2][row]*Bparg->value[col]/Bparg->value[2] ) );
			adjfuncreturn.tangent[row][col] = - ( Rhat[row]*Bparg->value[0]*adjvariable->tangent[0][col]/Bparg->value[2] + Rhat[row]*Bparg->value[1]*adjvariable->tangent[1][col]/Bparg->value[2] 
			+ (Xparg->loc[0]/Bparg->value[2]) * (Bparg->derivative[0][row]*adjvariable->tangent[0][col] + Bparg->derivative[1][row]*adjvariable->tangent[1][col] 
			 - Bparg->derivative[2][row]*Bparg->value[0]*adjvariable->tangent[0][col]/Bparg->value[2] - Bparg->derivative[2][row]*Bparg->value[1]*adjvariable->tangent[1][col]/Bparg->value[2] ) ) ;
		}
	}
	return adjfuncreturn;
}

struct position derivative_lambdatangent(struct field *Bp, struct position *Xp, struct position *lambdaQ, struct position *sperp, struct position *mu) {
	struct position adjfuncreturn;
	double Rhat[2];
	double object[2][2];
	int row, col, third;
	Rhat[0] =1.0; Rhat[1] =0.0;
	adjfuncreturn.tangent = malloc(2*sizeof(double*));
	adjfuncreturn.tangent[0] = malloc(2*sizeof(double)); adjfuncreturn.tangent[1] = malloc(2*sizeof(double));
	for (row=0;row<2;row++)
	{
		adjfuncreturn.loc[row] = 0.0;
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
			adjfuncreturn.loc[row] -=  ( ( lambdaQ->loc[col] / Bp->value[2] ) * ( Rhat[row]*Bp->value[col] + Xp->loc[0] * ( Bp->derivative[col][row] - Bp->derivative[2][row]*Bp->value[col]/Bp->value[2] ) ) + object[col][row]*sperp->loc[col] );
			adjfuncreturn.tangent[row][col] = - ( Rhat[row]*Bp->value[0]*lambdaQ->tangent[0][col]/Bp->value[2] + Rhat[row]*Bp->value[1]*lambdaQ->tangent[1][col]/Bp->value[2] 
			+ (Xp->loc[0]/Bp->value[2]) * (Bp->derivative[0][row]*lambdaQ->tangent[0][col] + Bp->derivative[1][row]*lambdaQ->tangent[1][col] 
			 - Bp->derivative[2][row]*Bp->value[0]*lambdaQ->tangent[0][col]/Bp->value[2] - Bp->derivative[2][row]*Bp->value[1]*lambdaQ->tangent[1][col]/Bp->value[2] ) ) ;
		}
	}
	return adjfuncreturn;
}

struct position derivative_lambdamu(struct field *Bp, struct position *Xp, struct position *lambdamu) {
   struct position adjfuncreturn;
   double Rhat[2];
   double object[2][2][2];
   int row, col, dep, third;
   Rhat[0] =1.0; Rhat[1] =0.0;
   adjfuncreturn.tangent = malloc(2*sizeof(double*));
   adjfuncreturn.tangent[0] = malloc(2*sizeof(double)); adjfuncreturn.tangent[1] = malloc(2*sizeof(double));
   for (row=0;row<2;row++) {
      adjfuncreturn.loc[row] = 0.0;
      for (col=0;col<2;col++) {	
         adjfuncreturn.tangent[row][col] = 0.0;
         for (dep=0;dep<2; dep++) {
	    object[col][row][dep] = 0.0;
	    for (third=0;third<2;third++) {
	       object[col][row][dep] 
	       += ( Rhat[col]*( Bp->derivative[third][row]/Bp->value[2] - Bp->value[third]*Bp->derivative[2][row]/ pow(Bp->value[2], 2.0) ) + Xp->loc[0]*( Bp->twoderivative[third][col][row]/(2.0*Bp->value[2]) - Bp->value[third]*Bp->twoderivative[2][row][col]/(2.0*pow(Bp->value[2], 2.0)) + Bp->value[third]*Bp->derivative[2][col]*Bp->derivative[2][row]/pow(Bp->value[2], 3.0) - Bp->derivative[third][col]*Bp->derivative[2][row]/pow(Bp->value[2], 2.0) )
	       +    Rhat[row]*( Bp->derivative[third][col]/Bp->value[2] - Bp->value[third]*Bp->derivative[2][col]/ pow(Bp->value[2], 2.0) ) + Xp->loc[0]*( Bp->twoderivative[third][row][col]/(2.0*Bp->value[2]) - Bp->value[third]*Bp->twoderivative[2][col][row]/(2.0*pow(Bp->value[2], 2.0)) + Bp->value[third]*Bp->derivative[2][row]*Bp->derivative[2][col]/pow(Bp->value[2], 3.0) - Bp->derivative[third][row]*Bp->derivative[2][col]/pow(Bp->value[2], 2.0) ) ) * lambdamu->tangent[third][dep]; 
	    }
	    adjfuncreturn.loc[row] -= ( object[col][row][dep]*Xp->tangent[col][dep] ); 
	 }
	 adjfuncreturn.loc[row] -=   ( lambdamu->loc[col] / Bp->value[2] ) * ( Rhat[row]*Bp->value[col] + Xp->loc[0] * ( Bp->derivative[col][row] - Bp->derivative[2][row]*Bp->value[col]/Bp->value[2] ) ) ;
	 adjfuncreturn.tangent[row][col] = - ( Rhat[row]*Bp->value[0]*lambdamu->tangent[0][col]/Bp->value[2] + Rhat[row]*Bp->value[1]*lambdamu->tangent[1][col]/Bp->value[2] + (Xp->loc[0]/Bp->value[2]) * (Bp->derivative[0][row]*lambdamu->tangent[0][col] + Bp->derivative[1][row]*lambdamu->tangent[1][col] - Bp->derivative[2][row]*Bp->value[0]*lambdamu->tangent[0][col]/Bp->value[2] - Bp->derivative[2][row]*Bp->value[1]*lambdamu->tangent[1][col]/Bp->value[2] ) ) ;
      }
   }
   return adjfuncreturn;
}

double *derivative_gradcirc(struct field *Bparg, struct field *gradBparg, struct position *Xparg, struct position *lambdaarg, int num_params)
{
	double *adjfunc = calloc(num_params,sizeof(double));
	int row, col;
	//printf("num_params=%d\n", num_params);
	for (row=0;row<num_params;row++) {
		adjfunc[row] = 0.0;
		for (col=0;col<2;col++) {
			adjfunc[row] += ( lambdaarg->loc[col] * Xparg->loc[0] * ( - gradBparg[row].value[col] + gradBparg[row].value[2]*Bparg->value[col] / Bparg->value[2] ) / Bparg->value[2] );
		}
	}
	return adjfunc;
}

double **derivative_gradXp(struct field *Bparg, struct field *gradBparg, struct position *Xparg, struct position *lambdaarg, int num_params)
{
	double **adjfunc = malloc(2*sizeof(double));
	int row, col;
	adjfunc[0] = calloc(num_params,sizeof(double));
	adjfunc[1] = calloc(num_params,sizeof(double));
	//printf("num_params=%d\n", num_params);
	for (row=0;row<num_params;row++) {
		adjfunc[0][row] = 0.0;
		adjfunc[1][row] = 0.0;
		for (col=0;col<2;col++) {
			adjfunc[0][row] += ( lambdaarg->tangent[col][0] * Xparg->loc[0] * ( - gradBparg[row].value[col] + gradBparg[row].value[2]*Bparg->value[col] / Bparg->value[2] ) / Bparg->value[2] );
			adjfunc[1][row] += ( lambdaarg->tangent[col][1] * Xparg->loc[0] * ( - gradBparg[row].value[col] + gradBparg[row].value[2]*Bparg->value[col] / Bparg->value[2] ) / Bparg->value[2] );
		}
	}
	return adjfunc;
}

double *derivative_gradtangent(struct field *Bparg, struct field *gradBparg, struct position *Xparg, struct position *lambdaarg, struct position *sperp, struct position *mu, int num_params)
{
	double *adjfunc=calloc(num_params, sizeof(double)), extraterm, Rhat[2];
	int pro, row, col;
	Rhat[0] = 1.0; Rhat[1] = 0.0;
	for (pro=0;pro<num_params;pro++) {
		for (row=0;row<2;row++) {
			extraterm = 0.0;
			for (col=0;col<2;col++) {
		    		extraterm += ( sperp->loc[col] / Bparg->value[2] ) * ( Rhat[col]*gradBparg[pro].value[row] - Rhat[col]*gradBparg[pro].value[2] * Bparg->value[row] / Bparg->value[2] + Xparg->loc[0] * ( gradBparg[pro].derivative[row][col] - gradBparg[pro].derivative[2][col]*Bparg->value[row]/Bparg->value[2] - Bparg->derivative[2][col]*gradBparg[pro].value[row]/Bparg->value[2] - Bparg->derivative[row][col]*gradBparg[pro].value[2]/Bparg->value[2] + 2.0*Bparg->derivative[2][col]*Bparg->value[row]*gradBparg[pro].value[2]/(Bparg->value[2]*Bparg->value[2]) ) ) ;
			}
			adjfunc[pro] += ( lambdaarg->loc[row] * Xparg->loc[0] 
			* ( - gradBparg[pro].value[row] + gradBparg[pro].value[2]* Bparg->value[row] / Bparg->value[2] ) 
		/ Bparg->value[2] - mu->loc[row] * extraterm );
		}
	}
	return adjfunc;
}

double *derivative_gradRes(struct field *Bparg, struct field *gradBparg, struct position *Xparg, struct position *lambdaarg, int num_params) {
    double *adjfunc=malloc(num_params*sizeof(double)), extraterm, Rhat[2];
    int pro, row, col, dep;
    Rhat[0] = 1.0; Rhat[1] = 0.0;
    for (pro=0;pro<num_params;pro++) {
   	adjfunc[pro] = 0.0;
   	for (col=0;col<2;col++) {
	    for (dep = 0; dep<2; dep++) {
	    	extraterm = 0.0;
	    	for (row=0;row<2;row++) {
		    extraterm += ( lambdaarg->tangent[row][dep] / Bparg->value[2] ) * ( Rhat[col]*gradBparg[pro].value[row] - Rhat[col]*gradBparg[pro].value[2] * Bparg->value[row] / Bparg->value[2] + Xparg->loc[0] * ( gradBparg[pro].derivative[row][col] - gradBparg[pro].derivative[2][col]*Bparg->value[row]/Bparg->value[2] - Bparg->derivative[2][col]*gradBparg[pro].value[row]/Bparg->value[2] - Bparg->derivative[row][col]*gradBparg[pro].value[2]/Bparg->value[2] + 2.0*Bparg->derivative[2][col]*Bparg->value[row]*gradBparg[pro].value[2]/(Bparg->value[2]*Bparg->value[2]) ) ) ;
	    	}
	    	adjfunc[pro] -= ( Xparg->tangent[col][dep] * extraterm ); 
	    }
	    adjfunc[pro] += ( lambdaarg->loc[col] * Xparg->loc[0]*(- gradBparg[pro].value[col] + gradBparg[pro].value[2]* Bparg->value[col] / Bparg->value[2] ) / Bparg->value[2] );
   	}
    }
    return adjfunc;
}

void RK4step(struct position *Xp, double varphi, double dvarphi, struct fieldparams allparams, struct field *Bfield_saved) {
	int row;
	struct position *dXp=malloc(4*sizeof(struct position)), Xpold;
	Xpold = *Xp;
	Bfield_saved[0] = Bfield(Xp->loc, varphi, allparams);
	dXp[0] = derivative_center(Bfield_saved, Xp);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp);
	Bfield_saved[1] = Bfield(Xp->loc, varphi + 0.5*dvarphi, allparams);
	dXp[1] = derivative_center(Bfield_saved+1, Xp);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp+1);
	Bfield_saved[2] = Bfield(Xp->loc, varphi + 0.5*dvarphi, allparams);	
	dXp[2] = derivative_center(Bfield_saved+2, Xp);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp+2);
	Bfield_saved[3] = Bfield(Xp->loc, varphi+dvarphi, allparams);
	dXp[3] = derivative_center(Bfield_saved+3, Xp);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp+1);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp+2);
	addstructsreassign(1.0, Xp,     (1.0/6.0)*dvarphi, dXp+3);
	for (row=0; row<4; row++) {
		free(dXp[row].tangent[0]); free(dXp[row].tangent[1]); 
		free(dXp[row].tangent);
	}
	free(dXp);
	return;
}

void RK4stepsave(struct position *Xp_adv, struct position *Xp, double varphi, double dvarphi, struct fieldparams allparams, struct field *Bfield_saved) {
	int row;
	struct position *dXp=malloc(4*sizeof(struct position));
	Xp[0] = *Xp_adv;
	Bfield_saved[0] = Bfield(Xp[0].loc, varphi, allparams);
	dXp[0] = derivative_center(Bfield_saved, Xp);
	Xp[1] = addstructs(1.0, Xp, 0.5*dvarphi, dXp);
	Bfield_saved[1] = Bfield(Xp[1].loc, varphi + 0.5*dvarphi, allparams);
	dXp[1] = derivative_center(Bfield_saved+1, Xp+1);
	Xp[2] = addstructs(1.0, Xp, 0.5*dvarphi, dXp+1);
	Bfield_saved[2] = Bfield(Xp[2].loc, varphi + 0.5*dvarphi, allparams);	
	dXp[2] = derivative_center(Bfield_saved+2, Xp+2);
	Xp[3] = addstructs(1.0, Xp, dvarphi, dXp+2);
	Bfield_saved[3] = Bfield(Xp[3].loc, varphi+dvarphi, allparams);
	dXp[3] = derivative_center(Bfield_saved+3, Xp+3);
	*Xp_adv = addstructs(1.0, Xp, (1.0/6.0)*dvarphi, dXp);
	addstructsreassign(1.0, Xp_adv,     (1.0/3.0)*dvarphi, dXp+1);
	addstructsreassign(1.0, Xp_adv,     (1.0/3.0)*dvarphi, dXp+2);
	addstructsreassign(1.0, Xp_adv,     (1.0/6.0)*dvarphi, dXp+3);
	for (row=0; row<4; row++) {
		free(dXp[row].tangent[0]); free(dXp[row].tangent[1]); 
		free(dXp[row].tangent);
	}
	free(dXp);
	return;
}

void RK4step_withsavedB(struct position *Xp, double varphi, double dvarphi, struct field *Bfield_saved)
{
	int row;
	struct position *dXp=malloc(4*sizeof(struct position)), Xpold;
	Xpold = *Xp;
	//printf("magnetic field before 1st segment of RK iteration: Bfield=%f\n", Bfield_saved[0].value[0]);
	dXp[0] = derivative_center(Bfield_saved, Xp);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp);
	//printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	dXp[1] = derivative_center(Bfield_saved+1, Xp);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp+1);
	//printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	dXp[2] = derivative_center(Bfield_saved+2, Xp);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp+2);
	//printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	dXp[3] = derivative_center(Bfield_saved+3, Xp);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp+1);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp+2);
	addstructsreassign(1.0, Xp,     (1.0/6.0)*dvarphi, dXp+3);
	//printf("after 4th segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	for (row=0; row<4; row++) {
		free(dXp[row].tangent[0]); free(dXp[row].tangent[1]); 
		free(dXp[row].tangent);
	}
	free(dXp);
	return;
}

void RK4step_lambdacirc_mutangent(struct position *Xp, struct position *adjvariable, double varphi, double dvarphi, struct field *Bfield_saved) {
	int row;
	struct position *dXp=malloc(4*sizeof(struct position)), *dadjvariable=malloc(4*sizeof(struct position)), Xpold, adjvariableold;
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
	dXp[0] = derivative_center(Bfield_saved, Xp);
	dadjvariable[0] = derivative_lambdacirc_mutangent(Bfield_saved, Xp, adjvariable);
	// I think this memory needs to be kept (because it corresponds to Xp
	// before the next iteration, which is preserved in the function that
	// calls this function. So the 2 lines below must remain commented out
	//free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	//free(adjvariable->tangent[0]); free(adjvariable->tangent[1]); free(adjvariable->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp);
	*adjvariable = addstructs(1.0, &adjvariableold, 0.5*dvarphi, dadjvariable);
	//printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	dXp[1] = derivative_center(Bfield_saved+1, Xp);
	dadjvariable[1] = derivative_lambdacirc_mutangent(Bfield_saved+1, Xp, adjvariable);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(adjvariable->tangent[0]); free(adjvariable->tangent[1]); free(adjvariable->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp+1);
	*adjvariable = addstructs(1.0, &adjvariableold, 0.5*dvarphi, dadjvariable+1);
	//printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	dXp[2] = derivative_center(Bfield_saved+2, Xp);
	dadjvariable[2] = derivative_lambdacirc_mutangent(Bfield_saved+2, Xp, adjvariable);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(adjvariable->tangent[0]); free(adjvariable->tangent[1]); free(adjvariable->tangent);
	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp+2);
	*adjvariable = addstructs(1.0, &adjvariableold, dvarphi, dadjvariable+2);
	//printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	dXp[3] = derivative_center(Bfield_saved+3, Xp);
	dadjvariable[3] = derivative_lambdacirc_mutangent(Bfield_saved+3, Xp, adjvariable);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(adjvariable->tangent[0]); free(adjvariable->tangent[1]); free(adjvariable->tangent);
	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp+1);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp+2);
	addstructsreassign(1.0, Xp,     (1.0/6.0)*dvarphi, dXp+3);
	*adjvariable = addstructs(1.0, &adjvariableold, (1.0/6.0)*dvarphi, dadjvariable);
	addstructsreassign(1.0, adjvariable,     (1.0/3.0)*dvarphi, dadjvariable+1);
	addstructsreassign(1.0, adjvariable,     (1.0/3.0)*dvarphi, dadjvariable+2);
	addstructsreassign(1.0, adjvariable,     (1.0/6.0)*dvarphi, dadjvariable+3);
	for (row=0; row<4; row++) {
		free(dXp[row].tangent[0]); free(dXp[row].tangent[1]); 
		free(dadjvariable[row].tangent[0]); free(dadjvariable[row].tangent[1]); 
		free(dXp[row].tangent);
		free(dadjvariable[row].tangent);
	}
	free(dXp);
	free(dadjvariable);
	return;
}

void RK4step_lambdatangent(struct position *Xp, struct position *lambda, struct position *sperp, struct position *mu, double varphi, double dvarphi, struct field *Bfield_saved) {
	int row;
	struct position *dXp=malloc(4*sizeof(struct position)), *dlambda=malloc(4*sizeof(struct position));
	struct position *dsperp= malloc(4*sizeof(struct position)), *dmu=malloc(4*sizeof(struct position));
	struct position Xpold, lambdaold, sperpold, muold;
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

	dXp[0] = derivative_center(Bfield_saved, Xp);
	dlambda[0] = derivative_lambdatangent(Bfield_saved, Xp, lambda, sperp, mu);
	dsperp[0] = derivative_sperp(Bfield_saved, Xp, sperp);
	dmu[0] = derivative_lambdacirc_mutangent(Bfield_saved, Xp, mu);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	free(sperp->tangent[0]); free(sperp->tangent[1]); free(sperp->tangent);
	//free(mu->tangent[0]); free(mu->tangent[1]); free(mu->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp);
	*lambda = addstructs(1.0, &lambdaold, 0.5*dvarphi, dlambda);
	*sperp = addstructs(1.0, &sperpold, 0.5*dvarphi, dsperp);
	*mu = addstructs(1.0, &muold, 0.5*dvarphi, dmu);
	//printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	dXp[1] = derivative_center(Bfield_saved+1, Xp);
	dlambda[1] = derivative_lambdatangent(Bfield_saved+1, Xp, lambda, sperp, mu);
	dsperp[1] = derivative_sperp(Bfield_saved+1, Xp, sperp);
	dmu[1] = derivative_lambdacirc_mutangent(Bfield_saved+1, Xp, mu);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	free(sperp->tangent[0]); free(sperp->tangent[1]); free(sperp->tangent);
	free(mu->tangent[0]); free(mu->tangent[1]); free(mu->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp+1);
	*lambda = addstructs(1.0, &lambdaold, 0.5*dvarphi, dlambda+1);
	*sperp = addstructs(1.0, &sperpold, 0.5*dvarphi, dsperp+1);
	*mu = addstructs(1.0, &muold, 0.5*dvarphi, dmu+1);
	//printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	dXp[2] = derivative_center(Bfield_saved+2, Xp);
	dlambda[2] = derivative_lambdatangent(Bfield_saved+2, Xp, lambda, sperp, mu);
	dsperp[2] = derivative_sperp(Bfield_saved+2, Xp, sperp);
	dmu[2] = derivative_lambdacirc_mutangent(Bfield_saved+2, Xp, mu);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	free(sperp->tangent[0]); free(sperp->tangent[1]); free(sperp->tangent);
	free(mu->tangent[0]); free(mu->tangent[1]); free(mu->tangent);
	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp+2);
	*lambda = addstructs(1.0, &lambdaold, dvarphi, dlambda+2);
	*sperp = addstructs(1.0, &sperpold, dvarphi, dsperp+2);
	*mu = addstructs(1.0, &muold, dvarphi, dmu+2);
	//printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	dXp[3] = derivative_center(Bfield_saved+3, Xp);
	dlambda[3] = derivative_lambdatangent(Bfield_saved+3, Xp, lambda, sperp, mu);
	dsperp[3] = derivative_sperp(Bfield_saved+3, Xp, sperp);
	dmu[3] = derivative_lambdacirc_mutangent(Bfield_saved+3, Xp, mu);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	free(sperp->tangent[0]); free(sperp->tangent[1]); free(sperp->tangent);
	free(mu->tangent[0]); free(mu->tangent[1]); free(mu->tangent);
	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp+1);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp+2);
	addstructsreassign(1.0, Xp,     (1.0/6.0)*dvarphi, dXp+3);
	*lambda = addstructs(1.0, &lambdaold, (1.0/6.0)*dvarphi, dlambda);
	addstructsreassign(1.0, lambda,     (1.0/3.0)*dvarphi, dlambda+1);
	addstructsreassign(1.0, lambda,     (1.0/3.0)*dvarphi, dlambda+2);
	addstructsreassign(1.0, lambda,     (1.0/6.0)*dvarphi, dlambda+3);
	*sperp = addstructs(1.0, &sperpold, (1.0/6.0)*dvarphi, dsperp);
	addstructsreassign(1.0, sperp,     (1.0/3.0)*dvarphi, dsperp+1);
	addstructsreassign(1.0, sperp,     (1.0/3.0)*dvarphi, dsperp+2);
	addstructsreassign(1.0, sperp,     (1.0/6.0)*dvarphi, dsperp+3);
	*mu = addstructs(1.0, &muold, (1.0/6.0)*dvarphi, dmu);
	addstructsreassign(1.0, mu,     (1.0/3.0)*dvarphi, dmu+1);
	addstructsreassign(1.0, mu,     (1.0/3.0)*dvarphi, dmu+2);
	addstructsreassign(1.0, mu,     (1.0/6.0)*dvarphi, dmu+3);
	for (row=0; row<4; row++) {
		free(dXp[row].tangent[0]); free(dXp[row].tangent[1]); 
		free(dlambda[row].tangent[0]); free(dlambda[row].tangent[1]); 
		free(dsperp[row].tangent[0]); free(dsperp[row].tangent[1]); 
		free(dmu[row].tangent[0]); free(dmu[row].tangent[1]); 
		free(dXp[row].tangent); 
		free(dlambda[row].tangent);
		free(dsperp[row].tangent);
		free(dmu[row].tangent);
	}
	free(dXp);
	free(dlambda);
	free(dsperp);
	free(dmu);
	return;
}

void RK4step_lambdaRes(struct position *lambdamu_out, struct position *lambdamu, struct position *Xp, double varphi, double dvarphi, struct field *Bfield_saved)
{
	int row;
	struct position *dlambdamu=malloc(4*sizeof(struct position)), lambdamuold;
	lambdamuold.tangent = set_identity();
	for (row = 0; row < 2; row++) {
		lambdamuold.loc[row] = lambdamu_out->loc[row];
		lambdamuold.tangent[row][0] = lambdamu_out->tangent[row][0];
		lambdamuold.tangent[row][1] = lambdamu_out->tangent[row][1];
	}
	lambdamu[0] = *lambdamu_out;
	dlambdamu[0] = derivative_lambdamu(Bfield_saved, Xp, lambdamu);
	lambdamu[1] = addstructs(1.0, lambdamu_out, 0.5*dvarphi, dlambdamu);
	//printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	dlambdamu[1] = derivative_lambdamu(Bfield_saved+1, Xp+1, lambdamu+1);
	lambdamu[2] = addstructs(1.0, lambdamu_out, 0.5*dvarphi, dlambdamu+1);
	//printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	dlambdamu[2] = derivative_lambdamu(Bfield_saved+2, Xp+2, lambdamu+2);
	lambdamu[3] = addstructs(1.0, lambdamu_out, dvarphi, dlambdamu+2);
	//printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	dlambdamu[3] = derivative_lambdamu(Bfield_saved+3, Xp+3, lambdamu+3);
	*lambdamu_out = addstructs(1.0, lambdamu_out, (1.0/6.0)*dvarphi, dlambdamu);
	addstructsreassign(1.0, lambdamu_out,     (1.0/3.0)*dvarphi, dlambdamu+1);
	addstructsreassign(1.0, lambdamu_out,     (1.0/3.0)*dvarphi, dlambdamu+2);
	addstructsreassign(1.0, lambdamu_out,     (1.0/6.0)*dvarphi, dlambdamu+3);
	for (row=0; row<4; row++) {
		free(dlambdamu[row].tangent[0]); free(dlambdamu[row].tangent[1]); 
		free(dlambdamu[row].tangent);
	}
	free(dlambdamu);
	return;
}

void RK4step_gradcirc(double *number, struct position *Xp, struct position *lambda, double varphi, double dvarphi, struct field *Bfield_saved, struct fieldparams allparams, int diffparam_ind1, int diffparam_ind2)
{
	int row, num_params = allparams.n_diff;
	struct field *Bpointgrad;
	struct position *dXp=malloc(4*sizeof(struct position)), Xpold;
	struct position *dlambda=malloc(4*sizeof(struct position)), lambdaold;
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
	Bpointgrad = gradBfield(Xp->loc, varphi, allparams, diffparam_ind1, diffparam_ind2);
	dXp[0] = derivative_center(Bfield_saved, Xp);
	//printf("dXp[0]=%f\n", dXp[0]->loc[0]);
	dlambda[0] = derivative_lambdacirc_mutangent(Bfield_saved, Xp, lambda);
	dnumber[0] = derivative_gradcirc(Bfield_saved, Bpointgrad, Xp, lambda, num_params);
	free(Bpointgrad);
	//free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	//free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	//printf("1st segment of RK iteratiion: BR=%f, GXBR=%f, dR=%f, dlambdaR=%f, dnumberX=%f\n", Bpoint->value[0], Bpointgrad[0].value[0], dXp[0]->loc[0], dlambda[0]->loc[0], dnumber[0][0]);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp);
	*lambda = addstructs(1.0, &lambdaold, 0.5*dvarphi, dlambda);
	//printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpointgrad = gradBfield(Xp->loc, varphi+0.5*dvarphi, allparams, diffparam_ind1, diffparam_ind2);
	dXp[1] = derivative_center(Bfield_saved+1, Xp);
	dlambda[1] = derivative_lambdacirc_mutangent(Bfield_saved+1, Xp, lambda);
	dnumber[1] = derivative_gradcirc(Bfield_saved+1, Bpointgrad, Xp, lambda, num_params);
	free(Bpointgrad);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp+1);
	*lambda = addstructs(1.0, &lambdaold, 0.5*dvarphi, dlambda+1);
	//printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpointgrad = gradBfield(Xp->loc, varphi+0.5*dvarphi, allparams, diffparam_ind1, diffparam_ind2);
	dXp[2] = derivative_center(Bfield_saved+2, Xp);
	dlambda[2] = derivative_lambdacirc_mutangent(Bfield_saved+2, Xp, lambda);
	dnumber[2] = derivative_gradcirc(Bfield_saved+2, Bpointgrad, Xp, lambda, num_params);
	free(Bpointgrad);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp+2);
	*lambda = addstructs(1.0, &lambdaold, dvarphi, dlambda+2);
	//printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpointgrad = gradBfield(Xp->loc, varphi+dvarphi, allparams, diffparam_ind1, diffparam_ind2);
	dXp[3] = derivative_center(Bfield_saved+3, Xp);
	dlambda[3] = derivative_lambdacirc_mutangent(Bfield_saved+3, Xp, lambda);
	dnumber[3] = derivative_gradcirc(Bfield_saved+3, Bpointgrad, Xp, lambda, num_params);
	free(Bpointgrad);
	free(Xp->tangent[0]); free(Xp->tangent[1]); free(Xp->tangent);
	free(lambda->tangent[0]); free(lambda->tangent[1]); free(lambda->tangent);
	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp+1);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp+2);
	addstructsreassign(1.0, Xp,     (1.0/6.0)*dvarphi, dXp+3);
	*lambda = addstructs(1.0, &lambdaold, (1.0/6.0)*dvarphi, dlambda);
	addstructsreassign(1.0, lambda,     (1.0/3.0)*dvarphi, dlambda+1);
	addstructsreassign(1.0, lambda,     (1.0/3.0)*dvarphi, dlambda+2);
	addstructsreassign(1.0, lambda,     (1.0/6.0)*dvarphi, dlambda+3);
	for (row=0;row<num_params;row++) {
		number[row] = numberold[row] + dvarphi*(dnumber[0][row]/6.0 + dnumber[1][row]/3.0 + dnumber[2][row]/3.0 + dnumber[3][row]/6.0);
	}

	for (row=0; row<4; row++) {
		free(dXp[row].tangent[0]); free(dXp[row].tangent[1]); 
		free(dXp[row].tangent); 
		free(dlambda[row].tangent[0]); free(dlambda[row].tangent[1]); 
		free(dlambda[row].tangent);
		free(dnumber[row]);
	}
	free(dXp);
	free(dlambda);
	return;
}

void RK4step_gradtangent(double *number, struct position *Xp, struct position *lambda, struct position *sperp, struct position *mu, double varphi, double dvarphi, struct field *Bfield_saved, struct fieldparams allparams, int diffparam_ind1, int diffparam_ind2) {
	int row, num_params=allparams.n_diff;
	struct field *Bpoint, *Bpointgrad;
	struct position *dXp=malloc(4*sizeof(struct position)), *dlambda=malloc(4*sizeof(struct position)), *dsperp=malloc(4*sizeof(struct position)), *dmu=malloc(4*sizeof(struct position)), Xpold, lambdaold, sperpold, muold;
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
	Bpoint = Bfield_saved;
	Bpointgrad = gradBfield(Xp->loc, varphi, allparams, diffparam_ind1, diffparam_ind2);
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
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp);
	*lambda = addstructs(1.0, &lambdaold, 0.5*dvarphi, dlambda);
	*sperp = addstructs(1.0, &sperpold, 0.5*dvarphi, dsperp);
	*mu = addstructs(1.0, &muold, 0.5*dvarphi, dmu);
	//printf("after 1st segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpoint = Bfield_saved+1;
	Bpointgrad = gradBfield(Xp->loc, varphi+0.5*dvarphi, allparams, diffparam_ind1, diffparam_ind2);
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
	*Xp = addstructs(1.0, &Xpold, 0.5*dvarphi, dXp+1);
	*lambda = addstructs(1.0, &lambdaold, 0.5*dvarphi, dlambda+1);
	*sperp = addstructs(1.0, &sperpold, 0.5*dvarphi, dsperp+1);
	*mu = addstructs(1.0, &muold, 0.5*dvarphi, dmu+1);
	//printf("after 2nd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpoint = Bfield_saved+2;
	Bpointgrad = gradBfield(Xp->loc, varphi+0.5*dvarphi, allparams, diffparam_ind1, diffparam_ind2);
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
	*Xp = addstructs(1.0, &Xpold, dvarphi, dXp+2);
	*lambda = addstructs(1.0, &lambdaold, dvarphi, dlambda+2);
	*sperp = addstructs(1.0, &sperpold, dvarphi, dsperp+2);
	*mu = addstructs(1.0, &muold, dvarphi, dmu+2);
	//printf("after 3rd segment of RK iteratiion: Xp->loc[0]=%f, Xp1->loc[1]=%f\n", Xp->loc[0], Xp->loc[1]);
	Bpoint = Bfield_saved+3;
	Bpointgrad = gradBfield(Xp->loc, varphi+dvarphi, allparams, diffparam_ind1, diffparam_ind2);
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
	*Xp = addstructs(1.0, &Xpold, (1.0/6.0)*dvarphi, dXp);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp+1);
	addstructsreassign(1.0, Xp,     (1.0/3.0)*dvarphi, dXp+2);
	addstructsreassign(1.0, Xp,     (1.0/6.0)*dvarphi, dXp+3);
	*lambda = addstructs(1.0, &lambdaold, (1.0/6.0)*dvarphi, dlambda);
	addstructsreassign(1.0, lambda,     (1.0/3.0)*dvarphi, dlambda+1);
	addstructsreassign(1.0, lambda,     (1.0/3.0)*dvarphi, dlambda+2);
	addstructsreassign(1.0, lambda,     (1.0/6.0)*dvarphi, dlambda+3);
	*sperp = addstructs(1.0, &sperpold, (1.0/6.0)*dvarphi, dsperp);
	addstructsreassign(1.0, sperp,     (1.0/3.0)*dvarphi, dsperp+1);
	addstructsreassign(1.0, sperp,     (1.0/3.0)*dvarphi, dsperp+2);
	addstructsreassign(1.0, sperp,     (1.0/6.0)*dvarphi, dsperp+3);
	*mu = addstructs(1.0, &muold, (1.0/6.0)*dvarphi, dmu);
	addstructsreassign(1.0, mu,     (1.0/3.0)*dvarphi, dmu+1);
	addstructsreassign(1.0, mu,     (1.0/3.0)*dvarphi, dmu+2);
	addstructsreassign(1.0, mu,     (1.0/6.0)*dvarphi, dmu+3);
	for (row=0;row<num_params;row++) {
		number[row] = numberold[row] + dvarphi*(dnumber[0][row]/6.0 + dnumber[1][row]/3.0 + dnumber[2][row]/3.0 + dnumber[3][row]/6.0);
	}
	for (row=0; row<4; row++) {
		free(dXp[row].tangent[0]); free(dXp[row].tangent[1]); 
		free(dlambda[row].tangent[0]); free(dlambda[row].tangent[1]); 
		free(dsperp[row].tangent[0]); free(dsperp[row].tangent[1]); 
		free(dmu[row].tangent[0]); free(dmu[row].tangent[1]); 
		free(dXp[row].tangent);
		free(dlambda[row].tangent);
		free(dsperp[row].tangent);
		free(dmu[row].tangent);
		free(dnumber[row]);
	}
	free(dXp);
	free(dlambda);
	free(dsperp);
	free(dmu);
	return;
}

void RK4step_gradRes(double *number, struct position *Xp, struct position *lambdamu, double varphi, double dvarphi, struct field *Bfield_saved, struct fieldparams allparams, int diffparam_ind1, int diffparam_ind2) {
   int row, num_params=allparams.n_diff;
   struct field **Bpointgrad = malloc(4*sizeof(struct field));
   double *dnumber[4], numberold[num_params];
   for (row = 0; row < num_params; row++) 
      numberold[row] = number[row];
   Bpointgrad[0] = gradBfield(Xp[0].loc, varphi, allparams, diffparam_ind1, diffparam_ind2);
   dnumber[0] = derivative_gradRes(Bfield_saved, Bpointgrad[0], Xp, lambdamu, num_params);
   Bpointgrad[1] = gradBfield(Xp[1].loc, varphi+0.5*dvarphi, allparams, diffparam_ind1, diffparam_ind2);
   dnumber[1] = derivative_gradRes(Bfield_saved+1, Bpointgrad[1], Xp+1, lambdamu+1, num_params);
   Bpointgrad[2] = gradBfield(Xp[2].loc, varphi+0.5*dvarphi, allparams, diffparam_ind1, diffparam_ind2);
   dnumber[2] = derivative_gradRes(Bfield_saved+2, Bpointgrad[2], Xp+2, lambdamu+2, num_params);
   Bpointgrad[3] = gradBfield(Xp[3].loc, varphi+dvarphi, allparams, diffparam_ind1, diffparam_ind2);
   dnumber[3] = derivative_gradRes(Bfield_saved+3, Bpointgrad[3], Xp+3, lambdamu+3, num_params);
   for (row=0;row<num_params;row++) 
      number[row] = numberold[row] + dvarphi*(dnumber[0][row]/6.0 + dnumber[1][row]/3.0 + dnumber[2][row]/3.0 + dnumber[3][row]/6.0);
   return;
}

void RK4step_gradXp(double **number, struct position *Xp, struct position *lambdaXp, double varphi, double dvarphi, struct field *Bfield_saved, struct fieldparams allparams, int diffparam_ind1, int diffparam_ind2) {
   int row, num_params=allparams.n_diff;
   struct field **Bpointgrad = malloc(4*sizeof(struct field));
   double **dnumber[4];
   double numberold[2][num_params];
   for (row = 0; row < num_params; row++) {
      numberold[0][row] = number[0][row];
      numberold[1][row] = number[1][row];
   }
   Bpointgrad[0] = gradBfield(Xp[0].loc, varphi, allparams, diffparam_ind1, diffparam_ind2);
   dnumber[0] = derivative_gradXp(Bfield_saved, Bpointgrad[0], Xp, lambdaXp, num_params);
   Bpointgrad[1] = gradBfield(Xp[1].loc, varphi+0.5*dvarphi, allparams, diffparam_ind1, diffparam_ind2);
   dnumber[1] = derivative_gradXp(Bfield_saved+1, Bpointgrad[1], Xp+1, lambdaXp+1, num_params);
   Bpointgrad[2] = gradBfield(Xp[2].loc, varphi+0.5*dvarphi, allparams, diffparam_ind1, diffparam_ind2);
   dnumber[2] = derivative_gradXp(Bfield_saved+2, Bpointgrad[2], Xp+2, lambdaXp+2, num_params);
   Bpointgrad[3] = gradBfield(Xp[3].loc, varphi+dvarphi, allparams, diffparam_ind1, diffparam_ind2);
   dnumber[3] = derivative_gradXp(Bfield_saved+3, Bpointgrad[3], Xp+3, lambdaXp+3, num_params);
   for (row=0;row<num_params;row++) {
      number[0][row] = numberold[0][row] + dvarphi*(dnumber[0][0][row]/6.0 + dnumber[1][0][row]/3.0 + dnumber[2][0][row]/3.0 + dnumber[3][0][row]/6.0);
      number[1][row] = numberold[1][row] + dvarphi*(dnumber[0][1][row]/6.0 + dnumber[1][1][row]/3.0 + dnumber[2][1][row]/3.0 + dnumber[3][1][row]/6.0);
   }
   return;
}
