//Currently under construction
//Author: Alessandro Geraldini
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "isc.h"
#include <gsl/gsl_sf_bessel.h>

int main()
{
	clock_t start = clock();
	int N_gridphi_fieldperiod, N_gridphi_tor, num_turns=1, m0, optimizeheli = 0;
	int *n_coils=NULL, **n_consts, qq=1, *qq_segs=NULL, qq_segsq=1, pol_mode, tor_mode, ind, row, col, row_largest, count;
	int find_axis = 1, deltaindmax = 20, deltaind, rowallowed[2][3];
	int L_fixedpoints, N_line, seg_step = 20;
	int coil_ind, seg_ind, max_coil_ind, num_segs;
	int sizeline = 100;
	int q0val = 0, q02val = 0, *q0_index=NULL, *q0_index2=NULL; // It is important to initialize the integer stored in q0_index to 0
	struct position *Xp=calloc(1,sizeof(struct position)), *axis, *island_center, **mu, **lambdaQ, lambda, **lambdamu_saved;
	struct position *island_center2, **island_center_saved2, **lambdamu_saved2;
	struct position *Xp2=calloc(1,sizeof(struct position));
	struct position *island_center_FD[2], **island_center_saved, **axis_saved;
	struct ext_position *ext_center, *ext_center2, *ext_center_FD[2];
	struct field **Bfield_island, **Bfield_island2, **Bfield_axis, Bpoint;
	struct fieldparams allparams;
	//struct field *BB, *gradBB; for field checks
	double ***coils=NULL, *shapeRes, *shapecirc, **shapetan, **shapewidth, **realshapetan, **realshapewidth, *realshapecirc, *realshapeRes, *Sigma;  
	double *shapeRes2, Res2;
	double Res, shapeRes_largest=0.0;// **adjXp;
	double powdeltain = -1.0, steppowdelta = -3.0, *doublepointer, **n_params=NULL;
	double phichec = 0.0, factor, dldown, dl, testpos[2] = {0.95, 0.0};
	struct field realBpointgrad[3];
	double **checkshapewidth, **checkshapeSigma, *checkshapecirc, *checkshapeRes, savedcoilseg;
	double dr[3], drdown[3], *checkparam = malloc(4*sizeof(double));
	double *iota=NULL, *rmin=NULL, *zmin=NULL, axis_phi0[2];
	double **evec=malloc(2*sizeof(double*)), *eval=malloc(2*sizeof(double)), trace, det, iota_axis;
	double delta, checkBpointgrad[3], checknablaBpointgrad[3][2];
	double normgradwidth[3], normgradcirc[3], normgradSigma[3], normgradRes[3], checknormgradRes[3][deltaindmax], checknormgradwidth[3][deltaindmax], checknormgradcirc[3][deltaindmax], checknormgradSigma[3][deltaindmax];
	struct field *Bpointgrad, Bpointdelta[3][2];
	char line[sizeline];
	char poiniota, findisla, calcshap, checshap, checfast, checkmag;
	FILE *OPTIONS;

	q0_index = &q0val;
	q0_index2 = &q02val;

	/* 
	Specify the options of the code:
	Poincare = ('y', 'n') --> want to produce a Poincare plot & iota profile?
	findisla = ('y', 'n') --> want to find magnetic islands?
	calcshap = ('y', 'n') --> want to calculate shape gradient of magnetic islands?
	checshap = ('y', 'n') --> want to check shape gradient using finite differences?
	*/
	if ( (OPTIONS = fopen("options.txt", "r")) == NULL ) {
		printf("Couldn't open options.txt\n");	
		exit(1);
	}
	while (fgets(line, sizeline, OPTIONS) != NULL) {
		if (strncmp("checkmag", line, 8) == 0)  checkmag = line[9];
		if (strncmp("Poincare", line, 8) == 0)  poiniota = line[9];
		if (strncmp("findisla", line, 8) == 0)  findisla = line[9];
		if (strncmp("calcshap", line, 8) == 0)  calcshap = line[9];
		if (strncmp("checshap", line, 8) == 0)  checshap = line[9];
		if (strncmp("checfast", line, 8) == 0)  checfast = line[9];
	}
	fclose(OPTIONS);
	printf("checkmag = %c\npoiniota = %c\nfindisla = %c\ncalcshap = %c\nchecshap = %c\nchecfast = %c\n", checkmag, poiniota, findisla, calcshap, checshap, checfast);

	/* 
	Read parameters of magnetic field:
	coils contains the magnetic field parameters
	for magnetic fields produced by coils these are coil segment positions
	for Reiman configurations coils contains the parameters 
	(iota0, iota1, epsilon_i) for i =(0, ..., N)
	*/
	n_coils = &qq; qq_segs = &qq_segsq; n_consts = &qq_segs; // Initializes pointers to parameters
	n_params = &doublepointer;
	allparams = fetchparams();
	m0 = allparams.m0_fieldperiods;
	//coils = fieldparams(type, &m0, n_coils, n_consts, n_params);
	N_gridphi_fieldperiod = 120/m0;

	for (ind=0; ind < *n_coils; ind++) {
		printf("%d ", n_consts[0][ind]);
	}
	printf("\n");
	//printf("coils[0][0][0] = %f", coils[0][0][0]);
	//printf("n_coils = %d\n", *n_coils);
	//printf("n_consts = %d\n", **n_consts);
	if (strncmp(allparams.type, "coil", 4) == 0) {
		pol_mode = 4; tor_mode = 1;
		delta = 0.005;
	}
	else if (strncmp(allparams.type, "heli", 4) == 0) {
		pol_mode = 8; tor_mode = 5;
		pol_mode = 3; tor_mode = 5; //for heli2
		delta = 0.00001;
	}
	else if (strncmp(allparams.type, "Reim", 4) == 0) {
		pol_mode = 6; tor_mode = 1;
		delta = 0.001;
	}
	N_gridphi_tor = N_gridphi_fieldperiod*m0;

	evec[0]=malloc(2*sizeof(double)); evec[1]=malloc(2*sizeof(double));
	/* make arrays of coils */
	printf("field periods = %d\n", m0);
	printf("N points per field period = %d\n", N_gridphi_fieldperiod);
	clock_t tcoils = clock();
	printf("... %f seconds ...\n", (double) (tcoils-start)/CLOCKS_PER_SEC);
	//Xp->tangent = malloc(2*sizeof(double*));
	//Xp->tangent[0] = malloc(2*sizeof(double)); Xp->tangent[1] = malloc(2*sizeof(double));
	Xp->tangent = set_identity();
	Xp2->tangent = set_identity();
	//printf("n_coils = %d, n_consts = %d\n", *n_coils, **n_consts);
	Bpoint = Bfield(testpos, 0.0, allparams);
	printf("(BR, BZ, Bphi) = (%f, %f, %f)\n", Bpoint.value[0], Bpoint.value[1], Bpoint.value[2]);

	coils = allparams.diffparams;
	printf("The magnetic field type is %s\n", allparams.type);
	if (checkmag == 'y') {
		if (strncmp("coil", allparams.type, 5) == 0 ) {
			Xp->loc[0] = 1.53; Xp->loc[1]= 0.0; //NCSX
			coil_ind = 5;
			seg_ind = 14;
			checkparam[0] = coils[coil_ind][seg_ind][0];
			checkparam[1] = coils[coil_ind][seg_ind][1];
			checkparam[2] = coils[coil_ind][seg_ind][2];
			checkparam[3] = coils[coil_ind][seg_ind][3];
			dr[0] = coils[coil_ind][seg_ind+1][0] - coils[coil_ind][seg_ind][0];
			dr[1] = coils[coil_ind][seg_ind+1][1] - coils[coil_ind][seg_ind][1];
			dr[2] = coils[coil_ind][seg_ind+1][2] - coils[coil_ind][seg_ind][2];
			//drlength = pow(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2], 0.5);
			Bpointgrad = gradBfield(Xp->loc, phichec, allparams, coil_ind, seg_ind);
			for (row=0;row<3;row++) {
				for (ind=0;ind<3; ind++) {
					realBpointgrad[ind].value[row] = (dr[(ind+1)%3]*Bpointgrad[(ind+2)%3].value[row] - dr[(ind+2)%3]*Bpointgrad[(ind+1)%3].value[row]); 
					/*
					to compare with more precisely with finite differences (below)
					it is necessary to change this slightly.
					work in progress
					*/
				}
			}
			coils[coil_ind][seg_ind][0] += delta; 
			Bpointdelta[0][1] = Bfield(Xp->loc, phichec, allparams); 
			coils[coil_ind][seg_ind][0] -= 2.0*delta; 
			Bpointdelta[0][0] = Bfield(Xp->loc, phichec, allparams); 
			coils[coil_ind][seg_ind][0] += delta; 
			coils[coil_ind][seg_ind][1] += delta; 
			Bpointdelta[1][1] = Bfield(Xp->loc, phichec, allparams); 
			coils[coil_ind][seg_ind][1] -= 2.0*delta; 
			Bpointdelta[1][0] = Bfield(Xp->loc, phichec, allparams); 
			coils[coil_ind][seg_ind][1] += delta; 
			coils[coil_ind][seg_ind][2] += delta; 
			Bpointdelta[2][1] = Bfield(Xp->loc, phichec, allparams); 
			coils[coil_ind][seg_ind][2] -= 2.0*delta; 
			Bpointdelta[2][0] = Bfield(Xp->loc, phichec, allparams); 
			coils[coil_ind][seg_ind][2] += delta; 
			printf("checking shape gradient\n");
			for (row = 0; row<3; row++) 
			checkBpointgrad[row] = ( Bpointdelta[0][1].value[row] - Bpointdelta[0][0].value[row] ) / (2.0*delta);
			printf("checkgradB = (%f, %f, %f) (direct), (%f, %f, %f) (FD) \n", realBpointgrad[0].value[0], realBpointgrad[0].value[1], realBpointgrad[0].value[2], checkBpointgrad[0], checkBpointgrad[1], checkBpointgrad[2]);
			for (row = 0; row<3; row++) 
			checkBpointgrad[row] = ( Bpointdelta[1][1].value[row] - Bpointdelta[1][0].value[row] ) / (2.0*delta);
			printf("checkgradB = (%f, %f, %f) (direct), (%f, %f, %f) (FD) \n", realBpointgrad[1].value[0], realBpointgrad[1].value[1], realBpointgrad[1].value[2], checkBpointgrad[0], checkBpointgrad[1], checkBpointgrad[2]);
			for (row = 0; row<3; row++) 
			checkBpointgrad[row] = ( Bpointdelta[2][1].value[row] - Bpointdelta[2][0].value[row] ) / (2.0*delta);
			printf("checkgradB = (%f, %f, %f) (direct), (%f, %f, %f) (FD) \n", realBpointgrad[2].value[0], realBpointgrad[2].value[1], realBpointgrad[2].value[2], checkBpointgrad[0], checkBpointgrad[1], checkBpointgrad[2]);
		}
		else if (strncmp(allparams.type, "Reim", 4) == 0) {
			Xp->loc[0] = 1.05; Xp->loc[1]= 0.0; //NCSX
			printf("Xp = (%f, %f)\n", Xp->loc[0], Xp->loc[1]);
			Bpoint = Bfield(Xp->loc, phichec, allparams); 
			printstructfield("Bpoint =", &Bpoint);
			//exit(0);
		}
		else if (strncmp("heli", allparams.type, 5) == 0 ) {
			Xp->loc[0] = 0.95; Xp->loc[1]= 0.0; 
			delta = 0.001;
			printf("Xp = (%f, %f)\n", Xp->loc[0], Xp->loc[1]);
			printf("n_coils[3] = %d\n", allparams.n_coils);
			Bpointgrad = gradBfield(Xp->loc, phichec, allparams, 0, 0);
			Bpoint = Bfield(Xp->loc, phichec, allparams); 
			for (row = 0; row<7; row++) {
				allparams.diffparams[0][0][row] += delta;
				Bpointdelta[1][1] = Bfield(Xp->loc, phichec, allparams); 
				allparams.diffparams[0][0][row] -= 2.0*delta;
				Bpointdelta[1][0] = Bfield(Xp->loc, phichec, allparams); 
				allparams.diffparams[0][0][row] += delta;
				for (col=0;col<3;col++) {
					checkBpointgrad[col] = ( Bpointdelta[1][1].value[col] - Bpointdelta[1][0].value[col] ) / (2.0*delta);
					for (ind=0; ind<2; ind++) {
						checknablaBpointgrad[col][ind] = ( Bpointdelta[1][1].derivative[col][ind] - Bpointdelta[1][0].derivative[col][ind] ) / (2.0*delta);
					}
					printf("checkgradnablaB = (%14.14f, %14.14f) (direct), (%14.14f, %14.14f) (FD) \n", Bpointgrad[row].derivative[col][0], Bpointgrad[row].derivative[col][1], checknablaBpointgrad[col][0], checknablaBpointgrad[col][1]);
				}
				printf("checkgradB = (%14.14f, %14.14f, %14.14f) (direct), (%14.14f, %14.14f, %14.14f) (FD) \n", Bpointgrad[row].value[0], Bpointgrad[row].value[1], Bpointgrad[row].value[2], checkBpointgrad[0], checkBpointgrad[1], checkBpointgrad[2]);
			}
			//Bpointdelta[2][1] = Bfield(Xp->loc, phichec, allparams); 
			//Bpointdelta[2][0] = Bfield(Xp->loc, phichec, allparams); 
		}
	}
	if (find_axis == 1) {
		if (strncmp(allparams.type, "coil", 4) == 0) {
			Xp->loc[0] = 1.55; Xp->loc[1]= 0.0; //NCSX
		}
		else if (strncmp(allparams.type, "heli", 4) == 0) {
			Xp->loc[0] = 0.97; Xp->loc[1]= 0.0; //helical
			Xp->loc[0] = 0.9834425342; Xp->loc[1]= 0.0; //helical2
		}
		else {
			Xp->loc[0] = 1.1; Xp->loc[1]= 0.0;
		}
		Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		Bfield_axis = malloc(N_gridphi_fieldperiod*sizeof(struct field));
		for (ind=0; ind<N_gridphi_fieldperiod; ind++) {
			Bfield_axis[ind] = malloc(4*sizeof(struct field));
		}
		printf("Entering function that solves for the magnetic axis-->\n");
		axis_saved = malloc(N_gridphi_fieldperiod*sizeof(struct position));
		for (ind=0; ind<N_gridphi_fieldperiod; ind++) axis_saved[ind] = malloc(4*sizeof(struct position));
		//axis = solve_magneticaxis(allparams, Bfield_axis, Xp, N_gridphi_fieldperiod);
		solve_magneticaxis_save(Xp->loc, axis_saved, allparams, Bfield_axis, N_gridphi_fieldperiod);
		axis = malloc(N_gridphi_fieldperiod*sizeof(struct position));
		for (ind=0; ind<N_gridphi_fieldperiod; ind++) axis[ind] = axis_saved[ind][0];

		//adjXp = adj2x2(Xp->tangent, &det);
		//lambdamu_saved = solve_lambdaRes(Bfield_axis, axis_saved, adjXp, m0, 1, N_gridphi_fieldperiod);
		printf("<--\n");
		clock_t t_afteraxis = clock();
		printf("... %f seconds ...\n", (double) (t_afteraxis-start)/CLOCKS_PER_SEC);
		printf("m0_fieldperiods = %d\n", m0);
		axis_phi0[0] = axis[0].loc[0];  axis_phi0[1] = axis[0].loc[1];
		printstructposition("axis", Xp);
		linalg2x2(Xp->tangent, evec, eval, &det, &trace);
		//printf("evalss=(%f, %f)\n", eval[0], eval[1]);
		printstructposition("axis", Xp);
		iota_axis = -eval[1]*m0/(2*M_PI);
		printf("iota_axis=%14.14f\n", iota_axis);
		//shapeRes = calloc(3,sizeof(double));
		//solve_gradRes(shapeRes, Bfield_axis, axis_saved, lambdamu_saved, 1, N_gridphi_fieldperiod, allparams, 0, 0) ;
		////grad_iota[i] = 0.5*m0*shapeRes[i]/(2.0*M_PI*sin(-eval[1]))
		//printf("iota_grad=%14.14f\n", 0.5*m0*shapeRes[0]/(2.0*M_PI*sin(-eval[1])));
		//printf("iota_grad=%14.14f\n", 0.5*m0*shapeRes[1]/(2.0*M_PI*sin(-eval[1])));
		//printf("iota_grad=%14.14f\n", 0.5*m0*shapeRes[2]/(2.0*M_PI*sin(-eval[1])));
		//printf("This number should be one: %f\n", eval[0]);
	}
	/* make Poincare plots (optional) */
	if (poiniota == 'y') {
		Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		//r_interval = 0.001; //n_points = 10;
		printf("About to enter Poincare module\n");
		printstructposition("axis", axis);
		iotaprofile(*axis, m0, N_gridphi_fieldperiod, rmin, zmin, iota, allparams);
		clock_t tPoincare = clock();
		printf("Time after filling Poincare plot file: %f\n", (double) (tPoincare-start)/CLOCKS_PER_SEC);
	}

	if (findisla == 'y') {
		//Xp->loc[0] = 1.3; Xp->loc[1]= -0.515; // island(?) is at(1.299008, -0.515194)
		//Xp->loc[0] = 0.945; Xp->loc[1]= 0.0;//Dommaschk (5,2) amp 1.73: island is at(1.299008, -0.515194)
		//Xp->loc[0] = 1.091; Xp->loc[1]= 0.0; // Dommaschk (5,2) amp 0.00001 on top of loop current
		//Xp->loc[0] = 1.0097; Xp->loc[1]= 0.0; // Dommaschk (5,2) amp 0.00001 on top of loop current
		if (strncmp(allparams.type, "coil", 4) == 0) {
			Xp->loc[0] = 1.299009; Xp->loc[1]= -0.515194; // NCSX: island is at(1.299008, -0.515194)
		}
		else if (strncmp(allparams.type, "Reim", 4) == 0) {
			Xp->loc[0] = 1.21; Xp->loc[1]= 0.0;  // Reiman 6-fold island chain
		}
		else if (strncmp(allparams.type, "heli", 4) == 0) {
			Xp->loc[0] = 0.87; Xp->loc[1]= 0.015;  // Reiman 6-fold island chain
			Xp->loc[0] = 1.00910608999481 ; Xp->loc[1]= -0.00000001285434;  // island for Cary Hanson case with I =0.0307
			Xp2->loc[0] = 0.8243366298; Xp2->loc[1] = 0.0000000608;  // X point for Cary Hanson case with I =0.0307
			printstructposition("Xp initial guess = \n", Xp);
		}
		Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		Xp2->tangent[0][0]=1.0; Xp2->tangent[0][1]=0.0; Xp2->tangent[1][0]=0.0; Xp2->tangent[1][1]=1.0;
		if (tor_mode % m0 == 0) 
			L_fixedpoints = pol_mode;	
		else 				 
			L_fixedpoints = m0*pol_mode;
		N_line = L_fixedpoints*N_gridphi_fieldperiod;
		Bfield_island = malloc(N_line*sizeof(struct field));
		Bfield_island2 = malloc(N_line*sizeof(struct field));
		island_center_saved = malloc(N_line*sizeof(struct position));
		island_center_saved2 = malloc(N_line*sizeof(struct position));
		for (ind=0; ind<N_line; ind++) {
			Bfield_island[ind] = malloc(4*sizeof(struct field));
			island_center_saved[ind] = malloc(4*sizeof(struct position));
			Bfield_island2[ind] = malloc(4*sizeof(struct field));
			island_center_saved2[ind] = malloc(4*sizeof(struct position));
		}
		printf("Entering function that solves for the magnetic island centre (O point)-->\n");
		//island_center = solve_islandcenter(allparams, Bfield_island, Xp, L_fixedpoints, N_gridphi_fieldperiod);
		//island_center2 = solve_islandcenter(allparams, Bfield_island2, Xp2, L_fixedpoints, N_gridphi_fieldperiod);
		island_center_saved[0][0].loc[0] = Xp->loc[0]; island_center_saved[0][0].loc[1] = Xp->loc[1];
		//solve_islandcenter_save(island_center_saved, &Res, allparams, Bfield_island, L_fixedpoints, N_gridphi_fieldperiod);
		solve_islandcenter_save(island_center_saved, &Res, allparams, Bfield_island, L_fixedpoints, N_gridphi_fieldperiod);
		island_center = malloc(L_fixedpoints*sizeof(struct position));
		for (ind=0; ind<L_fixedpoints; ind++) island_center[ind] = island_center_saved[ind*N_gridphi_fieldperiod][0];
		//solve_islandcenter_save(Xp2->loc, island_center_saved2, &Res2, allparams, Bfield_island2, L_fixedpoints, N_gridphi_fieldperiod);
		printf("island center at position (%14.14f, %14.14f)\n", island_center_saved[0][0].loc[0], island_center_saved[0][0].loc[1]);
		printf("island center at position (%14.14f, %14.14f)\n", island_center->loc[0], island_center->loc[1]);
		printf("<--\n");
		solve_islandcenter_save(island_center_saved, &Res, allparams, Bfield_island, L_fixedpoints, N_gridphi_fieldperiod);
		//exit(0);
		clock_t int_afterislandcentre = clock();
		printf("... %f seconds ...\n", (double) (int_afterislandcentre-start)/CLOCKS_PER_SEC);
		num_turns=1;
		ext_center = malloc(L_fixedpoints*sizeof(struct ext_position));
		ext_center2 = malloc(L_fixedpoints*sizeof(struct ext_position));
		ext_center_FD[0] = malloc(L_fixedpoints*sizeof(struct ext_position));
		ext_center_FD[1] = malloc(L_fixedpoints*sizeof(struct ext_position));
		ext_center[0].loc[0] = island_center[0].loc[0]; ext_center[0].loc[1] = island_center[0].loc[1];
		ext_center2[0].loc[0] = island_center2[0].loc[0]; ext_center2[0].loc[1] = island_center2[0].loc[1];
		//for (col=0; col< L_fixedpoints; col++) {
		//	ext_center[col].q0_index = 0;
		//}
		printf("Entering function that solves for the magnetic island centre (O point) in extended way\n");
		extsolve_periodicfieldline(ext_center, island_center, axis, m0, L_fixedpoints, pol_mode, q0_index, N_gridphi_fieldperiod);
		printf("Residue of fixed point = %f\n", ext_center[0].Res);
		//extsolve_periodicfieldline(ext_center2, axis, m0, L_fixedpoints, pol_mode, q0_index2, N_gridphi_fieldperiod, Bfield_island2);
		printf("Residue of fixed point = %f\n", ext_center2[0].Res);
		//storeq = ext_center->q0_index;
		clock_t int_afterextended = clock();
		printf("... %f seconds ...\n", (double) (int_afterextended-start)/CLOCKS_PER_SEC);
		//exit(0);
		printf("List of O point (R,Z) co-ordinates in island chain-->\n");
		for (row = 0; row<L_fixedpoints; row++) {
			printf("island_center[%d].loc = (%f, %f)\n", row, ext_center[row].loc[0], ext_center[row].loc[1]);
		}
		printf("<--\n");
		Sigma = malloc(L_fixedpoints*sizeof(double));
		//printf("Entering function that calculates the width of all islands in this chain\n");
		//width = calc_islandwidth(&circ, Sigma, ext_center, m0, L_fixedpoints, pol_mode);
		//printf("Left function that calculates the width of all islands in this chain\n");
		clock_t int_afterwidth = clock();
		printf("... %f seconds ...\n", (double) (int_afterwidth-start)/CLOCKS_PER_SEC);
		printf("For each of the %d islands, list:\nΣ_k = sum of tangent map matrix elements, w_k = island width of each individual island-->\n", L_fixedpoints);
		printf("<--\n");
		printf("Entering function that solves for the adjoint variable λ (circumference)\n");
		lambda = solve_lambda_circ(ext_center, m0, L_fixedpoints, N_gridphi_fieldperiod);
		printf("<--\n");
		clock_t int_afterlambda = clock();
		printf("... %f seconds ...\n", (double) (int_afterlambda-start)/CLOCKS_PER_SEC);
		//printf("%f\n", Bfield_island[0][0].value[0]);
		//printf("time after following adjoint variable for circumference: %f\n", (double) (int3-start)/CLOCKS_PER_SEC);
		//shapecirc = shapecircumference(coils, *n_coils, *n_consts, ext_center, &lambda, m0, N_gridphi_fieldperiod, tor_mode, pol_mode) ;
		//printf("time after calculating shape gradient of circumference: %f\n", (double) (int4-start)/CLOCKS_PER_SEC);
		printf("Entering function that solves for the adjoint variable μ\n");
		mu = solve_mu_tangent(ext_center, m0, L_fixedpoints, *q0_index, N_gridphi_fieldperiod);
		printf("<--\n");
		printf("Entering function that solves for the adjoint variable λ_Q -->\n");
		lambdaQ = solve_lambda_tangent(Bfield_island, ext_center, mu, m0, L_fixedpoints, *q0_index, N_gridphi_fieldperiod);
		printf("<--\n");
		printmat("adj_full_tangent", ext_center[0].adj_full_tangent, 2, 2);
		printmat("full_tangent", ext_center[0].full_tangent, 2, 2);
		lambdamu_saved = solve_lambdaRes(Bfield_island, island_center_saved, ext_center[0].adj_full_tangent, m0, L_fixedpoints, N_gridphi_fieldperiod);
		//lambdamu_saved2 = solve_lambdaRes(Bfield_island2, island_center_saved2, ext_center2[0].adj_full_tangent, m0, L_fixedpoints, N_gridphi_fieldperiod);
		clock_t int_afteralladjoints = clock();
		printf("... %f seconds ...\n", (double) (int_afteralladjoints-start)/CLOCKS_PER_SEC);

		if (checfast == 'y') max_coil_ind = 1;
		else max_coil_ind = allparams.n_coils;
		
		printf("N points per field period = %d\n", N_gridphi_fieldperiod);
	
		if (calcshap == 'y') {
			FILE *file_shapecirc, *file_shapetan, *file_shapewidth, *file_shapeRes;
			if ( (file_shapeRes = fopen("shape_Res.txt", "w") ) == NULL ) {	
				printf("Oops: couldn't open shape_Res.txt\n");
				exit(EXIT_FAILURE);
			}	
			if ( (file_shapecirc = fopen("shape_circ.txt", "w") ) == NULL ) {	
				printf("Oops: couldn't open shape_circ.txt\n");
				exit(EXIT_FAILURE);
			}	
			if ( (file_shapetan = fopen("shape_tan.txt", "w") ) == NULL ) {	
				printf("Oops: couldn't open shape_tan.txt\n");
				exit(EXIT_FAILURE);
			}	
			if ( (file_shapewidth = fopen("shape_width.txt", "w") ) == NULL ) {	
				printf("Oops: couldn't open shape_width.txt\n");
				exit(EXIT_FAILURE);
			}	
			shapecirc = calloc(allparams.n_diff,sizeof(double));
			shapeRes = calloc(allparams.n_diff,sizeof(double));
			shapetan = malloc(L_fixedpoints*sizeof(double));
			shapewidth = malloc(L_fixedpoints*sizeof(double));
			realshapecirc = calloc(allparams.n_diff,sizeof(double));
			realshapeRes = calloc(allparams.n_diff,sizeof(double));
			realshapetan = malloc(L_fixedpoints*sizeof(double));
			realshapewidth = malloc(L_fixedpoints*sizeof(double));
			for (ind=0; ind<L_fixedpoints; ind++) {
				shapetan[ind] = calloc(allparams.n_diff, sizeof(double));
				shapewidth[ind] = calloc(allparams.n_diff,sizeof(double));
				realshapetan[ind] = calloc(allparams.n_diff, sizeof(double));
				realshapewidth[ind] = calloc(allparams.n_diff,sizeof(double));
			}
			for (coil_ind=0; coil_ind<max_coil_ind; coil_ind++) {
				printf("coil_ind=%d/%d\n", coil_ind, *n_coils);
				if (strncmp(allparams.type, "coils", 4) == 0) num_segs = allparams.intparams[coil_ind];
				else num_segs = 1;
				for ( seg_ind=0; seg_ind < num_segs  ; seg_ind+=20) {
					printf("seg_ind=%d/%d\n", seg_ind, allparams.intparams[coil_ind]);
					solve_gradRes(shapeRes, Bfield_island, island_center_saved, lambdamu_saved, L_fixedpoints, N_gridphi_fieldperiod, allparams, coil_ind, seg_ind) ;
					printf("Entering function that solves for the gradient of the circumference\n");
					solve_gradcirc(shapecirc, Bfield_island, ext_center, &lambda, L_fixedpoints, N_gridphi_fieldperiod, allparams, coil_ind, seg_ind) ;
					printf("<--\n");
					printf("Entering function that solves for the gradient of Σ\n");
					solve_gradtangent(shapetan, Bfield_island, ext_center, lambdaQ, mu, L_fixedpoints, *q0_index, N_gridphi_fieldperiod, allparams, coil_ind, seg_ind) ;
					printf("<--\n");
					for (row = 0; row<allparams.n_diff; row++) {
						//shapeRes[row]*=(-1.0/4);
						for (ind=0; ind<L_fixedpoints; ind++) {
							shapewidth[ind][row] = ext_center[ind].width * ( shapecirc[row]/ext_center[0].circ - shapetan[ind][row]/ext_center[ind].Sigma );
							printf("shapeRes = %f\n", shapeRes[row]);
							printf("shapecirc = %f\n", shapecirc[row]);
							printf("shapetan[%d] = %f\n", ind, shapetan[ind][row]);
							printf("shapewidth[%d] = %f\n", ind, shapewidth[ind][row]);
						}
					}
					for (ind=0; ind<L_fixedpoints; ind++) {
						//printf("shapetan[%d] = %f\n", ind, shapetan[ind][0]);
						for (row = 0; row<allparams.n_diff; row++) {
							realshapeRes[row] = shapeRes[row];
							realshapecirc[row] = shapecirc[row];
							realshapewidth[ind][row] = shapewidth[ind][row];
							realshapetan[ind][row] = shapetan[ind][row];
						}
					}
					if (strncmp(allparams.type, "coil", 4) == 0) {
						for (row=0; row<3;row++) {
							if (seg_ind != allparams.intparams[coil_ind]-1) 
							dr[row] = coils[coil_ind][seg_ind+1][row] - coils[coil_ind][seg_ind][row];
							else
							dr[row] = coils[coil_ind][0][row] - coils[coil_ind][seg_ind][row];

							//printf("dr = %f\n", dr[row]);
							if (seg_ind != 0) 
							drdown[row] = - coils[coil_ind][seg_ind-1][row] + coils[coil_ind][seg_ind][row];
							else 
							drdown[row] = - coils[coil_ind][allparams.intparams[coil_ind]-1][row] + coils[coil_ind][seg_ind][row];
							//printf("drdown = %f\n", drdown[row]);
						}
						dldown = pow(drdown[0]*drdown[0] + drdown[1]*drdown[1] + drdown[2]*drdown[2], 0.5);
						dl = pow(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2], 0.5);
						factor = dldown/(dl+dldown);
						//printf("factor = %f\n", factor);
						//printf("seg_ind = %d, max_seg_ind = %d\n", seg_ind, (*n_consts)[coil_ind]);
						//printf("coils[%d][%d][%d] = %f, coils[%d][%d][%d] = %f\n", coil_ind, ((*n_consts)[coil_ind]+seg_ind-1)%(*n_consts)[coil_ind], 0, coils[coil_ind][((*n_consts)[coil_ind]+seg_ind-1)%(*n_consts)[coil_ind]][0], coil_ind, seg_ind, 0, coils[coil_ind][seg_ind][0]);
						for (row = 0; row<3; row++) {
							//realshapecirc[row] = factor*(dr[(row+1)%3]*shapecirc[(row+2)%3] - dr[(row+2)%3]*shapecirc[(row+1)%3]) + (1.0 - factor)*(drdown[(row+1)%3]*shapecirc[(row+2)%3] - drdown[(row+2)%3]*shapecirc[(row+1)%3]);
							realshapecirc[row] = dr[(row+1)%3]*shapecirc[(row+2)%3] - dr[(row+2)%3]*shapecirc[(row+1)%3];
							realshapeRes[row] = dr[(row+1)%3]*shapeRes[(row+2)%3] - dr[(row+2)%3]*shapeRes[(row+1)%3];
							printf("realshapecirc_%d = %f\n", row, realshapecirc[row]);
							printf("realshapeRes_%d = %f\n", row, realshapeRes[row]);
						}
						for (ind=0; ind<L_fixedpoints; ind++) {
							//printf("shapetan[%d] = %f\n", ind, shapetan[ind][0]);
							for (row = 0; row<3; row++) shapewidth[ind][row] = ext_center[ind].width* ( shapecirc[row]/ext_center[0].circ - shapetan[ind][row]/ext_center[ind].Sigma );
							//printf("shapetan[%d] = %f\n", ind, shapetan[ind][row]);
							//printf("shapewidth[%d] = %f\n", ind, shapewidth[ind][row]);
							for (row = 0; row<3; row++) {
								realshapetan[ind][row] = dr[(row+1)%3]*shapetan[ind][(row+2)%3] - dr[(row+2)%3]*shapetan[ind][(row+1)%3];
								realshapewidth[ind][row] = dr[(row+1)%3]*shapewidth[ind][(row+2)%3] - dr[(row+2)%3]*shapewidth[ind][(row+1)%3];
								//printf("realshapetan[%d][%d] = %f\n", ind, row, realshapetan[ind][row]);
								//printf("realshapewidth[%d][%d] = %f\n", ind, row, realshapewidth[ind][row]);
							}
							//printf("shapecirc, circ, shapetan, Sigma, width = %f %f %f %f %f\n", shapecirc[0], circ, shapetan[ind][0], Sigma[ind], width[ind]); 
						}
					}
						//free(shapecirc);
						//free(shapetan[centre_ind]
					fprintf(file_shapeRes, "%f %f %f\n", realshapeRes[0], realshapeRes[1], realshapeRes[2]);
					fprintf(file_shapecirc, "%f %f %f\n", realshapecirc[0], realshapecirc[1], realshapecirc[2]);
					fprintf(file_shapetan, "%f %f %f\n", realshapetan[0][0], realshapetan[0][1], realshapetan[0][2]);
					fprintf(file_shapewidth, "%f %f %f\n", realshapewidth[0][0], realshapewidth[0][1], realshapewidth[0][2]);
				}
				fprintf(file_shapeRes, "\n");
				fprintf(file_shapecirc, "\n");
				fprintf(file_shapetan, "\n");
				fprintf(file_shapewidth, "\n");
			}
		}
			//else if (strncmp(allparams.type, "coil", 4) == 0) {
			//	shapecirc = calloc(3,sizeof(double));
			//	shapeRes = calloc(3,sizeof(double));
			//	realshapecirc = calloc(3,sizeof(double));
			//	shapetan = malloc(L_fixedpoints*sizeof(double));
			//	shapewidth = malloc(L_fixedpoints*sizeof(double));
			//	realshapetan = malloc(L_fixedpoints*sizeof(double));
			//	realshapewidth = malloc(L_fixedpoints*sizeof(double));
			//	for (ind=0; ind<L_fixedpoints; ind++) {
			//		shapetan[ind] = calloc(3, sizeof(double));
			//		shapewidth[ind] = calloc(3,sizeof(double));
			//		realshapewidth[ind] = calloc(3,sizeof(double));
			//		realshapetan[ind] = calloc(3,sizeof(double));
			//	}
			//	FILE *file_shapecirc, *file_shapetan, *file_shapewidth;
			//	if ( (file_shapecirc = fopen("shape_circ.txt", "w") ) == NULL ) {	
			//		printf("Oups: couldn't open shape_circ.txt\n");
			//		exit(EXIT_FAILURE);
			//	}	
			//	if ( (file_shapetan = fopen("shape_tan.txt", "w") ) == NULL ) {	
			//		printf("Oups: couldn't open shape_tan.txt\n");
			//		exit(EXIT_FAILURE);
			//	}	
			//	if ( (file_shapewidth = fopen("shape_width.txt", "w") ) == NULL ) {	
			//		printf("Oups: couldn't open shape_width.txt\n");
			//		exit(EXIT_FAILURE);
			//	}	
			//	for (coil_ind=0; coil_ind<max_coil_ind; coil_ind++) {
			//		printf("coil_ind=%d/%d\n", coil_ind, *n_coils);
			//		num_segs = allparams.intparams[coil_ind];
			//		for ( seg_ind=0; seg_ind < num_segs  ; seg_ind+=20) {
			//			printf("seg_ind=%d/%d\n", seg_ind, allparams.intparams[coil_ind]);
			//			printf("Entering function that solves for the gradient of the circumference-->\n");
			//			printf("n_diffparams = %d\n", allparams.n_diff);
			//			solve_gradcirc(shapecirc, Bfield_island, ext_center, &lambda, L_fixedpoints, N_gridphi_fieldperiod, allparams, coil_ind, seg_ind) ;
			//			//printf("shapecirc_%d = %f\n", row, shapecirc[0]);
			//			printf("<--\n");
			//			printf("Entering function that solves for the gradient of Σ-->\n");
			//			solve_gradtangent(shapetan, Bfield_island, ext_center, lambdaQ, mu, L_fixedpoints, *q0_index, N_gridphi_fieldperiod, allparams, coil_ind, seg_ind) ;
			//			printf("<--\n");
			//			//printf("shapecirc = %f\n", *shapecirc);
			//			for (row=0; row<3;row++) {
			//				if (seg_ind != allparams.intparams[coil_ind]-1) 
			//				dr[row] = coils[coil_ind][seg_ind+1][row] - coils[coil_ind][seg_ind][row];
			//				else
			//				dr[row] = coils[coil_ind][0][row] - coils[coil_ind][seg_ind][row];

			//				//printf("coils = %f\n", coils[coil_ind][seg_ind][row]); printf("coils = %f\n", coils[coil_ind][seg_ind+1][row]);
			//				printf("dr = %f\n", dr[row]);
			//				if (seg_ind != 0) 
			//				drdown[row] = - coils[coil_ind][seg_ind-1][row] + coils[coil_ind][seg_ind][row];
			//				else 
			//				drdown[row] = - coils[coil_ind][allparams.intparams[coil_ind]-1][row] + coils[coil_ind][seg_ind][row];
			//				//printf("drdown = %f\n", drdown[row]);
			//			}
			//			dldown = pow(drdown[0]*drdown[0] + drdown[1]*drdown[1] + drdown[2]*drdown[2], 0.5);
			//			dl = pow(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2], 0.5);
			//			factor = dldown/(dl+dldown);
			//			//printf("factor = %f\n", factor);
			//			//printf("seg_ind = %d, max_seg_ind = %d\n", seg_ind, (*n_consts)[coil_ind]);
			//			//printf("coils[%d][%d][%d] = %f, coils[%d][%d][%d] = %f\n", coil_ind, ((*n_consts)[coil_ind]+seg_ind-1)%(*n_consts)[coil_ind], 0, coils[coil_ind][((*n_consts)[coil_ind]+seg_ind-1)%(*n_consts)[coil_ind]][0], coil_ind, seg_ind, 0, coils[coil_ind][seg_ind][0]);
			//			for (row = 0; row<3; row++) {
			//				realshapecirc[row] = factor*(dr[(row+1)%3]*shapecirc[(row+2)%3] - dr[(row+2)%3]*shapecirc[(row+1)%3]);
			//				realshapecirc[row] += (1.0 - factor)*(drdown[(row+1)%3]*shapecirc[(row+2)%3] - drdown[(row+2)%3]*shapecirc[(row+1)%3]);
			//				printf("realshapecirc_%d = %f\n", row, realshapecirc[row]);
			//			}
			//			for (ind=0; ind<L_fixedpoints; ind++) {
			//				//printf("shapetan[%d] = %f\n", ind, shapetan[ind][0]);
			//				for (row = 0; row<3; row++) shapewidth[ind][row] = ext_center[ind].width* ( shapecirc[row]/ext_center[0].circ - shapetan[ind][row]/ext_center[ind].Sigma );
			//				//printf("shapetan[%d] = %f\n", ind, shapetan[ind][row]);
			//				//printf("shapewidth[%d] = %f\n", ind, shapewidth[ind][row]);
			//				for (row = 0; row<3; row++) {
			//					realshapetan[ind][row] = dr[(row+1)%3]*shapetan[ind][(row+2)%3] - dr[(row+2)%3]*shapetan[ind][(row+1)%3];
			//					realshapewidth[ind][row] = dr[(row+1)%3]*shapewidth[ind][(row+2)%3] - dr[(row+2)%3]*shapewidth[ind][(row+1)%3];
			//					//printf("realshapetan[%d][%d] = %f\n", ind, row, realshapetan[ind][row]);
			//					//printf("realshapewidth[%d][%d] = %f\n", ind, row, realshapewidth[ind][row]);
			//				}
			//				//printf("shapecirc, circ, shapetan, Sigma, width = %f %f %f %f %f\n", shapecirc[0], circ, shapetan[ind][0], Sigma[ind], width[ind]); 
			//			}
			//			//free(shapecirc);
			//			//free(shapetan[centre_ind]
			//			fprintf(file_shapecirc, "%f %f %f\n", realshapecirc[0], realshapecirc[1], realshapecirc[2]);
			//			fprintf(file_shapetan, "%f %f %f\n", realshapetan[0][0], realshapetan[0][1], realshapetan[0][2]);
			//			fprintf(file_shapewidth, "%f %f %f\n", realshapewidth[0][0], realshapewidth[0][1], realshapewidth[0][2]);
			//		}
			//		fprintf(file_shapecirc, "\n");
			//		fprintf(file_shapetan, "\n");
			//		fprintf(file_shapewidth, "\n");
			//	}
			//}
		clock_t int5 = clock();
		printf("time after solving for shape gradient: %f\n", (double) (int5-start)/CLOCKS_PER_SEC);
		
		if (optimizeheli == 1) {
			rowallowed[0][0] = 2; rowallowed[0][1] = 4; rowallowed[0][2] = 6;
			rowallowed[1][0] = 1; rowallowed[1][1] = 4; rowallowed[1][2] = 5;
			rowallowed[0][0] = 2;
			rowallowed[1][0] = 1;
			for (coil_ind = 0; coil_ind < allparams.n_coils; coil_ind++) {
				for (row = 0; row<7; row++) printf("coil # %d, parameter = %f\n", coil_ind, allparams.diffparams[coil_ind][0][row]);
			}
			shapeRes = calloc(7,sizeof(double));
			shapeRes2 = calloc(7,sizeof(double));
			L_fixedpoints = 3;
			//for (L_fixedpoints = 4; L_fixedpoints<9; L_fixedpoints++) {
				coil_ind = 0;
				Xp->loc[0] = island_center_saved[0][0].loc[0];
				Xp->loc[1] = island_center_saved[0][0].loc[1];
				solve_islandcenter_save(island_center_saved, &Res, allparams, Bfield_island, L_fixedpoints, N_gridphi_fieldperiod);
				solve_islandcenter_save(island_center_saved2, &Res2, allparams, Bfield_island2, L_fixedpoints, N_gridphi_fieldperiod);
				printf("Res = (%f, %f)\n", Res, Res2);
				printf("Res^2 + Res2^2 = %f\n", Res*Res+Res2*Res2);
				lambdamu_saved = solve_lambdaRes(Bfield_island, island_center_saved, ext_center[0].adj_full_tangent, m0, L_fixedpoints, N_gridphi_fieldperiod);
				lambdamu_saved2 = solve_lambdaRes(Bfield_island2, island_center_saved2, ext_center2[0].adj_full_tangent, m0, L_fixedpoints, N_gridphi_fieldperiod);
				solve_gradRes(shapeRes, Bfield_island, island_center_saved, lambdamu_saved, L_fixedpoints, N_gridphi_fieldperiod, allparams, coil_ind, 0) ;
				solve_gradRes(shapeRes2, Bfield_island2, island_center_saved2, lambdamu_saved2, L_fixedpoints, N_gridphi_fieldperiod, allparams, coil_ind, 0) ;
				shapeRes_largest = 0.0;
				for (row = 0; row<1; row++) {
					printf("shapeRes = %f\n", shapeRes[rowallowed[0][row]]);
					if (fabs(shapeRes[rowallowed[0][row]]) > fabs(shapeRes_largest)) {
						shapeRes_largest = shapeRes[rowallowed[0][row]];
						row_largest = rowallowed[0][row];
					}
				}
				count = 1;
				while (Res*Res+Res2*Res2 > 0.0001){
					allparams.diffparams[coil_ind][0][row_largest] -= 1.0*((Res*Res+Res2*Res2)/(2.0*Res*0.25*shapeRes[row_largest] + 2.0*Res2*0.25*shapeRes2[row_largest]));
					coil_ind = count % allparams.n_coils;
					Xp->loc[0] = island_center_saved[0][0].loc[0];
					Xp->loc[1] = island_center_saved[0][0].loc[1];
					solve_islandcenter_save(island_center_saved, &Res, allparams, Bfield_island, L_fixedpoints, N_gridphi_fieldperiod);
					printf("Res = %f\n", Res);
					lambdamu_saved = solve_lambdaRes(Bfield_island, island_center_saved, ext_center[0].adj_full_tangent, m0, L_fixedpoints, N_gridphi_fieldperiod);
					for (row = 0; row<5; row++) shapeRes[row] = 0.0;
					solve_gradRes(shapeRes, Bfield_island, island_center_saved, lambdamu_saved, L_fixedpoints, N_gridphi_fieldperiod, allparams, coil_ind, 0) ;
					
					//solve_islandcenter_save(Xp2->loc, island_center_saved2, &Res2, allparams, Bfield_island2, L_fixedpoints, N_gridphi_fieldperiod);
					//printf("Res = %f\n", Res2);
					//for (row = 0; row<5; row++) shapeRes2[row] = 0.0;
					//solve_gradRes(shapeRes2, Bfield_island2, island_center_saved2, lambdamu_saved2, L_fixedpoints, N_gridphi_fieldperiod, allparams, coil_ind, 0) ;
					shapeRes_largest = 0.0;
					for (row = 0; row<1; row++) {
						printf("shapeRes = %f\n", shapeRes[rowallowed[coil_ind][row]]);
						if (fabs(shapeRes[rowallowed[coil_ind][row]]) > fabs(shapeRes_largest)) {
							shapeRes_largest = shapeRes[rowallowed[coil_ind][row]];
							row_largest = rowallowed[coil_ind][row];
						}
					}
					count += 1; 
					printf("Res = %f\n", Res);
				} 
				for (coil_ind = 0; coil_ind < allparams.n_coils; coil_ind++) {
					for (row = 0; row<5; row++) printf("coil # %d, parameter = %f\n", coil_ind, allparams.diffparams[coil_ind][0][row]);
				}
			//}
		}

			//FILE *file_checkshapewidth;
			//if ( (file_checkshapewidth = fopen("check_shape_width.txt", "w") ) == NULL ) {	
			//	printf("Oups: couldn't open check_shape_width.txt\n");
			//	exit(EXIT_FAILURE);
			//}	
			//	printf("coil_ind=%d/%d\n", coil_ind, *n_coils);
		if (strncmp(allparams.type, "coil", 4) == 0) num_segs = allparams.intparams[coil_ind];
		else num_segs = 1; // for anything but the "coil" type, the 2nd index of diff_params is always 0
		if (checshap == 'y') {
			checkshapeRes = calloc(allparams.n_diff,sizeof(double));
			checkshapecirc = calloc(allparams.n_diff,sizeof(double));
			checkshapewidth = malloc(L_fixedpoints*sizeof(double));
			checkshapeSigma = malloc(L_fixedpoints*sizeof(double));
			for (ind=0; ind<L_fixedpoints; ind++) {
				checkshapewidth[ind] = calloc(allparams.n_diff,sizeof(double));
				checkshapeSigma[ind] = calloc(allparams.n_diff,sizeof(double));
			}
			delta = pow(10.0, powdeltain);
			for (deltaind =0; deltaind<deltaindmax; deltaind++) {
				for (coil_ind=0; coil_ind<allparams.n_coils; coil_ind++) {
					for (seg_ind=0; seg_ind<num_segs; seg_ind+=seg_step) {
						for (row=0; row<allparams.n_diff; row++) {
							if (row == 0) printf("gradient with respect to 1st parameter\n");
							else if (row == 1) printf("gradient with respect to 2nd parameter\n");
							else if (row == 2) printf("gradient with respect to 3rd parameter\n");
							else printf("gradient with respect to %dth parameter\n", row-1);
							//savedcoilseg = allparams.diffparams[0][0][row];
							//printf("q0_index = %d\n", *q0_index);
							//savedcoilseg = 0.001;
							//savedcoilseg = allparams.diffparams[0][0][row];
							savedcoilseg = 0.0001; // pick a reference value for finite difference of parameter
							savedcoilseg = ext_center[0].width / shapewidth[0][row];
							//savedcoilseg = Res / (0.25*shapeRes[row]);
							printf("FD step = %f\n", delta*savedcoilseg);
							// begin the finite differencing procedure
							// step up for finite difference
							allparams.diffparams[coil_ind][seg_ind][row] += (delta*savedcoilseg);
							// having perturbed the parameter, find the new periodic field line position
							solve_islandcenter_save(island_center_saved, &Res, allparams, Bfield_island, L_fixedpoints, N_gridphi_fieldperiod);
							// analyze the periodic field line 
							// find the new residue, circumference, Sigma and width
							// q0_index already contains a non-zero value, so it won't be reassigned below
							for (ind=0; ind<L_fixedpoints; ind++) island_center[ind] = island_center_saved[ind*N_gridphi_fieldperiod][0];
							printf("periodic field line position = (%f, %f)\n", ext_center_FD[1][0].loc[0], ext_center_FD[1][0].loc[1]);
							extsolve_periodicfieldline(ext_center_FD[1], island_center, axis, m0, L_fixedpoints, pol_mode, q0_index, N_gridphi_fieldperiod);
							//printf("q0_index = %d\n", *q0_index);
							// step back down to original parameter
							allparams.diffparams[coil_ind][seg_ind][row] -= (delta*savedcoilseg);
							// step down for centred difference
							allparams.diffparams[coil_ind][seg_ind][row] -= (delta*savedcoilseg);
							//printf("%10.10f %10.10f\n", coils[0][0][row], savedcoilseg);
							//solve_islandcenter_save(allparams, Bfield_island, Xp, L_fixedpoints, N_gridphi_fieldperiod);
							solve_islandcenter_save(island_center_saved, &Res, allparams, Bfield_island, L_fixedpoints, N_gridphi_fieldperiod);
							for (ind=0; ind<L_fixedpoints; ind++) island_center[ind] = island_center_saved[ind*N_gridphi_fieldperiod][0];
							// q0_index already contains a non-zero value, so it won't be reassigned here
							extsolve_periodicfieldline(ext_center_FD[0], island_center, axis, m0, L_fixedpoints, pol_mode, q0_index, N_gridphi_fieldperiod);
							//printf("q0_index = %d\n", *q0_index);
							checkshapecirc[row] = ( log(ext_center_FD[1][0].circ) - log(ext_center_FD[0][0].circ) ) / (2.0*delta*savedcoilseg);
							checkshapeRes[row] = ( log(fabs(ext_center_FD[1][0].Res)) - log(fabs(ext_center_FD[0][0].Res)) ) / (2.0*delta*savedcoilseg);
							// step back up to original parameter
							allparams.diffparams[coil_ind][seg_ind][row] += (delta*savedcoilseg);
							for (ind=0; ind<L_fixedpoints; ind++) {
								checkshapeSigma[ind][row] = ( log(ext_center_FD[1][ind].Sigma) - log(ext_center_FD[0][ind].Sigma) ) / (2.0*delta*savedcoilseg);
								checkshapewidth[ind][row] = ( log(ext_center_FD[1][ind].width) - log(ext_center_FD[0][ind].width) ) / (2.0*delta*savedcoilseg);
								//checkshapeSigma[ind][row] = ( newSigma[1][ind] - Sigma[ind] ) / (delta*savedcoilseg);
								//checkshapewidth[ind][row] = ( newwidth[1][ind] - width[ind] ) / (delta*savedcoilseg);
								if (ind==0) {
									printf("shapeRes = %9.9f (direct) %9.9f (FD)\n", 0.25*shapeRes[row]/ext_center[0].Res, checkshapeRes[row]);
									printf("shapecirc = %9.9f (direct) %9.9f (FD)\n", shapecirc[row]/ext_center[0].circ, checkshapecirc[row]);
									printf("shapetan[%d] = %9.9f (direct) %9.9f (FD)\n", ind, shapetan[ind][row]/ext_center[ind].Sigma, checkshapeSigma[ind][row]);   
									printf("shapewidth[%d] = %9.9f (direct) %9.9f (FD)\n", ind, shapewidth[ind][row]/ext_center[ind].width, checkshapewidth[ind][row]);
								}
							}
							normgradwidth[row] = shapewidth[0][row]/ext_center[0].width;
							normgradcirc[row] = shapecirc[row]/ext_center[0].circ;
							normgradSigma[row] = shapetan[0][row]/ext_center[0].Sigma;
							normgradRes[row] = 0.25*shapeRes[row]/ext_center[0].Res; 
							checknormgradwidth[row][deltaind] = checkshapewidth[0][row];
							checknormgradcirc[row][deltaind] = checkshapecirc[row];
							checknormgradSigma[row][deltaind] = checkshapeSigma[0][row];
							checknormgradRes[row][deltaind] = checkshapeRes[row];
							// I don't yet understand why this is not -0.25
						}
					}
				}
				delta *= (pow(10.0, steppowdelta/deltaindmax));
			}
		}
		clock_t int6 = clock();
		printf("time after solving for everything: %f\n", (double) (int6-start)/CLOCKS_PER_SEC);
	}
			//else if (strncmp(allparams.type, "heli", 4) == 0) {
			//	newSigma[0] = malloc(L_fixedpoints*sizeof(double));
			//	newSigma[1] = malloc(L_fixedpoints*sizeof(double));
			//	checkshapeRes = calloc(7,sizeof(double));
			//	checkshapecirc = calloc(7,sizeof(double));
			//	checkshapewidth = malloc(L_fixedpoints*sizeof(double));
			//	checkshapeSigma = malloc(L_fixedpoints*sizeof(double));
			//	for (ind=0; ind<L_fixedpoints; ind++) {
			//		checkshapewidth[ind] = calloc(7,sizeof(double));
			//		checkshapeSigma[ind] = calloc(7,sizeof(double));
			//	}
			//	delta = pow(10.0, powdeltain);
			//	for (deltaind =0; deltaind<deltaindmax; deltaind++) {
			//		for (row=0; row<7; row++) {
			//			if (row == 0) printf("gradient with respect to 1st parameter\n");
			//			else if (row == 1) printf("gradient with respect to 2nd parameter\n");
			//			else if (row == 2) printf("gradient with respect to 3rd parameter\n");
			//			else printf("gradient with respect to %dth parameter\n", row-1);
			//			//savedcoilseg = allparams.diffparams[0][0][row];
			//			savedcoilseg = 0.01; // pick a reference value for finite difference of parameter
			//			//savedcoilseg = Res / (0.25*shapeRes[row]);
			//			printf("FD step = %f\n", delta*savedcoilseg);
			//			// begin the finite differencing procedure
			//			// step up for finite difference
			//			allparams.diffparams[0][0][row] += (delta*savedcoilseg);
			//			// having perturbed the parameter, find the new periodic field line position
			//			island_center_FD[1] = solve_islandcenter(allparams, Bfield_island, Xp, L_fixedpoints, N_gridphi_fieldperiod);
			//			// analyze the periodic field line 
			//			// find the new residue, circumference, Sigma and width
			//			// q0_index already contains a non-zero value, so it won't be reassigned below
			//			ext_center_FD[1] = extsolve_periodicfieldline(island_center_FD[1].loc, axis, m0, L_fixedpoints, pol_mode, q0_index, N_gridphi_fieldperiod, Bfield_island);
			//			//printf("q0_index = %d\n", *q0_index);
			//			// step back down to original parameter
			//			allparams.diffparams[0][0][row] -= (delta*savedcoilseg);
			//			// step down for centred difference
			//			allparams.diffparams[0][0][row] -= (delta*savedcoilseg);
			//			//printf("%10.10f %10.10f\n", coils[0][0][row], savedcoilseg);
			//			island_center_FD[0] = solve_islandcenter(allparams, Bfield_island, Xp, L_fixedpoints, N_gridphi_fieldperiod);
			//			// q0_index already contains a non-zero value, so it won't be reassigned here
			//			ext_center_FD[0] = extsolve_periodicfieldline(island_center_FD[0][0].loc, axis, m0, L_fixedpoints, pol_mode, q0_index, N_gridphi_fieldperiod, Bfield_island);
			//			//printf("q0_index = %d\n", *q0_index);
			//			checkshapecirc[row] = ( log(ext_center_FD[1][0].circ) - log(ext_center_FD[0][0].circ) ) / (2.0*delta*savedcoilseg);
			//			checkshapeRes[row] = ( log(ext_center_FD[1][0].Res) - log(ext_center_FD[0][0].Res) ) / (2.0*delta*savedcoilseg);
			//			// step back up to original parameter
			//			allparams.diffparams[0][0][row] += (delta*savedcoilseg);
			//			//checkshapecirc[row] = ( circup - circ ) / (delta*savedcoilseg);
			//			for (ind=0; ind<L_fixedpoints; ind++) {
			//				checkshapeSigma[ind][row] = ( ext_center_FD[1][ind] - ext_center_FD[0][ind] ) / (2.0*delta*savedcoilseg);
			//				checkshapewidth[ind][row] = ( ext_center_FD[1][ind] - ext_center_FD[0][ind] ) / (2.0*delta*savedcoilseg);
			//				//checkshapeSigma[ind][row] = ( newSigma[1][ind] - Sigma[ind] ) / (delta*savedcoilseg);
			//				//checkshapewidth[ind][row] = ( newwidth[1][ind] - width[ind] ) / (delta*savedcoilseg);
			//				if (ind==0) {
			//					printf("shapeRes = %9.9f (direct) %9.9f (FD)\n", 0.25*shapeRes[row]/Res, checkshapeRes[row]);
			//					printf("shapecirc = %9.9f (direct) %9.9f (FD)\n", shapecirc[row]/circ, checkshapecirc[row]);
			//					printf("shapetan[%d] = %9.9f (direct) %9.9f (FD)\n", ind, shapetan[ind][row]/Sigma[ind], checkshapeSigma[ind][row]);   
			//					printf("shapewidth[%d] = %9.9f (direct) %9.9f (FD)\n", ind, shapewidth[ind][row]/width[ind], checkshapewidth[ind][row]);
			//				}
			//			}
			//			normgradwidth[row] = shapewidth[0][row]/width[0];
			//			normgradcirc[row] = shapecirc[row]/circ;
			//			normgradSigma[row] = shapetan[0][row]/Sigma[0];
			//			normgradRes[row] = 0.25*shapeRes[row]/Res; // I don't yet understand why this is not -0.25
			//			//checknormgradwidth[row][deltaind] = checkshapewidth[0][row]/width[0];
			//			//checknormgradcirc[row][deltaind] = checkshapecirc[row]/circ;
			//			//checknormgradSigma[row][deltaind] = checkshapeSigma[0][row]/Sigma[0];
			//			checknormgradwidth[row][deltaind] = checkshapewidth[0][row];
			//			checknormgradcirc[row][deltaind] = checkshapecirc[row];
			//			checknormgradSigma[row][deltaind] = checkshapeSigma[0][row];
			//			checknormgradRes[row][deltaind] = checkshapeRes[row];
			//		}
			//		delta *= (pow(10.0, steppowdelta/deltaindmax));
			//	}
			//}

			//else if (strncmp(allparams.type, "coil", 4) == 0) {
			//	newSigma[0] = malloc(L_fixedpoints*sizeof(double));
			//	newSigma[1] = malloc(L_fixedpoints*sizeof(double));
			//	checkshapeRes = calloc(3,sizeof(double));
			//	checkshapecirc = calloc(3,sizeof(double));
			//	checkshapewidth = malloc(L_fixedpoints*sizeof(double));
			//	checkshapeSigma = malloc(L_fixedpoints*sizeof(double));
			//	for (ind=0; ind<L_fixedpoints; ind++) {
			//		checkshapewidth[ind] = calloc(3,sizeof(double));
			//		checkshapeSigma[ind] = calloc(3,sizeof(double));
			//	}
			//	FILE *file_checkshapewidth;
			//	if ( (file_checkshapewidth = fopen("check_shape_width.txt", "w") ) == NULL ) {	
			//		printf("Oups: couldn't open check_shape_width.txt\n");
			//		exit(EXIT_FAILURE);
			//	}	
			//	for (coil_ind=0; coil_ind<max_coil_ind; coil_ind++) {
			//		printf("coil_ind=%d/%d\n", coil_ind, *n_coils);
			//		max_seg_ind = allparams.intparams[coil_ind];
			//		for ( seg_ind=0; seg_ind < max_seg_ind; seg_ind+=20) {
			//			printf("seg_ind=%d/%d\n", seg_ind, (*n_consts)[coil_ind]);
			//			for (row=0; row<3; row++) {
			//				coils[coil_ind][seg_ind][row] += (delta) ;
			//				island_center_up = solve_islandcenter(allparams, Bfield_island, Xp, L_fixedpoints, N_gridphi_fieldperiod);
			//				ext_center_up = extsolve_periodicfieldline(island_center_up[0].loc, axis, m0, L_fixedpoints, pol_mode, q0_index, N_gridphi_fieldperiod, Bfield_island);
			//				//circup = 0.0;
			//				//newwidth[1] = calc_islandwidth(&circup, newSigma[1], ext_center_up, m0, L_fixedpoints, pol_mode);
			//				coils[coil_ind][seg_ind][row] -= (2.0*delta) ;
			//				island_center_down = solve_islandcenter(allparams, Bfield_island, Xp, L_fixedpoints, N_gridphi_fieldperiod);
			//				ext_center_down = extsolve_periodicfieldline(island_center_down[0].loc, axis, m0, L_fixedpoints, pol_mode, q0_index, N_gridphi_fieldperiod, Bfield_island);
			//				//circdown = 0.0;
			//				//newwidth[0] = calc_islandwidth(&circdown, newSigma[0], ext_center_down, m0, L_fixedpoints, pol_mode);
			//				coils[coil_ind][seg_ind][row] += (delta) ;
			//				checkshapecirc[row] = ( circup - circdown ) / (2.0*delta);
			//				printf("checkshapecirc[%d] = %f (FD)\n", row, checkshapecirc[row]);
			//				checkshapecirc[row] = ( circup - circ ) / (delta);
			//				printf("checkshapecirc[%d] = %f (FD)\n", row, checkshapecirc[row]);
			//				checkshapecirc[row] = ( circ - circdown ) / (delta);
			//				printf("checkshapecirc[%d] = %f (FD)\n", row, checkshapecirc[row]);
			//				for (ind=0; ind<1; ind++) {
			//					checkshapeSigma[ind][row] = ( newSigma[1][ind] - newSigma[0][ind] ) / (2.0*delta);
			//					checkshapewidth[ind][row] = ( newwidth[1][ind] - newwidth[0][ind] ) / (2.0*delta);
			//					//printf("checkshapetan[%d][%d] = %f (FD)\n", ind, row, checkshapeSigma[ind][row]);   
			//					printf("checkshapewidth[%d][%d] = %f (FD)\n", ind, row, checkshapewidth[ind][row]);
			//					//printf("shapecirc = %f (direct) %f (FD)\n", shapecirc[row], checkshapecirc[row]);
			//					//printf("shapetan[%d] = %f (direct) %f (FD)\n", ind, shapetan[ind][row], checkshapeSigma[ind][row]);   
			//					printf("shapewidth[%d] = %f (direct) %f (FD)\n", ind, realshapewidth[ind][row], checkshapewidth[ind][row]);
			//				}
			//			}
			//			//printf("shapecirc = %f\n", *shapecirc);
			//			//free(shapecirc);
			//			//free(shapetan[centre_ind]
			//			fprintf(file_checkshapewidth, "%f %f %f\n", checkshapewidth[0][0], checkshapewidth[0][1], checkshapewidth[0][2]);
			//		}
			//		fprintf(file_checkshapewidth, "\n");
			//	}
			//}
			

	if ( (checshap == 'y') && (strncmp(allparams.type, "Reim", 4) == 0) ) {
		printf("delta = [");
		delta = pow(10.0, powdeltain);
		for (deltaind = 0; deltaind<deltaindmax; deltaind++) {
			printf("%f ", delta);
			if (deltaind != deltaindmax - 1) printf(", ");
			delta *= (pow(10.0, steppowdelta/deltaindmax));
		}
		printf("]\n");
		for (row=0; row<3; row++) {
			printf("gw%d = %10.10f\ngw%d_delta = [", row, normgradwidth[row], row);
			for (deltaind = 0; deltaind<deltaindmax; deltaind++) {
				printf("%10.10f ", checknormgradwidth[row][deltaind]);
				if (deltaind != deltaindmax - 1) printf(", ");
			}
			printf("]\n");
			printf("gC%d = %10.10f\ngC%d_delta = [", row, normgradcirc[row], row);
			for (deltaind = 0; deltaind<deltaindmax; deltaind++) {
				printf("%10.10f ", checknormgradcirc[row][deltaind]);
				if (deltaind != deltaindmax - 1) printf(", ");
			}
			printf("]\n");
			printf("gSigma%d = %10.10f\ngSigma%d_delta = [", row, normgradSigma[row], row);
			for (deltaind = 0; deltaind<deltaindmax; deltaind++) {
				printf("%10.10f", checknormgradSigma[row][deltaind]);
				if (deltaind != deltaindmax - 1) printf(", ");
				
			}
			printf("]\n");
			printf("gRes%d = %10.10f\ngRes%d_delta = [", row, normgradRes[row], row);
			for (deltaind = 0; deltaind<deltaindmax; deltaind++) {
				printf("%10.10f", checknormgradRes[row][deltaind]);
				if (deltaind != deltaindmax - 1) printf(", ");
				
			}
			printf("]\n");
		}
	}

	if ( (checshap == 'y') && (strncmp(allparams.type, "heli", 4) == 0) ) {
		printf("delta = [");
		delta = pow(10.0, powdeltain);
		for (deltaind = 0; deltaind<deltaindmax; deltaind++) {
			printf("%f ", delta);
			if (deltaind != deltaindmax - 1) printf(", ");
			delta *= (pow(10.0, steppowdelta/deltaindmax));
		}
		printf("]\n");
		for (row=0; row<7; row++) {
			printf("gw[%d] = %10.10f\ngw_delta[%d] = [", row, normgradwidth[row], row);
			for (deltaind = 0; deltaind<deltaindmax; deltaind++) {
				printf("%10.10f ", checknormgradwidth[row][deltaind]);
				if (deltaind != deltaindmax - 1) printf(", ");
			}
			printf("]\n");
			printf("gC[%d] = %10.10f\ngC_delta[%d] = [", row, normgradcirc[row], row);
			for (deltaind = 0; deltaind<deltaindmax; deltaind++) {
				printf("%10.10f ", checknormgradcirc[row][deltaind]);
				if (deltaind != deltaindmax - 1) printf(", ");
			}
			printf("]\n");
			printf("gSigma[%d] = %10.10f\ngSigma_delta[%d] = [", row, normgradSigma[row], row);
			for (deltaind = 0; deltaind<deltaindmax; deltaind++) {
				printf("%10.10f", checknormgradSigma[row][deltaind]);
				if (deltaind != deltaindmax - 1) printf(", ");
				
			}
			printf("]\n");
			printf("gRes[%d] = %10.10f\ngRes_delta[%d] = [", row, normgradRes[row], row);
			for (deltaind = 0; deltaind<deltaindmax; deltaind++) {
				printf("%10.10f", checknormgradRes[row][deltaind]);
				if (deltaind != deltaindmax - 1) printf(", ");
				
			}
			printf("]\n");
		}
	}
	
	return 0;
}
