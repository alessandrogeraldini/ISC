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
	int N_gridphi_fieldperiod, N_gridphi_tor, num_turns=1, m0;
	int *n_coils=NULL, **n_consts, qq=1, *qq_segs=NULL, qq_segsq=1, pol_mode, tor_mode, ind, row, col;
	int find_axis = 1, deltaindmax = 20, deltaind;
	int L_fixedpoints, N_line;
	int coil_ind, seg_ind, max_coil_ind, max_seg_ind;
	int sizeline = 100;
	int q0val = 0, *q0_index=NULL; // It is important to initialize the integer stored in q0_index to 0
	struct position *Xp=calloc(1,sizeof(struct position)), *axis, *island_center, **mu, **lambdaQ, lambda, **lambdamu_saved;
	struct position *island_center_up, *island_center_down, **island_center_saved;
	struct ext_position *ext_center, *ext_center_up, *ext_center_down; //, *grad_ext_center;
	struct field **Bfield_island, **Bfield_axis, Bpoint;
	struct fieldparams allparams;
	//struct field *BB, *gradBB; for field checks
	double ***coils=NULL, *width, *shapeRes, *shapecirc, **shapetan, **shapewidth, **realshapetan, **realshapewidth, *realshapecirc, *Reimparam=NULL, circ, *Sigma;  
	double Res, Resup, Resdown;
	double powdeltain = -0.5, steppowdelta = -3.5, *doublepointer, **n_params=NULL;
	double phichec = 0.0, factor, dldown, dl, testpos[2] = {0.95, 0.0};
	struct field realBpointgrad[3];
	double *newwidth[2], **checkshapewidth, **checkshapeSigma, *checkshapecirc, *checkshapeRes, savedcoilseg;
	double circup, circdown, *newSigma[2], dr[3], drdown[3], *checkparam = malloc(4*sizeof(double));
	double *iota=NULL, *rmin=NULL, *zmin=NULL, axis_phi0[2];
	double **evec=malloc(2*sizeof(double*)), *eval=malloc(2*sizeof(double)), trace, det, iota_axis;
	double delta, checkBpointgrad[3];
	double normgradwidth[3], normgradcirc[3], normgradSigma[3], normgradRes[3], checknormgradRes[3][deltaindmax], checknormgradwidth[3][deltaindmax], checknormgradcirc[3][deltaindmax], checknormgradSigma[3][deltaindmax];
	struct field *Bpointgrad, Bpointdelta[3][2];
	char line[sizeline];
	char poiniota, findisla, calcshap, checshap, checfast, checkmag;
	FILE *OPTIONS;

	q0_index = &q0val;

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
	N_gridphi_fieldperiod = 60/m0;

	for (ind=0; ind < *n_coils; ind++) {
		printf("%d ", n_consts[0][ind]);
	}
	printf("\n");
	//printf("coils[0][0][0] = %f", coils[0][0][0]);
	//printf("n_coils = %d\n", *n_coils);
	//printf("n_consts = %d\n", **n_consts);
	if (strncmp(allparams.type, "coil", 4) == 0) {
		pol_mode = 4; tor_mode = 1;
		delta = 0.01;
	}
	else if (strncmp(allparams.type, "heli", 4) == 0) {
		pol_mode = 8; tor_mode = 5;
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
	//printf("n_coils = %d, n_consts = %d\n", *n_coils, **n_consts);
	Bpoint = Bfield(testpos, 0.0, allparams);
	printf("(BR, BZ, Bphi) = (%f, %f, %f)\n", Bpoint.value[0], Bpoint.value[1], Bpoint.value[2]);
	//exit(0);

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
				for (col=0;col<3;col++) 
				checkBpointgrad[col] = ( Bpointdelta[1][1].value[col] - Bpointdelta[1][0].value[col] ) / (2.0*delta);
				printf("checkgradB = (%14.14f, %14.14f, %14.14f) (direct), (%f, %f, %f) (FD) \n", Bpointgrad[row].value[0], Bpointgrad[row].value[1], Bpointgrad[row].value[2], checkBpointgrad[0], checkBpointgrad[1], checkBpointgrad[2]);
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
			Xp->loc[0] = 0.95; Xp->loc[1]= 0.0; //helical
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
		axis = solve_magneticaxis(allparams, Bfield_axis, Xp, N_gridphi_fieldperiod);
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
		printf("iota_axis=%f\n", iota_axis);
		//printf("This number should be one: %f\n", eval[0]);
	}
	/* make Poincare plots (optional) */
	if (poiniota == 'y') {
		Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		//r_interval = 0.001; //n_points = 10;
		printf("About to enter Poincare module\n");
		printstructposition("axis", axis);
		iotaprofile(allparams.type, *axis, m0, N_gridphi_fieldperiod, rmin, zmin, iota, coils, *n_coils, *n_consts);
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
			printstructposition("Xp initial guess = \n", Xp);
		}
		Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		if (tor_mode % m0 == 0) 
			L_fixedpoints = pol_mode;	
		else 				 
			L_fixedpoints = m0*pol_mode;
		N_line = L_fixedpoints*N_gridphi_fieldperiod;
		Bfield_island = malloc(N_line*sizeof(struct field));
		for (ind=0; ind<N_line; ind++) Bfield_island[ind] = malloc(4*sizeof(struct field));
		printf("Entering function that solves for the magnetic island centre (O point)-->\n");
		island_center = solve_islandcenter(allparams, Bfield_island, Xp, L_fixedpoints, N_gridphi_fieldperiod);
		//island_center_saved = solve_islandcenter_save(Xp->loc[0], Xp->loc[1], allparams, Bfield_island, L_fixedpoints, N_gridphi_fieldperiod);
		island_center_saved = solve_islandcenter_save(Xp->loc[0], Xp->loc[1], allparams, Bfield_island, L_fixedpoints, N_gridphi_fieldperiod);
		printf("island center at position (%14.14f, %14.14f)\n", island_center_saved[0][0].loc[0], island_center_saved[0][0].loc[1]);
		printf("island center at position (%14.14f, %14.14f)\n", island_center->loc[0], island_center->loc[1]);
		printf("<--\n");
		clock_t int_afterislandcentre = clock();
		printf("... %f seconds ...\n", (double) (int_afterislandcentre-start)/CLOCKS_PER_SEC);
		num_turns=1;
		printf("Entering function that solves for the magnetic island centre (O point) in extended way\n");
		ext_center = malloc(L_fixedpoints*sizeof(struct ext_position));
		//for (col=0; col< L_fixedpoints; col++) {
		//	ext_center[col].q0_index = 0;
		//}
		ext_center = solve_islandcenter_full(island_center[0].loc, axis, &num_turns, m0, L_fixedpoints, &Res, q0_index, N_gridphi_fieldperiod, Bfield_island);
		printf("Residue of fixed point = %f\n", Res);
		//storeq = ext_center->q0_index;
		printf("Left function that solves for the magnetic island centre (O point) in extended way\n");
		clock_t int_afterextended = clock();
		printf("... %f seconds ...\n", (double) (int_afterextended-start)/CLOCKS_PER_SEC);
		printf("List of O point (R,Z) co-ordinates in island chain-->\n");
		for (row = 0; row<L_fixedpoints; row++) {
			printf("island_center[%d].loc = (%f, %f)\n", row, ext_center[row].loc[0], ext_center[row].loc[1]);
		}
		printf("<--\n");
		Sigma = malloc(L_fixedpoints*sizeof(double));
		printf("Entering function that calculates the width of all islands in this chain\n");
		width = calc_islandwidth(&circ, Sigma, ext_center, m0, L_fixedpoints, pol_mode);
		printf("Left function that calculates the width of all islands in this chain\n");
		clock_t int_afterwidth = clock();
		printf("... %f seconds ...\n", (double) (int_afterwidth-start)/CLOCKS_PER_SEC);
		printf("Circumference of the island chain (sum of chords) = %10.10f\n", circ);
		printf("For each of the %d islands, list:\nΣ_k = sum of tangent map matrix elements, w_k = island width of each individual island-->\n", L_fixedpoints);
		for (row=0; row< L_fixedpoints; row++) {
			printf("Σ_%d = %10.10f, w_%d = %10.10f\n", row, Sigma[row], row, width[row]);
		}
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
		clock_t int_afteralladjoints = clock();
		printf("... %f seconds ...\n", (double) (int_afteralladjoints-start)/CLOCKS_PER_SEC);

		max_seg_ind = 3;
		if (checfast == 'y') max_coil_ind = 1;
		else max_coil_ind = *n_coils;
		
		printf("N points per field period = %d\n", N_gridphi_fieldperiod);
	
		if (calcshap == 'y') {
			if (strncmp(allparams.type, "Reim", 4) == 0) {
				shapecirc = calloc(3,sizeof(double));
				shapeRes = calloc(3,sizeof(double));
				shapetan = malloc(L_fixedpoints*sizeof(double));
				shapewidth = malloc(L_fixedpoints*sizeof(double));
				for (ind=0; ind<L_fixedpoints; ind++) {
					shapetan[ind] = calloc(3, sizeof(double));
					shapewidth[ind] = calloc(3,sizeof(double));
				}
				Reimparam = **coils;
				solve_gradRes(shapeRes, Bfield_island, island_center_saved, lambdamu_saved, L_fixedpoints, N_gridphi_fieldperiod, allparams, 0, 0) ;
				printf("Entering function that solves for the gradient of the circumference\n");
				solve_gradcirc(shapecirc, Bfield_island, ext_center, &lambda, L_fixedpoints, N_gridphi_fieldperiod, allparams, 0, 0) ;
				printf("<--\n");
				printf("Entering function that solves for the gradient of Σ\n");
				solve_gradtangent(shapetan, Bfield_island, ext_center, lambdaQ, mu, L_fixedpoints, *q0_index, N_gridphi_fieldperiod, allparams, 0, 0) ;
				printf("<--\n");
				for (row = 0; row<3; row++) {
					//shapeRes[row]*=(-1.0/4);
					for (ind=0; ind<L_fixedpoints; ind++) {
						shapewidth[ind][row] = width[ind]* ( shapecirc[row]/circ - shapetan[ind][row]/Sigma[ind] );
						printf("shapeRes = %f\n", shapeRes[row]);
						printf("shapecirc = %f\n", shapecirc[row]);
						printf("shapetan[%d] = %f\n", ind, shapetan[ind][row]);
						printf("shapewidth[%d] = %f\n", ind, shapewidth[ind][row]);
					}
				}
			}
			else if (strncmp(allparams.type, "heli", 4) == 0) {
				shapecirc = calloc(7,sizeof(double));
				shapeRes = calloc(7,sizeof(double));
				shapetan = malloc(L_fixedpoints*sizeof(double));
				shapewidth = malloc(L_fixedpoints*sizeof(double));
				for (ind=0; ind<L_fixedpoints; ind++) {
					shapetan[ind] = calloc(7, sizeof(double));
					shapewidth[ind] = calloc(7,sizeof(double));
				}
				solve_gradRes(shapeRes, Bfield_island, island_center_saved, lambdamu_saved, L_fixedpoints, N_gridphi_fieldperiod, allparams, 0, 0) ;
				printf("Entering function that solves for the gradient of the circumference\n");
				solve_gradcirc(shapecirc, Bfield_island, ext_center, &lambda, L_fixedpoints, N_gridphi_fieldperiod, allparams, 0, 0) ;
				printf("<--\n");
				printf("Entering function that solves for the gradient of Σ\n");
				solve_gradtangent(shapetan, Bfield_island, ext_center, lambdaQ, mu, L_fixedpoints, *q0_index, N_gridphi_fieldperiod, allparams, 0, 0) ;
				printf("<--\n");
				for (row = 0; row<7; row++) {
					//shapeRes[row]*=(-1.0/4);
					for (ind=0; ind<L_fixedpoints; ind++) {
						shapewidth[ind][row] = width[ind]* ( shapecirc[row]/circ - shapetan[ind][row]/Sigma[ind] );
						printf("shapeRes = %f\n", shapeRes[row]);
						printf("shapecirc = %f\n", shapecirc[row]);
						printf("shapetan[%d] = %f\n", ind, shapetan[ind][row]);
						printf("shapewidth[%d] = %f\n", ind, shapewidth[ind][row]);
					}
				}
			
			}
			else if (strncmp(allparams.type, "coil", 4) == 0) {
				shapecirc = calloc(3,sizeof(double));
				shapeRes = calloc(3,sizeof(double));
				realshapecirc = calloc(3,sizeof(double));
				shapetan = malloc(L_fixedpoints*sizeof(double));
				shapewidth = malloc(L_fixedpoints*sizeof(double));
				realshapetan = malloc(L_fixedpoints*sizeof(double));
				realshapewidth = malloc(L_fixedpoints*sizeof(double));
				for (ind=0; ind<L_fixedpoints; ind++) {
					shapetan[ind] = calloc(3, sizeof(double));
					shapewidth[ind] = calloc(3,sizeof(double));
					realshapewidth[ind] = calloc(3,sizeof(double));
					realshapetan[ind] = calloc(3,sizeof(double));
				}
				FILE *file_shapecirc, *file_shapetan, *file_shapewidth;
				if ( (file_shapecirc = fopen("shape_circ.txt", "w") ) == NULL ) {	
					printf("Oups: couldn't open shape_circ.txt\n");
					exit(EXIT_FAILURE);
				}	
				if ( (file_shapetan = fopen("shape_tan.txt", "w") ) == NULL ) {	
					printf("Oups: couldn't open shape_tan.txt\n");
					exit(EXIT_FAILURE);
				}	
				if ( (file_shapewidth = fopen("shape_width.txt", "w") ) == NULL ) {	
					printf("Oups: couldn't open shape_width.txt\n");
					exit(EXIT_FAILURE);
				}	
				for (coil_ind=0; coil_ind<max_coil_ind; coil_ind++) {
					printf("coil_ind=%d/%d\n", coil_ind, *n_coils);
					if (checfast != 'y') max_seg_ind = allparams.intparams[coil_ind];
					for ( seg_ind=0; seg_ind < max_seg_ind  ; seg_ind++) {
						printf("seg_ind=%d/%d\n", seg_ind, allparams.intparams[coil_ind]);
						printf("Entering function that solves for the gradient of the circumference-->\n");
						printf("n_diffparams = %d\n", allparams.n_diff);
						solve_gradcirc(shapecirc, Bfield_island, ext_center, &lambda, L_fixedpoints, N_gridphi_fieldperiod, allparams, coil_ind, seg_ind) ;
						//printf("shapecirc_%d = %f\n", row, shapecirc[0]);
						printf("<--\n");
						printf("Entering function that solves for the gradient of Σ-->\n");
						solve_gradtangent(shapetan, Bfield_island, ext_center, lambdaQ, mu, L_fixedpoints, *q0_index, N_gridphi_fieldperiod, allparams, coil_ind, seg_ind) ;
						printf("<--\n");
						//printf("shapecirc = %f\n", *shapecirc);
						for (row=0; row<3;row++) {
							if (seg_ind != allparams.intparams[coil_ind]-1) 
							dr[row] = coils[coil_ind][seg_ind+1][row] - coils[coil_ind][seg_ind][row];
							else
							dr[row] = coils[coil_ind][0][row] - coils[coil_ind][seg_ind][row];

							//printf("coils = %f\n", coils[coil_ind][seg_ind][row]); printf("coils = %f\n", coils[coil_ind][seg_ind+1][row]);
							printf("dr = %f\n", dr[row]);
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
							realshapecirc[row] = factor*(dr[(row+1)%3]*shapecirc[(row+2)%3] - dr[(row+2)%3]*shapecirc[(row+1)%3]);
							realshapecirc[row] += (1.0 - factor)*(drdown[(row+1)%3]*shapecirc[(row+2)%3] - drdown[(row+2)%3]*shapecirc[(row+1)%3]);
							printf("realshapecirc_%d = %f\n", row, realshapecirc[row]);
						}
						for (ind=0; ind<L_fixedpoints; ind++) {
							//printf("shapetan[%d] = %f\n", ind, shapetan[ind][0]);
							for (row = 0; row<3; row++) shapewidth[ind][row] = width[ind]* ( shapecirc[row]/circ - shapetan[ind][row]/Sigma[ind] );
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
						//free(shapecirc);
						//free(shapetan[centre_ind]
						fprintf(file_shapecirc, "%f %f %f\n", shapecirc[0], shapecirc[1], shapecirc[2]);
						fprintf(file_shapetan, "%f %f %f\n", shapetan[0][0], shapetan[0][1], shapetan[0][2]);
						fprintf(file_shapewidth, "%f %f %f\n", shapewidth[0][0], shapewidth[0][1], shapewidth[0][2]);
					}
					fprintf(file_shapecirc, "\n");
					fprintf(file_shapetan, "\n");
					fprintf(file_shapewidth, "\n");
				}
			}
		}
		clock_t int5 = clock();
		printf("time after solving for shape gradient: %f\n", (double) (int5-start)/CLOCKS_PER_SEC);
		
		if (checshap == 'y') {
			if (strncmp(allparams.type, "Reim", 4) == 0) {
				newSigma[0] = malloc(L_fixedpoints*sizeof(double));
				newSigma[1] = malloc(L_fixedpoints*sizeof(double));
				checkshapeRes = calloc(3,sizeof(double));
				checkshapecirc = calloc(3,sizeof(double));
				checkshapewidth = malloc(L_fixedpoints*sizeof(double));
				checkshapeSigma = malloc(L_fixedpoints*sizeof(double));
				for (ind=0; ind<L_fixedpoints; ind++) {
					checkshapewidth[ind] = calloc(3,sizeof(double));
					checkshapeSigma[ind] = calloc(3,sizeof(double));
				}
				delta = pow(10.0, powdeltain);
				for (deltaind =0; deltaind<deltaindmax; deltaind++) {
					for (row=0; row<3; row++) {
						printf("gradient with respect to parameter #%d\n", row+1);
						savedcoilseg = allparams.diffparams[0][0][row];
						savedcoilseg = width[0] / shapewidth[0][row];
						//savedcoilseg = 1.0;
						//coils[0][0][row] += (delta*savedcoilseg);
						allparams.diffparams[0][0][row] += (delta*savedcoilseg);
						//printf("%10.10f %10.10f\n", allparams.diffparams[0][0][row], savedcoilseg);
						//printf("Entering function that solves for the magnetic island centre (O point)-->\n");
						island_center_up = solve_islandcenter(allparams, Bfield_island, Xp, L_fixedpoints, N_gridphi_fieldperiod);
						//printf("<--\n");
						//printf("Entering function that solves for the magnetic island centre (O point) in extended way-->\n");
						//ext_center_up = malloc(L_fixedpoints*sizeof(struct ext_position));
						//ext_center_down = malloc(L_fixedpoints*sizeof(struct ext_position));
// since q0_index already contains a non-zero value, it won't be reassigned here
						ext_center_up = solve_islandcenter_full(island_center_up[0].loc, axis, &num_turns, m0, L_fixedpoints, &Resup, q0_index, N_gridphi_fieldperiod, Bfield_island);
						//printf("q0_index = %d\n", *q0_index);
						//printf("<--\n");
						newwidth[1] = calc_islandwidth(&circup, newSigma[1], ext_center_up, m0, L_fixedpoints, pol_mode);
						allparams.diffparams[0][0][row] -= (delta*savedcoilseg);
						allparams.diffparams[0][0][row] -= (delta*savedcoilseg);
						//printf("circup = %f\n", circup);
						//printf("%10.10f %10.10f\n", coils[0][0][row], savedcoilseg);
						island_center_down = solve_islandcenter(allparams, Bfield_island, Xp, L_fixedpoints, N_gridphi_fieldperiod);
// since q0_index already contains a non-zero value, it won't be reassigned here
						ext_center_down = solve_islandcenter_full(island_center_down[0].loc, axis, &num_turns, m0, L_fixedpoints, &Resdown, q0_index, N_gridphi_fieldperiod, Bfield_island);
						//printf("q0_index = %d\n", *q0_index);
						//ext_center_down->q0_index = storeq;
						newwidth[0] = calc_islandwidth(&circdown, newSigma[0], ext_center_down, m0, L_fixedpoints, pol_mode);
						circup =      log( circup);
						circdown =    log(circdown);
						printf("Resdown, Res, Resup = %f, %f, %f\n", Resdown, Res, Resup);
						printf("gradRes = %f, %f (num, dir)\n", (Resup - Resdown)/(2.0*delta*savedcoilseg), -0.25*shapeRes[row]) ;
						Resup = log(Resup);
						Resdown = log(Resdown);
						//printf("circdown = %f\n", circdown);
						//printf("circ = %f\n", circ);
						printf("delta = %f\nsavedcoilseg = %f\n", delta, savedcoilseg);
						//coils[0][0][row] += (delta*savedcoilseg);
						//printf("%10.10f %10.10f\n", coils[0][0][row], savedcoilseg);
						checkshapecirc[row] = ( circup - circdown ) / (2.0*delta*savedcoilseg);
						checkshapeRes[row] = ( Resup - Resdown ) / (2.0*delta*savedcoilseg);
						allparams.diffparams[0][0][row] += (delta*savedcoilseg);
						//checkshapecirc[row] = ( circup - circ ) / (delta*savedcoilseg);
						for (ind=0; ind<6; ind++) {
							newwidth[0][ind] = log(newwidth[0][ind]);
							newwidth[1][ind] = log(newwidth[1][ind]);
							newSigma[0][ind] = log(newSigma[0][ind]);
							newSigma[1][ind] = log(newSigma[1][ind]);
							checkshapeSigma[ind][row] = ( newSigma[1][ind] - newSigma[0][ind] ) / (2.0*delta*savedcoilseg);
							checkshapewidth[ind][row] = ( newwidth[1][ind] - newwidth[0][ind] ) / (2.0*delta*savedcoilseg);
							//checkshapeSigma[ind][row] = ( newSigma[1][ind] - Sigma[ind] ) / (delta*savedcoilseg);
							//checkshapewidth[ind][row] = ( newwidth[1][ind] - width[ind] ) / (delta*savedcoilseg);
							//printf("shapeRes = %9.9f (direct) %9.9f (FD)\n", shapeRes[row], checkshapeRes[row]);
							//printf("shapecirc = %9.9f (direct) %9.9f (FD)\n", shapecirc[row], checkshapecirc[row]);
							//printf("shapetan[%d] = %9.9f (direct) %9.9f (FD)\n", ind, shapetan[ind][row], checkshapeSigma[ind][row]);   
							//printf("shapewidth[%d] = %9.9f (direct) %9.9f (FD)\n", ind, shapewidth[ind][row], checkshapewidth[ind][row]);
						}
						normgradwidth[row] = shapewidth[0][row]/width[0];
						normgradcirc[row] = shapecirc[row]/circ;
						normgradSigma[row] = shapetan[0][row]/Sigma[0];
						normgradRes[row] = 0.25*shapeRes[row]/Res; // I don't yet understand why this is not -0.25
						//checknormgradwidth[row][deltaind] = checkshapewidth[0][row]/width[0];
						//checknormgradcirc[row][deltaind] = checkshapecirc[row]/circ;
						//checknormgradSigma[row][deltaind] = checkshapeSigma[0][row]/Sigma[0];
						checknormgradwidth[row][deltaind] = checkshapewidth[0][row];
						checknormgradcirc[row][deltaind] = checkshapecirc[row];
						checknormgradSigma[row][deltaind] = checkshapeSigma[0][row];
						checknormgradRes[row][deltaind] = checkshapeRes[row];
					}
					delta *= (pow(10.0, steppowdelta/deltaindmax));
				}
			}
			else if (strncmp(allparams.type, "heli", 4) == 0) {
				newSigma[0] = malloc(L_fixedpoints*sizeof(double));
				newSigma[1] = malloc(L_fixedpoints*sizeof(double));
				checkshapeRes = calloc(7,sizeof(double));
				checkshapecirc = calloc(7,sizeof(double));
				checkshapewidth = malloc(L_fixedpoints*sizeof(double));
				checkshapeSigma = malloc(L_fixedpoints*sizeof(double));
				for (ind=0; ind<L_fixedpoints; ind++) {
					checkshapewidth[ind] = calloc(7,sizeof(double));
					checkshapeSigma[ind] = calloc(7,sizeof(double));
				}
				delta = pow(10.0, powdeltain);
				for (deltaind =0; deltaind<deltaindmax; deltaind++) {
					for (row=0; row<7; row++) {
						printf("gradient with respect to parameter #%d\n", row+1);
						//savedcoilseg = allparams.diffparams[0][0][row];
						savedcoilseg = 0.001;
						//savedcoilseg = Res / (0.25*shapeRes[row]);
						allparams.diffparams[0][0][row] += (delta*savedcoilseg);
						island_center_up = solve_islandcenter(allparams, Bfield_island, Xp, L_fixedpoints, N_gridphi_fieldperiod);
// since q0_index already contains a non-zero value, it won't be reassigned here
						ext_center_up = solve_islandcenter_full(island_center_up[0].loc, axis, &num_turns, m0, L_fixedpoints, &Resup, q0_index, N_gridphi_fieldperiod, Bfield_island);
						//printf("q0_index = %d\n", *q0_index);
						//printf("<--\n");
						newwidth[1] = calc_islandwidth(&circup, newSigma[1], ext_center_up, m0, L_fixedpoints, pol_mode);
						allparams.diffparams[0][0][row] -= (delta*savedcoilseg);
						allparams.diffparams[0][0][row] -= (delta*savedcoilseg);
						//printf("circup = %f\n", circup);
						//printf("%10.10f %10.10f\n", coils[0][0][row], savedcoilseg);
						island_center_down = solve_islandcenter(allparams, Bfield_island, Xp, L_fixedpoints, N_gridphi_fieldperiod);
// since q0_index already contains a non-zero value, it won't be reassigned here
						ext_center_down = solve_islandcenter_full(island_center_down[0].loc, axis, &num_turns, m0, L_fixedpoints, &Resdown, q0_index, N_gridphi_fieldperiod, Bfield_island);
						//printf("q0_index = %d\n", *q0_index);
						//ext_center_down->q0_index = storeq;
						newwidth[0] = calc_islandwidth(&circdown, newSigma[0], ext_center_down, m0, L_fixedpoints, pol_mode);
						circup =      log( circup);
						circdown =    log(circdown);
						printf("Resdown, Res, Resup = %f, %f, %f\n", Resdown, Res, Resup);
						printf("gradRes = %f, %f (num, dir)\n", (Resup - Resdown)/(2.0*delta*savedcoilseg), -0.25*shapeRes[row]) ;
						Resup = log(Resup);
						Resdown = log(Resdown);
						//printf("circdown = %f\n", circdown);
						//printf("circ = %f\n", circ);
						printf("delta = %f\nsavedcoilseg = %f\n", delta, savedcoilseg);
						//coils[0][0][row] += (delta*savedcoilseg);
						//printf("%10.10f %10.10f\n", coils[0][0][row], savedcoilseg);
						checkshapecirc[row] = ( circup - circdown ) / (2.0*delta*savedcoilseg);
						checkshapeRes[row] = ( Resup - Resdown ) / (2.0*delta*savedcoilseg);
						allparams.diffparams[0][0][row] += (delta*savedcoilseg);
						//checkshapecirc[row] = ( circup - circ ) / (delta*savedcoilseg);
						for (ind=0; ind<L_fixedpoints; ind++) {
							newwidth[0][ind] = log(newwidth[0][ind]);
							newwidth[1][ind] = log(newwidth[1][ind]);
							newSigma[0][ind] = log(newSigma[0][ind]);
							newSigma[1][ind] = log(newSigma[1][ind]);
							checkshapeSigma[ind][row] = ( newSigma[1][ind] - newSigma[0][ind] ) / (2.0*delta*savedcoilseg);
							checkshapewidth[ind][row] = ( newwidth[1][ind] - newwidth[0][ind] ) / (2.0*delta*savedcoilseg);
							//checkshapeSigma[ind][row] = ( newSigma[1][ind] - Sigma[ind] ) / (delta*savedcoilseg);
							//checkshapewidth[ind][row] = ( newwidth[1][ind] - width[ind] ) / (delta*savedcoilseg);
							printf("yolo\n");
							printf("shapeRes = %9.9f (direct) %9.9f (FD)\n", 0.25*shapeRes[row]/Res, checkshapeRes[row]);
							printf("shapecirc = %9.9f (direct) %9.9f (FD)\n", shapecirc[row]/circ, checkshapecirc[row]);
							printf("shapetan[%d] = %9.9f (direct) %9.9f (FD)\n", ind, shapetan[ind][row]/Sigma[ind], checkshapeSigma[ind][row]);   
							printf("shapewidth[%d] = %9.9f (direct) %9.9f (FD)\n", ind, shapewidth[ind][row]/width[ind], checkshapewidth[ind][row]);
						}
						normgradwidth[row] = shapewidth[0][row]/width[0];
						normgradcirc[row] = shapecirc[row]/circ;
						normgradSigma[row] = shapetan[0][row]/Sigma[0];
						normgradRes[row] = 0.25*shapeRes[row]/Res; // I don't yet understand why this is not -0.25
						//checknormgradwidth[row][deltaind] = checkshapewidth[0][row]/width[0];
						//checknormgradcirc[row][deltaind] = checkshapecirc[row]/circ;
						//checknormgradSigma[row][deltaind] = checkshapeSigma[0][row]/Sigma[0];
						checknormgradwidth[row][deltaind] = checkshapewidth[0][row];
						checknormgradcirc[row][deltaind] = checkshapecirc[row];
						checknormgradSigma[row][deltaind] = checkshapeSigma[0][row];
						checknormgradRes[row][deltaind] = checkshapeRes[row];
					}
					delta *= (pow(10.0, steppowdelta/deltaindmax));
				}
			}
			else if (strncmp(allparams.type, "coil", 4) == 0) {
				newSigma[0] = malloc(L_fixedpoints*sizeof(double));
				newSigma[1] = malloc(L_fixedpoints*sizeof(double));
				checkshapeRes = calloc(3,sizeof(double));
				checkshapecirc = calloc(3,sizeof(double));
				checkshapewidth = malloc(L_fixedpoints*sizeof(double));
				checkshapeSigma = malloc(L_fixedpoints*sizeof(double));
				for (ind=0; ind<L_fixedpoints; ind++) {
					checkshapewidth[ind] = calloc(3,sizeof(double));
					checkshapeSigma[ind] = calloc(3,sizeof(double));
				}
				FILE *file_checkshapewidth;
				if ( (file_checkshapewidth = fopen("check_shape_width.txt", "w") ) == NULL ) {	
					printf("Oups: couldn't open check_shape_width.txt\n");
					exit(EXIT_FAILURE);
				}	
				//for (coil_ind=0; coil_ind<(*n_coils); coil_ind++) {
				for (coil_ind=0; coil_ind<max_coil_ind; coil_ind++) {
					printf("coil_ind=%d/%d\n", coil_ind, *n_coils);
					if (checfast != 'y') max_seg_ind = (*n_consts)[coil_ind];
					for ( seg_ind=0; seg_ind < max_seg_ind; seg_ind++) {
						printf("seg_ind=%d/%d\n", seg_ind, (*n_consts)[coil_ind]);
						for (row=0; row<3; row++) {
							coils[coil_ind][seg_ind][row] += (delta) ;
							island_center_up = solve_islandcenter(allparams, Bfield_island, Xp, L_fixedpoints, N_gridphi_fieldperiod);
							ext_center_up = solve_islandcenter_full(island_center_up[0].loc, axis, &num_turns, m0, L_fixedpoints, &Res, q0_index, N_gridphi_fieldperiod, Bfield_island);
							//circup = 0.0;
							newwidth[1] = calc_islandwidth(&circup, newSigma[1], ext_center_up, m0, L_fixedpoints, pol_mode);
							coils[coil_ind][seg_ind][row] -= (2.0*delta) ;
							island_center_down = solve_islandcenter(allparams, Bfield_island, Xp, L_fixedpoints, N_gridphi_fieldperiod);
							ext_center_down = solve_islandcenter_full(island_center_down[0].loc, axis, &num_turns, m0, L_fixedpoints, &Res, q0_index, N_gridphi_fieldperiod, Bfield_island);
							//circdown = 0.0;
							newwidth[0] = calc_islandwidth(&circdown, newSigma[0], ext_center_down, m0, L_fixedpoints, pol_mode);
							coils[coil_ind][seg_ind][row] += (delta) ;
							checkshapecirc[row] = ( circup - circdown ) / (2.0*delta);
							printf("checkshapecirc[%d] = %f (FD)\n", row, checkshapecirc[row]);
							checkshapecirc[row] = ( circup - circ ) / (delta);
							printf("checkshapecirc[%d] = %f (FD)\n", row, checkshapecirc[row]);
							checkshapecirc[row] = ( circ - circdown ) / (delta);
							printf("checkshapecirc[%d] = %f (FD)\n", row, checkshapecirc[row]);
							for (ind=0; ind<1; ind++) {
								checkshapeSigma[ind][row] = ( newSigma[1][ind] - newSigma[0][ind] ) / (2.0*delta);
								checkshapewidth[ind][row] = ( newwidth[1][ind] - newwidth[0][ind] ) / (2.0*delta);
								//printf("checkshapetan[%d][%d] = %f (FD)\n", ind, row, checkshapeSigma[ind][row]);   
								//printf("checkshapewidth[%d][%d] = %f (FD)\n", ind, row, checkshapewidth[ind][row]);
								//printf("shapecirc = %f (direct) %f (FD)\n", shapecirc[row], checkshapecirc[row]);
								//printf("shapetan[%d] = %f (direct) %f (FD)\n", ind, shapetan[ind][row], checkshapeSigma[ind][row]);   
								//printf("shapewidth[%d] = %f (direct) %f (FD)\n", ind, shapewidth[ind][row], checkshapewidth[ind][row]);
							}
						}
						//printf("shapecirc = %f\n", *shapecirc);
						//free(shapecirc);
						//free(shapetan[centre_ind]
						fprintf(file_checkshapewidth, "%f %f %f\n", checkshapewidth[0][0], checkshapewidth[0][1], checkshapewidth[0][2]);
					}
					fprintf(file_checkshapewidth, "\n");
				}
			}
			
		}
		clock_t int6 = clock();
		printf("time after solving for everything: %f\n", (double) (int6-start)/CLOCKS_PER_SEC);
	}

	if ( (checshap == 'y') && (strncmp(allparams.type, "Reim", 4) == 0) ) {
		printf("delta = [");
		delta = pow(10.0, powdeltain);
		for (deltaind = 0; deltaind<deltaindmax; deltaind++) {
			printf("%f ", delta);
			if (deltaind != deltaindmax - 1) printf(", ");
			delta *= (pow(10.0, steppowdelta/deltaindmax));
		}
		printf("\]\n");
		for (row=0; row<3; row++) {
			printf("gw%d = %10.10f\ngw%d_delta = \[", row, normgradwidth[row], row);
			for (deltaind = 0; deltaind<deltaindmax; deltaind++) {
				printf("%10.10f ", checknormgradwidth[row][deltaind]);
				if (deltaind != deltaindmax - 1) printf(", ");
			}
			printf("\]\n");
			printf("gC%d = %10.10f\ngC%d_delta = \[", row, normgradcirc[row], row);
			for (deltaind = 0; deltaind<deltaindmax; deltaind++) {
				printf("%10.10f ", checknormgradcirc[row][deltaind]);
				if (deltaind != deltaindmax - 1) printf(", ");
			}
			printf("\]\n");
			printf("gSigma%d = %10.10f\ngSigma%d_delta = \[", row, normgradSigma[row], row);
			for (deltaind = 0; deltaind<deltaindmax; deltaind++) {
				printf("%10.10f", checknormgradSigma[row][deltaind]);
				if (deltaind != deltaindmax - 1) printf(", ");
				
			}
			printf("\]\n");
			printf("gRes%d = %10.10f\ngRes%d_delta = \[", row, normgradRes[row], row);
			for (deltaind = 0; deltaind<deltaindmax; deltaind++) {
				printf("%10.10f", checknormgradRes[row][deltaind]);
				if (deltaind != deltaindmax - 1) printf(", ");
				
			}
			printf("\]\n");
		}
	}
	//free(Xp->tangent[0]); free(Xp->tangent[1]);
	//free(Xp->tangent);
	//free(evec[0]); free(evec[1]);
	//free(evec);
	//free(eval);
	//for (coil_ind = 0; coil_ind < *n_coils; coil_ind++) {
	//	for (seg_ind = 0; seg_ind < *n_coils; seg_ind++) {
	//		free(coils[coil_ind][seg_ind]); 
	//	}
	//	free(coils[coil_ind]); 
	//}
	//free(coils);
	//free(minor_radius);
	//free(iota);
	
	return 0;
}
