
//Currently under construction
//Author: Alessandro Geraldini
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "isc.h"
#include <gsl/gsl_sf_bessel.h>
#define error_max 1e-5

int main()
{
	clock_t start = clock();
	int search, search2;
	int N_gridphi_fieldperiod, N_gridphi_tor, num_turns=1, m0;
	int *n_coils=NULL, **n_consts, qq=1, *qq_segs=NULL, qq_segsq=1, pol_mode, tor_mode, ind, row, col;
	int stuck = 0, stuck_max = 15, count=0, allow = 0, number_minima = 0;
	int find_axis = 1, rowallowed[2][3];
	int L_fixedpoints, N_line;
	int coil_ind, seg_ind;
	int sizeline = 100;
	int q0val = 0, q02val = 0, q0val2, *q0_index=NULL, *q0_index2=NULL; // It is important to initialize the integer stored in q0_index to 0
	struct position *Xp=calloc(1,sizeof(struct position)), *axis, **mu, **lambdaQ, lambda, **lambdamu_saved, **lambdaXp, **lambdaXp2, **lambdaXp_axis;
	struct position **island_center_saved2, **lambdamu_saved2;
	struct position *Xp2=calloc(1,sizeof(struct position));
	struct position **island_center_saved, **axis_saved;
	struct ext_position *ext_center, *ext_center2, *ext_center_FD[2];
	struct field **Bfield_island, **Bfield_island2, **Bfield_axis, Bpoint;
	struct fieldparams allparams;
	//struct field *BB, *gradBB; for field checks
	double axis_optimized[2], island_center_optimized[2], island_center_optimized2[2];
	double gradXpmag, limit, saveR, dice;
	double goal = 0.0;
	double ***coils=NULL, *Sigma;  
	double **gradXp[2], **gradXp2[2], **gradXp_axis[2];
	double *shapeResheli[2], *shapeResheli2[2];
	double Res, Res2, Resopt[2], coilparams_opt[2], deltaparam[2], rmsqavRes_opt=1.0, rmsqavRes;
	double Xpguess[2], Xpguess2[2];
	double *doublepointer, **n_params=NULL;
	double phichec = 0.0, testpos[2] = {0.95, 0.0};
	double opfactor, addfactor = 1.0, tolstep;
	double **derivative_matrix;
	struct field realBpointgrad[3];
	double dr[3], *checkparam = malloc(4*sizeof(double));
	double *iota=NULL, *rmin=NULL, *zmin=NULL, axis_phi0[2];
	double **evec=malloc(2*sizeof(double*)), *eval=malloc(2*sizeof(double)), trace, det, iota_axis;
	double **matrix = malloc(2*sizeof(double)), **invertmatrix, **adjmatrix;
	double delta, checkBpointgrad[3], checknablaBpointgrad[3][2];
	double objective, objective_opt=100.0, gradobjective[2], objectiveold = 1000000000.1;
	double **Hessobjective;
	double tol_window[2];
	struct field *Bpointgrad, Bpointdelta[3][2];
	char line[sizeline];
	char poiniota, findisla, calcshap, checshap, checfast, checkmag;
	FILE *OPTIONS;

	gradXp[0]=malloc(2*sizeof(double));
	gradXp[1]=malloc(2*sizeof(double));
	gradXp2[0]=malloc(2*sizeof(double));
	gradXp2[1]=malloc(2*sizeof(double));
	gradXp_axis[0]=malloc(2*sizeof(double));
	gradXp_axis[1]=malloc(2*sizeof(double));

	matrix[0] = malloc(2*sizeof(double));
	matrix[1] = malloc(2*sizeof(double));

	q0_index = &q0val;
	q0_index2 = &q02val;

	printf("\n\n\n\n\n\n\n");

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
		// AXIS
		if (strncmp(allparams.type, "coil", 4) == 0) {
			Xp->loc[0] = 1.55; Xp->loc[1]= 0.0; //NCSX
		}
		else if (strncmp(allparams.type, "heli", 4) == 0) {
			Xp->loc[0] = 0.97; Xp->loc[1]= 0.0; //helical
			Xp->loc[0] = 0.9834425342; Xp->loc[1]= 0.0; //helical2
			//Xp->loc[0] = 0.943835; Xp->loc[1] = 0.000000; //optimized case of above
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
		printf("<--\n");
		solve_magneticaxis_save(Xp->loc, axis_saved, allparams, Bfield_axis, N_gridphi_fieldperiod, error_max);
		printf("<--\n");
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
		printstructposition("axis", axis_saved[0]);
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

	//Xp->loc[0] = 1.3; Xp->loc[1]= -0.515; // island(?) is at(1.299008, -0.515194)
	//Xp->loc[0] = 0.945; Xp->loc[1]= 0.0;//Dommaschk (5,2) amp 1.73: island is at(1.299008, -0.515194)
	//Xp->loc[0] = 1.091; Xp->loc[1]= 0.0; // Dommaschk (5,2) amp 0.00001 on top of loop current
	//Xp->loc[0] = 1.0097; Xp->loc[1]= 0.0; // Dommaschk (5,2) amp 0.00001 on top of loop current
	if (strncmp(allparams.type, "coil", 4) == 0) {
		Xpguess[0] = 1.299009; Xpguess[1] = -0.515194; // NCSX: island is at(1.299008, -0.515194)
	}
	else if (strncmp(allparams.type, "Reim", 4) == 0) {
		//Xp->loc[0] = 1.21; Xp->loc[1]= 0.0;  // Reiman 6-fold island chain
		Xpguess[0] = 1.21; Xpguess[1]= 0.0;  // Reiman 6-fold island chain
	}
	else if (strncmp(allparams.type, "heli", 4) == 0) {
		Xpguess[0] = 0.87; Xpguess[1]= 0.015; // island for I = 0.02 current 

		Xpguess[0] = 1.00910608999481 ; Xpguess[1]= -0.00000001285434;  // island for Cary Hanson case with I =0.0307
		Xpguess2[0] =  0.8243366300;  Xpguess2[1] = 0.0000000608;  // X point for Cary Hanson case with I =0.0307
		//Xpguess[0] = 1.0194035021; Xpguess[1] = 0.0000010816;// good optimization of above
		//Xpguess2[0]= 0.87173924; Xpguess2[1] = -0.00000182;// good optimization of above
		printstructposition("Xp initial guess = \n", Xp);
	}
	Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
	Xp2->tangent[0][0]=1.0; Xp2->tangent[0][1]=0.0; Xp2->tangent[1][0]=0.0; Xp2->tangent[1][1]=1.0;
	if (tor_mode % m0 == 0) 
		L_fixedpoints = pol_mode;	
	else 				 
		L_fixedpoints = m0*pol_mode;
	N_line = L_fixedpoints*N_gridphi_fieldperiod;
	Bfield_island = malloc((N_line)*sizeof(struct field));
	Bfield_island2 = malloc((N_line)*sizeof(struct field));
	// the field line is periodic but the tangent map is not
	// hence, to store both its initial and final value
	// the island_center_saved pointer is given space for one more structure
	// this is done for Bfield too just for convention (the last element here is undefined/useless)
	island_center_saved = malloc((N_line)*sizeof(struct position));
	island_center_saved2 = malloc((N_line)*sizeof(struct position));
	for (ind=0; ind<N_line; ind++) {
		Bfield_island[ind] = malloc(4*sizeof(struct field));
		island_center_saved[ind] = malloc(4*sizeof(struct position));
		Bfield_island2[ind] = malloc(4*sizeof(struct field));
		island_center_saved2[ind] = malloc(4*sizeof(struct position));
	}
	printf("Entering function that solves for the magnetic island centre (O point)-->\n");
	//island_center = solve_islandcenter(allparams, Bfield_island, Xp, L_fixedpoints, N_gridphi_fieldperiod);
	//island_center2 = solve_islandcenter(allparams, Bfield_island2, Xp2, L_fixedpoints, N_gridphi_fieldperiod);

	//for (ind=0; ind<L_fixedpoints; ind++) island_center[ind] = island_center_saved[ind*N_gridphi_fieldperiod][0];
	//solve_islandcenter_save(Xp2->loc, island_center_saved2, &Res2, allparams, Bfield_island2, L_fixedpoints, N_gridphi_fieldperiod);
	//printf("island center at position (%14.14f, %14.14f)\n", island_center_saved[0][0].loc[0], island_center_saved[0][0].loc[1]);
	//printf("<--\n");
	//solve_islandcenter_save(island_center_saved, &Res, allparams, Bfield_island, L_fixedpoints, N_gridphi_fieldperiod);
	//exit(0);
	clock_t int_afterislandcentre = clock();
	printf("... %f seconds ...\n", (double) (int_afterislandcentre-start)/CLOCKS_PER_SEC);
	num_turns=1;
	ext_center = malloc(L_fixedpoints*sizeof(struct ext_position));
	ext_center2 = malloc(L_fixedpoints*sizeof(struct ext_position));
	ext_center_FD[0] = malloc(L_fixedpoints*sizeof(struct ext_position));
	ext_center_FD[1] = malloc(L_fixedpoints*sizeof(struct ext_position));
	//ext_center[0].loc[0] = island_center[0].loc[0]; ext_center[0].loc[1] = island_center[0].loc[1];
	//ext_center2[0].loc[0] = island_center2[0].loc[0]; ext_center2[0].loc[1] = island_center2[0].loc[1];
	//for (col=0; col< L_fixedpoints; col++) {
	//	ext_center[col].q0_index = 0;
	//}
	printf("Entering function that solves for the magnetic island centre (O point) in extended way\n");
	search  = extsolve_periodicfieldline(NULL, Xpguess, ext_center, island_center_saved, Bfield_island, axis, allparams, m0, L_fixedpoints, pol_mode, q0_index, N_gridphi_fieldperiod, error_max);
	search2 = extsolve_periodicfieldline(NULL, Xpguess2, ext_center2, island_center_saved2, Bfield_island2, axis, allparams, m0, L_fixedpoints, pol_mode, q0_index, N_gridphi_fieldperiod, error_max);
	//island_center_saved[0][0].tangent[0][0] = 1.0;
	//island_center_saved[0][0].tangent[0][1] = 0.0;
	//island_center_saved[0][0].tangent[1][0] = 0.0;
	//island_center_saved[0][0].tangent[1][1] = 1.0;
	printf("Residue of fixed point = %f\n", ext_center[0].Res);
	printf("Residue of fixed point = %f\n", ext_center2[0].Res);
	//printf("Residue of fixed point = %f\n", ext_center2[0].Res);
	//exit(0);
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
	//printmat("adj_full_tangent", ext_center[0].adj_full_tangent, 2, 2);
	//printmat("full_tangent", ext_center[0].full_tangent, 2, 2);
	//printf("N_gridphi = %d, N_line = %d\n", N_gridphi_fieldperiod, N_line);
	lambdamu_saved = malloc(N_line*sizeof(struct position));
	for (ind=0; ind<N_line; ind++) lambdamu_saved[ind] = malloc(4*sizeof(struct position));
	lambdamu_saved[0][0].tangent = malloc(2*sizeof(double));
	lambdamu_saved[0][0].tangent[0] = malloc(2*sizeof(double));
	lambdamu_saved[0][0].tangent[1] = malloc(2*sizeof(double));
	lambdamu_saved[0][0].tangent[0][0] = ext_center[0].full_tangent[0][0];
	lambdamu_saved[0][0].tangent[0][1] = ext_center[0].full_tangent[1][0];
	lambdamu_saved[0][0].tangent[1][0] = ext_center[0].full_tangent[0][1];
	lambdamu_saved[0][0].tangent[1][1] = ext_center[0].full_tangent[1][1];
	solve_lambdaRes(lambdamu_saved, island_center_saved, m0, L_fixedpoints, N_gridphi_fieldperiod, Bfield_island);
	lambdamu_saved2 = malloc(N_line*sizeof(struct position));
	for (ind=0; ind<N_line; ind++) lambdamu_saved2[ind] = malloc(4*sizeof(struct position));
	lambdamu_saved2[0][0].tangent = malloc(2*sizeof(double));
	lambdamu_saved2[0][0].tangent[0] = malloc(2*sizeof(double));
	lambdamu_saved2[0][0].tangent[1] = malloc(2*sizeof(double));
	lambdamu_saved2[0][0].tangent[0][0] = ext_center2[0].full_tangent[0][0];
	lambdamu_saved2[0][0].tangent[0][1] = ext_center2[0].full_tangent[1][0];
	lambdamu_saved2[0][0].tangent[1][0] = ext_center2[0].full_tangent[0][1];
	lambdamu_saved2[0][0].tangent[1][1] = ext_center2[0].full_tangent[1][1];
	printmat("full_tangent", ext_center[0].full_tangent, 2, 2);
	solve_lambdaRes(lambdamu_saved2, island_center_saved2, m0, L_fixedpoints, N_gridphi_fieldperiod, Bfield_island2);
	clock_t int_afteralladjoints = clock();


	derivative_matrix = set_identity();
	Hessobjective = set_identity();
	gradobjective[0] = 0.0;
	gradobjective[1] = 0.0;
	for (coil_ind=0; coil_ind < allparams.n_coils; coil_ind++) {
		gradXp_axis[coil_ind][0] = calloc(allparams.n_diff,sizeof(double));
		gradXp_axis[coil_ind][1] = calloc(allparams.n_diff,sizeof(double));
		gradXp[coil_ind][0] = calloc(allparams.n_diff,sizeof(double));
		gradXp[coil_ind][1] = calloc(allparams.n_diff,sizeof(double));
		gradXp2[coil_ind][0] = calloc(allparams.n_diff,sizeof(double));
		gradXp2[coil_ind][1] = calloc(allparams.n_diff,sizeof(double));
	}
	lambdaXp = malloc(N_line*sizeof(struct position));
	lambdaXp2 = malloc(N_line*sizeof(struct position));
	lambdaXp_axis = malloc(N_gridphi_fieldperiod*sizeof(struct position));
	for (ind = 0; ind< N_line; ind++) {
		lambdaXp[ind] = malloc(4*sizeof(struct position));
		lambdaXp2[ind] = malloc(4*sizeof(struct position));
	}
	for (ind = 0; ind < N_gridphi_fieldperiod; ind++) lambdaXp_axis[ind] = malloc(4*sizeof(struct position));
	lambdaXp[0][0].tangent = set_identity();
	lambdaXp2[0][0].tangent = set_identity();
	lambdaXp_axis[0][0].tangent = set_identity();
	printf("... %f seconds ...\n", (double) (int_afteralladjoints-start)/CLOCKS_PER_SEC);
	//printf("Entering optimization loop\n");
	if (strncmp("heli", allparams.type, 4) == 0) {
		rowallowed[0][0] = 2; rowallowed[0][1] = 4; rowallowed[0][2] = 6;
		rowallowed[1][0] = 1; rowallowed[1][1] = 4; rowallowed[1][2] = 5;
		rowallowed[0][0] = 2;
		rowallowed[1][0] = 1;
		L_fixedpoints = 3;
	}
	shapeResheli[0] = calloc(allparams.n_diff,sizeof(double));
	shapeResheli[1] = calloc(allparams.n_diff,sizeof(double));
	shapeResheli2[0] = calloc(allparams.n_diff,sizeof(double));
	shapeResheli2[1] = calloc(allparams.n_diff,sizeof(double));
	coil_ind = 0;
	//tol_window[0] = 0.2*tolstep;
	tol_window[1] = 100.0;
	Xp->loc[0] = axis_saved[0][0].loc[0];
	Xp->loc[1] = axis_saved[0][0].loc[1];
	tolstep = 0.001;
	do {
		stuck = 0;
		do {
		q0val = 0;
		q0val2 = 0;
		printf("coil parameters = (%14.14f, %14.14f)\n", allparams.diffparams[0][0][rowallowed[0][0]], allparams.diffparams[1][0][rowallowed[1][0]]);
		printf("axis_guess = (%f, %f)\n", Xp->loc[0], Xp->loc[1]);
		printf("island_center_guess = (%10f, %10f)\n", Xpguess[0], Xpguess[1]);
		printf("island_center2_guess = (%10f, %10f)\n",Xpguess2[0], Xpguess2[1]);
		solve_magneticaxis_save(Xp->loc, axis_saved, allparams, Bfield_axis, N_gridphi_fieldperiod, error_max);
		ind = 0;
		saveR = Xp->loc[0];
		while ( fabs(axis_saved[0][0].loc[0] - Xp->loc[0]) > tolstep ) {
			Xp->loc[0] = saveR + (ind/tol_window[1])*tolstep;
			solve_magneticaxis_save(Xp->loc, axis_saved, allparams, Bfield_axis, N_gridphi_fieldperiod, error_max);
			if ( fabs(island_center_saved[0][0].loc[0] - Xpguess[0]) > tolstep ) {
				Xp->loc[0] = saveR - (ind/tol_window[1])*tolstep;
				solve_magneticaxis_save(Xp->loc, axis_saved, allparams, Bfield_axis, N_gridphi_fieldperiod, error_max);
			}
			ind++;
			if (ind>tol_window[1]) {
				printf("ERROR: magnetic axis lost during iteration\n");
				exit(1);
			}
		}
		printf("axis = (%f, %f)\n", axis_saved[0][0].loc[0], axis_saved[0][0].loc[1]);
		for (ind=0; ind<N_gridphi_fieldperiod; ind++) axis[ind] = axis_saved[ind][0];
		extsolve_periodicfieldline(NULL, Xpguess, ext_center, island_center_saved, Bfield_island, axis, allparams, m0, L_fixedpoints, pol_mode, &q0val, N_gridphi_fieldperiod, error_max);
		ind = 0;
		saveR = Xpguess[0];
		while ( fabs(island_center_saved[0][0].loc[0] - Xpguess[0]) > tolstep ) {
			q0val = 0;
			Xpguess[0] = saveR + (ind/tol_window[1])*tolstep;
			extsolve_periodicfieldline(NULL, Xpguess, ext_center, island_center_saved, Bfield_island, axis, allparams, m0, L_fixedpoints, pol_mode, &q0val, N_gridphi_fieldperiod, error_max);
			if ( fabs(island_center_saved[0][0].loc[0] - Xpguess[0]) > tolstep ) {
				Xpguess[0] = saveR - (ind/tol_window[1])*tolstep;
				extsolve_periodicfieldline(NULL, Xpguess, ext_center, island_center_saved, Bfield_island, axis, allparams, m0, L_fixedpoints, pol_mode, &q0val, N_gridphi_fieldperiod, error_max);
			}
			ind++;
			if (ind>tol_window[1]) {
				printf("ERROR: fixed point lost during iteration\n");
				exit(1);
			}
		}
		printf("island_center  = (%10f, %10f)\n", island_center_saved[0][0].loc[0], island_center_saved[0][0].loc[1]);
		extsolve_periodicfieldline(NULL, Xpguess2, ext_center2, island_center_saved2, Bfield_island2, axis, allparams, m0, L_fixedpoints, pol_mode, &q0val2, N_gridphi_fieldperiod, error_max);
		ind = 0;
		saveR = Xpguess2[0];
		while ( fabs(island_center_saved2[0][0].loc[0] - Xpguess2[0]) > tolstep ) {
			q0val2 = 0;
			Xpguess2[0] = saveR + (ind/tol_window[1])*tolstep;
			extsolve_periodicfieldline(NULL, Xpguess2, ext_center2, island_center_saved2, Bfield_island2, axis, allparams, m0, L_fixedpoints, pol_mode, &q0val2, N_gridphi_fieldperiod, error_max);
			if ( fabs(island_center_saved[0][0].loc[0] - Xpguess[0]) > tolstep ) {
				Xpguess[0] = saveR - (ind/tol_window[1])*tolstep;
				extsolve_periodicfieldline(NULL, Xpguess2, ext_center2, island_center_saved2, Bfield_island2, axis, allparams, m0, L_fixedpoints, pol_mode, &q0val2, N_gridphi_fieldperiod, error_max);
			}
			ind++;
			if (ind>tol_window[1]) {
				printf("ERROR: fixed point lost during iteration\n");
				exit(1);
			}
		}
		printf("island_center2 = (%10f, %10f)\n", island_center_saved2[0][0].loc[0], island_center_saved2[0][0].loc[1]);
		if (fabs(island_center_saved[0][0].loc[0] - axis_saved[0][0].loc[0]) < fabs(island_center_saved2[0][0].loc[0] - axis_saved[0][0].loc[0]) )
		tolstep = 0.3*fabs(island_center_saved[0][0].loc[0] - axis_saved[0][0].loc[0]);
		else 
		tolstep = 0.3*fabs(island_center_saved2[0][0].loc[0] - axis_saved[0][0].loc[0]);
		//tolstep = 0.005;
		tol_window[0] = tolstep;
		//printf("tolstep=%f\n", tolstep);
		//tolstep = 0.01*sqrt( pow(island_center_saved2[0][0].loc[0] - axis_saved[0][0].loc[0], 2.0) + pow(island_center_saved[0][0].loc[0] - axis_saved[0][0].loc[0], 2.0) ) ;
		//tol_window[0] = tolstep;
		Res = ext_center[0].Res;
		Res2 = ext_center2[0].Res;
		printf("Res = (%f, %f)\n", Res, Res2);
		lambdamu_saved[0][0].tangent[0][0] = ext_center[0].full_tangent[0][0];
		lambdamu_saved[0][0].tangent[0][1] = ext_center[0].full_tangent[1][0];
		lambdamu_saved[0][0].tangent[1][0] = ext_center[0].full_tangent[0][1];
		lambdamu_saved[0][0].tangent[1][1] = ext_center[0].full_tangent[1][1];
		solve_lambdaRes(lambdamu_saved, island_center_saved, m0, L_fixedpoints, N_gridphi_fieldperiod, Bfield_island);
		lambdamu_saved2[0][0].tangent[0][0] = ext_center2[0].full_tangent[0][0];
		lambdamu_saved2[0][0].tangent[0][1] = ext_center2[0].full_tangent[1][0];
		lambdamu_saved2[0][0].tangent[1][0] = ext_center2[0].full_tangent[0][1];
		lambdamu_saved2[0][0].tangent[1][1] = ext_center2[0].full_tangent[1][1];
		solve_lambdaRes(lambdamu_saved2, island_center_saved2, m0, L_fixedpoints, N_gridphi_fieldperiod, Bfield_island2);

		matrix[0][0] = 1.0 -  ext_center[0].adj_full_tangent[0][0];
		matrix[0][1] =     -  ext_center[0].adj_full_tangent[0][1];
		matrix[1][0] =     -  ext_center[0].adj_full_tangent[1][0];
		matrix[1][1] = 1.0 -  ext_center[0].adj_full_tangent[1][1];
		invertmatrix = invert2x2(matrix, &det);
		lambdaXp[0][0].tangent[0][0] = invertmatrix[0][0];
		lambdaXp[0][0].tangent[0][1] = invertmatrix[0][1];
		lambdaXp[0][0].tangent[1][0] = invertmatrix[1][0];
		lambdaXp[0][0].tangent[1][1] = invertmatrix[1][1];
		solve_lambdaXp(lambdaXp, island_center_saved, m0, L_fixedpoints, N_gridphi_fieldperiod, Bfield_island); 
		free(invertmatrix[0]); free(invertmatrix[1]);
		free(invertmatrix);

		matrix[0][0] = 1.0 -  ext_center2[0].adj_full_tangent[0][0];
		matrix[0][1] =     -  ext_center2[0].adj_full_tangent[0][1];
		matrix[1][0] =     -  ext_center2[0].adj_full_tangent[1][0];
		matrix[1][1] = 1.0 -  ext_center2[0].adj_full_tangent[1][1];
		invertmatrix = invert2x2(matrix, &det);
		lambdaXp2[0][0].tangent[0][0] = invertmatrix[0][0];
		lambdaXp2[0][0].tangent[0][1] = invertmatrix[0][1];
		lambdaXp2[0][0].tangent[1][0] = invertmatrix[1][0];
		lambdaXp2[0][0].tangent[1][1] = invertmatrix[1][1];
		solve_lambdaXp(lambdaXp2, island_center_saved2, m0, L_fixedpoints, N_gridphi_fieldperiod, Bfield_island2); 
		free(invertmatrix[0]); free(invertmatrix[1]);
		free(invertmatrix);

		adjmatrix = adj2x2(axis_saved[0][0].tangent, &det);
		matrix[0][0] = 1.0 -  adjmatrix[0][0];
		matrix[0][1] =     -  adjmatrix[0][1];
		matrix[1][0] =     -  adjmatrix[1][0];
		matrix[1][1] = 1.0 -  adjmatrix[1][1];
		invertmatrix = invert2x2(matrix, &det);
		lambdaXp_axis[0][0].tangent[0][0] = invertmatrix[0][0];
		lambdaXp_axis[0][0].tangent[0][1] = invertmatrix[0][1];
		lambdaXp_axis[0][0].tangent[1][0] = invertmatrix[1][0];
		lambdaXp_axis[0][0].tangent[1][1] = invertmatrix[1][1];
		solve_lambdaXp(lambdaXp_axis, axis_saved, m0, 1, N_gridphi_fieldperiod, Bfield_axis); 
		free(invertmatrix[0]); free(invertmatrix[1]);
		free(invertmatrix);
		free(adjmatrix[0]); free(adjmatrix[1]);
		free(adjmatrix);

		opfactor = 0.5;
		rmsqavRes = sqrt(0.5*Res*Res+0.5*Res2*Res2);
		objective = rmsqavRes;
		objective = (fabs(Res) + fabs(Res2))/2.0;
		objective = (1.0/(1.0+fabs(addfactor)))*(Res + addfactor*Res2);
		printf("objective = %f\n", objective);
		printf("objectiveold = %f\n", objectiveold);
		printf("rmsqavRes = %f\n", rmsqavRes);
		if ( rmsqavRes  < rmsqavRes_opt ) {
			Resopt[0] = Res;	
			Resopt[1] = Res2;	
			objective_opt = objective;
			coilparams_opt[0] = allparams.diffparams[0][0][rowallowed[0][0]];
			coilparams_opt[1] = allparams.diffparams[1][0][rowallowed[1][0]];
			rmsqavRes_opt = rmsqavRes;
			axis_optimized[0] = axis_saved[0][0].loc[0];
			axis_optimized[1] = axis_saved[0][0].loc[1];
			island_center_optimized[0] = island_center_saved[0][0].loc[0];
			island_center_optimized[1] = island_center_saved[0][0].loc[1];
			island_center_optimized2[0] = island_center_saved2[0][0].loc[0];
			island_center_optimized2[1] = island_center_saved2[0][0].loc[1];
		}

		if ( ( ( (count != 0) && (objectiveold  < objective) ) && (allow == 0) ) && (stuck < stuck_max) ) { //- gradobjective[0]*deltaparam[0] - gradobjective[1]*deltaparam[1]
			Xpguess[0] = island_center_saved[0][0].loc[0];			
			Xpguess2[0] = island_center_saved2[0][0].loc[0];
			Xp->loc[0] = axis_saved[0][0].loc[0];
			Xpguess[1] = 0.0;
			Xpguess2[1] = 0.0;
			Xp->loc[1] = 0.0;
			for (coil_ind=0;coil_ind<2; coil_ind++) {
				//printf("deltaparam[%d] = %f\n", coil_ind, deltaparam[coil_ind]);
				deltaparam[coil_ind] /= 2.0;
				allparams.diffparams[coil_ind][0][rowallowed[coil_ind][0]] += deltaparam[coil_ind];
				Xpguess[0] += gradXp[coil_ind][0][rowallowed[coil_ind][0]]*deltaparam[coil_ind];
				Xpguess2[0] += gradXp2[coil_ind][0][rowallowed[coil_ind][0]]*deltaparam[coil_ind];
				Xp->loc[0] += gradXp_axis[coil_ind][0][rowallowed[coil_ind][0]]*deltaparam[coil_ind];
			}
			stuck += 1;
			//printf("stuck = %d\n", stuck);
		}
		//////if (count == 0) { gradobjective[0] = 0.0; gradobjective[1]= 0.0;}
		//printf("allow = %d\n", allow);
		
		} while ( ( ( (count != 0) && (objectiveold  < objective) ) && (allow == 0) ) && (stuck < stuck_max) ); //gradobjective[0]*2.0*deltaparam[0] - gradobjective[1]*2.0*deltaparam[1]
		//
		if ( (allow == 1) || (stuck < stuck_max )  ) {
		stuck = 0;
		allow = 0;
		for (coil_ind=0;coil_ind<2; coil_ind++) {

			
			solve_gradRes(shapeResheli[coil_ind], Bfield_island, island_center_saved, lambdamu_saved, L_fixedpoints, N_gridphi_fieldperiod, allparams, coil_ind, 0) ;
			solve_gradRes(shapeResheli2[coil_ind], Bfield_island2, island_center_saved2, lambdamu_saved2, L_fixedpoints, N_gridphi_fieldperiod, allparams, coil_ind, 0) ;
			solve_gradXp(gradXp[coil_ind], Bfield_island, island_center_saved, lambdaXp, L_fixedpoints, N_gridphi_fieldperiod, allparams, coil_ind, 0) ;
			solve_gradXp(gradXp2[coil_ind], Bfield_island2, island_center_saved2, lambdaXp2, L_fixedpoints, N_gridphi_fieldperiod, allparams, coil_ind, 0) ;
			solve_gradXp(gradXp_axis[coil_ind], Bfield_axis, axis_saved, lambdaXp_axis, 1, N_gridphi_fieldperiod, allparams, coil_ind, 0) ;
			gradobjective[coil_ind] = 0.5*0.25*(Res*shapeResheli[coil_ind][rowallowed[coil_ind][0]] + Res2*shapeResheli2[coil_ind][rowallowed[coil_ind][0]])/objective;
			gradobjective[coil_ind] = 0.5*0.25*((fabs(Res)/Res)*shapeResheli[coil_ind][rowallowed[coil_ind][0]] + (fabs(Res2)/Res2)*shapeResheli2[coil_ind][rowallowed[coil_ind][0]]);
			gradobjective[coil_ind] = (1.0/(1.0+fabs(addfactor)))*0.25*(shapeResheli[coil_ind][rowallowed[coil_ind][0]] + addfactor*shapeResheli2[coil_ind][rowallowed[coil_ind][0]]);
			//gradobjective[coil_ind] = 4.0*0.25*(Res*Res*Res*shapeResheli[coil_ind][rowallowed[coil_ind][0]] + Res2*Res2*Res2*shapeResheli2[coil_ind][rowallowed[coil_ind][0]]);
			derivative_matrix[0][coil_ind] = shapeResheli[coil_ind][rowallowed[coil_ind][0]];
			derivative_matrix[1][coil_ind] = shapeResheli2[coil_ind][rowallowed[coil_ind][0]];
			//printf("gradobjective[%d] = %f\n", coil_ind, gradobjective[coil_ind]);
			//printf("shapeRes1[%d] = %f\n", coil_ind, shapeResheli[coil_ind][rowallowed[coil_ind][0]]);
			//printf("shapeRes1[%d] = %f\n", coil_ind, shapeResheli2[coil_ind][rowallowed[coil_ind][0]]);


			//dice = (double) rand()/RAND_MAX;
			//if ( count % 10 == 0 ) {
			//	randomnumber = (double) rand()/RAND_MAX;
			//	gradobjective[coil_ind] *= 2.0*(randomnumber - 0.5);
			//	printf("randomised gradobjective[%d] = %f\n", coil_ind, gradobjective[coil_ind]);
			//	//addfactor = (double) rand()/RAND_MAX;//float in range 0 to 1
			//	//printf("randomised weight = %f\n", addfactor);
			//}

			gradXpmag = sqrt( pow(gradXp[coil_ind][0][rowallowed[coil_ind][0]], 2.0) + pow(gradXp[coil_ind][1][rowallowed[coil_ind][0]], 2.0) );
			if ( gradXpmag < sqrt( pow(gradXp2[coil_ind][0][rowallowed[coil_ind][0]], 2.0) + pow(gradXp2[coil_ind][1][rowallowed[coil_ind][0]], 2.0) ) )  
				gradXpmag = sqrt( pow(gradXp2[coil_ind][0][rowallowed[coil_ind][0]], 2.0) + pow(gradXp2[coil_ind][1][rowallowed[coil_ind][0]], 2.0) );
			if ( gradXpmag < sqrt( pow(gradXp_axis[coil_ind][0][rowallowed[coil_ind][0]], 2.0) + pow(gradXp_axis[coil_ind][1][rowallowed[coil_ind][0]], 2.0) ) )  
				gradXpmag = sqrt( pow(gradXp_axis[coil_ind][0][rowallowed[coil_ind][0]], 2.0) + pow(gradXp_axis[coil_ind][1][rowallowed[coil_ind][0]], 2.0) );
			limit = (tolstep/gradXpmag)*fabs(gradobjective[coil_ind]/objective);

			if (opfactor > limit) opfactor = limit;
			//printf("opfactor= %f, limit = %f\n", opfactor, limit);
		}
		//stepdirection[0] = gradobjective[0]*random_direction[0]
		//opfactor = 0.5;
		//invertmatrix = invert2x2(derivative_matrix, &det);
		//if (fabs(gradXp[0][0][rowallowed[0][0]]) > fabs(gradXp[1][0][rowallowed[1][0]])) gradXpmag = gradXp[0][0][rowallowed[0][0]];
		//else  gradXpmag = gradXp[1][0][rowallowed[1][0]];
		//limit = fabs((tolstep/(gradXpmag*(invertmatrix[0][0]*Res + invertmatrix[0][1]*Res2))));
		//if (opfactor > limit) opfactor = limit;
		//printf("objfactor= %f, limit = %f\n", opfactor, limit);
		//opfactor = 1.0;
		//deltaparam[1] = opfactor*(invertmatrix[1][0]*Res + invertmatrix[1][1]*Res2);
		Xpguess[0] = island_center_saved[0][0].loc[0];			
		Xpguess2[0] = island_center_saved2[0][0].loc[0];
		Xp->loc[0] = axis_saved[0][0].loc[0];
		Xpguess[1] = 0.0;
		Xpguess2[1] = 0.0;
		Xp->loc[1] = 0.0;
		for (coil_ind=0;coil_ind<2; coil_ind++) {
			deltaparam[coil_ind] = opfactor*(objective/gradobjective[coil_ind]); 
			//deltaparam[coil_ind] = opfactor*(invertmatrix[coil_ind][0]*Res + invertmatrix[coil_ind][1]*Res2);
			//printf("deltaparam[%d] = %f\n", coil_ind, deltaparam[coil_ind]);
			allparams.diffparams[coil_ind][0][rowallowed[coil_ind][0]] -= deltaparam[coil_ind];
			Xpguess[0] -= gradXp[coil_ind][0][rowallowed[coil_ind][0]]*deltaparam[coil_ind];
			Xpguess2[0] -= gradXp2[coil_ind][0][rowallowed[coil_ind][0]]*deltaparam[coil_ind];
			Xp->loc[0] -= gradXp_axis[coil_ind][0][rowallowed[coil_ind][0]]*deltaparam[coil_ind];
		}
		objectiveold = objective;		
		}
		//free(invertmatrix);
		//if ( count % 50 == 0 ) {
		//	printf("optimized coil parameters = (%f, %f)\n", coilparams_opt[0], coilparams_opt[1]);
		//	printf("optimized Residues = (%f, %f)\n", Resopt[0], Resopt[1]);
		//	printf("optimized average of magnitude of residues = %f\n", objective_opt);
		//	printf("optimized rmsq average of residues = %f\n", rmsqavRes_opt);
		//}
		else {
			stuck = 0;
			allow = 1;
			number_minima += 1;
			if (number_minima == 30) count = 200;
			printf("# MINIMA = %d\n", number_minima);
			printf("optimized coil parameters = (%f, %f)\n", coilparams_opt[0], coilparams_opt[1]);
			printf("optimized Residues = (%f, %f)\n", Resopt[0], Resopt[1]);
			printf("optimized average of magnitude of residues = %f\n", objective_opt);
			printf("optimized rmsq average of residues = %f\n", rmsqavRes_opt);
		}
		count++;
	} while (fabs(objective) > goal && count < 200); 
	printf("optimized coil parameters = (%f, %f)\n", coilparams_opt[0], coilparams_opt[1]);
	printf("optimized Residues = (%f, %f)\n", Resopt[0], Resopt[1]);
	printf("optimized average of magnitude of residues = %f\n", objective_opt);
	printf("optimized rmsq average of residues = %f\n", rmsqavRes_opt);
	//printf("# iterations = %d\n", count);
	clock_t intopt = clock();
	printf("time after solving for everything: %f\n", (double) (intopt-start)/CLOCKS_PER_SEC);
		//exit(0);
	return 0;
}
