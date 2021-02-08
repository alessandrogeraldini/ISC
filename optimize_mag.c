//Author: Alessandro Geraldini
//Currently optimizer does not deal with explicit "coil" magnetic field type

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include "isc.h"
#include <gsl/gsl_sf_bessel.h>
#define error_max 1e-14
#define N_gridphi_fieldperiod 50

int main()
{
    clock_t start = clock();
    struct fieldparams allparams;
    int lenline = 1000, num_line_file, num_resonances;
    int *search, m0, *L_fixedpoints, N_line;
    int *q0_index; // It is important to initialize the integer stored in q0_index to 0
    int *pol_mode, *tor_mode, fixed_ind, res_ind, ind, coil_ind, param_ind;
    int stuck = 0, stuck_int = 4, stuck_max = 10, count=0, count_int = 0, count_max = 10000, count_since_optimum=0, allow = 0, number_minima = 0, rowallowed[2][3];
    int sizeline = 100, line_elements;
    int max_minima = 100;
    struct position *axis, ***lambdamu_saved, ***lambdaXp, **lambdaXp_axis;
    struct position ***periodicfield_saved, **axis_saved;
    struct ext_position **ext_periodicfield;
    struct field ***Bfield_island, **Bfield_axis, Bpoint;
    //struct field *BB, *gradBB; for field checks
    double stop_threshold = 1e-14;
    double **fixedpoint_RZ, *fixedpoint_info, axis_RZ[2], r_decrease = 2.0;
    double axis_optimized[2], **periodicfield_optimized;
    double *gradXpmag, numb, limit, saveR, det;
    double goal = 1e-10, reference_change=10000000.0, fracfpsep = 0.3;
    double ****gradXp, **gradXp_axis[2];
    double ***shapeResheli;
    double *Res, *Resopt, coilparams_opt[2], deltaparam[2], rmsqavRes_opt=100.0, rmsqavRes;
    double **Xpguess, axisguess[2];
    double testpos[2] = {0.95, 0.0};
    double opfactor_reference = 1.0, opfactor=opfactor_reference, tolstep;
    double **matrix = malloc(2*sizeof(double)), **invertmatrix, **adjmatrix;
    double objective=0.0, objective_opt=9e13, gradobjective[2], objectiveold = 9e13;
    double gradobjectivesquared;
    double tol_window[2];
    double **optparams, *factor, sumfactor = 0.0;
    char line[sizeline];

    FILE *fixed_points;
    if ( (fixed_points = fopen("fixedpoints.txt", "r")) == NULL) {
	printf("ERROR: file fixedpoints.txt didn't open\n");
	exit(1);
    }
    num_line_file = 0;
    while (fgets(line, lenline, fixed_points) != NULL) {	
	printf("num_line_file at start= %d\n", num_line_file);
	num_line_file+=1;
	printf("%s", line);
	printf("num_line_file at end = %d\n", num_line_file);
    }
    rewind(fixed_points);
    num_resonances = num_line_file-1;
    printf("number of resonances analyzed = %d\n", num_resonances);
    
    factor=malloc(num_resonances*sizeof(double));
    gradXp=malloc(num_resonances*sizeof(double));
    shapeResheli=malloc(num_resonances*sizeof(double));
    lambdaXp=malloc(num_resonances*sizeof(double));
    lambdamu_saved=malloc(num_resonances*sizeof(double));
    gradXp_axis[0]=malloc(2*sizeof(double));
    gradXp_axis[1]=malloc(2*sizeof(double));
    
    matrix[0] = malloc(2*sizeof(double));
    matrix[1] = malloc(2*sizeof(double));
    

    
    fixedpoint_RZ = malloc(num_resonances*sizeof(double));
    Xpguess = malloc(num_resonances*sizeof(double));
    ext_periodicfield = malloc(num_resonances*sizeof(struct ext_position));
    periodicfield_optimized = malloc(num_resonances*sizeof(double));
    Bfield_island = malloc(num_resonances*sizeof(struct field));
    periodicfield_saved = malloc(num_resonances*sizeof(struct position));
    q0_index = calloc(num_resonances, sizeof(int));
    tor_mode = malloc(num_resonances*sizeof(int));
    pol_mode = malloc(num_resonances*sizeof(int));
    L_fixedpoints = malloc(num_resonances*sizeof(int));
    search = malloc(num_resonances*sizeof(int));
    Res = malloc(num_resonances*sizeof(double));
    Resopt = malloc(num_resonances*sizeof(double));
    num_line_file = 0;
    line_elements = 4;
    while (fgets(line, lenline, fixed_points) != NULL) {	
	line_elements = 4;
	fixedpoint_info = linetodata(line, &line_elements);	
	if (num_line_file == 0) {
	    axis_RZ[0] = fixedpoint_info[0];
	    axis_RZ[1] = fixedpoint_info[1];
	    printf("(R, Z) = (%f, %f) for magnetic axis\n", axis_RZ[0], axis_RZ[1]);
	}
	else {
	    fixedpoint_RZ[num_line_file-1] = malloc(2*sizeof(double));
	    fixedpoint_RZ[num_line_file-1][0] = fixedpoint_info[0];
	    fixedpoint_RZ[num_line_file-1][1] = fixedpoint_info[1];
	    pol_mode[num_line_file-1] = (int) fixedpoint_info[2];
	    tor_mode[num_line_file-1] = (int) fixedpoint_info[3];
	    printf("(R, Z) = (%f, %f) for pol_mode = %d, tor_mode = %d\n", fixedpoint_RZ[num_line_file-1][0], fixedpoint_RZ[num_line_file-1][1], pol_mode[num_line_file-1], tor_mode[num_line_file-1]);
	}
	num_line_file += 1;
	free(fixedpoint_info);
    }
    fclose(fixed_points);
    

    /* 
    Read parameters of magnetic field:
    coils contains the magnetic field parameters
    for magnetic fields produced by coils these are coil segment positions
    for Reiman configurations coils contains the parameters 
    (iota0, iota1, epsilon_i) for i =(0, ..., N)
    */
    allparams = fetchparams();
    m0 = allparams.m0_fieldperiods;
    Bpoint = Bfield(testpos, 0.0, allparams);
    printf("(BR, BZ, Bphi) = (%f, %f, %f)\n", Bpoint.value[0], Bpoint.value[1], Bpoint.value[2]);
    

    Bfield_axis = malloc(N_gridphi_fieldperiod*sizeof(struct field));
    for (ind=0; ind<N_gridphi_fieldperiod; ind++) {
    	Bfield_axis[ind] = malloc(4*sizeof(struct field));
    }
    printf("Entering function that solves for the magnetic axis-->\n");
    axis_saved = malloc(N_gridphi_fieldperiod*sizeof(struct position));
    for (ind=0; ind<N_gridphi_fieldperiod; ind++) axis_saved[ind] = malloc(4*sizeof(struct position));
    solve_magneticaxis_save(axis_RZ, axis_saved, allparams, Bfield_axis, N_gridphi_fieldperiod, error_max);
    axis = malloc(N_gridphi_fieldperiod*sizeof(struct position));
    for (ind=0; ind<N_gridphi_fieldperiod; ind++) axis[ind] = axis_saved[ind][0];

    //adjXp = adj2x2(Xp->tangent, &det);
    //lambdamu_saved = solve_lambdaRes(Bfield_axis, axis_saved, adjXp, m0, 1, N_gridphi_fieldperiod);
    printf("<--\n");
    clock_t t_afteraxis = clock();
    printf("... %f seconds ...\n", (double) (t_afteraxis-start)/CLOCKS_PER_SEC);
    printf("m0_fieldperiods = %d\n", m0);
    //axis_phi0[0] = axis[0].loc[0];  axis_phi0[1] = axis[0].loc[1];
    printstructposition("axis", axis_saved[0]);
    //linalg2x2(axis_saved[0][0].tangent, evec, eval, &det, &trace);
    ////printf("evalss=(%f, %f)\n", eval[0], eval[1]);
    //printstructposition("axis", axis_saved[0]);
    //iota_axis = -eval[1]*m0/(2*M_PI);
    //printf("iota_axis=%14.14f\n", iota_axis);

    gradXpmag = malloc(allparams.n_coils*sizeof(double));
    for (coil_ind=0; coil_ind < allparams.n_coils; coil_ind++) {
	gradXp_axis[coil_ind] = malloc(2*sizeof(double));
	gradXp_axis[coil_ind][0] = calloc(allparams.n_diff,sizeof(double));
	gradXp_axis[coil_ind][1] = calloc(allparams.n_diff,sizeof(double));
    }
    for (res_ind=0; res_ind< num_resonances; res_ind++) {
	factor[res_ind] = 1.0;
	//if (res_ind == 1) factor[res_ind] = -1.0;
        periodicfield_optimized[res_ind] = malloc(2*sizeof(double));
        gradXp[res_ind]=malloc(allparams.n_coils*sizeof(double));
	shapeResheli[res_ind] = malloc(allparams.n_coils*sizeof(double));
	Xpguess[res_ind] = malloc(2*sizeof(double));
    	if (tor_mode[res_ind] % m0 == 0) L_fixedpoints[res_ind] = pol_mode[res_ind];	
    	else L_fixedpoints[res_ind] = m0*pol_mode[res_ind];
    	N_line = L_fixedpoints[res_ind]*N_gridphi_fieldperiod;
    	Bfield_island[res_ind] = malloc(N_line*sizeof(struct field));
    	periodicfield_saved[res_ind] = malloc(N_line*sizeof(struct position));
    	for (ind=0; ind<N_line; ind++) {
	    Bfield_island[res_ind][ind] = malloc(4*sizeof(struct field));
	    periodicfield_saved[res_ind][ind] = malloc(4*sizeof(struct position));
    	}
    	ext_periodicfield[res_ind] = malloc(L_fixedpoints[res_ind]*sizeof(struct ext_position));
    	printf("Entering function that solves for the periodic field line in extended way\n");
    	search[res_ind]  = extsolve_periodicfieldline(NULL, fixedpoint_RZ[res_ind], ext_periodicfield[res_ind], periodicfield_saved[res_ind], Bfield_island[res_ind], axis, allparams, m0, L_fixedpoints[res_ind], pol_mode[res_ind], q0_index+res_ind, N_gridphi_fieldperiod, error_max);
    	if (search[res_ind] == 0) {
	    //periodicfield_saved[0][0].tangent[0][0] = 1.0;
	    //periodicfield_saved[0][0].tangent[0][1] = 0.0;
	    //periodicfield_saved[0][0].tangent[1][0] = 0.0;
	    //periodicfield_saved[0][0].tangent[1][1] = 1.0;
	    printf("List of fixed points and their basic info: \n");
	    for (fixed_ind = 0; fixed_ind<L_fixedpoints[res_ind]; fixed_ind++) {
	    	printf("periodicfield[%d].(R,Z) = (%f, %f)\n", fixed_ind, ext_periodicfield[res_ind][fixed_ind].loc[0], ext_periodicfield[res_ind][fixed_ind].loc[1]);
	    	printf("periodicfield[%d].(Res, width, Sigma, circ) = (%f, %f, %f, %f)\n", fixed_ind, ext_periodicfield[res_ind][fixed_ind].Res, ext_periodicfield[res_ind][fixed_ind].width, ext_periodicfield[res_ind][fixed_ind].Sigma, ext_periodicfield[res_ind][fixed_ind].circ);
	    }
	    clock_t int_afterislandcentre = clock();
	    printf("... %f seconds ...\n", (double) (int_afterislandcentre-start)/CLOCKS_PER_SEC);



	    printf("Residue of fixed point = %f\n", ext_periodicfield[res_ind][0].Res);
	}

	lambdamu_saved[res_ind] = malloc(N_line*sizeof(struct position));
	for (ind=0; ind<N_line; ind++) lambdamu_saved[res_ind][ind] = malloc(4*sizeof(struct position));
	lambdamu_saved[res_ind][0][0].tangent = malloc(2*sizeof(double));
	lambdamu_saved[res_ind][0][0].tangent[0] = malloc(2*sizeof(double));
	lambdamu_saved[res_ind][0][0].tangent[1] = malloc(2*sizeof(double));
	lambdamu_saved[res_ind][0][0].tangent[0][0] = ext_periodicfield[res_ind][0].full_tangent[0][0];
	lambdamu_saved[res_ind][0][0].tangent[0][1] = ext_periodicfield[res_ind][0].full_tangent[1][0];
	lambdamu_saved[res_ind][0][0].tangent[1][0] = ext_periodicfield[res_ind][0].full_tangent[0][1];
	lambdamu_saved[res_ind][0][0].tangent[1][1] = ext_periodicfield[res_ind][0].full_tangent[1][1];
	solve_lambdaRes(lambdamu_saved[res_ind], periodicfield_saved[res_ind], m0, L_fixedpoints[res_ind], N_gridphi_fieldperiod, Bfield_island[res_ind]);
	//clock_t int_afteralladjoints = clock();


	//derivative_matrix = set_identity();
	for (coil_ind=0; coil_ind < allparams.n_coils; coil_ind++) {
		gradXp[res_ind][coil_ind] = malloc(2*sizeof(double));
		gradXp[res_ind][coil_ind][0] = calloc(allparams.n_diff,sizeof(double));
		gradXp[res_ind][coil_ind][1] = calloc(allparams.n_diff,sizeof(double));
		shapeResheli[res_ind][coil_ind] = calloc(allparams.n_diff,sizeof(double));
	}
	lambdaXp[res_ind] = malloc(N_line*sizeof(struct position));
	lambdaXp_axis = malloc(N_gridphi_fieldperiod*sizeof(struct position));
	for (ind = 0; ind< N_line; ind++) {
		lambdaXp[res_ind][ind] = malloc(4*sizeof(struct position));
	}
	for (ind = 0; ind < N_gridphi_fieldperiod; ind++) lambdaXp_axis[ind] = malloc(4*sizeof(struct position));
	lambdaXp[res_ind][0][0].tangent = set_identity();
	lambdaXp_axis[0][0].tangent = set_identity();
	Xpguess[res_ind][0] = ext_periodicfield[res_ind][0].loc[0];
	Xpguess[res_ind][1] = ext_periodicfield[res_ind][0].loc[1];
    }
	//if (strncmp("heli", allparams.type, 4) == 0) {
	//	rowallowed[0][0] = 2; rowallowed[0][1] = 4; rowallowed[0][2] = 6;
	//	rowallowed[1][0] = 1; rowallowed[1][1] = 4; rowallowed[1][2] = 5;
    optparams = malloc(allparams.n_coils*sizeof(double));
    for (coil_ind=0; coil_ind<allparams.n_coils; coil_ind++) {
    	//rowallowed[coil_ind] = opt_indices[coil_ind]; 
        optparams[coil_ind] = malloc(allparams.n_diff*sizeof(double));
	for (param_ind=0; param_ind<allparams.n_diff; param_ind++) {
            optparams[coil_ind] = allparams.diffparams[coil_ind][param_ind];
	}
    }

    rowallowed[0][0] = 2; // need to make this an input of the optimization
    rowallowed[1][0] = 1;
    gradobjective[0] = 0.0;
    gradobjective[1] = 0.0;
    /* 
    optimization part of code below
    */
    coil_ind = 0;
    //tol_window[0] = 0.2*tolstep;
    tol_window[1] = 100.0;
    axisguess[0] = axis_saved[0][0].loc[0];
    axisguess[1] = axis_saved[0][0].loc[1];
    for (res_ind=0; res_ind<num_resonances; res_ind++) {
	if (fabs(periodicfield_saved[res_ind][0][0].loc[0] - axis_saved[0][0].loc[0]) < reference_change) 
	    reference_change = fabs(periodicfield_saved[res_ind][0][0].loc[0] - axis_saved[0][0].loc[0]);
    }
    tolstep = fracfpsep*reference_change;
    do {
	//if (count < count_int) allow = 1;
	if (count_since_optimum < count_int)  allow = 1;
	else if ( (count_since_optimum == count_int) && (count_int != 0) ) {
        //if (count == count_int) {
            //opfactor_reference /= 2.0;
	    printf("???coil parameters = (%14.14f, %14.14f)\n", allparams.diffparams[0][0][rowallowed[0][0]], allparams.diffparams[1][0][rowallowed[1][0]]);
            objective = objective_opt ;
            objectiveold = 9e13;
            allparams.diffparams[0][0][rowallowed[0][0]] = coilparams_opt[0] ;
            allparams.diffparams[1][0][rowallowed[1][0]] = coilparams_opt[1] ;
            rmsqavRes = rmsqavRes_opt ;
            axis_saved[0][0].loc[0] = axis_optimized[0] ;
            axis_saved[0][0].loc[1] = axis_optimized[1] ;
            for (res_ind=0; res_ind<num_resonances; res_ind++) {
                Res[res_ind] = Resopt[res_ind];	
                periodicfield_saved[res_ind][0][0].loc[0]  = periodicfield_optimized[res_ind][0];
                periodicfield_saved[res_ind][0][0].loc[1]  = periodicfield_optimized[res_ind][1];
            }
        }
	// else do nothing
	printf("ITERATION # = %d \n", count);
	stuck = 0;
	do {
	    reference_change = 10000000.0;
	    solve_magneticaxis_save(axisguess, axis_saved, allparams, Bfield_axis, N_gridphi_fieldperiod, error_max);
	    ind = 0;
	    saveR = axisguess[0];
	    while ( fabs(axis_saved[0][0].loc[0] - axisguess[0]) > tolstep ) {
	    	axisguess[0] = saveR + (ind/tol_window[1])*tolstep;
	    	solve_magneticaxis_save(axisguess, axis_saved, allparams, Bfield_axis, N_gridphi_fieldperiod, error_max);
	    	if ( fabs(axis_saved[0][0].loc[0] - axisguess[0]) > tolstep ) {
	    		axisguess[0] = saveR - (ind/tol_window[1])*tolstep;
	    		solve_magneticaxis_save(axisguess, axis_saved, allparams, Bfield_axis, N_gridphi_fieldperiod, error_max);
	    	}
	    	ind++;
	    	if (ind>tol_window[1]) {
	    		printf("ERROR: magnetic axis lost during iteration\n");
	    		exit(1);
	    	}
	    }
	    //printf("axis_guess = (%f, %f)\n", axisguess[0], axisguess[1]);
	    printf("axis = (%f, %f)\n", axis_saved[0][0].loc[0], axis_saved[0][0].loc[1]);
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

	    for (ind=0; ind<N_gridphi_fieldperiod; ind++) axis[ind] = axis_saved[ind][0];
	    objective = rmsqavRes = sumfactor = 0.0;
	    printf("coil parameters = (%14.14f, %14.14f)\n", allparams.diffparams[0][0][rowallowed[0][0]], allparams.diffparams[1][0][rowallowed[1][0]]);
	    for (res_ind=0; res_ind<num_resonances; res_ind++) {
	        q0_index[res_ind] = 0;
	        //printf("periodicfield_guess = (%10f, %10f)\n", Xpguess[res_ind][0], Xpguess[res_ind][1]);
	        extsolve_periodicfieldline(NULL, Xpguess[res_ind], ext_periodicfield[res_ind], periodicfield_saved[res_ind], Bfield_island[res_ind], axis, allparams, m0, L_fixedpoints[res_ind], pol_mode[res_ind], q0_index+res_ind, N_gridphi_fieldperiod, error_max);
	        //printf("periodicfield = (%10f, %10f)\n", periodicfield_saved[res_ind][0][0].loc[0], periodicfield_saved[res_ind][0][0].loc[1]);
	        ind = 0;
	        saveR = Xpguess[res_ind][0];
	        while ( fabs(periodicfield_saved[res_ind][0][0].loc[0] - Xpguess[res_ind][0]) > tolstep ) {
		    q0_index[res_ind] = 0;
		    Xpguess[res_ind][0] = saveR + (ind/tol_window[1])*tolstep;
		    extsolve_periodicfieldline(NULL, Xpguess[res_ind], ext_periodicfield[res_ind], periodicfield_saved[res_ind], Bfield_island[res_ind], axis, allparams, m0, L_fixedpoints[res_ind], pol_mode[res_ind], q0_index+res_ind, N_gridphi_fieldperiod, error_max);
		    if ( fabs(periodicfield_saved[res_ind][0][0].loc[0] - Xpguess[res_ind][0]) > tolstep ) {
		    	Xpguess[res_ind][0] = saveR - (ind/tol_window[1])*tolstep;
		    	extsolve_periodicfieldline(NULL, Xpguess[res_ind], ext_periodicfield[res_ind], periodicfield_saved[res_ind], Bfield_island[res_ind], axis, allparams, m0, L_fixedpoints[res_ind], pol_mode[res_ind], q0_index+res_ind, N_gridphi_fieldperiod, error_max);
		    }
		    ind++;
		    if (ind>tol_window[1]) {
		    	printf("ERROR: fixed point lost during iteration\n");
		    	exit(1);
		    }
	        }
	        printf("periodicfield = (%10f, %10f)\n", periodicfield_saved[res_ind][0][0].loc[0], periodicfield_saved[res_ind][0][0].loc[1]);
	        Res[res_ind] = ext_periodicfield[res_ind][0].Res;
	        if (fabs(periodicfield_saved[res_ind][0][0].loc[0] - axis_saved[0][0].loc[0]) < reference_change) 
	        reference_change = fabs(periodicfield_saved[res_ind][0][0].loc[0] - axis_saved[0][0].loc[0]);

	        tolstep = fracfpsep*reference_change;
	        tol_window[0] = tolstep;
	        lambdamu_saved[res_ind][0][0].tangent[0][0] = ext_periodicfield[res_ind][0].full_tangent[0][0];
	        lambdamu_saved[res_ind][0][0].tangent[0][1] = ext_periodicfield[res_ind][0].full_tangent[1][0];
	        lambdamu_saved[res_ind][0][0].tangent[1][0] = ext_periodicfield[res_ind][0].full_tangent[0][1];
	        lambdamu_saved[res_ind][0][0].tangent[1][1] = ext_periodicfield[res_ind][0].full_tangent[1][1];
	        solve_lambdaRes(lambdamu_saved[res_ind], periodicfield_saved[res_ind], m0, L_fixedpoints[res_ind], N_gridphi_fieldperiod, Bfield_island[res_ind]);
	        matrix[0][0] = 1.0 -  ext_periodicfield[res_ind][0].adj_full_tangent[0][0];
	        matrix[0][1] =     -  ext_periodicfield[res_ind][0].adj_full_tangent[0][1];
	        matrix[1][0] =     -  ext_periodicfield[res_ind][0].adj_full_tangent[1][0];
	        matrix[1][1] = 1.0 -  ext_periodicfield[res_ind][0].adj_full_tangent[1][1];
	        invertmatrix = invert2x2(matrix, &det);
	        lambdaXp[res_ind][0][0].tangent[0][0] = invertmatrix[0][0];
	        lambdaXp[res_ind][0][0].tangent[0][1] = invertmatrix[0][1];
	        lambdaXp[res_ind][0][0].tangent[1][0] = invertmatrix[1][0];
	        lambdaXp[res_ind][0][0].tangent[1][1] = invertmatrix[1][1];
	        solve_lambdaXp(lambdaXp[res_ind], periodicfield_saved[res_ind], m0, L_fixedpoints[res_ind], N_gridphi_fieldperiod, Bfield_island[res_ind]); 
	        free(invertmatrix[0]); free(invertmatrix[1]);
	        free(invertmatrix);
	        objective += factor[res_ind]*(Res[res_ind]);
	        rmsqavRes += Res[res_ind]*Res[res_ind];
		sumfactor += fabs(factor[res_ind]);
	    }
	    objective /= sumfactor;
	    rmsqavRes /= num_resonances;
	    rmsqavRes = sqrt(rmsqavRes);
	    //objective = rmsqavRes;

	    if ( ( ( fabs(objectiveold) < fabs(objective) ) && (allow == 0) ) && (stuck < stuck_max) ) { //- gradobjective[0]*deltaparam[0] - gradobjective[1]*deltaparam[1]
	        axisguess[0] = axis_saved[0][0].loc[0];
	        axisguess[1] = 0.0;
	        for (res_ind=0; res_ind<num_resonances; res_ind++) {
	            Xpguess[res_ind][0] = periodicfield_saved[res_ind][0][0].loc[0];			
	            Xpguess[res_ind][1] = 0.0;
	        }
	        for (coil_ind=0;coil_ind<2; coil_ind++) {
	            //printf("deltaparam[%d] = %f\n", coil_ind, deltaparam[coil_ind]);
		    //if (stuck>10) deltaparam[coil_ind] *= 4.0; 
	            allparams.diffparams[coil_ind][0][rowallowed[coil_ind][0]] += deltaparam[coil_ind];
	            axisguess[0] += gradXp_axis[coil_ind][0][rowallowed[coil_ind][0]]*deltaparam[coil_ind];
	            for (res_ind=0; res_ind<num_resonances; res_ind++) {
	                Xpguess[res_ind][0] += gradXp[res_ind][coil_ind][0][rowallowed[coil_ind][0]]*deltaparam[coil_ind];
	            }
	            deltaparam[coil_ind] /= r_decrease;
	            allparams.diffparams[coil_ind][0][rowallowed[coil_ind][0]] -= deltaparam[coil_ind];
	            axisguess[0] -= gradXp_axis[coil_ind][0][rowallowed[coil_ind][0]]*deltaparam[coil_ind];
	            for (res_ind=0; res_ind<num_resonances; res_ind++) {
	                Xpguess[res_ind][0] -= gradXp[res_ind][coil_ind][0][rowallowed[coil_ind][0]]*deltaparam[coil_ind];
	            }
	        }
		
		if (stuck == stuck_int) opfactor_reference /= r_decrease;
	        printf("stuck = %d\n", stuck);
	        stuck += 1;
	    }

	    printf("objective = %10.10f (previous iteration %10.10f)\n", objective, objectiveold);
	    printf("rmsqavRes = %10.10f\n", rmsqavRes);
	    opfactor = opfactor_reference;
	} while ( ( ( fabs(objectiveold) < fabs(objective) ) && (allow == 0) ) && (stuck < stuck_max) ); //gradobjective[0]*2.0*deltaparam[0] - gradobjective[1]*2.0*deltaparam[1]
	if ( ( fabs(objective)  < fabs(objective_opt) ) && (allow == 1) ) {
	    printf("STORE OPTIMUM\n");
	    objective_opt = objective;
	    coilparams_opt[0] = allparams.diffparams[0][0][rowallowed[0][0]];
	    coilparams_opt[1] = allparams.diffparams[1][0][rowallowed[1][0]];
	    for (coil_ind=0; coil_ind< allparams.n_coils; coil_ind++) 
		optparams[coil_ind][rowallowed[coil_ind][0]];
	    rmsqavRes_opt = rmsqavRes;
	    axis_optimized[0] = axis_saved[0][0].loc[0];
	    axis_optimized[1] = axis_saved[0][0].loc[1];
	    for (res_ind=0; res_ind<num_resonances; res_ind++) {
	    Resopt[res_ind] = Res[res_ind];	
	    periodicfield_optimized[res_ind][0] = periodicfield_saved[res_ind][0][0].loc[0];
	    periodicfield_optimized[res_ind][1] = periodicfield_saved[res_ind][0][0].loc[1];
	    }
	    //if (count_since_optimum < 1) opfactor_reference *= 1.1;
	    count_since_optimum = 0;
	}
	//if (stuck == stuck_max) printf("\n\n\n\n\n\n\n\nstuck = %d\n\n\n\n\n\n\n\n\n", stuck);
	//} while (find_local == 1);
	if ( (allow == 1) || (stuck < stuck_max )  ) {
	    printf("Take a new step in parameter space\n");
	    stuck = 0;
	    allow = 0;
	    gradobjectivesquared = 0.0;
	    for (coil_ind=0;coil_ind<2; coil_ind++) {
	        solve_gradXp(gradXp_axis[coil_ind], Bfield_axis, axis_saved, lambdaXp_axis, 1, N_gridphi_fieldperiod, allparams, coil_ind, 0) ;
	        gradXpmag[coil_ind] = sqrt( pow(gradXp_axis[coil_ind][0][rowallowed[coil_ind][0]], 2.0) + pow(gradXp_axis[coil_ind][1][rowallowed[coil_ind][0]], 2.0) );
	        axisguess[0] = axis_saved[0][0].loc[0];
	        axisguess[1] = 0.0;
	        gradobjective[coil_ind] = 0.0;
	        for (res_ind=0; res_ind<num_resonances; res_ind++) {	
	            solve_gradRes(shapeResheli[res_ind][coil_ind], Bfield_island[res_ind], periodicfield_saved[res_ind], lambdamu_saved[res_ind], L_fixedpoints[res_ind], N_gridphi_fieldperiod, allparams, coil_ind, 0) ;
	            solve_gradXp(gradXp[res_ind][coil_ind], Bfield_island[res_ind], periodicfield_saved[res_ind], lambdaXp[res_ind], L_fixedpoints[res_ind], N_gridphi_fieldperiod, allparams, coil_ind, 0) ;
	            //gradobjective[coil_ind] += 0.25*2.0*Res[res_ind]*shapeResheli[res_ind][coil_ind][rowallowed[coil_ind][0]];
	            gradobjective[coil_ind] += 0.25*factor[res_ind]*shapeResheli[res_ind][coil_ind][rowallowed[coil_ind][0]];
	            //gradobjective[coil_ind] += 0.25*(Res[res_ind]/fabs(Res[res_ind]))*shapeResheli[res_ind][coil_ind][rowallowed[coil_ind][0]];
	            if ( gradXpmag[coil_ind] < sqrt( pow(gradXp[res_ind][coil_ind][0][rowallowed[coil_ind][0]], 2.0) + pow(gradXp[res_ind][coil_ind][1][rowallowed[coil_ind][0]], 2.0) ) )  
	                gradXpmag[coil_ind] = sqrt( pow(gradXp[res_ind][coil_ind][0][rowallowed[coil_ind][0]], 2.0) + pow(gradXp[res_ind][coil_ind][1][rowallowed[coil_ind][0]], 2.0) );
	        }
	        gradobjective[coil_ind] /= sumfactor;
	        //gradobjective[coil_ind] /= (num_resonances*objective);
		gradobjectivesquared += (gradobjective[coil_ind]*gradobjective[coil_ind]);
	    }
	    numb = 0.0;
	    for (coil_ind=0;coil_ind<2; coil_ind++) {
	        //limit = (tolstep/gradXpmag)*fabs(gradobjective[coil_ind]/objective);
	        numb += fabs( (objective/gradobjectivesquared)*gradobjective[coil_ind]*gradXpmag[coil_ind] );
	        //limit = (tolstep/gradXpmag)*fabs(gradobjectivesquared/(objective*gradobjective[coil_ind]));
	        printf("gradXpmag = %f, gradobjective[%d] = %f, objective = %f\n", gradXpmag[coil_ind], coil_ind, gradobjective[coil_ind], objective);
	    }
	    limit = fabs(tolstep/numb);
	    if (opfactor > limit) opfactor = limit;
	    printf("opfactor= %f, limit = %f, tolstep = %f\n", opfactor, limit, tolstep);
	    for (res_ind=0; res_ind < num_resonances; res_ind++) {
	    	Xpguess[res_ind][0] = periodicfield_saved[res_ind][0][0].loc[0];			
	    	Xpguess[res_ind][1] = 0.0;
	    }
	    for (coil_ind=0;coil_ind<2; coil_ind++) {
	        //deltaparam[coil_ind] = opfactor*(objective/gradobjective[coil_ind]); 
	        deltaparam[coil_ind] = opfactor*(objective*gradobjective[coil_ind]/gradobjectivesquared); 
	        //printf("deltaparam[%d] = %f\n", coil_ind, deltaparam[coil_ind]);
	        allparams.diffparams[coil_ind][0][rowallowed[coil_ind][0]] -= deltaparam[coil_ind];
	        axisguess[0] -= gradXp_axis[coil_ind][0][rowallowed[coil_ind][0]]*deltaparam[coil_ind];
	        for (res_ind=0; res_ind < num_resonances; res_ind++) 
	            Xpguess[res_ind][0] -= gradXp[res_ind][coil_ind][0][rowallowed[coil_ind][0]]*deltaparam[coil_ind];
	    }
	    //if (count%2 == 0) coil_ind = 0;
	    //else coil_ind = 1;
	    //    deltaparam[coil_ind] = opfactor*(objective/gradobjective[coil_ind]);
	    //    //printf("deltaparam[%d] = %f\n", coil_ind, deltaparam[coil_ind]);
	    //    allparams.diffparams[coil_ind][0][rowallowed[coil_ind][0]] -= deltaparam[coil_ind];
	    //    axisguess[0] -= gradXp_axis[coil_ind][0][rowallowed[coil_ind][0]]*deltaparam[coil_ind];
	    //    for (res_ind=0; res_ind < num_resonances; res_ind++) 
	    //        Xpguess[res_ind][0] -= gradXp[res_ind][coil_ind][0][rowallowed[coil_ind][0]]*deltaparam[coil_ind];
	}
	else {
	    stuck = 0;
	    allow = 1;
	    number_minima += 1;
	    if (number_minima == max_minima) count = count_max;
	    //printf("# MINIMA = %d\n", number_minima);
	    //printf("optimized coil parameters = (%f, %f)\n", coilparams_opt[0], coilparams_opt[1]);
	    //printf("optimized Residues = (%f, %f)\n", Resopt[0], Resopt[1]);
	    //printf("optimized average of magnitude of residues = %f\n", objective_opt);
	    //printf("optimized rmsq average of residues = %f\n", rmsqavRes_opt);
	}
	if (fabs( fabs(objectiveold) - fabs(objective) ) < stop_threshold ) {
	    count = count_max;
	}
	objectiveold = objective;		
	count++;
	count_since_optimum++;
	printf("count = %d\n", count);
    } while (fabs(objective) > goal && count < count_max); 
    printf("optimized coil parameters = (%f, %f)\n", coilparams_opt[0], coilparams_opt[1]);
    printf("optimized Residues = (%f, %f)\n", Resopt[0], Resopt[1]);
    printf("optimized rmsq average of residues = %f\n", rmsqavRes_opt);
    printf("optimized objective function = %f\n", objective_opt);

    FILE *fixedpoints_optimized, *magfieldparams_optimized;
    if ( (fixedpoints_optimized = fopen("fixedpoints_optimized.txt", "w")) == NULL) {
	printf("ERROR: file fixedpoints_optimized.txt didn't open\n");
	exit(1);
    }
    fprintf(fixedpoints_optimized, "%f %f\n", axis_optimized[0], axis_optimized[1]);
    for (res_ind=0; res_ind < num_resonances; res_ind++) {
        fprintf(fixedpoints_optimized, "%f %f %d %d\n", periodicfield_optimized[res_ind][0], periodicfield_optimized[res_ind][1], pol_mode[res_ind], tor_mode[res_ind]);
    }
    fclose(fixedpoints_optimized);
    //if ( (magfieldparams_optimized = fopen("magfieldparams.txt", "r+")) == NULL) {
    //    printf("ERROR: file magfieldparams.txt didn't open\n");
    //    exit(1);
    //}
    //while ( (fgets(line, lenline, fixed_points) != NULL) || (isdigit(line[0]) == 0) ) printf("%s\n", line);	
    //for (coil_ind=0; coil_ind < allparams.n_coils; coil_ind++) {
    //    for (param_ind=0; param_ind < allparams.n_diff; param_ind++) 
    //        fprintf(magfieldparams_optimized, "%f ", optparams[coil_ind][param_ind]);
    //    fprintf(magfieldparams_optimized, "\n");
    //}
    //fclose(magfieldparams_optimized);

    clock_t intopt = clock();
    printf("time after solving for everything: %f\n", (double) (intopt-start)/CLOCKS_PER_SEC);
    return 0;
}
