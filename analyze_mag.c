//Currently under construction
//Author: Alessandro Geraldini
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "isc.h"
#include <gsl/gsl_sf_bessel.h>
#define N_gridphi_fieldperiod 30
#define error_max 1e-13
#define poiniota 'n'
#define calcshap 'y'
#define convergence_test 0
#define testwidth_Reim 0

int main()
{
// store starting time to keep track of time
    int study_Reim = 0;
    clock_t start = clock();
// provide values of options to the code
    char checshap = 'y', checfast = 'y', checkmag = 'y';
// variables below are just there to test the magnetic field evaluation
    double test_RZ[2] = {0.95, 0.0}, test_varphi = 0.0;
    struct field Bpoint;
    int *search, L_fixedpoints;
    int *pol_mode, *tor_mode;
    double *fixedpoint_info, axis_RZ[2], **fixedpoint_RZ;
    struct position *axis, **mu, **lambdaQ, lambda, **lambdamu_saved;
    double RZ[2];
    int m0, N_gridphi_tor;
    int num_coils, num_resonances;
    int ind, row, col, coil_ind, seg_ind, fixed_ind, res_ind, param_ind;
    int delta_ind_max = 20, delta_ind;
    double powdeltain = -0.5, steppowdelta = -3.5;
    int N_line;
    int seg_step = 1;
    int max_coil_ind, num_segs, max_num_segs = 80, seg_start = 0;
    int line_elements, *q0_index; 
    struct position ***periodicfield_saved, **axis_saved;
    struct ext_position **ext_periodicfield, *ext_periodicfield_FD[2];
    struct field ***Bfield_island, **Bfield_axis;
// structure below contains all the magnetic field parameters and its type --> see function fetchparams() in magfield.c
    struct fieldparams allparams;
// pointer coils is just equivalent to allparams.diffparams but a convenient shorthand for when allparams.type == "coil"
    double ***coils;
/* 
the grad pointers store the gradient of the Residue, circumference, Sigma and width with respect to all the parameters of the magnetic field
 first pointer dimension corresponds to number of resonances studied (# lines - 1 in the file fixedpoints.txt)
 second pointer dimension corresponds to number of coils in magnetic field parameter file (one for allparams.type = "Reim")
 third pointer dimension corresponds to number of coil segments (this is always one except for when allparams.type = "coil")
 fourth pointer dimension corresponds to number of variables per coil segment (equal to allparams.n_diff)
 ONLY for gradSigma and gradwidth, fifth pointer dimension corresponds to number of fixed points in island chain (L_fixedpoints)
this order may not be the most intuitive (especially the first and last pointer dimensions would have been better next to each other)
*/
    double ****gradRes, ****gradcirc, *****gradSigma, *****gradwidth;
// below, additional pointer dimension corresponds to number of finite difference step sizes used to check the gradient values
    double *****checknormgradRes, *****checknormgradcirc, ******checknormgradSigma, ******checknormgradwidth;
/* 
Pointers below only differ from above ones if allparams.type == "coil". In this case realgrad evaluates the gradient with respect to variation of one coil segment from the shape gradient. A posteriori, I realize this might be confusing: the real shape gradient (in cross product form, from the paper) is grad realgrad is what you get when you take a set of coil positions, infinitesimally displace one of them, and divide by the displacement
*/
    double *****realgradSigma, *****realgradwidth, ****realgradcirc, ****realgradRes;  
/* 
 pointers below are like the ones above but the gradient of the quantity is normalised by the value of the quantity i.e. normgradquantity = gradlogquantity this pointer has lower dimensionality because it is not stored, but it is always reevaluated replacing the previous value
*/
    double normgradwidth, normgradcirc, normgradSigma, normgradRes;
    double phichec = 0.0, factor, dldown, dl, savedcoilseg;
    struct field realBpointgrad[3];
    double dr[3], drdown[3], *checkparam = malloc(4*sizeof(double));
    double *iota=NULL, *rmin=NULL, *zmin=NULL, axis_phi0[2];
    double **evec=malloc(2*sizeof(double*)), *eval=malloc(2*sizeof(double)), trace, det, iota_axis;
    double delta = 0.000001, checkBpointgrad[3], checknablaBpointgrad[3][2];
    struct field *Bpointgrad, Bpointdelta[3][2];
// variables below are just for reading the file fixedpoints.txt
    int lenline = 1000, num_line_file=0; 
    char line[lenline];
/* 
 the file fixedpoints.txt contains the (R, Z) coordinate of the magnetic axis on the first line on the following lines, in no particular order, the approximate fixed point (R, Z) coordinates are listed the last two integers on each line (but the first) are the poloidal and toroidal mode number of the fixed point
*/
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
    fixedpoint_RZ = malloc(num_resonances*sizeof(double));
    ext_periodicfield = malloc(num_resonances*sizeof(struct ext_position));
    Bfield_island = malloc(num_resonances*sizeof(struct field));
    periodicfield_saved = malloc(num_resonances*sizeof(struct position));
    q0_index = calloc(num_resonances, sizeof(int));
    tor_mode = malloc(num_resonances*sizeof(int));
    pol_mode = malloc(num_resonances*sizeof(int));
    search = malloc(num_resonances*sizeof(int));
    num_line_file = 0;
    line_elements = 4;
    while (fgets(line, lenline, fixed_points) != NULL) {	
	line_elements = 4; // has to be larger than what the number will actually turn out to be
	fixedpoint_info = linetodata(line, &line_elements);	
	if (num_line_file == 0) {
	    //printf("line_elements = %d\n", line_elements);
	    axis_RZ[0] = fixedpoint_info[0];
	    axis_RZ[1] = fixedpoint_info[1];
	    printf("(R, Z) = (%f, %f) for magnetic axis\n", axis_RZ[0], axis_RZ[1]);
	}
	else {
	    //printf("line_elements = %d\n", line_elements);
	    fixedpoint_RZ[num_line_file-1] = malloc(2*sizeof(double));
	    fixedpoint_RZ[num_line_file-1][0] = fixedpoint_info[0];
	    fixedpoint_RZ[num_line_file-1][1] = fixedpoint_info[1];
	    if (line_elements == 4) {
	        pol_mode[num_line_file-1] = (int) fixedpoint_info[2];
	        tor_mode[num_line_file-1] = (int) fixedpoint_info[3];
	    }
	    else {
	        pol_mode[num_line_file-1] = 0;
	        tor_mode[num_line_file-1] = 0;
	    }
	    printf("(R, Z) = (%f, %f) for pol_mode = %d, tor_mode = %d\n", fixedpoint_RZ[num_line_file-1][0], fixedpoint_RZ[num_line_file-1][1], pol_mode[num_line_file-1], tor_mode[num_line_file-1]);
	    printf("N.B. if pol_mode and tor_mode are unknown they should be zero\n");
	}
	num_line_file += 1;
	free(fixedpoint_info);
    }
    fclose(fixed_points);
// insert explanation of q0_index here
/* 
Specify the options of the code:
Poincare = ('y', 'n') --> want to produce a Poincare plot & iota profile?
checshap = ('y', 'n') --> want to check shape gradient using finite differences?
*/
    printf("checkmag = %c\npoiniota = %c\n\nchecshap = %c\nchecfast = %c\n", checkmag, poiniota, checshap, checfast);
/* 
 the function fetchparams is ESSENTIAL (see more detail on the function in the file magfield.c)
 it returns the structure allparams which contains all the magnetic field parameters 
 in principle it can describe any magnetic field configuration type 
 of course different types require existing functions to evaluate the magnetic field using the parameters in allparams
 currently supported types are: coil, Reim, heli
 coil provides a list of position coordinates along a set of coils
 Reim provides the list of parameters of the analytical Reiman model:
 (iota0, iota1, epsilon_i) for i =(0, ..., N)
 heli provides a list of parameters describing a pair of helical coils
 (if time allows could implement Dommaschk field as I have it on a separate old script)
*/
    allparams = fetchparams(); // assign magnetic field parameters to the structure allparams
    m0 = allparams.m0_fieldperiods; // assign an element of the structure to the integer m0 (for convenience)
    N_gridphi_tor = m0*N_gridphi_fieldperiod; // number of grid points in a toroidal field period 2*pi/m0
    //N_gridphi_fieldperiod = N_gridphi_tor/m0; // number of grid points in a toroidal field period 2*pi/m0
    num_coils = allparams.n_coils; // assign element of structure to integer num_coils (for convenience) 
    printf("The magnetic field type is %s\n", allparams.type);
    evec[0]=malloc(2*sizeof(double)); evec[1]=malloc(2*sizeof(double));
    printf("field periods = %d\nN points per field period = %d", m0, N_gridphi_fieldperiod);
    clock_t tcoils = clock(); printf("... %f seconds ...\n", (double) (tcoils-start)/CLOCKS_PER_SEC);
    // a test of the magnetic field evaluation at a point specified by test_RZ = (R, Z) and toroidal angle test_varphi
    Bpoint = Bfield(test_RZ, test_varphi, allparams);
    printf("At test position (R, Z, φ) = (%f, %f, %f) magnetic field is (BR, BZ, Bphi) = (%f, %f, %f)\n", test_RZ[0], test_RZ[1], test_varphi, Bpoint.value[0], Bpoint.value[1], Bpoint.value[2]);
/* 
 this just introduces a shorthand notation for coil positions
 useful for coil type as otherwise allparams.diffparams gets written often
*/
    coils = allparams.diffparams;
// the code in the loop below might need to be modified depending on which field one wants to test
    if (checkmag == 'y') { 
	if (strncmp("coil", allparams.type, 5) == 0 ) {
	    RZ[0] = 1.53; RZ[1]= 0.0; //NCSX
	    coil_ind = 5;
	    seg_ind = 14;
	    checkparam[0] = coils[coil_ind][seg_ind][0];
	    checkparam[1] = coils[coil_ind][seg_ind][1];
	    checkparam[2] = coils[coil_ind][seg_ind][2];
	    checkparam[3] = coils[coil_ind][seg_ind][3];
	    dr[0] = coils[coil_ind][seg_ind+1][0] - coils[coil_ind][seg_ind][0];
	    dr[1] = coils[coil_ind][seg_ind+1][1] - coils[coil_ind][seg_ind][1];
	    dr[2] = coils[coil_ind][seg_ind+1][2] - coils[coil_ind][seg_ind][2];
	    Bpointgrad = gradBfield(RZ, phichec, allparams, coil_ind, seg_ind);
	    for (row=0;row<3;row++) {
	    	for (ind=0;ind<3; ind++) realBpointgrad[ind].value[row] = (dr[(ind+1)%3]*Bpointgrad[(ind+2)%3].value[row] - dr[(ind+2)%3]*Bpointgrad[(ind+1)%3].value[row]); 
	    }
	    coils[coil_ind][seg_ind][0] += delta; 
	    Bpointdelta[0][1] = Bfield(RZ, phichec, allparams); 
	    coils[coil_ind][seg_ind][0] -= 2.0*delta; 
	    Bpointdelta[0][0] = Bfield(RZ, phichec, allparams); 
	    coils[coil_ind][seg_ind][0] += delta; 
	    coils[coil_ind][seg_ind][1] += delta; 
	    Bpointdelta[1][1] = Bfield(RZ, phichec, allparams); 
	    coils[coil_ind][seg_ind][1] -= 2.0*delta; 
	    Bpointdelta[1][0] = Bfield(RZ, phichec, allparams); 
	    coils[coil_ind][seg_ind][1] += delta; 
	    coils[coil_ind][seg_ind][2] += delta; 
	    Bpointdelta[2][1] = Bfield(RZ, phichec, allparams); 
	    coils[coil_ind][seg_ind][2] -= 2.0*delta; 
	    Bpointdelta[2][0] = Bfield(RZ, phichec, allparams); 
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
	    RZ[0] = 1.05; RZ[1]= 0.0; //NCSX
	    printf("(R, Z) = (%f, %f)\n", RZ[0], RZ[1]);
	    Bpoint = Bfield(RZ, phichec, allparams); 
	    printstructfield("Bpoint =", &Bpoint);
	}
	else if (strncmp("heli", allparams.type, 5) == 0 ) {
	    RZ[0] = 0.95; RZ[1]= 0.0; 
	    coil_ind = 1;
	    delta = 0.0001;
	    printf("(R, Z) = (%f, %f)\n", RZ[0], RZ[1]);
	    printf("num_coils = %d\n", num_coils);
	    Bpointgrad = gradBfield(RZ, phichec, allparams, coil_ind, 0);
	    Bpoint = Bfield(RZ, phichec, allparams); 
	    for (row = 0; row<allparams.n_diff; row++) {
		allparams.diffparams[coil_ind][0][row] += delta;
		Bpointdelta[1][1] = Bfield(RZ, phichec, allparams); 
		allparams.diffparams[coil_ind][0][row] -= 2.0*delta;
		Bpointdelta[1][0] = Bfield(RZ, phichec, allparams); 
		allparams.diffparams[coil_ind][0][row] += delta;
		for (col=0;col<3;col++) {
		    checkBpointgrad[col] = ( Bpointdelta[1][1].value[col] - Bpointdelta[1][0].value[col] ) / (2.0*delta);
		    for (ind=0; ind<2; ind++) {
		    	checknablaBpointgrad[col][ind] = ( Bpointdelta[1][1].derivative[col][ind] - Bpointdelta[1][0].derivative[col][ind] ) / (2.0*delta);
		    }
		    printf("checkgradnablaB = (%14.14f, %14.14f) (direct), (%14.14f, %14.14f) (FD) \n", Bpointgrad[row].derivative[col][0], Bpointgrad[row].derivative[col][1], checknablaBpointgrad[col][0], checknablaBpointgrad[col][1]);
		}
		printf("checkgradB = (%14.14f, %14.14f, %14.14f) (direct), (%14.14f, %14.14f, %14.14f) (FD) \n", Bpointgrad[row].value[0], Bpointgrad[row].value[1], Bpointgrad[row].value[2], checkBpointgrad[0], checkBpointgrad[1], checkBpointgrad[2]);
	    }
	}
    }
/* 
below the magnetic axis solution is obtained 
the first guess of the axis position in the plane varphi = 0 is in the file fixedpoints.txt
the first line in the file contains the R, Z coordinates of a guess for the magnetic axis
the guess has to be good enough for the Newton iteration in the function solve_magneticaxis_save to converge
*/
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
    axis_phi0[0] = axis[0].loc[0];  axis_phi0[1] = axis[0].loc[1];
    printstructposition("axis", axis_saved[0]);
    linalg2x2(axis_saved[0][0].tangent, evec, eval, &det, &trace);
    //printf("evalss=(%f, %f)\n", eval[0], eval[1]);
    printstructposition("axis", axis_saved[0]);
    iota_axis = -eval[1]*m0/(2*M_PI);
    printf("iota_axis=%14.14f\n", iota_axis);
    //gradRes = calloc(3,sizeof(double));
    //solve_gradRes(gradRes, Bfield_axis, axis_saved, lambdamu_saved, 1, N_gridphi_fieldperiod, allparams, 0, 0) ;
    ////grad_iota[i] = 0.5*m0*gradRes[i]/(2.0*M_PI*sin(-eval[1]))
    //printf("iota_grad=%14.14f\n", 0.5*m0*gradRes[0]/(2.0*M_PI*sin(-eval[1])));
    //printf("iota_grad=%14.14f\n", 0.5*m0*gradRes[1]/(2.0*M_PI*sin(-eval[1])));
    //printf("iota_grad=%14.14f\n", 0.5*m0*gradRes[2]/(2.0*M_PI*sin(-eval[1])));
    //printf("This number should be one: %f\n", eval[0]);
/* 
 the following is an option to make a Poincare plot and obtain the iota profile
 perhaps it is useful to have a separate script for this, as one normally chooses to make a Poincare plot 
 in order to examine it and gain more information 
 alternatively one can only obtain an iota profile (which is faster) 
 the iota profile can give a sense of where to expect resonances
*/
    if (poiniota == 'y') {
    	//r_interval = 0.001; //n_points = 10;
    	printf("About to enter Poincare module\n");
    	printstructposition("axis", axis);
    	iotaprofile(*axis, m0, N_gridphi_fieldperiod, rmin, zmin, iota, allparams);
    	clock_t tPoincare = clock();
    	printf("Time after filling Poincare plot file: %f\n", (double) (tPoincare-start)/CLOCKS_PER_SEC);
    	//exit(0);
    }

    if (calcshap == 'y') {
	gradSigma = malloc(num_resonances*sizeof(double));
	gradwidth = malloc(num_resonances*sizeof(double));
	gradcirc = malloc(num_resonances*sizeof(double));
	gradRes = malloc(num_resonances*sizeof(double));
	realgradSigma = malloc(num_resonances*sizeof(double));
	realgradwidth = malloc(num_resonances*sizeof(double));
	realgradcirc = malloc(num_resonances*sizeof(double));
	realgradRes = malloc(num_resonances*sizeof(double));
	if (checshap == 'y') {
	    checknormgradSigma = malloc(num_resonances*sizeof(double));
	    checknormgradwidth = malloc(num_resonances*sizeof(double));
	    checknormgradcirc = malloc(num_resonances*sizeof(double));
	    checknormgradRes = malloc(num_resonances*sizeof(double));
	}
    }
//Xp->loc[0] = 1.3; Xp->loc[1]= -0.515; // island(?) is at(1.299008, -0.515194)
//Xp->loc[0] = 0.945; Xp->loc[1]= 0.0;//Dommaschk (5,2) amp 1.73: island is at(1.299008, -0.515194)
//Xp->loc[0] = 1.091; Xp->loc[1]= 0.0; // Dommaschk (5,2) amp 0.00001 on top of loop current
//Xp->loc[0] = 1.0097; Xp->loc[1]= 0.0; // Dommaschk (5,2) amp 0.00001 on top of loop current

/* 
Below the periodic field line solutions are obtained.
The first guess of each periodic field line position is in fixedpoints.txt
Each line (except the first one) in the file contains the following: 
R, Z, pol_mode, tor_mode
*/
	printf("YOLO\n\n\n\n\n\n\n\n");
    FILE *island_Reiman;
    if ( (island_Reiman = fopen("island_Reiman.txt", "w")) == NULL ) {
	printf("Error opening file island_Reiman.txt");
    }
    if (testwidth_Reim == 1) study_Reim = 1; 
    do {
    for (res_ind=0; res_ind< num_resonances; res_ind++) {
	printf("Resonance number %d\n", res_ind+1);
	q0_index[res_ind] = 0;
    	if (tor_mode[res_ind] % m0 == 0) L_fixedpoints = pol_mode[res_ind];	
    	else L_fixedpoints = m0*pol_mode[res_ind];
    	N_line = L_fixedpoints*N_gridphi_fieldperiod;
    	Bfield_island[res_ind] = malloc(N_line*sizeof(struct field));
    	periodicfield_saved[res_ind] = malloc(N_line*sizeof(struct position));
    	for (ind=0; ind<N_line; ind++) {
	    Bfield_island[res_ind][ind] = malloc(4*sizeof(struct field));
	    periodicfield_saved[res_ind][ind] = malloc(4*sizeof(struct position));
    	}
    	ext_periodicfield[res_ind] = malloc(L_fixedpoints*sizeof(struct ext_position));
    	ext_periodicfield_FD[0] = malloc(L_fixedpoints*sizeof(struct ext_position));
    	ext_periodicfield_FD[1] = malloc(L_fixedpoints*sizeof(struct ext_position));
    	printf("Entering function that solves for the periodic field line in extended way\n");
    	search[res_ind]  = extsolve_periodicfieldline(fixedpoint_RZ[res_ind], ext_periodicfield[res_ind], periodicfield_saved[res_ind], Bfield_island[res_ind], axis, allparams, m0, L_fixedpoints, pol_mode+res_ind, q0_index+res_ind, N_gridphi_fieldperiod, error_max);
	printf("eperp = (%f, %f)\nepar = (%f, %f)\n", ext_periodicfield[0][0].eperp[0], ext_periodicfield[0][0].eperp[1], ext_periodicfield[0][0].epar[0], ext_periodicfield[0][0].epar[1]);
	printf("q0 = %d\n", q0_index[0]);
    	if (search[res_ind] == 0) {
	    //periodicfield_saved[0][0].tangent[0][0] = 1.0;
	    //periodicfield_saved[0][0].tangent[0][1] = 0.0;
	    //periodicfield_saved[0][0].tangent[1][0] = 0.0;
	    //periodicfield_saved[0][0].tangent[1][1] = 1.0;
	    printf("List of fixed points and their basic info: \n");
	    for (fixed_ind = 0; fixed_ind<L_fixedpoints; fixed_ind++) {
	    	printf("island_center[%d].(R,Z) = (%f, %f)\n", fixed_ind, ext_periodicfield[res_ind][fixed_ind].loc[0], ext_periodicfield[res_ind][fixed_ind].loc[1]);
	    	printf("island_center[%d].(Res, width, Sigma, circ) = (%f, %f, %f, %f)\n", fixed_ind, ext_periodicfield[res_ind][fixed_ind].Res, ext_periodicfield[res_ind][fixed_ind].width, ext_periodicfield[res_ind][fixed_ind].Sigma, ext_periodicfield[res_ind][fixed_ind].circ);
	    }
	    printf("List of fixed points and their basic info: \n");
	    clock_t int_afterislandcentre = clock();
	    printf("... %f seconds ...\n", (double) (int_afterislandcentre-start)/CLOCKS_PER_SEC);
	    printf("calcshap = %c\n", calcshap);
	    if (calcshap == 'y') {
 	        //printf("Entering function that calculates the width of all islands in this chain\n");
 	        //printf("Left function that calculates the width of all islands in this chain\n");
 	        printf("Entering function that solves for the adjoint variable λ (circumference)\n");
 	        lambda = solve_lambda_circ(ext_periodicfield[res_ind], m0, L_fixedpoints, N_gridphi_fieldperiod);
 	        clock_t int_afterlambda = clock();
 	        printf("... %f seconds ...\n", (double) (int_afterlambda-start)/CLOCKS_PER_SEC);
 	        //printf("time after following adjoint variable for circumference: %f\n", (double) (int3-start)/CLOCKS_PER_SEC);
 	        //printf("time after calculating shape gradient of circumference: %f\n", (double) (int4-start)/CLOCKS_PER_SEC);
 	        printf("Entering function that solves for the adjoint variable μ\n");
 	        mu = solve_mu_tangent(ext_periodicfield[res_ind], m0, L_fixedpoints, q0_index[res_ind], N_gridphi_fieldperiod);
 	        printf("Entering function that solves for the adjoint variable λ_Q -->\n");
 	        lambdaQ = solve_lambda_tangent(Bfield_island[res_ind], ext_periodicfield[res_ind], mu, m0, L_fixedpoints, q0_index[res_ind], N_gridphi_fieldperiod);
 	        //printmat("adj_full_tangent", ext_periodicfield[0].adj_full_tangent, 2, 2);
 	        //printmat("full_tangent", ext_periodicfield[0].full_tangent, 2, 2);
 	        //printf("N_gridphi = %d, N_line = %d\n", N_gridphi_fieldperiod, N_line);
 	        lambdamu_saved = malloc(N_line*sizeof(struct position));
 	        for (ind=0; ind<N_line; ind++) lambdamu_saved[ind] = malloc(4*sizeof(struct position));
 	        lambdamu_saved[0][0].tangent = malloc(2*sizeof(double));
 	        lambdamu_saved[0][0].tangent[0] = malloc(2*sizeof(double));
 	        lambdamu_saved[0][0].tangent[1] = malloc(2*sizeof(double));
 	        lambdamu_saved[0][0].tangent[0][0] = ext_periodicfield[res_ind][0].full_tangent[0][0];
 	        lambdamu_saved[0][0].tangent[0][1] = ext_periodicfield[res_ind][0].full_tangent[1][0];
 	        lambdamu_saved[0][0].tangent[1][0] = ext_periodicfield[res_ind][0].full_tangent[0][1];
 	        lambdamu_saved[0][0].tangent[1][1] = ext_periodicfield[res_ind][0].full_tangent[1][1];
 	        solve_lambdaRes(lambdamu_saved, periodicfield_saved[res_ind], m0, L_fixedpoints, N_gridphi_fieldperiod, Bfield_island[res_ind]);
 	        //printmat("full_tangent", ext_periodicfield[res_ind][0].full_tangent, 2, 2);
 	        clock_t int_afteralladjoints = clock();
 	        printf("... %f seconds ...\n", (double) (int_afteralladjoints-start)/CLOCKS_PER_SEC);
 
 	        if ( (checfast == 'y') && (strncmp(allparams.type, "coil", 4) == 0) ) max_coil_ind = 1;
 	        else max_coil_ind = num_coils;
 	        
 	        //FILE *file_gradcirc, *file_gradSigma, *file_gradwidth, *file_gradRes;
		FILE *file_gradwidth, *file_gradRes;
 	        //if ( (file_gradcirc = fopen("grad_circ.txt", "w") ) == NULL ) {	
 	        //	printf("ERROR: couldn't open grad_circ.txt\n");
 	        //	exit(EXIT_FAILURE);
 	        //}	
 	        //if ( (file_gradSigma = fopen("grad_tan.txt", "w") ) == NULL ) {	
 	        //	printf("ERROR: couldn't open grad_tan.txt\n");
 	        //	exit(EXIT_FAILURE);
 	        //}	
 	        if ( (file_gradwidth = fopen("grad_width.txt", "w") ) == NULL ) {	
 	        	printf("ERROR: couldn't open grad_width.txt\n");
 	        	exit(EXIT_FAILURE);
 	        }	
 	        if ( (file_gradRes = fopen("grad_Res.txt", "w") ) == NULL ) {	
 	        	printf("ERROR: couldn't open grad_Res.txt\n");
 	        	exit(EXIT_FAILURE);
 	        }	
		max_coil_ind = 1; // TEMPORARY delete this line 
 	        gradcirc[res_ind] = malloc(max_coil_ind*sizeof(double));
 	        gradRes[res_ind] = malloc(max_coil_ind*sizeof(double));
 	        gradSigma[res_ind] = malloc(max_coil_ind*sizeof(double));
 	        gradwidth[res_ind] = malloc(max_coil_ind*sizeof(double));
 	        realgradcirc[res_ind] =malloc(max_coil_ind*sizeof(double));
 	        realgradRes[res_ind] = malloc(max_coil_ind*sizeof(double));
 	        realgradSigma[res_ind] = malloc(max_coil_ind*sizeof(double));
 	        realgradwidth[res_ind]=malloc(max_coil_ind*sizeof(double));
 	        for (coil_ind=0; coil_ind<max_coil_ind; coil_ind++) {
 	            printf("coil_ind=%d/%d\n", coil_ind, num_coils);
 	            if (strncmp(allparams.type, "coils", 4) == 0) num_segs = allparams.intparams[coil_ind];
 	            else num_segs = 1; // for anything but the "coil" type, the 2nd index of diff_params is always 0
 	            gradcirc[res_ind][coil_ind] = malloc(num_segs*sizeof(double));
 	            gradRes[res_ind][coil_ind] = malloc(num_segs*sizeof(double));
 	            gradSigma[res_ind][coil_ind] = malloc(num_segs*sizeof(double));
 	            gradwidth[res_ind][coil_ind] = malloc(num_segs*sizeof(double));
 	            realgradcirc[res_ind][coil_ind] =malloc(num_segs*sizeof(double));
 	            realgradRes[res_ind][coil_ind] = malloc(num_segs*sizeof(double));
 	            realgradSigma[res_ind][coil_ind] = malloc(num_segs*sizeof(double));
 	            realgradwidth[res_ind][coil_ind] =malloc(num_segs*sizeof(double));
// This is temporary, remove!
//num_segs = max_num_segs + seg_start;
 	            for ( seg_ind=seg_start; seg_ind < num_segs  ; seg_ind+=seg_step) {
 	            	gradcirc[res_ind][coil_ind][seg_ind] = calloc(allparams.n_diff,sizeof(double));
 	            	gradRes[res_ind][coil_ind][seg_ind] = calloc(allparams.n_diff,sizeof(double));
 	            	gradSigma[res_ind][coil_ind][seg_ind] = malloc(L_fixedpoints*sizeof(double));
 	            	gradwidth[res_ind][coil_ind][seg_ind] = malloc(L_fixedpoints*sizeof(double));
 	            	realgradcirc[res_ind][coil_ind][seg_ind] =calloc(allparams.n_diff,sizeof(double));
 	            	realgradRes[res_ind][coil_ind][seg_ind] = calloc(allparams.n_diff,sizeof(double));
 	            	realgradSigma[res_ind][coil_ind][seg_ind] = malloc(L_fixedpoints*sizeof(double));
 	            	realgradwidth[res_ind][coil_ind][seg_ind] =malloc(L_fixedpoints*sizeof(double));
 	            	for (fixed_ind=0; fixed_ind<L_fixedpoints; fixed_ind++) {
			    gradSigma[res_ind][coil_ind][seg_ind][fixed_ind] = calloc(allparams.n_diff, sizeof(double));
			    gradwidth[res_ind][coil_ind][seg_ind][fixed_ind] = calloc(allparams.n_diff,sizeof(double));
			    realgradSigma[res_ind][coil_ind][seg_ind][fixed_ind] = calloc(allparams.n_diff, sizeof(double));
			    realgradwidth[res_ind][coil_ind][seg_ind][fixed_ind] = calloc(allparams.n_diff,sizeof(double));
 	            	}
 	            	printf("seg_ind=%d/%d\n", seg_ind, allparams.intparams[coil_ind]);
 	            	solve_gradRes(gradRes[res_ind][coil_ind][seg_ind], Bfield_island[res_ind], periodicfield_saved[res_ind], lambdamu_saved, L_fixedpoints, N_gridphi_fieldperiod, allparams, coil_ind, seg_ind) ;
 	            	printf("Entering function that solves for the gradient of the circumference\n");
 	            	solve_gradcirc(gradcirc[res_ind][coil_ind][seg_ind], Bfield_island[res_ind], ext_periodicfield[res_ind], &lambda, L_fixedpoints, N_gridphi_fieldperiod, allparams, coil_ind, seg_ind) ;
 	            	printf("<--\n");
 	            	printf("Entering function that solves for the gradient of Σ\n");
 	            	solve_gradtangent(gradSigma[res_ind][coil_ind][seg_ind], Bfield_island[res_ind], ext_periodicfield[res_ind], lambdaQ, mu, L_fixedpoints, q0_index[res_ind], N_gridphi_fieldperiod, allparams, coil_ind, seg_ind) ;
 	            	printf("<--\n");
 	            	for (param_ind = 0; param_ind<allparams.n_diff; param_ind++) {
 	            	    //gradRes[param_ind]*=(-1.0/4);
 	            	    for (fixed_ind=0; fixed_ind<L_fixedpoints; fixed_ind++) {
 	            	    	gradwidth[res_ind][coil_ind][seg_ind][fixed_ind][param_ind] = ext_periodicfield[res_ind][fixed_ind].width * ( gradcirc[res_ind][coil_ind][seg_ind][param_ind]/ext_periodicfield[res_ind][0].circ - gradSigma[res_ind][coil_ind][seg_ind][fixed_ind][param_ind]/ext_periodicfield[res_ind][fixed_ind].Sigma );
 	            	    	printf("gradRes[%d][%d][%d][%d] = %f\n", res_ind, coil_ind, seg_ind, param_ind, gradRes[res_ind][coil_ind][seg_ind][param_ind]);
 	            	    	printf("gradcirc[%d][%d][%d][%d] = %f\n", res_ind, coil_ind, seg_ind, param_ind, gradcirc[res_ind][coil_ind][seg_ind][param_ind]);
 	            	    	printf("gradSigma[%d][%d][%d][%d][%d] = %f\n", res_ind, coil_ind, seg_ind, param_ind, fixed_ind, gradSigma[res_ind][coil_ind][seg_ind][fixed_ind][param_ind]);
 	            	    	printf("gradwidth[%d][%d][%d][%d][%d] = %f\n", res_ind, coil_ind, seg_ind, param_ind, fixed_ind, gradwidth[res_ind][coil_ind][seg_ind][fixed_ind][param_ind]);
 	            	    	realgradRes[res_ind][coil_ind][seg_ind][param_ind] = gradRes[res_ind][coil_ind][seg_ind][param_ind];
 	            	    	realgradcirc[res_ind][coil_ind][seg_ind][param_ind] = gradcirc[res_ind][coil_ind][seg_ind][param_ind];
 	            	    	realgradwidth[res_ind][coil_ind][seg_ind][fixed_ind][param_ind] = gradwidth[res_ind][coil_ind][seg_ind][fixed_ind][param_ind];
 	            	    	realgradSigma[res_ind][coil_ind][seg_ind][fixed_ind][param_ind] = gradSigma[res_ind][coil_ind][seg_ind][fixed_ind][param_ind];
 	            	    }
 	            	}
 	            	if (strncmp(allparams.type, "coil", 4) == 0) {
			    for (row=0; row<3;row++) {
			    	if ( (seg_ind != allparams.intparams[coil_ind]-1) && (seg_ind != 0) ) 
			    	dr[row] = coils[coil_ind][seg_ind+1][row] - coils[coil_ind][seg_ind-1][row];
			    	else if (seg_ind == 0) 
			    	dr[row] = coils[coil_ind][1][row] - coils[coil_ind][allparams.intparams[coil_ind]-1][row];
				else if (seg_ind == allparams.intparams[coil_ind]-1)	
			    	dr[row] = coils[coil_ind][0][row] - coils[coil_ind][allparams.intparams[coil_ind]-2][row];

				dr[row]*=0.5;

			    	//printf("dr = %f\n", dr[row]);
			    	if (seg_ind != 0) 
			    	drdown[row] = - coils[coil_ind][seg_ind-1][row] + coils[coil_ind][seg_ind][row];
			    	else 
			    	drdown[row] = - coils[coil_ind][allparams.intparams[coil_ind]-1][row] + coils[coil_ind][seg_ind][row];
			    	//printf("drdown = %f\n", drdown[row]);
			    }
			    dldown = pow(drdown[0]*drdown[0] + drdown[1]*drdown[1] + drdown[2]*drdown[2], 0.5);
			    dl = pow(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2], 0.5);
			    //factor = dldown/(dl+dldown);
			    //factor = 2.0/allparams.intparams[coil_ind];
			    printf("dl = %f, dlc = %f\n", dl, dldown);
			    for (row = 0; row<3; row++) {
			    	//realgradcirc[row] = factor*(dr[(row+1)%3]*gradcirc[(row+2)%3] - dr[(row+2)%3]*gradcirc[(row+1)%3]) + (1.0 - factor)*(drdown[(row+1)%3]*gradcirc[(row+2)%3] - drdown[(row+2)%3]*gradcirc[(row+1)%3]);
			    	realgradcirc[res_ind][coil_ind][seg_ind][row] = dr[(row+1)%3]*gradcirc[res_ind][coil_ind][seg_ind][(row+2)%3] - dr[(row+2)%3]*gradcirc[res_ind][coil_ind][seg_ind][(row+1)%3];
			    	realgradRes[res_ind][coil_ind][seg_ind][row] = dr[(row+1)%3]*gradRes[res_ind][coil_ind][seg_ind][(row+2)%3] - dr[(row+2)%3]*gradRes[res_ind][coil_ind][seg_ind][(row+1)%3];
			    	printf("realgradcirc_%d = %f\n", row, realgradcirc[res_ind][coil_ind][seg_ind][row]);
			    	printf("realgradRes_%d = %f\n", row, realgradRes[res_ind][coil_ind][seg_ind][row]);
			    }
			    for (fixed_ind=0; fixed_ind<L_fixedpoints; fixed_ind++) {
			    	for (row = 0; row<3; row++) {
				    realgradSigma[res_ind][coil_ind][seg_ind][fixed_ind][row] = dr[(row+1)%3]*gradSigma[res_ind][coil_ind][seg_ind][fixed_ind][(row+2)%3] - dr[(row+2)%3]*gradSigma[res_ind][coil_ind][seg_ind][fixed_ind][(row+1)%3];
				    realgradwidth[res_ind][coil_ind][seg_ind][fixed_ind][row] = dr[(row+1)%3]*gradwidth[res_ind][coil_ind][seg_ind][fixed_ind][(row+2)%3] - dr[(row+2)%3]*gradwidth[res_ind][coil_ind][seg_ind][fixed_ind][(row+1)%3];
				    //printf("realgradSigma[%d][%d][%d][%d][%d] = %f\n", res_ind, coil_ind, seg_ind, fixed_ind, param_ind, realgradSigma[res_ind][coil_ind][seg_ind][fixed_ind][row]);
				    //printf("realgradwidth[%d][%d][%d][%d][%d] = %f\n", res_ind, coil_ind, seg_ind, fixed_ind, param_ind, realgradwidth[fixed_ind][row]);
			    	}
			    }
 	            	}
 	            	//fprintf(file_gradRes, "%f %f %f\n", realgradRes[res_ind][coil_ind][seg_ind][0], realgradRes[res_ind][coil_ind][seg_ind][1], realgradRes[res_ind][coil_ind][seg_ind][2]);
 	            	//fprintf(file_gradcirc, "%f %f %f\n", realgradcirc[res_ind][coil_ind][seg_ind][0], realgradcirc[res_ind][coil_ind][seg_ind][1], realgradcirc[res_ind][coil_ind][seg_ind][2]);
 /* 
Note for now the index 0 below because I only print to file gradients with respect to one magnetix island. Perhaps zero should be replaced with the index corresponding to the largest island? Or maybe it should be the gradient of the sum of island widths?
 */
 	            	//fprintf(file_gradSigma, "%f %f %f\n", realgradSigma[res_ind][coil_ind][seg_ind][0][0], realgradSigma[res_ind][coil_ind][seg_ind][0][1], realgradSigma[res_ind][coil_ind][seg_ind][0][2]);
			for (param_ind=0; param_ind< allparams.n_diff; param_ind++) 
 	            	fprintf(file_gradwidth, "%f ", realgradwidth[res_ind][coil_ind][seg_ind][0][param_ind]/ext_periodicfield[res_ind][0].width);
 	            	fprintf(file_gradwidth, "\n");
 	            	//fprintf(file_gradwidth, "%f %f %f\n", realgradwidth[res_ind][coil_ind][seg_ind][0][0]/ext_periodicfield[res_ind][0].width, realgradwidth[res_ind][coil_ind][seg_ind][0][1]/ext_periodicfield[res_ind][0].width, realgradwidth[res_ind][coil_ind][seg_ind][0][2]/ext_periodicfield[res_ind][0].width);
 	            }
 	            //fprintf(file_gradRes, "\n");
 	            //fprintf(file_gradcirc, "\n");
 	            //fprintf(file_gradSigma, "\n");
 	            fprintf(file_gradwidth, "\n");
 	        }
 	        clock_t int5 = clock();
 	        printf("time after solving for shape gradient: %f\n", (double) (int5-start)/CLOCKS_PER_SEC);
 	        
 	        if (checshap == 'y') {
		    FILE *file_checkgradwidth;
		    if ( (file_checkgradwidth = fopen("checkgrad_width.txt", "w") ) == NULL ) {	
		        printf("Error: couldn't open checkgrad_width.txt\n");
		        exit(EXIT_FAILURE);
		    }	
 	            checknormgradSigma[res_ind] = malloc(max_coil_ind*sizeof(double));
		    checknormgradwidth[res_ind] = malloc(max_coil_ind*sizeof(double));
		    checknormgradcirc[res_ind] = malloc(max_coil_ind*sizeof(double));
		    checknormgradRes[res_ind] = malloc(max_coil_ind*sizeof(double));
		    for (coil_ind=0; coil_ind<max_coil_ind; coil_ind++) {
		        if (strncmp(allparams.type, "coil", 4) == 0) num_segs = allparams.intparams[coil_ind];
		        else num_segs = 1; // for anything but the "coil" type, the 2nd index of diff_params is always 0
			printf("num_segs = %d\n", num_segs);
		        checknormgradSigma[res_ind][coil_ind] = malloc(num_segs*sizeof(double));
		        checknormgradwidth[res_ind][coil_ind] = malloc(num_segs*sizeof(double));
		        checknormgradcirc[res_ind][coil_ind] = malloc(num_segs*sizeof(double));
		        checknormgradRes[res_ind][coil_ind] = malloc(num_segs*sizeof(double));
//
//num_segs = max_num_segs;
//num_segs = max_num_segs + seg_start;
		        for (seg_ind=seg_start; seg_ind<num_segs; seg_ind+=seg_step) {
		    	    checknormgradSigma[res_ind][coil_ind][seg_ind] = malloc(L_fixedpoints*sizeof(double));
		    	    checknormgradwidth[res_ind][coil_ind][seg_ind] = malloc(L_fixedpoints*sizeof(double));
		    	    for (fixed_ind = 0; fixed_ind < L_fixedpoints; fixed_ind++) {
		    	        checknormgradSigma[res_ind][coil_ind][seg_ind][fixed_ind] = malloc(allparams.n_diff*sizeof(double));
		    	        checknormgradwidth[res_ind][coil_ind][seg_ind][fixed_ind] = malloc(allparams.n_diff*sizeof(double));
		    	    }
		    	    checknormgradcirc[res_ind][coil_ind][seg_ind] = malloc(allparams.n_diff*sizeof(double));
		    	    checknormgradRes[res_ind][coil_ind][seg_ind] = malloc(allparams.n_diff*sizeof(double));
		    	    for (param_ind=0; param_ind<allparams.n_diff; param_ind++) {
		    	    	for (fixed_ind = 0; fixed_ind < L_fixedpoints; fixed_ind++) {
		    	            checknormgradSigma[res_ind][coil_ind][seg_ind][fixed_ind][param_ind] = malloc(delta_ind_max*sizeof(double));
		    	            checknormgradwidth[res_ind][coil_ind][seg_ind][fixed_ind][param_ind] = malloc(delta_ind_max*sizeof(double));
		    	        }
		    	        checknormgradcirc[res_ind][coil_ind][seg_ind][param_ind] = malloc(delta_ind_max*sizeof(double));
		    	        checknormgradRes[res_ind][coil_ind][seg_ind][param_ind] = malloc(delta_ind_max*sizeof(double));
		    	    	delta = pow(10.0, powdeltain);
				if (convergence_test == 0) {
				    delta_ind_max = 1;
				    delta = 0.0001;
				}
		    	        for (delta_ind =0; delta_ind < delta_ind_max; delta_ind++) {
		    	    	    if (param_ind == 0) printf("gradient with respect to 1st parameter\n");
		    	    	    else if (param_ind == 1) printf("gradient with respect to 2nd parameter\n");
		    	    	    else if (param_ind == 2) printf("gradient with respect to 3rd parameter\n");
		    	    	    else printf("gradient with respect to %dth parameter\n", param_ind+1);
		    	    	    //savedcoilseg = allparams.diffparams[0][0][param_ind];
		    	    	    //savedcoilseg = Res / (gradRes[param_ind]);
		    	    	    savedcoilseg = ext_periodicfield[res_ind][0].Res / gradRes[res_ind][coil_ind][seg_ind][param_ind];
		    	    	    //savedcoilseg = 1.0;
		    	    	    if (fabs(savedcoilseg) > 100.0)  savedcoilseg = 1.0; // pick a reference value for finite difference of parameter
		    	    	    printf("FD step = %f\n", delta*savedcoilseg);
// begin the finite differencing procedure by stepping up 
		                    allparams.diffparams[coil_ind][seg_ind][param_ind] += (delta*savedcoilseg);
		                    printf("q0_index = %d\n", q0_index[res_ind]);
// having perturbed the parameter, find the new periodic field line position, the new residue, circumference, Sigma and width
// q0_index+res_ind already points to a non-zero value, so it won't be reassigned below
		    	    	    extsolve_periodicfieldline(ext_periodicfield[res_ind][0].loc, ext_periodicfield_FD[1], periodicfield_saved[res_ind], Bfield_island[res_ind], axis, allparams, m0, L_fixedpoints, pol_mode+res_ind, q0_index+res_ind, N_gridphi_fieldperiod, error_max);
		    	    	    printf("q0_index = %d\n", q0_index[res_ind]);
		    	    	    //printf("q0_index = %d\n", q0_index[res_ind]);
// step back down to original parameter
			            allparams.diffparams[coil_ind][seg_ind][param_ind] -= (delta*savedcoilseg);
			            printf("q0_index = %d\n", q0_index[res_ind]);
// step down for centred difference
			            allparams.diffparams[coil_ind][seg_ind][param_ind] -= (delta*savedcoilseg);
// q0_index+res_ind already points to a non-zero value, so it won't be reassigned here
		    	    	    extsolve_periodicfieldline(ext_periodicfield[res_ind][0].loc, ext_periodicfield_FD[0], periodicfield_saved[res_ind], Bfield_island[res_ind], axis, allparams, m0, L_fixedpoints, pol_mode+res_ind, q0_index+res_ind, N_gridphi_fieldperiod, error_max);
		    	    	    checknormgradcirc[res_ind][coil_ind][seg_ind][param_ind][delta_ind] = ( log(ext_periodicfield_FD[1][0].circ) - log(ext_periodicfield_FD[0][0].circ) ) / (2.0*delta*savedcoilseg);
		    	    	    printf("q0_index = %d\n", q0_index[res_ind]);
		    	    	    checknormgradRes[res_ind][coil_ind][seg_ind][param_ind][delta_ind] = ( log(fabs(ext_periodicfield_FD[1][0].Res)) - log(fabs(ext_periodicfield_FD[0][0].Res)) ) / (2.0*delta*savedcoilseg);
// step back up to original parameter
		    	    	    allparams.diffparams[coil_ind][seg_ind][param_ind] += (delta*savedcoilseg);
//printf("res_ind=%d, coil_ind=%d, seg_ind = %d, param_ind = %d\n", res_ind, coil_ind, seg_ind, param_ind);
		    	    	    printf("gradRes[%d][%d][%d][%d] = %9.9f (direct) %9.9f (FD)\n", res_ind, coil_ind, seg_ind, param_ind, realgradRes[res_ind][coil_ind][seg_ind][param_ind]/ext_periodicfield[res_ind][0].Res, checknormgradRes[res_ind][coil_ind][seg_ind][param_ind][delta_ind]);
		    	    	    printf("gradcirc[%d][%d][%d][%d] = %9.9f (direct) %9.9f (FD)\n", res_ind, coil_ind, seg_ind, param_ind, realgradcirc[res_ind][coil_ind][seg_ind][param_ind]/ext_periodicfield[res_ind][0].circ, checknormgradcirc[res_ind][coil_ind][seg_ind][param_ind][delta_ind]);
		    	    	    for (fixed_ind=0; fixed_ind<L_fixedpoints; fixed_ind++) {
		    	    	        checknormgradSigma[res_ind][coil_ind][seg_ind][fixed_ind][param_ind][delta_ind] = ( log(ext_periodicfield_FD[1][fixed_ind].Sigma) - log(ext_periodicfield_FD[0][fixed_ind].Sigma) ) / (2.0*delta*savedcoilseg);
		    	    	        checknormgradwidth[res_ind][coil_ind][seg_ind][fixed_ind][param_ind][delta_ind] = ( log(ext_periodicfield_FD[1][fixed_ind].width) - log(ext_periodicfield_FD[0][fixed_ind].width) ) / (2.0*delta*savedcoilseg);
		    	    	        //checknormgradSigma[fixed_ind][param_ind] = ( newSigma[1][fixed_ind] - Sigma[fixed_ind] ) / (delta*savedcoilseg);
		    	    	        //checknormgradwidth[fixed_ind][param_ind] = ( newwidth[1][fixed_ind] - width[fixed_ind] ) / (delta*savedcoilseg);
		    	    	        if (fixed_ind==0) {
		    	    	    	printf("gradSigma[%d][%d]%d]%d[%d] = %9.9f (direct) %9.9f (FD)\n", res_ind, coil_ind, seg_ind, param_ind, fixed_ind, realgradSigma[res_ind][coil_ind][seg_ind][fixed_ind][param_ind]/ext_periodicfield[res_ind][fixed_ind].Sigma, checknormgradSigma[res_ind][coil_ind][seg_ind][fixed_ind][param_ind][delta_ind]);   
		    	    	    	printf("gradwidth[%d][%d]%d[%d][%d] = %9.9f (direct) %9.9f (FD)\n", res_ind, coil_ind, seg_ind, param_ind, fixed_ind, realgradwidth[res_ind][coil_ind][seg_ind][fixed_ind][param_ind]/ext_periodicfield[res_ind][fixed_ind].width, checknormgradwidth[res_ind][coil_ind][seg_ind][fixed_ind][param_ind][delta_ind]);
		    	    	        }
		    	            
		    	    	    }
				    delta *= (pow(10.0, steppowdelta/delta_ind_max));
		    	        }
			    	fprintf(file_checkgradwidth, "%f ", checknormgradwidth[res_ind][coil_ind][seg_ind][0][param_ind][delta_ind_max-1]);
			    	//fprintf(file_gradwidth, "%f %f %f\n", checknormgradwidth[res_ind][coil_ind][seg_ind][0][0], checknormgradwidth[res_ind][coil_ind][seg_ind][0][1], checknormgradwidth[res_ind][coil_ind][seg_ind][0][2]);
			    }
			    fprintf(file_checkgradwidth, "\n");
			    if (strncmp(allparams.type, "coil", 4) != 0) {
		                printf("delta = [");
		                delta = pow(10.0, powdeltain);
		                for (delta_ind = 0; delta_ind<delta_ind_max; delta_ind++) {
		                	printf("%f ", delta);
		                	if (delta_ind != delta_ind_max - 1) printf(", ");
		                	delta *= (pow(10.0, steppowdelta/delta_ind_max));
		                }
		                printf("]\n");
		                for (param_ind=0; param_ind<allparams.n_diff; param_ind++) {
		                	normgradwidth = gradwidth[res_ind][coil_ind][seg_ind][0][param_ind]/ext_periodicfield[res_ind][0].width;
		                	normgradcirc = gradcirc[res_ind][coil_ind][seg_ind][param_ind]/ext_periodicfield[res_ind][0].circ;
/* 
note aga                        in the index 0
*/
		                	normgradSigma = gradSigma[res_ind][coil_ind][seg_ind][0][param_ind]/ext_periodicfield[res_ind][0].Sigma;
		                	normgradRes = gradRes[res_ind][coil_ind][seg_ind][param_ind]/ext_periodicfield[res_ind][0].Res; 
		                	printf("gw[%d] = %10.10f\ngw_delta[%d] = [", param_ind, normgradwidth, param_ind);
		                	for (delta_ind = 0; delta_ind<delta_ind_max; delta_ind++) {
		                	    printf("%10.10f ", checknormgradwidth[res_ind][coil_ind][seg_ind][0][param_ind][delta_ind]);
		                	    if (delta_ind != delta_ind_max - 1) printf(", ");
		                	}
		                	printf("]\n");
		                	printf("gC[%d] = %10.10f\ngC_delta[%d] = [", param_ind, normgradcirc, param_ind);
		                	for (delta_ind = 0; delta_ind<delta_ind_max; delta_ind++) {
		                        printf("%10.10f ", checknormgradcirc[res_ind][coil_ind][seg_ind][param_ind][delta_ind]);
		                        if (delta_ind != delta_ind_max - 1) printf(", ");
		                	}
		                	printf("]\n");
		                	printf("gSigma[%d] = %10.10f\ngSigma_delta[%d] = [", param_ind, normgradSigma, param_ind);
		                	for (delta_ind = 0; delta_ind<delta_ind_max; delta_ind++) {
		                       printf("%10.10f", checknormgradSigma[res_ind][coil_ind][seg_ind][0][param_ind][delta_ind]);
		                       if (delta_ind != delta_ind_max - 1) printf(", ");
		                    
		                	}
		                	printf("]\n");
		                	printf("gRes[%d] = %10.10f\ngRes_delta[%d] = [", param_ind, normgradRes, param_ind);
		                	for (delta_ind = 0; delta_ind<delta_ind_max; delta_ind++) {
		                        printf("%10.10f", checknormgradRes[res_ind][coil_ind][seg_ind][param_ind][delta_ind]);
		                        if (delta_ind != delta_ind_max - 1) printf(", ");
		                		
		                	}
		                	printf("]\n");
		                }
			    }
		        }
			fprintf(file_checkgradwidth, "\n");
		    }
 	        }
	    }
	}
    	else { 
    		printf("WARNING: periodic field line indexed %d not found\nthis corresponds to the guess in line %d of the file fixedpoint.txt (convention: first line = 1)\n", res_ind, res_ind+2); 
    	}
    }
    fprintf(island_Reiman, "%.14e %.14e %.14e %.14e %.14e\n", allparams.diffparams[0][0][2], ext_periodicfield[0][0].Res, ext_periodicfield[0][0].width, ext_periodicfield[0][0].Sigma, ext_periodicfield[0][0].circ);
    if (allparams.diffparams[0][0][2] > 1e-10) allparams.diffparams[0][0][2]*= pow(10.0,-0.05);
    else study_Reim = 0; 
    printf("epsilon = %f\n", allparams.diffparams[0][0][2]); 
    } while (study_Reim == 1); 
    fclose(island_Reiman);
    clock_t int6 = clock();
    printf("... %f ...\n", (double) (int6-start)/CLOCKS_PER_SEC);
    return 0;
}
