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
	int N_gridphi_fieldperiod=100, m0_fieldperiods=1, N_gridphi_tor, num_turns=1;
	int *n_coils=NULL, **n_segs=NULL, qq=1, *qq_segs=NULL, qq_segsq=1, pol_mode = 6, tor_mode = 1, ind;
	int make_Poincare_iota=2, find_axis = 1, find_islands=2, n_points=25; 
	int L_fixedpoints, N_line;
	int coil_ind, seg_ind;
	struct position *Xp=calloc(1,sizeof(struct position)), *axis, *island_center, **mu, **lambdaQ, lambda;
	struct ext_position *ext_center; //, *grad_ext_center;
	struct field **Bfield_island, **Bfield_axis;
	//struct field *BB, *gradBB; for field checks
	double ***coils=NULL, *width, *number, ***shapecirc, *gradtangent, *Reimparam=NULL;  //*gradwidth,
	double r_interval = 0.1, *iota=NULL, *minor_radius=NULL, axis_phi0[2];
	double **evec=malloc(2*sizeof(double*)), *eval=malloc(2*sizeof(double)), trace, det, iota_axis;
	char *type = "Reim";

	if (strncmp(type, "coil", 4) == 0) {
		m0_fieldperiods = 3; //NCSX
		N_gridphi_fieldperiod = 20;
		pol_mode = 4; tor_mode = 1;
	}
	N_gridphi_tor = N_gridphi_fieldperiod*m0_fieldperiods,

	evec[0]=malloc(2*sizeof(double)); evec[1]=malloc(2*sizeof(double));
	/* make arrays of coils */
	printf("field periods = %d\n", m0_fieldperiods);
	printf("N points per field period = %d\n", N_gridphi_fieldperiod);
	n_coils = &qq; qq_segs = &qq_segsq; n_segs = &qq_segs; // Initializes all pointers
	coils = coil_grid(n_coils, n_segs);
	clock_t tcoils = clock();
	printf("Time after creating arrays of coils: %f\n", (double) (tcoils-start)/CLOCKS_PER_SEC);
	Xp->loc[0] = 1.002; Xp->loc[1]= 0.0;
	Xp->tangent = malloc(2*sizeof(double*));
	Xp->tangent[0] = malloc(2*sizeof(double)); Xp->tangent[1] = malloc(2*sizeof(double));
	/* find magnetic axis */
	//Bpoint = Bfield(Xp->loc, varphi, coils, n_coils, n_segs);
	//printf("B=(%f, %f, %f)\n", Bpoint->value[0], Bpoint->value[1], Bpoint->value[2]);
	printf("Do you want to produce a Poincare plot and an iota profile? (1=yes; 0=no) ");
	scanf("%d", &make_Poincare_iota);
	printf("Do you want to find magnetic islands? (1=yes; 0=no) ");
	scanf("%d", &find_islands);
	if (find_axis == 1) {
		//Xp->loc[0] = 1.55; Xp->loc[1]= 0.0;
		Xp->loc[0] = 1.1; Xp->loc[1]= 0.0;
		Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		Bfield_axis = malloc(N_gridphi_tor*sizeof(struct field));
		for (ind=0; ind<N_gridphi_tor; ind++) {
			Bfield_axis[ind] = malloc(4*sizeof(struct field));
		}
		axis = solve_magneticaxis(coils, *n_coils, *n_segs, Bfield_axis, Xp, N_gridphi_tor);
		axis_phi0[0] = axis[0].loc[0];  axis_phi0[1] = axis[0].loc[1];
		linalg2x2(Xp->tangent, evec, eval, &det, &trace);
		//printf("evec=%f\n", evec[0][0]);
		iota_axis = eval[1]/(2*M_PI);
		printf("iota_axis=%f\n", iota_axis);
		printf("unity?%f\n", eval[0]);
	}
	/* make Poincare plots (optional) */
	if (make_Poincare_iota==1) {
		//Xp->loc[0] = 1.09; Xp->loc[1]=0.0;
		Xp->loc[0] = 1.13; Xp->loc[1]=0.1;
		//Xp->loc[0] = 1.206; Xp->loc[1]=0.0;
		Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		//r_interval = 0.001;
		//n_points = 10;
		//printf("center is (%f, %f)\n", Xp->loc[0], Xp->loc[1]);
		minor_radius = malloc(n_points*sizeof(double));
		iota = malloc(n_points*sizeof(double));
		printf("About to enter Poincare module\n");
		iotaprofile(Xp->loc[0], r_interval, n_points, m0_fieldperiods, N_gridphi_fieldperiod, minor_radius, iota, coils, *n_coils, *n_segs);
		clock_t tPoincare = clock();
		printf("Time after filling Poincare plot file: %f\n", (double) (tPoincare-start)/CLOCKS_PER_SEC);
	}

	if (find_islands == 1) {
		//Xp->loc[0] = 1.3; Xp->loc[1]= -0.515; // island is at(1.299008, -0.515194)
		Xp->loc[0] = 1.299009; Xp->loc[1]= -0.515194; // NCSX: island is at(1.299008, -0.515194)
		//Xp->loc[0] = 0.945; Xp->loc[1]= 0.0;//Dommaschk (5,2) amp 1.73: island is at(1.299008, -0.515194)
		//Xp->loc[0] = 1.091; Xp->loc[1]= 0.0; // Dommaschk (5,2) amp 0.00001 on top of loop current
		//Xp->loc[0] = 1.0097; Xp->loc[1]= 0.0; // Dommaschk (5,2) amp 0.00001 on top of loop current
		Xp->loc[0] = 1.21; Xp->loc[1]= 0.0;  // Reiman 6-fold island chain
		//Xp->loc[0] = 1.21016105; Xp->loc[1]= 0.0; ???
		//Xp->loc[0] = 1.0; Xp->loc[1]= 0.21;  ???
		Xp->tangent[0][0]=1.0; Xp->tangent[0][1]=0.0; Xp->tangent[1][0]=0.0; Xp->tangent[1][1]=1.0;
		if (tor_mode % m0_fieldperiods == 0) 
			L_fixedpoints = pol_mode;	
		else 				 
			L_fixedpoints = m0_fieldperiods*pol_mode;
		N_line = L_fixedpoints*N_gridphi_fieldperiod;
		Bfield_island = malloc(N_line*sizeof(struct field));
		for (ind=0; ind<N_line; ind++) {
			Bfield_island[ind] = malloc(4*sizeof(struct field));
		}
		island_center = solve_islandcenter(coils, *n_coils, *n_segs, Bfield_island, Xp, m0_fieldperiods, L_fixedpoints, N_gridphi_fieldperiod);
		clock_t int1 = clock();
		printf("time after finding island center: %f\n", (double) (int1-start)/CLOCKS_PER_SEC);
		//printstruct("Xp", Xp);
		//printstruct("lambda", lambda);
		num_turns=1;
		ext_center = solve_islandcenter_full(island_center[0].loc, axis_phi0, &num_turns, m0_fieldperiods, L_fixedpoints, N_gridphi_fieldperiod, Bfield_island);
		clock_t int2 = clock();
		printf("time after following field line: %f\n", (double) (int2-start)/CLOCKS_PER_SEC);
		width = calc_islandwidth(ext_center, m0_fieldperiods, pol_mode, L_fixedpoints, N_gridphi_fieldperiod);
		printf("time after calculating width: %f\n", (double) (int2-start)/CLOCKS_PER_SEC);
		lambda = solve_lambda_circ(ext_center, m0_fieldperiods, L_fixedpoints, N_gridphi_fieldperiod);
		printf("time after following lambdacirc: %f\n", (double) (int2-start)/CLOCKS_PER_SEC);
		printf("%f\n", Bfield_island[0][0].value[0]);
		number = solve_gradcirc(Bfield_island, ext_center, &lambda, m0_fieldperiods, L_fixedpoints, N_gridphi_fieldperiod, Reimparam) ;
		//clock_t int3 = clock();
		//printf("time after following adjoint variable for circumference: %f\n", (double) (int3-start)/CLOCKS_PER_SEC);
		//shapecirc = shapecircumference(coils, *n_coils, *n_segs, ext_center, &lambda, m0_fieldperiods, N_gridphi_fieldperiod, tor_mode, pol_mode) ;
		//printf("time after calculating shape gradient of circumference: %f\n", (double) (int4-start)/CLOCKS_PER_SEC);
	
		if (strncmp(type, "coil", 4) == 0) {
			FILE *file_shapecirc;
			if ( (file_shapecirc = fopen("shape_circ.txt", "w") ) == NULL ) {	
				printf("Oups: couldn't open shape_circ.txt\n");
				exit(EXIT_FAILURE);
			}	
			shapecirc = malloc((*n_coils)*sizeof(double));

			for (coil_ind=0; coil_ind<(*n_coils); coil_ind++) {
				shapecirc[coil_ind] = malloc((*n_segs)[coil_ind]*sizeof(double));
				printf("coil_ind=%d/%d\n", coil_ind, *n_coils);
				for ( seg_ind=0; seg_ind < (*n_segs)[coil_ind] ; seg_ind++) {
					printf("seg_ind=%d/%d\n", seg_ind, (*n_segs)[coil_ind]);
					shapecirc[coil_ind][seg_ind] = solve_gradcirc(Bfield_island, ext_center, &lambda, m0_fieldperiods, L_fixedpoints, N_gridphi_fieldperiod, coils[coil_ind][seg_ind]) ;
					fprintf(file_shapecirc, "%f %f %f\n", shapecirc[coil_ind][seg_ind][0], shapecirc[coil_ind][seg_ind][1], shapecirc[coil_ind][seg_ind][2]);
				}
				fprintf(file_shapecirc, "\n");
			}

		}

		mu = solve_mu_tangent(ext_center, m0_fieldperiods, L_fixedpoints, N_gridphi_fieldperiod);
		lambdaQ = solve_lambda_tangent(Bfield_island, ext_center, mu, m0_fieldperiods, L_fixedpoints, N_gridphi_fieldperiod);
		clock_t int4 = clock();
		printf("time after following adjoint variables for tangent map: %f\n", (double) (int4-start)/CLOCKS_PER_SEC);
		gradtangent = solve_gradtangent(Bfield_island, ext_center, lambdaQ, mu, m0_fieldperiods, L_fixedpoints, N_gridphi_fieldperiod, Reimparam) ;
		clock_t int5 = clock();
		printf("time after solving for gradient of tangent map: %f\n", (double) (int5-start)/CLOCKS_PER_SEC);
		//// The module below contains everything: island width and island width adjoint gradient
		//// However, I think it ruins modularity
		coil_ind = seg_ind = 0;
		//width = islandwidthnew(Bfield_island, ext_center, &lambda, lambdaQ, mu, m0_fieldperiods, pol_mode, L_fixedpoints, N_gridphi_fieldperiod, coils[coil_ind][seg_ind]);
		//clock_t int6 = clock();
		//printf("time after solving for everything: %f\n", (double) (int6-start)/CLOCKS_PER_SEC);
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
