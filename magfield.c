//Author: Alessandro Geraldini

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include "isc.h"
#define STANDSIZE 200

struct fieldparams fetchparams() {
	struct fieldparams allparams;
	int lenline = 1000, row=0, rowst=0, ctr = 0, ii=0;
	char line[lenline];
	int nn, coil_index=0, coilseg_index=0, startstore = 0, count=0;
	char *intstr = malloc(6*sizeof(char));
	FILE *file_field;

	file_field = fopen("magfieldparams.txt", "r");
	if (file_field == NULL) {	
		printf("Cannot open magfieldparams.txt\n");
		exit(1);
	}
	while (fgets(line, lenline, file_field) != NULL) {	
		if (strncmp(line, "type", 4) == 0) {
			row = 4;
			while (line[row] == 32) {
				row++;
			}
			rowst = row;
			while (line[row] != 10) {
				allparams.type[row-rowst] = line[row];
				row++;
			}
			printf("%s\n", line);
			printf("%s\n", allparams.type);
		}
		else if (strncmp(line, "periods", 7) == 0) {
			printf("%s\n", line);
			row = 7;
			while (line[row] == 32) {
				row++;
			}
			sscanf(line+row, "%2d", &(allparams.m0_fieldperiods));
			printf("m0 = %d\n", allparams.m0_fieldperiods);
		}
		else if (strncmp(line, "ncoils", 6) == 0) {
			printf("%s\n", line);
			row = 6;
			while (line[row] == 32) {
				row++;
			}
			sscanf(line+row, "%2d", &(allparams.n_coils));
			//*num_coils = m0_temp;
			allparams.diffparams =  malloc((allparams.n_coils)*sizeof(double));
			printf("n_coils = %d\n", allparams.n_coils);
		}
		else if (strncmp(line, "ndiff", 5) == 0) {
			printf("%s\n", line);
			row = 5;
			while (line[row] == 32) {
				row++;
			}
			sscanf(line+row, "%2d", &(allparams.n_diff));
			//*num_coils = m0_temp;
			printf("n_diffparams = %d\n", allparams.n_diff);
		}
		else if ( strncmp(line, "const", 5) == 0 ) {
			ctr = 0;
			//n_constparams = malloc(STANDSIZE*sizeof(int));
			allparams.intparams = malloc(STANDSIZE*sizeof(int));
			row= 6;
			do {
				if ( (line[row] != 32) && (line[row] != '\n') ) {
					//printf("line[row]=%c\n", line[row]);
					intstr[count] = line[row];
					startstore = 1;
					count++;
				}
				else {
					if (startstore == 1) {
						//printf("%s\n", &intstr);
						//sscanf(intstr, "%4d", (allparams.intparams)+ctr);
						//printf("%s\n", intstr);
						allparams.intparams[ctr] = atoi(intstr);
						ctr++;
						startstore = 0;
						count=0;
					}
				}
				row++;
			} while (line[row-1] != '\n');
			//n_constparams[1] = 0;
			for (coil_index=0; coil_index<allparams.n_coils; coil_index++) {
				if ( strncmp(allparams.type, "coil", 4) == 0 ) {
					// for a coil type all segments of the coil are considered
					allparams.diffparams[coil_index] = malloc(allparams.intparams[coil_index]*sizeof(double));
				}
				else {
					// for any other type we either have no coils 
					// or coils parameterized in some way
					allparams.diffparams[coil_index] = malloc(sizeof(double));
				}	
			}
		}
		else if ( strncmp(line, "param", 5) == 0 ) {
			ctr = 0;
			//(*num_params) = malloc(2*sizeof(int));
			row= 5;
			nn = 0;
			allparams.constparams = linetodata(line, &nn); 
			//printf("constparams[0][1]=%f\n", allparams.constparams[1]);
			coil_index = coilseg_index = 0;
			//n_constparams[1] = 0;
		}
		else {
			ii=0;
			while (line[ii] != '\n') ii+=1;
			if (ii>0) {
				nn=0;
				if ( strncmp(allparams.type, "coil", 4) == 0 ) {
					// for a coil type all segments of the coil are considered
					nn = 4;
					allparams.diffparams[coil_index][coilseg_index] = linetodata(line, &nn);
					coilseg_index += 1;
					if (coilseg_index ==  allparams.intparams[coil_index]) {
						coilseg_index = 0;
						coil_index++;
					}
				}
				else {
					// for any other type we either have no coils 
					// or coils parameterized in some way
					nn = allparams.n_diff;
					allparams.diffparams[coil_index][0] = linetodata(line, &nn);
					coil_index++;
				}	
			}
		}
	}
	rewind(file_field);
	printf("diffparams=%f\n", allparams.diffparams[0][0][0]);
	return allparams;
}

//double ***fieldparams(char *type, int *m0_symmetry, int *num_coils, int **num_consts, double **num_params) {
//	int lenline = 1000, row=0, rowst=0, ctr = 0;
//	char line[lenline];
//	int nn, coil_index=0, coilseg_index=0, startstore = 0, count=0;
//	char intstr[5];
//	double ***XXvect_coil= NULL;
//	FILE *file_field;
//
//	file_field = fopen("magfieldparams.txt", "r");
//	if (file_field == NULL) {	
//		printf("Cannot open magfieldparams.txt\n");
//		exit(1);
//	}
//	while (fgets(line, lenline, file_field) != NULL) {	
//		if (strncmp(line, "type", 4) == 0) {
//			row = 4;
//			while (line[row] == 32) {
//				row++;
//			}
//			rowst = row;
//			//printf("%s\n", line);
//			//printf("%s\n", type);
//			while (line[row] != 10) {
//				type[row-rowst] = line[row];
//				row++;
//			}
//		}
//		else if (strncmp(line, "periods", 7) == 0) {
//			row = 7;
//			while (line[row] == 32) {
//				row++;
//			}
//			sscanf(line+row, "%2d", m0_symmetry);
//			printf("m0 = %d\n", *m0_symmetry);
//		}
//		else if (strncmp(line, "ncoils", 6) == 0) {
//			row = 6;
//			while (line[row] == 32) {
//				row++;
//			}
//			sscanf(line+row, "%2d", num_coils);
//			//*num_coils = m0_temp;
//			XXvect_coil =  malloc((*num_coils)*sizeof(double));
//			printf("num_coils = %d\n", *num_coils);
//		}
//		else if ( strncmp(line, "const", 5) == 0 ) {
//			ctr = 0;
//			//n_constparams = malloc(STANDSIZE*sizeof(int));
//			*num_consts = malloc(STANDSIZE*sizeof(int));
//			row= 6;
//			do {
//				if ( (line[row] != 32) && (line[row] != '\n') ) {
//					//printf("line[row]=%c\n", line[row]);
//					intstr[count] = line[row];
//					startstore = 1;
//					count++;
//				}
//				else {
//					if (startstore == 1) {
//						//printf("%s\n", &intstr);
//						//sscanf(intstr, "%4d", (*num_consts)+ctr);
//						(*num_consts)[ctr] = atoi(intstr);
//						ctr++;
//						startstore = 0;
//						count=0;
//					}
//				}
//				row++;
//			} while (line[row-1] != '\n');
//			//n_constparams[1] = 0;
//			for (coil_index=0; coil_index<*num_coils; coil_index++) {
//				if ( strncmp(type, "coil", 4) == 0 ) {
//					// for a coil type all segments of the coil are considered
//					XXvect_coil[coil_index] = malloc((*num_consts)[coil_index]*sizeof(double));
//				}
//				else {
//					// for any other type we either have no coils 
//					// or coils parameterized in some way
//					XXvect_coil[coil_index] = malloc(sizeof(double));
//				}	
//			}
//		}
//		else if ( strncmp(line, "param", 5) == 0 ) {
//			ctr = 0;
//			//(*num_params) = malloc(2*sizeof(int));
//			row= 5;
//			*num_params = linetodata(line, &nn); 
//			printf("num_params[0][1]=%f\n", num_params[0][1]);
//			coil_index = coilseg_index = 0;
//			//n_constparams[1] = 0;
//		}
//		else {
//			nn=0;
//			while (line[nn] != '\n') nn+=1;
//			if (nn>0) {
//				nn=0;
//				XXvect_coil[coil_index][coilseg_index] = linetodata(line, &nn);
//				if (coilseg_index ==  (*num_consts)[coil_index]) {
//					coilseg_index = 0;
//					coil_index++;
//				}
//			}
//		}
//	}
//	rewind(file_field);
//				
//			//if ( strncmp(type, "coil", 4) == 0 ) {
//			//	if ((*num_coils) == 0) {
//			//		n_params = 4;
//			//		if (nn > n_params)
//			//		{
//			//			coil_index += 1;
//			//		}
//			//	}
//			//	else {
//			//		XXvect_coil[coil_index][coilseg_index] = malloc(4*sizeof(double));
//			//		for (row=0; row<4; row++) 
//			//		XXvect_coil[coil_index][coilseg_index][row] = *(storevals+row); 
//			//		coilseg_index+=1;
//			//		if (coilseg_index == (*num_consts)[coil_index]) {
//			//			coil_index+=1;
//			//			coilseg_index=0; 
//			//		}
//			//	}
//			//}
//			//else if ( strncmp(type, "Reim", 4) == 0 ) {
//			//	if (isdigit(line[0]) != 0) {
//			//		coil_index += 1;
//			//		n_params = nn;
//			//		XXvect_coil = calloc(1,sizeof(double));
//			//		*XXvect_coil = calloc(1,sizeof(double));
//			//		**XXvect_coil = calloc(n_params,sizeof(double));
//			//		for (row=0;row<n_params;row++) {
//			//			XXvect_coil[0][0][row] = *(storevals+row); 
//			//			//printf("storevals= %f\n", *(storevals+row));
//			//		}
//			//	}
//			//}
//			//else if ( strncmp(type, "heli", 4) == 0 ) {
//			//	//printf("%s\n", line);
//			//	//printf("coil_index=%d\n", coil_index);
//			//	if (isdigit(line[0]) != 0) {
//			//		n_params = nn;
//			//		printf("n_params = %d\n", n_params);
//			//	}
//			//}
//			//}
//			//free(storevals);
//	//num_coils_temp = coil_index;
//	//num_segs_temp = calloc(coil_index,sizeof(int));
//	//max_num_segs = 0;
//	//coil_index = 0;
//	//if (strncmp(type, "coil", 4) == 0) {
//	//	*num_consts = calloc(coil_index,sizeof(int));
//	//	while (fgets(line, lenline, file_field) != NULL)
//	//	{	
//	//		storevals = linetodata(line, &nn);
//	//		if (nn == 4)
//	//		{
//	//			coilseg_index += 1; 
//	//		}
//	//		else if (nn > 4)
//	//		{
//	//			num_segs_temp[coil_index] = coilseg_index;
//	//			coil_index += 1;
//	//			if (coilseg_index > max_num_segs)
//	//			{
//	//				max_num_segs = coilseg_index;
//	//			}
//	//			coilseg_index = 0;
//	//		}
//	//	}
//	//}
//	//else if (strncmp(type, "heli", 4) == 0) {
//	//	printf("n_params = %d\n", n_params);
//	//	*num_consts = calloc(n_params,sizeof(int));
//	//	coil_index = 0;
//	//	while (fgets(line, lenline, file_field) != NULL) {
//	//		storevals = linetodata(line, &nn);
//	//		printf("nn = %d, n_params=%d\n", nn, n_params);
//	//		if (isdigit(line[0]) != 0) {
//	//			XXvect_coil[0][coil_index] = malloc(nn*sizeof(double));
//	//			for (row=0; row< nn; row++) {
//	//				printf("row = %d/%d\n", row, nn);
//	//				XXvect_coil[0][coil_index][row] = storevals[row];	
//	//			}
//	//			coil_index += 1;
//	//		}
//	//	}
//	//}
//	//else {
//	//	max_num_segs = 1;
//	//	//*num_consts = (*num_params);
//	//}
//	//printf("num_coils_temp = %d\n", num_coils_temp);
//	////printf("max_num_segs = %d\n", max_num_segs);
//	//if (strncmp(type, "coil", 4) == 0) {
//	//	coil_index = 0;
//	//	coilseg_index = 0;
//	//	XXvect_coil = calloc(num_coils_temp,sizeof(double));
//	//	for (coil_index = 0; coil_index < num_coils_temp; coil_index ++)
//	//	{
//	//		XXvect_coil[coil_index] = calloc(max_num_segs,sizeof(double));
//	//		for (coilseg_index = 0; coilseg_index < max_num_segs; coilseg_index ++)
//	//		{
//	//				XXvect_coil[coil_index][coilseg_index] = calloc(n_params,sizeof(double));
//	//		}
//	//	}
//	//	coil_index = 0;
//	//	coilseg_index = 0;
//	//	rewind(file_field);
//	//	while (fgets(line, lenline, file_field) != NULL)
//	//	{	
//	//		//printf("coil_index=%d/%d\n", coil_index, num_coils_temp-1);
//	//		storevals = linetodata(line, &nn);
//	//		//printf("nn=%d\n", nn);
//	//		if (nn == n_params)
//	//		{
//	//			//printf("coilseg_index=%d/%d\n", coilseg_index, max_num_segs-1);
//	//			for (row=0;row<n_params;row++) {
//	//				XXvect_coil[coil_index][coilseg_index][row] = *(storevals+row); 
//	//				//printf("XX[%d] = %f\n", row, XXvect_coil[coil_index][coilseg_index][row]);
//	//			}
//	//			coilseg_index += 1; 
//	//		}
//	//		if (nn > n_params)
//	//		{
//	//			coil_index += 1;
//	//			coilseg_index = 0;
//	//		}
//	//		free(storevals);
//	//	}
//	//	fclose(file_field);
//	//	*num_coils = num_coils_temp;
//	//	*num_consts = num_segs_temp;
//	//}
//	//else {
//	//	*num_coils = n_params;
//	//}
//	////printf("num_coils = %d\n", *num_coils);
//	////for (coil_index=0; coil_index< *num_coils; coil_index++) {
//	////	printf("%d\n", (*num_segs)[coil_index]);
//	////}
//	//XXvect_coil[0][0][2] += 0.00;// check shape gradient (first entry)
//	//clock_t int2 = clock();
//	//printf("params = %f %f %f\n", XXvect_coil[0][0][0], XXvect_coil[0][0][1],XXvect_coil[0][0][2]);
//	printf("type = %s\n", type);
//	//printf("Time out of coil module: %f\n", (double) (int2-start)/CLOCKS_PER_SEC);
//	return XXvect_coil;
//}
		
struct field BReim(double RR, double ZZ, double varphi, int m0_symmetry, double iota0, double iota1, double *epsilon, int *k_theta, int size_epsilon) { // all working
		int ind;
		double R_axis = 1.0, theta, rmin, combo, combo1, dcombodR, dcombodZ, dcombo1dR, dcombo1dZ; // check_q;
		double d2combodRR, d2combodRZ, d2combodZZ, d2combo1dRR, d2combo1dRZ, d2combo1dZZ; // check_q;
		//double check, checkdlen, dlen=0.000001; // uncomment when wishing to check gradients
		double quantity, dquantitydZ=0.0;
		struct field mag;
		//printf("epsilon=%f, iota0=%f, iota1=%f\n", epsilon[0], iota0, iota1);
		//printf("k_theta=%d, size_epsilon=%d, m0_symmetry=%d\n", k_theta[0], size_epsilon, m0_symmetry);
		//printf("(RR,ZZ) = (%f, %f)\n", RR, ZZ);
		theta = atan2(ZZ, RR - R_axis);
		rmin = sqrt(pow((RR-R_axis), 2.0) + pow(ZZ, 2.0));
		combo = iota0 + iota1*rmin*rmin;
		//printf("rmin = %f\n", rmin);
		//printf("theta = %f\n", theta);
		combo1 = 0.0;
		dcombodR = 2.0*iota1*(RR - R_axis);
		dcombodZ = 2.0*iota1*ZZ;
		dcombo1dR = 0.0;
		dcombo1dZ = 0.0;
		d2combodRR = 2.0*iota1;
		d2combodRZ = 0.0;
		d2combodZZ = 2.0*iota1;
		d2combo1dRR = 0.0;
		d2combo1dRZ = 0.0;
		d2combo1dZZ = 0.0;
		//epsilon[0] += 0.00001; 
		//this was a test to verify that the tangent map determines the width
		// the small change of the position of the fixed point has no effect
		quantity = 0.0;
		//printf("combo = %f\n", combo);
		//printf("size_epsilon = %d\n", size_epsilon);
		for (ind=0; ind < size_epsilon; ind++) {
			combo  -= k_theta[ind] * epsilon[ind] * pow(rmin, k_theta[ind] - 2) * cos(k_theta[ind]*theta - m0_symmetry*varphi);
			combo1 += k_theta[ind] * epsilon[ind] * pow(rmin, k_theta[ind] - 2) * sin(k_theta[ind]*theta - m0_symmetry*varphi);

			dcombodR -= k_theta[ind] * pow( rmin, k_theta[ind] - 4 ) * epsilon[ind] * ( k_theta[ind] * ZZ * sin( k_theta[ind]*theta - m0_symmetry*varphi )  +  (k_theta[ind] - 2) *  (RR - R_axis) * cos(k_theta[ind]*theta - m0_symmetry*varphi) ) ;
			dcombodZ += pow( rmin, k_theta[ind] - 4 ) * epsilon[ind] * k_theta[ind] * ( k_theta[ind] * sin( k_theta[ind]*theta - m0_symmetry*varphi ) * ( RR - R_axis ) -  (k_theta[ind] - 2) * ZZ * cos(k_theta[ind]*theta - m0_symmetry * varphi) );
			quantity += pow( rmin, k_theta[ind] - 4 ) * epsilon[ind] * k_theta[ind] * ( k_theta[ind] * sin( k_theta[ind]*theta - m0_symmetry*varphi ) * ( RR - R_axis ) );
			dcombo1dR += k_theta[ind]*pow(rmin, k_theta[ind] - 4) * epsilon[ind] * ( - k_theta[ind]* ZZ * cos(k_theta[ind]*theta - m0_symmetry*varphi) + ( k_theta[ind] -2 ) * sin(k_theta[ind]*theta - m0_symmetry*varphi) * (RR - R_axis) ) ;
			dcombo1dZ += k_theta[ind] * pow(rmin, k_theta[ind] - 4) * epsilon[ind] * ( k_theta[ind] * cos(k_theta[ind]*theta - m0_symmetry *varphi) * (RR - R_axis) +  (k_theta[ind] - 2) * sin(k_theta[ind]*theta - m0_symmetry*varphi) * ZZ ) ;

			//right
			d2combodRR -= k_theta[ind]*epsilon[ind]*pow(rmin, k_theta[ind]-6) * ( ( (k_theta[ind]-2)*(rmin*rmin + (k_theta[ind]-4)*pow(RR-R_axis, 2.0)) - pow(k_theta[ind], 2.0)*pow(ZZ, 2.0) ) *cos(k_theta[ind]*theta - m0_symmetry*varphi) + k_theta[ind]*ZZ*(RR-R_axis)*( 2.0*k_theta[ind]-6.0 )*sin(k_theta[ind]*theta - m0_symmetry*varphi) );
			//right 
			d2combodRZ -= k_theta[ind]*epsilon[ind]*pow(rmin, k_theta[ind]-6) * ( k_theta[ind]*( rmin*rmin - (k_theta[ind]-2)*pow(RR-R_axis, 2.0) + (k_theta[ind] - 4)*pow(ZZ, 2.0) ) *sin(k_theta[ind]*theta - m0_symmetry*varphi) + ZZ*(RR-R_axis)*(2*k_theta[ind]*k_theta[ind]-6*k_theta[ind] + 8) * cos(k_theta[ind]*theta - m0_symmetry*varphi)  );
			//right
			d2combodZZ -= k_theta[ind]*epsilon[ind]*pow(rmin, k_theta[ind]-6) * ( ( (k_theta[ind]-2)*(rmin*rmin + (k_theta[ind]-4)*ZZ*ZZ) - pow(k_theta[ind], 2.0)*pow(RR-R_axis, 2.0) ) *cos(k_theta[ind]*theta - m0_symmetry*varphi) - k_theta[ind]*ZZ*(RR-R_axis)*(2*k_theta[ind]-6) * sin(k_theta[ind]*theta - m0_symmetry*varphi)  );

			dquantitydZ += k_theta[ind]*epsilon[ind]*pow(rmin, k_theta[ind]-6) * ( k_theta[ind]*((RR-R_axis)*ZZ*(k_theta[ind]-4)*sin(k_theta[ind]*theta - m0_symmetry*varphi) + k_theta[ind]*pow(RR-R_axis, 2.0)* cos(k_theta[ind]*theta - m0_symmetry*varphi)  ) );

			//right
			d2combo1dRR += k_theta[ind]*epsilon[ind]*pow(rmin, k_theta[ind]-6) * ( ( (k_theta[ind]-2)*(rmin*rmin + (k_theta[ind]-4)*pow(RR-R_axis, 2.0)) - pow(k_theta[ind], 2.0)*pow(ZZ, 2.0) ) *sin(k_theta[ind]*theta - m0_symmetry*varphi) - k_theta[ind]*ZZ*(RR-R_axis)*( 2.0*k_theta[ind]-6.0 )*cos(k_theta[ind]*theta - m0_symmetry*varphi) );
			//right
			d2combo1dRZ += k_theta[ind]*epsilon[ind]*pow(rmin, k_theta[ind]-6) * ( -k_theta[ind]*(rmin*rmin - (k_theta[ind]-2)*pow(RR-R_axis, 2.0) + (k_theta[ind] - 4)*pow(ZZ, 2.0) ) *cos(k_theta[ind]*theta - m0_symmetry*varphi) + ZZ*(RR-R_axis)*(2*k_theta[ind]*k_theta[ind]-6*k_theta[ind] + 8) * sin(k_theta[ind]*theta - m0_symmetry*varphi)  );
			//right
			d2combo1dZZ += k_theta[ind]*epsilon[ind]*pow(rmin, k_theta[ind]-6) * ( ( (k_theta[ind]-2)*(rmin*rmin + (k_theta[ind]-4)*ZZ*ZZ) - pow(k_theta[ind], 2.0)*pow(RR-R_axis, 2.0) ) *sin(k_theta[ind]*theta - m0_symmetry*varphi) + k_theta[ind]*ZZ*(RR-R_axis)*(2*k_theta[ind]-6) * cos(k_theta[ind]*theta - m0_symmetry*varphi)  );
			//printf("combo1 = %f\n", combo1);
			//printf("combo = %f\n", combo);

		}
		//ZZ+=dlen;
		//theta = atan2(ZZ, RR - R_axis);
		//rmin = sqrt(pow((RR-R_axis), 2.0) + pow(ZZ, 2.0));
		//check = dcombo1dR;
		////checkdlen = 2.0*iota1*(RR-R_axis);
		////checkdlen = 2.0*iota1*ZZ;
		//checkdlen = 0.0;
		//for (ind=0; ind < size_epsilon; ind++) {
		//	//checkdlen -= k_theta[ind] * pow( rmin, k_theta[ind] - 4 ) * epsilon[ind] * ( k_theta[ind] * ZZ * sin( k_theta[ind]*theta - m0_symmetry*varphi )  +  (k_theta[ind] - 2) *  (RR - R_axis) * cos(k_theta[ind]*theta - m0_symmetry*varphi) ) ; //dcombodR
		//	checkdlen += k_theta[ind]*pow(rmin, k_theta[ind] - 4) * epsilon[ind] * ( - k_theta[ind]* ZZ * cos(k_theta[ind]*theta - m0_symmetry*varphi) + ( k_theta[ind] -2 ) * sin(k_theta[ind]*theta - m0_symmetry*varphi) * (RR - R_axis) ) ; //dcombo1dR
		//	//checkdlen += pow( rmin, k_theta[ind] - 4 ) * epsilon[ind] * k_theta[ind] * ( k_theta[ind] * sin( k_theta[ind]*theta - m0_symmetry*varphi ) * ( RR - R_axis ) -  (k_theta[ind] - 2) * ZZ * cos(k_theta[ind]*theta - m0_symmetry * varphi) );//dcombodZ
		//	//checkdlen += k_theta[ind] * pow(rmin, k_theta[ind] - 4) * epsilon[ind] * ( k_theta[ind] * cos(k_theta[ind]*theta - m0_symmetry *varphi) * (RR - R_axis) +  (k_theta[ind] - 2) * sin(k_theta[ind]*theta - m0_symmetry*varphi) * ZZ ) ;//dcombo1dZ
		//	//checkdlen += pow( rmin, k_theta[ind] - 4 ) * epsilon[ind] * k_theta[ind] * ( k_theta[ind] * sin( k_theta[ind]*theta - m0_symmetry*varphi ) * ( RR - R_axis ) );//quantity
		//}
		//printf("numerical = %f, exact = %f\n", (checkdlen - check)/dlen, d2combo1dRZ);
		//ZZ-=dlen;
		//theta = atan2(ZZ, RR - R_axis);
		//rmin = sqrt(pow((RR-R_axis), 2.0) + pow(ZZ, 2.0));

		mag.value[0] = ( ZZ / RR ) *  combo + ( ( RR - R_axis ) / RR ) * combo1 ;
		mag.value[1] = - ( (RR - R_axis) / RR ) * combo + ( ZZ / RR ) * combo1  ; 
		mag.value[2] = -1.0;
		mag.derivative[0][0] = ( - ZZ / pow( RR, 2.0) ) *  combo + (ZZ / RR) * dcombodR + combo1 * R_axis / pow(RR, 2.0) + dcombo1dR * (RR - R_axis) / RR;
 //+ gradR*( (2.0*Z.pow(RR, 3.0))*combo0 - (2.0*ZZ/pow(RR, 2.0))*
		mag.derivative[0][1] = ( 1.0 / RR ) * combo + (ZZ / RR) * dcombodZ + dcombo1dZ * (RR - R_axis) / RR;
		mag.derivative[1][0] = ( - R_axis / pow( RR, 2.0) ) *  combo - ((RR - R_axis) / RR) * dcombodR - combo1 * ZZ / pow(RR, 2.0) + dcombo1dR * ZZ / RR;
		mag.derivative[1][1] = - ((RR - R_axis) / RR) * dcombodZ + combo1 * ( 1.0 / RR ) + dcombo1dZ * ZZ / RR ;
		mag.derivative[2][0] = 0.0;
		mag.derivative[2][1] = 0.0;
		mag.twoderivative[0][0][0] = (  2*ZZ / pow( RR, 3.0) ) *combo - (2.0*R_axis/pow(RR, 3.0))*combo1 - (2*ZZ / pow(RR, 2.0)) * dcombodR +  (2*R_axis / pow(RR, 2.0) )*dcombo1dR  + (ZZ/RR)*d2combodRR + ((RR-R_axis)/RR)*d2combo1dRR; 
		mag.twoderivative[0][0][1] = - (  1.0 / pow( RR, 2.0) ) *combo - (ZZ/pow(RR, 2.0))*dcombodZ + (1.0/ RR) * dcombodR +  (R_axis / pow(RR, 2.0) )*dcombo1dZ  + (ZZ/RR)*d2combodRZ + ((RR-R_axis)/RR)*d2combo1dRZ; 
		mag.twoderivative[0][1][0] = - (  1.0 / pow( RR, 2.0) ) *combo - (ZZ/pow(RR, 2.0))*dcombodZ + (1.0/ RR) * dcombodR +  (R_axis / pow(RR, 2.0) )*dcombo1dZ  + (ZZ/RR)*d2combodRZ + ((RR-R_axis)/RR)*d2combo1dRZ; 
		mag.twoderivative[0][1][1] = ( 2.0 / RR ) *dcombodZ +  (ZZ/RR)*d2combodZZ + ((RR-R_axis)/RR)*d2combo1dZZ;
		mag.twoderivative[1][0][0] = (  2*R_axis / pow( RR, 3.0) ) *combo + (2.0*ZZ/pow(RR, 3.0))*combo1 - (2*ZZ / pow(RR, 2.0)) * dcombo1dR -  (2*R_axis / pow(RR, 2.0) )*dcombodR  + (ZZ/RR)*d2combo1dRR - ((RR-R_axis)/RR)*d2combodRR; 
		mag.twoderivative[1][1][0] = mag.twoderivative[1][0][1] = - (  1.0 / pow( RR, 2.0) ) *combo1 - (R_axis/pow(RR, 2.0))*dcombodZ - (ZZ/pow(RR, 2.0))*dcombo1dZ + (1.0/ RR)*dcombo1dR - ((RR - R_axis) / RR)*d2combodRZ  + (ZZ/RR)*d2combo1dRZ;
		mag.twoderivative[1][1][1] = ( 2.0 / RR ) *dcombo1dZ +  (ZZ/RR)*d2combo1dZZ - ((RR-R_axis)/RR)*d2combodZZ;
		mag.twoderivative[2][0][0] = 0.0;
		mag.twoderivative[2][1][0] = mag.twoderivative[2][0][1] = 0.0;
		mag.twoderivative[2][1][1] = 0.0;

		//below checks gradients
		//RR+=dlen;
		//theta = atan2(ZZ, RR - R_axis);
		//rmin = sqrt(pow((RR-R_axis), 2.0) + pow(ZZ, 2.0));
		//check = mag->derivative[1][1];
		//combo = iota0 + iota1*rmin*rmin;
		//combo1 = 0.0;
		//dcombodR = 2.0*iota1*(RR - R_axis);
		//dcombodZ = 2.0*iota1*ZZ;
		//dcombo1dR = 0.0;
		//dcombo1dZ = 0.0;
		//for (ind=0; ind < size_epsilon; ind++) {
		//	combo  -= k_theta[ind] * epsilon[ind] * pow(rmin, k_theta[ind] - 2) * cos(k_theta[ind]*theta - m0_symmetry*varphi);
		//	combo1 += k_theta[ind] * epsilon[ind] * pow(rmin, k_theta[ind] - 2) * sin(k_theta[ind]*theta - m0_symmetry*varphi);

		//	dcombodR -= k_theta[ind] * pow( rmin, k_theta[ind] - 4 ) * epsilon[ind] * ( k_theta[ind] * ZZ * sin( k_theta[ind]*theta - m0_symmetry*varphi )  +  (k_theta[ind] - 2) *  (RR - R_axis) * cos(k_theta[ind]*theta - m0_symmetry*varphi) ) ;
		//	dcombodZ += pow( rmin, k_theta[ind] - 4 ) * epsilon[ind] * k_theta[ind] * ( k_theta[ind] * sin( k_theta[ind]*theta - m0_symmetry*varphi ) * ( RR - R_axis ) -  (k_theta[ind] - 2) * ZZ * cos(k_theta[ind]*theta - m0_symmetry * varphi) );
		//	dcombo1dR += k_theta[ind]*pow(rmin, k_theta[ind] - 4) * epsilon[ind] * ( - k_theta[ind]* ZZ * cos(k_theta[ind]*theta - m0_symmetry*varphi) + ( k_theta[ind] -2 ) * sin(k_theta[ind]*theta - m0_symmetry*varphi) * (RR - R_axis) ) ;
		//	dcombo1dZ += k_theta[ind] * pow(rmin, k_theta[ind] - 4) * epsilon[ind] * ( k_theta[ind] * cos(k_theta[ind]*theta - m0_symmetry *varphi) * (RR - R_axis) +  (k_theta[ind] - 2) * sin(k_theta[ind]*theta - m0_symmetry*varphi) * ZZ ) ;
		//}
		////checkdlen = ( - ZZ / pow( RR, 2.0) ) *  combo + (ZZ / RR) * dcombodR + combo1 * R_axis / pow(RR, 2.0) + dcombo1dR * (RR - R_axis) / RR; // dBRdR
		////checkdlen = ( 1.0 / RR ) * combo + (ZZ / RR) * dcombodZ + dcombo1dZ * (RR - R_axis) / RR; // dBRdZ
		//checkdlen = - ((RR - R_axis) / RR) * dcombodZ + combo1 * ( 1.0 / RR ) + dcombo1dZ * ZZ / RR ; // dBZdZ
		////checkdlen = ( - R_axis / pow( RR, 2.0) ) *  combo - ((RR - R_axis) / RR) * dcombodR - combo1 * ZZ / pow(RR, 2.0) + dcombo1dR * ZZ / RR; // dBZdR
		//printf("numerical = %f, exact = %f\n", (checkdlen - check)/dlen, mag->twoderivative[1][1][0]);
		//RR-=dlen;

		//printf("Exit Reiman magnetic field evaluation: OK\n");
		return mag;
}

struct field Bcoil(double RR, double ZZ, double varphi, double ***coils, int n_params, int *num_segs) {
	/*
		Declarations
	*/
	double Xminxc=0.0, Yminyc=0.0, Zminzc=0.0, Rminrc=0.0, dxc=0.0, dyc=0.0, dzc=0.0, drc, lencoil = 0.0;
	double Xminxcp=0.0, Yminycp=0.0, Zminzcp=0.0, Rminrcp=0.0;
	double XX=RR*cos(varphi), YY=RR*sin(varphi), tancoil[3], tancoilx[3], tancoil_length;
	struct field magfield;// Bvect;// gradB, gradgradBcheck, twogradB;
	//struct field *magfield, *magfieldpdR = calloc(1,sizeof(struct field)), *magfieldpdZ = calloc(1,sizeof(struct field)), *magfieldpdvarphi = calloc(1, sizeof(struct field));
	//struct field *magfieldmdR = calloc(1,sizeof(struct field)), *magfieldmZ = calloc(1,sizeof(struct field));
	//struct field *magfieldmdRpdZ = calloc(1,sizeof(struct field)), *magfieldmdZpdR = calloc(1,sizeof(struct field));
	//struct field *magfieldmdRmdZ = calloc(1,sizeof(struct field)), *magfieldpdZpdR = calloc(1,sizeof(struct field));
	//struct field dmagfielddR, dmagfielddZ, dmagfielddvarphi;
	//struct field d2magfielddRR, d2magfieldZZ;
	//struct field d2magfielddRZ, d2magfieldZR;
	int coil_index=0, coilseg_index=0, row, col, dep;
	double Rminrcvect[3], drcvect[3], Rhat[3], Zhat[3], phihat[3], current;
	double Rminrcvectp[3];
	//double drcvectminus[3], dxcminus=0.0, dycminus=0.0, dzcminus=0.0, drcminus; 
	double object;
	/* 
		Initialize magnetic field elements
	*/
	Rhat[0] = 1.0; Rhat[1] = 0.0; Rhat[2] = 0.0;
	Zhat[0] = 0.0; Zhat[1] = 1.0; Zhat[2] = 0.0;
	phihat[0] = 0.0; phihat[1] = 0.0; phihat[2] = 1.0;
	for (row=0;row<3;row++) {
		magfield.value[row] = 0.0; 
		//Bvect[row] = 0.0; 
		for (col=0;col<2;col++) {
			magfield.derivative[row][col] = 0.0;
			//gradB[row][col] = 0.0;
			//gradBcheck[row][col] = 0.0;
			for (dep=0;dep<2;dep++) {
				magfield.twoderivative[row][col][dep] = 0.0;
				//gradgradB[row][col][dep] = 0.0;
				//gradgradBcheck[row][col][dep] = 0.0;
			}
		}
	}
	/* 
		Loop over coils and perform Biot-Savart integral
	*/
	for (coil_index=0;coil_index<n_params;coil_index++) {
		//printf("coil_index=%d/%d\n", coil_index, n_params);
		lencoil = 0.0;
		for (coilseg_index=0;coilseg_index<num_segs[coil_index];coilseg_index++) {
			//printf("coilseg_index=%d/%d\n", coilseg_index, num_segs[coil_index]);
			current = coils[coil_index][coilseg_index][3];

			//printf("current = %f\n", current);
			////printf("num_segs[%d] = %d\n", coil_index, num_segs[coil_index]);
			////printf("coilseg_index = %d\n", coilseg_index);
			//printf("coilsegment = %f,", coils[coil_index][coilseg_index][0]);
			//printf(              "%f,", coils[coil_index][coilseg_index][1]);
			//printf(              "%f\n", coils[coil_index][coilseg_index][2]);
			//printf("current = %10.10f\n", current);
			//Xminxc = XX - 0.5*coils[coil_index][coilseg_index][0] - 0.5*coils[coil_index][(coilseg_index+1)%num_segs[coil_index]][0];
			//Yminyc = YY - 0.5*coils[coil_index][coilseg_index][1] - 0.5*coils[coil_index][(coilseg_index+1)%num_segs[coil_index]][1]; 
			//Zminzc = ZZ - 0.5*coils[coil_index][coilseg_index][2] - 0.5*coils[coil_index][(coilseg_index+1)%num_segs[coil_index]][2]; 
			//Rminrc = sqrt(pow(Xminxc,  2.0) + pow(Yminyc, 2.0) + pow(Zminzc, 2.0)); 
			//dxc = coils[coil_index][(coilseg_index+1)%num_segs[coil_index]][0] - coils[coil_index][coilseg_index][0];
			//dyc = coils[coil_index][(coilseg_index+1)%num_segs[coil_index]][1] - coils[coil_index][coilseg_index][1];
			//dzc = coils[coil_index][(coilseg_index+1)%num_segs[coil_index]][2] - coils[coil_index][coilseg_index][2];
			//drc = sqrt(pow(dxc, 2.0) + pow(dyc, 2.0) + pow(dzc, 2.0));
			////if (coil_index == 0)  printf("drc = (%f, %f, %f)\n", dxc, dyc, dzc);
			////if (coil_index == 0) {
			////	printf("coil_index = %d/%d\n", coil_index, n_params);
			////	printf("coilseg_index = %d/%d, dxc = %10.10f\n", coilseg_index, num_segs[coil_index], dxc); 
			////	printf("Xcoil = %f\n", coils[coil_index][coilseg_index][0]);
			////}
			///* project into cartesian axes aligned with local poloidal plane 
			//Rhat dot (phihat cross Zhat) = 1 prescription -> (R,phi,Z) right-handed */
			//Rminrcvect[0] = cos(varphi)*Xminxc + sin(varphi)*Yminyc;
			//Rminrcvect[1] = Zminzc;
			//Rminrcvect[2] = - sin(varphi)*Xminxc + cos(varphi)*Yminyc;
			//drcvect[0] = cos(varphi)*dxc + sin(varphi)*dyc;
			//drcvect[1] = dzc;
			//drcvect[2] = - sin(varphi)*dxc + cos(varphi)*dyc; 


			/////////////////////////////////////////////
			/////////////////////////////////////////////
			/////////////////////////////////////////////

			
			Xminxc = XX - coils[coil_index][coilseg_index][0];
			Yminyc = YY - coils[coil_index][coilseg_index][1];
			Zminzc = ZZ - coils[coil_index][coilseg_index][2];
			Rminrc = sqrt(pow(Xminxc,  2.0) + pow(Yminyc, 2.0) + pow(Zminzc, 2.0)); 
			Rminrcvect[0] = cos(varphi)*Xminxc + sin(varphi)*Yminyc;
			Rminrcvect[1] = Zminzc;
			Rminrcvect[2] = - sin(varphi)*Xminxc + cos(varphi)*Yminyc;

			Xminxcp = XX - coils[coil_index][(coilseg_index+1)%num_segs[coil_index]][0];
			Yminycp = YY - coils[coil_index][(coilseg_index+1)%num_segs[coil_index]][1];
			Zminzcp = ZZ - coils[coil_index][(coilseg_index+1)%num_segs[coil_index]][2];
			Rminrcp = sqrt(pow(Xminxcp,  2.0) + pow(Yminycp, 2.0) + pow(Zminzcp, 2.0)); 
			Rminrcvectp[0] = cos(varphi)*Xminxcp + sin(varphi)*Yminycp;
			Rminrcvectp[1] = Zminzcp;
			Rminrcvectp[2] = - sin(varphi)*Xminxcp + cos(varphi)*Yminycp;


			dxc = coils[coil_index][(coilseg_index+1)%num_segs[coil_index]][0] - coils[coil_index][coilseg_index][0];
			dyc = coils[coil_index][(coilseg_index+1)%num_segs[coil_index]][1] - coils[coil_index][coilseg_index][1];
			dzc = coils[coil_index][(coilseg_index+1)%num_segs[coil_index]][2] - coils[coil_index][coilseg_index][2];
			drc = sqrt(pow(dxc, 2.0) + pow(dyc, 2.0) + pow(dzc, 2.0));
			lencoil +=drc;
			//printf("drc = %f\n", drc);
			drcvect[0] = cos(varphi)*dxc + sin(varphi)*dyc;
			drcvect[1] = dzc;
			drcvect[2] = - sin(varphi)*dxc + cos(varphi)*dyc; 

			if (coilseg_index == 0) {
				tancoilx[0] = coils[coil_index][1][0] - coils[coil_index][num_segs[coil_index] - 1][0];
				tancoilx[1] = coils[coil_index][1][1] - coils[coil_index][num_segs[coil_index] - 1][1];
				tancoilx[2] = coils[coil_index][1][2] - coils[coil_index][num_segs[coil_index] - 1][2];
			}
			else if (coilseg_index == num_segs[coil_index]-1) {
				tancoilx[0] = coils[coil_index][0][0] - coils[coil_index][num_segs[coil_index] - 2][0];
				tancoilx[1] = coils[coil_index][0][1] - coils[coil_index][num_segs[coil_index] - 2][1];
				tancoilx[2] = coils[coil_index][0][2] - coils[coil_index][num_segs[coil_index] - 2][2];
			}
			else {
				tancoilx[0] = coils[coil_index][coilseg_index+1][0] - coils[coil_index][coilseg_index-1][0];
				tancoilx[1] = coils[coil_index][coilseg_index+1][1] - coils[coil_index][coilseg_index-1][1];
				tancoilx[2] = coils[coil_index][coilseg_index+1][2] - coils[coil_index][coilseg_index-1][2];
			}
			tancoil[0] = cos(varphi)*tancoilx[0] + sin(varphi)*tancoilx[1];
			tancoil[1] = tancoilx[2];
			tancoil[2] = - sin(varphi)*tancoilx[0] + cos(varphi)*tancoilx[1]; 
			tancoil_length = sqrt(tancoil[0]*tancoil[0] + tancoil[1]*tancoil[1] + tancoil[2]*tancoil[2]);	
			//lencoil += 0.5*tancoil_length*
			//printf("x drc vs tancoil %f %f\n", drcvect[0], tancoil[0]);
			//printf("z drc vs tancoil %f %f\n", drcvect[1], tancoil[1]);
			//printf("y drc vs tancoil %f %f\n", drcvect[2], tancoil[2]);
			//tancoil[0] /= tancoil_length;
			//tancoil[1] /= tancoil_length;
			//tancoil[2] /= tancoil_length;


			////////////////
			for (row=0;row<3;row++) {
				magfield.value[row] +=  (0.5)*pow(10.0,-7.0) * ( current / pow(Rminrc, 3.0) ) * ( Rhat[row] * ( -tancoil[1]*Rminrcvect[2] + tancoil[2]*Rminrcvect[1] ) + Zhat[row] * (-tancoil[2]*Rminrcvect[0] + tancoil[0]*Rminrcvect[2] ) + phihat[row] * (-tancoil[0]*Rminrcvect[1] + tancoil[1]*Rminrcvect[0] ) ) ;
				// different method below (below works best with shape gradient)
				//magfield.value[row] +=  0.5*pow(10.0,-7.0) * ( current / pow(Rminrc, 3.0) ) * ( Rhat[row] * ( -drcvect[1]*Rminrcvect[2] + drcvect[2]*Rminrcvect[1] ) + Zhat[row] * (-drcvect[2]*Rminrcvect[0] + drcvect[0]*Rminrcvect[2] ) + phihat[row] * (-drcvect[0]*Rminrcvect[1] + drcvect[1]*Rminrcvect[0] ) );
				//magfield.value[row] +=  0.5*pow(10.0,-7.0) * ( current / pow(Rminrcp, 3.0) ) * ( Rhat[row] * ( -drcvect[1]*Rminrcvectp[2] + drcvect[2]*Rminrcvectp[1] ) + Zhat[row] * (-drcvect[2]*Rminrcvectp[0] + drcvect[0]*Rminrcvectp[2] ) + phihat[row] * (-drcvect[0]*Rminrcvectp[1] + drcvect[1]*Rminrcvectp[0] ) );
				//printf("%10.10f\n",  -drcvect[1]*Rminrcvect[2] + drcvect[2]*Rminrcvect[1] );
				for (col=0;col<2;col++) {
					magfield.derivative[row][col] += (0.5)*( pow(10.0,-7.0) * ( current / pow(Rminrc, 3.0) ) * ( - 3.0 * Rminrcvect[col] * ( Rhat[row] * (-tancoil[1]*Rminrcvect[2] + tancoil[2]*Rminrcvect[1] ) + Zhat[row] * (-tancoil[2]*Rminrcvect[0] + tancoil[0]*Rminrcvect[2] ) + phihat[row] * (-tancoil[0]*Rminrcvect[1] + tancoil[1]*Rminrcvect[0] ) ) / pow(Rminrc, 2.0) + ( tancoil[2]*( Rhat[row]*Zhat[col] - Rhat[col]*Zhat[row])  + tancoil[1] * phihat[row] * Rhat[col] - tancoil[0] * phihat[row] *Zhat[col] ) ) ) ; 
					// different method below
					//magfield.derivative[row][col] += 0.5*( pow(10.0,-7.0) * ( current / pow(Rminrc, 3.0) ) * ( - 3.0 * Rminrcvect[col] * ( Rhat[row] * (-drcvect[1]*Rminrcvect[2] + drcvect[2]*Rminrcvect[1] ) + Zhat[row] * (-drcvect[2]*Rminrcvect[0] + drcvect[0]*Rminrcvect[2] ) + phihat[row] * (-drcvect[0]*Rminrcvect[1] + drcvect[1]*Rminrcvect[0] ) ) / pow(Rminrc, 2.0) + ( drcvect[2]*( Rhat[row]*Zhat[col] - Rhat[col]*Zhat[row])  + drcvect[1] * phihat[row] * Rhat[col] - drcvect[0] * phihat[row] *Zhat[col] ) ) ) ; 
					//magfield.derivative[row][col] += 0.5*( pow(10.0,-7.0) * ( current / pow(Rminrcp, 3.0) ) * ( - 3.0 * Rminrcvectp[col] * ( Rhat[row] * (-drcvect[1]*Rminrcvectp[2] + drcvect[2]*Rminrcvectp[1] ) + Zhat[row] * (-drcvect[2]*Rminrcvectp[0] + drcvect[0]*Rminrcvectp[2] ) + phihat[row] * (-drcvect[0]*Rminrcvectp[1] + drcvect[1]*Rminrcvectp[0] ) ) / pow(Rminrcp, 2.0) + ( drcvect[2]*( Rhat[row]*Zhat[col] - Rhat[col]*Zhat[row])  + drcvect[1] * phihat[row] * Rhat[col] - drcvect[0] * phihat[row] *Zhat[col] ) ) ) ; 
					for (dep=0; dep<2; dep++) {
						object = 3.0 * ( - Rhat[dep] *Rhat[col] - Zhat[dep]* Zhat[col] + 5.0 * Rminrcvect[dep] * Rminrcvect[col] / pow(Rminrc, 2.0) ) * ( Rhat[row] * ( - drcvect[1]*Rminrcvect[2] + drcvect[2]*Rminrcvect[1] ) + Zhat[row] * (-drcvect[2]*Rminrcvect[0] + drcvect[0]*Rminrcvect[2] ) + phihat[row] * (-drcvect[0]*Rminrcvect[1] + drcvect[1]*Rminrcvect[0] ) );
						object += 3.0 * Rminrcvect[dep]* ( drcvect[2]*( Rhat[col]*Zhat[row] - Rhat[row]*Zhat[col])  - drcvect[1] * phihat[row] * Rhat[col] + drcvect[0] * phihat[row] *Zhat[col] ) ; 
						object += 3.0 * Rhat[dep] * Rminrcvect[col] * (   drcvect[2] * Zhat[row] - drcvect[1] * phihat[row] );
						object += 3.0 * Zhat[dep] * Rminrcvect[col] * ( - drcvect[2]* Rhat[row] + drcvect[0] * phihat[row]  );
						//gradgradB[row][col][dep] += ( pow(10.0,-7.0) * ( current / pow(Rminrc, 5.0) ) * object ) ;
						magfield.twoderivative[row][col][dep] += ( pow(10.0,-7.0) * ( current / pow(Rminrc, 5.0) ) * object ) ;
						// different method below
						//object = 3.0 * ( - Rhat[dep] *Rhat[col] - Zhat[dep]* Zhat[col] + 5.0 * Rminrcvect[dep] * Rminrcvect[col] / pow(Rminrc, 2.0) ) * ( Rhat[row] * ( - drcvect[1]*Rminrcvect[2] + drcvect[2]*Rminrcvect[1] ) + Zhat[row] * (-drcvect[2]*Rminrcvect[0] + drcvect[0]*Rminrcvect[2] ) + phihat[row] * (-drcvect[0]*Rminrcvect[1] + drcvect[1]*Rminrcvect[0] ) );
						//object += 3.0 * Rminrcvect[dep]* ( drcvect[2]*( Rhat[col]*Zhat[row] - Rhat[row]*Zhat[col])  - drcvect[1] * phihat[row] * Rhat[col] + drcvect[0] * phihat[row] *Zhat[col] ) ; 
						//object += 3.0 * Rhat[dep] * Rminrcvect[col] * (   drcvect[2] * Zhat[row] - drcvect[1] * phihat[row] );
						//object += 3.0 * Zhat[dep] * Rminrcvect[col] * ( - drcvect[2]* Rhat[row] + drcvect[0] * phihat[row]  );
						////gradgradB[row][col][dep] += ( pow(10.0,-7.0) * ( current / pow(Rminrc, 5.0) ) * object ) ;
						//magfield.twoderivative[row][col][dep] += 0.5* ( pow(10.0,-7.0) * ( current / pow(Rminrc, 5.0) ) * object ) ;
						//object = 3.0 * ( - Rhat[dep] *Rhat[col] - Zhat[dep]* Zhat[col] + 5.0 * Rminrcvectp[dep] * Rminrcvectp[col] / pow(Rminrcp, 2.0) ) * ( Rhat[row] * ( - drcvect[1]*Rminrcvectp[2] + drcvect[2]*Rminrcvectp[1] ) + Zhat[row] * (-drcvect[2]*Rminrcvectp[0] + drcvect[0]*Rminrcvectp[2] ) + phihat[row] * (-drcvect[0]*Rminrcvectp[1] + drcvect[1]*Rminrcvectp[0] ) );
						//object += 3.0 * Rminrcvectp[dep]* ( drcvect[2]*( Rhat[col]*Zhat[row] - Rhat[row]*Zhat[col])  - drcvect[1] * phihat[row] * Rhat[col] + drcvect[0] * phihat[row] *Zhat[col] ) ; 
						//object += 3.0 * Rhat[dep] * Rminrcvectp[col] * (   drcvect[2] * Zhat[row] - drcvect[1] * phihat[row] );
						//object += 3.0 * Zhat[dep] * Rminrcvectp[col] * ( - drcvect[2]* Rhat[row] + drcvect[0] * phihat[row]  );
						////gradgradB[row][col][dep] += ( pow(10.0,-7.0) * ( current / pow(Rminrc, 5.0) ) * object ) ;
						//magfield.twoderivative[row][col][dep] += 0.5* ( pow(10.0,-7.0) * ( current / pow(Rminrcp, 5.0) ) * object ) ;
					}
				}
			}
		}
	//if (coil_index == 0) printf("lencoil = %f\n", lencoil);
	}
	return magfield;
}

double delta(int wantszero) {
	double answer;
	if (wantszero == 0) answer = 1.0;
	else answer = 0.0;
	return answer;
}


/* 
 evaluates magnetic field, its gradient (first derivative with respect to position) and Hessian (second derivative...)
*/
struct field Bfield(double *Xp, double varphi, struct fieldparams allparams) {
	double divB;
	double helicoil[4], helitangent[3], helitangentx[3];
	double **AA, **BB, thetacoil, phicoil, Rmaj, rmin, deltaphi, dthetadphi;
	int ind, l0=2, n_diff;
	struct field magfield, magfielddelta[3];
	struct field dmagfield[3];
	int coil_index=0, coilseg_index=0; 
	int check = 0, row, col, dep;
	int N_gridphi_toroidalturn;
	int m0_symmetry, n_coils, *intparams;
	double delta = 0.0001;
	double current;
	double gradB[3][2], gradgradB[3][2][2], gradBcheck[3][2], gradgradBcheck[3][2][2]; 
	double *epsilon, iota[2]; 
	//int amp_Domm[1], pol_Domm[1], tor_Domm[1];
	double Xminxc, Zminzc, Yminyc, XX=Xp[0]*cos(varphi), YY=Xp[0]*sin(varphi), ZZ=Xp[1], Rminrcvect[3], Rminrc, Zhat[3]={0.0, 1.0, 0.0}, Rhat[3]={1.0, 0.0, 0.0}, phihat[3] = {0.0, 0.0, 1.0};
	double object;
	double ***coils;
	m0_symmetry = allparams.m0_fieldperiods;
	//printf("m0_symmetry=%d\n", m0_symmetry);
	//printf("type %s\n", allparams.type);
	coils = allparams.diffparams;
	n_coils = allparams.n_coils;
	n_diff = allparams.n_diff;
	intparams = allparams.intparams;
	for (row=0;row<3;row++) {
		for (col=0;col<2;col++) {
			gradB[row][col] = 0.0;
			gradBcheck[row][col] = 0.0;
			for (dep=0;dep<2;dep++) {
				gradgradB[row][col][dep] = 0.0;
				gradgradBcheck[row][col][dep] = 0.0;
			}
		}
	}
	//printf("Sanity check: print line of coil file\n");
	//printf("coils[0][0][0]=%Le\n", coils[0][0][0]); printf("coils[0][0][1]=%Le\n", coils[0][0][1]);
	//printf("coils[0][0][2]=%Le\n", coils[0][0][2]); printf("coils[0][0][3]=%Le\n", coils[0][0][3]);
	if (strncmp(allparams.type, "coil", 4) == 0) {
		magfield = Bcoil(Xp[0], Xp[1], varphi, coils, n_coils, intparams); 
		if (check == 1) {
			magfielddelta[0] = Bcoil(Xp[0]+delta, Xp[1], varphi, coils, n_coils, intparams);
			magfielddelta[1] = Bcoil(Xp[0], Xp[1]+delta, varphi, coils, n_coils, intparams);
			magfielddelta[2] = Bcoil(Xp[0], Xp[1], varphi + delta, coils, n_coils, intparams);
			dmagfield[0] = addstructsfield(1/delta, magfielddelta, -1/delta, &magfield);
			dmagfield[1] = addstructsfield(1/delta, magfielddelta, -1/delta, &magfield);
			dmagfield[2] = addstructsfield(1/delta, magfielddelta, -1/delta, &magfield);
			printf("dBRdR   %f %f\n", dmagfield[0].value[0], magfield.derivative[0][0]);
			printf("dBZdR   %f %f\n", dmagfield[0].value[1], magfield.derivative[1][0]);
			printf("dBphidR %f %f\n", dmagfield[0].value[2], magfield.derivative[2][0]);
			printf("dBRdZ   %f %f\n", dmagfield[1].value[0], magfield.derivative[0][1]);
			printf("dBZdZ   %f %f\n", dmagfield[1].value[1], magfield.derivative[1][1]);
			printf("dBphidZ %f %f\n", dmagfield[1].value[2], magfield.derivative[2][1]);
			divB = ( dmagfield[0].value[0] + magfield.value[0] / Xp[0] + dmagfield[1].value[1] + dmagfield[2].value[2] / Xp[0] ) / sqrt( pow( dmagfield[0].value[0] + magfield.value[0] / Xp[0] , 2.0 ) + pow( dmagfield[1].value[1], 2.0 ) + pow( dmagfield[2].value[2] / Xp[0] , 2.0 ) );
			printf("div(B) = %f\n", divB);
		}
	}
	else if (strncmp(allparams.type, "heli", 4) == 0) {
		for (row=0;row<3;row++) {
			magfield.value[row] = 0.0; 
			for (col=0;col<2;col++) {
				magfield.derivative[row][col] = 0.0;
				for (dep=0;dep<2;dep++) {
					magfield.twoderivative[row][col][dep] = 0.0;
				}
			}
		}
		l0 = allparams.intparams[0];
		N_gridphi_toroidalturn = allparams.intparams[1];
		deltaphi = 2.0*M_PI/N_gridphi_toroidalturn;
		AA = malloc(n_coils*sizeof(double));
		BB = malloc(n_coils*sizeof(double));
		Rmaj = allparams.constparams[0];
		rmin = allparams.constparams[1];
		//printf("Rmaj = %f, rmin = %f\n", Rmaj, rmin);
		for (coil_index=0;coil_index<n_coils;coil_index++)
		{
			//printf("coil_index=%d/%d\n", coil_index, n_coil);
			AA[coil_index] = malloc((n_diff/2)*sizeof(double));
			BB[coil_index] = malloc((n_diff/2)*sizeof(double));
			AA[coil_index][0] = coils[coil_index][0][0];
			BB[coil_index][0] = 0.0;
			for (row=1; row<n_diff-1; row++) {
				if (row%2 == 1)  AA[coil_index][(row+1)/2] = coils[coil_index][0][row];
				else BB[coil_index][row/2] = coils[coil_index][0][row];
			}
			current = coils[coil_index][0][n_diff-1];
			for (coilseg_index=0;coilseg_index<l0*N_gridphi_toroidalturn;coilseg_index++)
			{
				//printf("coilseg_index = %d/%d\n", coilseg_index, l0*N_gridphi_toroidalturn);
				phicoil = coilseg_index*deltaphi;
				thetacoil = m0_symmetry*phicoil/l0;
				dthetadphi = (double) m0_symmetry/l0;
				for (ind=0;ind<n_diff/2;ind++) {
					thetacoil  += (AA[coil_index][ind]*cos(ind*m0_symmetry*phicoil/l0) + BB[coil_index][ind]*sin(ind*m0_symmetry*phicoil/l0));
					dthetadphi += ( ind*m0_symmetry*(-AA[coil_index][ind]*sin(ind*m0_symmetry*phicoil/l0) + BB[coil_index][ind]*cos(ind*m0_symmetry*phicoil/l0))/l0 );
				}
				helicoil[0] = ( Rmaj + rmin*cos(thetacoil) )* cos(phicoil);
				helicoil[1] = ( Rmaj + rmin*cos(thetacoil) )* sin(phicoil);
				helicoil[2] =   -rmin*sin(thetacoil) ;
				helicoil[3] =   current; 
				helitangentx[0] = - ( Rmaj + rmin*cos(thetacoil) )* sin(phicoil) - rmin*sin(thetacoil)*cos(phicoil)*dthetadphi;
				helitangentx[1] = ( Rmaj + rmin*cos(thetacoil) )* cos(phicoil) - rmin*sin(thetacoil)*sin(phicoil)*dthetadphi;
				helitangentx[2] =   -rmin*cos(thetacoil)*dthetadphi ;
				//printf("helicoil = (%f %f %f) with current = %f\n", helicoil[coil_index][coilseg_index][0], helicoil[coil_index][coilseg_index][1], helicoil[coil_index][coilseg_index][2], helicoil[coil_index][coilseg_index][3] );
				//printf("helitangentx = (%f %f %f)\n", helitangentx[0], helitangentx[1], helitangentx[2]);
				Xminxc = XX - helicoil[0];
				Yminyc = YY - helicoil[1]; 
				Zminzc = ZZ - helicoil[2]; 
				Rminrc = sqrt(Xminxc*Xminxc + Yminyc*Yminyc + Zminzc*Zminzc); 

				/* project into cartesian axes aligned with local poloidal plane 
				Rhat dot (phihat cross Zhat) = 1 prescription -> (R,phi,Z) right-handed */
				Rminrcvect[0] = cos(varphi)*Xminxc + sin(varphi)*Yminyc;
				Rminrcvect[1] = Zminzc;
				Rminrcvect[2] = - sin(varphi)*Xminxc + cos(varphi)*Yminyc;

				helitangent[0] = cos(varphi)*helitangentx[0] + sin(varphi)*helitangentx[1];
				helitangent[1] = helitangentx[2];
				helitangent[2] = - sin(varphi)*helitangentx[0] + cos(varphi)*helitangentx[1];
				//printf("helitangent = (%f %f %f)\n", helitangent[0], helitangent[1], helitangent[2]);
				//printf("dthetadphi = %f\n", dthetadphi); //printf("phicoil = %f; thetacoil = %f\n", phicoil, thetacoil);
				for (row=0;row<3;row++) {
					magfield.value[row] +=   deltaphi*( current / pow(Rminrc, 3.0) ) * ( Rhat[row] * ( -helitangent[1]*Rminrcvect[2] + helitangent[2]*Rminrcvect[1] ) + Zhat[row] * (-helitangent[2]*Rminrcvect[0] + helitangent[0]*Rminrcvect[2] ) + phihat[row] * (-helitangent[0]*Rminrcvect[1] + helitangent[1]*Rminrcvect[0] ) );
					for (col=0;col<2;col++) {
						magfield.derivative[row][col] += deltaphi*( ( current / pow(Rminrc, 3.0) ) * ( - 3.0 * Rminrcvect[col] * ( Rhat[row] * (-helitangent[1]*Rminrcvect[2] + helitangent[2]*Rminrcvect[1] ) + Zhat[row] * (-helitangent[2]*Rminrcvect[0] + helitangent[0]*Rminrcvect[2] ) + phihat[row] * (-helitangent[0]*Rminrcvect[1] + helitangent[1]*Rminrcvect[0] ) ) / pow(Rminrc, 2.0) + ( helitangent[2]*( Rhat[row]*Zhat[col] - Rhat[col]*Zhat[row])  + helitangent[1] * phihat[row] * Rhat[col] - helitangent[0] * phihat[row] *Zhat[col] ) ) ) ; 
						for (dep=0; dep<2; dep++) {
							object = 3.0 * ( - Rhat[dep] *Rhat[col] - Zhat[dep]* Zhat[col] + 5.0 * Rminrcvect[dep] * Rminrcvect[col] / pow(Rminrc, 2.0) ) * ( Rhat[row] * ( - helitangent[1]*Rminrcvect[2] + helitangent[2]*Rminrcvect[1] ) + Zhat[row] * (-helitangent[2]*Rminrcvect[0] + helitangent[0]*Rminrcvect[2] ) + phihat[row] * (-helitangent[0]*Rminrcvect[1] + helitangent[1]*Rminrcvect[0] ) );
							object += 3.0 * Rminrcvect[dep]* ( helitangent[2]*( Rhat[col]*Zhat[row] - Rhat[row]*Zhat[col])  - helitangent[1] * phihat[row] * Rhat[col] + helitangent[0] * phihat[row] *Zhat[col] ) ; 
							object += 3.0 * Rhat[dep] * Rminrcvect[col] * (   helitangent[2] * Zhat[row] - helitangent[1] * phihat[row] );
							object += 3.0 * Zhat[dep] * Rminrcvect[col] * ( - helitangent[2]* Rhat[row] + helitangent[0] * phihat[row]  );
							magfield.twoderivative[row][col][dep] += deltaphi*(  ( current / pow(Rminrc, 5.0) ) * object ) ;
						}
					}
				}
			}
		}
		magfield.value[2] += 1.0/Xp[0];
		magfield.derivative[2][0] -= 1.0/pow(Xp[0], 2.0);
		magfield.twoderivative[2][0][0] += 2.0/pow(Xp[0], 3.0);
		//printstructfield("magfield", magfield);
	}
	else if (strncmp(allparams.type, "Reim", 4) == 0) {
		iota[0] = coils[0][0][0]; // 0.15; 
		iota[1] = coils[0][0][1]; //0.38;
		epsilon = coils[0][0]+2; //0.001* 
		//printf("n_params=%d\n", n_params);
		//printf("intparams = %d\n", intparams[0]);
		//printf("m0_symmetry = %d\n", m0_symmetry);
		//printf("iota0 = %f, iota1 = %f epsilon = %f\n", iota[0], iota[1], epsilon[0]);
		magfield = BReim(Xp[0], Xp[1], varphi, m0_symmetry, iota[0], iota[1], epsilon, intparams+1, intparams[0]-2);
	}
	else {
		printf("ERROR: No identifiable magnetic field type was found\n");
		exit(1);
	}
	return magfield;
}

/* 
 the function below evaluates the derivative of the Reiman model magnetic field with respect to parameters:
 ι_ax (or ι_0) , ι'_ax ( or ι_1), and ε_n for all allowed n
*/
struct field *gradBReim(int m0_symmetry, double iota0, double iota1, double *epsilon, int *k_theta, int size_epsilon, double RR, double ZZ, double varphi) {
		int ind;
		double R_axis = 1.0, theta, rmin; // check_q;
		double dcomboiota, dcomboiota1, dcomboepsilon[size_epsilon], dcombo1epsilon[size_epsilon];
		double d_dcombodR_epsilon[size_epsilon], d_dcombodZ_epsilon[size_epsilon], d_dcombo1dR_epsilon[size_epsilon], d_dcombo1dZ_epsilon[size_epsilon];
		double d_dcombodR_iota1, d_dcombodZ_iota1;
		struct field *mag = malloc((2+size_epsilon)*sizeof(struct field));
		//printf("k_theta[0] = %d\n", k_theta[0]);
		//printf("m0 = %d\n", m0_symmetry);
		//printf("iota0 = %f\n", iota0);
		//printf("iota1 = %f\n", iota1);
		//printf("epsilon = %f\n", epsilon[0]);
		theta = atan2(ZZ, RR - R_axis);
		rmin = sqrt(pow((RR-R_axis), 2.0) + pow(ZZ, 2.0));
		//combo = iota0 + iota1*rmin*rmin - 2.0*epsilon[0]*cos(2.0*theta - varphi) - 3.0*epsilon[1]*rmin*cos(3.0*theta - varphi);
		//combo1 = 2.0*epsilon[0]*rmin*rmin*sin(2.0*theta - varphi) + 3.0*epsilon[1]*pow(rmin, 3.0)*sin(3.0*theta - varphi);

		dcomboiota = 1.0;
		dcomboiota1 = rmin*rmin;
		mag[0].value[0] = (ZZ / RR) * dcomboiota; 
		mag[0].value[1] = - ( (RR - R_axis) / RR ) * dcomboiota; 
		mag[0].value[2] = 0.0;
		mag[1].value[0] = (ZZ / RR) * dcomboiota1; 
		mag[1].value[1] = - ( (RR - R_axis) / RR ) * dcomboiota1; 
		mag[1].value[2] = 0.0;
		d_dcombodR_iota1 = 2.0*(RR - R_axis);
		d_dcombodZ_iota1 = 2.0*ZZ;
		mag[0].derivative[0][0] = ( - ZZ / pow( RR, 2.0) ) *  dcomboiota ;
		mag[0].derivative[1][0] = - ( R_axis / pow( RR, 2.0) ) *  dcomboiota;
		mag[0].derivative[2][0] = 0.0;
		mag[0].derivative[0][1] = ( 1.0 / RR ) * dcomboiota;
		mag[0].derivative[1][1] = 0.0;
		mag[0].derivative[2][1] = 0.0;
		mag[1].derivative[0][0] = ( - ZZ / pow( RR, 2.0) ) *  dcomboiota1 + (ZZ / RR) * d_dcombodR_iota1;
		mag[1].derivative[1][0] = - ( R_axis / pow( RR, 2.0) ) *  dcomboiota1 - ((RR - R_axis) / RR) * d_dcombodR_iota1;
		mag[1].derivative[2][0] = 0.0;
		mag[1].derivative[0][1] = ( 1.0 / RR ) * dcomboiota1 + (ZZ / RR) * d_dcombodZ_iota1;
		mag[1].derivative[1][1] = - ((RR - R_axis) / RR) * d_dcombodZ_iota1;
		mag[1].derivative[2][1] = 0.0;
		for (ind=0; ind < size_epsilon; ind++) {
			dcomboepsilon[ind]  = - k_theta[ind] * pow(rmin, k_theta[ind] - 2) * cos(k_theta[ind]*theta - m0_symmetry*varphi);
			dcombo1epsilon[ind] =   k_theta[ind] * pow(rmin, k_theta[ind] - 2) * sin(k_theta[ind]*theta - m0_symmetry*varphi);
			mag[ind+2].value[0] = (ZZ / RR) * dcomboepsilon[ind] + dcombo1epsilon[ind] * (RR - R_axis) / RR;
			mag[ind+2].value[1] = - ( (RR - R_axis) / RR ) * dcomboepsilon[ind] + dcombo1epsilon[ind] * ZZ / RR;
			mag[ind+2].value[2] = 0.0;
			d_dcombodR_epsilon[ind] = - k_theta[ind] * pow( rmin, k_theta[ind] - 4 ) * ( k_theta[ind] * sin( k_theta[ind]*theta - m0_symmetry*varphi ) * ZZ + (k_theta[ind] - 2) *  (RR - R_axis) * cos(k_theta[ind]*theta - m0_symmetry * varphi) );
			d_dcombodZ_epsilon[ind] = - pow( rmin, k_theta[ind] - 4 ) * k_theta[ind] * (  - k_theta[ind] * sin( k_theta[ind]*theta - m0_symmetry*varphi ) * ( RR - R_axis ) + (k_theta[ind] - 2) * ZZ * cos(k_theta[ind]*theta - m0_symmetry * varphi) );
			d_dcombo1dR_epsilon[ind] = k_theta[ind] * pow(rmin, k_theta[ind] - 4) * ( - k_theta[ind] * cos(k_theta[ind]*theta - m0_symmetry*varphi) * ZZ  +  (k_theta[ind] - 2) * sin(k_theta[ind]*theta - m0_symmetry*varphi) * (RR - R_axis) ) ;
			d_dcombo1dZ_epsilon[ind] = k_theta[ind] * pow(rmin, k_theta[ind] - 4) * ( k_theta[ind] * cos(k_theta[ind]*theta - m0_symmetry *varphi) * (RR - R_axis) + (k_theta[ind] - 2) * sin(k_theta[ind]*theta - m0_symmetry*varphi) * ZZ ) ;
			mag[ind+2].derivative[0][0] = ( - ZZ / pow( RR, 2.0) ) *  dcomboepsilon[ind] + (ZZ / RR) * d_dcombodR_epsilon[ind] + dcombo1epsilon[ind] * R_axis / pow(RR, 2.0) + d_dcombo1dR_epsilon[ind] * (RR - R_axis) / RR ;
			mag[ind+2].derivative[1][0] = - ( R_axis / pow( RR, 2.0) ) *  dcomboepsilon[ind] - ((RR - R_axis) / RR) * d_dcombodR_epsilon[ind] - dcombo1epsilon[ind] * ( ZZ / pow(RR, 2.0) ) + d_dcombo1dR_epsilon[ind] * ZZ / RR;
			mag[ind+2].derivative[2][0] = 0.0;
			mag[ind+2].derivative[0][1] = ( 1.0 / RR ) * dcomboepsilon[ind] + (ZZ / RR) * d_dcombodZ_epsilon[ind] + d_dcombo1dZ_epsilon[ind] * ( RR - R_axis ) / RR;
			mag[ind+2].derivative[1][1] = - ((RR - R_axis) / RR) * d_dcombodZ_epsilon[ind] + dcombo1epsilon[ind] * ( 1.0 / RR ) + d_dcombo1dZ_epsilon[ind] * ZZ / RR;
			mag[ind+2].derivative[2][1] = 0.0;
		}
		return mag;
}

/* 
 the function below evaluates the shape gradient at a coil location with respect to that coil location
*/
struct field *shapeBcoil(double *Xp, double varphi, double *coilseg) {
	int row, col, dep;
	double Rminrc, Rminrcvectp[3], Rminrcvectx[3], Rvectp[3], rvectp[3], current;
	double rotation[3][3];
	struct field *gradmagfield = malloc(3*sizeof(struct field));
		current = coilseg[3];
		Rminrc = 0.0;
		rvectp[0] = coilseg[0]*cos(varphi) + coilseg[1]*sin(varphi);
		rvectp[1] = coilseg[2];
		rvectp[2] = coilseg[1]*cos(varphi) - coilseg[0]*sin(varphi);
		Rvectp[0] = Xp[0];
		Rvectp[1] = Xp[1];
		Rvectp[2] = 0.0;  // coordinate defined in locally cartesian plane at point
		rotation[0][0] = rotation[1][2] = cos(varphi); 		
		rotation[1][0] = sin(varphi); 		
		rotation[0][2] = -sin(varphi); 		
		rotation[2][1] = 1.0;
		rotation[2][0] = rotation[0][1] = rotation[2][2] = rotation[1][1] = 0.0;
		for (row=0;row<3;row++) {
			Rminrcvectp[row] = Rvectp[row] - rvectp[row];
			Rminrc += pow(Rminrcvectp[row], 2.0);
		}
		for (row=0;row<3;row++) {
			//Rminrcvectx[row] = rotation[0][row]*Rminrcvectp[0] + rotation[1][row]*Rminrcvectp[1] + rotation[2][row]*Rminrcvectp[2];
			Rminrcvectx[row] = rotation[row][0]*Rminrcvectp[0] + rotation[row][1]*Rminrcvectp[1] + rotation[row][2]*Rminrcvectp[2]; // should be correct
		}
		Rminrc = pow(Rminrc, 0.5);
		//printf("current = %f, Rminrc = %f\n", current, Rminrc);
		for (row=0;row<3;row++) {
			//printf("Rminrcvectx[%d]=%f, Rminrcvectp[%d]=%f\n", row, Rminrcvectx[row], row, Rminrcvectp[row]);
			for (col=0; col<3; col++) {
				gradmagfield[row].value[col] 
				= pow(10.0, -7.0)*(current/pow(Rminrc, 3.0))*( - rotation[row][col] + 3.0* Rminrcvectx[row] * Rminrcvectp[col] / pow(Rminrc, 2.0) );			
				//= pow(10.0, -7.0)*(current/pow(Rminrc, 3.0))*( - rotation[row][col] + Rminrcvectx[row] * Rminrcvectp[col] / pow(Rminrc, 2.0) );			
				//printf("rotation=%f\n", rotation[row][col]);
				//printf("grad=%f\n", gradmagfield[row].value[col]);
				for (dep=0; dep<2;dep++) {
					//printf("param = (%f, %f, %f, %f)\n", param[0], param[1], param[2], param[3]);
					gradmagfield[row].derivative[col][dep] 
					= pow(10.0, -7.0)*(current/pow(Rminrc, 5.0))*( 3.0*Rminrcvectp[dep]*( rotation[row][col] - 5.0 *Rminrcvectx[row] *Rminrcvectp[col] / pow(Rminrc, 2.0) ) + 3.0 *delta(dep - col) *Rminrcvectx[row] + 3.0 * rotation[row][dep] *Rminrcvectp[col] ) ;	
				}
			}
		}
	return gradmagfield;
}

/* 
 the function below evaluates the gradient of the Cary-Hanson helical magnetic field model with respect to the Fourier mode amplitudes of the perturbations to the given helical coil and with respect to the current in the coil. One coil at a time is done: Cary & Hanson's configurations have a pair of coils so coil_index = 0 or 1
*/
struct field *gradBheli(double *Xp, double varphi, int m0_symmetry, struct fieldparams allparams, int coil_index) {
	int ind_param, row, col, dep, ind, coilseg_index; 
	int l0, N_gridphi_toroidalturn, n_segs, n_diff;
	double XX = Xp[0]*cos(varphi), YY = Xp[0]*sin(varphi), ZZ = Xp[1];
	double Xminxc, Yminyc, Zminzc, Rminrcvect[3], drcdkappa[3], crossprdct[3];
	double Rminrc, Rminrcvectp[3], Rminrcvectx[3], Rvectp[3], rvectp[3], current;
	double helicoil[3], helitanx[3], helitanp[3];
	double Rmaj, rmin;
	double rotation[3][3], *AA, *BB;
	double thetacoil, phicoil, deltaphi, dthetadphi;
	struct field *gradmagfield;
	n_diff = allparams.n_diff;
	gradmagfield = calloc(n_diff, sizeof(struct field));
	for (row=0;row<n_diff;row++) {
		for (col=0;col<3;col++) {
			gradmagfield[row].value[col] = 0.0; 
			for (dep=0;dep<2;dep++) {
				gradmagfield[row].derivative[col][dep] = 0.0;
			}
		}
	}
	l0 = allparams.intparams[0];
	//printf("l0=%d\n", l0);
	//printf("m0=%d\n", m0_symmetry);
	//printf("n_coil=%d\n", n_coil);
	//printf("n_diff=%d\n", n_diff);
	N_gridphi_toroidalturn = allparams.intparams[1];
	deltaphi = 2.0*M_PI/N_gridphi_toroidalturn;
	Rmaj = allparams.constparams[0];
	rmin = allparams.constparams[1];
	AA = malloc((n_diff/2)*sizeof(double));
	BB = malloc((n_diff/2)*sizeof(double));
	//printf("coil_index=%d\n", coil_index);
	//printf("n_diff=%d\n", n_diff);
	//printf("diffparams=%f\n", allparams.diffparams[0][0][0]);
	for (ind_param=0; ind_param< n_diff-1; ind_param++) {
		//printf("diffparams=%f\n", allparams.diffparams[0][0][0]);
		//printf("diffparams=%f\n", allparams.diffparams[1][0][0]);
		//printf("ind_param=%d/%d\n", ind_param, n_diff-1);
		AA[0] = allparams.diffparams[coil_index][0][0];
		BB[0] = 0.0;
		for (row=1; row<n_diff-1; row++) {
			if (row%2 == 1)  AA[(row+1)/2] = allparams.diffparams[coil_index][0][row];
			else             BB[(row)/2] =   allparams.diffparams[coil_index][0][row];
		}
		n_segs = N_gridphi_toroidalturn*l0; 
		//printf("%d=%d\n", 3 + 2*n_modes, n_params);
		current = allparams.diffparams[coil_index][0][n_diff-1];
		//printf("coil_index = %d\n", coil_index);
		//printf("current = %f\n", current);
		//printf("n_params = %d\n", num_coils);
		//printf("coil_index = %d\n", coil_index);
		//for (coilseg_index=0;coilseg_index<*(*num_segs+coil_index)-1;coilseg_index++)
		for (coilseg_index=0;coilseg_index<n_segs;coilseg_index++)
		{
			phicoil = coilseg_index*deltaphi;
			//thetacoil = M_PI*coil_index + m0_symmetry*phicoil/l0;
			thetacoil = m0_symmetry*phicoil/l0;
			dthetadphi = m0_symmetry/(1.0*l0);
			//printf("m0_symmetry = %d\n", m0_symmetry);
			for (ind=0;ind<n_diff/2;ind++) {
				thetacoil  += (AA[ind]*cos(ind*m0_symmetry*phicoil/l0) + BB[ind]*sin(ind*m0_symmetry*phicoil/l0));
				dthetadphi += ( ind*m0_symmetry*(-AA[ind]*sin(ind*m0_symmetry*phicoil/l0) + BB[ind]*cos(ind*m0_symmetry*phicoil/l0))/l0 );
			}
			//printf("dthetadphi = %f\n", dthetadphi);

			helicoil[0] = ( Rmaj + rmin*cos(thetacoil) )* cos(phicoil);
			helicoil[1] = ( Rmaj + rmin*cos(thetacoil) )* sin(phicoil);
			helicoil[2] =   -rmin*sin(thetacoil) ;
			helitanx[0] = - ( Rmaj + rmin*cos(thetacoil) )* sin(phicoil) - rmin*sin(thetacoil)*cos(phicoil)*dthetadphi;
			helitanx[1] = ( Rmaj + rmin*cos(thetacoil) )* cos(phicoil) - rmin*sin(thetacoil)*sin(phicoil)*dthetadphi;
			helitanx[2] =   -rmin*cos(thetacoil)*dthetadphi ;
			Xminxc = XX - helicoil[0];
			Yminyc = YY - helicoil[1]; 
			Zminzc = ZZ - helicoil[2]; 
			Rminrc = sqrt(Xminxc*Xminxc + Yminyc*Yminyc + Zminzc*Zminzc); 

			/* project into cartesian axes aligned with local poloidal plane 
			Rhat dot (phihat cross Zhat) = 1 prescription -> (R,phi,Z) right-handed */
			Rminrcvect[0] = cos(varphi)*Xminxc + sin(varphi)*Yminyc;
			Rminrcvect[1] = Zminzc;
			Rminrcvect[2] = - sin(varphi)*Xminxc + cos(varphi)*Yminyc;

			helitanp[0] = cos(varphi)*helitanx[0] + sin(varphi)*helitanx[1];
			helitanp[1] = helitanx[2];
			helitanp[2] = - sin(varphi)*helitanx[0] + cos(varphi)*helitanx[1];
			ind = (ind_param+1)/2;
			if (ind_param == 0) {
				drcdkappa[0] = -rmin*sin(thetacoil)*cos(phicoil);
				drcdkappa[1] = -rmin*sin(thetacoil)*sin(phicoil);
				drcdkappa[2] = -rmin*cos(thetacoil);
			}
			else if (ind_param%2 == 1) {
				drcdkappa[0] = -rmin*cos(ind*m0_symmetry*phicoil/l0)*sin(thetacoil)*cos(phicoil);
				drcdkappa[1] = -rmin*cos(ind*m0_symmetry*phicoil/l0)*sin(thetacoil)*sin(phicoil);
				drcdkappa[2] = -rmin*cos(ind*m0_symmetry*phicoil/l0)*cos(thetacoil);
			}
			else {
				drcdkappa[0] = -rmin*sin(ind*m0_symmetry*phicoil/l0)*sin(thetacoil)*cos(phicoil);
				drcdkappa[1] = -rmin*sin(ind*m0_symmetry*phicoil/l0)*sin(thetacoil)*sin(phicoil);
				drcdkappa[2] = -rmin*sin(ind*m0_symmetry*phicoil/l0)*cos(thetacoil);
			}

			for (row=0; row<3; row++) crossprdct[row] = drcdkappa[(row+1)%3]*helitanx[(row+2)%3] -  drcdkappa[(row+2)%3]*helitanx[(row+1)%3]; 
				//printf("drcdkappa = %f, helitanx = %f\n", drcdkappa[(row+1)%3], helitanx[(row+2)%3]);
			
			Rminrc = 0.0;
			rvectp[0] = helicoil[0]*cos(varphi) + helicoil[1]*sin(varphi);
			rvectp[1] = helicoil[2];
			rvectp[2] = helicoil[1]*cos(varphi) - helicoil[0]*sin(varphi);
			Rvectp[0] = Xp[0];
			Rvectp[1] = Xp[1];
			Rvectp[2] = 0.0;  // coordinate defined in locally cartesian plane at point
			rotation[0][0] = rotation[1][2] = cos(varphi); 		
			rotation[1][0] = sin(varphi); 		
			rotation[0][2] = -sin(varphi); 		
			rotation[2][1] = 1.0;
			rotation[2][0] = rotation[0][1] = rotation[2][2] = rotation[1][1] = 0.0;
			for (row=0;row<3;row++) {
				Rminrcvectp[row] = Rvectp[row] - rvectp[row];
				Rminrc += pow(Rminrcvectp[row], 2.0);
			}
			for (row=0;row<3;row++) {
				//Rminrcvectx[row] = rotation[0][row]*Rminrcvectp[0] + rotation[1][row]*Rminrcvectp[1] + rotation[2][row]*Rminrcvectp[2];
				Rminrcvectx[row] = rotation[row][0]*Rminrcvectp[0] + rotation[row][1]*Rminrcvectp[1] + rotation[row][2]*Rminrcvectp[2]; // should be correct
			}
			Rminrc = pow(Rminrc, 0.5);
			//printf("current = %f, Rminrc = %f\n", current, Rminrc);
			for (row=0;row<3;row++) {
				for (col=0; col<3; col++) {
					gradmagfield[ind_param].value[col] 
					+= deltaphi* ( crossprdct[row]*(current/pow(Rminrc, 3.0))*( - rotation[row][col] + 3.0* Rminrcvectx[row] * Rminrcvectp[col] / pow(Rminrc, 2.0) ) );			
					//printf("crossprdct=%f, Rminrcvectx=%f, Rminrcvectp=%f, current= %f, Rminrc=%f, rotation = %f\n", crossprdct[row], Rminrcvectx[row], Rminrcvectp[col], current, Rminrc, rotation[row][col]);
					//printf("what is being added: %14.14f\n", crossprdct[row]*(current/pow(Rminrc, 3.0))*( - rotation[row][col] + 3.0* Rminrcvectx[row] * Rminrcvectp[col] / pow(Rminrc, 2.0) ) );
					for (dep=0; dep<2;dep++) {
						//printf("%d, %d, %f\n", dep, col, delta(dep - col));
						gradmagfield[ind_param].derivative[col][dep] 
						+= deltaphi*crossprdct[row]*(current/pow(Rminrc, 5.0))*( 3.0*Rminrcvectp[dep]*( rotation[row][col] - 5.0 *Rminrcvectx[row] *Rminrcvectp[col] / pow(Rminrc, 2.0) ) + 3.0 *delta(dep - col) *Rminrcvectx[row] + 3.0 * rotation[row][dep] *Rminrcvectp[col] ) ;	
					}
				}
			}
		}
		//printf("gradB[%d] = (%f, %f, %f)\n", ind_param, gradmagfield[ind_param].value[0], gradmagfield[ind_param].value[1], gradmagfield[ind_param].value[2]);
	}
/* 
 The last one is the gradient with respect to coil current, so derivative of current with respect to the last parameter is just one. the derivative is therefore equivalent to evaluating the magnetic field with current = 1.0
*/
	current = 1.0; 
	ind_param = n_diff - 1;
	AA[0] = allparams.diffparams[coil_index][0][0];
	BB[0] = 0.0;
	for (row=1; row<n_diff-1; row++) {
		if (row%2 == 1)  AA[(row+1)/2] = allparams.diffparams[coil_index][0][row];
		else             BB[(row)/2] =   allparams.diffparams[coil_index][0][row];
	}
	n_segs = N_gridphi_toroidalturn*l0; 
	//printf("%d=%d\n", 3 + 2*n_modes, n_params); 
	//printf("coil_index = %d\n", coil_index);
	//printf("current = %f\n", current);
	//printf("n_params = %d\n", num_coils);
	//printf("coil_index = %d\n", coil_index);
	//for (coilseg_index=0;coilseg_index<*(*num_segs+coil_index)-1;coilseg_index++)
	for (coilseg_index=0;coilseg_index<n_segs;coilseg_index++)
	{
		phicoil = coilseg_index*deltaphi;
		//thetacoil = M_PI*coil_index + m0_symmetry*phicoil/l0;
		thetacoil = m0_symmetry*phicoil/l0;
		dthetadphi = m0_symmetry/(1.0*l0);
		//printf("m0_symmetry = %d\n", m0_symmetry);
		for (ind=0;ind<n_diff/2;ind++) {
			thetacoil  += (AA[ind]*cos(ind*m0_symmetry*phicoil/l0) + BB[ind]*sin(ind*m0_symmetry*phicoil/l0));
			dthetadphi += ( ind*m0_symmetry*(-AA[ind]*sin(ind*m0_symmetry*phicoil/l0) + BB[ind]*cos(ind*m0_symmetry*phicoil/l0))/l0 );
		}
		//printf("dthetadphi = %f\n", dthetadphi);

		helicoil[0] = ( Rmaj + rmin*cos(thetacoil) )* cos(phicoil);
		helicoil[1] = ( Rmaj + rmin*cos(thetacoil) )* sin(phicoil);
		helicoil[2] =   -rmin*sin(thetacoil) ;
		helitanx[0] = - ( Rmaj + rmin*cos(thetacoil) )* sin(phicoil) - rmin*sin(thetacoil)*cos(phicoil)*dthetadphi;
		helitanx[1] = ( Rmaj + rmin*cos(thetacoil) )* cos(phicoil) - rmin*sin(thetacoil)*sin(phicoil)*dthetadphi;
		helitanx[2] =   -rmin*cos(thetacoil)*dthetadphi ;
		Xminxc = XX - helicoil[0];
		Yminyc = YY - helicoil[1]; 
		Zminzc = ZZ - helicoil[2]; 
		Rminrc = sqrt(Xminxc*Xminxc + Yminyc*Yminyc + Zminzc*Zminzc); 

		/* project into cartesian axes aligned with local poloidal plane 
		Rhat dot (phihat cross Zhat) = 1 prescription -> (R,phi,Z) right-handed */
		Rminrcvect[0] = cos(varphi)*Xminxc + sin(varphi)*Yminyc;
		Rminrcvect[1] = Zminzc;
		Rminrcvect[2] = - sin(varphi)*Xminxc + cos(varphi)*Yminyc;

		helitanp[0] = cos(varphi)*helitanx[0] + sin(varphi)*helitanx[1];
		helitanp[1] = helitanx[2];
		helitanp[2] = - sin(varphi)*helitanx[0] + cos(varphi)*helitanx[1];
		ind = (ind_param+1)/2;
		if (ind_param == 0) {
			drcdkappa[0] = -rmin*sin(thetacoil)*cos(phicoil);
			drcdkappa[1] = -rmin*sin(thetacoil)*sin(phicoil);
			drcdkappa[2] = -rmin*cos(thetacoil);
		}
		else if (ind_param%2 == 1) {
			drcdkappa[0] = -rmin*cos(ind*m0_symmetry*phicoil/l0)*sin(thetacoil)*cos(phicoil);
			drcdkappa[1] = -rmin*cos(ind*m0_symmetry*phicoil/l0)*sin(thetacoil)*sin(phicoil);
			drcdkappa[2] = -rmin*cos(ind*m0_symmetry*phicoil/l0)*cos(thetacoil);
		}
		else {
			drcdkappa[0] = -rmin*sin(ind*m0_symmetry*phicoil/l0)*sin(thetacoil)*cos(phicoil);
			drcdkappa[1] = -rmin*sin(ind*m0_symmetry*phicoil/l0)*sin(thetacoil)*sin(phicoil);
			drcdkappa[2] = -rmin*sin(ind*m0_symmetry*phicoil/l0)*cos(thetacoil);
		}

		for (row=0; row<3; row++) {
			crossprdct[row] = drcdkappa[(row+1)%3]*helitanx[(row+2)%3] -  drcdkappa[(row+2)%3]*helitanx[(row+1)%3]; 
			//printf("drcdkappa = %f, helitanx = %f\n", drcdkappa[(row+1)%3], helitanx[(row+2)%3]);
		}
		
		Rminrc = 0.0;
		rvectp[0] = helicoil[0]*cos(varphi) + helicoil[1]*sin(varphi);
		rvectp[1] = helicoil[2];
		rvectp[2] = helicoil[1]*cos(varphi) - helicoil[0]*sin(varphi);
		Rvectp[0] = Xp[0];
		Rvectp[1] = Xp[1];
		Rvectp[2] = 0.0;  // coordinate defined in locally cartesian plane at point
		rotation[0][0] = rotation[1][2] = cos(varphi); 		
		rotation[1][0] = sin(varphi); 		
		rotation[0][2] = -sin(varphi); 		
		rotation[2][1] = 1.0;
		rotation[2][0] = rotation[0][1] = rotation[2][2] = rotation[1][1] = 0.0;
		for (row=0;row<3;row++) {
			Rminrcvectp[row] = Rvectp[row] - rvectp[row];
			Rminrc += pow(Rminrcvectp[row], 2.0);
		}
		for (row=0;row<3;row++) {
			//Rminrcvectx[row] = rotation[0][row]*Rminrcvectp[0] + rotation[1][row]*Rminrcvectp[1] + rotation[2][row]*Rminrcvectp[2];
			Rminrcvectx[row] = rotation[row][0]*Rminrcvectp[0] + rotation[row][1]*Rminrcvectp[1] + rotation[row][2]*Rminrcvectp[2]; // should be correct
		}
		Rminrc = pow(Rminrc, 0.5);
		//printf("current = %f, Rminrc = %f\n", current, Rminrc);
		for (row=0;row<3;row++) {
			gradmagfield[ind_param].value[row] += 
			( current* deltaphi*( 1.0 / pow(Rminrc, 3.0) ) * ( Rminrcvect[(row+1)%3]*helitanp[(row+2)%3] - Rminrcvect[(row+2)%3]*helitanp[(row+1)%3] ) );
			for (col=0;col<2;col++) {
				gradmagfield[ind_param].derivative[row][col] += 
				( current* deltaphi* ( ( current / pow(Rminrc, 3.0) ) * ( - 3.0 * Rminrcvect[col] * ( -helitanp[(row+1)%3]*Rminrcvect[(row+2)%3] + helitanp[(row+2)%3]*Rminrcvect[(row+1)%3] ) ) / pow(Rminrc, 2.0) + helitanp[(row+1)%3]*delta((row+2)%3 - col) ) ) ;
//( helitanp[2]*( Rhat[row]*Zhat[col] - Rhat[col]*Zhat[row])  + helitanp[1] * phihat[row] * Rhat[col] - helitanp[0] * phihat[row] *Zhat[col] ) ) ) ; 
			}
		}
	}
	//printf("gradB[%d] = (%f, %f, %f)\n", ind_param, gradmagfield[ind_param].value[0], gradmagfield[ind_param].value[1], gradmagfield[ind_param].value[2]);
	free(AA); free(BB);
	return gradmagfield;
}	
/* 
 the function below evaluates the derivative of the magnetic field with respect to parameters
 it just redirects to other functions that deal with each magnetic field type separately
*/
struct field *gradBfield(double *Xp, double varphi, struct fieldparams allparams, int diffparam_ind1, int diffparam_ind2) {
	int m0 = allparams.m0_fieldperiods;
	double *epsilon, iota[2];
	struct field *gradmagfield = NULL;
	if (strncmp(allparams.type, "Reim", 4) == 0) {
		iota[0] = allparams.diffparams[0][0][0]; 
		iota[1] = allparams.diffparams[0][0][1];
		epsilon = allparams.diffparams[0][0]+2; 
		//printf("iota0 = %f, iota1 = %f epsilon = %f\n", iota[0], iota[1], epsilon[0]);
		gradmagfield = gradBReim(m0, iota[0], iota[1], epsilon, allparams.intparams+1, allparams.intparams[0]-2, Xp[0], Xp[1], varphi);
	}
	else if (strncmp(allparams.type, "coil", 4) == 0) {
		gradmagfield = shapeBcoil(Xp, varphi, allparams.diffparams[diffparam_ind1][diffparam_ind2]);
	}
	else if (strncmp(allparams.type, "heli", 4) == 0) {
		gradmagfield = gradBheli(Xp, varphi, m0, allparams, diffparam_ind1);
	}
	return gradmagfield;
}
