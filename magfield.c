//Author: Alessandro Geraldini


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "isc.h"
#include <string.h>

double ***coil_grid(int *num_coils, int **num_segs) {
	clock_t start = clock();
	int lenline = 1000, line_num, ii=0;
	double *storevals;
	char line[lenline];
	int nn, max_num_segs, num_coils_temp=0, coil_index=0, coilseg_index=0, *num_segs_temp;
	double ***XXvect_coil= NULL;
	XXvect_coil = calloc(4,sizeof(double));
	FILE *file_coils;

	clock_t int1 = clock();
	printf("Time after declarations: %f\n", (double) (int1-start)/CLOCKS_PER_SEC);
	file_coils = fopen("coils.c09r00", "r");
	if (file_coils == NULL)
	{	
		printf("Cannot open coils.c09r00\n");
		exit(1);
	}
	line_num = 0;
	while (fgets(line, lenline, file_coils) != NULL)
	{	
		line_num += 1;
		storevals = linetodata(line, &nn);
		if (nn > 4)
		{
			coil_index += 1;
		}
	}
	rewind(file_coils);
	num_coils_temp = coil_index;
	num_segs_temp = calloc(coil_index,sizeof(int));
	*num_segs = calloc(coil_index,sizeof(int));
	max_num_segs = 0;
	coil_index = 0;
	while (fgets(line, lenline, file_coils) != NULL)
	{	
		storevals = linetodata(line, &nn);
		if (nn == 4)
		{
			coilseg_index += 1; 
		}
		else if (nn > 4)
		{
			num_segs_temp[coil_index] = coilseg_index;
			coil_index += 1;
			if (coilseg_index > max_num_segs)
			{
				max_num_segs = coilseg_index;
			}
			coilseg_index = 0;
		}
	}
	printf("num_coils_temp = %d\n", num_coils_temp);
	coil_index = 0;
	coilseg_index = 0;
	XXvect_coil = calloc(num_coils_temp,sizeof(double));
	for (ii=0;ii<4;ii+=1)
	{
		XXvect_coil[ii] = calloc(num_coils_temp,sizeof(double));
		for (coil_index = 0; coil_index < num_coils_temp; coil_index ++)
		{
			XXvect_coil[ii][coil_index] = calloc(max_num_segs,sizeof(double));
		}
	}
	coil_index = 0;
	rewind(file_coils);
	while (fgets(line, lenline, file_coils) != NULL)
	{	
		storevals = linetodata(line, &nn);
		if (nn > 3)
		{
			//printf("hello world, number of coils = %d, coil_index = %d, coilseg_index = %d\n", num_coils, coil_index, coilseg_index);
			XXvect_coil[0][coil_index][coilseg_index] = *storevals; 
			XXvect_coil[1][coil_index][coilseg_index] = *(storevals+1);
			XXvect_coil[2][coil_index][coilseg_index] = *(storevals+2);
			XXvect_coil[3][coil_index][coilseg_index] = *(storevals+3);
			//printf("coil_index=%d, coilseg_index = %d\n", coil_index, coilseg_index);
			//printf("storevals[0] = %Le\n", *storevals);
			//printf("storevals[1] = %Le\n", *(storevals+1));
			//printf("storevals[2] = %Le\n", *(storevals+2));
			//printf("storevals[3] = %Le\n", *(storevals+3));
			coilseg_index += 1; 
		}
		if (nn > 4)
		{
			coil_index += 1;
			coilseg_index = 0;
		}
		free(storevals);
	}
	fclose(file_coils);
	*num_coils = num_coils_temp;
	*num_segs = num_segs_temp;
	return XXvect_coil;
}
		
struct field *BReim(int m0_symmetry, double iota0, double iota1, double *epsilon, int *k_theta, int size_epsilon, double RR, double ZZ, double varphi) { // all working
		int ind;
		double R_axis = 1.0, theta, rmin, combo, combo1, dcombodR, dcombodZ, dcombo1dR, dcombo1dZ; // check_q;
		double d2combodRR, d2combodRZ, d2combodZZ, d2combo1dRR, d2combo1dRZ, d2combo1dZZ; // check_q;
		//double check, checkdlen, dlen=0.000001; // uncomment when wishing to check gradients
		double quantity, dquantitydZ=0.0;
		struct field *mag = calloc(1,sizeof(struct field));
		theta = atan2(ZZ, RR - R_axis);
		rmin = sqrt(pow((RR-R_axis), 2.0) + pow(ZZ, 2.0));
		combo = iota0 + iota1*rmin*rmin;
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

		mag->value[0] = ( ZZ / RR ) *  combo + ( ( RR - R_axis ) / RR ) * combo1 ;
		mag->value[1] = - ( (RR - R_axis) / RR ) * combo + ( ZZ / RR ) * combo1  ; 
		mag->value[2] = -1.0;
		mag->derivative[0][0] = ( - ZZ / pow( RR, 2.0) ) *  combo + (ZZ / RR) * dcombodR + combo1 * R_axis / pow(RR, 2.0) + dcombo1dR * (RR - R_axis) / RR;
 //+ gradR*( (2.0*ZZ/pow(RR, 3.0))*combo0 - (2.0*ZZ/pow(RR, 2.0))*
		mag->derivative[0][1] = ( 1.0 / RR ) * combo + (ZZ / RR) * dcombodZ + dcombo1dZ * (RR - R_axis) / RR;
		mag->derivative[1][0] = ( - R_axis / pow( RR, 2.0) ) *  combo - ((RR - R_axis) / RR) * dcombodR - combo1 * ZZ / pow(RR, 2.0) + dcombo1dR * ZZ / RR;
		mag->derivative[1][1] = - ((RR - R_axis) / RR) * dcombodZ + combo1 * ( 1.0 / RR ) + dcombo1dZ * ZZ / RR ;
		mag->derivative[2][0] = 0.0;
		mag->derivative[2][1] = 0.0;
		// right
		mag->twoderivative[0][0][0] = (  2*ZZ / pow( RR, 3.0) ) *combo - (2.0*R_axis/pow(RR, 3.0))*combo1 - (2*ZZ / pow(RR, 2.0)) * dcombodR +  (2*R_axis / pow(RR, 2.0) )*dcombo1dR  + (ZZ/RR)*d2combodRR + ((RR-R_axis)/RR)*d2combo1dRR; 
		// right
		mag->twoderivative[0][0][1] = - (  1.0 / pow( RR, 2.0) ) *combo - (ZZ/pow(RR, 2.0))*dcombodZ + (1.0/ RR) * dcombodR +  (R_axis / pow(RR, 2.0) )*dcombo1dZ  + (ZZ/RR)*d2combodRZ + ((RR-R_axis)/RR)*d2combo1dRZ; 
		// right
		mag->twoderivative[0][1][0] = - (  1.0 / pow( RR, 2.0) ) *combo - (ZZ/pow(RR, 2.0))*dcombodZ + (1.0/ RR) * dcombodR +  (R_axis / pow(RR, 2.0) )*dcombo1dZ  + (ZZ/RR)*d2combodRZ + ((RR-R_axis)/RR)*d2combo1dRZ; 
		// right
		mag->twoderivative[0][1][1] = ( 2.0 / RR ) *dcombodZ +  (ZZ/RR)*d2combodZZ + ((RR-R_axis)/RR)*d2combo1dZZ;
		// right
		mag->twoderivative[1][0][0] = (  2*R_axis / pow( RR, 3.0) ) *combo + (2.0*ZZ/pow(RR, 3.0))*combo1 - (2*ZZ / pow(RR, 2.0)) * dcombo1dR -  (2*R_axis / pow(RR, 2.0) )*dcombodR  + (ZZ/RR)*d2combo1dRR - ((RR-R_axis)/RR)*d2combodRR; 
		//right
		mag->twoderivative[1][1][0] = mag->twoderivative[1][0][1] = - (  1.0 / pow( RR, 2.0) ) *combo1 - (R_axis/pow(RR, 2.0))*dcombodZ - (ZZ/pow(RR, 2.0))*dcombo1dZ + (1.0/ RR)*dcombo1dR - ((RR - R_axis) / RR)*d2combodRZ  + (ZZ/RR)*d2combo1dRZ;
		//right
		mag->twoderivative[1][1][1] = ( 2.0 / RR ) *dcombo1dZ +  (ZZ/RR)*d2combo1dZZ - ((RR-R_axis)/RR)*d2combodZZ;
		mag->twoderivative[2][0][0] = 0.0;
		mag->twoderivative[2][1][0] = mag->twoderivative[2][0][1] = 0.0;
		mag->twoderivative[2][1][1] = 0.0;

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
		return mag;
}

struct field *gradBReim(int m0_symmetry, double iota0, double iota1, double *epsilon, int *k_theta, int size_epsilon, double RR, double ZZ, double varphi) {
		int ind;
		double R_axis = 1.0, theta, rmin; // check_q;
		double dcomboiota, dcomboiota1, dcomboepsilon[size_epsilon], dcombo1epsilon[size_epsilon];
		double d_dcombodR_epsilon[size_epsilon], d_dcombodZ_epsilon[size_epsilon], d_dcombo1dR_epsilon[size_epsilon], d_dcombo1dZ_epsilon[size_epsilon];
		struct field *mag = malloc(size_epsilon*sizeof(struct field));
		theta = atan2(ZZ, RR - R_axis);
		rmin = sqrt(pow((RR-R_axis), 2.0) + pow(ZZ, 2.0));
		//combo = iota0 + iota1*rmin*rmin - 2.0*epsilon[0]*cos(2.0*theta - varphi) - 3.0*epsilon[1]*rmin*cos(3.0*theta - varphi);
		//combo1 = 2.0*epsilon[0]*rmin*rmin*sin(2.0*theta - varphi) + 3.0*epsilon[1]*pow(rmin, 3.0)*sin(3.0*theta - varphi);

		dcomboiota = 1.0;
		dcomboiota1 = rmin*rmin;
		//mag[0].value[0] = ( ZZ / RR ) *  dcomboiota ;
		//mag[0].value[1] = - ( (RR - R_axis) / RR ) * dcomboiota;
		//mag[0].value[2] = 0.0;
		//mag[1].value[0] = ( ZZ / RR ) *  dcomboiota1 ;
		//mag[1].value[1] = - ( (RR - R_axis) / RR ) * dcomboiota1;
		//mag[1].value[2] = 0.0;
		for (ind=0; ind < size_epsilon; ind++) {
			dcomboepsilon[ind]  = - k_theta[ind] * pow(rmin, k_theta[ind] - 2) * cos(k_theta[ind]*theta - m0_symmetry*varphi);
			dcombo1epsilon[ind] =   k_theta[ind] * pow(rmin, k_theta[ind] - 2) * sin(k_theta[ind]*theta - m0_symmetry*varphi);
			mag[ind].value[0] = (ZZ / RR) * dcomboepsilon[ind] + dcombo1epsilon[ind] * (RR - R_axis) / RR;
			mag[ind].value[1] = - ( (RR - R_axis) / RR ) * dcomboepsilon[ind] + dcombo1epsilon[ind] * ZZ / RR;
			mag[ind].value[2] = 0.0;
			d_dcombodR_epsilon[ind] = - k_theta[ind] * pow( rmin, k_theta[ind] - 4 ) * ( k_theta[ind] * sin( k_theta[ind]*theta - m0_symmetry*varphi ) * ZZ + (k_theta[ind] - 2) *  (RR - R_axis) * cos(k_theta[ind]*theta - m0_symmetry * varphi) );
			d_dcombodZ_epsilon[ind] = - pow( rmin, k_theta[ind] - 4 ) * k_theta[ind] * (  - k_theta[ind] * sin( k_theta[ind]*theta - m0_symmetry*varphi ) * ( RR - R_axis ) + (k_theta[ind] - 2) * ZZ * cos(k_theta[ind]*theta - m0_symmetry * varphi) );
			d_dcombo1dR_epsilon[ind] = k_theta[ind] * pow(rmin, k_theta[ind] - 4) * ( - k_theta[ind] * cos(k_theta[ind]*theta - m0_symmetry*varphi) * ZZ  +  (k_theta[ind] - 2) * sin(k_theta[ind]*theta - m0_symmetry*varphi) * (RR - R_axis) ) ;
			d_dcombo1dZ_epsilon[ind] = k_theta[ind] * pow(rmin, k_theta[ind] - 4) * ( k_theta[ind] * cos(k_theta[ind]*theta - m0_symmetry *varphi) * (RR - R_axis) + (k_theta[ind] - 2) * sin(k_theta[ind]*theta - m0_symmetry*varphi) * ZZ ) ;
			mag[ind].derivative[0][0] = ( - ZZ / pow( RR, 2.0) ) *  dcomboepsilon[ind] + (ZZ / RR) * d_dcombodR_epsilon[ind] + dcombo1epsilon[ind] * R_axis / pow(RR, 2.0) + d_dcombo1dR_epsilon[ind] * (RR - R_axis) / RR ;
			mag[ind].derivative[1][0] = - ( R_axis / pow( RR, 2.0) ) *  dcomboepsilon[ind] - ((RR - R_axis) / RR) * d_dcombodR_epsilon[ind] - dcombo1epsilon[ind] * ( ZZ / pow(RR, 2.0) ) + d_dcombo1dR_epsilon[ind] * ZZ / RR;
			mag[ind].derivative[2][0] = 0.0;
			mag[ind].derivative[0][1] = ( 1.0 / RR ) * dcomboepsilon[ind] + (ZZ / RR) * d_dcombodZ_epsilon[ind] + d_dcombo1dZ_epsilon[ind] * ( RR - R_axis ) / RR;
			mag[ind].derivative[1][1] = - ((RR - R_axis) / RR) * d_dcombodZ_epsilon[ind] + dcombo1epsilon[ind] * ( 1.0 / RR ) + d_dcombo1dZ_epsilon[ind] * ZZ / RR;
			mag[ind].derivative[2][1] = 0.0;
	//		mag[ind].derivative[0][0] =  (ZZ / RR) * d_dcombodR_epsilon[ind] + dcombo1epsilon[ind] * R_axis / pow(RR, 2.0) + d_dcombo1dR_epsilon[ind] * (RR - R_axis) / RR ;
	//		mag[ind].derivative[1][0] =  - ((RR - R_axis) / RR) * d_dcombodR_epsilon[ind] - dcombo1epsilon[ind] * ( ZZ / pow(RR, 2.0) ) + d_dcombo1dR_epsilon[ind] * ZZ / RR;
	//		mag[ind].derivative[2][0] = 0.0;
	//		mag[ind].derivative[0][1] = (ZZ / RR) * d_dcombodZ_epsilon[ind] + d_dcombo1dZ_epsilon[ind] * ( RR - R_axis ) / RR;
	//		mag[ind].derivative[1][1] = - ((RR - R_axis) / RR) * d_dcombodZ_epsilon[ind]  + d_dcombo1dZ_epsilon[ind] * ZZ / RR;
	//		mag[ind].derivative[2][1] = 0.0;
		}
		//mag[0].value[0] = (ZZ / RR) * dcomboiota; mag[0].value[1] = - ( (RR - R_axis) / RR ) * dcomboiota; mag[0].value[2] = 0.0;
		//mag[0].value[0] = (ZZ / RR) * dcomboiota1; mag[0].value[1] = - ( (RR - R_axis) / RR ) * dcomboiota1; mag[0].value[2] = 0.0;
		return mag;
}

struct field *Bfield(double *Xp, double varphi, double ***coils, int *num_coils, int **num_segs) {
	//clock_t start = clock();
	/* below: declare B and gradB calculated with the old way */
	//double BZ=0.0, BX=0.0, BY=0.0, BR=0.0, Bvarphi=0.0, dBXdR=0.0, dBYdR=0.0, dBZdR=0.0, dBXdZ=0.0, dBYdZ=0.0, dBZdZ=0.0, dBRdR=0.0, dBvarphidR=0.0, dBRdZ=0.0, dBvarphidZ=0.0;
	/* above: declare B and gradB calculated with the old way */
	double Xminxc=0.0, Yminyc=0.0, Zminzc=0.0, Rminrc=0.0, dxc=0.0, dyc=0.0, dzc=0.0;
	double XX=Xp[0]*cos(varphi), YY=Xp[0]*sin(varphi), ZZ = Xp[1];
	double divB;
	struct field *magfield = calloc(1,sizeof(struct field)), *magfieldpdR = calloc(1,sizeof(struct field)), *magfieldpdZ = calloc(1,sizeof(struct field)), *magfieldpdvarphi = calloc(1, sizeof(struct field));
	//struct field *magfieldmdR = calloc(1,sizeof(struct field)), *magfieldmZ = calloc(1,sizeof(struct field));
	//struct field *magfieldmdRpdZ = calloc(1,sizeof(struct field)), *magfieldmdZpdR = calloc(1,sizeof(struct field));
	//struct field *magfieldmdRmdZ = calloc(1,sizeof(struct field)), *magfieldpdZpdR = calloc(1,sizeof(struct field));
	struct field dmagfielddR, dmagfielddZ, dmagfielddvarphi;
	//struct field d2magfielddRR, d2magfieldZZ;
	//struct field d2magfielddRZ, d2magfieldZR;
	int coil_index=0, coilseg_index=0; 
	int check = 0, m0_symmetry = 1, row, col, dep;
	//double XpdR, YpdR, RpdRminrc, XpdRminxc, YpdRminyc, BXpdR, BYpdR, BZpdR, 
	double dR[2], dZ[2], delta = 0.0000001, dvarphi = 0.00001;
	double Rminrcvect[3], drcvect[3], Rhat[3], Zhat[3], phihat[3], Bvect[3];
	double gradB[3][2], gradgradB[3][2][2], object, gradBcheck[3][2], gradgradBcheck[3][2][2]; 
	//double RpdZminrc, ZpdZminzc, BXpdZ, BYpdZ, BZpdZ;
	double amp_Domm[1], epsilon[1], iota[2]; 
	int pol_Domm[1], tor_Domm[1], k_theta[1];
	char *type = "coil";
	for (row=0;row<3;row++) {
		Bvect[row] = 0.0; 
		for (col=0;col<2;col++) {
			gradB[row][col] = 0.0;
			gradBcheck[row][col] = 0.0;
			for (dep=0;dep<2;dep++) {
				gradgradB[row][col][dep] = 0.0;
				gradgradBcheck[row][col][dep] = 0.0;
			}
		}
	}
	//printf("n_segs[1]=%d\n", n_segs[1]);
	//printf("Sanity check: print line of coil file\n");
	//printf("coils[0][%d][%d]=%f\n", j, k, coils[0][j][k]); printf("coils[1][%d][%d]=%f\n", j, k, coils[1][j][k]);
	//printf("coils[2][%d][%d]=%f\n", j, k, coils[2][j][k]); printf("coils[3][%d][%d]=%f\n", j, k, coils[3][j][k]);
	//printf("Sanity check: print another line of coil file\n");
	//printf("coils[0][0][1]=%Le\n", coils[0][0][1]); printf("coils[1][0][1]=%Le\n", coils[1][0][1]);
	//printf("coils[2][0][1]=%Le\n", coils[2][0][1]); printf("coils[3][0][1]=%Le\n", coils[3][0][1]);
	if (strncmp(type, "coil", 4) == 0) {
		Rhat[0] = 1.0; Rhat[1] = 0.0; Rhat[2] = 0.0;
		Zhat[0] = 0.0; Zhat[1] = 1.0; Zhat[2] = 0.0;
		phihat[0] = 0.0; phihat[1] = 0.0; phihat[2] = 1.0;
		for (coil_index=0;coil_index<*num_coils;coil_index++)
		{
			//printf("num_coils = %d\n", *num_coils);
			//printf("coil_index = %d\n", coil_index);
			for (coilseg_index=0;coilseg_index<*(*num_segs+coil_index)-1;coilseg_index++)
			{
				//printf("num_segs[%d] = %d\n", coil_index, *(*num_segs+coil_index));
				//printf("coilseg_index = %d\n", coilseg_index);
				/* extract useful info */
				Xminxc = XX - 0.5*coils[0][coil_index][coilseg_index] - 0.5*coils[0][coil_index][coilseg_index+1];
				Yminyc = YY - 0.5*coils[1][coil_index][coilseg_index] - 0.5*coils[1][coil_index][coilseg_index+1]; 
				Zminzc = ZZ - 0.5*coils[2][coil_index][coilseg_index] - 0.5*coils[2][coil_index][coilseg_index+1]; 
				Rminrc = sqrt(pow(Xminxc,  2.0) + pow(Yminyc, 2.0) + pow(Zminzc, 2.0)); 

				dxc = coils[0][coil_index][coilseg_index+1] - coils[0][coil_index][coilseg_index];
				dyc = coils[1][coil_index][coilseg_index+1] - coils[1][coil_index][coilseg_index];
				dzc = coils[2][coil_index][coilseg_index+1] - coils[2][coil_index][coilseg_index];

				/* project into cartesian axes aligned with local poloidal plane 
				Rhat dot (phihat cross Zhat) = 1 prescription -> (R,phi,Z) right-handed */
				Rminrcvect[0] = cos(varphi)*Xminxc + sin(varphi)*Yminyc;
				Rminrcvect[1] = Zminzc;
				Rminrcvect[2] = - sin(varphi)*Xminxc + cos(varphi)*Yminyc;
				drcvect[0] = cos(varphi)*dxc + sin(varphi)*dyc;
				drcvect[1] = dzc;
				drcvect[2] = - sin(varphi)*dxc + cos(varphi)*dyc; 

				/* below: old way to calculate B and grad B (can soon cancel it) */
				//BX += (coils[3][coil_index][coilseg_index]*(dyc*Zminzc - dzc*Yminyc)/pow(Rminrc, 3.0));
				//BY += (coils[3][coil_index][coilseg_index]*(dzc*Xminxc - dxc*Zminzc)/pow(Rminrc, 3.0));
				//BZ += (coils[3][coil_index][coilseg_index]*(dxc*Yminyc - dyc*Xminxc)/pow(Rminrc, 3.0));

				//
				//dBXdR += (coils[3][coil_index][coilseg_index]*(dyc*0.0         - dzc*sin(varphi))/pow(Rminrc, 3.0));
				//dBYdR += (coils[3][coil_index][coilseg_index]*(dzc*cos(varphi) - dxc*0.0        )/pow(Rminrc, 3.0));
				//dBZdR += (coils[3][coil_index][coilseg_index]*(dxc*sin(varphi) - dyc*cos(varphi))/pow(Rminrc, 3.0));

				//dBXdR -= 3*(coils[3][coil_index][coilseg_index]*(Xminxc*cos(varphi) + Yminyc*sin(varphi))*(dyc*Zminzc - dzc*Yminyc)/pow(Rminrc, 5.0));
				//dBYdR -= 3*(coils[3][coil_index][coilseg_index]*(Xminxc*cos(varphi) + Yminyc*sin(varphi))*(dzc*Xminxc - dxc*Zminzc)/pow(Rminrc, 5.0));
				//dBZdR -= 3*(coils[3][coil_index][coilseg_index]*(Xminxc*cos(varphi) + Yminyc*sin(varphi))*(dxc*Yminyc - dyc*Xminxc)/pow(Rminrc, 5.0));

				//dBXdZ += (coils[3][coil_index][coilseg_index]*( dyc)/pow(Rminrc, 3.0));
				//dBYdZ += (coils[3][coil_index][coilseg_index]*(-dxc)/pow(Rminrc, 3.0));

				//dBXdZ -= 3*(coils[3][coil_index][coilseg_index]*(Zminzc)*(dyc*Zminzc - dzc*Yminyc)/pow(Rminrc, 5.0));
				//dBYdZ -= 3*(coils[3][coil_index][coilseg_index]*(Zminzc)*(dzc*Xminxc - dxc*Zminzc)/pow(Rminrc, 5.0));
				//dBZdZ -= 3*(coils[3][coil_index][coilseg_index]*(Zminzc)*(dxc*Yminyc - dyc*Xminxc)/pow(Rminrc, 5.0));
				// above: old way to calculate B and grad B (can soon cancel it)

				for (row=0;row<3;row++) {
					Bvect[row] +=  pow(10.0,-7.0) * ( coils[3][coil_index][coilseg_index] / pow(Rminrc, 3.0) ) * ( Rhat[row] * (-drcvect[1]*Rminrcvect[2] + drcvect[2]*Rminrcvect[1] ) + Zhat[row] * (-drcvect[2]*Rminrcvect[0] + drcvect[0]*Rminrcvect[2] ) + phihat[row] * (-drcvect[0]*Rminrcvect[1] + drcvect[1]*Rminrcvect[0] ) );
					for (col=0;col<2;col++) {
						gradB[row][col] += ( pow(10.0,-7.0) * ( coils[3][coil_index][coilseg_index] / pow(Rminrc, 3.0) ) * ( - 3.0 * Rminrcvect[col] * ( Rhat[row] * (-drcvect[1]*Rminrcvect[2] + drcvect[2]*Rminrcvect[1] ) + Zhat[row] * (-drcvect[2]*Rminrcvect[0] + drcvect[0]*Rminrcvect[2] ) + phihat[row] * (-drcvect[0]*Rminrcvect[1] + drcvect[1]*Rminrcvect[0] ) ) / pow(Rminrc, 2.0) + ( drcvect[2]*( Rhat[row]*Zhat[col] - Rhat[col]*Zhat[row])  + drcvect[1] * phihat[row] * Rhat[col] - drcvect[0] * phihat[row] *Zhat[col] ) ) ) ; 
						for (dep=0; dep<2; dep++) {
							object = 3.0 * ( - Rhat[dep] *Rhat[col] - Zhat[dep]* Zhat[col] + 5.0 * Rminrcvect[dep] * Rminrcvect[col] / pow(Rminrc, 2.0) ) * ( Rhat[row] * ( - drcvect[1]*Rminrcvect[2] + drcvect[2]*Rminrcvect[1] ) + Zhat[row] * (-drcvect[2]*Rminrcvect[0] + drcvect[0]*Rminrcvect[2] ) + phihat[row] * (-drcvect[0]*Rminrcvect[1] + drcvect[1]*Rminrcvect[0] ) );
							object += 3.0 * Rminrcvect[dep]* ( drcvect[2]*( Rhat[col]*Zhat[row] - Rhat[row]*Zhat[col])  - drcvect[1] * phihat[row] * Rhat[col] + drcvect[0] * phihat[row] *Zhat[col] ) ; 
							object += 3.0 * Rhat[dep] * Rminrcvect[col] * (   drcvect[2] * Zhat[row] - drcvect[1] * phihat[row] );
							object += 3.0 * Zhat[dep] * Rminrcvect[col] * ( - drcvect[2]* Rhat[row] + drcvect[0] * phihat[row]  );
							gradgradB[row][col][dep] += ( pow(10.0,-7.0) * ( coils[3][coil_index][coilseg_index] / pow(Rminrc, 5.0) ) * object ) ;
						}
					}
				}
				
//				for (col=0;col<2;col++) {
//					gradB[2][col] += ( pow(10.0,-7.0) * coils[3][coil_index][coilseg_index] * ( - 3.0 * Rminrcvect[col] * ( - drcvect[0]*Rminrcvect[1] + drcvect[1]*Rminrcvect[0] ) / pow(Rminrc, 5.0) + ( drcvect[1] * Rhat[col] - drcvect[0] * Zhat[col] ) / pow(Rminrc, 3.0) ) ) ; 
//				}
//
				if (check==1) {
					for (dep=0;dep<2;dep++) {
						if (dep==0) {
							dR[dep] = delta;
							dZ[dep] = 0.0;
						}
						else if (dep==1) {
							dR[dep] = 0.0;
							dZ[dep] = delta;
						}

						Rminrcvect[0] = cos(varphi)*Xminxc + sin(varphi)*Yminyc + dR[dep];
						Rminrcvect[1] = Zminzc + dZ[dep];
						Rminrcvect[2] = - sin(varphi)*Xminxc + cos(varphi)*Yminyc;
						Rminrc = pow(pow(Rminrcvect[0], 2.0) + pow(Rminrcvect[1], 2.0) + pow(Rminrcvect[2], 2.0), 0.5);
						for (row=0;row<3;row++) {
							gradBcheck[row][dep] +=  pow(10.0,-7.0) * ( coils[3][coil_index][coilseg_index] / pow(Rminrc, 3.0) ) * ( Rhat[row] * (-drcvect[1]*Rminrcvect[2] + drcvect[2]*Rminrcvect[1] ) + Zhat[row] * (-drcvect[2]*Rminrcvect[0] + drcvect[0]*Rminrcvect[2] ) + phihat[row] * (-drcvect[0]*Rminrcvect[1] + drcvect[1]*Rminrcvect[0] ) );
							for (col=0;col<2;col++) {
								gradgradBcheck[row][col][dep] += ( pow(10.0,-7.0) * ( coils[3][coil_index][coilseg_index] / pow(Rminrc, 3.0) ) * ( - 3.0 * Rminrcvect[col] * ( Rhat[row] * (-drcvect[1]*Rminrcvect[2] + drcvect[2]*Rminrcvect[1] ) + Zhat[row] * (-drcvect[2]*Rminrcvect[0] + drcvect[0]*Rminrcvect[2] ) + phihat[row] * (-drcvect[0]*Rminrcvect[1] + drcvect[1]*Rminrcvect[0] ) ) / pow(Rminrc, 2.0) + ( drcvect[2]*( Rhat[row]*Zhat[col] - Rhat[col]*Zhat[row])  + drcvect[1] * phihat[row] * Rhat[col] - drcvect[0] * phihat[row] *Zhat[col] ) ) ) ; 
							}
						}
					}
				}
			}
		}
		/* below: final calculation of B and gradB using old way */
		//BX*=pow(10.0,-7.0); BY*=pow(10.0,-7.0); BZ*=pow(10.0,-7.0); // multiply by mu_0/4*PI
		//BR = cos(varphi)*BX + sin(varphi)*BY; Bvarphi = -sin(varphi)*BX + cos(varphi)*BY; // BZ = BZ;
		//dBXdR*=pow(10.0,-7.0); dBYdR*=pow(10.0,-7.0); dBZdR*=pow(10.0,-7.0); // multiply by mu_0/4*PI
		//dBXdZ*=pow(10.0,-7.0); dBYdZ*=pow(10.0,-7.0); dBZdZ*=pow(10.0,-7.0); // multiply by mu_0/4*PI
		//dBRdR = dBXdR*cos(varphi) + dBYdR*sin(varphi); dBvarphidR = -dBXdR*sin(varphi) + dBYdR*cos(varphi);
		//dBRdZ = dBXdZ*cos(varphi) + dBYdZ*sin(varphi); dBvarphidZ = -dBXdZ*sin(varphi) + dBYdZ*cos(varphi);
		//magfield->value[0] = BR; magfield->value[1] = BZ; magfield->value[2] = Bvarphi;
		//magfield->derivative[0][0] = dBRdR; 
		//magfield->derivative[1][0] = dBZdR; 
		//magfield->derivative[2][0] = dBvarphidR;
		//magfield->derivative[0][1] = dBRdZ; 
		//magfield->derivative[1][1] = dBZdZ; 
		//magfield->derivative[2][1] = dBvarphidZ;
		/* above: final calculation of B and gradB using old way */
		/* below: comparison with old way */
		//printf("dBRdR=%f and %f\n", dBRdR, gradB[0][0]);
		//printf("dBRdZ=%f and %f\n", dBRdZ, gradB[0][1]);
		//printf("dBZdR=%f and %f\n", dBZdR, gradB[1][0]);
		//printf("dBZdZ=%f and %f\n", dBZdZ, gradB[1][1]);
		//printf("dBphidR=%f and %f\n", dBvarphidR, gradB[2][0]);
		//printf("dBphidZ=%f and %f\n", dBvarphidZ, gradB[2][1]);
		//printf("BR=%f and %f\n", BR, Bvect[0]);
		//printf("BZ=%f and %f\n", BZ, Bvect[1]);
		//printf("Bvarphi=%f and %f\n", Bvarphi, Bvect[2]);
		/* above: comparison with old way */

		/* assign Bfield structure the correct values of B, gradB and gradgradB */
		/* if check == 1 then check with numerical derivatives */
		for (row=0;row<3;row++) {
			magfield->value[row] = Bvect[row]; 
			for (col=0;col<2;col++) {
				magfield->derivative[row][col] = gradB[row][col]; 
				if (check==1) {
					gradBcheck[row][col] = (gradBcheck[row][col] - Bvect[row])/delta;
					printf("gradB[%d][%d] = %f (check against %f)\n", row, col, gradB[row][col], gradBcheck[row][col]);
				}
				for (dep=0;dep<2;dep++) {
					magfield->twoderivative[row][col][dep] = gradgradB[row][col][dep]; 
					if (check==1) {
					
						gradgradBcheck[row][col][dep] = (gradgradBcheck[row][col][dep] - gradB[row][col])/delta;
						printf("gradgradB[%d][%d][%d] = %f (check against %f)\n", row, col, dep,  gradgradB[row][col][dep], gradgradBcheck[row][col][dep]);
					}
				}
			}
		}
			//dBXdR = (BXpdR - BX)/dR;
			//dBYdR = (BYpdR - BY)/dR;
			//dBZdR = (BZpdR - BZ)/dR;
			//dBXdZ = (BXpdZ - BX)/dZ;
			//dBYdZ = (BYpdZ - BY)/dZ;
			//dBZdZ = (BZpdZ - BZ)/dZ;
			//dBXdR*=pow(10.0,-7.0); dBYdR*=pow(10.0,-7.0); dBZdR*=pow(10.0,-7.0); // multiply by mu_0/4*PI
			//dBXdZ*=pow(10.0,-7.0); dBYdZ*=pow(10.0,-7.0); dBZdZ*=pow(10.0,-7.0); // multiply by mu_0/4*PI
			//dBRdR = dBXdR*cos(varphi) + dBYdR*sin(varphi); dBvarphidR = -dBXdR*sin(varphi) + dBYdR*cos(varphi);
			//dBRdZ = dBXdZ*cos(varphi) + dBYdZ*sin(varphi); dBvarphidZ = -dBXdZ*sin(varphi) + dBYdZ*cos(varphi);

			//printf("%f %f\n", magfield->derivative[0][0], dBRdR);
			//printf("%f %f\n", magfield->derivative[1][0], dBZdR);
			//printf("%f %f\n", magfield->derivative[2][0], dBvarphidR);
			//printf("%f %f\n", magfield->derivative[0][1], dBRdZ);
			//printf("%f %f\n", magfield->derivative[1][1], dBZdZ);
			//printf("%f %f\n", magfield->derivative[2][1], dBvarphidZ);
	}
	//clock_t int3 = clock();
	//printf("Time after evaluating B field at a point: %f\n", (double) (int3-start)/CLOCKS_PER_SEC);
	else if (strncmp(type, "Reim", 4) == 0) {
		epsilon[0] = 0.1; //epsilon[1] = 0.0;
		//epsilon[0] = 0.0; epsilon[1] = 0.0;
		k_theta[0] = 6; //k_theta[1] = 3;
		//epsilon[0] = 0.00; epsilon[1] = 0.00;
		iota[0] = 0.15; iota[1] = 0.38;
		magfield = BReim(m0_symmetry, iota[0], iota[1], epsilon, k_theta, 1, Xp[0], Xp[1], varphi);
		//theta = atan2(Xp[1], Xp[0] - R_axis);
		//rmin = sqrt(pow((Xp[0]-R_axis), 2.0) + pow(Xp[1], 2.0));
		//combo = iota[0] + iota[1]*rmin*rmin - 2.0*epsilon[0]*cos(2.0*theta - varphi) - 3.0*epsilon[1]*rmin*cos(3.0*theta - varphi);
		//combo1 = 2.0*epsilon[0]*rmin*rmin*sin(2.0*theta - varphi) + 3.0*epsilon[1]*pow(rmin, 3.0)*sin(3.0*theta - varphi);
		//BR = ( Xp[1] / Xp[0] ) * ( combo + combo1 / (pow(Xp[1], 2.0) + pow(Xp[0] - R_axis, 2.0)) ) ;
		//BZ = - ( (Xp[0] - R_axis) / Xp[0] ) * ( combo + combo1 / (pow(Xp[1], 2.0) + pow(Xp[0] - R_axis, 2.0)) ) ; 
		//Bvarphi = -1.0; 
		//printf("theta = %f\t rmin = %f\n", theta, rmin);
		//printf("BR = %f\tBZ = %f\tBvarphi = %f\n", BR, BZ, Bvarphi);

		//BR = ( Xp[1] / ( Xp[0] + dR ) ) * ( combo + combo1 / (pow(Xp[1], 2.0) + pow(Xp[0] - R_axis, 2.0)) ) ;
		//BZ = - ( (Xp[0] - R_axis) / Xp[0] ) * ( combo + combo1 / (pow(Xp[1], 2.0) + pow(Xp[0] - R_axis, 2.0)) ) ; 
		//dBRdR = 0.0;
		//dBZdR = 0.0;
		//dBvarphidR = 0.0;
		//dBRdZ = 0.0;
		//dBZdZ = 0.0;
		//dBvarphidZ = 0.0;
		//printf("work in progress\n");	
		//magfield->value[0] = BR; 
		//magfield->value[1] = BZ; 
		//magfield->value[2] = Bvarphi;
		//magfield->derivative[0][0] = dBRdR; 
		//magfield->derivative[1][0] = dBZdR; 
		//magfield->derivative[2][0] = dBvarphidR;
		//magfield->derivative[0][1] = dBRdZ; 
		//magfield->derivative[1][1] = dBZdZ; 
		//magfield->derivative[2][1] = dBvarphidZ;
		if (check == 1) {
			magfieldpdR      = BReim(m0_symmetry, iota[0], iota[1], epsilon, k_theta, 1, Xp[0]+delta, Xp[1], varphi);
			magfieldpdZ      = BReim(m0_symmetry, iota[0], iota[1], epsilon, k_theta, 1, Xp[0], Xp[1]+delta, varphi);
			magfieldpdvarphi = BReim(m0_symmetry, iota[0], iota[1], epsilon, k_theta, 1, Xp[0], Xp[1], varphi + dvarphi);
			//BX = 0.0; BY = 0.0; BZ = 0.0;
			//BXpdR = 0.0; BYpdR = 0.0; BZpdR = 0.0;
			//BXpdZ = 0.0; BYpdZ = 0.0; BZpdZ = 0.0;
			//dBXdR = 0.0; dBYdR = 0.0; dBZdR = 0.0;
			//dBXdZ = 0.0; dBYdZ = 0.0; dBZdZ = 0.0;
			//BXpdR = 0.0; BYpdR = 0.0; BZpdR = 0.0;
			dmagfielddR = addstructsfield(1/delta, magfieldpdR, -1/delta, magfield);
			dmagfielddZ = addstructsfield(1/delta, magfieldpdZ, -1/delta, magfield);
			dmagfielddvarphi = addstructsfield(1/dvarphi, magfieldpdvarphi, -1/dvarphi, magfield);
			printf("dBRdR   %f %f\n", dmagfielddR.value[0], magfield->derivative[0][0]);
			printf("dBZdR   %f %f\n", dmagfielddR.value[1], magfield->derivative[1][0]);
			printf("dBphidR %f %f\n", dmagfielddR.value[2], magfield->derivative[2][0]);
			printf("dBRdZ   %f %f\n", dmagfielddZ.value[0], magfield->derivative[0][1]);
			printf("dBZdZ   %f %f\n", dmagfielddZ.value[1], magfield->derivative[1][1]);
			printf("dBphidZ %f %f\n", dmagfielddZ.value[2], magfield->derivative[2][1]);
			//printf("%f %f\n", magfield->derivative[1][0], dBZdR);
			//printf("%f %f\n", magfield->derivative[2][0], dBvarphidR);
			//printf("%f %f\n", magfield->derivative[0][1], dBRdZ);
			//printf("%f %f\n", magfield->derivative[1][1], dBZdZ);
			//printf("%f %f\n", magfield->derivative[2][1], dBvarphidZ);
			//divb = ( dmagfielddr.value[0] + magfield->value[0] / xp[0] + dmagfielddz.value[1] ) / sqrt( pow(magfield->value[0], 2.0) + pow(magfield->value[1], 2.0) );
			divB = ( dmagfielddR.value[0] + magfield->value[0] / Xp[0] + dmagfielddZ.value[1] + dmagfielddvarphi.value[2] / Xp[0] ) / sqrt( pow( dmagfielddR.value[0] + magfield->value[0] / Xp[0] , 2.0 ) + pow( dmagfielddZ.value[1], 2.0 ) + pow( dmagfielddvarphi.value[2] / Xp[0] , 2.0 ) );
			printf("div(B) = %f\n", divB);
			//printf("%f %f\n", magfield->derivative[1][0], dBZdR);
			//printf("%f %f\n", magfield->derivative[2][0], dBvarphidR);
			//printf("%f %f\n", magfield->derivative[0][1], dBRdZ);
			//printf("%f %f\n", magfield->derivative[1][1], dBZdZ);
			//printf("%f %f\n", magfield->derivative[2][1], dBvarphidZ);
		}
	}
	else if (strncmp(type, "Domm", 4) == 0) {
		amp_Domm[0]=1.73; pol_Domm[0] = 2; tor_Domm[0] = 5;
		//amp_Domm[0]=2.0; pol_Domm[0] = 5; tor_Domm[0] = 5;
		//amp_Domm[0]=0.001; pol_Domm[0] = 4; tor_Domm[0] = 8;
		//amp_Domm[0]=0.00000; pol_Domm[0] = 2; tor_Domm[0] = 5;
		magfield = DommBfield(1, amp_Domm, tor_Domm, pol_Domm, Xp[0], Xp[1], varphi);
		if (check==1) {
			magfieldpdR = DommBfield(1, amp_Domm, tor_Domm, pol_Domm, Xp[0]+delta, Xp[1], varphi);
			magfieldpdZ = DommBfield(1, amp_Domm, tor_Domm, pol_Domm, Xp[0], Xp[1]+delta, varphi);
			magfieldpdvarphi = DommBfield(1, amp_Domm, tor_Domm, pol_Domm, Xp[0], Xp[1], varphi + dvarphi);
			//BX = 0.0; BY = 0.0; BZ = 0.0;
			//BXpdR = 0.0; BYpdR = 0.0; BZpdR = 0.0;
			//BXpdZ = 0.0; BYpdZ = 0.0; BZpdZ = 0.0;
			//dBXdR = 0.0; dBYdR = 0.0; dBZdR = 0.0;
			//dBXdZ = 0.0; dBYdZ = 0.0; dBZdZ = 0.0;
			//BXpdR = 0.0; BYpdR = 0.0; BZpdR = 0.0;
			dmagfielddR = addstructsfield(1/delta, magfieldpdR, -1/delta, magfield);
			dmagfielddZ = addstructsfield(1/delta, magfieldpdZ, -1/delta, magfield);
			dmagfielddvarphi = addstructsfield(1/dvarphi, magfieldpdvarphi, -1/dvarphi, magfield);
			printf("dBRdR   %f %f\n", dmagfielddR.value[0], magfield->derivative[0][0]);
			printf("dBZdR   %f %f\n", dmagfielddR.value[1], magfield->derivative[1][0]);
			printf("dBphidR %f %f\n", dmagfielddR.value[2], magfield->derivative[2][0]);
			printf("dBRdZ   %f %f\n", dmagfielddZ.value[0], magfield->derivative[0][1]);
			printf("dBZdZ   %f %f\n", dmagfielddZ.value[1], magfield->derivative[1][1]);
			printf("dBphidZ %f %f\n", dmagfielddZ.value[2], magfield->derivative[2][1]);
			//printf("%f %f\n", magfield->derivative[1][0], dBZdR);
			//printf("%f %f\n", magfield->derivative[2][0], dBvarphidR);
			//printf("%f %f\n", magfield->derivative[0][1], dBRdZ);
			//printf("%f %f\n", magfield->derivative[1][1], dBZdZ);
			//printf("%f %f\n", magfield->derivative[2][1], dBvarphidZ);
			//divB = ( dmagfielddR.value[0] + magfield->value[0] / Xp[0] + dmagfielddZ.value[1] + dmagfielddvarphi.value[2] / Xp[0] ) / sqrt( pow(magfield->value[0], 2.0) + pow(magfield->value[1], 2.0) + pow(magfield->value[2], 2.0) );
			divB = ( dmagfielddR.value[0] + magfield->value[0] / Xp[0] + dmagfielddZ.value[1] + dmagfielddvarphi.value[2] / Xp[0] ) / sqrt( pow( dmagfielddR.value[0] + magfield->value[0] / Xp[0] , 2.0 ) + pow( dmagfielddZ.value[1], 2.0 ) + pow( dmagfielddvarphi.value[2] / Xp[0] , 2.0 ) );
			printf("div(B) = %f\n", divB);
		}
	}

	if (check==1) { 
	free(magfieldpdR); free(magfieldpdZ);
	}
	return magfield;
}

struct field *gradBfield(double *Xp, double varphi, double ***coils, int *num_coils, int **num_segs) {
	struct field *gradmagfield = calloc(1,sizeof(struct field));
	int m0_symmetry = 1;
	double epsilon[1], iota[2]; 
	int k_theta[1];
	epsilon[0] = 0.1; 
	//epsilon[0] = 0.0; epsilon[1] = 0.0;
	k_theta[0] = 6; 
	//epsilon[0] = 0.00; epsilon[1] = 0.00;
	iota[0] = 0.15; iota[1] = 0.38;
	gradmagfield = gradBReim(m0_symmetry, iota[0], iota[1], epsilon, k_theta, 1, Xp[0], Xp[1], varphi);
	return gradmagfield;
}
