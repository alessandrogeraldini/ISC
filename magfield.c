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
			XXvect_coil[0][coil_index][coilseg_index] = *storevals; XXvect_coil[1][coil_index][coilseg_index] = *(storevals+1);
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
		
struct field *BReim(int m0_symmetry, double iota0, double iota1, double *epsilon, int *k_theta, int size_epsilon, double RR, double ZZ, double varphi) {
		int ind;
		double R_axis = 1.0, theta, rmin, combo, combo1, dcombodR, dcombodZ, dcombo1dR, dcombo1dZ; // check_q;
		struct field *mag = calloc(1,sizeof(struct field));
		theta = atan2(ZZ, RR - R_axis);
		rmin = sqrt(pow((RR-R_axis), 2.0) + pow(ZZ, 2.0));
		combo = iota0 + iota1*rmin*rmin;
		dcombodR = 2.0*iota1*(RR - R_axis);
		dcombodZ = 2.0*iota1*ZZ;
		combo1 = 0.0;
		dcombo1dR = 0.0;
		dcombo1dZ = 0.0;
		for (ind=0; ind < size_epsilon; ind++) {
			combo  -= k_theta[ind] * epsilon[ind] * pow(rmin, k_theta[ind] - 2) * cos(k_theta[ind]*theta - m0_symmetry*varphi);
			combo1 += k_theta[ind] * epsilon[ind] * pow(rmin, k_theta[ind] - 2) * sin(k_theta[ind]*theta - m0_symmetry*varphi);
		}
		//epsilon[0] += 0.001; 
		//this was a test to verify that the tangent map determines the width
		// the small change of the position of the fixed point has no effect
		for (ind=0; ind < size_epsilon; ind++) {
			dcombodR -= k_theta[ind] * pow( rmin, k_theta[ind] - 4 ) * epsilon[ind] * ( k_theta[ind] * ZZ * sin( k_theta[ind]*theta - m0_symmetry*varphi )  +  (k_theta[ind] - 2) *  (RR - R_axis) * cos(k_theta[ind]*theta - m0_symmetry*varphi) ) ;
			//check_q = ( k_theta[ind] * pow( rmin, k_theta[ind] - 4 ) * epsilon[ind] * ( k_theta[ind] * pow( rmin, 2.0 ) * sin( k_theta[ind]*theta - m0_symmetry*varphi ) * ( - ZZ / pow(rmin, 2.0) ) - (k_theta[ind] - 2)  *  (RR - R_axis) * cos(k_theta[ind]*theta - m0_symmetry * varphi) ) );
			dcombodZ += pow( rmin, k_theta[ind] - 4 ) * epsilon[ind] * k_theta[ind] * ( k_theta[ind] * sin( k_theta[ind]*theta - m0_symmetry*varphi ) * ( RR - R_axis ) -  (k_theta[ind] - 2) * ZZ * cos(k_theta[ind]*theta - m0_symmetry * varphi) );
			dcombo1dR += k_theta[ind]*pow(rmin, k_theta[ind] - 4) * epsilon[ind] * ( - k_theta[ind]* ZZ * cos(k_theta[ind]*theta - m0_symmetry*varphi) + ( k_theta[ind] -2 ) * sin(k_theta[ind]*theta - m0_symmetry*varphi) * (RR - R_axis) ) ;
			dcombo1dZ += k_theta[ind] * pow(rmin, k_theta[ind] - 4) * epsilon[ind] * ( k_theta[ind] * cos(k_theta[ind]*theta - m0_symmetry *varphi) * (RR - R_axis) +  (k_theta[ind] - 2) * sin(k_theta[ind]*theta - m0_symmetry*varphi) * ZZ ) ;
		}
		mag->value[0] = ( ZZ / RR ) *  combo + ( ( RR - R_axis ) / RR ) * combo1 ;
		mag->value[1] = - ( (RR - R_axis) / RR ) * combo + ( ZZ / RR ) * combo1  ; 
		mag->value[2] = -1.0;
		mag->derivative[0][0] = ( - ZZ / pow( RR, 2.0) ) *  combo + (ZZ / RR) * dcombodR + combo1 * R_axis / pow(RR, 2.0) + dcombo1dR * (RR - R_axis) / RR;
		mag->derivative[1][0] = ( - R_axis / pow( RR, 2.0) ) *  combo - ((RR - R_axis) / RR) * dcombodR - combo1 * ZZ / pow(RR, 2.0) + dcombo1dR * ZZ / RR;
		mag->derivative[2][0] = 0.0;
		mag->derivative[0][1] = ( 1.0 / RR ) * combo + (ZZ / RR) * dcombodZ + dcombo1dZ * (RR - R_axis) / RR;
		mag->derivative[1][1] = - ((RR - R_axis) / RR) * dcombodZ + combo1 * ( 1.0 / RR ) + dcombo1dZ * ZZ / RR ;
		mag->derivative[2][1] = 0.0;
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
		}
		//mag[0].value[0] = (ZZ / RR) * dcomboiota; mag[0].value[1] = - ( (RR - R_axis) / RR ) * dcomboiota; mag[0].value[2] = 0.0;
		//mag[0].value[0] = (ZZ / RR) * dcomboiota1; mag[0].value[1] = - ( (RR - R_axis) / RR ) * dcomboiota1; mag[0].value[2] = 0.0;
		return mag;
}

struct field *Bfield(double *Xp, double varphi, double ***coils, int *num_coils, int **num_segs)
{
	//clock_t start = clock();
	double BZ=0.0, BX=0.0, BY=0.0, BR=0.0, Bvarphi=0.0, dBXdR=0.0, dBYdR=0.0, dBZdR=0.0, dBXdZ=0.0, dBYdZ=0.0, dBZdZ=0.0, dBRdR=0.0, dBvarphidR=0.0, dBRdZ=0.0, dBvarphidZ=0.0;
	double Xminxc=0.0, Yminyc=0.0, Zminzc=0.0, Rminrc=0.0, dl=0.0, dxc=0.0, dyc=0.0, dzc=0.0;
	double XX=Xp[0]*cos(varphi), YY=Xp[0]*sin(varphi), ZZ = Xp[1];
	double divB;
	struct field *magfield = calloc(1,sizeof(struct field)), *magfielddR = calloc(1,sizeof(struct field)), *magfielddZ = calloc(1,sizeof(struct field)), *magfielddvarphi = calloc(1, sizeof(struct field));
	struct field dmagfielddR, dmagfielddZ, dmagfielddvarphi;
	int coil_index=0, coilseg_index=0; 
	int check_coil = 0, check = 0, m0_symmetry = 1;
	double XpdR, YpdR, RpdRminrc, XpdRminxc, YpdRminyc, BXpdR, BYpdR, BZpdR, dR = 0.00001, dvarphi = 0.00001;
	double RpdZminrc, ZpdZminzc, BXpdZ, BYpdZ, BZpdZ, dZ = 0.00001;
	double amp_Domm[1], epsilon[1], iota[2]; 
	int pol_Domm[1], tor_Domm[1], k_theta[2];
	char *type = "Reim";
	//amp_Domm[0]=1.73; pol_Domm[0] = 2; tor_Domm[0] = 5;
	amp_Domm[0]=2.0; pol_Domm[0] = 5; tor_Domm[0] = 5;
	//amp_Domm[0]=0.001; pol_Domm[0] = 4; tor_Domm[0] = 8;
	//amp_Domm[0]=0.00000; pol_Domm[0] = 2; tor_Domm[0] = 5;
	//printf("n_segs[1]=%d\n", n_segs[1]);
	//printf("Sanity check: print line of coil file\n");
	//printf("coils[0][%d][%d]=%f\n", j, k, coils[0][j][k]);
	//printf("coils[1][%d][%d]=%f\n", j, k, coils[1][j][k]);
	//printf("coils[2][%d][%d]=%f\n", j, k, coils[2][j][k]);
	//printf("coils[3][%d][%d]=%f\n", j, k, coils[3][j][k]);
	//printf("Sanity check: print another line of coil file\n");
	//printf("coils[0][0][1]=%Le\n", coils[0][0][1]);
	//printf("coils[1][0][1]=%Le\n", coils[1][0][1]);
	//printf("coils[2][0][1]=%Le\n", coils[2][0][1]);
	//printf("coils[3][0][1]=%Le\n", coils[3][0][1]);
	if (strncmp(type, "coil", 4) == 0) {
		for (coil_index=0;coil_index<*num_coils;coil_index++)
		{
			//printf("num_coils = %d\n", *num_coils);
			//printf("coil_index = %d\n", coil_index);
			for (coilseg_index=0;coilseg_index<*(*num_segs+coil_index)-1;coilseg_index++)
			{
				//printf("num_segs[%d] = %d\n", coil_index, *(*num_segs+coil_index));
				//printf("coilseg_index = %d\n", coilseg_index);
				Xminxc = XX - 0.5*coils[0][coil_index][coilseg_index] - 0.5*coils[0][coil_index][coilseg_index+1];
				Yminyc = YY - 0.5*coils[1][coil_index][coilseg_index] - 0.5*coils[1][coil_index][coilseg_index+1]; 
				Zminzc = ZZ - 0.5*coils[2][coil_index][coilseg_index] - 0.5*coils[2][coil_index][coilseg_index+1]; 
				Rminrc = sqrt(pow(Xminxc,  2.0) + pow(Yminyc, 2.0) + pow(Zminzc, 2.0)); 

				dxc = coils[0][coil_index][coilseg_index+1] - coils[0][coil_index][coilseg_index];
				dyc = coils[1][coil_index][coilseg_index+1] - coils[1][coil_index][coilseg_index];
				dzc = coils[2][coil_index][coilseg_index+1] - coils[2][coil_index][coilseg_index];
				dl = sqrt(pow(dxc, 2.0) + pow(dyc, 2.0) + pow(dzc, 2.0));

				BX += (coils[3][coil_index][coilseg_index]*(dyc*Zminzc - dzc*Yminyc)/pow(Rminrc, 3.0));
				BY += (coils[3][coil_index][coilseg_index]*(dzc*Xminxc - dxc*Zminzc)/pow(Rminrc, 3.0));
				BZ += (coils[3][coil_index][coilseg_index]*(dxc*Yminyc - dyc*Xminxc)/pow(Rminrc, 3.0));
				
				dBXdR += (coils[3][coil_index][coilseg_index]*(dyc*0.0         - dzc*sin(varphi))/pow(Rminrc, 3.0));
				dBYdR += (coils[3][coil_index][coilseg_index]*(dzc*cos(varphi) - dxc*0.0        )/pow(Rminrc, 3.0));
				dBZdR += (coils[3][coil_index][coilseg_index]*(dxc*sin(varphi) - dyc*cos(varphi))/pow(Rminrc, 3.0));

				dBXdR -= 3*(coils[3][coil_index][coilseg_index]*(Xminxc*cos(varphi) + Yminyc*sin(varphi))*(dyc*Zminzc - dzc*Yminyc)/pow(Rminrc, 5.0));
				dBYdR -= 3*(coils[3][coil_index][coilseg_index]*(Xminxc*cos(varphi) + Yminyc*sin(varphi))*(dzc*Xminxc - dxc*Zminzc)/pow(Rminrc, 5.0));
				dBZdR -= 3*(coils[3][coil_index][coilseg_index]*(Xminxc*cos(varphi) + Yminyc*sin(varphi))*(dxc*Yminyc - dyc*Xminxc)/pow(Rminrc, 5.0));

				dBXdZ += (coils[3][coil_index][coilseg_index]*( dyc)/pow(Rminrc, 3.0));
				dBYdZ += (coils[3][coil_index][coilseg_index]*(-dxc)/pow(Rminrc, 3.0));
				//dBZdZ += (coils[3][coil_index][coilseg_index]*( 0.0)/pow(Rminrc, 3.0));

				dBXdZ -= 3*(coils[3][coil_index][coilseg_index]*(Zminzc)*(dyc*Zminzc - dzc*Yminyc)/pow(Rminrc, 5.0));
				dBYdZ -= 3*(coils[3][coil_index][coilseg_index]*(Zminzc)*(dzc*Xminxc - dxc*Zminzc)/pow(Rminrc, 5.0));
				dBZdZ -= 3*(coils[3][coil_index][coilseg_index]*(Zminzc)*(dxc*Yminyc - dyc*Xminxc)/pow(Rminrc, 5.0));
			}
		}

		//magfield = calloc(9,sizeof(double)); // first triplet B; if derivatives[0] ==1: second triplet dB/dR, third triplet dB/dZ 
		BX*=pow(10.0,-7.0); BY*=pow(10.0,-7.0); BZ*=pow(10.0,-7.0); // multiply by mu_0/4*PI
		BR = cos(varphi)*BX + sin(varphi)*BY; Bvarphi = -sin(varphi)*BX + cos(varphi)*BY; // BZ = BZ;
		magfield->value[0] = BR; magfield->value[1] = BZ; magfield->value[2] = Bvarphi;
		dBXdR*=pow(10.0,-7.0); dBYdR*=pow(10.0,-7.0); dBZdR*=pow(10.0,-7.0); // multiply by mu_0/4*PI
		dBXdZ*=pow(10.0,-7.0); dBYdZ*=pow(10.0,-7.0); dBZdZ*=pow(10.0,-7.0); // multiply by mu_0/4*PI
		dBRdR = dBXdR*cos(varphi) + dBYdR*sin(varphi); dBvarphidR = -dBXdR*sin(varphi) + dBYdR*cos(varphi);
		dBRdZ = dBXdZ*cos(varphi) + dBYdZ*sin(varphi); dBvarphidZ = -dBXdZ*sin(varphi) + dBYdZ*cos(varphi);
		//printf("dBXdR=%f attempt 2\n", dBXdR);
		magfield->derivative[0][0] = dBRdR; 
		magfield->derivative[1][0] = dBZdR; 
		magfield->derivative[2][0] = dBvarphidR;
		magfield->derivative[0][1] = dBRdZ; 
		magfield->derivative[1][1] = dBZdZ; 
		magfield->derivative[2][1] = dBvarphidZ;
	}
	else if (strncmp(type, "Reim", 4) == 0) {
		epsilon[0] = 0.005; //epsilon[1] = 0.0;
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
		dBRdR = 0.0;
		dBZdR = 0.0;
		dBvarphidR = 0.0;
		dBRdZ = 0.0;
		dBZdZ = 0.0;
		dBvarphidZ = 0.0;
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
			magfielddR      = BReim(m0_symmetry, iota[0], iota[1], epsilon, k_theta, 2, Xp[0]+dR, Xp[1], varphi);
			magfielddZ      = BReim(m0_symmetry, iota[0], iota[1], epsilon, k_theta, 2, Xp[0], Xp[1]+dZ, varphi);
			magfielddvarphi = BReim(m0_symmetry, iota[0], iota[1], epsilon, k_theta, 2, Xp[0], Xp[1], varphi + dvarphi);
			BX = 0.0; BY = 0.0; BZ = 0.0;
			BXpdR = 0.0; BYpdR = 0.0; BZpdR = 0.0;
			BXpdZ = 0.0; BYpdZ = 0.0; BZpdZ = 0.0;
			dBXdR = 0.0; dBYdR = 0.0; dBZdR = 0.0;
			dBXdZ = 0.0; dBYdZ = 0.0; dBZdZ = 0.0;
			BXpdR = 0.0; BYpdR = 0.0; BZpdR = 0.0;
			dmagfielddR = addstructsfield(1/dR, magfielddR, -1/dR, magfield);
			dmagfielddZ = addstructsfield(1/dZ, magfielddZ, -1/dZ, magfield);
			dmagfielddvarphi = addstructsfield(1/dvarphi, magfielddvarphi, -1/dvarphi, magfield);
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
		magfield = DommBfield(1, amp_Domm, tor_Domm, pol_Domm, Xp[0], Xp[1], varphi);
		if (check==1) {
			magfielddR = DommBfield(1, amp_Domm, tor_Domm, pol_Domm, Xp[0]+dR, Xp[1], varphi);
			magfielddZ = DommBfield(1, amp_Domm, tor_Domm, pol_Domm, Xp[0], Xp[1]+dZ, varphi);
			magfielddvarphi = DommBfield(1, amp_Domm, tor_Domm, pol_Domm, Xp[0], Xp[1], varphi + dvarphi);
			BX = 0.0; BY = 0.0; BZ = 0.0;
			BXpdR = 0.0; BYpdR = 0.0; BZpdR = 0.0;
			BXpdZ = 0.0; BYpdZ = 0.0; BZpdZ = 0.0;
			dBXdR = 0.0; dBYdR = 0.0; dBZdR = 0.0;
			dBXdZ = 0.0; dBYdZ = 0.0; dBZdZ = 0.0;
			BXpdR = 0.0; BYpdR = 0.0; BZpdR = 0.0;
			dmagfielddR = addstructsfield(1/dR, magfielddR, -1/dR, magfield);
			dmagfielddZ = addstructsfield(1/dZ, magfielddZ, -1/dZ, magfield);
			dmagfielddvarphi = addstructsfield(1/dvarphi, magfielddvarphi, -1/dvarphi, magfield);
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

	if (check_coil==1) {
		BX = 0.0; BY = 0.0; BZ = 0.0;
		BXpdR = 0.0; BYpdR = 0.0, BZpdR = 0.0;
		BXpdZ = 0.0; BYpdZ = 0.0, BZpdZ = 0.0;
		dBXdR = 0.0; dBYdR = 0.0; dBZdR = 0.0;
		dBXdZ = 0.0; dBYdZ = 0.0; dBZdZ = 0.0;
		BXpdR = 0.0; BYpdR = 0.0; BZpdR = 0.0;
		for (coil_index=0;coil_index<*num_coils;coil_index++) {
			for (coilseg_index=0;coilseg_index<*(*num_segs+coil_index)-1;coilseg_index++) {
				Xminxc = XX - 0.5*coils[0][coil_index][coilseg_index] - 0.5*coils[0][coil_index][coilseg_index+1];
				Yminyc = YY - 0.5*coils[1][coil_index][coilseg_index] - 0.5*coils[1][coil_index][coilseg_index+1]; 
				Zminzc = ZZ - 0.5*coils[2][coil_index][coilseg_index] - 0.5*coils[2][coil_index][coilseg_index+1]; 
				Rminrc = sqrt(pow(Xminxc,  2.0) + pow(Yminyc, 2.0) + pow(Zminzc, 2.0)); 
				
				dxc = coils[0][coil_index][coilseg_index+1] - coils[0][coil_index][coilseg_index];
				dyc = coils[1][coil_index][coilseg_index+1] - coils[1][coil_index][coilseg_index];
				dzc = coils[2][coil_index][coilseg_index+1] - coils[2][coil_index][coilseg_index];
				dl = sqrt(pow(dxc, 2.0) + pow(dyc, 2.0) + pow(dzc, 2.0));

				BX += (coils[3][coil_index][coilseg_index]*(dyc*Zminzc - dzc*Yminyc)/pow(Rminrc, 3.0));
				BY += (coils[3][coil_index][coilseg_index]*(dzc*Xminxc - dxc*Zminzc)/pow(Rminrc, 3.0));
				BZ += (coils[3][coil_index][coilseg_index]*(dxc*Yminyc - dyc*Xminxc)/pow(Rminrc, 3.0));

				XpdR = (Xp[0]+dR)*cos(varphi); YpdR = (Xp[0]+dR)*sin(varphi); 
				XpdRminxc = XpdR - 0.5*coils[0][coil_index][coilseg_index] - 0.5*coils[0][coil_index][coilseg_index+1];
				YpdRminyc = YpdR - 0.5*coils[1][coil_index][coilseg_index] - 0.5*coils[1][coil_index][coilseg_index+1]; 
				ZpdZminzc = Xp[1] + dZ - 0.5*coils[2][coil_index][coilseg_index] - 0.5*coils[2][coil_index][coilseg_index+1]; 
				RpdRminrc = sqrt(pow(XpdRminxc,  2.0) + pow(YpdRminyc, 2.0) + pow(Zminzc, 2.0)); 
				RpdZminrc = sqrt(pow(Xminxc,  2.0) + pow(Yminyc, 2.0) + pow(ZpdZminzc, 2.0)); 

				BXpdR += (coils[3][coil_index][coilseg_index]*(dyc*Zminzc - dzc*YpdRminyc)/pow(RpdRminrc, 3.0));
				BYpdR += (coils[3][coil_index][coilseg_index]*(dzc*XpdRminxc - dxc*Zminzc)/pow(RpdRminrc, 3.0));
				BZpdR += (coils[3][coil_index][coilseg_index]*(dxc*YpdRminyc - dyc*XpdRminxc)/pow(RpdRminrc, 3.0));
				BXpdZ += (coils[3][coil_index][coilseg_index]*(dyc*ZpdZminzc - dzc*Yminyc)/pow(RpdZminrc, 3.0));
				BYpdZ += (coils[3][coil_index][coilseg_index]*(dzc*Xminxc - dxc*ZpdZminzc)/pow(RpdZminrc, 3.0));
				BZpdZ += (coils[3][coil_index][coilseg_index]*(dxc*Yminyc - dyc*Xminxc)/pow(RpdZminrc, 3.0));
				}
			}
		dBXdR = (BXpdR - BX)/dR;
		dBYdR = (BYpdR - BY)/dR;
		dBZdR = (BZpdR - BZ)/dR;
		dBXdZ = (BXpdZ - BX)/dZ;
		dBYdZ = (BYpdZ - BY)/dZ;
		dBZdZ = (BZpdZ - BZ)/dZ;
		dBXdR*=pow(10.0,-7.0); dBYdR*=pow(10.0,-7.0); dBZdR*=pow(10.0,-7.0); // multiply by mu_0/4*PI
		dBXdZ*=pow(10.0,-7.0); dBYdZ*=pow(10.0,-7.0); dBZdZ*=pow(10.0,-7.0); // multiply by mu_0/4*PI
		dBRdR = dBXdR*cos(varphi) + dBYdR*sin(varphi); dBvarphidR = -dBXdR*sin(varphi) + dBYdR*cos(varphi);
		dBRdZ = dBXdZ*cos(varphi) + dBYdZ*sin(varphi); dBvarphidZ = -dBXdZ*sin(varphi) + dBYdZ*cos(varphi);

		//printf("%f %f\n", magfield->derivative[0][0], dBRdR);
		//printf("%f %f\n", magfield->derivative[1][0], dBZdR);
		//printf("%f %f\n", magfield->derivative[2][0], dBvarphidR);
		//printf("%f %f\n", magfield->derivative[0][1], dBRdZ);
		//printf("%f %f\n", magfield->derivative[1][1], dBZdZ);
		//printf("%f %f\n", magfield->derivative[2][1], dBvarphidZ);
	}
	//clock_t int3 = clock();
	//printf("Time after evaluating B field at a point: %f\n", (double) (int3-start)/CLOCKS_PER_SEC);
	free(magfielddR); free(magfielddZ);
	return magfield;
}

struct field *gradBfield(double *Xp, double varphi, double ***coils, int *num_coils, int **num_segs) {
	struct field *gradmagfield = calloc(1,sizeof(struct field));
	int m0_symmetry = 1;
	double epsilon[1], iota[2]; 
	int k_theta[2];
	epsilon[0] = 0.005; 
	//epsilon[0] = 0.0; epsilon[1] = 0.0;
	k_theta[0] = 6; 
	//epsilon[0] = 0.00; epsilon[1] = 0.00;
	iota[0] = 0.15; iota[1] = 0.38;
	gradmagfield = gradBReim(m0_symmetry, iota[0], iota[1], epsilon, k_theta, 1, Xp[0], Xp[1], varphi);
	return gradmagfield;
}
