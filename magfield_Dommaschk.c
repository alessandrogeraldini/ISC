//Author: Alessandro Geraldini
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_sf_ellint.h>
#include "isc.h"

#define aloop 1.0 
//#define Itor 0.0158*1000.0
//#define Itor 0.0002*1000.0
#define Itor 0.00008*1000.0
#define IZ 1000.0
#define mu0 4*M_PI*pow(10.0, -7.0)

double BRloop(double R, double Z) {
	double alphasq, betasq, beta, ksq, kk;
	double ellipe, ellipk, y, C;
	C = mu0*Itor/M_PI;
	alphasq = aloop*aloop + R*R + Z*Z - 2*aloop*R;
	betasq = aloop*aloop + R*R + Z*Z + 2*aloop*R;
	beta = sqrt(betasq);
	ksq = 1.0 - alphasq / betasq;
	kk = sqrt(ksq);
	//printf("R = %f, Z = %f\n", R, Z);
	//printf("is it here (before1)? k = %f\n", kk);
	ellipk = gsl_sf_ellint_Kcomp(kk, GSL_PREC_DOUBLE);
	ellipe = gsl_sf_ellint_Ecomp(kk, GSL_PREC_DOUBLE);
	y = (C*Z/(2.0*alphasq*beta*R))*((aloop*aloop + R*R+ Z*Z)*ellipe- alphasq*ellipk);
	//printf("k = %f\n", kk);
	return y;
}

double BZloop(double R, double Z) {
	double alphasq, betasq, beta, ksq, kk;
	double ellipe, ellipk, y, C;
	C = mu0*Itor/M_PI;
	alphasq = aloop*aloop + R*R + Z*Z - 2*aloop*R;
	betasq = aloop*aloop + R*R + Z*Z + 2*aloop*R;
	beta = sqrt(betasq);
	ksq = 1.0 - alphasq / betasq;
	kk = sqrt(ksq);
	//printf("is it here (before2)? k = %f\n", kk);
	ellipk = gsl_sf_ellint_Kcomp(kk, GSL_PREC_DOUBLE);
	ellipe = gsl_sf_ellint_Ecomp(kk, GSL_PREC_DOUBLE);
	y = (C/(2.0*alphasq*beta))*((aloop*aloop - R*R - Z*Z)*ellipe + alphasq*ellipk);
	//printf("k = %f\n", kk);
	return y;
}

double dRBRloop(double R, double Z) {
	double alphasq, betasq, beta, ksq, kk;
	double ellipe, ellipk, y, C;
	C = mu0*Itor/M_PI;
	alphasq = aloop*aloop + R*R + Z*Z - 2*aloop*R;
	betasq = aloop*aloop + R*R + Z*Z + 2*aloop*R;
	beta = sqrt(betasq);
	ksq = 1.0 - alphasq / betasq;
	kk = sqrt(ksq);
	//printf("is it here (before3)? k = %f\n", kk);
	ellipk = gsl_sf_ellint_Kcomp(kk, GSL_PREC_DOUBLE);
	ellipe = gsl_sf_ellint_Ecomp(kk, GSL_PREC_DOUBLE);
	y = -C*Z/(2.0*R*R*(pow(alphasq, 2.0)*pow(beta, 3.0)))*((pow(aloop,6.0) + pow((R*R + Z*Z),2.0)*(2*R*R +Z*Z) + pow(aloop,4.0)*(3*Z*Z-8*R*R) + aloop*aloop*(5.0*pow(R,4.0) -4*R*R*Z*Z + 3*pow(Z,4.0)))*ellipe - alphasq*(pow(aloop,4.0) - 3*aloop*aloop*R*R + 2*pow(R,4.0) + (2*aloop*aloop + 3*R*R)*Z*Z + pow(Z,4.0))*ellipk);
	//printf("k = %f\n", kk);
	return y;
}

double dZBRloop(double R, double Z) {
	double alphasq, betasq, beta, ksq, kk;
	double ellipe, ellipk, y;
	double C = mu0*Itor/M_PI;
	alphasq = aloop*aloop + R*R + Z*Z - 2*aloop*R;
	betasq = aloop*aloop + R*R + Z*Z + 2*aloop*R;
	beta = sqrt(betasq);
	ksq = 1.0 - alphasq / betasq;
	kk = sqrt(ksq);
	//printf("is it here (before4)? k = %f\n", kk);
	ellipk = gsl_sf_ellint_Kcomp(kk, GSL_PREC_DOUBLE);
	ellipe = gsl_sf_ellint_Ecomp(kk, GSL_PREC_DOUBLE);
	//y = C/(2.0*R*alphasq*alphasq*pow(beta,3.0))*(((aloop*aloop+R*R)*(pow(Z, 4.0) + pow((aloop*aloop - R*R),2.0)) + 2*Z*Z*(pow(aloop,4.0) - 6*aloop*aloop*R*R +R*R))*ellipe - alphasq*(pow((aloop*aloop - R*R),2.0) + (aloop*aloop+R*R)*Z*Z)*ellipk);
	y = C/(2.0*R*alphasq*alphasq*pow(beta,3.0))*(((aloop*aloop+R*R)*(pow(Z, 4.0) + pow((aloop*aloop - R*R),2.0)) + 2*Z*Z*(pow(aloop,4.0) - 6*aloop*aloop*R*R + pow(R,4.0)))*ellipe - alphasq*(pow((aloop*aloop - R*R),2.0) + (aloop*aloop+R*R)*Z*Z)*ellipk);
	//printf("k = %f\n", kk);
	return y;
}

double  dZBZloop(double R, double Z) {
	double alphasq, betasq, beta, ksq, kk;
	double ellipe, ellipk, y;
	double C = mu0*Itor/M_PI;
	alphasq = aloop*aloop + R*R + Z*Z - 2*aloop*R;
	betasq =  aloop*aloop + R*R + Z*Z + 2*aloop*R;
	beta = sqrt(betasq);
	//printf("alphasq=%f\n", alphasq);
	//printf("betasq=%f\n", betasq);
	ksq = 1.0 - alphasq / betasq;
	kk = sqrt(ksq);
	//printf("is it here (before5)? k = %f\n", kk);
	ellipk = gsl_sf_ellint_Kcomp(kk, GSL_PREC_DOUBLE);
	ellipe = gsl_sf_ellint_Ecomp(kk, GSL_PREC_DOUBLE);
	y = C*Z/(2.0*alphasq*alphasq*pow(beta,3.0))*((6*aloop*aloop*(R*R-Z*Z) - 7*aloop*aloop + pow((R*R + Z*Z),2.0))*ellipe + alphasq*(aloop*aloop - R*R - Z*Z)*ellipk);
	//printf("k = %f\n", kk);
	return y;
}

double dRBZloop(double R, double Z) {
	double alphasq, betasq, beta, ksq, kk;
	double ellipe, ellipk, y;
	double C = mu0*Itor/M_PI;
	alphasq = aloop*aloop + R*R + Z*Z - 2*aloop*R;
	betasq =  aloop*aloop + R*R + Z*Z + 2*aloop*R;
	beta = sqrt(betasq);
	ksq = 1.0 - alphasq / betasq;
	kk = sqrt(ksq);
	//printf("is it here (before6)? k = %f\n", kk);
	ellipk = gsl_sf_ellint_Kcomp(kk, GSL_PREC_DOUBLE);
	ellipe = gsl_sf_ellint_Ecomp(kk, GSL_PREC_DOUBLE);
	y = C/(2.0*R*alphasq*alphasq*pow(beta,3.0))*(((aloop*aloop+R*R)*(pow(Z, 4.0) + pow((aloop*aloop - R*R),2.0)) + 2*Z*Z*(pow(aloop,4.0) - 6*aloop*aloop*R*R + pow(R,4.0)))*ellipe - alphasq*(pow((aloop*aloop - R*R),2.0) + (aloop*aloop+R*R)*Z*Z)*ellipk);
	//printf("is it here? k = %f\n", kk);
	return y;
}

//int main() {
//	double BR, BZ;
//	BR = BRloop(0.0000, 0.0);
//	BZ = BZloop(0.0000, 0.0);
//	printf("BR = %f\nBZ = %f\n", BR, BZ);
//	return 0;
//}

int fac(int number) {
	int ind, result=1;
	for (ind=0;ind<number;ind++) {
		result = result*(ind+1);
		//printf("IN fac: output = %d\n, ind = %d", result, ind);
	}	
	//printf("IN fac: output = %d\n", result);
	return result;
}

double alpha(int m,int l) {
	double y;
	y = pow(-1.0,l)/(fac(m+l)*fac(l)*pow(2.0,2*l+m));
	//printf("alpha(%d, %d) = %.9f\n", m, l, y); //checked
	return y;
}

double beta(int m,int l) {
	double y;
	y = fac(m-l-1)/(fac(l)*pow(2.0,2*l-m+1));
	//printf("beta(%d, %d) = %.9f\n", m, l, y); //checked
	return y;
}

double Dmn(int m, int n, double R, double Z) {
	double sumD=0.0, y=0.0;
	int j, k;
	for (k=0; k< n/2 + 1; k++) {
		sumD = 0.0;
		for (j=0;j<k+1;j++) {
			//printf(alpha(m,j), beta(m,k-j), alpha(m, k-j), beta(m,j))
			//printf("number is %f\n", R);
			sumD += ( - (2*k-2*j - m)*alpha(m,j)*beta(m,k-j)*pow(R,(2*j+m)) + (2*k-2*j + m)*alpha(m,k-j)*beta(m,j)*pow(R,2*j-m) );
			//printf("sumD is %f\n", sumD);
		}

		y += (pow(Z,n-2*k)/fac(n-2*k))*sumD;
		//printf("Dmn is %f\n", y);
	}
	//printf("IN Dmn: output = %f\n", y);
	//printf("Dmn(%d, %d, %f, %f) = %.10f\n", m, n, R, Z, y);
	return y;
}

double Nmn(int m, int n, double R, double Z) {
	double sumN = 0.0, y = 0.0;
	int j, k;
	for (k=0;k<n/2 + 1; k++) {
		sumN = 0.0;
		for (j=0;j<k+1;j++) {
			//print(alpha(m,j), beta(m,k-j), alpha(m, k-j), beta(m,j))
			sumN += ( alpha(m,j)*beta(m,k-j)*pow(R,2*j+m) - alpha(m,k-j)*beta(m,j)*pow(R,2*j-m) );
		}
		y += (pow(Z,n-2*k)/fac(n-2*k))*sumN;
	}
	//printf("IN Nmn: output = %f\n", y);
	//printf("Nmn is %f", y);
	//printf("Nmn(%d, %d, %f, %f) = %.10f\n", m, n, R, Z, y);
	return y;
}

double dRDmn(int m, int n, double R, double Z) {
	double sumD= 0.0, y=0.0, intval;
	int j, k;
	for (k=0; k<n/2 + 1; k++) {
		sumD = 0.0;
		for (j=0;j<k+1;j++) {
			intval = ( - (2*j+m)*(2*k-2*j - m)*alpha(m,j)*beta(m,k-j)*pow(R,2*j+m-1) + (2*j-m)*(2*k-2*j + m)*alpha(m,k-j)*beta(m,j)*pow(R,2*j-m-1) );
			sumD += intval;
			//printf("sumD=%f\n", sumD);
			//printf("intval=%f\n", intval);
		}
		y += (pow(Z,n-2*k)/fac(n-2*k))*sumD;
		//printf("y=%f\n", y);
	}
	//printf("dRDmn(%d, %d, %f, %f) = %.10f\n", m, n, R, Z, y);
	//printf("returned y=%f\n", y);
	return y;
}

double dZDmn(int m, int n, double R, double Z) {
	double sumD= 0.0, y=0.0;
	int j, k;
	for (k=0;k<(n-1)/2 + 1; k++) {
		sumD = 0.0;
		for (j=0; j<k+1;j++) {
			//printf("pow = %f\n", pow(R,2*j+m));
			sumD += ( - (2*k-2*j - m)*alpha(m,j)*beta(m,k-j)*pow(R,2*j+m) + (2*k-2*j + m)*alpha(m,k-j)*beta(m,j)*pow(R,2*j-m) );
		}
		y += (pow(Z,n-2*k-1)/fac(n-2*k-1))*sumD;
	}
	//printf("dZDmn(%d, %d, %f, %f) = %.10f\n", m, n, R, Z, y);
	return y;
}

double dRRDmn(int m, int n, double R, double Z) {
	double sumD= 0.0, y=0.0;
	int j, k;
	for (k=0;k<n/2 + 1;k++) {
		sumD = 0.0;
		for (j=0; j<k+1; j++) {
			sumD += ( - (2*j+m-1)*(2*j+m)*(2*k-2*j - m)*alpha(m,j)*beta(m,k-j)*pow(R,2*j+m-2) + (2*j-m-1)*(2*j-m)*(2*k-2*j + m)*alpha(m,k-j)*beta(m,j)*pow(R,2*j-m-2) );
		}
		y += (pow(Z,n-2*k)/fac(n-2*k))*sumD;
	}
	//printf("dRRDmn(%d, %d, %f, %f) = %.10f\n", m, n, R, Z, y);
	return y;
}

double dZZDmn(int m, int n, double R, double Z) {
	double sumD= 0.0, y=0.0;
	int j, k;
	for (k=0;k<n/2;k++) {
		sumD = 0.0;
		for (j=0; j<k+1; j++) {
			sumD += ( - (2*k-2*j - m)*alpha(m,j)*beta(m,k-j)*pow(R,2*j+m) + (2*k-2*j + m)*alpha(m,k-j)*beta(m,j)*pow(R,2*j-m) );
		}
		y += (pow(Z,n-2*k-2)/fac(n-2*k-2))*sumD;
	}
	//printf("dZZDmn(%d, %d, %f, %f) = %.10f\n", m, n, R, Z, y);
	return y;
}

double dRZDmn(int m, int n, double R, double Z) {
	double sumD= 0.0, y=0.0;
	int j, k;
	for (k=0;k<(n-1)/2 + 1;k++) {
		sumD = 0.0;
		for (j=0; j<k+1; j++) {
			sumD += ( - (2*j+m)*(2*k-2*j - m)*alpha(m,j)*beta(m,k-j)*pow(R,2*j+m) + (2*j-m)*(2*k-2*j + m)*alpha(m,k-j)*beta(m,j)*pow(R,2*j-m) );
		}
		y += (pow(Z,n-2*k-1)/fac(n-2*k-1))*sumD;
	}
	//printf("IN dRZDmn: output = %f\n", y);
	//printf("dRZDmn(%d, %d, %f, %f) = %.10f\n", m, n, R, Z, y);
	return y;
}

double dRNmn(int m, int n, double R, double Z) {
	double sumN= 0.0, y=0.0;
	int j, k;
	for (k=0;k<n/2 + 1; k++) {
		sumN = 0.0;
		for (j=0; j<k+1; j++) {
			sumN += ( (2*j+m)*alpha(m,j)*beta(m,k-j)*pow(R,2*j+m-1) - (2*j-m)*alpha(m,k-j)*beta(m,j)*pow(R,2*j-m-1) );
		}	
	y += (pow(Z,n-2*k)/fac(n-2*k))*sumN;
	}
	//printf("dRNmn(%d, %d, %f, %f) = %.10f\n", m, n, R, Z, y);
	//printf("IN dRNmn: output = %f\n", y);
	return y;
}

double dZNmn(int m, int n, double R, double Z) {
	double sumN= 0.0, y=0.0;
	int j, k;
	for (k=0;k<(n-1)/2 + 1;k++) {
		sumN = 0.0;
		for (j=0; j<k+1; j++) {
			sumN += ( alpha(m,j)*beta(m,k-j)*pow(R,2*j+m) - alpha(m,k-j)*beta(m,j)*pow(R,2*j-m) );
		}
		y += (pow(Z,n-2*k-1)/fac(n-2*k-1))*sumN;
	}
	//printf("dZNmn(%d, %d, %f, %f) = %.10f\n", m, n, R, Z, y);
	//printf("IN dZNmn: output = %f\n", y);
	return y;
}
		
double dRRNmn(int m, int n, double R, double Z) {
	double sumN= 0.0, y=0.0;
	int j, k;
	for (k=0; k< n/2 + 1; k++) {
		sumN = 0.0;
		for (j=0; j< k+1; j++) {
			sumN += ( (2*j+m-1)*(2*j+m)*alpha(m,j)*beta(m,k-j)*pow(R,2*j+m-2) - (2*j-m-1)*(2*j-m)*alpha(m,k-j)*beta(m,j)*pow(R,2*j-m-2) );
		}
		y += (pow(Z,n-2*k)/fac(n-2*k))*sumN;
	}
	//printf("dRRNmn(%d, %d, %f, %f) = %.10f\n", m, n, R, Z, y);
	//printf("IN dRRNmn: output = %f\n", y);
	return y;
}

double dZZNmn(int m, int n, double R, double Z) {
	double sumN= 0.0, y=0.0;
	int j, k;
	for (k=0; k< n/2; k++) {
		sumN = 0.0;
		for (j =0; j<k+1; j++) {
			sumN += ( alpha(m,j)*beta(m,k-j)*pow(R,2*j+m) - alpha(m,k-j)*beta(m,j)*pow(R,2*j-m) );
		}
		y += (pow(Z,n-2*k-2)/fac(n-2*k-2))*sumN;
	}
	//printf("IN dZZNmn: output = %f\n", y);
	//printf("dZZNmn(%d, %d, %f, %f) = %.10f\n", m, n, R, Z, y);
	return y;
}
		
double dRZNmn(int m, int n, double R, double Z) {
	double sumN= 0.0, y=0.0;
	int j, k;
	for (k=0;k<(n-1)/2 + 1; k++) {
		sumN = 0.0;
		for (j=0; j< k+1; j++) {
			sumN += ( (2*j+m)*alpha(m,j)*beta(m,k-j)*pow(R,2*j+m-1) - (2*j-m)*alpha(m,k-j)*beta(m,j)*pow(R,2*j-m-1) );
		}
		y += (pow(Z,n-2*k-1)/fac(n-2*k-1))*sumN;
	}
	//printf("dRZNmn(%d, %d, %f, %f) = %.10f\n", m, n, R, Z, y);
	return y;
}

double Phi(int m, int nn, double R, double Z, double phi) {
	double y;
	int a, b, c, d;
	if (nn%2 == 0) {
		a = d = 0;
		b = 1;
		c = 1;
	}
	else {
		a = 1; 
		d = -1; 
		b = c = 0;
	}
	y = (a*cos(m*phi) + b*sin(m*phi))*Dmn(m,nn,R,Z) + (c*cos(m*phi) + d*sin(m*phi))*Nmn(m,nn-1,R,Z);
	return y;
}

double BR(int m, int n, double R, double Z, double phi) {
	double y;
	int a, b, c, d;
	if (n%2 == 0) {
		a = d = 0;
		b = 1;
		c = 1;
	}
	else {
		a = 1; 
		d = -1; 
		b = c = 0;
	}
	y = (a*cos(m*phi) + b*sin(m*phi))*dRDmn(m,n,R,Z) + (c*cos(m*phi) + d*sin(m*phi))*dRNmn(m,n-1,R,Z);
	//printf("dRDmn = %f\n", dRDmn(m,n,R,Z));
	//printf("dRNmn = %f\n", dRNmn(m,n,R,Z));
	//printf("B = %f\n", y);
	//printf("BR(%d, %d, %f, %f, %f) = %.10f\n", m, n, R, Z, phi, y);
	return y;
}

double BZ(int m, int n, double R, double Z, double phi) {
	double y;
	int a, b, c, d;
	if (n%2 == 0) {
		a = d = 0;
		b = 1;
		c = 1;
	}
	else {
		a = 1;
		d = -1; 
		b = c = 0;
	}
	y = (a*cos(m*phi) + b*sin(m*phi))*dZDmn(m,n,R,Z) + (c*cos(m*phi) + d*sin(m*phi))*dZNmn(m,n-1,R,Z);
	//printf("BZ(%d, %d, %f, %f, %f) = %.10f\n", m, n, R, Z, phi, y);
	return y;
}

double Bphi(int m, int n, double R, double Z, double phi) {
	double y;
	int a, b, c, d;
	if (n%2 == 0) {
		a = d = 0;
		b = 1;
		c = 1;
	}
	else {
		a = 1;
		d = -1; 
		b = c = 0;
	}
	y = m*(-a*sin(m*phi) + b*cos(m*phi))*Dmn(m,n,R,Z)/R + m*(-c*sin(m*phi) + d*cos(m*phi))*Nmn(m,n-1,R,Z)/R;
	//printf("Bphi(%d, %d, %f, %f, %f) = %.10f\n", m, n, R, Z, phi, y);
	return y;
}

double dRBR(int m, int n, double R, double Z, double phi) {
	double y;
	int a, b, c, d;
	if (n%2 == 0) {
		a = d = 0;
		b = 1;
		c = 1;
	}
	else {
		a = 1;
		d = -1; 
		b = c = 0;
	}
	y = (a*cos(m*phi) + b*sin(m*phi))*dRRDmn(m,n,R,Z) + (c*cos(m*phi) + d*sin(m*phi))*dRRNmn(m,n-1,R,Z);
	return y;
}

double dZBZ(int m, int n, double R, double Z, double phi) {
	double y;
	int a, b, c, d;
	if (n%2 == 0) {
		a = d = 0;
		b = 1;
		c = 1;
	}
	else {
		a = 1 ;
		d = -1; 
		b = c = 0;
	}
	y = (a*cos(m*phi) + b*sin(m*phi))*dZZDmn(m,n,R,Z) + (c*cos(m*phi) + d*sin(m*phi))*dZZNmn(m,n-1,R,Z);
	return y;
}

double dRBZ(int m, int n, double R, double Z, double phi) {
	double y;
	int a, b, c, d;
	if (n%2 == 0) {
		a = d = 0;
		b = 1;
		c = 1;
	}
	else {
		a = 1; 
		d = -1; 
		b = c = 0;
	}
	y = (a*cos(m*phi) + b*sin(m*phi))*dRZDmn(m,n,R,Z) + (c*cos(m*phi) + d*sin(m*phi))*dRZNmn(m,n-1,R,Z);
	return y;
}

double dZBR(int m, int n, double R, double Z, double phi) {
	double y;
	int a, b, c, d;
	if (n%2 == 0) {
		a = d = 0; 
		b = 1;
		c = 1;
	}
	else {
		a = 1; 
		d = -1; 
		b = c = 0;
	}
	y = (a*cos(m*phi) + b*sin(m*phi))*dRZDmn(m,n,R,Z) + (c*cos(m*phi) + d*sin(m*phi))*dRZNmn(m,n-1,R,Z);
	return y;
}

double dRBphi(int m, int n, double R, double Z, double phi) {
	double y;
	int a, b, c, d;
	if (n%2 == 0) {
		a = d = 0;
		b = 1;
		c = 1;
	}
	else {
		a = 1 ; 
		d = -1; 
		b = c = 0;
	}
	y =  m*(-a*sin(m*phi) + b*cos(m*phi))*dRDmn(m,n,R,Z)/R + m*(-c*sin(m*phi) + d*cos(m*phi))*dRNmn(m,n-1,R,Z)/R - m*(-a*sin(m*phi) + b*cos(m*phi))*Dmn(m,n,R,Z)/pow(R,2) - m*(-c*sin(m*phi) + d*cos(m*phi))*Nmn(m,n-1,R,Z)/pow(R,2);
	return y;
}

double dZBphi(int m, int n, double R, double Z, double phi) {
	double y;
	int a, b, c, d;
	if (n%2 == 0){
		a = d = 0;
		b = 1;
		c = 1;
	}
	else {
		a = 1 ;
		d = -1;
		b = c = 0;
	}
	y = m*(-a*sin(m*phi) + b*cos(m*phi))*dZDmn(m,n,R,Z)/R + m*(-c*sin(m*phi) + d*cos(m*phi))*dZNmn(m,n-1,R,Z)/R;
	return y;
}

struct field *DommBfield(int N_modes, double *amp, int *tor_mode, int *pol_mode, double R, double Z, double phi) {
	struct field *magfield=malloc(sizeof(struct field));
	int ind;
	//printf("YOLOO\n");
	magfield->value[0] = BRloop(R, Z);
	magfield->value[1] = BZloop(R, Z);
	magfield->value[2] = mu0*IZ/(2.0*M_PI*R);
	magfield->derivative[0][0] = dRBRloop(R, Z);
	magfield->derivative[1][0] = dRBZloop(R, Z);
	magfield->derivative[2][0] = -mu0*IZ/(2.0*M_PI*R*R);
	magfield->derivative[0][1] = dZBRloop(R, Z);
	magfield->derivative[1][1] = dZBZloop(R, Z);
	magfield->derivative[2][1] = 0.0;
	for (ind=0; ind<N_modes; ind++) {
		//printf("ind=%d\n", ind);
		magfield->value[0] += amp[ind]*BR(tor_mode[ind], pol_mode[ind], R, Z, phi);
		magfield->value[1] += amp[ind]*BZ(tor_mode[ind], pol_mode[ind], R, Z, phi);
		magfield->value[2] += amp[ind]*Bphi(tor_mode[ind], pol_mode[ind], R, Z, phi);
		magfield->derivative[0][0] += amp[ind]*dRBR(tor_mode[ind], pol_mode[ind], R, Z, phi);
		magfield->derivative[1][0] += amp[ind]*dRBZ(tor_mode[ind], pol_mode[ind], R, Z, phi);
		magfield->derivative[2][0] += amp[ind]*dRBphi(tor_mode[ind], pol_mode[ind], R, Z, phi);
		magfield->derivative[0][1] += amp[ind]*dZBR(tor_mode[ind], pol_mode[ind], R, Z, phi);
		magfield->derivative[1][1] += amp[ind]*dZBZ(tor_mode[ind], pol_mode[ind], R, Z, phi);
		magfield->derivative[2][1] += amp[ind]*dZBphi(tor_mode[ind], pol_mode[ind], R, Z, phi);
		//printf("amp[ind]=%f\n", amp[ind]);
	}
	//printf("R=%f\n", R);
	//printf("Z=%f\n", Z);
	//printf("phi=%f\n", phi);
	//printf("BR=%f\n", magfield->value[0]);
	//printf("BZ=%f\n", magfield->value[1]);
	//printf("Bphi%f\n", magfield->value[2]);
	return magfield;
}

