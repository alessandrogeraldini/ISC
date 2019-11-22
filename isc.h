#define small 0.0000000001

double *linetodata(char *line, int *size);

struct field {
	double value[3];
	double derivative[3][2];
};

struct position {
	double loc[2];
	double tangent[2][2];
};

struct ext_position {
	double loc[2];
	double **part_tangent;
	double **full_tangent;
	double **symmfull_eig;
	double **long_tangent;
	double angle;
	int q0_index;
};

struct field *Bfield(double *Xp, double varphi, double ***coils, int *num_coils, int **num_segs);

double ***coil_grid();
 
void RK4(struct position *Xp, double varphi, double dvarphi, double ***coils, int *num_coils, int **num_segs);

double **invert2x2(double input2x2[2][2], double* det_tangent); 

double **multiply2x2(double **input2x2, double **input2xn, int dims);

struct position *findcentre(double ***coils, int *n_coils, int **n_segs, struct position *fieldline);

struct position *findisland(double ***coils, int *n_coils, int **n_segs, struct position *fieldline, int tor_mode, int pol_mode);

struct ext_position *alongcentre(double RR, double ZZ, int tor_mode, int pol_mode, double ***coils, int *n_coils, int **n_segs);

struct position addstructs(double num1, struct position *struct1, double num2, struct position *struct2);

void iotaprofile(double *rmin, double *iota, double ***coils, int *n_coils, int **n_segs);

void linalg2x2(double input[2][2], double evecs[2][2], double evals[2], double *determinant, double *trace);
