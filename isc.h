#define small 0.0000000001

/* defines a function that extracts data from lines of a file */
double *linetodata(char *line, int *size);

/* defines a structure containing the magnetic field value and its derivatives (assigned by Bfield function below) */
struct field {
	double value[3];
	double derivative[3][2];
};

/* defines a structure containing the position and tangent maps on a point along a magnetic field line */
struct position {
	double loc[2];
	double **tangent;
};

/* defines a structure containing position and several different tangent maps of a point along a magnetic field line 
   contains the eigenvector of the symmetrized full orbit tangent map (which are the island axes) 
   contains the angle neighbouring points rotate around the island centre after a full orbit 
   contans q0 index = number of field periods necessary for the rotation of neighbouring points to be pi/2 */
struct ext_position {
	double loc[2];
	double circumference;
	double **part_tangent;
	double **full_tangent;
	double angle;
	double *epar;
	double *eperp;
	int q0_index;
	double ***long_tangent;
	double *eparkdotlong_tangentdoteperp0;
};

/* prints the nrows x ncols matrix input (and the name of the matrix) */
void printmat(char *name, double **input, int nrows, int ncols);

/* returns the inverse of the 2x2 matrix input2x2, and assigns the determinant to det_tangent */
double **invert2x2(double **input2x2, double* det_tangent); 

/* returns the 2xdims matrix product of 2x2 matrix input1 and 2xdims matrix input2 */
double **multiply2x2(double **input2x2, double **input2xn, int dims);

/* mutliply 2x2 matrix input1 by 2x2 matrix input2 and assign the result to inputargument_number */
void multiply2x2reassign(double **input1, double **input2, int argument_number);

void linalg2x2(double **input, double **evecs, double *evals, double *determinant, double *trace);

void **symmeigs(double **input, double *evec_largeeval, double *evec_smalleval);

double inner(double *left_vect, double **matrix, double *right_vect);

double **set_identity();

struct position addstructs(double num1, struct position *struct1, double num2, struct position *struct2);

double ***coil_grid();

struct field *Bfield(double *Xp, double varphi, double ***coils, int *num_coils, int **num_segs);
 
void RK4(struct position *Xp, double varphi, double dvarphi, double ***coils, int *num_coils, int **num_segs);

struct position *findcentre(double ***coils, int *n_coils, int **n_segs, struct position *fieldline);

void iotaprofile(double *rmin, double *iota, double ***coils, int *n_coils, int **n_segs);

struct position *findisland(double ***coils, int *n_coils, int **n_segs, struct position *fieldline, int tor_mode, int pol_mode);

struct ext_position *alongcentre(double RR, double ZZ, int tor_mode, int pol_mode, double ***coils, int *n_coils, int **n_segs);

double *islandwidth(struct ext_position *ext_fieldline, int tor_mode, int pol_mode);
