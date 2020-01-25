#define small 0.0000000001

struct field {
	double value[3];
	double derivative[3][2];
};
/* contains the magnetic field value and its derivatives (assigned by Bfield function below) */

struct position {
	double loc[2];
	double **tangent;
};
/* contains the position and tangent maps on a point along a magnetic field line */

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
/* contains position and several different tangent maps of a point along a magnetic field line 
   contains the eigenvector of the symmetrized full orbit tangent map (which are the island axes) 
   contains the angle neighbouring points rotate around the island centre after a full orbit 
   contains q0 index = number of field periods necessary for the rotation of neighbouring points to be pi/2 */

//struct field_grad {
//	double value[3];
//	double derivative[3][2];
//};
///* contains the shape gradients of the magnetic field and its derivative
//   for coils, this structure will be an array of the same configuration as the coils */
//
//struct position_grad {
//	double loc_grad[2];
//	double **tangent_grad;
//};
///* contains the shape gradients of the position and tangent map on a point along a magnetic field line 
//   for coils, this structure will be an array of the same configuration as the coils */
//
//   //there are three gradients per parameter to take into account vector parameters (coil positions)
//   //for scalar parameters, only the first component is non-zero

/* contains the gradients of several quantities needed to evaluate the island width
   also contains the gradient of the island width
   to be used with analytical configurations where the island can be tuned by varying a single (few) parameter(s) */ 

// FUNCTIONS

double *linetodata(char *line, int *size);
/* extracts data from lines of a file */

void printmat(char *name, double **input, int nrows, int ncols);
/* prints the nrows x ncols matrix input (and the name of the matrix) */

double **invert2x2(double **input2x2, double* det_tangent); 
/* returns the inverse of the 2x2 matrix input2x2, and assigns the determinant to det_tangent */

double **multiply2x2(double **input2x2, double **input2xn, int dims);
/* returns the 2xdims matrix product of 2x2 matrix input1 and 2xdims matrix input2 */

void multiply2x2reassign(double **input1, double **input2, int argument_number);
/* mutliply 2x2 matrix input1 by 2x2 matrix input2 and assign the result to inputargument_number */

double **add2x2(double c1, double **input2x2, double c2, double **input2xn, int dims);
/* returns the 2xdims matrix product of 2x2 matrix input1 and 2xdims matrix input2 */

void add2x2reassign(double c1, double **input1, double c2, double **input2, int argument_number);
/* mutliply 2x2 matrix input1 by 2x2 matrix input2 and assign the result to inputargument_number */

void linalg2x2(double **input, double **evecs, double *evals, double *determinant, double *trace);

void **symmeigs(double **input, double *evec_largeeval, double *evec_smalleval);

double inner(double *left_vect, double **matrix, double *right_vect);

double **set_identity();

double **set_zeros();

struct position addstructs(double num1, struct position *struct1, double num2, struct position *struct2);

struct field addstructsfield(double num1, struct field *struct1, double num2, struct field *struct2);

double ***coil_grid();

struct field *Bfield(double *Xp, double varphi, double ***coils, int *num_coils, int **num_segs);

struct field *gradBfield(double *Xp, double varphi, double ***coils, int *num_coils, int **num_segs);

struct field *DommBfield(int N_modes, double *amp, int *tor_mode, int *pol_mode, double R, double Z, double phi); 
 
void RK4(struct position *Xp, double varphi, double dvarphi, double ***coils, int *num_coils, int **num_segs);
/* takes a Runge Kutta 4th order step (in toroidal angle) 
   follows field line position and tangent map */

void RK4_wgrad(struct position *Xp, struct position *Xpgrad, double varphi, double dvarphi, double ***coils, int *num_coils, int **num_segs);
/* takes a Runge Kutta 4th order step (in toroidal angle) 
   follows shape gradients of field line position and tangent map */

struct position *findcentre(double ***coils, int *n_coils, int **n_segs, struct position *fieldline, int N_gridphi_toroidal);

void iotaprofile(double RRin, double r_interval, int n_points, int m0_symmetry, int N_gridphi_per_field_period, double *minor_radius, double *iota, double ***coils, int *n_coils, int **n_segs);

struct position *findisland(double ***coils, int *n_coils, int **n_segs, struct position *fieldline, int field_periods, int N_gridphi_per_field_period, int tor_mode, int pol_mode);
/* finds the centre of an island from a nearby guess */

struct ext_position *alongcentre(double RR, double ZZ, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode, double ***coils, int *n_coils, int **n_segs);
/* evaluates all quantities necessary to calculate the island width */

struct ext_position *gradalongcentre(double RR, double ZZ, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode, double ***coils, int *n_coils, int **n_segs);
/* evaluates all quantities necessary to calculate the island width gradient */

double *islandwidth(struct ext_position *ext_fieldline, int field_periods, int N_gridphi_per_field_period, int tor_mode, int pol_mode);
/* evaluates island width */

double *gradislandwidth(struct ext_position *ext_fieldline, struct ext_position *grad_ext_fieldline, int field_periods, int N_gridphi_per_field_period, int tor_mode, int pol_mode);
/* evaluates island width gradient */

struct position *gradcentre(double RR, double ZZ, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode, double ***coils, int *n_coils, int **n_segs);
/* evaluates gradient of O point position and tangent map */

int fac(int number);
/* factorial function */
