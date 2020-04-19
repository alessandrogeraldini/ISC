#define small 0.0000000001

//STRUCTURES

struct field {
	double value[3];
	double derivative[3][2];
	double twoderivative[3][2][2];
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
	double **adj_part_tangent;
	double **adj_full_tangent;
	double angle;
	double *epar;
	double *eperp;
	int q0_index;
	double ***long_tangent;
	double **sperp;
	double sperp_final[2];
	double chord[2];
	double chordplus[2];
/* contains position and several different tangent maps of a point along a magnetic field line 
   contains the eigenvector of the symmetrized full orbit tangent map (which are the island axes) 
   contains the angle neighbouring points rotate around the island centre after a full orbit 
   contains q0 index = number of field periods necessary for the rotation of neighbouring points to be pi/2 */
};

//struct field_grad {
//	double value[3];
//	double derivative[3][2];
//};
///* contains the shape gradients of the magnetic field and its derivative
//   for coils, this structure will be an array of the same configuration as the coils */
//   there are three gradients per parameter to take into account vector parameters (coil positions)
//   for scalar parameters, only the first component is non-zero

/* contains the gradients of several quantities needed to evaluate the island width
   also contains the gradient of the island width
   to be used with analytical configurations where the island can be tuned by varying a single (few) parameter(s) */ 

/******         *******/

// FUNCTIONS

int fac(int number);
/* factorial function */

double *linetodata(char *line, int *size);
/* extracts data from lines of a file */

void printmat(char *name, double **input, int nrows, int ncols);
/* prints the nrows x ncols matrix input (and the name of the matrix) */

void printstructposition(char *name, struct position *input);

double **invert2x2(double **input2x2, double* det_tangent); 
/* returns the inverse of the 2x2 matrix input2x2, and assigns the determinant to det_tangent */

double **multiply2x2(double **input2x2, double **input2xn, int dims);
/* returns the 2xdims matrix product of 2x2 matrix input1 and 2xdims matrix input2 */

double *multiply(double **input2x2, double *inputvec);
/* returns the product of 2x2 matrix input1 and 2-vector input2 */

void multiply2x2reassign(double **input1, double **input2, int argument_number);
/* mutliply 2x2 matrix input1 by 2x2 matrix input2 and assign the result to inputargument_number */

double **add2x2(double c1, double **input2x2, double c2, double **input2xn, int dims);
/* returns the 2xdims matrix product of 2x2 matrix input1 and 2xdims matrix input2 */

void add2x2reassign(double c1, double **input1, double c2, double **input2, int argument_number);
/* mutliply 2x2 matrix input1 by 2x2 matrix input2 and assign the result to inputargument_number */

void linalg2x2(double **input, double **evecs, double *evals, double *determinant, double *trace);

void **symmeigs(double **input, double *evec_largeeval, double *evec_smalleval, double *evals);

double inner(double *left_vect, double **matrix, double *right_vect);

double **set_identity();

double **set_zeros();

struct position addstructs(double num1, struct position *struct1, double num2, struct position *struct2);

struct field addstructsfield(double num1, struct field *struct1, double num2, struct field *struct2);

double ***coil_grid();
/* forms an array with the coil information
   the array can be fed into the function Bfield
   to evaluate the magnetic field structure at each point */

struct field *Bfield(double *Xp, double varphi, double ***coils, int num_coils, int *num_segs);
/* magnetic field evaluation at each point Xp->loc */

struct field *gradBfield(double *Xp, double varphi, double *param);
/* magnetic field (eventually shape) gradient evaluation at each point Xp->loc */

struct field ***shapeBfield(double *Xp, double varphi, double ***coils, int num_coils, int *num_segs);

struct field *DommBfield(int N_modes, double *amp, int *tor_mode, int *pol_mode, double R, double Z, double phi); 
/* magnetic field evaluation using Dommaschk potentials */
 
void RK4step(struct position *Xp, double varphi, double dvarphi, double ***coils, int num_coils, int *num_segs, struct field *Bfieldsaved);
/* takes a Runge Kutta 4th order step (in toroidal angle) 
   follows field line position and tangent map */

void RK4step_lambdacirc_mutangent(struct position *Xp, struct position *adjvariable, double varphi, double dvarphi, struct field *Bfield_saved);
/* takes a Runge Kutta 4th order step (in toroidal angle) 
   follows adjoint variables of field line position and tangent map 
   includes only homogeneous terms (linear in adjvariable), jumps from source terms taken care of separately */

void RK4step_lambdatangent(struct position *Xp, struct position *lambda, struct position *sperp, struct position *mu, double varphi, double dvarphi, struct field *Bfieldsaved);
/* takes a Runge Kutta 4th order step (in toroidal angle) 
   follows adjoint variable of field line position with additional constraint from tangent map 
   includes only homogeneous terms (linear in mu or lambda), jumps from source terms taken care of separately */

void RK4step_gradcirc(double *number, struct position *centre, struct position *lambda_circ, double varphi, double dvarphi, struct field *Bfieldsaved, double *param);
/* takes a Runge Kutta 4th order step (in toroidal angle) 
   follows shape gradients of field line circumference, tangent map (eventually), and island width (eventually) */

void RK4step_gradtangent(double *number, struct position *Xp, struct position *lambda, struct position *sperp, struct position *mu, double varphi, double dvarphi, struct field *Bfield_saved) ;

void RK4_adjshapecirc(double ***shapecirc, struct position *Xp, struct position *lambda, double varphi, double dvarphi, double ***coils, int num_coils, int *num_segs);
/* takes a Runge Kutta 4th order step (in toroidal angle) 
   follows shape gradients of field line circumference */

struct position *solve_magneticaxis(double ***coils, int n_coils, int *n_segs, struct field **Bfield_saved, struct position *fieldline, int N_gridphi_toroidal);
/* finds the magnetic axis using a Newton iteration */

void iotaprofile(double RRin, double r_interval, int n_points, int m0_symmetry, int N_gridphi_per_field_period, double *minor_radius, double *iota, double ***coils, int n_coils, int *n_segs);
/* computes iota as a function of distance from axis in a specified direction */

struct position *solve_islandcenter(double ***coils, int n_coils, int *n_segs, struct field **Bfield_saved, struct position *fieldline, int m0_fieldperiods, int L_fixedpoints, int N_gridphi_fieldperiod);
/* finds the centre of an island from a nearby guess using a Newton iteration */

struct ext_position *solve_islandcenter_full(double *islandcenter, double *axis, int *n_turns, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod, struct field **Bfield_saved);
/* evaluates all quantities necessary to calculate the island width */

double *calc_islandwidth(struct ext_position *ext_fieldline, int m0_symmetry, int L_fixedpoints, int pol_mode, int N_gridphi_fieldperiod);
/* evaluates island width */

double *islandwidthnew(struct field **Bfield_saved, struct ext_position *ext_centre, struct position *lambda_circ, struct position **lambda_tangent, struct position **mu_tangent, int m0_fieldperiods, int N_polmode, int L_fixedpoints, int N_gridphi_fieldperiod, double *param) ;

struct position solve_lambda_circ(struct ext_position *ext_centre, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod);
/* finds the periodic adjoint solution lambda (its initial value) that constrains the field line to be periodic */

struct position **solve_mu_tangent(struct ext_position *ext_centre, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod);

struct position **solve_lambda_tangent(struct field **Bfield_saved, struct ext_position *ext_centre, struct position **mu, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod);
/* finds the adjoint solutions that constrain the field line to be periodic (lambdaQ) 
   and constrain the initial tangent map variation to be zero (mu) */

double *solve_gradcirc(struct field **Bfield_island, struct ext_position *ext_centre, struct position *lambda_circ, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod, double *param) ;
/* evaluates gradient of circumference */

double *solve_gradtangent(struct field **Bfield_island, struct ext_position *ext_centre, struct position **lambdaQ, struct position **muQ, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod, double *param);
/* evaluates gradient of Sigma (sum of tangent map matrix elements) */


//double ***shapecircumference(double ***coils, int n_coils, int *n_segs, struct ext_position *ext_centre, struct position *lambda_circ, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode) ;
/* evaluates shape gradient of circumference */

//can delete this eventually
//struct ext_position *gradalongcentrealt(double RR, double ZZ, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode, double ***coils, int n_coils, int *n_segs);
/* evaluates all quantities necessary to calculate the island width gradient */

//can delete this eventually
//double *gradislandwidth(struct ext_position *ext_fieldline, struct ext_position *grad_ext_fieldline, int field_periods, int N_gridphi_per_field_period, int tor_mode, int pol_mode);
/* evaluates island width gradient */

//can delete this eventually
//struct position *gradcentre(double RR, double ZZ, int m0_symmetry, int N_gridphi_per_field_period, int tor_mode, int pol_mode, double ***coils, int n_coils, int *n_segs);
/* evaluates gradient of O point position and tangent map */

