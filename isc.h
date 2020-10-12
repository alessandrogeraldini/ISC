#define small 0.0000000001

//STRUCTURES

struct fieldparams {
	char type[5];
	int m0_fieldperiods;
	int n_coils;
	int *intparams;
	double *constparams;
	double ***diffparams;
	int n_diff;
};
/* contains the magnetic field parameters */

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
	double **part_tangent;
	double **full_tangent;
	double **adj_part_tangent;
	double **adj_full_tangent;
	double *epar;
	double *eperp;
	int *sign;
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

struct islandchain {
   int L;
   double circumference;
   double Res;
   double *Sigma;
   double *width;
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
 
struct fieldparams fetchparams();

void printmat(char *name, double **input, int nrows, int ncols);
/* prints the nrows x ncols matrix input (and the name of the matrix) */

void printstructposition(char *name, struct position *input);

void printstructfield(char *name, struct field *input);

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

double ***fieldparams(char *type, int *m0_fieldperiods, int *n_coils, int **n_consts, double **n_params);
/* forms an array with the coil information
   the array can be fed into the function Bfield
   to evaluate the magnetic field structure at each point */

//struct field Bfield(double *Xp, double varphi, char *type, int m0_fieldperiods, double ***diffparams, int n_coils, int *intparams, double *constparams);
struct field Bfield(double *Xp, double varphi, struct fieldparams allparams) ;
/* magnetic field evaluation at each point Xp->loc */

struct field *gradBfield(double *Xp, double varphi, struct fieldparams allparams, int diffparam_ind1, int diffparam_ind2);
/* magnetic field (eventually shape) gradient evaluation at each point Xp->loc */

//struct field ***shapeBfield(double *Xp, double varphi, double ***coils, int num_coils, int *num_segs);

struct field *DommBfield(int N_modes, double *amp, int *tor_mode, int *pol_mode, double R, double Z, double phi); 
/* magnetic field evaluation using Dommaschk potentials */
 
void RK4step(struct position *Xp, double varphi, double dvarphi, struct fieldparams allparams, struct field *Bfieldsaved);
/* takes a Runge Kutta 4th order step (in toroidal angle) 
   follows field line position and tangent map */

//struct position RK4stepsave(struct position *Xp, double varphi, double dvarphi, struct fieldparams allparams, struct field *Bfield_saved);
void RK4stepsave(struct position *Xp_adv, struct position *Xp, double varphi, double dvarphi, struct fieldparams allparams, struct field *Bfield_saved) ;

void RK4step_lambdacirc_mutangent(struct position *Xp, struct position *adjvariable, double varphi, double dvarphi, struct field *Bfield_saved);
/* takes a Runge Kutta 4th order step (in toroidal angle) 
   follows adjoint variables of field line position and tangent map 
   includes only homogeneous terms (linear in adjvariable), jumps from source terms taken care of separately */

void RK4step_lambdatangent(struct position *Xp, struct position *lambda, struct position *sperp, struct position *mu, double varphi, double dvarphi, struct field *Bfieldsaved);
/* takes a Runge Kutta 4th order step (in toroidal angle) 
   follows adjoint variable of field line position with additional constraint from tangent map 
   includes only homogeneous terms (linear in mu or lambda), jumps from source terms taken care of separately */

void RK4step_lambdaRes(struct position *lambdamu_out, struct position *lambdamu, struct position *Xp, double varphi, double dvarphi, struct field *Bfield_saved);

void RK4step_gradcirc(double *number, struct position *Xp, struct position *lambda, double varphi, double dvarphi, struct field *Bfield_saved, struct fieldparams allparams, int diffparams_ind1, int diffparams_ind2) ;
/* takes a Runge Kutta 4th order step (in toroidal angle) 
   follows shape gradients of field line circumference, tangent map (eventually), and island width (eventually) */

void RK4step_gradtangent(double *number, struct position *Xp, struct position *lambda, struct position *sperp, struct position *mu, double varphi, double dvarphi, struct field *Bfield_saved, struct fieldparams allparams, int diffparams_ind1, int diffparams_ind2) ;

//void RK4step_gradRes(double number, &fieldline, &lambdamu, varphi, -dvarphi, Bfield_island[i % N_line], allparams, diffparams_ind1, diffparams_ind2);
void RK4step_gradRes(double *number, struct position *Xp, struct position *lambdamu, double varphi, double dvarphi, struct field *Bfield_saved, struct fieldparams allparams, int diffparam_ind1, int diffparam_ind2);

void RK4_adjshapecirc(double ***shapecirc, struct position *Xp, struct position *lambda, double varphi, double dvarphi, double ***coils, int num_coils, int *num_segs);
/* takes a Runge Kutta 4th order step (in toroidal angle) 
   follows shape gradients of field line circumference */

struct position *solve_magneticaxis(struct fieldparams allparams, struct field **Bfield_saved, struct position *fieldline, int N_gridphi_fieldperiod);
/* finds the magnetic axis using a Newton iteration */

struct position *solve_magneticaxis_save(struct fieldparams allparams, struct field **Bfield_saved, struct position *fieldline, int N_gridphi_fieldperiod);
/* finds the magnetic axis using a Newton iteration */

void iotaprofile(char *type, struct position axis, int m0_symmetry, int N_gridphi_per_field_period, double *rmin, double *zmin, double *iota, double ***coils, int n_coils, int *n_segs) ;
//void iotaprofile(double RRin, double r_interval, int n_points, int m0_symmetry, int N_gridphi_per_field_period, double *minor_radius, double *iota, double ***coils, int n_coils, int *n_segs);
/* computes iota as a function of distance from axis in a specified direction */

struct position *solve_islandcenter(struct fieldparams allparams, struct field **Bfield_saved, struct position *fieldline, int L_fixedpoints, int N_gridphi_fieldperiod);
/* finds the centre of an island from a nearby guess using a Newton iteration */

struct position **solve_islandcenter_save(double Rguess, double Zguess, struct fieldparams allparams, struct field **Bfieldsaved, int L_fixedpoints, int N_gridphi_fieldperiod) ;

struct ext_position *solve_islandcenter_full(double *islandcenter, struct position *axis, int *n_turns, int m0_symmetry, int L_fixedpoints, double *Res, int *q0_index, int N_gridphi_fieldperiod, struct field **Bfield_saved);
/* evaluates all quantities necessary to calculate the island width */

double *calc_islandwidth(double *circ, double *Sigma, struct ext_position *ext_fieldline, int m0_symmetry, int L_fixedpoints, int pol_mode);
/* evaluates island width */

double *islandwidthnew(struct field **Bfield_saved, struct ext_position *ext_centre, struct position *lambda_circ, struct position **lambda_tangent, struct position **mu_tangent, int m0_fieldperiods, int N_polmode, int L_fixedpoints, int N_gridphi_fieldperiod, char *type, double *param, int num_params) ;

struct position solve_lambda_circ(struct ext_position *ext_centre, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod);
/* finds the periodic adjoint solution lambda (its initial value) that constrains the field line to be periodic */

struct position **solve_mu_tangent(struct ext_position *ext_centre, int m0_symmetry, int L_fixedpoints, int q0_index, int N_gridphi_fieldperiod);

struct position **solve_lambda_tangent(struct field **Bfield_saved, struct ext_position *ext_centre, struct position **mu, int m0_symmetry, int L_fixedpoints, int q0_index, int N_gridphi_fieldperiod);
/* finds the adjoint solutions that constrain the field line to be periodic (lambdaQ) 
   and constrain the initial tangent map variation to be zero (mu) */
struct position **solve_lambdaRes(struct field **Bfield_saved, struct position **Xp_saved, double **adjfulltangent, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod) ;

//void solve_gradcirc(double *gradcircumference, struct field **Bfield_island, struct ext_position *ext_centre, struct position *lambda_circ, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod, char *type, double *param, int num_params, int *params2) ;
void solve_gradcirc(double *gradcircumference, struct field **Bfield_island, struct ext_position *ext_centre, struct position *lambda_circ, int L_fixedpoints, int N_gridphi_fieldperiod, struct fieldparams allparams, int diffparams_ind1, int diffparams_ind2) ;
/* evaluates gradient of circumference */

//void solve_gradtangent(double **gradtangent, struct field **Bfield_island, struct ext_position *ext_centre, struct position **lambdaQ, struct position **muQ, int m0_symmetry, int L_fixedpoints, int q0_index, int N_gridphi_fieldperiod, char *type, double *param, int num_params, int *params2);
void solve_gradtangent(double **number, struct field **Bfield_island, struct ext_position *ext_centre, struct position **lambdaQ, struct position **muQ, int L_fixedpoints, int q0_ind, int N_gridphi_fieldperiod, struct fieldparams allparams, int diffparams_ind1, int diffparams_ind2) ;
/* evaluates gradient of Sigma (sum of tangent map matrix elements) */

//void solve_gradRes(double *number, struct field **Bfield_island, struct position fieldline, struct position lambdamu, int L_fixedpoints, int N_gridphi_fieldperiod, struct fieldparams allparams, int diffparams_ind1, int diffparams_ind2) ;
void solve_gradRes(double *number, struct field **Bfield_island, struct position **fieldline, struct position **lambdamu, int L_fixedpoints, int N_gridphi_fieldperiod, struct fieldparams allparams, int diffparams_ind1, int diffparams_ind2) ;

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

