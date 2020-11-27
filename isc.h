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
// this is evaluated as a pointer for every periodic field line position in the set of symmetry phi-planes of choice
	double loc[2]; // (R,Z) coordinates of the periodic field line at the symmetry phi-planes
	double **part_tangent; // partial tangent map from current to next symmetry plane
	double **full_tangent; // full orbit tangent map from current symmetry plane to the Lth one after
	double **adj_part_tangent; // adjoint of above
	double **adj_full_tangent; // adjoint of above 
	// the above two are actually not necessary: the adjoint tangent map is the transpose of the inverse tangent map
	// however, they are kept and evaluated separately (mostly because I initially did not realise this)
	double *epar; // unit eigenvector parallel to flux surface 
	double *eperp; // unit eigenvector perpendicular to flux surface
	// the above are obtained from the W matrix which is obtained from the full orbit tangent map
	double ***long_tangent; // the tangent map from the current symmetry plane to the (q0+q)th one after, for 0<=q<L
	double **sperp; // probably unnecessary but for now I need to keep it
	double sperp_final[2]; // ditto
	double chord[2]; // chord joining fixed point position at previous poloidally reordered position with current one
	double chordplus[2]; // chord joining fixed point position at current location with next poloidally reordered one
	double width; // island width
	double Sigma; // Sigma quantity entering in island width Sigma = epar_k+q . S_k^q . eperp_k
	int *sign; // this is probably not necessary, I will get rid of it (DELETE unless good reason found)
	// these final numbers are equivalent for all fixed points in the island chain
	// this introduces a redundancy in the pointer to values of ext_position, as all values below repeat
	//int q0_index; // number of symmetry planes one needs to follow the linearized field line for it to rotate ~=90deg
	double circ; // circumference
	double Res; // Residue
};


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

void transpose2x2reassign(double **input2x2);

double **transpose2x2(double **input2x2);

double **adj2x2(double **input2x2, double* det_tangent); 
/* returns the transpose of the inverse of the 2x2 matrix input2x2, and assigns the determinant to det_tangent */

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

//struct position *solve_magneticaxis_save(struct fieldparams allparams, struct field **Bfield_saved, struct position *fieldline, int N_gridphi_fieldperiod);
void solve_magneticaxis_save(double *axisguess, struct position **centre, struct fieldparams allparams, struct field **Bfieldsaved, int N_gridphi_fieldperiod);
/* finds the magnetic axis using a Newton iteration */

void iotaprofile(struct position axis, int m0_symmetry, int N_gridphi_per_field_period, double *rmin, double *zmin, double *iota, struct fieldparams allparams) ;
//void iotaprofile(double RRin, double r_interval, int n_points, int m0_symmetry, int N_gridphi_per_field_period, double *minor_radius, double *iota, double ***coils, int n_coils, int *n_segs);
/* computes iota as a function of distance from axis in a specified direction */

struct position *solve_islandcenter(struct fieldparams allparams, struct field **Bfield_saved, struct position *fieldline, int L_fixedpoints, int N_gridphi_fieldperiod);
/* finds the centre of an island from a nearby guess using a Newton iteration */

void solve_islandcenter_save(struct position **centre, double *Res, struct fieldparams allparams, struct field **Bfieldsaved, int L_fixedpoints, int N_gridphi_fieldperiod) ;

void extsolve_periodicfieldline(struct ext_position *ext_centre, struct position *fieldline, struct position *axis, int m0_symmetry, int L_fixedpoints, int pol_mode, int *q0_index, int N_gridphi_fieldperiod);
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
void solve_lambdaRes(struct position **lambdaRes, struct position **Xp_saved, int m0_symmetry, int L_fixedpoints, int N_gridphi_fieldperiod, struct field **Bfield_saved) ;

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

