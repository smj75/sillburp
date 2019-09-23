
#ifndef _SILLBURP_H
#define _SILLBURP_H

#define ABSOLUTE_ZERO -273.15

/* From Maclennan & Lovell, Geology, 2002 */

#define LATENT_HEAT 700e3
#define SPECIFIC_HEAT 1.7e3

#define LABILE 1
#define REFRACTORY 2
#define VITRINITE 3
#define OIL 4

/*
	Parameters to control the 4 organic reactions.
  reaction 1: CH gas from labile kerogen
  reaction 2: CH gas direct from refractory kerogen
  reaction 3: Vitrinite reflectance
  reaction 4: CH gas cracked from oil from labile kerogen
  1, 2 & 4 PARAMETERS FROM QUIGLEY & MACKENZIE, NATURE, 1988
  3 PARAMETERS FROM ROWLEY & WHITE, EPSL, 1997
 */

const int n_reactions = 4;
double mean_pre_exp_const[n_reactions+1] = {0.0, 
	1.58e13, 
	1.83e18, 
	4.00e10, 
	1.00e13	}; // Unit: /s
double mean_activation_energy[n_reactions+1] = {0.0, 
	208.0e3, 
	279.0e3, 
	242.0e3, 
	230.0e3	}; // Unit: J
double standard_deviation[n_reactions+1] = {0.0, 
	5.0e3, 
	13.0e3, 
	41.0e3, 
	5.0e3	}; // Unit: J
const double mass_frac_labile_to_gas = 0.2;
const double gas_constant = 8.314510;
const double carbon_to_methane = 1.34;
const double def_edge_reaction_zone[n_reactions+1] = {0.0,
	0.01,
	0.01,
	0.01,
	0.01};

/* 
	Internal function definitions
 */

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
static double minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

/* FUNCTION DEFINITIONS */

double arrhenius_reaction_rate_coeff( 
	double energy, 
	double temp_c );

void assign_aureole_uncertainties (
	void );

void avevar(
	double data[], 
	int n, 
	double *ave, 
	double *var	);

double brent(
	double ax, 
	double bx, 
	double cx, 
	double (*f)(), 
	double tol, 
	double *xmin);

void calculate_and_print_profile_single_step (
	double end_time );

void calculate_reaction_profiles(
	double end_time );

void calculate_temperature_profile(
	double end_time );

double distance_integral ( 
	double dist   );
	
double distance_wrt_sheet_margin ( 
	double distance );

double energy_integral ( 
	double energy   );
	
double expint(
	int n, 
	double x);

double f1dim(
	double x	);

void fdjac( 
	int n, 
	double x[], 
	double fvec[], 
	double **df, 
	void (*vecfunc)()   );
	
void find_aureole_edges (  
	int i_grid_top, 
	int i_grid_base, 
	int i_integrate_down_to,
	double time_step  );

void find_aureole_misfit ( void );
	
double find_dless_solidification_boundary( 
	double dless_solidification_boundary   );	
double left_side;

void find_solidification_time ( void );

double find_time_step (
	double time	);

double fmin(
	double x[]  );

double gaussian_distribution( 
	double energy   );

double get_progress ( 
	int i_reaction,
	int i_grid  );

double get_rate ( 
	int i_reaction,
	int i_grid  );

double ierf(
	double xin  );

void integrate_aureole ( 
	int i_grid_top, 
	int i_grid_base, 
	int i_integrate_down_to,
	double time, 
	double time_step );

void integrate_dless_aureole ( 
	int i_grid_top, 
	int i_grid_base, 
	int i_integrate_down_to );

void linmin(
	double p[], 
	double xi[], 
	int n, 
	double *fret, 
	double (*func)(double [])	);

void load(
	double x1,
	double v[],
	double T[]  );

void lnsrch(
	int n, 
	double xold[], 
	double fold, 
	double g[], 
	double p[], 
	double x[], 
	double *f, 
	double stpmax, 
	int *check, 
	double (*func)()	);

void lubksb(
	double **a, 
	int n, 
	int *indx, 
	double b[]  );

void ludcmp(
	double **a, 
	int n, 
	int *indx, 
	double *d,
	int *check   );

void machar(
	int *ibeta,
	int *it,
	int *irnd,
	int *ngrd,
	int *machep,
	int *negep,
	int *iexp,
	int *minexp,
	int *maxexp,
	double *eps,
	double *epsneg,
	double *xmin,
	double *xmax	);

double midexp(
	double (*funk)(double),
	double aa,
	double bb,
	int n   );

double midpnt_energy(
	double (*func)(double),
	double a,
	double b,
	int n   );

double midpnt_time(
	double (*func)(double),
	double a,
	double b,
	int n   );

void mnbrak(
	double *ax, 
	double *bx, 
	double *cx, 
	double *fa, 
	double *fb, 
	double *fc, 
	double (*func)()	);

void newt(
	double x[], 
	int n, 
	int *check, 
	void (*vecfunc)(int, double*, double*)   );

void odeint(
	double ystart[], 
	int nvar, 
	double x1, 
	double x2, 
	double eps,
	double h1,
	double hmin,
	int *nok,
	int *nbad,
	void (*derivs)(double, double*, double*)   );

void organic_reaction_derivs(
	double t,
	double C[],
	double dCdt[]	);

void powell (
	double p[], 
	double **xi, 
	int n, 
	double ftol, 
	int *iter, 
	double *fret, 
	double (*func)( double [] )	);

double qtrap_energy(
	double (*func)(double),
	double a,
	double b   );

double qtrap_time(
	double (*func)(double),
	double a,
	double b	);

double qtrap_time_2(
	double (*func)(double),
	double a,
	double b	);

void recreate_distance_grid ( void ); 

void rkck(
	double y[],
	double dydx[],
	int n,
	double x,
	double h,
	double yout[],
	double yerr[],
	void (*derivs)(double, double*, double*)	);

void rkqs(
	double y[],
	double dydx[],
	int n,
	double *x,
	double htry,
	double eps,
	double yscal[],
	double *hdid,
	double *hnext,
	void (*derivs)(double, double*, double*)	);

double rtbis( 
	double x1, 
	double x2, 
	double xacc,
	double (*func)(double)	);

void scale_reaction_products( void );

void scale_reaction_products_old( void );

void score(
	double x2,
	double T[],
	double f[]  );

void set_up_distance_grid( void );

void set_up_reaction_energies( void );

void shoot(
	int n,
	double v[],
	double f[]  );

void spline(double x[], 
	double y[], 
	int n, 
	double yp1, 
	double ypn, 
	double y2[]	);
	
void splint(
	double xa[], 
	double ya[], 
	double y2a[], 
	int n, 
	double x, 
	double *y	);

void sum_reactions_find_max_rate ( 
	double time, 
	double time_step, 
	int i_grid_base );

void take_reaction_step ( 
	double time,
	double time_step );

void take_temperature_step ( 
	double time,
	double time_step );

void thermal_parameters (
	double temperature );

double time_integral ( 
	double time   );
	
double time_integral_2 ( 
	double time   );
	
void write_out_profile (
	double end_time );

void write_out_time_series (
	double time,
	bool print_label );

#endif /* _SILLBURP_H */
