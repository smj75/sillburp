/*
 * sillburp.cc
 *
 * Stephen M Jones, Birmingham
 * January 2010 v1.0
 *   Solidification time and thermal history after solidification.
 *   Integration of maturation scheme using Numerical Recipes routines.
 * June 2010 v1.1
 *   Fixed bug in 3D integration of reaction products over sill aureole.
 *   Simple venting model. 
 * October 2011 v2.0 
 *   Finite difference integration using analytical solution for 
 *   constant temperature within each timestep.
 *   Burial history pre-intrusion.
 *   Intrusion can be sill or dyke.  
 *   Fixed bug in individual reaction energy parameterization.
 * March 2012 v2.1
 *   Can read in pre-existing temperature history.  
 * October 2012 v3.0
 *   Compares observed and calculated vitrinite reflectance profiles.
 * Autumn 2016 v4.0
 *   Time stepping improved.
 *   Tracks peak reaction rate and records temperature etc there.
 * July 2017 v4.1
 *   Improve reaction rate tracking.
 *   Improve aureole thicknesses tracking.
 * September 2017 v4.2
 *   Read kinetic parameters from command line.
 *   
 */

/* 
	Calculate history of methane generated in and expelled from
	the thermal aureole of a magma sheet using a thermal model 
	and a kinetic maturation model.
*/

/* Parameter and funtion definitions common to sillburp & lipburp */

#include <burp.h>

/* Parameter and funtion definitions for sillburp only */

#include <sillburp.h>

/* Variable definitions and initial values where appropriate */

bool compare_vitrinite_data = false;
bool constant_temperature_reaction = false;
bool do_labile_reaction = false;
bool do_refractory_reaction = false;
bool do_vitrinite_reaction = false;
bool print_aureole = false;
bool print_flux = false;
bool print_mass_expelled = false;
bool print_mass_generated = false;
bool print_misfit = false;
bool print_profile = false;
bool print_reactions = false;
bool print_solidification_time = false;
bool print_temperature = false;
bool print_time_series = false;
bool print_time_series_at_point = false;
bool print_vitrinite = false;
bool scale_to_whole_sill = false;
bool sheet_is_dyke = false;
bool subtract_mass_at_injection = true;
bool use_analytical_solution = false;
bool vent_methane_instantaneously = false;
bool use_output_time_step = false;
bool use_output_time_step_power = false;
bool use_window = false;
int grid_point_centre;
int grid_point_wall;
int grid_points_in_sheet = 5;
int grid_points_in_country_rock = 50;
int integral_grid_point;
int i_integrate_down_to;
int i_max_rate[n_reactions+1] = {0, 0, 0, 0, 0};
int kmax;   // FOR ODEINT
int kount;   // FOR ODEINT
int latent_heat_scheme;
int n_approx_reactions[n_reactions+1] = {0, 7, 21, 55, 7};
int n_approx_max = 0;
int n_aureole;
int n_cols_starting = 7; 
int N_MAX_TIME_SERIES_ARRAY;
int n_grid_points;
int n_print = 0;
int reaction;
double ambient_temperature = 0.0;   //  DEG C
double **approx_reaction_energies;
double *aureole_dless;
double **aureole_model;
double **aureole_model_prev;
double **aureole_rates;
double *aureole_vitrinite;
double aureole_window = 1.0; // 
double begin_constant_temperature;
double burial_rate = 0.01e-3;	// M/YR
//double burial_heating_rate = 1.0e-6;	// degC/YR
double CLOSE_TO_ZERO = 1.0e-30;
double constant_temperature;
double **data_aureole;
double **data_starting;
double density_country_rock = 2400.0 * 1.0e-9; //  MT/M^3
double depth_of_sheet = 1000.0; // M
double *dist;
double dist_increment;
//double dist_increment_sq;
double dless_solidification_boundary;
double dxsav;   // STOREAGE SPACE FOR ODEINT
double end_time = 0.0;
double energy_bound_factor = 2.0;
double energy_in_time_integral;
double energy_interval = 2000.0;
double erfc_dless_solidification_boundary;
double first_time_to_output = 0.0;
	double flux_from_labile;
	double flux_from_refractory;
	double flux_from_oil;
	double flux_methane_generated;
	double flux_methane_vented;
	double flux_from_labile_new;
	double flux_from_refractory_new;
	double flux_from_oil_new;
	double flux_methane_generated_new;
	double flux_methane_vented_new;
double geothermal_gradient = 30.0e-3;	// C/M
double heat_ratio;
double heating_rate_from_burial;
double integral_time;
double latent_heat; //   J/kg
double magma_temperature = 1100; // deg C
double mass_fraction_labile = 0.04;
double mass_fraction_refractory = 0.01;
	double mass_from_labile;
	double mass_from_refractory;
	double mass_from_oil;
	double mass_from_labile_dless_aureole;
	double mass_from_refractory_dless_aureole;
	double mass_from_oil_dless_aureole;
	double mass_methane_generated;
	double mass_from_labile_at_injection = 0.0;
	double mass_from_refractory_at_injection = 0.0;
	double mass_from_oil_at_injection = 0.0;
	double mass_methane_generated_at_injection = 0.0;
	double mass_from_labile_previously = 0.0;
	double mass_from_refractory_previously = 0.0;
	double mass_from_oil_previously = 0.0;
	double mass_methane_generated_previously = 0.0;
	double mass_methane_vented;
	double max_mass_methane_vented;
	double max_flux_methane_vented;
double new_sheet_thickness = 0.0;
double next_time_to_output;
double output_time_power;
double output_time_step;
double output_time_step_power;
double *pos_outer_edge_lower_aureole;
double *pos_inner_edge_lower_aureole;
double *pos_max_rate_lower_aureole;
double *pos_outer_edge_upper_aureole;
double *pos_inner_edge_upper_aureole;
double *pos_max_rate_upper_aureole;
double ***progress_of_reactions;
double **progress_at_injection;
double ***rate_of_reactions;
double reflectance_min = 0.2; // Corresponds to EasyRo min
double reflectance_max = 4.7; // Corresponds to EasyRo max
double sheet_half_thickness = 0.0; // m
double sheet_thickness = 0.0; // m
double sill_radius;  
double sill_radius_squared;
double solidification_time;
double solid_crust;
double specific_heat; //  /kg/K
double stable_burial_time_step = 0.0;
double stable_conduction_time_step = 0.0;
double start_time;
double surface_temperature = 0.0;
double *temp_max_rate_lower_aureole;
double *temp_max_rate_upper_aureole;
double *temperature;
double *temperature_increment;
double thermal_diffusivity; // m^2/s
double time_scale = S_PER_YR;
double *time_series_flux_upper;
double *time_series_flux;
double *time_series_mass_upper;
double *time_series_mass;
double time_so_far;
double time_step_print = 50.0;
double time_step_after = 50.0;
double time_step_before = 1.0e6;
double vent_formation_time;
double vent_decay_time;
double vent_initial_flux;
double vent_background_flux;
double vent_methane_saturation;
double *vit_outer_edge_lower_aureole;
double *vit_inner_edge_lower_aureole;
double *vit_max_rate_lower_aureole;
double *vit_outer_edge_upper_aureole;
double *vit_inner_edge_upper_aureole;
double *vit_max_rate_upper_aureole;
double window_dless_width = 0.0;
double *xp;   // STOREAGE SPACE FOR ODEINT
double **yp;   // STOREAGE SPACE FOR ODEINT
const char *time_unit = "yr";
char file_aureole[BUFSIZ];
char file_starting_temperature[BUFSIZ];
FILE *fp_aureole = NULL; // Vitrinite reflectance profile 
FILE *fp_starting_temperature = NULL; // Starting temperature history 




/* Main program */

int main(int argc, char **argv)
{

/* Local variables */

	bool reactions_finished = false;
	int i_grid_base;
	int i_grid_top;
	double prev_time_to_output;
	double time;
	double time_end_starting;
	double time_step;
	double time_step_starting;
	double time_step_for_integrate_aureole;

/* Banner */

	cerr << "\n";
	cerr << " * * * sillburp * * *\n";
	cerr << " (c) Stephen M Jones\n";
	cerr << " v4.2 September 2017\n\n";

/* Read controlling parameters from command line 
   Write them to first line of output file
 */

	get_user_parameters( argc, argv );
	for (int i=0; i<argc; i++) cout << argv[i] << " ";
	cout << "\n";

/* SCALE TIME */

	time_step_after *= time_scale;
	time_step_before *= time_scale;
	time_step_print *= time_scale;
	end_time *= time_scale;
	output_time_step *= time_scale;

/* SET UP ENERGY GRID */

	set_up_reaction_energies();

/* READ IN PRE-EXISTING TEMPERATURE HISTORY IF SUPPLIED */

	if ( fp_starting_temperature ) {
		data_starting = dmatrix(0, MAX_FIELDS, 0, MAX_RECORDS);
		cerr << " Reading STARTING TEMPERATURE & REACTION data file...\n";
		if ( read_data_file (fp_starting_temperature, &n_cols_starting, &n_grid_points, data_starting, &time_end_starting, &time_step_starting, 1 ) ) {
			cerr << " Fatal error whilst reading data\n";
			exit (0);
		}
		if ( n_cols_starting != 7 ) {
			cerr << " ERROR: Data file not produced by sillburp -Pr: wrong number of columns\n";
			exit(0);
		}

/* SET UP DISTANCE GRID */

	  recreate_distance_grid();
	} else {
	  set_up_distance_grid();
	}
	sill_radius_squared = sill_radius * sill_radius;

/* READ IN AUREOLE FILE IF SUPPLIED */

	double dummy1, dummy2;
	int n_cols_aureole = 3;
	if ( fp_aureole ) {
		data_aureole = dmatrix(0, MAX_FIELDS, 0, MAX_RECORDS);
		cerr << " Reading AUREOLE (Vitrinite Reflectance) data file...\n";
		if ( read_data_file (fp_aureole, &n_cols_aureole, &n_aureole, data_aureole, &dummy1, &dummy2, 0 ) ) {
			cerr << "   Fatal error whilst reading data\n";
			exit (0);
		}
		if ( n_cols_aureole != 3 ) {
			cerr << "   ERROR: Data file has wrong number of columns\n";
			cerr << "          " << n_aureole << " cols read, 3 expected\n";
			exit(0);
		}
		cerr << "    " << n_aureole << " vitrinite reflectance data successfully read\n";
	}

/* Find how long it will take for sheet to solidify 
	 AND SET UP CONSTANTS FOR ANALYTICAL SOLIDIFICATION CALCULATION 
 */

	if ( !fp_starting_temperature ) 
		find_solidification_time();
	else
		start_time = 0.0;	
	if ( print_solidification_time ) exit(0);

// Set up storeage arrays.  

		aureole_dless = dvector(1,n_reactions);
		pos_inner_edge_lower_aureole = dvector(1,n_reactions);
		pos_max_rate_lower_aureole = dvector(1,n_reactions);
		pos_outer_edge_lower_aureole = dvector(1,n_reactions);
		pos_inner_edge_upper_aureole = dvector(1,n_reactions);
		pos_max_rate_upper_aureole = dvector(1,n_reactions);
		pos_outer_edge_upper_aureole = dvector(1,n_reactions);
		progress_at_injection = dmatrix(1,n_reactions,1,n_grid_points);
		temp_max_rate_lower_aureole = dvector(1,n_reactions);
		temp_max_rate_upper_aureole = dvector(1,n_reactions);
		vit_inner_edge_lower_aureole = dvector(1,n_reactions);
		vit_max_rate_lower_aureole = dvector(1,n_reactions);
		vit_outer_edge_lower_aureole = dvector(1,n_reactions);
		vit_inner_edge_upper_aureole = dvector(1,n_reactions);
		vit_max_rate_upper_aureole = dvector(1,n_reactions);
		vit_outer_edge_upper_aureole = dvector(1,n_reactions);
	
//	THESE ALREADY SET UP IF PREVIOUS TEMPERATURE HISTORY WAS READ IN

	if ( !fp_starting_temperature ) {
		progress_of_reactions = d3tensor(1,n_approx_max,1,n_reactions,1,n_grid_points);
		rate_of_reactions = d3tensor(1,n_approx_max,1,n_reactions,1,n_grid_points);		
		aureole_model = dmatrix(1,n_reactions,1,n_grid_points);
		aureole_model_prev = dmatrix(1,n_reactions,1,n_grid_points);
		aureole_rates = dmatrix(1,n_reactions,1,n_grid_points);
		aureole_vitrinite = dvector(1,n_grid_points);
		temperature = dvector(1,n_grid_points);
		temperature_increment = dvector(1,n_grid_points);
	}

//	Initiate reaction progress (= concentration of product) at zero.
//	Initiate temperature at absolute zero. 
	
		for ( int i_grid=1; i_grid<=n_grid_points; i_grid++ ) {
			temperature[i_grid] = ABSOLUTE_ZERO;
			temperature_increment[i_grid] = 0.0;
			for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ )
				for ( int i_approx=1; i_approx<=n_approx_max; i_approx++ )
					if ( abs(dist[i_grid]) <= sheet_half_thickness ) {
						progress_of_reactions[i_approx][i_reaction][i_grid] = 0.0;
						rate_of_reactions[i_approx][i_reaction][i_grid] = 0.0;
					}
		}

/* GRID POINT AT TOP OF SILL */

	if (sheet_is_dyke) {
			i_grid_top = grid_points_in_sheet+1;
	} else {
		for ( int i=1; i<=n_grid_points; i++) {
//		double distance = (sheet_is_dyke) ? dist[i] : dist[i]-depth_of_sheet;
			double distance = dist[i]-depth_of_sheet;
			if ( distance == -sheet_half_thickness ) {
				i_grid_top = i;
				break;
			}
		}
	}

/* GRID POINT AT BOTTOM OF SILL */

	for ( int i=i_grid_top; i<=n_grid_points; i++) {
		double distance = (sheet_is_dyke) ? dist[i] : dist[i]-depth_of_sheet;
		if ( distance == sheet_half_thickness ) {
			i_grid_base = i;
			break;
		}
	}

/* INTEGRATING ONLY DOWN TO A GIVEN DEPTH  
   IN THE CASE OF A SILL */

	if ( use_window && sheet_is_dyke ) {
		cerr << "   -J window not used for a dyke\n";
		use_window = false;
	}
	if ( use_window && aureole_window >= dist[n_grid_points] ) {
		cerr << "   -J window not used: max grid depth less than top window\n";
		use_window = false;
	}
	if ( use_window && !compare_vitrinite_data) {
		i_integrate_down_to = 0;
		do {
			i_integrate_down_to++;
		} while ( dist[i_integrate_down_to] < aureole_window );
	} else {
		i_integrate_down_to = n_grid_points;
	}

/* GRID NODES AT SHEET EDGES */

	if (sheet_is_dyke) {
		cerr << "   Grid node " << i_grid_top << " at dyke edge, " << n_grid_points << " max\n";
	} else {
		cerr << "   Grid node " << i_grid_top << " at sill top, " 
				<< i_grid_base << " at base, " << n_grid_points << " max\n";
		cerr << "   Max depth is " << dist[n_grid_points] << " m, integrating to " << dist[i_integrate_down_to] << " m\n";
	}	

/* Print list of what is in each output column */

	if ( print_time_series ) {
		write_out_time_series ( 0.0, true );
	}

/* Prepare for time loop */

	int n_steps = 0;
	double background_aureole_generation_rate = 0.0;
	time = start_time;
	cerr << " Beginning thermal maturation calculation at " << start_time/time_scale << " yr, aiming for " << end_time/time_scale << " yr\n";
//  time_step = find_time_step( time );
	next_time_to_output = first_time_to_output;
	output_time_power = 0.0;
	output_time_power = log10( time_scale )-output_time_step_power;
	
/* Start of time loop */
	
	while ( time < end_time ) {
		n_steps ++;
		
/* Time stepping */

  	time_step = find_time_step( time );
		time += time_step;
		
/* Time of magma injection: Copy 1d reaction profile across the whole grid. */

		if ( time == 0.0 && sheet_is_dyke ) {
			for ( int i_grid=2; i_grid<=n_grid_points; i_grid++ )
				for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ )
					for ( int i_approx=1; i_approx<=n_approx_max; i_approx++ )
						if ( abs(dist[i_grid]) >= sheet_half_thickness )
							progress_of_reactions[i_approx][i_reaction][i_grid] = progress_of_reactions[i_approx][i_reaction][1];				
			cerr << "   Copied 1d profile to 2d profile at injection time, t=0\n";
		}

/* Calculate temperature profile at this time step */

    take_temperature_step ( time, time_step );

/* Advance reactions */

    take_reaction_step ( time, time_step );
    
/* Sum the individual reactions; find points of max reaction rates */

		sum_reactions_find_max_rate ( time, time_step, i_grid_base );

/* Print time series at a specific point on the distance grid */

    if ( print_time_series ) {
			if ( time >= next_time_to_output ) {

/* Identify reaction aureole edges and position of max reaction rates */

				find_aureole_edges( i_grid_top, i_grid_base, i_integrate_down_to, time_step );

/* Find d'less aureole thickness (scaled by sheet thickness) from total mass */ 

				integrate_dless_aureole ( i_grid_top, i_grid_base, i_integrate_down_to );

/* Integrate over approximate reactions and scale */

		    scale_reaction_products();
		
/* Integrate over aureole */
/* CAN WE DO WITHOUT THIS IF PRINTING OUT CERTAIN SORTS OF TIME SERIES? */

				time_step_for_integrate_aureole = next_time_to_output - prev_time_to_output;
				integrate_aureole ( i_grid_top, i_grid_base, i_integrate_down_to, time, time_step_for_integrate_aureole );
				prev_time_to_output = next_time_to_output;

/* Write out */

				write_out_time_series ( time, false );

/* Next time step */
				
    		if ( use_output_time_step ) {
    			if ( use_output_time_step_power ) {
						output_time_power += output_time_step_power;
						next_time_to_output = pow(10.0, output_time_power);
					} else {
						next_time_to_output += output_time_step;
					}
					if ( time > 0.0 ) cerr << "written " << "\n";
					if ( time < end_time ) cerr << "   Next time to write out is " << next_time_to_output/time_scale << " " << time_unit << " ... ";
				} else {
					next_time_to_output = time;
				} 
			}
    }
        
/* Test whether reaction has finished (i.e. generation is within
   and order of magnitude of the background rate)
 */

/*	
		if ( time < 0.0 ) background_aureole_generation_rate = (mass_methane_generated - mass_methane_generated_previously) / time_step;
		if ( time > 0.0 && abs(mass_methane_generated-mass_methane_generated_previously) < 10.0*background_aureole_generation_rate*time_step ) {
			cerr << "   Aureole reactions finished at " << time/time_scale << " yr\n";
			time = end_time;
		}
*/


/* End of TIME LOOP */

	} 

/* PROGRAM ENDS HERE IF WRITING OUT A TIME SERIES */

	if ( print_time_series ) {
		cerr << "   Sheet radius: " << sill_radius*1.0e-3 << " km"<<  "\n";
		cerr << "   Total methane generated: " << mass_methane_generated << " Mt"<<  "\n";
		cerr << "\n";
		exit(0);
	}

/* WRITE OUT PROFILE */

	if ( print_profile ) {
		if ( !print_misfit ) write_out_profile ( end_time );
		
/* CALCULATE DATA MISFIT IF AUREOLE DATA SUPPLIED */

		if ( fp_aureole ) find_aureole_misfit (); 

/* PROGRAM ENDS HERE IF PRINTING OUT A PROFILE */

		cerr << "\n";
		exit (0);
	}

/* INTEGRATE OVER DISTANCE FROM SILL */
/*
	double width = 250.0;
	cerr << " Testing mass integration at aureole " << width << " m from centre line\n";
	for ( int i=1; i<=n_grid_points; i++) {
		double distance = (sheet_is_dyke) ? dist[i] : dist[i]-depth_of_sheet;
		if ( abs(distance) <= width ) {
			progress_of_reactions[1][LABILE][i] = 1.0;
			progress_of_reactions[1][REFRACTORY][i] = 1.0;
			progress_of_reactions[1][OIL][i] = 1.0;
		} else {
			progress_of_reactions[1][LABILE][i] = 0.0;
			progress_of_reactions[1][REFRACTORY][i] = 0.0;
			progress_of_reactions[1][OIL][i] = 0.0;
		}
	}
*/


/* Integrate methane mass generated throughout aureole */

	integrate_aureole ( i_grid_top, i_grid_base, i_integrate_down_to, time, time_step );

/* WRITE METHANE VENTED */


/* END OF MAIN */

  exit(0);
}


void integrate_dless_aureole ( int i_grid_top, int i_grid_base, int i_integrate_down_to )
{

// INTEGRATE AUREOLE DEFINED BY REACTION PROGRESSES
// TO FIND TOTAL WITH OF DIMENSIONLESS AUREOLE
// THE AUREOLE IS DEFINED RELATIVE TO REACTION PROGRESS AT INJECTION
// WINDOWING IS NOT DONE

  int i_end;
  int i_start;
  double symmetry_factor;

	symmetry_factor = sheet_is_dyke ? 2.0 : 1.0;
	for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ ) aureole_dless[i_reaction] = 0;

/* Integrate methane produced in aureole above sill, or adjacent to dyke */

	if ( sheet_is_dyke ) {
		i_start = i_grid_top + 1;
		i_end = n_grid_points; 
	} else {
		i_start = 2;
		i_end = i_grid_top;
	}
	for ( int i=i_start; i<=i_end; i++) {
		for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ ) {
			aureole_dless[i_reaction] += symmetry_factor * 0.5*( (aureole_model[i_reaction][i] - progress_at_injection[i_reaction][i]) + (aureole_model[i_reaction][i-1] - progress_at_injection[i_reaction][i])) * ( dist[i] - dist[i-1] );
		}
	}

	
/* Integrate methane produced in aureole below sill only */

	if ( !sheet_is_dyke ) {
		for ( int i=i_grid_base+1; i<=n_grid_points; i++) {
			for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ ) {
			aureole_dless[i_reaction] += symmetry_factor * 0.5*( (aureole_model[i_reaction][i] - progress_at_injection[i_reaction][i]) + (aureole_model[i_reaction][i-1] - progress_at_injection[i_reaction][i])) * ( dist[i] - dist[i-1] );
			}
		}
	}

/* Make dimensionless by scaling using sheet thickness */

	for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ ) aureole_dless[i_reaction] /= sheet_thickness;
	
/* Total gas masses from this method */
	
	mass_from_oil_dless_aureole = aureole_dless[OIL] * sheet_thickness * (1.0-mass_frac_labile_to_gas) * mass_fraction_labile * carbon_to_methane * density_country_rock * PI;
	mass_from_labile_dless_aureole = aureole_dless[LABILE] * sheet_thickness  * mass_frac_labile_to_gas * mass_fraction_labile * carbon_to_methane * density_country_rock * PI;
	mass_from_refractory_dless_aureole = aureole_dless[REFRACTORY] * sheet_thickness  * mass_fraction_refractory * carbon_to_methane * density_country_rock * PI;


}



void integrate_aureole ( int i_grid_top, int i_grid_base, int i_integrate_down_to, double time, double time_step )
{

  int i_end;
  int i_start;
  double symmetry_factor;

	flux_from_labile_new = 0.0;
	flux_from_refractory_new = 0.0;
	flux_from_oil_new = 0.0;
	flux_methane_generated_new = 0.0;
	mass_from_labile = 0.0;
	mass_from_refractory = 0.0;
	mass_from_oil = 0.0;
	mass_methane_generated = 0.0;
	symmetry_factor = sheet_is_dyke ? 2.0 : 1.0;

/* Integrate methane produced in aureole above sill, or adjacent to dyke */

	if ( sheet_is_dyke ) {
		i_start = i_grid_top + 1;
		i_end = n_grid_points; 
	} else {
		i_start = 2;
		i_end = ( use_window && i_integrate_down_to<i_grid_top ) ? i_integrate_down_to : i_grid_top;
	}
	for ( int i=i_start; i<=i_end; i++) {
			flux_from_labile_new += symmetry_factor * 0.5*( aureole_rates[LABILE][i] + aureole_rates[LABILE][i-1]) * ( dist[i] - dist[i-1] );
			flux_from_refractory_new += symmetry_factor * 0.5*( aureole_rates[REFRACTORY][i] + aureole_rates[REFRACTORY][i-1]) * ( dist[i] - dist[i-1] );
			flux_from_oil_new += symmetry_factor * 0.5*( aureole_rates[OIL][i] + aureole_rates[OIL][i-1]) * ( dist[i] - dist[i-1] );
			mass_from_labile += symmetry_factor * 0.5*( aureole_model[LABILE][i] + aureole_model[LABILE][i-1]) * ( dist[i] - dist[i-1] );
			mass_from_refractory += symmetry_factor * 0.5*( aureole_model[REFRACTORY][i] + aureole_model[REFRACTORY][i-1]) * ( dist[i] - dist[i-1] );
			mass_from_oil += symmetry_factor * 0.5*( aureole_model[OIL][i] + aureole_model[OIL][i-1]) * ( dist[i] - dist[i-1] );
	}
	
/* Integrate methane produced in aureole below sill only */

	if ( !sheet_is_dyke && i_integrate_down_to > i_grid_base+1 )
	for ( int i=i_grid_base+1; i<=i_integrate_down_to; i++) {
		flux_from_labile_new += 0.5*( aureole_rates[LABILE][i] + aureole_rates[LABILE][i-1]) * ( dist[i] - dist[i-1] );
		flux_from_refractory_new += 0.5*( aureole_rates[REFRACTORY][i] + aureole_rates[REFRACTORY][i-1]) * ( dist[i] - dist[i-1] );
		flux_from_oil_new += 0.5*( aureole_rates[OIL][i] + aureole_rates[OIL][i-1]) * ( dist[i] - dist[i-1] );
		mass_from_labile += 0.5*( aureole_model[LABILE][i] + aureole_model[LABILE][i-1]) * ( dist[i] - dist[i-1] );
		mass_from_refractory += 0.5*( aureole_model[REFRACTORY][i] + aureole_model[REFRACTORY][i-1]) * ( dist[i] - dist[i-1] );
		mass_from_oil += 0.5*( aureole_model[OIL][i] + aureole_model[OIL][i-1]) * ( dist[i] - dist[i-1] );
	}
	
	flux_from_labile_new *= PI;
	flux_from_refractory_new *= PI;
	flux_from_oil_new *= PI;
	mass_from_labile *= PI;
	mass_from_refractory *= PI;
	mass_from_oil *= PI;
	
	if ( scale_to_whole_sill ) {
		flux_from_labile_new *= sill_radius_squared;
		flux_from_refractory_new *= sill_radius_squared;
		flux_from_oil_new *= sill_radius_squared;
		mass_from_labile *= sill_radius_squared;
		mass_from_refractory *= sill_radius_squared;
		mass_from_oil *= sill_radius_squared;
	}
	
	flux_methane_generated_new = flux_from_labile_new + flux_from_refractory_new + flux_from_oil_new;
	mass_methane_generated = mass_from_labile + mass_from_refractory + mass_from_oil;
	
/* Store masses at time of injection */
	
	if ( time == 0 ) {
		mass_from_labile_at_injection = mass_from_labile;
		mass_from_refractory_at_injection = mass_from_refractory;
		mass_from_oil_at_injection = mass_from_oil;
		mass_methane_generated_at_injection = mass_methane_generated;
		if ( print_time_series )
			cerr << "   Reaction progresses at injection time: labile, " <<
			mass_from_labile_at_injection << "; refractory, " <<
			mass_from_refractory_at_injection << "; oil, " <<
			mass_from_oil_at_injection << "\n";
	}
	
/* Calculate production rates */

	flux_from_labile = (mass_from_labile - mass_from_labile_previously) / (time_step/time_scale);
	flux_from_refractory = (mass_from_refractory - mass_from_refractory_previously) / (time_step/time_scale);
	flux_from_oil = (mass_from_oil - mass_from_oil_previously) / (time_step/time_scale);
	flux_methane_generated = (mass_methane_generated - mass_methane_generated_previously) / (time_step/time_scale);


/* Store masses to calculate rate at next time step */

	mass_from_labile_previously = mass_from_labile;
	mass_from_refractory_previously = mass_from_refractory;
	mass_from_oil_previously = mass_from_oil;
	mass_methane_generated_previously = mass_methane_generated;
/*
	cerr 	<< "Fluxes from labile: " <<	flux_from_labile << "\t" << flux_from_labile_new << "\n"
				<< "Fluxes from refractory: " <<	flux_from_refractory << "\t" << flux_from_refractory_new << "\n"
				<< "Fluxes from oil: " <<	flux_from_oil << "\t" << flux_from_oil_new << "\n"
				<< "Total flux: " <<	flux_methane_generated << "\t" << flux_methane_generated_new << "\n";
*/
flux_from_labile = flux_from_labile_new;
flux_from_refractory = flux_from_refractory_new;
flux_from_oil = flux_from_oil_new;
flux_methane_generated = flux_methane_generated_new;



/* MAXIMUM METHANE THAT COULD HAVE BEEN VENTED DURING THIS TIME PERIOD */

	if ( print_mass_expelled ) {
		if ( time/time_scale <= vent_formation_time ) {
			max_mass_methane_vented = 0.0;
		} else {
			double vent_duration = end_time/time_scale - vent_formation_time;
			max_mass_methane_vented = vent_methane_saturation * ( vent_background_flux * vent_duration + ( vent_initial_flux - vent_background_flux ) * vent_decay_time * ( 1.0 - exp( -vent_duration / vent_decay_time ) ) );
			max_flux_methane_vented = vent_methane_saturation * ( vent_background_flux + ( vent_initial_flux - vent_background_flux ) * exp( -vent_duration / vent_decay_time ) );
		}
		mass_methane_vented = FMIN( max_mass_methane_vented, mass_methane_generated );
	}

}


void find_aureole_edges ( int i_grid_top, int i_grid_base, int i_integrate_down_to, double time_step )
{

/* 
	Find edge of reaction zone with reference to progress of reactions
	at time of injection (stored).  
	Uses constant "def_edge_reaction_zone" (=del) defined in sillburp.h
	Outer reaction zone edges taken to be when progress = progress at injection + del.
	Inner reaction zone edges taken to be when progress = 1 - del.
*/

	bool found_upper_inner;
	bool found_upper_outer;
	bool found_upper_rate;
	bool found_lower_inner;
	bool found_lower_outer;
	bool found_lower_rate;
	int i_end;
	int i_grid;
	int i_outer;
	int i_reaction;
	double **check_rates;
	double distance;
	double interpolation_factor;
	double max_rate;
	double progress;
	double rate;

	check_rates = dmatrix(1, n_reactions, 1, n_grid_points);
	for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ ) {
		found_upper_outer = false;
		found_upper_inner = false;
		found_lower_outer = false;
		found_lower_inner = false;
		pos_outer_edge_upper_aureole[i_reaction] = 0.0;
		pos_max_rate_upper_aureole[i_reaction] = 0.0;
		pos_inner_edge_upper_aureole[i_reaction] = 0.0;
		pos_inner_edge_lower_aureole[i_reaction] = 0.0;
		pos_max_rate_lower_aureole[i_reaction] = 0.0;
		pos_outer_edge_lower_aureole[i_reaction] = 0.0;
		temp_max_rate_upper_aureole[i_reaction] = 0.0;
		temp_max_rate_lower_aureole[i_reaction] = 0.0;
		vit_outer_edge_upper_aureole[i_reaction] = 0.0;
		vit_max_rate_upper_aureole[i_reaction] = 0.0;
		vit_inner_edge_upper_aureole[i_reaction] = 0.0;
		vit_inner_edge_lower_aureole[i_reaction] = 0.0;
		vit_max_rate_lower_aureole[i_reaction] = 0.0;
		vit_outer_edge_lower_aureole[i_reaction] = 0.0;
		
/* Sum approximate reactions to find overall reaction rate */

		for (int i_grid=1; i_grid<=n_grid_points; i_grid++) {
//			check_rates[i_reaction][i_grid] = time_scale * get_rate( i_reaction, i_grid );
//			aureole_rates[i_reaction][i_grid] = (aureole_model[i_reaction][i_grid] - aureole_model_prev[i_reaction][i_grid])/(time_step/time_scale);
			aureole_rates[i_reaction][i_grid] = time_scale * get_rate( i_reaction, i_grid );

/* Scale vitrinite reaction to %Ro */

			if ( i_reaction == VITRINITE )
//				progress = get_progress(  i_reaction, i_grid );
				aureole_vitrinite[i_grid] = reflectance_min + aureole_model[i_reaction][i_grid] * (reflectance_max-reflectance_min);
		}
		
/* Upper aureole (only for sill) */

		if ( !sheet_is_dyke ) {

			i_end = ( use_window && i_integrate_down_to<i_grid_top ) ? i_integrate_down_to+1 : i_grid_top;
			max_rate = 0.0;
			for ( i_grid=i_end-1; i_grid>=1; i_grid-- ) {
/*
				double check = aureole_model[i_reaction][i_grid];
				distance = dist[i_grid]-depth_of_sheet;
				cerr << distance << "\t" << i_reaction << "\t" << check << "\t" << 
					progress_at_injection[i_reaction][i_grid]+def_edge_reaction_zone[i_reaction] << "\t" <<
					aureole_rates[i_reaction][i_grid]  << "\t" <<
					check_rates[i_reaction][i_grid] << "\t" <<
					aureole_vitrinite[i_grid] << "\n";
*/				
/* Find inner edge of upper aureole reaction zone */

				if ( !found_upper_inner && aureole_model[i_reaction][i_grid] < (1.0-def_edge_reaction_zone[i_reaction]) ) {
					distance = dist[i_grid]-depth_of_sheet;
					interpolation_factor = (1.0-def_edge_reaction_zone[i_reaction] - aureole_model[i_reaction][i_grid]) / (aureole_model[i_reaction][i_grid+1] - aureole_model[i_reaction][i_grid]);
					pos_inner_edge_upper_aureole[i_reaction] = distance + (dist[i_grid+1]-dist[i_grid]) * interpolation_factor + sheet_half_thickness;
					vit_inner_edge_upper_aureole[i_reaction] = aureole_vitrinite[i_grid] + (aureole_vitrinite[i_grid+1]-aureole_vitrinite[i_grid]) * interpolation_factor;
					found_upper_inner = true;
//					cerr << "pos_inner_edge_upper_aureole: " << pos_inner_edge_upper_aureole[i_reaction] << "\t" << vit_inner_edge_upper_aureole[i_reaction] << "\n";
				}

/* Find maximum reaction rate in upper aureole */

				if ( aureole_rates[i_reaction][i_grid] >= max_rate ) {
					distance = dist[i_grid]-depth_of_sheet;
					max_rate = aureole_rates[i_reaction][i_grid];
					pos_max_rate_upper_aureole[i_reaction] = distance + sheet_half_thickness;
					temp_max_rate_upper_aureole[i_reaction] = temperature[i_grid];
					vit_max_rate_upper_aureole[i_reaction] = aureole_vitrinite[i_grid];
				}
				
/* Find outer edge of upper aureole reaction zone */

				if ( !found_upper_outer && aureole_model[i_reaction][i_grid] < (progress_at_injection[i_reaction][i_grid] + def_edge_reaction_zone[i_reaction]) ) {
					distance = dist[i_grid]-depth_of_sheet;
					interpolation_factor = (progress_at_injection[i_reaction][i_grid+1] + def_edge_reaction_zone[i_reaction] - aureole_model[i_reaction][i_grid]) / (aureole_model[i_reaction][i_grid+1] - aureole_model[i_reaction][i_grid]);
					pos_outer_edge_upper_aureole[i_reaction] = distance + (dist[i_grid+1]-dist[i_grid]) * interpolation_factor + sheet_half_thickness;
					vit_outer_edge_upper_aureole[i_reaction] = aureole_vitrinite[i_grid] + (aureole_vitrinite[i_grid+1]-aureole_vitrinite[i_grid]) * interpolation_factor;
					found_upper_outer = true;
					i_outer = i_grid;
//					cerr << "pos_outer_edge_upper_aureole: " << pos_outer_edge_upper_aureole[i_reaction] << "\t" << vit_outer_edge_upper_aureole[i_reaction] << "\n";
//					break;
				}

			} // End of grid loop
//			cerr << "pos, temp and vit_max_rate_upper_aureole: " << pos_max_rate_upper_aureole[i_reaction] << "\t" << temp_max_rate_upper_aureole[i_reaction] << "\t" << vit_max_rate_upper_aureole[i_reaction] << "\n";
		}
						
/* Aureole of dyke, or lower aureole of sill */

		if ( i_integrate_down_to > i_grid_base ) {
			max_rate = 0.0;
			for ( i_grid=i_grid_base+1; i_grid<=i_integrate_down_to; i_grid++ ) {
/*
				double check = aureole_model[i_reaction][i_grid];
				distance = dist[i_grid]-depth_of_sheet;
				cerr << distance << "\t" << i_reaction << "\t" << check << "\t" << progress_at_injection[i_reaction][i_grid]+def_edge_reaction_zone[i_reaction] << "\t" <<
					aureole_rates[i_reaction][i_grid]  << "\t" <<
					check_rates[i_reaction][i_grid] << "\t" <<
					aureole_vitrinite[i_grid] << "\n";
*/
/* Find inner edge of lower aureole reaction zone */

				if ( !found_lower_inner && (aureole_model[i_reaction][i_grid] < (1.0-def_edge_reaction_zone[i_reaction])) ) {
					distance = (sheet_is_dyke) ? dist[i_grid-1] : dist[i_grid-1]-depth_of_sheet;
					interpolation_factor = (1.0-def_edge_reaction_zone[i_reaction] - aureole_model[i_reaction][i_grid-1]) / (aureole_model[i_reaction][i_grid] - aureole_model[i_reaction][i_grid-1]);
					pos_inner_edge_lower_aureole[i_reaction] = distance + (dist[i_grid]-dist[i_grid-1]) * interpolation_factor - sheet_half_thickness;
					vit_inner_edge_lower_aureole[i_reaction] = aureole_vitrinite[i_grid-1] + (aureole_vitrinite[i_grid]-aureole_vitrinite[i_grid-1]) * interpolation_factor;
					found_lower_inner = true;
//					cerr << "pos_inner_edge_lower_aureole: " << pos_inner_edge_lower_aureole[i_reaction] << "\t" << vit_inner_edge_lower_aureole[i_reaction] << "\n";
				}

/* Find maximum reaction rate in lower aureole */

				if ( aureole_rates[i_reaction][i_grid] >= max_rate ) {
					distance = (sheet_is_dyke) ? dist[i_grid] : dist[i_grid]-depth_of_sheet;
					max_rate = aureole_rates[i_reaction][i_grid];
					pos_max_rate_lower_aureole[i_reaction] = distance - sheet_half_thickness;
					temp_max_rate_lower_aureole[i_reaction] = temperature[i_grid];
					vit_max_rate_lower_aureole[i_reaction] = aureole_vitrinite[i_grid];
				}
	
/* Find outer edge of lower aureole reaction zone */

				if ( !found_lower_outer && (aureole_model[i_reaction][i_grid] < (progress_at_injection[i_reaction][i_grid] + def_edge_reaction_zone[i_reaction])) ) {
					distance = (sheet_is_dyke) ? dist[i_grid-1] : dist[i_grid-1]-depth_of_sheet;
					interpolation_factor = (progress_at_injection[i_reaction][i_grid] + def_edge_reaction_zone[i_reaction] - aureole_model[i_reaction][i_grid-1]) / (aureole_model[i_reaction][i_grid] - aureole_model[i_reaction][i_grid-1]);
					pos_outer_edge_lower_aureole[i_reaction] = distance + (dist[i_grid]-dist[i_grid-1]) * interpolation_factor - sheet_half_thickness;
					vit_outer_edge_lower_aureole[i_reaction] = aureole_vitrinite[i_grid-1] + (aureole_vitrinite[i_grid]-aureole_vitrinite[i_grid-1]) * interpolation_factor;
					found_lower_outer = true;
//					cerr << "pos_outer_edge_lower_aureole: " << pos_outer_edge_lower_aureole[i_reaction] << "\t" << vit_outer_edge_lower_aureole[i_reaction] << "\n";
//					break;
				}
				
				if ( found_lower_inner && found_lower_outer && found_lower_rate ) break;
			} // End of grid loop
//			cerr << "pos, temp and vit_max_rate_lower_aureole: " << pos_max_rate_lower_aureole[i_reaction] << "\t" << temp_max_rate_lower_aureole[i_reaction] << "\t" << vit_max_rate_lower_aureole[i_reaction] << "\n";
		}
	} // End of reaction loop
}

void find_aureole_misfit ( void )
{

/* 
	Calculate misfit for each side of sheet separately
	because of the corners in the model profile at the sheet boundary
 */

	cerr << " Calculating data misfit\n";
	
/* 
	Assume that dyke observations are relative to sheet edge 
	and sill observations are relative to sheet centre
*/

	bool distance_rel_to_sheet_edge = false;
	if ( sheet_is_dyke ) distance_rel_to_sheet_edge = true;

/* Aureole model */
	
	int n_model;
	double *model_aureole;
	double *model_distance;
	model_aureole = dvector(1,n_grid_points);
	model_distance = dvector(1,n_grid_points);
	if ( sheet_is_dyke ) {
	
/* Dyke, distance positive and relative to sheet egde */

		if ( distance_rel_to_sheet_edge ) {
			for (int i=grid_points_in_sheet; i<=n_grid_points; i++) {
				model_aureole[i-grid_points_in_sheet] = progress_of_reactions[1][VITRINITE][i];
				model_distance[i-grid_points_in_sheet] = dist[i] - sheet_half_thickness;
			}
			n_model = n_grid_points - grid_points_in_sheet;
		}

/* Sill, distance positive and relative to sheet centre */

	} else {

		if ( !distance_rel_to_sheet_edge ) {
			for (int i=grid_point_centre+grid_points_in_sheet; i<=n_grid_points; i++) {
				model_aureole[i] = progress_of_reactions[1][VITRINITE][i];
				model_distance[i] = dist[i];
//				model_distance[i] = dist[i] - depth_of_sheet - sheet_half_thickness;
			}
			n_model = n_grid_points - grid_point_centre - grid_points_in_sheet;
		}
	}

/* Get model 2nd derivative */

	double *model_2deriv;
	model_2deriv = dvector(1,n_grid_points);
	spline(	model_distance, model_aureole, n_model, 0.0, 0.0, model_2deriv	);
	
/* Determine misfit for dyke or lower side of sill */

	double model_value;
	double misfit = 0.0;
	int n_aureole_used = 0;
//	cerr << "Distance, Dless distance, Obs, Model, Misfit\n";
 cerr << "Window: " << window_dless_width << "\n";
	for (int i=1; i<=n_aureole; i++ ) {
		if ( data_aureole[0][i] < model_distance[1] || data_aureole[0][i] > model_distance[n_model] ) continue;
		if ( distance_rel_to_sheet_edge ) {
			if ( data_aureole[0][i] < window_dless_width*sheet_thickness ) continue;
		} else {
			if ( data_aureole[0][i] < (window_dless_width+0.5)*sheet_thickness ) continue;
		}
		n_aureole_used++;
		splint( model_distance, model_aureole, model_2deriv, n_model, data_aureole[0][i], &model_value	);
		double dif = (data_aureole[1][i] - model_value)/data_aureole[2][i];
		misfit += dif*dif;
		cerr	<< data_aureole[0][i] << "\t"
					<< data_aureole[0][i]*window_dless_width << "\t"
					<< data_aureole[1][i] << "\t"
					<< model_value << "\t"
					<< misfit << "\n";

	}
	
/* Determine misfit from upper side of sill */	

	if ( !sheet_is_dyke ) {
		for (int i=1; i<=grid_point_centre-grid_points_in_sheet; i++) {
			model_aureole[i] = progress_of_reactions[1][VITRINITE][i];
			model_distance[i] = dist[i];
		}
		n_model = n_grid_points - grid_point_centre - grid_points_in_sheet;

/* Get model 2nd derivative */

		spline(	model_distance, model_aureole, n_model, 0.0, 0.0, model_2deriv	);
	
/* Determine misfit */

		for (int i=1; i<=n_aureole; i++ ) {
			if ( data_aureole[0][i] < model_distance[1] || data_aureole[0][i] > model_distance[n_model] ) continue;
			if ( distance_rel_to_sheet_edge ) {
				if ( data_aureole[0][i] > -window_dless_width*sheet_thickness ) continue;
			} else {
				if ( data_aureole[0][i] > -(window_dless_width+0.5)*sheet_thickness ) continue;
			}
			n_aureole_used++;
			splint( model_distance, model_aureole, model_2deriv, n_model, data_aureole[0][i], &model_value	);
			double dif = (data_aureole[1][i] - model_value)/data_aureole[2][i];
			misfit += dif*dif;
			cerr	<< data_aureole[0][i] << "\t"
						<< data_aureole[1][i] << "\t"
						<< model_value << "\t"
						<< misfit << "\n";
		}
	}

/* RMS misfit for both sides of sheet */
	
	misfit /= n_aureole_used;
	if (misfit > 0.0) {
		misfit = pow(misfit, 0.5);
		cerr << "   RMS misfit is " << misfit << " from " << n_aureole_used << " data points for " << sheet_thickness << " m thick sheet\n";
		if ( print_misfit ) {
			cout	<< misfit << "\t"
						<< sheet_thickness << "\t"
						<< depth_of_sheet << "\t"
						<< heating_rate_from_burial << "\n";
			cerr << "   Writing misfit\n";
			cerr << "   column 1: misfit (normalized)\n";
			cerr << "   column 2: sheet thickness (m)\n";
			cerr << "   column 3: depth of sheet (m)\n";
			cerr << "   column 4: heating_rate_from_burial (C/yr)\n";
		}
	} else {
		cerr << "   Not enough useable data points to calculate misfit\n";
	}  
	
	
}



void find_solidification_time( void ) 
{

	double activation_energy;
	double max_flux;
	double reaction_rate;
	double sum;

	cerr << " Prior to magma injection\n";
	heating_rate_from_burial = burial_rate*geothermal_gradient;

/* STARTING TIME SUPPLIED */

	if ( constant_temperature_reaction ) {
		start_time = -begin_constant_temperature*time_scale;
		cerr << "   Constant temperature " << constant_temperature;
		cerr << " C beginning at " << start_time/time_scale << " yr\n";

/* STARTING TIME DEPENDS ON DEEPEST GRID POINT AND BURIAL RATE */

	} else {
		cerr << "   Burial rate " << burial_rate*time_scale*1.0e-3*1.0e6 
				<< " km/Myr, Geothermal gradient " << geothermal_gradient*1.0e3 << " C/km\n";
		cerr << "   Heating rate " << heating_rate_from_burial*time_scale*1.0e6 << " C/Myr\n";
	  start_time = sheet_is_dyke ? -depth_of_sheet : -dist[n_grid_points];
  	start_time /= burial_rate;
  }
	cerr << "   Starting time " << start_time/time_scale*1.0e-6 << " Myr relative to magma injection\n";

/* 
   FIND DIMENSIONLESS SOLIDIFICATION BOUNDARY
   USING EQ 4.149 OF TURCOTTE & SCHUBERT 2002 
*/

	cerr << " After magma injection\n";

/* 
	Find ambient temperature at time of injection from geothermal gradient.
	If the sheet is a sill, assume the ambient temperature is the average of
  the temperatures at the top and the base.	
 */

	ambient_temperature = surface_temperature + geothermal_gradient*depth_of_sheet;
	cerr << "   Ambient temperature " << ambient_temperature 
			<< " C, Magma temperature " << magma_temperature << " C\n";

/* GET THERMAL PARAMETERS */

  thermal_parameters ( magma_temperature );

/* Solves Eq 4.149 of Turcotte & Schubert */

	left_side = heat_ratio * sqrt(PI) / (magma_temperature-ambient_temperature);
	dless_solidification_boundary = rtbis( 1.0e-3, 10.0 , 0.0001, find_dless_solidification_boundary );
	erfc_dless_solidification_boundary = erfc(-dless_solidification_boundary);
	double wall_temperature = ambient_temperature + (magma_temperature-ambient_temperature) / (1.0+erf(dless_solidification_boundary));
	solidification_time = pow(sheet_half_thickness, 2.0) / 4.0 / thermal_diffusivity / dless_solidification_boundary/dless_solidification_boundary;
	double boil_zone = 2.0 * sqrt( solidification_time * thermal_diffusivity );
	solid_crust = -2.0 * dless_solidification_boundary * sqrt( thermal_diffusivity * end_time*time_scale);
	cerr << "   Time for sheet to solidify: " << solidification_time/time_scale << " " << time_unit <<  "\n";
	cerr << "   Wall temperature during solidification: " << wall_temperature << " C"<<  "\n";
	cerr << "   Denominator of d'less temperature equation (T&S 4-143): " << erfc_dless_solidification_boundary << "\n";

/* Maximum reaction rates 
	(when concentration of reactant = 1 at wall temperature during solidification) */

	for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ ) {
		if ( i_reaction == VITRINITE ) continue;
		reaction=i_reaction;
		sum = 0.0;
//		for ( int i_approx=1; i_approx<=n_approx_reactions[i_reaction]; i_approx++ ) {
			activation_energy = approx_reaction_energies[i_reaction][n_approx_reactions[i_reaction]];
			reaction_rate = arrhenius_reaction_rate_coeff( activation_energy, wall_temperature );
			sum += reaction_rate;
/*
			aureole_model[OIL][i_grid] *= (1.0-mass_frac_labile_to_gas) * mass_fraction_labile * carbon_to_methane * density_country_rock;
			aureole_model[LABILE][i_grid] *= mass_frac_labile_to_gas * mass_fraction_labile * carbon_to_methane * density_country_rock;
			aureole_model[REFRACTORY][i_grid] *= mass_fraction_refractory * carbon_to_methane * density_country_rock;
*/
//		}
		if ( i_reaction == OIL  ) {
			max_flux = sum * (1.0-mass_frac_labile_to_gas) * mass_fraction_labile * carbon_to_methane * density_country_rock  * PI * sill_radius_squared;
		};
		if ( i_reaction == LABILE  ) {
			max_flux = sum * mass_fraction_labile * carbon_to_methane * density_country_rock  * PI * sill_radius_squared;
		};
		if ( i_reaction == REFRACTORY  ) {
			max_flux = sum * mass_fraction_refractory * carbon_to_methane * density_country_rock  * PI * sill_radius_squared;
		};
		cerr << "   Reaction " << i_reaction << ";  max rate: " << reaction_rate*time_scale << " per " << time_unit <<";  max flux: " << max_flux*time_scale << " per " << time_unit << "\n";
	}
/*				
				temperature_0 = temperature[i_grid];
				temperature_1 = temperature[i_grid]+temperature_increment[i_grid];
				reaction_rate_0 = arrhenius_reaction_rate_coeff( activation_energy, temperature_0 );
				reaction_rate_1 = arrhenius_reaction_rate_coeff( activation_energy, temperature_1 );
				reaction_rate_harm = 1.0 / ( 1.0/reaction_rate_0 + 1.0/reaction_rate_1);
				reaction_rate_mean = 0.5*(reaction_rate_0 + reaction_rate_1);
*/

/* WRITE OUT IF REQUIRED */

	if ( print_solidification_time ) {
	  cerr << " Writing out sheet width (m), solidification time (yr), wall temp. (C)\n\n";
		cout	<< sheet_thickness << "\t" 
				<< solidification_time/time_scale << "\t"
//				<< solid_crust << "\t"
				<< wall_temperature << "\n";
	}
}


double find_time_step (	double time	)
{
	double t;

/* FIND STABLE TIME STEPS ON FIRST TIME STEP AND STORE THEM */

  if ( stable_burial_time_step == 0.0 && stable_conduction_time_step == 0.0 ) {
  
/* FIND MINIMUM DISTANCE INCREMENT */

		double dz_min = 1.0e30;
		for ( int i=2; i<=n_grid_points; i++)
			if ( dz_min > abs(dist[i]-dist[i-1]) ) dz_min = abs(dist[i]-dist[i-1]);
		double dz_min_sq = pow(dz_min,2.0);

/* MAXIMUM TIME STEP FOR WHEN HEAT MOVES BY BURIAL */

    stable_burial_time_step = dz_min / burial_rate;

/* MAXIMUM TIME STEP FOR WHEN HEAT MOVES BY CONDUCTION */

		thermal_parameters(surface_temperature);
		stable_conduction_time_step = dz_min_sq / thermal_diffusivity / 2.0;
		cerr << "   Stable conduction time step: " << stable_conduction_time_step/time_scale << " yr\n";

/* SET UP ARRAY FOR STORING TIME SERIES */

		t = FMIN( stable_conduction_time_step, time_step_after );
		N_MAX_TIME_SERIES_ARRAY = end_time / t + 10;
		time_series_flux = dvector(1,N_MAX_TIME_SERIES_ARRAY);
		time_series_flux_upper = dvector(1,N_MAX_TIME_SERIES_ARRAY);
		time_series_mass = dvector(1,N_MAX_TIME_SERIES_ARRAY);
		time_series_mass_upper = dvector(1,N_MAX_TIME_SERIES_ARRAY);

	}

/* BEFORE TIME OF MAGMA INJECTION */

	if ( time < 0.0 ) {
		if ( constant_temperature_reaction ) {
//			t = ; 
		} else {
			t = FMIN( stable_burial_time_step, time_step_before );
//		t = time_step_before;
  	  if ( time+t > 0.0 ) t = 0.0 - time;
		}

/* AFTER TIME OF MAGMA INJECTION */

	} else {
		t = FMIN( stable_conduction_time_step, time_step_after );
	}

	return t;
}




void take_temperature_step(double time, double time_step)
{
	int grid_point_wall;
	double depth;
	double distance;
	double dist_increment_sq;
	double dless_temperature;

/*
	Before injection of magma sheet, the temperature is controlled by the
	burial rate and the geothermal gradient.  
	Also applies to time of injection because analytical 
	solution not valid for time=0.
 */

	if ( time <= 0.0 && !fp_starting_temperature ) {
		for (int i=1; i<=n_grid_points; i++) {
			if ( constant_temperature_reaction ) {
				temperature[i] = constant_temperature; 
				temperature_increment[i] = 0.0;
			} else {
				depth = (sheet_is_dyke) ? depth_of_sheet : dist[i];
				double temp = surface_temperature + geothermal_gradient*( depth + burial_rate*time );
				if ( temp >= surface_temperature ) {
					temperature[i] = temp; 
					temperature_increment[i] = geothermal_gradient * burial_rate * time_step;
				}
			}
//			cerr << dist[i] << "\t" << temperature[i] << "\t" << temperature_increment[i] << "\n";
		}

/* 
  ANALYTICAL SOLUTION IF SILL HAS NOT FULLY SOLIDIFIED 
  TURCOTTE & SCHUBERT EQS 4.142 - 4.145
*/ 

	} else if ( time > 0.0 && time <= solidification_time && !fp_starting_temperature ) {

		thermal_parameters (magma_temperature);
		double solid_crust = -2.0 * dless_solidification_boundary * sqrt( thermal_diffusivity * time);

		double dist_scale = 2.0 * sqrt( thermal_diffusivity * time);
		for (int i=1; i<=n_grid_points; i++) {
			double temperature_start = temperature[i];

/* 
  Set local ambient temperature to temperature sheet margin would
  have if sheet present.
  This method is good for dykes and some sills.
  It will introduce uncertainty for sills intruded shallow enough
  to feel the land surface boundary.  These will be tackled using
  finite differencing with the modified thermal diffusivity method
  in the future.  
  A smaller uncertainty is introduced because the analytical solution
  assumes background temperature, whereas a sill is intruded within
  a temperature gradient.  This effect is not expected to be
  significant for gas expulsion because most gas generation occurs
  during and just after sill solidification, and the temperature 
  curvature is dominated by the sill and not the geotherm at these
  times.
*/

		if ( sheet_is_dyke ) {
			depth = depth_of_sheet;
		} else {
			int n_point_centre = ( n_grid_points + 1 ) / 2;
			if ( i <= n_point_centre ) {
				grid_point_wall = n_point_centre - grid_points_in_sheet;
			} else {
				grid_point_wall = n_point_centre + grid_points_in_sheet;
			}
			depth = dist[grid_point_wall];
			depth = dist[i];
		}
		double local_ambient_temperature = surface_temperature + geothermal_gradient*( depth + burial_rate*time );
		double temperature_difference = magma_temperature - local_ambient_temperature;

/* 
 FOR THE FOLLOWING ANALYTICAL SOLUTION, DISTANCE IS RELATIVE TO THE THE SHEET
 MARGIN.  IN THE MAIN PROGRAM, DISTANCE IS RELATIVE TO THE SHEET CENTRE AND
 CAN BE POSTIVE OR NEGATIVE.  
 */

			distance =  distance_wrt_sheet_margin ( dist[i] );

/* ANALYTICAL SOLUTION */

			if ( distance < solid_crust ) 
				temperature[i] = magma_temperature;
			else {
				double dless_dist = distance / dist_scale;
				dless_temperature = erfc(dless_dist) / erfc_dless_solidification_boundary;
				temperature[i] = local_ambient_temperature + dless_temperature * temperature_difference;
			}

/* Must also work out temperature increment */

			temperature_increment[i] = temperature[i] - temperature_start;

//cerr << dist[i] << "\t" << dless_temperature << "\t" << temperature[i] << "\n";
		}
/* 
  FINITE DIFFERENCE SOLUTION AFTER SILL HAS SOLIDIFIED 
*/ 
	} else {

/* If the temperature is symmetrical about this point */

		if ( sheet_is_dyke ) {
			thermal_parameters (temperature[2]);
			dist_increment_sq = pow( (dist[2]-dist[1]) , 2.0);
			temperature_increment[1] = thermal_diffusivity*2.0*(temperature[2]-temperature[1])/dist_increment_sq;

/* FTCS Explicit differencing for the interior nodes */

    	for (int i=2; i<n_grid_points; i++) {
				thermal_parameters (temperature[i]);
				dist_increment_sq = pow( (dist[i]-dist[i-1]) , 2.0);
				temperature_increment[i] = thermal_diffusivity*(temperature[i+1]-2.0*temperature[i]+temperature[i-1])/dist_increment_sq;
			}

/* Add temperature increment to find final temperature */

      for (int i=1; i<n_grid_points; i++)
        temperature[i] += time_step*temperature_increment[i];

/* Exterior point far from dyke */

			temperature[n_grid_points] = 2.0*temperature[n_grid_points-1] - temperature[n_grid_points-2];

/* FINITE DIFFERENCE SOLUTION FOR SILL */

		} else {

/* FTCS Explicit differencing for the interior nodes */

    	for (int i=2; i<n_grid_points; i++) {
				thermal_parameters (temperature[i]);
				dist_increment_sq = pow( (dist[i]-dist[i-1]) , 2.0);
				temperature_increment[i] = thermal_diffusivity*(temperature[i+1]-2.0*temperature[i]+temperature[i-1])/dist_increment_sq;
			}
			
/* Assume constant curvature for deepest node */

//				temperature_increment[n_grid_points] = temperature_increment[n_grid_points-1];
				
//				temperature_increment[n_grid_points] = ( temperature_increment[n_grid_points-3] - 2.0*temperature_increment[n_grid_points-2] + temperature_increment[n_grid_points-1] ) * ( dist[n_grid_points-2] - 2.0*dist[n_grid_points-1] + dist[n_grid_points] ) / ( dist[n_grid_points-3] - 2.0*dist[n_grid_points-2] + dist[n_grid_points-1] ) + 2.0*temperature_increment[n_grid_points-1] - temperature_increment[n_grid_points-2];

				temperature_increment[n_grid_points] = ( temperature_increment[n_grid_points-3] - 2.0*temperature_increment[n_grid_points-2] + temperature_increment[n_grid_points-1] ) + 2.0*temperature_increment[n_grid_points-1] - temperature_increment[n_grid_points-2];

				
//			thermal_parameters (temperature[n_grid_points]);
//			dist_increment_sq = pow( (dist[n_grid_points]-dist[n_grid_points-1]) , 2.0);
//			temperature_increment[n_grid_points] = thermal_diffusivity*2.0*(temperature[n_grid_points]-temperature[n_grid_points-1])/dist_increment_sq;


/* Add temperature increment to find final temperature */

      for (int i=2; i<=n_grid_points; i++)
        temperature[i] += time_step*temperature_increment[i];

		}
  }

	return;
}


double distance_wrt_sheet_margin ( double distance )
{
return (sheet_is_dyke) ? abs(distance)-sheet_half_thickness : abs(distance-depth_of_sheet)-sheet_half_thickness;
}



void thermal_parameters ( double temperature )
{
 
/*  THERMAL DIFFUSIVITY ALTERED TO ACCOUNT FOR LATENT HEAT OF CRYSTALLIZATION
 */

	const double liquidus_temperature = 1350.0; 
	const double magma_density = 2720.0;
	double thermal_conductivity;

/* THESE VALUES USED IN THE ORIGINAL VERSIONS OF THE CODE V1 AND V2
   FOR 2 NERC PROPOSALS
		latent_heat = 320.0e3;
		specific_heat = 1200.0;
		thermal_conductivity = 1.632;
 */

/* THESE VALUES FROM MACLENNAN & LOVELL, GEOLOGY, 2002
 */

	if ( latent_heat_scheme==0 ) {
		latent_heat = 700.0e3;
		specific_heat = 1700.0;
		thermal_conductivity = 1.632;


	} else if ( latent_heat_scheme==1 ) {
		latent_heat = 700.0e3;
		specific_heat = 1700.0;
		if ( temperature > magma_temperature ) specific_heat += latent_heat*( liquidus_temperature - temperature );
		thermal_conductivity = 1.18 + (474.0/(350.0 + temperature));

/* INCORRECT -L VALUE */

	} else {
		cerr << " ERROR: Unknown value in -L\n";
		exit(0);
	}
//	thermal_diffusivity = thermal_conductivity / magma_density / specific_heat;
	thermal_diffusivity = 1.0e-6;
	heat_ratio = latent_heat / specific_heat;
}




void take_reaction_step ( double time, double time_step )
{
/*
 	The breakdown of labile and refractory kerogen and the maturation of vitrinite
  is described by irreversible reactions.
  
  Each reaction is governed by 
  
  	dR/dt = -k.R	[1]
  
  where R(t) is the concentration of the reactant.
  The concentration of the product (= progress of the reaction) is P = 1-R, so
  
  	d(1-P)/dt = -k.(1-P)	[2]

  In each time step we assume the reaction rate k is constant, so the
  solution to [2] is
  
  	P = 1 - (1 - P0).exp(-k.t)		[3]
  	
  where P0 is the concentration of the product at the start of the time step.  
  The other reaction considered is the oil cracking reaction, governed by

  	d(1-P)/dt = -k.(1-P) - S	[4]

  where S(t) is the rate of oil generation from breakdown of labile kerogen.  
  If we assume that S is constant over each time step then the solution to [4] is
  
  	P = 1 - S/k - (1 - P0 - S/k).exp(-k.t)		[5]
  
 */
 
 // (k - k P - S + E^(k T) (S + k (-1 + P + k T - S T)))/(E^(k T) k^2)

	int n_max;
	double activation_energy;
	double final_product_conc;
	double initial_product_conc;
	double oil_production_rate;
	double progress;
	double reaction_rate;
	double reaction_rate_mean;
	double reaction_rate_harm;
	double reaction_rate_0;
	double reaction_rate_1;
	double S_over_k;
	double temperature_average;
	double temperature_0;
	double temperature_1;

/* 
  If the sheet is a dyke, only need organic reactions
  at one node because the temperature structure is 1D
  prior to injection of dyke 
*/

	if ( sheet_is_dyke ) {
		n_max = (time < 0.0) ? 1 : n_grid_points;
	} else {
		n_max = n_grid_points;
	}

/* Loop over grid nodes to carry out reaction step */

	for ( int i_grid=1; i_grid<=n_max; i_grid++ ) {

/* No organic reactions within the magma sheet */

		if ( abs(dist[i_grid]) < sheet_half_thickness ) continue;

// No need to do reactions below region where products will be integrated

//	for ( int i_grid=1; i_grid<=i_integrate_down_to; i_grid++ ) {


/* 
	Temperature array holds temperature at the end of the time step
	during which the temperature changed by temperature_increment.
	The mean of the initial and final temperatures will be used to
	estimate the mean reaction rate during the time step.  
*/

		temperature_average = temperature[i_grid] - 0.5*temperature_increment[i_grid];

/* Loop over the 4 kerogen reactions */

		for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ ) {
		
/* No calculation if reaction is not required */

			if ( !do_labile_reaction && i_reaction==LABILE ) continue;
			if ( !do_labile_reaction && i_reaction==OIL ) continue;
			if ( !do_refractory_reaction && i_reaction==REFRACTORY ) continue;
			if ( !do_vitrinite_reaction && !print_vitrinite && !print_time_series && i_reaction==VITRINITE ) continue;
			if ( print_vitrinite && i_reaction!=VITRINITE ) continue;
			
/* Communication of reaction type to activation_energy and reaction_rate functions */

			reaction = i_reaction;
			
/* Loop over the individual reactions approximating each kerogen distribution */

			for ( int i_approx=1; i_approx<=n_approx_reactions[i_reaction]; i_approx++ ) {
				initial_product_conc = progress_of_reactions[i_approx][i_reaction][i_grid];

/* Find the reaction rate at the mean temperature */
				
				activation_energy = approx_reaction_energies[i_reaction][i_approx];
				reaction_rate = arrhenius_reaction_rate_coeff( activation_energy, temperature_average );
				
				temperature_0 = temperature[i_grid];
				temperature_1 = temperature[i_grid]+temperature_increment[i_grid];
				reaction_rate_0 = arrhenius_reaction_rate_coeff( activation_energy, temperature_0 );
				reaction_rate_1 = arrhenius_reaction_rate_coeff( activation_energy, temperature_1 );
				reaction_rate_harm = 1.0 / ( 1.0/reaction_rate_0 + 1.0/reaction_rate_1);
				reaction_rate_mean = 0.5*(reaction_rate_0 + reaction_rate_1);
//				if ( i_grid == n_print ) cerr << "t, k = " << time/time_scale << "\t" << reaction_rate << "\t" << reaction_rate_mean << "\t" << reaction_rate_harm << "\n";
				reaction_rate = reaction_rate_harm;

/* Take reaction step using equation [3] above */

				if ( i_reaction != OIL ) {
					progress_of_reactions[i_approx][i_reaction][i_grid] = (1.0 - ( 1.0 - initial_product_conc ) * exp( -reaction_rate * time_step ));

/* Save reaction rate */

					rate_of_reactions[i_approx][i_reaction][i_grid] = ( 1.0 - initial_product_conc ) * ( 1.0 - exp( -reaction_rate * time_step )) / time_step;

/*	Oil cracking reaction */

				} else {

/* The oil production rate is spread evenly over the oil cracking reactions, hence divide by n_approx_reactions[OIL] */

					S_over_k = (reaction_rate==0) ? 0.0 : oil_production_rate / reaction_rate / n_approx_reactions[OIL]; 

/* Take reaction step using equation [5] above */

					progress_of_reactions[i_approx][i_reaction][i_grid] = 1.0 - S_over_k - ( 1.0 - initial_product_conc - S_over_k ) * exp( -reaction_rate * time_step );

/* Save average reaction rate in time step */

//					rate_of_reactions[i_approx][i_reaction][i_grid] = (1.0 - S_over_k) + ( -1.0 - exp( -reaction_rate * time_step )) * ( ( initial_product_conc - 1.0 )/reaction_rate + S_over_k ) / time_step;
					rate_of_reactions[i_approx][i_reaction][i_grid] = time_step * ( 1.0 - S_over_k) - exp( -reaction_rate * time_step ) * ( initial_product_conc/reaction_rate - 1.0/reaction_rate + S_over_k );
					rate_of_reactions[i_approx][i_reaction][i_grid] = -1.0 * ( initial_product_conc/reaction_rate - 1.0/reaction_rate + S_over_k );
					rate_of_reactions[i_approx][i_reaction][i_grid] = -1.0 * ( initial_product_conc*reaction_rate - reaction_rate + oil_production_rate/n_approx_reactions[OIL] ) / reaction_rate / reaction_rate;
					rate_of_reactions[i_approx][i_reaction][i_grid] = ( 1.0 - initial_product_conc - S_over_k ) * ( 1.0 - exp( -reaction_rate * time_step )) / time_step;
				}

/* Store average oil production rate over the time step for use in oil cracking calculation */

				if ( i_reaction==LABILE ) {
					if (i_approx == 1) oil_production_rate = 0.0;
					
//					oil_production_rate += -reaction_rate * 0.5 * ( (1.0-initial_product_conc) + (1.0-progress_of_reactions[i_approx][i_reaction][i_grid]) ) / n_approx_reactions[LABILE];

					oil_production_rate += -reaction_rate * (1.0-progress_of_reactions[i_approx][i_reaction][i_grid]) / n_approx_reactions[LABILE];

/*
					oil_production_rate += -reaction_rate * (
						(1.0-progress_of_reactions[i_approx][i_reaction][i_grid]) * time -
						(1.0-initial_product_conc) * (time-time_step) ) 
						/ time_step / n_approx_reactions[LABILE];
*/

					if (i_approx == n_approx_reactions[LABILE]) oil_production_rate *= ( 1.0 - mass_frac_labile_to_gas );
				}
//				if ( i_grid == n_print ) cerr << "S = " << oil_production_rate << "\n";

			} // End of approximate reaction loop
						
		} // End of kerogen type loop

	} // End of grid point loop

}



void scale_reaction_products_old ( void )
{
	double progress;

	for (int i_grid=1; i_grid<=n_grid_points; i_grid++) {

/* No organic reactions within the sheet */
/*
		for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ )
			if ( abs(dist[i_grid]-depth_of_sheet) < sheet_half_thickness ) {
				progress_of_reactions[1][i_reaction][i_grid] = -999;
				continue;
			}
*/
/* Sum the approximate reactions */

		for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ ) {
			progress = get_progress( i_reaction, i_grid );
			progress_of_reactions[1][i_reaction][i_grid] = progress;
		}

/* SCALE PROFILES ACCORDING TO REQUIRED OUTPUT */

		if ( print_vitrinite ) 
			progress_of_reactions[1][VITRINITE][i_grid] = reflectance_min + progress_of_reactions[1][VITRINITE][i_grid] * (reflectance_max-reflectance_min);

		if ( print_flux || print_mass_generated ) {
//			progress_of_reactions[1][OIL][i_grid] = (progress_of_reactions[1][LABILE][i_grid]-progress_of_reactions[1][OIL][i_grid]) * (1.0-mass_frac_labile_to_gas) * mass_fraction_labile * carbon_to_methane * density_country_rock;
			progress_of_reactions[1][OIL][i_grid] *= (1.0-mass_frac_labile_to_gas) * mass_fraction_labile * carbon_to_methane * density_country_rock;
			progress_of_reactions[1][LABILE][i_grid] *= mass_frac_labile_to_gas * mass_fraction_labile * carbon_to_methane * density_country_rock;
			progress_of_reactions[1][REFRACTORY][i_grid] *= mass_fraction_refractory * carbon_to_methane * density_country_rock;
		}

	} // End of grid point loop
}


void sum_reactions_find_max_rate ( double time, double time_step, int i_grid_base )
{

/* 
	Sums individual reactions to find progress of breakdown of each kerogen type.
	Scales progress to find amount of oil/gas generated.
	Also searches for grid node with highest reaction rate.
*/

	double max_rate;
	double previous_progress;
	double progress;
	double rate;

/* Initialize counters for finding maximum reaction rate */

	for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ ) i_max_rate[i_reaction] = n_grid_points;

/* Loop over reaction types */

	for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ ) {
		max_rate = 0.0;

/* Loop over grid */

		for (int i_grid=1; i_grid<=n_grid_points; i_grid++) {

/* Sum the approximate reactions.  Save aureole shape this time step */

			previous_progress = aureole_model[i_reaction][i_grid];
			aureole_model_prev[i_reaction][i_grid] = previous_progress;
			progress = get_progress( i_reaction, i_grid );
			aureole_model[i_reaction][i_grid] = progress;
			
/* Save progress at time of injection */

			if ( time == 0.0 ) {
				progress_at_injection[i_reaction][i_grid] = previous_progress;
			}

/* Test to find maximum reaction rate */

//			if ( i_grid > i_grid_base ) {
			if ( i_grid > i_grid_base ) {
				rate = (progress - previous_progress) / (time_step/time_scale);
//				rate = (progress - previous_progress) ;
				if ( rate > max_rate ) {
					i_max_rate[i_reaction] = i_grid;
					max_rate = rate;
				}
//				if (time/time_scale > 1.258) cerr << i_grid << "\t" << previous_progress << "\t" << progress << "\t" << rate << "\n";
			}

		} // End of grid loop
	} // End of reaction loop
	if ( time == 0.0 ) cerr << "   Saved progress of reactions at injection time, t=0\n";
}


void scale_reaction_products ( void )
{

/* 
	Scales progress to find amount of oil/gas generated.
*/

/* Loop over grid */

	for (int i_grid=1; i_grid<=n_grid_points; i_grid++) {

/* SCALE PROFILES ACCORDING TO REQUIRED OUTPUT */

		if ( print_vitrinite ) 
			aureole_model[VITRINITE][i_grid] = reflectance_min + progress_of_reactions[1][VITRINITE][i_grid] * (reflectance_max-reflectance_min);

		if ( print_flux || print_mass_generated ) {
			aureole_model[OIL][i_grid] *= (1.0-mass_frac_labile_to_gas) * mass_fraction_labile * carbon_to_methane * density_country_rock;
			aureole_rates[OIL][i_grid] *= (1.0-mass_frac_labile_to_gas) * mass_fraction_labile * carbon_to_methane * density_country_rock;
			aureole_model[LABILE][i_grid] *= mass_frac_labile_to_gas * mass_fraction_labile * carbon_to_methane * density_country_rock;
			aureole_rates[LABILE][i_grid] *= mass_frac_labile_to_gas * mass_fraction_labile * carbon_to_methane * density_country_rock;
			aureole_model[REFRACTORY][i_grid] *= mass_fraction_refractory * carbon_to_methane * density_country_rock;
			aureole_rates[REFRACTORY][i_grid] *= mass_fraction_refractory * carbon_to_methane * density_country_rock;
		}

	} // End of grid point loop
}



void write_out_time_series ( double time, bool print_labels )
{

	double depth = (sheet_is_dyke) ? depth_of_sheet : dist[n_print];
	double distance = (sheet_is_dyke) ? dist[n_print] : dist[n_print]-depth_of_sheet;
	double dis;
	double progress;
		
	if ( print_labels ) {
		if ( print_time_series_at_point ) 
			cerr << " Writing out time series at node " << n_print << "\n";
		else
			cerr << " Writing out time series\n";
			cerr << "   column 1, time after magma injection (yr)\n";

		if ( print_temperature ) {
			cerr << "   column 2, temperature (C)\n";
			cerr << "   column 3, distance from sheet centre (m)\n";
			cerr << "   column 4, depth from surface (m)\n";

		} else if ( print_reactions ) {
			cerr << "   column 2, labile kerogen reaction\n";
			cerr << "   column 3, refractory kerogen reaction\n";
			cerr << "   column 4, vitrinite reaction\n";
			cerr << "   column 5, oil cracking reaction\n";
			cerr << "   column 6, distance from sheet centre (m)\n";
			cerr << "   column 7, depth from surface (m)\n";

		} else if ( print_flux ) {
			cerr << "   column 2, methane flux density (Mt/kg/yr)\n";
			cerr << "   column 3, distance from sheet centre (m)\n";
			cerr << "   column 4, depth from surface (m)\n";

		} else if ( print_vitrinite ) {
			cerr << "   column 2, vitrinite reflectance (%Ro)\n";
			cerr << "   column 3, distance from sheet centre (m)\n";
			cerr << "   column 4, depth from surface (m)\n";

		} else if ( print_mass_generated ) {
			if ( scale_to_whole_sill ) {
				cerr << "   column 2, methane mass from labile kerogen reaction (Mt)\n";
				cerr << "   column 3, methane mass from refractory kerogen reaction (Mt)\n";
				cerr << "   column 4, methane mass from oil cracking reaction (Mt)\n";
				cerr << "   column 5, total methane mass (Mt)\n";
				if ( subtract_mass_at_injection )
					cerr << "      NB methane already generated at time of injection subtracted from all masses\n"; 
				cerr << "   column 6, methane flux from labile kerogen reaction (Mt/yr)\n";
				cerr << "   column 7, methane flux from refractory kerogen reaction (Mt/yr)\n";
				cerr << "   column 8, methane flux from oil cracking reaction (Mt/yr)\n";
				cerr << "   column 9, total methane flux (Mt/yr)\n";
				cerr << "      NB all area densities have been multiplied by pi\n";
			} else {
				cerr << "   column 2, methane area density from labile kerogen reaction (Mt/m^2)\n";
				cerr << "   column 3, methane area density from refractory kerogen reaction (Mt/m^2)\n";
				cerr << "   column 4, methane area density from oil cracking reaction (Mt/m^2)\n";
				cerr << "   column 5, total methane area density (Mt/m^2)\n";
				if ( subtract_mass_at_injection )
					cerr << "      NB methane already generated at time of injection subtracted from all area densities\n"; 
				cerr << "   column 6, methane area density flux from labile kerogen reaction (Mt/m^2/yr)\n";
				cerr << "   column 7, methane area density flux from refractory kerogen reaction (Mt/m^2/yr)\n";
				cerr << "   column 8, methane area density flux from oil cracking reaction (Mt/m^2/yr)\n";
				cerr << "   column 9, total methane area density flux (Mt/m^2/yr)\n";
				cerr << "      NB all area densities have been multiplied by pi\n"; 
			}
			if ( sheet_is_dyke) {
				cerr << "   column 10, position of outer edge of aureole, labile (km)\n";
				cerr << "   column 11, position of inner edge of aureole, labile (km)\n";
				cerr << "   column 12, position of max reaction rate, labile (km)\n";
				cerr << "   column 13, temperature at max reaction rate, labile (C)\n";
				cerr << "   column 14, position of outer edge of aureole, refractory (km)\n";
				cerr << "   column 15, position of inner edge of aureole, refractory (km)\n";
				cerr << "   column 16, position of max reaction rate, refractory (km)\n";
				cerr << "   column 17, temperature at max reaction rate, refractory (C)\n";
				cerr << "   column 18, position of outer edge of aureole, cracking (km)\n";
				cerr << "   column 19, position of inner edge of aureole, cracking (km)\n";
				cerr << "   column 20, position of max reaction rate, cracking (km)\n";
				cerr << "   column 21, temperature at max reaction rate, cracking (C)\n";
				cerr << "   column 22, position of outer edge of aureole, vitrinite (km)\n";
				cerr << "   column 23, position of inner edge of aureole, vitrinite (km)\n";
				cerr << "   column 24, position of max reaction rate, vitrinite (km)\n";
				cerr << "   column 25, temperature at max reaction rate, vitrinite (C)\n";
				cerr << "   column 26, %Ro at outer edge of aureole, labile\n";
				cerr << "   column 27, %Ro at inner edge of aureole, labile\n";
				cerr << "   column 28, %Ro at max reaction rate in aureole, labile\n";				
				cerr << "   column 29, %Ro at outer edge of aureole, refractory\n";
				cerr << "   column 30, %Ro at inner edge of aureole, refractory\n";
				cerr << "   column 31, %Ro at max reaction rate in aureole, refractory\n";
				cerr << "   column 32, %Ro at outer edge of aureole, cracking\n";
				cerr << "   column 33, %Ro at inner edge of aureole, cracking\n";
				cerr << "   column 34, %Ro at max reaction rate in aureole, cracking\n";
				cerr << "   column 35, total dless aureole width (upper & lower), labile\n";
				cerr << "   column 36, total dless aureole width (upper & lower), refractory\n";
				cerr << "   column 37, total dless aureole width (upper & lower), cracking\n";
				cerr << "   column 38, total dless aureole width (upper & lower), vitrinite\n";

			} else {
				cerr << "   column 10, position of outer edge of upper aureole, labile (km)\n";
				cerr << "   column 11, position of inner edge of upper aureole, labile (km)\n";
				cerr << "   column 12, position of max reaction rate in upper aureole, labile (km)\n";
				cerr << "   column 13, temperature at max reaction rate in upper aureole, labile (C)\n";
				cerr << "   column 14, position of outer edge of lower aureole, labile (km)\n";
				cerr << "   column 15, position of inner edge of lower aureole, labile (km)\n";
				cerr << "   column 16, position of max reaction rate in lower aureole, labile (km)\n";
				cerr << "   column 17, temperature at max reaction rate in lower aureole, labile (C)\n";

				cerr << "   column 18, position of outer edge of upper aureole, refractory (km)\n";
				cerr << "   column 19, position of inner edge of upper aureole, refractory (km)\n";
				cerr << "   column 20, position of max reaction rate in upper aureole, refractory (km)\n";
				cerr << "   column 21, temperature at max reaction rate in upper aureole, refractory (C)\n";
				cerr << "   column 22, position of outer edge of lower aureole, refractory (km)\n";
				cerr << "   column 23, position of inner edge of lower aureole, refractory (km)\n";
				cerr << "   column 24, position of max reaction rate in lower aureole, refractory (km)\n";
				cerr << "   column 25, temperature at max reaction rate in lower aureole, refractory (C)\n";

				cerr << "   column 26, position of outer edge of upper aureole, cracking (km)\n";
				cerr << "   column 27, position of inner edge of upper aureole, cracking (km)\n";
				cerr << "   column 28, position of max reaction rate in upper aureole, cracking (km)\n";
				cerr << "   column 29, temperature at max reaction rate in upper aureole, cracking (C)\n";
				cerr << "   column 30, position of outer edge of lower aureole, cracking (km)\n";
				cerr << "   column 31, position of inner edge of lower aureole, cracking (km)\n";
				cerr << "   column 32, position of max reaction rate in lower aureole, cracking (km)\n";
				cerr << "   column 33, temperature at max reaction rate in lower aureole, cracking (C)\n";

				cerr << "   column 34, position of outer edge of upper aureole, vitrinite (km)\n";
				cerr << "   column 35, position of inner edge of upper aureole, vitrinite (km)\n";
				cerr << "   column 36, position of max reaction rate in upper aureole, vitrinite (km)\n";
				cerr << "   column 37, temperature at max reaction rate in upper aureole, vitrinite (C)\n";
				cerr << "   column 38, position of outer edge of lower aureole, vitrinite (km)\n";
				cerr << "   column 39, position of inner edge of lower aureole, vitrinite (km)\n";
				cerr << "   column 40, position of max reaction rate in lower aureole, vitrinite (km)\n";
				cerr << "   column 41, temperature at max reaction rate in lower aureole, vitrinite (C)\n";			

				cerr << "   column 42, %Ro at outer edge of upper aureole, labile\n";
				cerr << "   column 43, %Ro at inner edge of upper aureole, labile\n";
				cerr << "   column 44, %Ro at max reaction rate in upper aureole, labile\n";
				cerr << "   column 45, %Ro at outer edge of lower aureole, labile\n";
				cerr << "   column 46, %Ro at inner edge of lower aureole, labile\n";
				cerr << "   column 47, %Ro at max reaction rate in lower aureole, labile\n";
				
				cerr << "   column 48, %Ro at outer edge of upper aureole, refractory\n";
				cerr << "   column 49, %Ro at inner edge of upper aureole, refractory\n";
				cerr << "   column 50, %Ro at max reaction rate in upper aureole, refractory\n";
				cerr << "   column 51, %Ro at outer edge of lower aureole, refractory\n";
				cerr << "   column 52, %Ro at inner edge of lower aureole, refractory\n";
				cerr << "   column 53, %Ro at max reaction rate in lower aureole, refractory\n";

				cerr << "   column 54, %Ro at outer edge of upper aureole, cracking\n";
				cerr << "   column 55, %Ro at inner edge of upper aureole, cracking\n";
				cerr << "   column 56, %Ro at max reaction rate in upper aureole, cracking\n";
				cerr << "   column 57, %Ro at outer edge of lower aureole, cracking\n";
				cerr << "   column 58, %Ro at inner edge of lower aureole, cracking\n";
				cerr << "   column 59, %Ro at max reaction rate in lower aureole, cracking\n";
				cerr << "   column 60, total dless aureole width (upper & lower), labile\n";
				cerr << "   column 61, total dless aureole width (upper & lower), refractory\n";
				cerr << "   column 62, total dless aureole width (upper & lower), cracking\n";
				cerr << "   column 63, total dless aureole width (upper & lower), vitrinite\n";
				cerr << "   column 64, methane area density from labile kerogen reaction, dless aureole method\n";
				cerr << "   column 65, methane area density from refractory kerogen reaction, dless aureole method\n";
				cerr << "   column 66, methane area density from oil cracking reaction, dless aureole method\n";
			}

		} else if ( print_mass_expelled ) {
			cerr << "   column 2, methane mass vented (Mt)\n";
			cerr << "   column 3, methane mass generated (Mt)\n";
			cerr << "   column 4, maximum methane mass vented (Mt)\n";
			cerr << "   column 5, methane venting flux (Mt/yr)\n";
			cerr << "   column 6, methane generation flux (Mt/yr)\n";
			cerr << "   column 7, maximum venting flux (Mt/yr)\n";

		}

	} else {
		if ( print_temperature ) {
			cout	<< time/time_scale << "\t"
					<< temperature[n_print] << "\t"
					<< distance << "\t" 
					<< depth << "\n"; 

		} else if ( print_reactions ) {
			cout << fixed << showpoint;
			cout << setprecision(5) << time/time_scale << "\t";
			cout.unsetf(ios_base::floatfield);
			if (progress_of_reactions[1][LABILE][n_print] < 0.0)	
				cout << "nan\t";
			else {
				progress = get_progress( LABILE, n_print );
				cout << progress << "\t";
			}
			if (progress_of_reactions[1][REFRACTORY][n_print] < 0.0)	
				cout << "nan\t";
			else {
				progress = get_progress( REFRACTORY, n_print );
				cout << progress << "\t";
			}
			if (progress_of_reactions[1][VITRINITE][n_print] < 0.0)	
				cout << "nan\t";
			else {
				progress = get_progress( VITRINITE, n_print );
				cout << progress << "\t";
			}
			if (progress_of_reactions[1][OIL][n_print] < 0.0)	
				cout << "nan\t";
			else {
				progress = get_progress( OIL, n_print );
				cout << progress << "\t";
			}
			cout	<< distance << "\t" 
					<< depth << "\n"; 

		} else if ( print_flux ) {
			cout << fixed << showpoint;
			cout << setprecision(5) << time/time_scale << "\t";
			cout.unsetf(ios_base::floatfield);
			cout	<< progress_of_reactions[1][LABILE][n_print] << "\t"
					<< progress_of_reactions[1][REFRACTORY][n_print] << "\t"
					<< progress_of_reactions[1][OIL][n_print] << "\t";
			cout	<< distance << "\t" 
					<< depth << "\n";

		} else if ( print_vitrinite ) {
			cout << fixed << showpoint;
			cout << setprecision(5) << time/time_scale << "\t";
			cout.unsetf(ios_base::floatfield);
			if (progress_of_reactions[1][VITRINITE][n_print] < 0.0)	
				cout << "nan\t";
			else {
				progress = get_progress( VITRINITE, n_print );
				progress = reflectance_min + progress * (reflectance_max-reflectance_min);
				cout << progress << "\t";
			}
			cout	<< distance << "\t" 
					<< depth << "\n"; 

		} else if ( print_mass_generated ) {
			cout << fixed << showpoint;
			cout << setprecision(5) << time/time_scale << "\t";
			cout.unsetf(ios_base::floatfield);
			if ( subtract_mass_at_injection ) {
				cout	<< mass_from_labile - mass_from_labile_at_injection << "\t"
							<< mass_from_refractory - mass_from_refractory_at_injection << "\t"
							<< mass_from_oil - mass_from_oil_at_injection << "\t"
							<< mass_methane_generated - mass_methane_generated_at_injection << "\t";
			} else {
				cout	<< mass_from_labile << "\t"
							<< mass_from_refractory << "\t"
							<< mass_from_oil << "\t"
							<< mass_methane_generated << "\t";
			}
			cout	<< flux_from_labile << "\t"
						<< flux_from_refractory << "\t"
						<< flux_from_oil << "\t"
						<< flux_methane_generated << "\t";

			if ( !sheet_is_dyke ) {
				cout << pos_outer_edge_upper_aureole[LABILE] << "\t";
				cout << pos_inner_edge_upper_aureole[LABILE] << "\t";
				cout << pos_max_rate_upper_aureole[LABILE] << "\t";
				cout << temp_max_rate_upper_aureole[LABILE] << "\t";
			}
			cout << pos_outer_edge_lower_aureole[LABILE] << "\t";
			cout << pos_inner_edge_lower_aureole[LABILE] << "\t";
			cout << pos_max_rate_lower_aureole[LABILE] << "\t";
			cout << temp_max_rate_lower_aureole[LABILE] << "\t";
			if ( !sheet_is_dyke ) {
				cout << pos_outer_edge_upper_aureole[REFRACTORY] << "\t";
				cout << pos_inner_edge_upper_aureole[REFRACTORY] << "\t";
				cout << pos_max_rate_upper_aureole[REFRACTORY] << "\t";
				cout << temp_max_rate_upper_aureole[REFRACTORY] << "\t";
			}
			cout << pos_outer_edge_lower_aureole[REFRACTORY] << "\t";
			cout << pos_inner_edge_lower_aureole[REFRACTORY] << "\t";
			cout << pos_max_rate_lower_aureole[REFRACTORY] << "\t";
			cout << temp_max_rate_lower_aureole[REFRACTORY] << "\t";
			if ( !sheet_is_dyke ) {
				cout << pos_outer_edge_upper_aureole[OIL] << "\t";
				cout << pos_inner_edge_upper_aureole[OIL] << "\t";
				cout << pos_max_rate_upper_aureole[OIL] << "\t";
				cout << temp_max_rate_upper_aureole[OIL] << "\t";
			}
			cout << pos_outer_edge_lower_aureole[OIL] << "\t";
			cout << pos_inner_edge_lower_aureole[OIL] << "\t";
			cout << pos_max_rate_lower_aureole[OIL] << "\t";
			cout << temp_max_rate_lower_aureole[OIL] << "\t";
			if ( !sheet_is_dyke ) {
				cout << pos_outer_edge_upper_aureole[VITRINITE] << "\t";
				cout << pos_inner_edge_upper_aureole[VITRINITE] << "\t";
				cout << pos_max_rate_upper_aureole[VITRINITE] << "\t";
				cout << temp_max_rate_upper_aureole[VITRINITE] << "\t";
			}
			cout << pos_outer_edge_lower_aureole[VITRINITE] << "\t";
			cout << pos_inner_edge_lower_aureole[VITRINITE] << "\t";
			cout << pos_max_rate_lower_aureole[VITRINITE] << "\t";
			cout << temp_max_rate_lower_aureole[VITRINITE] << "\t";
			if ( !sheet_is_dyke ) {
				cout << vit_outer_edge_upper_aureole[LABILE] << "\t";
				cout << vit_inner_edge_upper_aureole[LABILE] << "\t";
				cout << vit_max_rate_upper_aureole[LABILE] << "\t";
			}
			cout << vit_outer_edge_lower_aureole[LABILE] << "\t";
			cout << vit_inner_edge_lower_aureole[LABILE] << "\t";
			cout << vit_max_rate_lower_aureole[LABILE] << "\t";
			if ( !sheet_is_dyke ) {
				cout << vit_outer_edge_upper_aureole[REFRACTORY] << "\t";
				cout << vit_inner_edge_upper_aureole[REFRACTORY] << "\t";
				cout << vit_max_rate_upper_aureole[REFRACTORY] << "\t";
			}
			cout << vit_outer_edge_lower_aureole[REFRACTORY] << "\t";
			cout << vit_inner_edge_lower_aureole[REFRACTORY] << "\t";
			cout << vit_max_rate_lower_aureole[REFRACTORY] << "\t";
			if ( !sheet_is_dyke ) {
				cout << vit_outer_edge_upper_aureole[OIL] << "\t";
				cout << vit_inner_edge_upper_aureole[OIL] << "\t";
				cout << vit_max_rate_upper_aureole[OIL] << "\t";
			}
			cout << vit_outer_edge_lower_aureole[OIL] << "\t";
			cout << vit_inner_edge_lower_aureole[OIL] << "\t";
			cout << vit_max_rate_lower_aureole[OIL] << "\t";
			cout << aureole_dless[LABILE] << "\t";
			cout << aureole_dless[REFRACTORY] << "\t";
			cout << aureole_dless[OIL] << "\t";
			cout << aureole_dless[VITRINITE] << "\t";
			cout << mass_from_labile_dless_aureole  << "\t";
			cout << mass_from_refractory_dless_aureole  << "\t";
			cout << mass_from_oil_dless_aureole  << "\n";

// NB These not yet corrected for mass of methane at time of injection

		} else if ( print_mass_expelled ) {
			cout << fixed << showpoint;
			cout << setprecision(5) << time/time_scale << "\t";
			cout.unsetf(ios_base::floatfield);
			cout <<	 mass_methane_vented << "\t"
						<< mass_methane_generated << "\t"
						<< max_mass_methane_vented << "\t"
						<< flux_methane_vented << "\t"
						<< flux_methane_generated << "\t"
						<< max_flux_methane_vented << "\n";

		}
	}

}




double get_progress ( int i_reaction, int i_grid )
{
/*
	Sums the approximate reactions to get the mean reaction progress.
*/
	double progress = 0.0;
	for ( int i_approx=1; i_approx<=n_approx_reactions[i_reaction]; i_approx++ )
		progress += progress_of_reactions[i_approx][i_reaction][i_grid] / n_approx_reactions[i_reaction];
	return progress;
}

double get_rate ( int i_reaction, int i_grid )
{
/*
	Sums the approximate reaction rates to get the mean reaction rate.
*/
	double rate = 0.0;
	for ( int i_approx=1; i_approx<=n_approx_reactions[i_reaction]; i_approx++ )
		rate += rate_of_reactions[i_approx][i_reaction][i_grid] / n_approx_reactions[i_reaction];
	return rate;
}




void write_out_profile ( double profile_time )
{

	for (int i=1; i<=n_grid_points; i++) {

/* PRINT PROFILE IF NECESSARY */

		double depth = (sheet_is_dyke) ? depth_of_sheet : dist[i];
		double distance = (sheet_is_dyke) ? dist[i] : dist[i]-depth_of_sheet;
		double distance_from_sheet_edge = ( distance<0.0 ) ? distance+sheet_half_thickness : distance-sheet_half_thickness;
		
		if (i==1) cerr << " Writing profile perpendicular to sheet\n";

		if ( print_temperature ) {
			if (i==1) cerr << "   column 1, distance from sheet centre (m)\n";
			if (i==1) cerr << "   column 2, temperature (C)\n";
			if (i==1) cerr << "   column 3, distance from surface (m)\n";
			if (i==1) cerr << "   column 4, time after magma injection (yr)\n";
			if (i==1) cerr << "   column 5, distance from sheet margin (m)\n\n";
			cout	<< distance << "\t"
					<< temperature[i] << "\t"
					<< depth << "\t" 
					<< profile_time/time_scale << "\t"
					<< distance_from_sheet_edge << "\n"; 

		} else if ( print_reactions ) {
			if (i==1) {
				cerr << "   column 1, distance from sheet centre (m)\n";
				cerr << "   column 2, labile kerogen reaction\n";
				cerr << "   column 3, refractory kerogen reaction\n";
				cerr << "   column 4, vitrinite reaction\n";
				cerr << "   column 5, oil cracking reaction\n";
				cerr << "   column 6, distance from surface (m)\n";
				cerr << "   column 7, temperature (C)\n";
				cerr << "   column 8, distance from sheet margin (m)\n\n";
			}
			cout << distance << "\t";
			if ( abs(distance) < sheet_half_thickness )	
				cout << "nan\t";
			else
				cout << aureole_model[LABILE][i] << "\t";
			if ( abs(distance) < sheet_half_thickness )	
				cout << "nan\t";
			else
				cout << aureole_model[REFRACTORY][i] << "\t";
			if ( abs(distance) < sheet_half_thickness )	
				cout << "nan\t";
			else
				cout << aureole_model[VITRINITE][i] << "\t";
			if ( abs(distance) < sheet_half_thickness )	
				cout << "nan\t";
			else
				cout << aureole_model[OIL][i] << "\t";
			cout	<< depth << "\t" 
					<< temperature[i] << "\t"
					<< distance_from_sheet_edge << "\n"; 

		} else if ( print_flux ) {
			cout	<< distance << "\t" 
					<< aureole_model[LABILE][i] << "\t"
					<< aureole_model[REFRACTORY][i] << "\t"
					<< aureole_model[OIL][i] << "\n";

		} else if ( print_vitrinite ) {
			if (i==1) cerr << "   column 1, distance from sheet centre (m)\n";
			if (i==1) cerr << "   column 2, vitrinite reflectance (%)\n";
			if (i==1) cerr << "   column 3, distance from surface (m)\n";
			if (i==1) cerr << "   column 4, time after magma injection (yr)\n";
			if (i==1) cerr << "   column 5, distance from sheet margin (m)\n";
			cout	<< distance << "\t";
			if ( abs(distance) < sheet_half_thickness )	
				cout << "nan\t";
			else
				cout << aureole_model[VITRINITE][i] << "\t";
			cout	<< depth << "\t" 
					<< profile_time/time_scale << "\t" 
					<< distance_from_sheet_edge << "\n"; 
		}

	}
}





void set_up_distance_grid ( void ) 
{
  cerr << " Setting up distance grid\n";
/* 
	DISTANCES ARE IN METRES RELATIVE TO THE CENTRE LINE OF THE MAGMA SHEET.
	CALCULATIONS ARE EASIER WITH NODES EXACTLY ON THE CENTRE LINE AND THE
	SHEET MARGINS.  THEREFORE THE NODE SPACING IS DEFINED BY STATING THE 
	NUMBER OF NODES IN THE SHEET, AND THE SAME SPACING IS USED FOR THE
	PART OF THE GRID IN THE COUNTRY ROCK.  
 */

//	grid_point_wall = grid_points_in_sheet + 1;
	dist_increment = double(sheet_half_thickness/grid_points_in_sheet);
//	double dist_increment_sq = pow(dist_increment,2.0);
	n_grid_points = grid_points_in_sheet + grid_points_in_country_rock + 1;
  cerr << "   Grid point spacing is " << dist_increment << " m\n";

/* 
	IF THE SHEET IS A DYKE, THEN THE CALCULATION GRID IS HORIZONTAL.  
	THE CALCULATION IS PERFORMED IN ONE OF THE TWO WALLS BECAUSE THE
	THERMAL AUREOLE IS SYMMETRICAL.
 */

	if ( sheet_is_dyke ) {
		grid_point_centre = 1;
		dist = dvector(1,n_grid_points);
		for (int i=1; i<=n_grid_points; i++) 
			dist[i] = double(i-1)*dist_increment;
  	cerr << "   Sheet is a dyke of thickness " << sheet_thickness << " m at depth " << depth_of_sheet << "m\n";
  	cerr << "   Centreline at 0 m, margins at ±"<<dist[grid_points_in_sheet+1]<<" m (though only +ve half modelled)\n";

/* 
	IF THE SHEET IS A SILL, THEN THE CALCULATION GRID IS VERTICAL
	AND THE TEMPERATURE WITHIN THE GRID DEPENDS ON THE GEOTHERMAL
	GRADIENT AS WELL AS THE EFFECT OF THE SILL.  THEREFORE THE
	THERMAL AUREOLE IS NOT SYMMETRICAL AND THE CALCULATION MUST
	BE PERFORMED IN BOTH WALLS: NEGATIVE DEPTHS ARE ABOVE THE
	CENTRE OF THE SILL AND POSITIVE DEPTHS BELOW IT.
 */

	} else {
		grid_point_centre = n_grid_points;
		n_grid_points += n_grid_points-1;
		dist = dvector(1,n_grid_points);
		while ( (depth_of_sheet+(1-grid_point_centre)*dist_increment) < 0.0 ) {
			grid_point_centre--;
			n_grid_points--;
		}
		for (int i=1; i<=n_grid_points; i++) 
			dist[i] = depth_of_sheet + double(i-grid_point_centre)*dist_increment;
			
//		dist[n_grid_points] += 1000.0;
			
  	cerr << "   Sheet is a sill of thickness " << sheet_thickness << " m, centreline at depth " << depth_of_sheet << " m\n";
  	cerr << "   Grid goes from " << dist[1] <<" m to " << dist[n_grid_points] << " m\n";

	}
}


void recreate_distance_grid ( void ) 
{
  cerr << " Re-creating existing distance grid\n";
/* 
	THE EXISTING FILE MUST BE IN THE FORMAT PRODUCED BY SILLBURP 
  WITH -Pr
   column 1, distance from sheet centre (m)
   column 2, labile kerogen reaction
   column 3, refractory kerogen reaction
   column 4, vitrinite reaction
   column 5, oil cracking reaction
   column 6, distance from surface (m)
   column 7, temperature (C)
 */

	double starting_sheet_thickness = 0.0;
  double new_sheet_thickness = sheet_thickness;
	int n_starting_grid_points = n_grid_points;
  int grid_point_upper_margin;
  int grid_point_lower_margin;

/*  */

	dist_increment = data_starting[0][2] - data_starting[0][1];
//	double dist_increment_sq = pow(dist_increment,2.0);
  cerr << "   Grid point spacing is " << dist_increment << "m\n";

/* 
	SEE WHETHER NEW SHEET THICKNESS IS MULIPLE OF GRID SPACING
  STOP IF NOT
 */

	double remainder = (sheet_thickness/dist_increment) - int(sheet_thickness/dist_increment);
	if ( remainder != 0.0 ) {
		cerr << " ERROR: New sheet thickness is not a multiple of the grid spacing\n";
		cerr << "        Try with new sheet thickness of " << int(sheet_thickness/dist_increment)*dist_increment <<"\n";
		cerr << "         or with new sheet thickness of " << (int(sheet_thickness/dist_increment)+1)*dist_increment <<"\n";
		cerr << "         or change grid spacing in original sheet\n";
		exit(0);
	}
	int grid_points_in_new_half_sheet = new_sheet_thickness / dist_increment / 2;
	if (grid_points_in_new_half_sheet%2 == 0) {
		cerr << " ERROR: Even number of grid points in new sheet\n";
		exit(0);
	}

/* 	DECIDE WHETHER SHEET IS SILL OR DYKE
		IF THE SHEET IS A DYKE, THE FIRST DISTANCE VALUE WILL BE ZERO 
 */

	int grid_points_in_starting_sheet = 0;
	if ( data_starting[0][1] == 0.0 ) {
		if ( !sheet_is_dyke ) cerr << "   WARNING: Re-setting sheet orientation to dyke\n";
		sheet_is_dyke = true;

/*
		COUNT NODES THE STARTING SHEET & READ OFF STARTING SHEET THICKNESS
 */

		grid_point_centre = 1;
		grid_points_in_starting_sheet = 0;
		while ( is_dnan( data_starting[1][grid_points_in_starting_sheet+1] ) ) grid_points_in_starting_sheet++;
cerr << grid_points_in_starting_sheet << "\n";
		starting_sheet_thickness = 2.0 * data_starting[0][grid_points_in_starting_sheet+1];
		depth_of_sheet = data_starting[5][1];
  	cerr << "   Starting sheet is a dyke of thickness " << starting_sheet_thickness << " m at depth " << depth_of_sheet << "m\n";
  	cerr << "   Starting centreline at 0 m, margins at ±" << data_starting[0][grid_points_in_starting_sheet+1] << " m\n";

/* NUMBER OF NODES IN COUNTRY ROCK */

		grid_points_in_country_rock = n_grid_points - grid_points_in_starting_sheet - 1;
		n_grid_points += grid_points_in_sheet;
		grid_points_in_sheet = grid_points_in_starting_sheet + grid_points_in_sheet;
		sheet_thickness += starting_sheet_thickness;
cerr << sheet_thickness << "\n";
		sheet_half_thickness = 0.5 * sheet_thickness;
cerr << sheet_half_thickness << "\n";

/* SET UP STOREAGE ARRAYS */

		dist = dvector(1,n_grid_points);
		progress_of_reactions = d3tensor(1,n_approx_max,1,n_reactions,1,n_grid_points);
		aureole_model = dmatrix(1,n_reactions,1,n_grid_points);
		aureole_model_prev = dmatrix(1,n_reactions,1,n_grid_points);
		aureole_rates = dmatrix(1,n_reactions,1,n_grid_points);
		aureole_vitrinite = dvector(1,n_grid_points);
		progress_at_injection = dmatrix(1,n_reactions,1,n_grid_points);
		temperature = dvector(1,n_grid_points);
		temperature_increment = dvector(1,n_grid_points);

/* COPY EXISTING DATA INTO ARRAYS */ 

		for ( int i_grid=1; i_grid<=n_grid_points; i_grid++ ) {

/* COPY EXISTING DISTANCE ARRAY AND ADD THICKNESS OF NEW SHEET AT END */
			
			if (i_grid <= n_starting_grid_points)
				dist[i_grid] = data_starting[0][i_grid];
			else
				dist[i_grid] = data_starting[0][n_starting_grid_points] + dist_increment 
												* (i_grid-n_starting_grid_points);

/* TEMPERATURE IS MAGMA TEMPERATURE WITHIN NEW SHEET, COPY TEMPERATURE OUTSIDE IT */

//			if (i_grid <= grid_points_in_sheet) {
			if (dist[i_grid] <= 0.5*new_sheet_thickness) {
cerr << dist[i_grid] << "\t" << sheet_half_thickness << "\n";
				temperature[i_grid] = magma_temperature;
			} else {
				temperature[i_grid] = data_starting[6][i_grid-grid_points_in_new_half_sheet];
			}
			temperature_increment[i_grid] = 0.0;

/* NO REACTIONS INSIDE NEW SHEET, COPY REACTION PROGRESS OUTSIDE IT */ 

			for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ )
				for ( int i_approx=1; i_approx<=n_approx_max; i_approx++ )
//					if (i_grid < grid_points_in_sheet)
					if (dist[i_grid] <= 0.5*new_sheet_thickness)
						progress_of_reactions[i_approx][i_reaction][i_grid] = 0.0;
					else
						progress_of_reactions[i_approx][i_reaction][i_grid] =
							data_starting[i_reaction][i_grid-grid_points_in_new_half_sheet];

		}

  	cerr << "   New sheet is a composite dyke of thickness " << sheet_thickness << " m at depth " << depth_of_sheet << "m\n";
  	cerr << "   New centreline at 0 m, margins at ±" << data_starting[0][grid_points_in_sheet+1] << " m\n";

/* 	SHEET IS A SILL */

	} else { 
		if ( sheet_is_dyke ) cerr << "   WARNING: Re-setting sheet orientation to sill\n";
		sheet_is_dyke = false;

/* CENTRE OF STARTING SHEET */

		grid_point_centre = 1;
		while ( data_starting[0][grid_point_centre] != 0.0) grid_point_centre++;
		depth_of_sheet = data_starting[5][grid_point_centre];
//cerr << grid_point_centre << "\n";

/*
		COUNT HOW MANY NODES THERE ARE IN THE STARTING SHEET
		HENCE WORK OUT NUMBER OF NODES IN COUNTRY ROCK
 */

//cerr << depth_of_sheet << "\n";
//exit(0);
		grid_point_upper_margin = 0;
		while ( !is_dnan( data_starting[1][grid_point_upper_margin] ) ) grid_point_upper_margin++;
		grid_point_lower_margin = grid_point_upper_margin;
		while ( is_dnan( data_starting[1][grid_point_lower_margin] ) ) grid_point_lower_margin++;
		grid_points_in_starting_sheet = grid_point_lower_margin - grid_point_upper_margin + 1;
//cerr << grid_points_in_starting_sheet << "\n";
		starting_sheet_thickness = 2.0 * data_starting[0][grid_point_lower_margin];
  	cerr << "   Starting sheet is a sill of thickness " << starting_sheet_thickness << " m at depth " << depth_of_sheet << "m\n";
  	cerr << "   Starting centreline at 0 m, margins at ±" << data_starting[0][grid_point_lower_margin] << " m\n";

/* NUMBER OF NODES IN COUNTRY ROCK */

		grid_points_in_country_rock = n_grid_points - grid_points_in_starting_sheet - 1;
		int grid_points_in_new_sheet = grid_points_in_sheet;
		n_grid_points += grid_points_in_sheet;
		grid_points_in_sheet += grid_points_in_sheet;
		sheet_thickness += starting_sheet_thickness;
		sheet_half_thickness = 0.5 * sheet_thickness;

/* SET UP STOREAGE ARRAYS */

		dist = dvector(1,n_grid_points);
		progress_of_reactions = d3tensor(1,n_approx_max,1,n_reactions,1,n_grid_points);
		aureole_model = dmatrix(1,n_reactions,1,n_grid_points);
		progress_at_injection = dmatrix(1,n_reactions,1,n_grid_points);
		temperature = dvector(1,n_grid_points);
		temperature_increment = dvector(1,n_grid_points);

/* COPY EXISTING DATA INTO ARRAYS */ 

		for ( int i_grid=1; i_grid<=n_grid_points; i_grid++ ) {

/* COPY EXISTING DISTANCE ARRAY AND ADD THICKNESS OF NEW SHEET AT END */
			
			if (i_grid <= n_starting_grid_points)
				dist[i_grid] = data_starting[5][i_grid] - 0.5*new_sheet_thickness;
			else
				dist[i_grid] = data_starting[5][n_starting_grid_points] + dist_increment 
												* (i_grid-n_starting_grid_points);

/* TEMPERATURE IS MAGMA TEMPERATURE WITHIN NEW SHEET, COPY TEMPERATURE OUTSIDE IT */
/* NO REACTIONS INSIDE NEW SHEET, COPY REACTION PROGRESS OUTSIDE IT */ 
/* NB GRID POINT CENTRE IS OLD GRID POINT CENTRE */

			if (i_grid < grid_point_centre) {
				temperature[i_grid] = data_starting[6][i_grid];
				for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ )
					for ( int i_approx=1; i_approx<=n_approx_max; i_approx++ )
						progress_of_reactions[i_approx][i_reaction][i_grid] =
							data_starting[i_reaction][i_grid];
			} else if (i_grid > grid_point_centre+2*grid_points_in_new_sheet) {
				temperature[i_grid] = data_starting[6][i_grid-2*grid_points_in_new_sheet];
				for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ )
					for ( int i_approx=1; i_approx<=n_approx_max; i_approx++ )
						progress_of_reactions[i_approx][i_reaction][i_grid] =
							data_starting[i_reaction][i_grid-2*grid_points_in_new_sheet];
			} else {
				temperature[i_grid] = magma_temperature;
				for ( int i_reaction=LABILE; i_reaction<=OIL; i_reaction++ )
					for ( int i_approx=1; i_approx<=n_approx_max; i_approx++ )
						progress_of_reactions[i_approx][i_reaction][i_grid] = 0.0;
			}
			temperature_increment[i_grid] = 0.0;

		}

  	cerr << "   New sheet is a composite sill of thickness " << sheet_thickness << " m at depth " << depth_of_sheet << "m\n";
//  	cerr << "   New centreline at 0 m, margins at ±" << data_starting[0][grid_point_centre+2*grid_points_in_new_sheet] << " m\n";
	}
  
}



void set_up_reaction_energies ( void ) {

	double s;
	double E;
	double sigma_root_2;
	int N;
	double fraction;
	double E_0;
	double E_1;
	double right_side;
	double erf_inv;

/*
	The kerogen types are assumed to have a normal distribution of activation energies
  given by the probability distribution P(E)  [eq 1]

                        1              /  -( E - E_input )^2  \
    P(E) dE  =  ------------------ exp | -------------------- | dE
                S_input sqrt[2.PI]     \       2 S_input^2    /

	where E_input and S_input are the activation energy and the standard deviation
	supplied in the parameterization: mean_activation_energy[] and standard_deviation[]
	in sillburp.h

	This distribution is approximated using N discrete reactions, where N is an odd
	number input by the user through -N on the command line.  Each of the N reactions
  covers a range of activation energies between E_0 and E_1 given by  [eq 2]

    / E_1
    |                  1
    |     P(E) dE  =  ---
    |                  N
    / E_0


	The integral in eq 2 evaluates to  [eq 3]

   -1      /   E_av - E  \ | E_1        1
   --- erf | ----------- | |        =   -
    2      \  s.sqrt[2]  / | E_0        N

  Starting with a value for E_0, E_1 can be found using the inverse error function since  [eq 4]

                                   /     /  E_av - E_0  \     2  \
    E_1  =   E  -  s.sqrt[2]  ierf | erf | ------------ |  -  -  |
                                   \     \   s.sqrt[2]  /     N  /

  For the reaction in the centre of the probability distribution, eq 2 becomes  [eq 2a]


    / E_1
    |                  1
    |     P(E) dE  =  ---
    |                 2.N
    / E_av

  and eq 4 becomes  [eq 4a]


    E_1  =  E_av  -  s.sqrt[2]  ierf[ -1/N ]


  Then the mean activation energy for each approximate reaction is found by  [eq 5]


            / E_1
            |              
    E  =  N |     E P(E) dE
            |              
            / E_0


  which evaluates to  [eq 6]


            /     -s         / -(E_av - E)^2 \     E_av     /  E_av - E  \ \ | E_1
   E  =  N  | ---------- exp |---------------|  -  ---- erf | ---------- | | |    
            \ sqrt[2.PI]     \    s sqrt[2]  /      2       \ s.sqrt[2]  / / | E_0


*/

/* MAX NUMBER OF APPROX REACTIONS */

	for ( int i_reaction=1; i_reaction<=n_reactions; i_reaction++ )
		if ( n_approx_reactions[i_reaction] > n_approx_max ) n_approx_max = n_approx_reactions[i_reaction];

/* SET UP MATRIX TO STORE REACTION ENERGIES */

	approx_reaction_energies = dmatrix(1,n_reactions,1,n_approx_max);
  cerr << " Organic maturation reactions\n";

// LOOP OVER THE 4 KEROGEN TYPES

	for ( int i_reaction=1; i_reaction<=n_reactions; i_reaction++ ) {

		s = standard_deviation[i_reaction];
		E = mean_activation_energy[i_reaction];
		sigma_root_2 = s * sqrt(2.0);
		N = n_approx_reactions[i_reaction];
		fraction = 2.0 / N;
		E_0 = 0.0;
		E_1 = 0.0;

/* NUMBER OF THE REACTION IN THE CENTRE OF THE PROBABILITY DISTRIBUTION */

		int n_middle = N/2 + 1;

/* LOOP OVER APPROXIMATE REACTIONS */

		for ( int i_approx=n_middle; i_approx<=N; i_approx++ ) {

/* CENTRAL REACTION: SET THE ACTIVATION ENERGY TO THE MEAN ACTIVATION ENERGY (EQ 4A) */

			if ( i_approx==n_middle ) {
				approx_reaction_energies[i_reaction][i_approx] = E;
				if ( N != 1 ) E_0 = E - sigma_root_2 * ierf( -1.0 / N );
				continue;
			}
			
/* REACTION ON THE EDGE OF THE DISTRIBUTION */

			if ( i_approx==N ) {
				approx_reaction_energies[i_reaction][i_approx] = N * (s / SQRT_TWOPI *
					exp( -pow((E-E_0),2.0)/(2.0*s*s) ) +
					E / 2.0 * ( 1.0 + erf( (E-E_0)/sigma_root_2 ) ) );				

/* OTHER REACTIONS (EQ 4) */

			} else {
				right_side = erf( (E-E_0)/sigma_root_2 ) - fraction;
				erf_inv = ierf( right_side );
				E_1 = E - erf_inv * sigma_root_2;

/* FIND ACTIVATION ENERGY FOR EACH APPROXIMATE REACTION (EQ 6) */

				approx_reaction_energies[i_reaction][i_approx] = N * (-s / SQRT_TWOPI *
					( exp( -pow((E-E_1),2.0)/(2.0*s*s) ) - exp( -pow((E-E_0),2.0)/(2.0*s*s) ) ) -
					E / 2.0 * ( erf( (E-E_1)/sigma_root_2 ) - erf( (E-E_0)/sigma_root_2 ) ) );
			}

/* PROBABILITY DISTRIBUTION IS SYMMETRICAL */

			approx_reaction_energies[i_reaction][N-i_approx+1] = 2.0*E - approx_reaction_energies[i_reaction][i_approx];

/* REMEMBER ACTIVATION ENERGY FOR NEXT APPROXIMATE REACTION */

			E_0 = E_1;

		} // END OF APPROXIMATE REACTION LOOP

	if ( i_reaction == 1 )
	  cerr << "   1 Labile kerogen,";
	else if ( i_reaction == 2 )
	  cerr << "   2 Refractory kerogen,";
	else if ( i_reaction == 3 )
	  cerr << "   3 Vitrinite,";
	else if ( i_reaction == 4 )
	  cerr << "   4 Oil cracking,";
  cerr << " activation energy " << E*1.0e-3
			<< "kJ, standard deviation " << s*1.0e-3 << "kJ\n";
	cerr << "     approximated by " << N 
			<< " reactions from " << approx_reaction_energies[i_reaction][1]*1.0e-3
			<< " to " << approx_reaction_energies[i_reaction][N]*1.0e-3 << "kJ\n";

	} // END OF KEROGEN TYPE LOOP

/*
	for ( int i_reaction=1; i_reaction<=n_reactions; i_reaction++ )
		for ( int i_approx=1; i_approx<=n_approx_reactions[i_reaction]; i_approx++ )
			cerr << "Kerogen type " << i_reaction 
					<< ", Reaction " << i_approx
					<< ", Energy " << approx_reaction_energies[i_reaction][i_approx]*1.0e-3 << "kJ\n";
*/

}






double arrhenius_reaction_rate_coeff( 
	double energy, 
	double temp_c )
{
	double temp_k = temp_c + 273.15;
	double rate = mean_pre_exp_const[reaction] *
		exp( -energy / gas_constant / temp_k );
	return rate;
}


/*
double average_reaction_rate_coeff( 
	double energy, 
	double temp_c_0,
	double temp_c_1)
{
	double temp_k = temp_c + 273.15;
	double rate = mean_pre_exp_const[reaction] *
		exp( -energy / gas_constant / temp_k );
	return rate;
}
*/



double find_dless_solidification_boundary( double dless_solidification_boundary )
{
	double right_side = exp(-dless_solidification_boundary*dless_solidification_boundary) / dless_solidification_boundary / ( 1.0 + erf( dless_solidification_boundary ) );
	return right_side - left_side;
}



double ierf( double xin )
{
	const int JMAX = 100;
	const double x1 = 0.0; 
	const double x2 = 100.0; 
	const double xacc = 1.0e-10; 

	double sign = (xin<0.0) ? (xin*=-1.0,-1.0) : 1.0;

	if ( xin == 0.0 ) return 0.0;
	if ( xin >= 1.0 ) {
		cerr << "Argument > 1.0 in ierf\n";
		return xin*sign;
	}

	int j;
	double dx,f,fmid,xmid,rtb;

	f = erf(x1)-xin;
	fmid = erf(x2)-xin;
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++) {
		fmid = erf(xmid=rtb+(dx *= 0.5))-xin;
		if (fmid <= 0.0) rtb=xmid;
		if (abs(dx) < xacc || fmid == 0.0) return rtb*sign;
	}
	cerr << "Too many root bisections in ierf\n";
	exit (0);
}




#define JMAX 100
double rtbis( 
	double x1, 
	double x2, 
	double xacc, 
	double (func)(double) )
{
        int j;
        double dx,f,fmid,xmid,rtb;

        f=func(x1);
        fmid=func(x2);
        if (f*fmid >= 0.0) {
			cerr << "Root must be bracketed for bisection in rtbis\n";
			exit (0);
		}
        rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
        for (j=1;j<=JMAX;j++) {
                fmid=func(xmid=rtb+(dx *= 0.5));
                if (fmid <= 0.0) rtb=xmid;
                if (abs(dx) < xacc || fmid == 0.0) return rtb;
        }
        cerr << "Too many bisections in rtbis\n";
        return 0.0;
}
#undef JMAX





void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
	int i,k;
	double p,qn,sig,un,*u;

	u=dvector(1,n-1);
	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free_dvector(u,1,n-1);
}



void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
	int klo,khi,k;
	double h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0)
		cerr << " ERROR: Bad xa input to routine splint\n\n";
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}


#define MAXIT 100
#define EULER 0.5772156649
#define FPMIN 1.0e-30
#define EPS 1.0e-7
double expint(
	int n, 
	double x)
{
	int i,ii,nm1;
	double a,b,c,d,del,fact,h,psi,ans;
	
	nm1=n-1;
	if (n < 0 || x < 0.0 || (x==0.0 && (n==0 || n==1))) {
		cerr << "bad arguments in expint\n";
		exit(0);
	} else {
		if (n == 0) ans=exp(-x)/x;
		else {
			if (x == 0.0) ans=1.0/nm1;

			else {
				if (x > 1.0) {
					b=x+n;
					c=1.0/FPMIN;
					d=1.0/b;
					h=d;
					for (i=1;i<=MAXIT;i++) {
						a = -i*(nm1+i);
						b += 2.0;
						d=1.0/(a*d+b);
						c=b+a/c;
						del=c*d;
						h *= del;
						if (fabs(del-1.0) < EPS) {
							ans=h*exp(-x);
							return ans;
						}
					}
					cerr << "continued fraction failed in expint\n";
					exit(0);
				} else {
					ans = (nm1!=0 ? 1.0/nm1 : -log(x)-EULER);
					fact=1.0;
					for (i=1;i<=MAXIT;i++) {
						fact *= -x/i;
						if (i != nm1) del = -fact/(i-nm1);
						else {
							psi = -EULER;
							for (ii=1;ii<=nm1;ii++) psi += 1.0/ii;
							del=fact*(-log(x)+psi);
						}
						ans += del;
						if (fabs(del) < fabs(ans)*EPS) return ans;
					}
					cerr << "series failed in expint\n";
					exit(0);
				}
			}
		}
	}
	return ans;
}
#undef MAXIT
#undef EPS
#undef FPMIN
#undef EULER


/*
void avevar(double data[], int n, double *ave, double *var)
{
        int j;
        double s,ep;

        for (*ave=0.0,j=1;j<=n;j++) *ave += data[j];
        *ave /= n;
        *var=ep=0.0;
        for (j=1;j<=n;j++) {
                s=data[j]-(*ave);
                ep += s;
                *var += s*s;
        }
        *var=(*var-ep*ep/n)/(n-1);
}

*/




void get_user_parameters( int argc, char **argv )
{
	if ( argc == 1 ) {
cerr << "                                                               \n";
cerr << " sillburp - temperature and gas production history when a      \n";
cerr << "            magma sheet intrudes mudrock                       \n";
cerr << "                                                               \n";
cerr << "                                                               \n";
cerr << " reaction 1: CH gas direct from labile kerogen \n";
cerr << " reaction 2: CH gas from refractory kerogen \n";
cerr << " reaction 3: vitrinite reflectance \n";
cerr << " reaction 4: CH gas cracked from oil from labile kerogen \n";
cerr << "                                                               \n";
cerr << " -A<temp>:  magma temperature ["<<magma_temperature<<" C]      \n";
cerr << "                                                               \n";
cerr << " -C<T>/<t>:  constant temperature <T> beginning at time <t>    \n";
cerr << "             before magma injection                            \n";
cerr << "                                                               \n";
cerr << " -d: all input and output times in days [years]                \n";
cerr << "                                                               \n";
cerr << " -D<density>: density of country rock ["<<density_country_rock<<" Mt/m^3]\n";
cerr << "                                                               \n";
//cerr << " -e[<form_time>/<decay_time>/<init_flux>/<bkgrnd_flux>/<CH4_conc>]:\n";
cerr << " -E[<form_time>/<decay_time>/<init_flux>/<bkgrnd_flux>/<CH4_conc>]:\n";
cerr << "   Calculate total mass of methane expelled (using -E)         \n";
//cerr << "   or rate of methane expulsion (using -e).                    \n";
//cerr << "   Use -e to generate input to 'lipburp'.                     \n";
cerr << "   Parameters:                                                 \n";
cerr << "    <form_time>, time of vent formation after sill intrusion   \n";
cerr << "    <decay_time>, characteristic time for decay of vent flux   \n";
cerr << "    <init_flux>, initial Darcy flux                            \n";
cerr << "    <bkgrnd_flux>, background Darcy flux                       \n";
cerr << "    <CH4_conc>, concentration of methane in venting fluid      \n";
cerr << "   Omit all 5 parameters to expel the methane instantaneously. \n";
cerr << "                                                               \n";
cerr << " -F  With -I, only print misfit to standard output.            \n";
cerr << "                                                               \n";
cerr << " -G<dTdt>[/<T0>]: Geothermal gradient controls                 \n";
cerr << "   <dTdz>: burial heating rate (C/yr) [";
cerr << heating_rate_from_burial;
cerr <<                                         "]\n";
cerr << "   <T0>: temperature at the sediment surface (C) [";
cerr << surface_temperature;
cerr <<                                                   "]\n";
cerr << "   For constant temperature T use -G0/T                        \n";
cerr << "                                                               \n";
cerr << " -g<T0>/<dTdz>/<dzdt>: Geothermal gradient controls            \n";
cerr << "   <T0>: temperature at the sediment surface (C) [";
cerr << surface_temperature;
cerr <<                                                   "]\n";
cerr << "   <dTdz>: geothermal gradient (C/m) [";
cerr << geothermal_gradient;
cerr <<                                                   "]\n";
cerr << "   <dzdt>: burial rate (m/yr) [";
cerr << burial_rate;
cerr <<                                "]\n";
cerr << "                                                               \n";
cerr << " -H<file>: Starting temperature history supplied in <file>     \n";
cerr << "   in same format as output using -Pr.  The grid spacing       \n";
cerr << "   in -X must be compatible with the input file or it will     \n";
cerr << "   be reset.  Cannot use -H with -G.                           \n";
cerr << "                                                               \n";
cerr << " -I<file>: Vitrinite reflectance profile supplied in <file>.   \n";
cerr << "   Report misfit between model and data.                       \n";
cerr << "                                                               \n";
cerr << " -J<z>: Window surrounding sheet.  [";
cerr << aureole_window;
cerr <<                                			"]\n";
cerr << "   If -I is supplied:\n";
cerr << "     Window surrounding sheet in multiples of sheet thickness.\n";
cerr << "     Profile data outside this window will be used to find the   \n";
cerr << "     geothermal gradient (if a sill) or background maturity \n";
cerr << "     (if a dyke).\n";
cerr << "   If -I is not supplied and sheet is a sill:\n";
cerr << "     Integrate total methane produced above depth <z>\n";
cerr << "                                                               \n";
cerr << " -K<A>/<E>/<sigma>: Give alternative reaction kinetic parameters.\n";
cerr << "   <A> is the frequency factor                                 \n";
cerr << "   <E> is the activation energy                                \n";
cerr << "   <sigma> is the standard deviation of activation energies    \n";
cerr << "   These parameters are used to replace the pre-defined        \n";
cerr << "   kinetic parameters for the vitrinite reaction.              \n";
cerr << "                                                               \n";
cerr << " -L<n>: Scheme for accounting for latent heat of crystallization\n";
cerr << "   <n> can be 0: [default] analytical scheme in Turcotte &     \n";
cerr << "                 Schubert with L=320,000J/kg and Cp=1200/kg/C  \n";
cerr << "              1: Scheme of Maclennan & Lovell (Geology, 2002)  \n";
cerr << "                                                               \n";
cerr << " -M<time>[/<radius>]: mass of methane produced by circular sill\n";
cerr << "   of given radius (in km) up to given time                    \n";
cerr << "   output: time (yr), mass from reactions 1, 2, 4 (kg), total mass (kg)\n";
cerr << "   If <radius> is not supplied, the area density (kg/m^2) is   \n";
cerr << "   printed; NB this value has been multiplied by PI.           \n";
cerr << "   The mass of methane already generated at the time of sheet  \n";
cerr << "   injection is subtracted from the printed value; choose -m   \n";
cerr << "   to print the total mass generated                           \n";
cerr << "                                                               \n";
cerr << " -N<nl>[/nr/nv]: Accuracy of organic reaction calculation      \n";
cerr << "   Each <n> is an odd number of individial reactions used to   \n";
cerr << "   approximate the probability distribution of activation      \n";
cerr << "   energies                                                    \n";
cerr << "   <nl> ["<<n_approx_reactions[LABILE]<<
                  "] is the number of reactions for labile kerogen, and \n";
cerr << "   this number is used by default for the other reactions      \n";
cerr << "   if <nr> and <nv> are not specified                          \n";
cerr << "   <nr> ["<<n_approx_reactions[REFRACTORY]<<
                  "] is the number of reactions for refractory kerogen  \n";
cerr << "   <nv> ["<<n_approx_reactions[VITRINITE]<<
                  "] is the number of reactions for vitrinite           \n";
cerr << "                                                               \n";
cerr << " -o<time_step>:       Output time step                         \n";
cerr << " -O<time_step_power>: Output time step as 10^<time_step_power> \n";
cerr << "   Output is printed at the calculation time step after        \n";
cerr << "   the given output time step has elapsed (i.e. not normally   \n"; 
cerr << "   at intervals of <time_step>)                                \n";
cerr << "                                                               \n";
cerr << " -P<type><time>[/<node>]:  Print profile or time series        \n";
cerr << "   <type> can be a: aureole width (m), temp (C) and reaction rate (/yr)\n";
cerr << "                 q: methane flux density (Mt/kg/yr)            \n";
cerr << "                 r: concentration of reaction product          \n";
cerr << "                 t: temperature (C)                            \n";
cerr << "                 v: vitrinite reflectance (%)                  \n";
cerr << "   <time> is the age after sill intrusion                      \n";
cerr << "   If <node> is not supplied, prints profile perpendicular to  \n";
cerr << "   sheet at the end <time>, otherwise prints the value(s) of   \n";
cerr << "   <type> at node number <node> for every time step.  <node>   \n";
cerr << "   is ignored for option 'a' because the appropriate node(s)   \n";
cerr << "   are to be determined from the calculated reaction profiles  \n";
cerr << "                                                               \n";
cerr << " -R<starting_reflectance>/<max_reflectance>:  ["<<reflectance_min<<"/"<<reflectance_max<<"]\n";
cerr << "   Used for calculating vitrinite reflectance.\n";
cerr << "													\n";
cerr << " -S: Output sheet width (m), solidification time (yr), wall temp. (C)\n";
cerr << "                                                               \n";
cerr << " -t<first_time>: First time to write out (yr) ["<<first_time_to_output<<"]\n";
cerr << "                                                               \n";
cerr << " -T<ta>/<tb>:  time steps (yr)                                 \n";
//cerr << "   <tp>: time step to print out                                \n";
cerr << "   <ta>: max time step for calculation after magma injection [";
cerr << time_step_after;
cerr <<                                                               "]\n";
cerr << "   <tb>: max time step for calculation before magma injection [";
cerr << time_step_before;
cerr <<                                                            "]   \n";
cerr << "   NB These timesteps can be reduced automatically to ensure   \n";
cerr << "   stability                                                   \n";
cerr << "                                                               \n";
cerr << " -V: Do vitrinite reflection reaction                          \n";
cerr << "                                                               \n";
cerr << " -W<mass_fraction_labile>/<mass_fraction_refractory>:  \n";
cerr << "   Fraction of carbon available for conversion to methane\n";
cerr << "   in the form of labile & refractory kerogen ["<<mass_fraction_labile<<"/"<<mass_fraction_refractory<<"]	\n";
cerr << "                                                               \n";
cerr << " -X<n1>/<n2>: Calculation grid                                 \n";
cerr << "   <n1> is the number of grid points in HALF the magma sheet [";
cerr << grid_points_in_sheet;
cerr <<                                                                "]\n";
cerr << "   <n2> is the number of grid points in the country rock [";
cerr << grid_points_in_country_rock;
cerr <<                                                           "]    \n";
cerr << "                                                               \n";
cerr << " -Y: Magma sheet is a dyke [default is a sill]                \n";
cerr << "                                                               \n";
cerr << " -Z<thickness>/<depth>: full thickness and intrusion depth of  \n";
cerr << "   the magma sheet (m)                                         \n";
cerr << "                                                               \n";
    exit(0);			
	}

	bool gotC = false;
	bool gotg = false;
	bool gotG = false;
	bool gotH = false;
  bool gotI = false;
  bool gotJ = false;
  bool gotK = false;
  bool gotM = false;
  bool got_o = false;
  bool gotO = false;
  bool gotP = false;
  bool gott = false;
  bool gotT = false;
	bool gotX = false;
	bool gotZ = false;
	bool input_error = false;
	char input[32];
	char *ptr;
	int entry;
	
  for (int i=1; i<argc; i++){
		if (argv[i][0] == '-'){
			switch(argv[i][1]){
			case 'a':
				use_analytical_solution = true;
				break;

			case 'A':
				if (sscanf(&argv[i][2], "%lf", &magma_temperature) != 1){
					cerr << " Invalid -A option\n";
					input_error = true;
				}
				break;

			case 'C':
				gotC = true;
				if (sscanf(&argv[i][2], "%lf/%lf", &constant_temperature, &begin_constant_temperature) != 2){
					cerr << " Invalid -C option\n";
					input_error = true;
				}
				constant_temperature_reaction = true;
				break;

			case 'd':
				time_scale = S_PER_DAY;
				time_unit = "day";
				break;

			case 'D':
				if (sscanf(&argv[i][2], "%lf", &density_country_rock) != 1){
					cerr << " Invalid -D option\n";
					input_error = true;
				}
				break;

			case 'E':
				print_mass_expelled = true;
				if (sscanf(&argv[i][2], "%lf/%lf/%lf/%lf/%lf", &vent_formation_time, &vent_decay_time, &vent_initial_flux, &vent_background_flux, &vent_methane_saturation) == EOF) {
					vent_methane_instantaneously = true;
					break;
				}
				if (sscanf(&argv[i][2], "%lf/%lf/%lf/%lf/%lf", &vent_formation_time, &vent_decay_time, &vent_initial_flux, &vent_background_flux, &vent_methane_saturation) != 5){
					cerr << " Invalid -E option\n";
					input_error = true;
				}
				if ( print_mass_generated ) print_mass_generated = false;
				break;

			case 'F':
				print_misfit = true;
				break;

			case 'g':
				gotg = true;
				if (sscanf(&argv[i][2], "%lf/%lf/%lf", &surface_temperature, &geothermal_gradient, &burial_rate) != 3){
					cerr << " Invalid -g option\n";
					input_error = true;
				}
				if ( burial_rate <= 0.0 ) {
					cerr << " Invalid -g option: burial rate must be > zero\n";
					input_error = true;
				}
				heating_rate_from_burial = geothermal_gradient*burial_rate;
				break;

			case 'G':
				gotG = true;
				if (sscanf(&argv[i][2], "%lf/%lf", &heating_rate_from_burial, &surface_temperature) > 2){
					cerr << " Invalid -G option\n";
					input_error = true;
				}
				if ( heating_rate_from_burial < 0.0 ) {
					cerr << " Invalid -G option: heating rate must be > zero\n";
					input_error = true;
				}
				break;

			case 'H':
				gotH = true;
				strncpy( file_starting_temperature, &argv[i][2], BUFSIZ );
			  if ((fp_starting_temperature = fopen (&argv[i][2], "r")) == NULL) {
					cerr << " Invalid -H option: Cannot open file " << &argv[i][2] << "\n";
					input_error = true;
				}
				break;

			case 'I':
				gotI = true;
				compare_vitrinite_data = true;
				strncpy( file_aureole, &argv[i][2], BUFSIZ );
			  if ((fp_aureole = fopen (&argv[i][2], "r")) == NULL) {
					cerr << " Invalid -I option: Cannot open file " << &argv[i][2] << "\n";
					input_error = true;
				}
				break;

			case 'J':
				gotJ = true;
				use_window = true;
				if (sscanf(&argv[i][2], "%lf", &aureole_window) != 1) {
					cerr << " Invalid -J option\n";
					input_error = true;
				}
				break;

			case 'K':
				gotK = true;
				double dumb1,dumb2,dumb3;
				if (sscanf(&argv[i][2], "%lf/%lf/%lf", &mean_pre_exp_const[3], &mean_activation_energy[3], &standard_deviation[3]) != 3){
					cerr << " Invalid -K option\n";
					input_error = true;
				}
				if ( burial_rate <= 0.0 ) {
					cerr << " Invalid -K option: all parameters must be positive\n";
					input_error = true;
				}
				do_vitrinite_reaction = true;
				do_labile_reaction = false;
				do_refractory_reaction = false;
				break;

			case 'L':
				if (sscanf(&argv[i][2], "%i", &latent_heat_scheme) != 1){
					cerr << " Invalid -L option\n";
					input_error = true;
				}
				break;

			case 'm':
				subtract_mass_at_injection = false;
			case 'M':
				gotM = true;
				print_time_series = true;
				print_aureole = true;
				if (sscanf(&argv[i][2], "%lf/%lf", &end_time, &sill_radius) == 1) {
					print_mass_generated = true;
				} else if (sscanf(&argv[i][2], "%lf/%lf", &end_time, &sill_radius) == 2) {
					print_mass_generated = true;
					scale_to_whole_sill = true;
					sill_radius *= 1000.0;
				} else {
					cerr << " Invalid -M option\n";
					input_error = true;
				}
				if ( print_mass_expelled ) print_mass_generated = false;
				break;

			case 'N':
				if (sscanf(&argv[i][2], "%i/%i/%i", &n_approx_reactions[LABILE], &n_approx_reactions[REFRACTORY], &n_approx_reactions[VITRINITE]) != 3){
					cerr << " Invalid -N option\n";
					input_error = true;
				}
				if ( n_approx_reactions[LABILE]%2 != 1 ) {
					cerr << " -N option: increasing number by 1 to make it odd\n";
					n_approx_reactions[LABILE]++;
				}
				n_approx_reactions[OIL] = n_approx_reactions[LABILE];
				break;

			case 'o':
				got_o = true;
				use_output_time_step = true;
				if (sscanf(&argv[i][2], "%lf", &output_time_step) != 1){
					cerr << " Invalid -o option\n";
					input_error = true;
				}
				if ( output_time_step <= 0 ) {
					cerr << " Invalid -o option: <time_step> must be positive\n";
					input_error = true;
				}
				break;
			case 'O':
				gotO = true;
				use_output_time_step = true;
				use_output_time_step_power = true;
				if (sscanf(&argv[i][2], "%lf", &output_time_step_power) != 1){
					cerr << " Invalid -O option\n";
					input_error = true;
				}
				if ( output_time_step_power <= 0 ) {
					cerr << " Invalid -O option: <time_step_power> must be positive\n";
					input_error = true;
				}
				break;

			case 'P':
				gotP = true;
				print_profile = true;
				switch(argv[i][2]){
				case 'A':
				case 'a':
					print_aureole = true;
					break;
				case 'Q':
				case 'q':
					print_flux = true;
					break;
				case 'R':
				case 'r':
					print_reactions = true;
					break;
				case 'T':
				case 't':
					print_temperature = true;
					break;
				case 'V':
				case 'v':
					print_vitrinite = true;
					do_vitrinite_reaction = true;
					break;
				default:
					cerr << " Invalid -P option: <type> must be one of a,m,r,t,v\n";
					input_error = true;
					break;
				}
				input[0] = 0;
				strcpy (input, &argv[i][3]); 
				ptr = strtok (input, "/");
				entry = 1;
				while (ptr) {
				 	if (entry==1) end_time = atof (ptr);
				 	if (entry==2) {
						print_time_series_at_point = true;
						print_time_series = true;
						n_print = atoi (ptr);
					}
					ptr = strtok (CNULL, "/");
					entry++;
				}
				if ( entry > 3 ) {
					cerr << " Invalid -P option: too many numbers supplied\n";
					input_error = true;
				}
				break;

			case 'R':
				if (sscanf(&argv[i][2], "%lf/%lf", &reflectance_min, &reflectance_max) != 2){
					cerr << " Invalid -R option\n";
					input_error = true;
				}
				break;

			case 'S':
				print_solidification_time = true;
				break;

			case 't':
				gott = true;
				if (sscanf(&argv[i][2], "%lf", &first_time_to_output) != 1){
					cerr << " Invalid -t option\n";
					input_error = true;
				}
				break;

			case 'T':
				gotT = true;
				if (sscanf(&argv[i][2], "%lf/%lf", &time_step_after, &time_step_before) != 2){
					cerr << " Invalid -T option\n";
					input_error = true;
				}

			case 'V':
				do_vitrinite_reaction = true;
				break;

			case 'W':
				if (sscanf(&argv[i][2], "%lf/%lf", &mass_fraction_labile, &mass_fraction_refractory) != 2){
					cerr << " Invalid -W option\n";
					input_error = true;
				}
				if ( mass_fraction_labile > 0.0) do_labile_reaction = true;
				if ( mass_fraction_refractory > 0.0) do_refractory_reaction = true;
				break;

			case 'X':
				gotX = true;
				if (sscanf(&argv[i][2], "%i/%i", &grid_points_in_sheet, &grid_points_in_country_rock) != 2){
					cerr << " Invalid -X option\n";
					input_error = true;
				}
				break;

			case 'Y':
				sheet_is_dyke = true;
				break;

			case 'Z':
				gotZ = true;
				if (sscanf(&argv[i][2], "%lf/%lf", &sheet_thickness, &depth_of_sheet) != 2){
					cerr << " Invalid -Z option\n";
					input_error = true;
				}
				sheet_half_thickness = 0.5 * sheet_thickness;
				if ( sheet_half_thickness <= 0.0 ) {
					cerr << " Invalid -Z option: Sheet thickness must be > 0 km\n";
					input_error = true;
				}
				break;

			default:
				cerr << " Command line argument " << argv[i] << " not recognized\n";
				input_error = true;
				break;
			}
		}
  }

/*
	if ( gotC && !gotT ) {
		cerr << " Error: Must use -T with -C\n";
		input_error = true;
	}
*/
	if ( gotJ && !gotI && sheet_is_dyke ) {
		cerr << " Warning: -J supplied but will not be used for a dyke\n";
		use_window = false;
	}

	if ( print_misfit && !gotI ) {
		cerr << " Warning: -K supplied but will not be used without with -I\n";
	}

	if ( !gotZ ) {
		cerr << " Error: Must supply -Z option\n";
		input_error = true;
	}

	if ( gotH && gotG ) {
		cerr << " Error: Cannot use -G with -H\n";
		input_error = true;
	}

  if ( print_mass_expelled && !scale_to_whole_sill ) {
		cerr << " Error: Must use -M with sill radius with -E or -e\n";
		input_error = true;
	}
	
	if ( scale_to_whole_sill && sheet_is_dyke ) {
		cerr << " Error: Mass calculation not available for dyke\n";
		input_error = true;
	}

	if ( gotM && gotP ) {
		cerr << " Error: Use only one of -M and -P\n";
		input_error = true;
	}


  if ( input_error ) {
		cerr << "\n";
		exit(0);
	}



/* CHANGE TO WORKING UNITS: LENGTH IN M, TIME IN S */

	first_time_to_output *= time_scale;
 	burial_rate /= time_scale;
	heating_rate_from_burial /= time_scale;

}
