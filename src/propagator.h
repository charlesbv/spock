//
// Includes
#include <math.h>
#include "prop_math.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nrlmsise-00.h"
#include "SpiceUsr.h"


// Constants
#define N_SATS 160
#define N_STORM 20
#define N_GROUND_STATION 100
#define N_ENSEMBLES 600// not used anymore
#define N_FILES_ENSEMBLES 30
#define N_SURFACES 100
#define MAX_DEGREE 100
#define MAX_ORDER 100

#ifndef __PROPAGATOR_H__
#define __PROPAGATOR_H__


typedef struct {

  char            name_ground_station[N_GROUND_STATION][256];
  double          latitude_ground_station[N_GROUND_STATION];
  double          longitude_ground_station[N_GROUND_STATION];
  double          altitude_ground_station[N_GROUND_STATION];
  double          min_elevation_angle_ground_station[N_GROUND_STATION];
  double          ecef_ground_station[N_GROUND_STATION][3];
  int             nb_ground_stations;
  
} GROUND_STATION_T;

typedef struct {

  GROUND_STATION_T   GROUND_STATION;
  //  STORM_T            STORM;
    
} COVERAGE_T;

typedef struct {

  double mu;
  double j2;
  double Clm[MAX_DEGREE+1][MAX_ORDER+1];
  double Slm[MAX_DEGREE+1][MAX_ORDER+1];
  double radius;
  double dlat_map; // size of latitude bins for the map (in degrees) 
  double dlon_map; // size of longitude bins for the map (in degrees)
  double dradius_map; // size of radius bins for the map (in km);
  double  min_lat_map; // min latitude for the map
  double max_lat_map;
  double max_radius_map; // max radius for the map
  double min_radius_map;
  double earth_max_radius_map; // max radius for the map
  double earth_min_radius_map;
  
  int nlat_map;
  int nlon_map;
  int nradius_map;
  
    double dzenith_map; // size of zenith bins for the map (in km);
    int nzenith_map;
    double max_zenith_map; // max zenith for the map
  double min_zenith_map;
    double delev_elt_map; // size of elev_elt bins for the map (in km);
    int nelev_elt_map;
    double max_elev_elt_map; // max elev_elt for the map
  double min_elev_elt_map;
    double dazim_elt_map; // size of azim_elt bins for the map (in km);
    int nazim_elt_map;
    double max_azim_elt_map; // max azim_elt for the map
  double min_azim_elt_map;

    double delev_surf_map; // size of elev_surf bins for the map (in km);
    int nelev_surf_map;
    double max_elev_surf_map; // max elev_surf for the map
  double min_elev_surf_map;
    double dazim_surf_map; // size of azim_surf bins for the map (in km);
    int nazim_surf_map;
    double max_azim_surf_map; // max azim_surf for the map
  double min_azim_surf_map;

  
  double ****gravity_map;
  double *****earth_pressure_map;
  double *radius_map;
  double *lat_map;
  double *lon_map;

  double *zenith_map;
  double *elev_surf_map;
  double *azim_surf_map;
  FILE *file_gravity_map;
    char filename_gravity_map[200];
  FILE *file_earth_pressure_map;
    char filename_earth_pressure_map[200];

}  GRAVITY_T;

typedef struct {

  GRAVITY_T   GRAVITY;
  double      flattening;
  double      radius;
  char            earth_fixed_frame[100]; // Earth body-fixed rotating frame for conversion ECI to ECEF  
} PLANET_T;

typedef struct {
  
  struct nrlmsise_flags flags;


 double nb_gitm_file;
  char array_gitm_file[700][70]; // Their should not be more than 700 files in the GITM directory you're looking at. 
  double array_gitm_date[700][2]; // column 0 is the time in Julian date, second column is the index in array_gitm_file of the file which date is this one. 
  int gitm_index_lower;             // represents the time (in Julian date) of the most recent GITM file
  int gitm_index_higher;            // represents the time (in Julian date) of the subsequent GITM file  
  int is_first_step_of_run;          // matters if using GITM for the computation of the density is drag is used. 1 if it is the first step of the run, 0 otherwise
  double ***density_gitm_right_before; // density of the GITM file which date is right before the epoch of the sc
    double ***density_gitm_right_after;// density of the GITM file which date is right after the epoch of the sc
  double ***longitude_gitm, ***latitude_gitm, ***altitude_gitm; // lon/alt/lat of the GITM files
    int nLons_gitm, nLats_gitm, nAlts_gitm, nVars_gitm; // number of lat/lon/alt/variables in the GITM file  
  int index_altitude_right_below_perigee; // index in altitude_gitm of the altitude that is right below the perigee altitude (calculated at the initialization). Note: if the duration of the run is long (so that the sc looses a good amount of altitude, so like 6 months) then this index is set to 0 so that we go over all the altitudes and not only the ones that start below the perigee calcualted at the initialization.
 
}  ATMOSPHERE_T;


typedef struct {

    
  PLANET_T        EARTH;
  PLANET_T        MOON;
  PLANET_T        SUN;
  ATMOSPHERE_T    ATMOSPHERE;
  SpiceDouble geophs[8];
} PARAMS_T;


typedef struct {

  double sma;             // semi-major axis
  double eccentricity;    // eccentricity 
  double inclination;     // inclination
  double long_an;         // longitude of the ascending node OR RAAN (see output file -> called RAAN)? As used in the code, it's the RAAN (see in propagator.c the function cart2kep with eq (2-82) in Vallado3)
  double w;               // argument of periapsis
  double f;               // true anomaly
  double tp;              // time of periapsis
  double E;               // eccentric anomaly
  double ra;              // right ascension (J2000)
  double w_ave_temp;             // obrit average argument of periapsis
  double w_ave;             // obrit average argument of periapsis

  double sma_ave_temp;             // temporary var obrit average sma
  double sma_ave;             // obrit average sma

  double ecc_ave_temp;             // temporary var obrit average ecc
  double ecc_ave;             // obrit average ecc

  double period;          // orbital period
  double an_to_sc;       // angle an to sc (= w+f)
  double initial_an_to_sc ; // initial an_to_sc (for orbit average computation)
  int ave_increm;            // nb of iteration in one orbit (to compute orbit average elements)
  double zenith;             // zenith angle (angle earth to sc, sc to sun) varies form 0 to 360. 0/360 deg if Sun at zenith of satellite, 90 deg if at horizon of satellite and 'behind’ sc deg, 180 deg if Sun at nadir of sc, 270 deg if Sun at horizon of sc and 'in front' of sc 
} ORBITAL_ELEMENTS_T;

typedef struct {

  double latitude;
  double longitude;
  double altitude;


} GEODETIC_T;

// Attitude options
typedef struct {
  double quaternion_current[4];
  double pitch_current; // pitch at the current time step of the propagation 
  double roll_current; // roll at the current time step of the propagation
  double yaw_current; // yaw at the current time step of the propagation
  int order_pitch_current; // order pitch at the current time step of the propagation 
  int order_roll_current; // order roll at the current time step of the propagation
  int order_yaw_current; // order yaw at the current time step of the propagation

  double *pitch; // pitch for all the time steps of the propagation
  double *roll; // roll for all the time steps of the propagation
  double *yaw; // yaw for all the time steps of the propagation
  double pitch_ini_ensemble;
  double roll_ini_ensemble;
  double yaw_ini_ensemble;
  double pitch_angular_velocity_ensemble;
  double roll_angular_velocity_ensemble;
  double yaw_angular_velocity_ensemble;
  double pitch_sigma_angular_velocity_ensemble;
  double roll_sigma_angular_velocity_ensemble;
  double yaw_sigma_angular_velocity_ensemble;  
  double pitch_sigma_ensemble;
  double roll_sigma_ensemble;
  double yaw_sigma_ensemble;  
  double pitch_for_attitude_ensemble;
  double roll_for_attitude_ensemble;
  double yaw_for_attitude_ensemble;
  double pitch_angular_velocity_constant;
  double roll_angular_velocity_constant;
  double yaw_angular_velocity_constant;
  int    *order_pitch;
  int    *order_roll;
  int    *order_yaw;
  double lvlh_alongtrack_in_body_cartesian[3]; // component of the along-track axis of the LVLH frame in the SC reference system under the cartesian representation
  double lvlh_crosstrack_in_body_cartesian[3]; // component of the cross-track axis of the LVLH frame in the SC reference system under the cartesian representation
  char   attitude_profile[256]; // nadir, sun_pointed, other...
  double **quaternion;
} ATTITUDE_T;


// Surface options
typedef struct {

  double Cd;        // Coefficient of drag
  double Cd_sigma;        // Sigma for the ensembles on the coefficient of drag
  double area;      // area of the surface
  double area_solar_panel; // solar cell area
  double normal[3]; // normal vector to the surface in the SC reference system
  double specular_reflectivity; // specular reflectivity
  double diffuse_reflectivity; // diffusive reflectivity
  double solar_radiation_coefficient; // coefficient for solar radiation (theory used in STK)
  char name_of_surface[100]; // name of the surface
  double power_per_surface; // solar power per area
  double acco_coeff; // accommodation coefficient of surface
  int ieff;
} SURFACE_T;

// Integrator options
typedef struct {

  double Ap_static, f107_static, f107A_static;
  double Ap_hist_static[7];
  double *Ap;                     // magnetic index(daily)
  double **Ap_hist;                     // magnetic index (historical)
  double *f107;                   // Daily average of F10.7 flux
  double *f107A;                  // 81 day average of F10.7 flux
  double *density;                // density from the input density file if the user chooses to input the density like this
  double dt;                     // Integrator timestep
  double dt_pos_neg;                     // Integrator timestep (negative if backward propagation)
  double mass;                   // Mass of spacecraft
  double solar_cell_efficiency;  // Efficiency of the solar cells
  char opengl_filename_solar_power[1000]; // name of file with corners of solar panels, read by opengl
  SURFACE_T surface[N_SURFACES]; // Properties of the surface
  SURFACE_T surface_eff[N_SURFACES]; // Properties of the surface effective
  double nb_surfaces;            // number of surfaces
  int nb_surfaces_eff;
  ATTITUDE_T attitude;           // Attitude of the satellite
  double degree;                 // Gravity degree
  double order;                  // Gravity order
  int include_drag;              // include drag
  int include_moon;              // include moon perturbations
  int include_sun;               // include sun perturbations
  int include_solar_pressure;    // include solar pressure
    int include_earth_pressure;    // include earth pressure
  char shadow[256];               // in umbra/penumbra or light of Earth
  char shadow_moon[256];               // in umbra/penumbra or light of Moon
  char format_density_driver[256]; // dynamic, static or density_file
  double bstar;                  // if the orbit is initialized with a tle then bstar is read from the tle
  int initialize_geo_with_bstar; // 1 if the geometry is initialized reading bstar from the tle and not from a geometry file; 0 othersize
  int index_in_attitude_interpolated;
  int index_in_driver_interpolated;
  int index_in_attitude_interpolated_first; // the first value in the propagation of index_in_attitude_interpolated
  int index_in_driver_interpolated_first; // // the first value in the propagation of index_in_driver_interpolated

  int isGPS;
  double sun_elevation; // value is -999 if the sc is not in light
  double beta_angle;
  FILE                *file_given_output; // this is a file used to output things in particular. Not used unless hard coded
  char                filename_given_output[1000]; // this is a file used to output things in particular. Not used unless hard coded
  int sc_main_nb;
  int sc_ensemble_nb;
  int write_given_output;
  double cd_tot_norm;
  double Ta; // atmospheric temperature in K, from NRLSMSIS if the user didn't chose "density_file" or "gitm", in which case atmo_temperature = 800K
  double density_mod; // the desnity given by msis is multiplied by density_mod + density_amp * sin(2*pi*t/T + density_phase*T) where T is the orbital period  
  double density_mod_amp;
  double density_mod_phase;

  int file_is_quaternion;
  double A_ref_tot;// see compute_drag (at the time I add this variable, I think this is the total cross section area (area with respect to relative velocity vector) but I haven't verified that yet (11-12-2017)
  double sum_cd_a_cos;//only used in kalman filter for now
  int last_compute_dxdt; // 1 if last time call compute_dxdt in propagate_spacecraft, 0 otherwise
  int opengl;
  int opengl_power; // 1 if name of file after solar cell efficiency value (4th line of main input file) -> will compute power from opengl_filename_solar_power. 0 otherwise -> will compute solar power from the geometry file, not the 3d model in opengl. if 0, then doesn't take into account the shadow effect (a surface shades another surface).
  int coll_vcm; // 1 is computing the probability of collision and the input file has the format of a VCM. 0 otherwise
    double bc_vcm; // if computing the probability of collision and the input file has the fomat of a VCM, this is the value of the ballistic coefficient reported in the VCM
  double srp_vcm; // if computing the probability of collision and the input file has the fomat of a VCM, this is the value of the solar radiation pressure coefficient reported in the VCM
  double et_vcm;// if computing the probability of collision and the input file has the format of a VCM, et_vcm is the epoch of the VCM
  double              a_i2cg_INRTL_ir_earth[3];
  int thrust; // 1 is section #THRUST exists in the main input file. 0 otherwise
  SpiceDouble elems[10];
}   INTEGRATOR_T;

typedef struct {
  
  double              density_here; // density at position of sc
  char                filename_kalman_init[1000];
  double              rho; // kalman
  double              et;                     // Seconds
  double              r_i2cg_INRTL[3];        // Postition of spacecraft in inertial coordinates
  double              v_i2cg_INRTL[3];        // Velocity of spacecraft in inertial coordinates
  double              a_i2cg_INRTL[3];        // Acceleration of spacecraft in inertial coordinates

  double              a_i2cg_LVLH[3]; // kalman
  double              a_i2cg_LVLH_drag[3];
  double              a_i2cg_LVLH_gravity[3];
  double              a_i2cg_INRTL_drag[3];
  double              a_i2cg_INRTL_gravity[3];

  double              r_ecef2cg_ECEF[3];      // Position of spacecraft in ECEF coordinates
  double              v_ecef2cg_ECEF[3];      // Velocity of spacecraft in ECEF coordinates
  double              number_of_collisions;
  ORBITAL_ELEMENTS_T  OE;
  GEODETIC_T          GEODETIC;
  INTEGRATOR_T        INTEGRATOR;
  int                 see_storm[N_STORM]; // one if the SC flies over the storm, 0 otherwise
  FILE                *fp;
  char                filename[1000];
  FILE                *fprho;
  char                filenamerho[1000];

  FILE                *fpout;
  char                filenameout[1000];
  FILE                *fpower;
  char                filenamepower[1000];
  FILE                *fpeclipse;
  char                filenameeclipse[1000];
  char                filenameatt[1000];
  FILE                *fpatt;
  FILE                *fpecef;
  char                filenameecef[1000];
  FILE                *fptle; // this contains the positions of the satellite from the TLE epoch to the constellation epoch starts
  char                filenametle[1000];
  char                name_sat[1000];
  FILE                *fpiproc[N_FILES_ENSEMBLES]; 
  char                filenameiproc[200][N_FILES_ENSEMBLES];
  FILE                *fpiproc_attitude[N_FILES_ENSEMBLES]; 
  char                filenameiproc_attitude[200][N_FILES_ENSEMBLES];
  FILE                *fp_coverage_ground_station[N_GROUND_STATION];
  char                filename_coverage_ground_station[N_GROUND_STATION][1000];
  int                 ispan; //used to save position of sc in the time span of the TCA if collisions assessment is on
  char                filenamekalman[1000]; // kalman
  char                filenamekalman_meas[1000]; // kalman measurement output file: the measruement output here is the same unit as the state. So it's the conversion of the measuremnt input file to the unit of the state. Ex: if the measurements are r/v in ECEF but the state r/v in ECI then this file is the measurements r/v in ECI. 
  double et_sc_initial; // initial epoch of each sc - can be different from a sc to another if TLE initialization or collision with VCM input file. Otherwize all equal to CONSTELLATION.et
  double a_i2cg_kalman[3];
  double et_next_time_step;// used in the kalman filter
  int already_output_cd_ensemble; // only used if collision with VCM. to output only the cd ensmebles onece ( not a t all time steps since it doesnt change with time if VCM collision_). 2 for 2 sc.
  int already_output_srp_ensemble; // only used if collision with VCM. to output only the srp ensmebles onece ( not a t all time steps since it doesnt change with time if VCM collision_). 2 for 2 sc.
  int orbit_number; // number of orbits travelled by the s/c since the start of the simu
  double et_last_orbit; // last time a full orbit was travelled (reference: start of simu)
} SPACECRAFT_T;

typedef struct {

  SPACECRAFT_T    **spacecraft;//[N_SATS][N_ENSEMBLES];
  double          et;
  FILE            *file_CYGNSS_constellation_for_specular;
  char            filename_CYGNSS_constellation_for_specular[1000];
  FILE            *file_GPS_constellation_for_specular;
  char            filename_GPS_constellation_for_specular[1000];
  char            filename_collision[1000];
  double          collision_time_span;
  int             *aaa_sigma;
  int             *aaa_mod;
  double **ensemble_array_per_iproc_f107_at_given_time;
  double **ensemble_array_per_iproc_f107_at_given_time_sorted;
  double **ensemble_array_per_iproc_ap_at_given_time;
  double **ensemble_array_per_iproc_ap_at_given_time_sorted;
  double **sum_sigma_for_f107_average;

  int nb_faces;
  int **nb_faces_theta_phi;
  double **normal_face;
  int ***which_face_theta_phi; // which face is seen at a given theta and phi
  double ***area_attitude_opengl; // cross section area computed by opengl as a function of the face, theta and ph i
  double **area_attitude_opengl_total; // total cross section area computed by opengl as a function theta and phi 

  double area_attitude_opengl_phi0;
  double area_attitude_opengl_theta0;
  double area_attitude_opengl_dtheta;
  double area_attitude_opengl_dphi;
  double **area_solar_panel_attitude_opengl;


} CONSTELLATION_T;


/* typedef struct { */

/*   double et_storm[N_STORM]; // epoch of the storm */
/*   double latitude[N_STORM]; // latitude of the storm */
/*   double longitude[N_STORM]; // longitude of the storm */
/*   double storm_ECEF[N_STORM][3]; // position of the storm in ECEF coordinates */
/*   double radius_uncertainty[N_STORM]; // radius of uncertainty of the storm */
/*   FILE *file_storm[N_STORM]; */
/*   FILE *file_storm_interpolated[N_STORM]; */
/*   char filename_storm[100][N_STORM]; */
/*   char filename_storm_interpolated[100][N_STORM]; */
/*   char storm_name[100][N_STORM]; */
/* } STORM_T; */

#endif

// Prototypes
int orbit_ave(  SPACECRAFT_T *SC,
			    CONSTELLATION_T *CONSTELLATION,
			    int iProc,
			    int iDebugLevel,
		double previous_an_to_sc
		    );

int eci2lla(double pos[3], double et,  double geodetic[3]  );
double  gstime (double jdut1 );
int jday (int year, int mon, int day, int hr, int minute, double sec, double *jd );
int coverage_ground_station(   SPACECRAFT_T *SC, // in: spacecraft (position and time). out: elevation, azimuth, range 
			       GROUND_STATION_T *GROUND_STATION, // in: list of ground stations
			       PARAMS_T *PARAMS, // in: parameters for the propagation   
			       int index_in_attitude_interpolated, // in: index of the current time step (take into account RK4)
			       INTEGRATOR_T    *INTEGRATOR, // in: paramters for the propagation (attitude of the sc here)
			       double          et_initial_epoch, // in: time of the inital epoch of the constellation
			       double          et_oldest_tle_epoch,
			       double sc_ecef_previous_time_step[3], // in: ECEF position of the spacecraft at the previous time step of the propagation (used to linear interpolate the position between the previous position and the current position)
			       double sc_eci_previous_time_step[3], //  in: ECI position of the spacecraft at the previous time step of the propagation (used to linear interpolate the position between the previous position and the current position)
			       double sc_eci_v_previous_time_step[3], //  in: ECI velocity of the spacecraft at the previous time step of the propagation (used to linear interpolate the position between the previous position and the current position)
			       double time_step_interpolation // in: time step of the linear interpolation of the sc position (in seconds)
			       );

int geodetic_to_geocentric(double flattening,            /* IN:     flattening parameter     */
			   double h,                     /* IN:  M  height above ellipsoid  */
			   double lat,                   /* IN:  r  geodetic latitude       */
			   double longitude,              /* IN:  r  longitude               */
			   double equatorial_radius,       /* IN: equatorial radius */
			   double  R_ecef_2_cg_ECEF[3]);   /* OUT:     vector in ECEF           */


int geocentric_to_geodetic(
			   double  R_pt_wrt_ecef_ECEF[3], /* vector in ECEF           */
			   double *semimajor_axis,        /* planetary radius         */
			   double *flattening,            /* flattening parameter     */
			   double *h,                     /* height above ellipsoid  */
			   double *lat,                   /* geodetic latitude       */
			   double *longitude   );         /* longitude               */


int cart2kep( ORBITAL_ELEMENTS_T *OE, double r[3], double v[3], double time, double u);

int kep2cart(   double              r_i2cg_INRTL[3],
                double              v_i2cg_INRTL[3],
                double             *mu,
                ORBITAL_ELEMENTS_T *OE);



/* int propagate_spacecraft(   SPACECRAFT_T *SC, */
/*                             PARAMS_T     *PARAMS, */
/* 			    double et_initial_epoch,  */
/* 			    double *density, */
/* 			    GROUND_STATION_T *GROUND_STATION); */
/* 			    //			       	    OPTIONS_T *OPTIONS); */
//newstructure
//int load_gravity(   GRAVITY_T *GRAVITY, char main_directory_location[256]);
int load_gravity(   GRAVITY_T *GRAVITY, char path_to_spice[256]);
//newstructure 

int compute_gravity(    double      a_i2cg_INRTL[3],
                        double      r_i2cg_INRTL[3],
                        double      et,
                        GRAVITY_T   *Gravity,
                        int         degree,
			char earth_fixed_frame[100],  double earth_flattening, double earth_radius, SPACECRAFT_T *SC, int gravity_map, CONSTELLATION_T *CONSTELLATION);
			//          int         order);

	
//newstructure 		
//int load_params( PARAMS_T *PARAMS, char main_directory_location[256], int iDebugLevel, char earth_fixed_frame[100] , double use_ap_hist, int iProc);
int load_params( PARAMS_T *PARAMS,  int iDebugLevel, char earth_fixed_frame[100] , double use_ap_hist, int iProc, char path_to_spice[256], int degree, int gravity_map_use, int earth_pressure);
//newstructure 

int compute_solar_pressure(double          a_solar_pressure_INRTL[3],
			   double          r_i2cg_INRTL[3],
			   double          v_i2cg_INRTL[3],
			   double          et,
			   PARAMS_T        *PARAMS,
			   INTEGRATOR_T    *INTEGRATOR,
			   CONSTELLATION_T *CONSTELLATION,
			   double          et_initial_epoch,
			       double          et_oldest_tle_epoch,
			   int             index_in_attitude_interpolated);

int compute_earth_pressure(double          a_solar_pressure_INRTL[3],
			   double          r_i2cg_INRTL[3],
			   double          v_i2cg_INRTL[3],
			   double          et,
			   PARAMS_T        *PARAMS,
			   INTEGRATOR_T    *INTEGRATOR,
			   CONSTELLATION_T *CONSTELLATION,
			   double          et_initial_epoch,
			       double          et_oldest_tle_epoch,
			   int             index_in_attitude_interpolated);


int shadow_light( char      shadow[256],
		  double    r_i2cg_INRTL[3],
		  double    et,
		  PARAMS_T  *PARAMS);

int shadow_light_moon( char      shadow[256],
		  double    r_i2cg_INRTL[3],
		  double    et,
		  PARAMS_T  *PARAMS);


int compute_power(INTEGRATOR_T    *INTEGRATOR,
		  CONSTELLATION_T *CONSTELLATION,
		  double          r_i2cg_INRTL[3],
		  double          v_i2cg_INRTL[3],
		  double          *et,
		  PARAMS_T        *PARAMS,
		  double          et_initial_epoch, 
			       double          et_oldest_tle_epoch,
		  int             index_in_attitude_interpolated);
		  
int print_test();
int print_teste();

int gitm_density( double *density, double et, double altitude, double latitude, double longitude, PARAMS_T *PARAMS);

int test_print( char to_print[1]  );
//int read_gitm_bin( double small_array[5], double ***longitude_gitm, double ***latitude_gitm, double ***altitude_gitm, double ***density_gitm, char time_gitm[256], int nLons, int nLats, int nAlts, char name_gitm_file_to_read[256]);
int test_print_iproc( int iProc , char to_print[1] );

int calculate_cd(double *cd,
		 double acco_coeff, 
		 double v_sc_mag, //in km/s
		 double Ta, // atmospheric temperature in K, from NRLSMSIS
		 double surface_area, // in m^2
		 double gamma,
		 double Ma // in kg/mol!!!!!from NRLMSIS?
		 );
int calculate_cd_opengl(double *cd,
		 double acco_coeff, 
		 double v_sc_mag, //in km/s
		 double Ta, // atmospheric temperature in K, from NRLSMSIS
		 double surface_area_projected, // in km^2
		 double gamma,
		 double Ma // in kg/mol!!!!!from NRLMSIS?
			);

double factorial(unsigned long f);
int build_gravity_map(GRAVITY_T  *Gravity, int degree,  int iProc);
int read_gravity_map(GRAVITY_T  *Gravity, int degree,  int iProc);
int build_earth_pressure_map(GRAVITY_T  *Gravity, int iProc);
int read_earth_pressure_map(GRAVITY_T  *Gravity, int iProc);

double compute_earth_albedo();
double compute_earth_emissivity();
int polynomial_interpo(double *val, int  n, double *x, double *y, double to_inter);
int compute_iradius_gravity_map(int iradius_arr[2], GRAVITY_T *Gravity, double rmag);
int compute_ilat_gravity_map(int ilat_arr[4], GRAVITY_T *Gravity, double lat_gc);
int compute_ilon_gravity_map(int ilon_arr[4], GRAVITY_T *Gravity, double lon_gc_corr);
int gravity_map_yinter_lat(double *yinter, int *order_interpo_map, double lat_gc, GRAVITY_T *Gravity, double y_lat0, double y_lat1, double y_lat2, double y_lat3);
int gravity_map_yinter_radius(double *yinter, int *order_interpo_map, double rmag, GRAVITY_T *Gravity, double y_radius0, double y_radius1, double y_radius2, double y_radius3);
int gravity_map_xinter(double *xinter_lon, double *xinter_lat, double *xinter_radius,  double long_gc_corr, double lat_gc, double rmag, GRAVITY_T *Gravity, int *ilon_arr, int *ilat_arr, int *iradius_arr);
int gravity_map_yinter_lon(double *yinter, int *order_interpo_map, double long_gc_corr, GRAVITY_T *Gravity, double y_lon0, double y_lon1, double y_lon2, double y_lon3);



int compute_iradius_earth_pressure_map(double *xinter_radius, int iradius_arr[2], int *nradius_interpo, GRAVITY_T *Gravity, double rmag);
int compute_izenith_earth_pressure_map(double *xinter_zenith, int izenith_arr[2], int *nzenith_interpo, GRAVITY_T *Gravity, double zenith);
int compute_ielev_surf_earth_pressure_map(double *xinter_elev_surf, int ielev_surf_arr[2], int *nelev_interpo, GRAVITY_T *Gravity, double elev_surf);
int compute_iazim_surf_earth_pressure_map(double *xinter_azim_surf, int iazim_surf_arr[2], int *nazim_interpo, GRAVITY_T *Gravity, double azim_surf_corr);

int earth_pressure_map_yinter_azim_surf(double **yinter, int *order_interpo_map, double azim_surf_corr, GRAVITY_T *Gravity, double **y_az, int iazim_surf_arr[2]);
int earth_pressure_map_yinter_elev_surf(double **yinter, int *order_interpo_map, double elev_surf, GRAVITY_T *Gravity, double **y_el);
int earth_pressure_map_yinter_zenith(double **yinter, int *order_interpo_map, double zenith, GRAVITY_T *Gravity, double **y_zen);
int gravity_map_yinter_radius_earth(double **yinter, int *order_interpo_map, double rmag, GRAVITY_T *Gravity, double **y_rad);
int polynomial_interpo_earth(double *val, int n, double *x, double **y, double to_inter);
