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
#define N_SATS 60
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
  double mass;                   // Mass of spacecraft
  double solar_cell_efficiency;  // Efficiency of the solar cells
  SURFACE_T surface[N_SURFACES]; // Properties of the surface
  double nb_surfaces;            // number of surfaces
  ATTITUDE_T attitude;           // Attitude of the satellite
  double degree;                 // Gravity degree
  double order;                  // Gravity order
  int include_drag;              // include drag
  int include_moon;              // include moon perturbations
  int include_sun;               // include sun perturbations
  int include_solar_pressure;    // include solar pressure
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
  FILE                *file_given_output; // this is a file used to output things in particular. Not used unless hard coded
  char                filename_given_output[1000]; // this is a file used to output things in particular. Not used unless hard coded
  int sc_main_nb;
  int sc_ensemble_nb;
  int write_given_output;
  double cd_tot_norm;
  double Ta; // atmospheric temperature in K, from NRLSMSIS if the user didn't chose "density_file" or "gitm", in which case atmo_temperature = 800K
  double density_mod; //factor to apply on density at position of satellite (calculated by NRLMSIS, gitm or from density file)
  int file_is_quaternion;
  double A_ref_tot;// see compute_drag (at the time I add this variable, I think this is the total cross section area (area with respect to relative velocity vector) but I haven't verified that yet (11-12-2017)
  double sum_cd_a_cos;//only used in kalman filter for now
  int last_compute_dxdt; // 1 if last time call compute_dxdt in propagate_spacecraft, 0 otherwise
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
  double et_sc_initial; // initial epoch of each sc - can be different from a sc to another if TLE initialization. Otherwize all equal to CONSTELLATION.et
  double a_i2cg_kalman[3];
  double et_next_time_step;// used in the kalman filter

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
			char earth_fixed_frame[100],  double earth_flattening, double earth_radius, SPACECRAFT_T *SC);
			//          int         order);

	
//newstructure 		
//int load_params( PARAMS_T *PARAMS, char main_directory_location[256], int iDebugLevel, char earth_fixed_frame[100] , double use_ap_hist, int iProc);
int load_params( PARAMS_T *PARAMS,  int iDebugLevel, char earth_fixed_frame[100] , double use_ap_hist, int iProc, char path_to_spice[256]);
//newstructure 

int compute_solar_pressure(double          a_solar_pressure_INRTL[3],
			   double          r_i2cg_INRTL[3],
			   double          v_i2cg_INRTL[3],
			   double          et,
			   PARAMS_T        *PARAMS,
			   INTEGRATOR_T    *INTEGRATOR,
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
		  double          r_i2cg_INRTL[3],
		  double          v_i2cg_INRTL[3],
		  double          *et,
		  PARAMS_T        *PARAMS,
		  double          et_initial_epoch, 
			       double          et_oldest_tle_epoch,
		  int             index_in_attitude_interpolated);
		  
int print_test();

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
