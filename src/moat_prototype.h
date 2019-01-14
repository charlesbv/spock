#include <mpi.h>
#include "propagator.h"
#include "options.h"

int initialize_constellation( CONSTELLATION_T  *CONSTELLATION,
                              OPTIONS_T        *OPTIONS,
                              PARAMS_T         *PARAMS,
			      GROUND_STATION_T *GROUND_STATION,
			      int              iDebugLevel,
			      int iProc,
			      int nProcs);

int generate_ephemerides(   CONSTELLATION_T  *CONSTELLATION,
                            OPTIONS_T        *OPTIONS,
                            PARAMS_T         *PARAMS, 
			    GROUND_STATION_T *GROUND_STATION,
			    int              iProc,
			    int              nProcs,
			    int              iDebugLevel);

int write_output(   SPACECRAFT_T    *SC,
                    int             init_flag,
                    int             choose_tle_to_initialise_orbit,
                    int             sc_index,
                    int             n_sc,
		    int             n_gps,
		    int             nb_ensembles_min,
		    //                    int             nb_ensembles_attitude,
		    int             nb_ensemble_min_per_proc,
		    //                     int             nb_ensemble_attitude_per_proc,
		    int             iProc,
		    int             nb_ensembles_output_files,
		    char            name_file_ensembles[30][1000],// !!! the number of columns has to correpond to n_files_ensembles
		    //		    char            name_file_ensembles_attitude[1000][1], // !!! the number of columns has to correpond to n_files_ensembles_attitude
		    double previous_lat,
		    OPTIONS_T       *OPTIONS,
		    char earth_fixed_frame[100],
		    int              write_reference_sc,
		    int              write_ensembles,
		    double constellation_et,
		    int nProcs,
		    int iDebugLevel,
		    int compute_collisions,
		    int *start_ensemble,
		    int *array_sc,
		    CONSTELLATION_T *CONSTELLATION,
		    PARAMS_T *PARAMS);


int interpolate_position(   char filenameout[1000],
			    char filenameecef[1000],
			    //		    char filename[1000],
                            double          dt_before,
                            double          dt_after, 
			    int             sat_number,
			    OPTIONS_T *OPTIONS,
			    PARAMS_T *PARAMS);

int ancas(double min_dist_close_approach, double et_end, double dt_interval, double r1_start[3], double v1_start[3], double a1_start[3], double r1_end[3],double  v1_end[3], double  a1_end[3],double r2_start[3], double v2_start[3], double a2_start[3], double r2_end[3],double  v2_end[3], double  a2_end[3], double gamma0, double gamma1, double gamma2, double gamma3, double *tca1, double *dca1, double *tca2, double *dca2, double *tca3, double *dca3);

int ancas_existence_of_min_using_two_values_of_function_and_two_values_of_its_derivative( double et_start, double dt_interval, double r1_start[3], double v1_start[3], double a1_start[3], double r1_end[3],double  v1_end[3], double  a1_end[3],double r2_start[3], double v2_start[3], double a2_start[3], double r2_end[3],double  v2_end[3], double  a2_end[3], int *min_exists,  double *gamma0, double *gamma1, double *gamma2, double *gamma3);

int close_approach_ensemble( double *eci_x_primary_sc_in_span, double *eci_y_primary_sc_in_span, double *eci_z_primary_sc_in_span, double  *eci_vx_primary_sc_in_span, double *eci_vy_primary_sc_in_span, double *eci_vz_primary_sc_in_span, double *eci_ax_primary_sc_in_span, double *eci_ay_primary_sc_in_span, double *eci_az_primary_sc_in_span,  
double *eci_x_secondary_sc_in_span, double *eci_y_secondary_sc_in_span, double *eci_z_secondary_sc_in_span, double  *eci_vx_secondary_sc_in_span, double *eci_vy_secondary_sc_in_span, double *eci_vz_secondary_sc_in_span, double *eci_ax_secondary_sc_in_span, double *eci_ay_secondary_sc_in_span, double *eci_az_secondary_sc_in_span,  
			     int time_step_of_tca,double *gamma0, double *gamma1, double *gamma2, double *gamma3,  int *min_exists, double r1_start[3], double v1_start[3], double a1_start[3], double r1_end[3],double  v1_end[3], double  a1_end[3],double r2_start[3], double v2_start[3], double a2_start[3], double r2_end[3],double  v2_end[3], double a2_end[3], int *time_step_start_interval, double *min_distance_in_time_spanning_tca, double *direction_distance, int initial_epoch_time_step_in_span, int final_epoch_time_step_in_span,  int iProc, double et_time_step_of_save_tca, OPTIONS_T *OPTIONS);

int print_progress(double min_end_time, double et , double starttime, int iProc, int nb_gps);

int print_progress_epoch_sc_to_epoch_constellation(double min_end_time, double et , double starttime, int iProc, int nb_gps);
//int send_r_v_a_in_tca_span( double ****save_r_i2cg_INRTL, double ****save_v_i2cg_INRTL, double ****save_a_i2cg_INRTL, int iProc, int nProcs, int ii, int eee, int nb_ensemble_min_per_proc, CONSTELLATION_T *CONSTELLATION);

int receive_r_v_a_in_tca_span( double ****save_r_i2cg_INRTL, double ****save_v_i2cg_INRTL, double ****save_a_i2cg_INRTL, int iProc, int nProcs, int ii, int eee, int nb_ensemble_min_per_proc, CONSTELLATION_T *CONSTELLATION, int nProcs_that_are_gonna_run_ensembles);
/* int interpolate_position(   SPACECRAFT_T    *SC, */
/*                             double          dt_before, */
/*                             double          dt_after,  */
/* 			    OPTIONS_T       *OPTIONS, */
/* 			    int             sat_number, */
/* 			    PARAMS_T        *PARAMS); */
double distance_between_two_sc( double eci_r_sc1[3], double eci_r_sc2[3] );

int allocate_memory_r_a_v_in_span(double ****save_x_i2cg_INRTL, double ****save_y_i2cg_INRTL, double ****save_z_i2cg_INRTL, double ****save_vx_i2cg_INRTL, double ****save_vy_i2cg_INRTL, double ****save_vz_i2cg_INRTL, double ****save_ax_i2cg_INRTL, double ****save_ay_i2cg_INRTL, double ****save_az_i2cg_INRTL, int nb_tca, int nb_satellites_not_including_gps, int total_ensemble_final_with_ref, int nb_time_steps_in_tca_time_span, int iProc, int iDebugLevel);

int compute_collision_between_one_secondary_and_all_primary(double *save_x_i2cg_INRTL_sec, double *save_y_i2cg_INRTL_sec, double *save_z_i2cg_INRTL_sec, double *save_vx_i2cg_INRTL_sec, double *save_vy_i2cg_INRTL_sec, double *save_vz_i2cg_INRTL_sec,double *save_ax_i2cg_INRTL_sec, double *save_ay_i2cg_INRTL_sec, double *save_az_i2cg_INRTL_sec, double *save_x_i2cg_INRTL_prim, double *save_y_i2cg_INRTL_prim, double *save_z_i2cg_INRTL_prim, double *save_vx_i2cg_INRTL_prim, double *save_vy_i2cg_INRTL_prim, double *save_vz_i2cg_INRTL_prim,double *save_ax_i2cg_INRTL_prim, double *save_ay_i2cg_INRTL_prim, double *save_az_i2cg_INRTL_prim,OPTIONS_T *OPTIONS, int iProc,  int *nb_coll_per_step_per_iproc_in_tca, double *et_time_step_of_save_tca, int nb_time_steps_in_tca_time_span, int iiitca, int eee_prim, int eee_sec, FILE *tca_file, FILE *dca_file,FILE *sample_file, int write_collision_files, int *eee_prim_that_collide);

int print_progress_collision(int eee_sec, int iProc, int nb_ensemble_min_per_proc, int nb_tca);


int compute_heading(double *heading, double v_ecef[3], double lon, double lat, double earth_flattening);



