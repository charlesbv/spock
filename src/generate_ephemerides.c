// Licensed to the Apache Software Foundation (ASF) under one
// or more contributor license agreements.  See the NOTICE file
// distributed with this work for additional information
// regarding copyright ownership.  The ASF licenses this file
// to you under the Apache License, Version 2.0 (the
// "License"); you may not use this file except in compliance
// with the License.  You may obtain a copy of the License at

//   http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.

#include "options.h"
#include "moat_prototype.h"
#include "gsl/gsl_poly.h"

/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           generate_ephemerides
//  Purpose:        Generates ephemerides for constellation
//  Assumptions:    None
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | J. Getchius   | 05/20/2015    |   ---     | Initial Implementation
//      | C. Bussy-Virat| 08/04/2015    |   ---     | The end time is read from the input file
//
/////////////////////////////////////////////////////////////////////////////////////////
int generate_ephemerides(   CONSTELLATION_T  *CONSTELLATION,
                            OPTIONS_T        *OPTIONS,
                            PARAMS_T         *PARAMS, 
			    GROUND_STATION_T *GROUND_STATION,
			    int iProc, 
			    int nProcs,
			    int iDebugLevel)
{
 
 
  int eee_prim_that_collide = -1;
  // if set to 0 then don't write in the given files. if set to 1 the write in given files
  int write_density = 0;

  int output_only_at_tca = 1; // set output_only_at_tca to 1 to output only at tca  (if computing collisions)
  if ( ( iDebugLevel >= 1 ) ){
    printf("-- (generate_ephemerides) Just got in generate_ephemerides. (iProcc %d)\n", iProc);
  }


  /* print_oe(CONSTELLATION->spacecraft[0][0].OE, PARAMS); */
  /* print_oe(CONSTELLATION->spacecraft[1][0].OE, PARAMS); */

  /* // ECI */
  /* 	  /\* printf("(%.10f; %.10f; %.10f) (%.10f; %.10f; %.10f)\n", CONSTELLATION->spacecraft[0][0].r_i2cg_INRTL[0]*1000., CONSTELLATION->spacecraft[0][0].r_i2cg_INRTL[1]*1000., CONSTELLATION->spacecraft[0][0].r_i2cg_INRTL[2]*1000., CONSTELLATION->spacecraft[0][0].v_i2cg_INRTL[0]*1000., CONSTELLATION->spacecraft[0][0].v_i2cg_INRTL[1]*1000., CONSTELLATION->spacecraft[0][0].v_i2cg_INRTL[2]*1000.); *\/ */
  /* 	  /\* printf("(%.10f; %.10f; %.10f) (%.10f; %.10f; %.10f)\n", CONSTELLATION->spacecraft[1][0].r_i2cg_INRTL[0]*1000., CONSTELLATION->spacecraft[1][0].r_i2cg_INRTL[1]*1000., CONSTELLATION->spacecraft[1][0].r_i2cg_INRTL[2]*1000., CONSTELLATION->spacecraft[1][0].v_i2cg_INRTL[0]*1000., CONSTELLATION->spacecraft[1][0].v_i2cg_INRTL[1]*1000., CONSTELLATION->spacecraft[1][0].v_i2cg_INRTL[2]*1000.); *\/ */
  /* 	  // ECEF / */
  /* double fac_ecef = 1.; */
  /* 	  printf("(%.10f; %.10f; %.10f) (%.10f; %.10f; %.10f)\n", CONSTELLATION->spacecraft[0][0].r_ecef2cg_ECEF[0]*fac_ecef, CONSTELLATION->spacecraft[0][0].r_ecef2cg_ECEF[1]*fac_ecef, CONSTELLATION->spacecraft[0][0].r_ecef2cg_ECEF[2]*fac_ecef, CONSTELLATION->spacecraft[0][0].v_ecef2cg_ECEF[0]*fac_ecef, CONSTELLATION->spacecraft[0][0].v_ecef2cg_ECEF[1]*fac_ecef, CONSTELLATION->spacecraft[0][0].v_ecef2cg_ECEF[2]*fac_ecef); */
  /* 	  printf("(%.10f; %.10f; %.10f) (%.10f; %.10f; %.10f)\n", CONSTELLATION->spacecraft[1][0].r_ecef2cg_ECEF[0]*fac_ecef, CONSTELLATION->spacecraft[1][0].r_ecef2cg_ECEF[1]*fac_ecef, CONSTELLATION->spacecraft[1][0].r_ecef2cg_ECEF[2]*fac_ecef, CONSTELLATION->spacecraft[1][0].v_ecef2cg_ECEF[0]*fac_ecef, CONSTELLATION->spacecraft[1][0].v_ecef2cg_ECEF[1]*fac_ecef, CONSTELLATION->spacecraft[1][0].v_ecef2cg_ECEF[2]*fac_ecef); */

  /* exitall(); */
  /* double v_angle[3]; */
  /* v_angle[0] = 0;   v_angle[1] = 0;  v_angle[2] = 0; */

  // ensembles on attitude
  /* int n_files_ensembles_attitude = 1; */
  /* char name_file_ensembles_attitude[256][1]; // !!! the number of columns has to correpond to n_files_ensembles_attitude. CHANGE ALSO in declaration of write_output (in .c and .h files) */
  /* strcpy(name_file_ensembles_attitude[0], "power"); */
  //  char temp_iproc_file_attitude[256], temp_nb_proc_attitude[256];


  /* int n_files_ensembles = 6; */
  /* char name_file_ensembles[256][6]; // !!! the number of columns has to correpond to n_files_ensembles. CHANGE ALSO in declaration of write_output (in .c and .h files) */
  /* strcpy(name_file_ensembles[0], "x"); */
  /* strcpy(name_file_ensembles[1], "y"); */
  /* strcpy(name_file_ensembles[2], "z"); */
  /* strcpy(name_file_ensembles[3], "vx"); */
  /* strcpy(name_file_ensembles[4], "vy"); */
  /* strcpy(name_file_ensembles[5], "vz"); */

  /* Declaration */
  //!!!!!!!!!! COMMENT THIS BLOCK !!!!!!!!!!
    /* char remove_this_var_temp[256]; */
    /* strcpy(remove_this_var_temp, "2016-11-29T00:00:00"); */
    /* double remove_this_var; */
    /* str2et_c(remove_this_var_temp, &remove_this_var); */
    //!!!!!!!!!! END OF COMMENT THIS BLOCK !!!!!!!!!!

  char temp_iproc_file[256], temp_nb_proc[256];
  double et_current_tca;
    int start_itca=0;
	char time_itca2[256];
      double **save_first_distance_unpertubed_orbit = NULL;
      double **save_last_distance_unpertubed_orbit = NULL;
    int nb_tca_without_collisions;
    int nb_tca_with_collisions;



  int compute_collisions_was_on_but_no_tca_found = 0;
      double *total_collision_each_dt;
  int already_propagated_ref_sc = 0;
    FILE *file_collision = NULL;
	double et_step_collision;

	//  int eee_sec;
  int eee_prim;

  int ispan_constellation;
      int eee_all_other_iproc;
  int total_ensemble_final = OPTIONS->nb_ensemble_min_per_proc*nProcs;
  int total_ensemble_final_with_ref = total_ensemble_final + 1;
  double et_start_of_span, et_end_of_span;
    int **total_nb_collisions = NULL;
      char time_itca[256];
  int ppp;
  int ***nb_coll_per_step_in_TCA = NULL;
  int iiitca;
  int ****nb_coll_per_step_per_iproc_in_tca = NULL;

  //  int istepspan;
  float min_altitude_constellation = 1.0e32;
  double **previous_eci_r = NULL;
  double **previous_eci_v = NULL;
  double **previous_eci_a = NULL;
   int isc_ref;
   // int  isc_ens;

  //  int time_step_start_interval;
  //  double r_primary_start[3], v_primary_start[3], a_primary_start[3], r_primary_end[3], v_primary_end[3],  a_primary_end[3],r_secondary_start[3], v_secondary_start[3], a_secondary_start[3], r_secondary_end[3], v_secondary_end[3], a_secondary_end[3];
  int close_approach_exists;
  double gamma0,  gamma1, gamma2, gamma3;
  int time_step_of_tca;
  int *done_with_tca = NULL;
  int itca;
  //  int eee_primary_sc;
  int iground;
/*   FILE *density_file; */
/*   char density_file_name[500]; */
/*   strcpy(density_file_name, OPTIONS->dir_output_run_name_sat_name[0]); */
/*   strcat(density_file_name, "/density_"); */
/*   strcat(density_file_name, OPTIONS->filename_output[0]); */
/*   density_file = fopen(density_file_name, "w+"); */



  FILE *tca_file;
  FILE *dca_file;
  FILE *sample_file;

  // ensembles on COE
  /* FILE *fp_temp_att1; */
  /*   fp_temp_att1 = fopen("attitude_132_1.a", "w+"); */
  /* /\* FILE *fp_temp_att2; *\/ */
  /* /\*   fp_temp_att2 = fopen("attitude_646.a", "w+"); *\/ */

  /* FILE *fp_temp_att_us; */
  /* fp_temp_att_us = fopen("attitude_us_132.a", "w+"); */
  //	      float frac1, frac2;
  double density;
  char times_att[300]; 
  double previous_lat = 0;
  double tca1 = -1, dca1 = 1e6, tca2 = -1, dca2 = 1e6, tca3 = -1, dca3 = 1e6;
  //  double tca_ensemble_1 = -1, dca_ensemble_1 = 1e6, tca_ensemble_2 = -1, dca_ensemble_2 = 1e6, tca_ensemble_3 = -1, dca_ensemble_3 = 1e6;
  int nb_tca = 0;
  double *save_tca = NULL; // !!!!!!!!!! FOR NOW WORKS ONLY IF TWO REFERENCE SATELLIES ONLY  (add a dimension if more than 2 sc)
  double *et_time_step_of_save_tca = NULL; // !!!!!!!!!! FOR NOW WORKS ONLY IF TWO REFERENCE SATELLIES ONLY  (add a dimension if more than 2 sc)
  double *save_dca; // !!!!!!!!!! FOR NOW WORKS ONLY IF TWO REFERENCE SATELLIES ONLY  (add a dimension if more than 2 sc)
  double  max_tca = -1;
  double dca_at_max_tca = 1.e6;
  double ***save_x_i2cg_INRTL = NULL, ***save_y_i2cg_INRTL = NULL, ***save_z_i2cg_INRTL = NULL;
  double ***save_vx_i2cg_INRTL = NULL, ***save_vy_i2cg_INRTL = NULL, ***save_vz_i2cg_INRTL = NULL;
  double ***save_ax_i2cg_INRTL = NULL, ***save_ay_i2cg_INRTL = NULL, ***save_az_i2cg_INRTL = NULL;
  int nb_time_steps_in_tca_time_span = (int) (nearbyint( CONSTELLATION->collision_time_span / OPTIONS->dt )) + 1 ; // the cast here is not necessary as we made sure in initialize_constellation that CONSTELLATION->collision_time_span is an even multiple of OPTIONS->dt
  //pti(nb_time_steps_in_tca_time_span, "nb_time_steps_in_tca_time_span");exitall();
  double min_end_time;
  // other variables
  int ii, eee, fff;
  double starttime, endtime;
  double twrite = 0.0;
  double time_between_last_gps_epoch_and_constellation_epoch_starttime, new_dt_for_gps;
  int save_include_drag, save_include_solar_pressure, save_include_earth_pressure;
  double save_solar_cell_efficiency;
  int choose_tle_to_initialise_orbit = 0;
  int ccc;
  int compute_collisions = 0;
  //  int start_ensemble;
  // COLLISIONS
  if ( ( OPTIONS->nb_satellites_not_including_gps > 1 ) && ( OPTIONS->nb_ensembles > 0 ) && ( (strcmp(OPTIONS->type_orbit_initialisation, "collision" ) == 0 ) || (strcmp(OPTIONS->type_orbit_initialisation, "collision_vcm" ) == 0 ))){
    compute_collisions = 1;
  }
  // end of COLLISIONS

  /*************************************************************************************/
  /*************************************************************************************/
  /*********************** SET THINGS UP FOR PARALLEL PROGRAMMING **********************/
  /*************************************************************************************/
  /*************************************************************************************/
/*   Things we set up here: */
/*   - which iProc runs which main sc / which sc is run by which iProc */

  int nProcs_that_are_gonna_run_ensembles;
  nProcs_that_are_gonna_run_ensembles = nProcs;
  if ( nProcs > OPTIONS->nb_ensembles_min ){
    nProcs_that_are_gonna_run_ensembles = OPTIONS->nb_ensembles_min;
  }


  // For each iProc, set up the first main sc (iStart_save[iProcf]) and the last main sc (iEnd_save[iProcf]) that it's going to run. iStart_save and iEnd_save are two arrays that have the same values for all procs (so if you're iProc 0 or iProc 1, you have value recorded for iStart_save[0] and iEnd_save[0] and the same value recorded for iStart_save[1] and iEnd_save[1]) -> they are not "iProc-dependent"
    int *iStart_save, *iEnd_save;
    int nscEachPe, nscLeft;
    int iProcf;
    int *iProc_run_main_sc;
    int i;
    iProc_run_main_sc = malloc(nProcs * sizeof(int));
    
    nscEachPe = (OPTIONS->n_satellites)/nProcs;
    nscLeft = (OPTIONS->n_satellites) - (nscEachPe * nProcs);

    iStart_save = malloc( nProcs * sizeof(int));
    iEnd_save = malloc( nProcs  * sizeof(int));
    for (iProcf = 0; iProcf < nProcs; iProcf++){
      iStart_save[iProcf] = 0;
      iEnd_save[iProcf] = 0;
    }
    for (iProcf = 0; iProcf < nProcs; iProcf++){
      for (i=0; i<iProcf; i++) {
	iStart_save[iProcf] += nscEachPe;
	if (i < nscLeft && iProcf > 0) iStart_save[iProcf]++;
      }
      iEnd_save[iProcf] = iStart_save[iProcf]+nscEachPe;
      if (iProcf  < nscLeft) iEnd_save[iProcf]++;
      iStart_save[iProcf] = iStart_save[iProcf];
      iEnd_save[iProcf] = iEnd_save[iProcf];
      if ( iStart_save[iProcf] >= iEnd_save[iProcf] ){ // this happens if the iProc iProcf does not run any main sc
	iProc_run_main_sc[iProcf] = 0;
      }
      else{
	iProc_run_main_sc[iProcf] = 1;
      }
/*       if ( iProc == 0){ */
/*       printf("iProc_run_main_sc[%d] = %d\n", iProcf, iProc_run_main_sc[iProcf]); */
/*       } */

    }
    
     /*   if (iProc == 0){ */
/*     for (iProcf = 0; iProcf < nProcs; iProcf++){ */
/*       printf("%d - %d\n", iStart_save[iProcf], iEnd_save[iProcf] - 1) ; */
/*     } */
/*     } */

    // For each main sc, start_ensemble is 0 if the iProc runs this main sc. start_ensemble is 1 is the iProc does not run this main sc. (so each iProc has a different array start_ensemble -> start_ensemble is "iProc-dependent") 
  int *start_ensemble;
  start_ensemble = malloc(OPTIONS->n_satellites * sizeof(int));
  for (ii = 0; ii < OPTIONS->n_satellites; ii++){
    if ( (ii >= iStart_save[iProc]) & ( ii < iEnd_save[iProc]) ){
      start_ensemble[ii] = 0;
    }
    else{
      start_ensemble[ii] = 1;
    }

    //        printf("iProc %d | start_ensemble[%d] = %d\n", iProc, ii, start_ensemble[ii]);

  }

  // array_sc is the array of sc (main and ensemble sc) run by this iProc. So array_sc is "iProc-dependent": each iProc has a different array array_sc. What array_sc has: the first elt is 0, whatever iProc it is (because it represents main sc. However, it does not mean that all iProc will run a main sc, this is decided later in the code (using start_ensemble)). The next elts (1, 2, ..., OPTIONS->nb_ensemble_min_per_proc) represent the ensemble sc run by this iProc. So for example if there are 20 ensembles and 2 iProc, array_sc is:
  // for iProc 0: [0, 1, 2, 3, ..., 9, 10]
  // for iProc 1: [0, 11, 12, 13, ..., 19, 20]
  int ielt;
  int *array_sc;
  array_sc = malloc((OPTIONS->nb_ensemble_min_per_proc + 2) * sizeof(int)); // + 2 (and not + 1) because if there is no ensemble, we still want array_sc[1] to exist (because later in the code we call array_sc[start_ensemble[ii]], and start_ensemble[ii] = 1 if the iProc does not run main sc ii). 

  array_sc[0] = 0;
  array_sc[1] = -1; // if there is no ensemble, we still want array_sc[1] to exist (because later in the code we call array_sc[start_ensemble[ii]], and start_ensemble[ii] = 1 if the iProc does not run main sc ii). We set to -1 because we want to make sure that if there is no ensemble, eee will never be equal to array_sc[1]. If there are ensembles, array_sc[1] is overwritten right below anyway

  for ( ielt = 1; ielt < OPTIONS->nb_ensemble_min_per_proc + 1; ielt ++ ){
    if (iProc < nProcs_that_are_gonna_run_ensembles){
    array_sc[ielt] = iProc * OPTIONS->nb_ensemble_min_per_proc + ielt;
    }
  }
/*   nscEachPe = (OPTIONS->n_satellites)/nProcs; */
/*   nscLeft = (OPTIONS->n_satellites) - (nscEachPe * nProcs); */

/*   if (iProc == 0){ */
/*     printf("XXXXXXXXX\nnscEachPe = %d | nscLeft = %d\nXXXXXXXXX\n", nscEachPe, nscLeft); */
/*   } */

  // for each main, which_iproc_is_running_main_sc is the iProc that runs it. For example, which_iproc_is_running_main_sc[3] is equal to the iProc that runs the main sc 3. So which_iproc_is_running_main_sc is the same array for all iProc -> which_iproc_is_running_main_sc is not "iProc-dependent"
  int *which_iproc_is_running_main_sc;
  which_iproc_is_running_main_sc = malloc( OPTIONS->n_satellites * sizeof(int));
  for (ii = 0; ii < OPTIONS->n_satellites; ii++){
    for (ccc = 0; ccc < nProcs; ccc++){
      if ( ( ii >= iStart_save[ccc] ) && ( ii < iEnd_save[ccc]  ) ){
	which_iproc_is_running_main_sc[ii] = ccc;
      }  
    }
    
  }

/*   for (ii = 0; ii < OPTIONS->n_satellites; ii++){ */
/*     if (iProc == 0){ */
/*     printf("iProc %d | which_iproc_is_running_main_sc[%d] = %d\n", iProc, ii, which_iproc_is_running_main_sc[ii]); */
/*     } */
/*   } */

  /*************************************************************************************/
  /*************************************************************************************/
  /****************** end of SET THINGS UP FOR PARALLEL PROGRAMMING ********************/
  /*************************************************************************************/
  /*************************************************************************************/

/*     if (iProc == 0){ */
/*       start_ensemble = 0; */
/*   } */
/*     else{ */
/*       start_ensemble = 1; */
/*     } */


  // Compute the starttime and endtime
  str2et_c(OPTIONS->initial_epoch, &starttime);
  str2et_c(OPTIONS->final_epoch, &endtime);
  min_end_time = endtime;

  //  endtime = endtime - OPTIONS->dt; // just remove the last time step
  // endtime = CONSTELLATION->et + 86400.0 * OPTIONS->ddays; // previous version (v2)
  twrite = CONSTELLATION->et;
  //  et2utc_c(CONSTELLATION->spacecraft[0][0].et, "ISOC" ,0 ,255 , times_att);



    // Check the TLE epochs of the satellites are older than the simulation epoch (if the user chose to initilaize the orbits with TLEs)
    // // If running TLE (GPS or other satellites), the epochs of the last TLEs have be older than the epoch of the start time of the constellation (the current version does not allow for propagating satellite going back in time)
  if ( (strcmp( OPTIONS->type_orbit_initialisation, "tle" ) == 0 ) || (strcmp(OPTIONS->type_orbit_initialisation, "tle_sgp4" ) == 0 ) ){
      for (ii = 0; ii < OPTIONS->n_satellites - OPTIONS->nb_gps; ii++){
	if ( start_ensemble[ii] == 0 ){ // if main sc ii is run by this iProc
	if (CONSTELLATION->spacecraft[ii][0].et > CONSTELLATION->et){
	  printf("The epochs of the TLEs of the satellites have to be older than the epoch of the start time of the constellation. The program will stop. \n");
	  printf("\nYou can choose 'now n' as the first line of the #TIME section in the input file input.d to run the constellation as from the current time. This will guarantee that the TLEs of the satellites are older than the epoch of the start time of the constellation ('n' is the number of hours to propagate the constellation for (n can be a decimal value)). \n");
	  exit(0);
	}
      }
      } 
      choose_tle_to_initialise_orbit = 1; // this means that the user chose to initialize the orbits with TLEs
    }


    // Check the TLE epochs of the GPS are older than the simulation epoch (if the user chose to run GPS)
    if (OPTIONS->nb_gps > 0){
      for (ii = OPTIONS->n_satellites - OPTIONS->nb_gps; ii < OPTIONS->n_satellites; ii++){
        if ( start_ensemble[ii] == 0 ){ // if main sc ii is run by this iProc                                        

	if (CONSTELLATION->spacecraft[ii][0].et > CONSTELLATION->et){
	  printf("The epochs of the last TLEs of the GPS satellites have to be older than the epoch of the start time of the constellation. The program will stop. \n");
	  printf("\nYou can choose 'now n' as the first line of the #TIME section in the input file input.d to run the constellation as from the current time. This will guarantee that the last TLEs of the GPS satellites are older than the epoch of the start time of the constellation ('n' is the number of hours to propagate the constellation for (n can be a decimal value)). \n");
	  exit(0);
	}
	}
      }
    }


    if ( ( iDebugLevel >= 2 ) ){
      printf("-- (generate_ephemerides) Opening output files... (iProc %d)\n", iProc);
    }

    // Open SpOCK's output files (other than ensemble output files)
    for (ii = 0; ii < OPTIONS->n_satellites; ii++) { // go through all main sc


	// Now CONSTELLATION->spacecraft[ii][0].filenameecef is moved to generate_ephemerides because we want iProc 0 to know CONSTELLATION->spacecraft[ii][0].filenameecef even for the main sc ii that it does not run (this is because iProc 0 will gather all ECEF files at the end of the propagation)
/*       if (iProc == 1){ */
/*       printf("<%d %d %d>\n", ii , OPTIONS->n_satellites, iProc); */
/*       } */
		strcpy(CONSTELLATION->spacecraft[ii][0].filenameecef, OPTIONS->dir_output_run_name_sat_name[ii]);
	strcat(CONSTELLATION->spacecraft[ii][0].filenameecef, "/");
	strcat(CONSTELLATION->spacecraft[ii][0].filenameecef, "ECEF_");
	strcat(CONSTELLATION->spacecraft[ii][0].filenameecef, OPTIONS->filename_output[ii]);
/* 	if (iProc == 0){ */
/* 	  printf("<%s> |  %d - %d\n", CONSTELLATION->spacecraft[ii][0].filenameecef, ii, OPTIONS->n_satellites-1); */
/* 	} */
	//CONSTELLATION->spacecraft[ii][0].fpecef = fopen( CONSTELLATION->spacecraft[ii][0].filenameecef, "w+");   // (ECEF format) (CBV)

      if ( start_ensemble[ii] == 0 ){ // if main sc ii is run by this iProc
	///	if ( iProc == 1 ){
	//	printf("iProc %d opens file <%s> for sc %d: \n", iProc, CONSTELLATION->spacecraft[ii][0].filename, ii);
	//	}

	//    if (ii > 7){

      //    }
	//      printf("%d %d <%s> <%s>\n", ii, OPTIONS->n_satellites,  CONSTELLATION->spacecraft[ii][0].filename, CONSTELLATION->spacecraft[ii][0].fp );
      if ( CONSTELLATION->spacecraft[ii][0].INTEGRATOR.isGPS == 0 ){ 
      CONSTELLATION->spacecraft[ii][0].fp     = fopen( CONSTELLATION->spacecraft[ii][0].filename, "w+");      // (longer format) (CBV)
	}

CONSTELLATION->spacecraft[ii][0].fpecef = fopen( CONSTELLATION->spacecraft[ii][0].filenameecef, "w+");   // (ECEF format) (CBV)

      if ( CONSTELLATION->spacecraft[ii][0].INTEGRATOR.isGPS == 0 ){ 
      CONSTELLATION->spacecraft[ii][0].fpout  = fopen( CONSTELLATION->spacecraft[ii][0].filenameout, "w+"); // (lla format) (CBV)
	}
      if ( CONSTELLATION->spacecraft[ii][0].INTEGRATOR.isGPS == 0 ){ 
      CONSTELLATION->spacecraft[ii][0].fprho  = fopen( CONSTELLATION->spacecraft[ii][0].filenamerho, "w+"); // (lla format) (CBV)
	}
      if ( CONSTELLATION->spacecraft[ii][0].INTEGRATOR.isGPS == 0 ){
            CONSTELLATION->spacecraft[ii][0].fpatt  = fopen( CONSTELLATION->spacecraft[ii][0].filenameatt, "w+"); // (lla format) (CBV)
	}
      // line below ot comment if you don't want to use it
      /* if ( ii <  OPTIONS->n_satellites - OPTIONS->nb_gps){ */
      /* 	CONSTELLATION->spacecraft[ii][0].INTEGRATOR.file_given_output = fopen( CONSTELLATION->spacecraft[ii][0].INTEGRATOR.filename_given_output, "w+" ); */
      /* } */

      if (GROUND_STATION->nb_ground_stations > 0){
	if ( ii <  OPTIONS->n_satellites - OPTIONS->nb_gps){ // do not compute the ground station coverage for gps
	  for (iground = 0; iground < GROUND_STATION->nb_ground_stations; iground++){
	    CONSTELLATION->spacecraft[ii][0].fp_coverage_ground_station[iground] =fopen( CONSTELLATION->spacecraft[ii][0].filename_coverage_ground_station[iground], "w+");
	  }
	}
      }
      if (( strcmp( OPTIONS->type_orbit_initialisation, "tle" ) == 0) || (strcmp(OPTIONS->type_orbit_initialisation, "tle_sgp4" ) == 0 )){
	if ( CONSTELLATION->spacecraft[ii][0].INTEGRATOR.isGPS == 0 ){
	CONSTELLATION->spacecraft[ii][0].fptle = fopen( CONSTELLATION->spacecraft[ii][0].filenametle, "w+");      // (longer format) (CBV)
	}
      }
      if (ii < OPTIONS->n_satellites - OPTIONS-> nb_gps){
	if (CONSTELLATION->spacecraft[ii][0].INTEGRATOR.solar_cell_efficiency != -1){
	  CONSTELLATION->spacecraft[ii][0].fpower =fopen( CONSTELLATION->spacecraft[ii][0].filenamepower, "w+"); // power file
	  CONSTELLATION->spacecraft[ii][0].fpeclipse =fopen( CONSTELLATION->spacecraft[ii][0].filenameeclipse, "w+"); // power file
	}
      }
      } // end of if main sc ii is run by this iProc
    } // end of go through all main sc
    // end of Open SpOCK's output files (other than ensemble output files)
  
    if ( ( iDebugLevel >= 2 ) ){
      printf("-- (generate_ephemerides) Done opening output files... (iProc %d)\n", iProc);
    }



      if ( ( iDebugLevel >= 5 ) ){
	printf("-- (generate_ephemerides) Sending choose_tle_to_initialise_orbit to all proc... (iProc %d)\n", iProc);
      }

      // if this iProc does not run any main sc then it is still needs to receive choose_tle_to_initialise_orbit from the iProc that run main sc (because it needs choose_tle_to_initialise_orbit for ensemble sc)
/*       if (iProc == 0){// iProc 0 always runs at least one main sc so it sends it to the iProc that do not know choose_tle_to_initialise_orbit */
/* 	for (ccc = 0; ccc < nProcs; ccc++){ */
/* 	  if (iProc_run_main_sc[ccc] == 0){ // if iProc ccc does not run any main sc (so it does not know choose_tle_to_initialise_orbit) */
/* 	    printf("XXXXXXXXX\n%d %d\nXXXXXXX\n", iProc, choose_tle_to_initialise_orbit); */
/* 	    MPI_Send(&choose_tle_to_initialise_orbit, 1, MPI_INT, ccc, 0, MPI_COMM_WORLD); */
/* 	  } */
/* 	} */
       
/* 	if ( ( iDebugLevel >= 5 ) ){ */
/* 	  printf("-- (generate_ephemerides) Done sending choose_tle_to_initialise_orbit to all proc... (iProc %d)\n", iProc); */
/* 	} */
/*       } // end of if iproc == 0 */
/*       else{ */
/* 	if (nProcs > 1){ */
/* 	  if ( ( iDebugLevel >= 5 ) ){ */
/* 	    //    printf("-- (generate_ephemerides) All proc receiving choose_tle_to_initialise_orbit from iProc 0... (iProc %d)\n", iProc); */
/* 	  } */
/* 	  if (iProc_run_main_sc[iProc] == 0){ // if this iProc does not run any main sc (so it does not know choose_tle_to_initialise_orbit)  */
/* 	    printf("YYYYYYYYY\n%d\nYYYYYYY\n", iProc); */
/* 	    MPI_Recv(&choose_tle_to_initialise_orbit, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); */
/* 	    //            printf("Iproc %d has value %d\n", iProc, choose_tle_to_initialise_orbit); */
/* 	  } */
/* 	  if ( ( iDebugLevel >= 5 ) ){ */
/* 	    //    printf("-- (generate_ephemerides) All proc DONE receiving choose_tle_to_initialise_orbit from iProc 0... (iProc %d)\n", iProc); */
/* 	  } */

/* 	}  // end of if nProcs > 1 */

/*       } // end of else (so if iProc > 0) */

  // Open SpOCK's ensemble output files
  for (ii = 0; ii < OPTIONS->n_satellites; ii++) { // go over all sc
    if (ii < OPTIONS->n_satellites - OPTIONS-> nb_gps){ // if sc is not a GPS
      //  if ( start_ensemble[ii] == 0 ){ // if main sc ii is run by this iProc
	if ( OPTIONS->nb_ensembles_min > 0 ){ // if running ensembles

	if ( array_sc[1] > 0 )  { // if this iProc runs ensembles (otherwise array_sc[1] = - 1)
	// names of files per proc
	sprintf(temp_iproc_file, "%d", iProc+1);
	sprintf(temp_nb_proc, "%d", nProcs_that_are_gonna_run_ensembles);
	for (fff = 0; fff < OPTIONS->nb_ensembles_output_files ; fff ++){
	  if ( (strcmp(OPTIONS->filename_output_ensemble[fff], "tca" )!= 0 ) && (strcmp(OPTIONS->filename_output_ensemble[fff], "dca" )!= 0 ) && (strcmp(OPTIONS->filename_output_ensemble[fff], "sample" )!= 0 )){
			 
	    strcpy(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff], OPTIONS->dir_output_run_name_sat_name[ii]);
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff], "/");
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff], "/ensemble/iproc_");
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],temp_iproc_file);
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],"-");
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],temp_nb_proc);
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],"_");
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],OPTIONS->filename_output_ensemble[fff]);
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],"_");
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],OPTIONS->filename_output[ii]);
/* 	    if (iProc == 1){ */
/* 	      printf("iProc %d opening file <%s> for sc %d\n", iProc, CONSTELLATION->spacecraft[ii][0].filenameiproc[fff], ii); */
/* 	    } */

	    CONSTELLATION->spacecraft[ii][0].fpiproc[fff] = fopen(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],"w+");
	  }
	  else if (ii == 0){
	    if (strcmp(OPTIONS->filename_output_ensemble[fff], "tca" ) == 0){
	      strcpy(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff], OPTIONS->dir_output_run_name_collision_tca);
	    }
	    if (strcmp(OPTIONS->filename_output_ensemble[fff], "dca" ) == 0){
	      strcpy(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff], OPTIONS->dir_output_run_name_collision_dca);
	    }

	    if (strcmp(OPTIONS->filename_output_ensemble[fff], "sample" ) == 0){
	      strcpy(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff], OPTIONS->dir_output_run_name_collision_sample);
	    }

	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff], "/");
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff], "iproc_");
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],temp_iproc_file);
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],"-");
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],temp_nb_proc);
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],"_");
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],OPTIONS->filename_output_ensemble[fff]);
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],"_");
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],OPTIONS->dir_output_run_name_temp);
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],".txt");

	    if (strcmp(OPTIONS->filename_output_ensemble[fff], "tca" ) == 0){
	      tca_file = fopen(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],"w+");
	    }
	    if (strcmp(OPTIONS->filename_output_ensemble[fff], "dca" ) == 0){
	      dca_file = fopen(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],"w+");
	    }
/* 	    if (strcmp(OPTIONS->filename_output_ensemble[fff], "sample" ) == 0){ */
/* 	      sample_file = fopen(CONSTELLATION->spacecraft[ii][0].filenameiproc[fff],"w+"); */
/* 	    } */


	  }
	 
	}
	} // end of if this iProc runs ensembles (otherwise array_sc[1] = - 1)
	} // end of if running ensembles

      /* if ( OPTIONS->nb_ensembles_attitude > 0 ){ */
      /* 	if (ii == 0){ // ensemble on attitude only works if there is one satellite only in the constellation for now */
      /* 	  // names of files per proc */
      /* 	  sprintf(temp_iproc_file_attitude, "%d", iProc); */
      /* 	  sprintf(temp_nb_proc_attitude, "%d", nProcs-1); */
      /* 	  fff = 0; */
      /* 	  strcpy(CONSTELLATION->spacecraft[ii][0].filenameiproc_attitude[iProc][fff],"./output/ensemble/iproc_"); */
      /* 	  strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc_attitude[iProc][fff],temp_iproc_file_attitude); */
      /* 	  strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc_attitude[iProc][fff],"-"); */
      /* 	  strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc_attitude[iProc][fff],temp_nb_proc_attitude); */
      /* 	  strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc_attitude[iProc][fff],"_"); */
      /* 	  strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc_attitude[iProc][fff],name_file_ensembles_attitude[fff]); */
      /* 	  strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc_attitude[iProc][fff],"_"); */
      /* 	  strcat(CONSTELLATION->spacecraft[ii][0].filenameiproc_attitude[iProc][fff],OPTIONS->filename_output[ii]); */
      /* 	  CONSTELLATION->spacecraft[ii][0].fpiproc_attitude[iProc][fff] = fopen(CONSTELLATION->spacecraft[ii][0].filenameiproc_attitude[iProc][fff],"w+"); */
      /* 	} */
      /* } */
	// } // end of if main sc ii is run by this iProc
    } // end of sc is not a GPS


  } // end of go trhough all sc
  // end of Open SpOCK's ensemble output files

  // !!!!!!!!!!!!!!!!!!!!!!!!!! THE BLOCK BELOW IS TO INITIALIZE SATELLLITE 2 WITH THE CORRECT SPACING WITH RESPECT TO SATELLITE 1. IT IS JUST FOR A TRY, AND SHOULD BE REMOVED AND IMPLEMENTED IN THE CODE ITSELF. SO IF IT'S STILL HERE AFTER JUNE 10TH, 2016 THEN REMOVE IT!
  /* if (iProc == 0){ */
  /* if ( strcmp( OPTIONS->type_orbit_initialisation, "oe" ) == 0 ){ */
  /*   for (ii = 1; ii < OPTIONS->n_satellites - OPTIONS->nb_gps; ii++){ */
  /*     while ( fmod( CONSTELLATION->spacecraft[ii][0].OE.w + CONSTELLATION->spacecraft[ii][0].OE.f - ( CONSTELLATION->spacecraft[0][0].OE.w + CONSTELLATION->spacecraft[0][0].OE.f ), 2*M_PI ) < ( OPTIONS->w[ii] + OPTIONS->f[ii] ) * DEG2RAD){ */
  /* 	propagate_spacecraft( &CONSTELLATION->spacecraft[ii][0], PARAMS, starttime, &density, GROUND_STATION, OPTIONS, CONSTELLATION, iProc ); */

  /*     } */
  /*     CONSTELLATION->spacecraft[ii][0].et = starttime; */
  /*   } */
  /* } */
  /* } */
  // !!!!!!!!!!!!!!!!!!!!!!!!!! END OF THE BLOCK BELOW IS TO INITIALIZE SATELLLITE 2 WITH THE CORRECT SPACING WITH RESPECT TO SATELLITE 1. IT IS JUST FOR A TRY, AND SHOULD BE REMOVED AND IMPLEMENTED IN THE CODE ITSELF. SO IF IT'S STILL HERE AFTER JUNE 10TH, 2016 THEN REMOVE IT!



  // Write header in SpOCK's output files
  for (ii = 0; ii < OPTIONS->n_satellites; ii++) {

    //    if ( start_ensemble[ii] == 0 ){ // if main sc ii is run by this iProc
    //      print_test();
/*     if (ii > 7){ */
/*     printf("%d %d <%s>\n", ii, OPTIONS->n_satellites, OPTIONS->gps_file_name[ii]); */
/*     } */
    write_output( CONSTELLATION->spacecraft[ii], 1, choose_tle_to_initialise_orbit, ii, OPTIONS->n_satellites,OPTIONS->nb_gps,  OPTIONS->nb_ensembles_min, OPTIONS->nb_ensemble_min_per_proc,  iProc, OPTIONS->nb_ensembles_output_files, OPTIONS->filename_output_ensemble, previous_lat,OPTIONS, PARAMS->EARTH.earth_fixed_frame,1,1, ( CONSTELLATION->et ), nProcs, iDebugLevel, compute_collisions, start_ensemble, array_sc, CONSTELLATION, PARAMS);
    //print_test();
      //    }
  }
  //  print_test(); exitf();
  // end of Write header in SpOCK's output files

  // If the orbit initialization of the satellites was done with TLEs, then propagate each satellite until their epoch reaches the same as the constellation's epoch start time
  if ( (strcmp( OPTIONS->type_orbit_initialisation, "tle" ) == 0) || (strcmp(OPTIONS->type_orbit_initialisation, "tle_sgp4" ) == 0 ) || (strcmp( OPTIONS->type_orbit_initialisation, "collision_vcm" ) == 0)){
    if (iProc == 0){
      //    printf("Propagating the satellites until the epoch start time of the constellation (with no drag)...\n");
    printf("Propagating the satellites until the epoch start time of the constellation...\n");
    }
    for (ii = 0; ii < OPTIONS->n_satellites - OPTIONS->nb_gps; ii++){ // go through all main sc except gps
      //      if (iProc == 0){

// write the initial position (at the TLE epoch start)
      if ((strcmp( OPTIONS->type_orbit_initialisation, "tle" ) == 0) || (strcmp(OPTIONS->type_orbit_initialisation, "tle_sgp4" ) == 0 )){
      write_output( CONSTELLATION->spacecraft[ii] , 2, choose_tle_to_initialise_orbit, ii, OPTIONS->n_satellites,OPTIONS->nb_gps, OPTIONS->nb_ensembles_min,  OPTIONS->nb_ensemble_min_per_proc,  iProc,  OPTIONS->nb_ensembles_output_files, OPTIONS->filename_output_ensemble,previous_lat,OPTIONS,PARAMS->EARTH.earth_fixed_frame,1,1,( CONSTELLATION->et ), nProcs, iDebugLevel, compute_collisions, start_ensemble, array_sc, CONSTELLATION, PARAMS);
      }
      if ( start_ensemble[ii] == 0){ // if this iProc runs main sc ii

      save_include_drag = CONSTELLATION->spacecraft[ii][0].INTEGRATOR.include_drag;
      save_include_solar_pressure = CONSTELLATION->spacecraft[ii][0].INTEGRATOR.include_solar_pressure;
            save_include_earth_pressure = CONSTELLATION->spacecraft[ii][0].INTEGRATOR.include_earth_pressure;
      save_solar_cell_efficiency = CONSTELLATION->spacecraft[ii][0].INTEGRATOR.solar_cell_efficiency;
/*       CONSTELLATION->spacecraft[ii][0].INTEGRATOR.include_drag  = 0; // do not include drag until the epoch becomes the epoch start time of the constellation (for now! (this is for interpolation of the thermospheric drivers reasons)) */
/*       CONSTELLATION->spacecraft[ii][0].INTEGRATOR.include_solar_pressure  = 0; // do not include solar_pressure until the epoch becomes the epoch start time of the constellation (for now! (this is for interpolation of the thermospheric drivers reasons)) */
/*       CONSTELLATION->spacecraft[ii][0].INTEGRATOR.solar_cell_efficiency  = -1; // do not include solar_power until the epoch becomes the epoch start time of the constellation */

      //    } end of if iProc == 0
      } // end of if this iProc runs main sc ii
      if (ii < OPTIONS->n_satellites - OPTIONS-> nb_gps){ // if main sc is not a GPS
	if ( OPTIONS->nb_ensembles_min > 0 ) { // if running ensembles
	  for (eee = 1L + iProc * OPTIONS->nb_ensemble_min_per_proc; eee< 1 + iProc * OPTIONS->nb_ensemble_min_per_proc + OPTIONS->nb_ensemble_min_per_proc ; eee++){ // go through all ensembles
	    save_include_drag = CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_drag;
	    save_include_solar_pressure = CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_solar_pressure;
	    	    save_include_earth_pressure = CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_earth_pressure;
	    save_solar_cell_efficiency = CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.solar_cell_efficiency;

	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.dt = OPTIONS->dt;
	    	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.dt_pos_neg = OPTIONS->dt_pos_neg;
/* 	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_drag  = 0; // do not include drag until the epoch becomes the epoch start time of the constellation (for now! (this is for interpolation of the thermospheric drivers reasons)) */
/* 	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_solar_pressure  = 0; // do not include solar_pressure until the epoch becomes the epoch start time of the constellation (for now! (this is for interpolation of the thermospheric drivers reasons)) */
/* 	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.solar_cell_efficiency  = -1; // do not include solar_power until the epoch becomes t */

	  } // end of go through all ensembles
	} // end of if running ensembles
      } // end of if main sc is not a GPS
      //      if (iProc == 0){
      if (strcmp(OPTIONS->type_orbit_initialisation, "tle_sgp4" ) != 0 ){ // if tle_sgp4 then we directly jump from the tle epoch to the constellation start time (ie we don't calculate the r/v at every time step between those two times)
      if ( start_ensemble[ii] == 0){ // if this iProc runs main sc ii
	while ( ( ( CONSTELLATION->spacecraft[ii][0].et - twrite ) < 0 ) && ( ( twrite - CONSTELLATION->spacecraft[ii][0].et ) > CONSTELLATION->spacecraft[ii][0].INTEGRATOR.dt_pos_neg ) ){
	  //	  etprint(CONSTELLATION->spacecraft[ii][0].et, "");
	  	  	  if (ii == OPTIONS->which_sc_oldest_tle_epoch){
			    //	    	    	    print_progress_epoch_sc_to_epoch_constellation(starttime, CONSTELLATION->spacecraft[ii][0].et, OPTIONS->et_oldest_tle_epoch, iProc, OPTIONS->nb_gps);

			    	  }

	  propagate_spacecraft( &CONSTELLATION->spacecraft[ii][0], PARAMS, starttime, OPTIONS->et_oldest_tle_epoch, &density, GROUND_STATION, OPTIONS, CONSTELLATION, iProc, iDebugLevel,  start_ensemble, array_sc ); // don't care about starttime here because this is used for linear interpolation with the density drivers (F10.7, Ap, ...) and the attitude, which we do not care for the satellites from their TLE epoch time to the constellation epoch start time (for now!) (this is for interpolation of the thermospheric drivers reasons)
	  if ((strcmp( OPTIONS->type_orbit_initialisation, "tle" ) == 0) || (strcmp(OPTIONS->type_orbit_initialisation, "tle_sgp4" ) == 0 ) ){
	    write_output( CONSTELLATION->spacecraft[ii] , 2, choose_tle_to_initialise_orbit, ii, OPTIONS->n_satellites,OPTIONS->nb_gps, OPTIONS->nb_ensembles_min,  OPTIONS->nb_ensemble_min_per_proc,  iProc,  OPTIONS->nb_ensembles_output_files, OPTIONS->filename_output_ensemble,previous_lat,OPTIONS,PARAMS->EARTH.earth_fixed_frame,1,1, ( CONSTELLATION->et + OPTIONS->dt_pos_neg ), nProcs, iDebugLevel, compute_collisions, start_ensemble, array_sc, CONSTELLATION, PARAMS);
      }
	} // end while propagation
      } // end of if this iProc runs main sc ii
      //      } // end of if iProc == 0

      if ( OPTIONS->nb_ensembles_min > 0 ) { // if running ensembles

	while ( ( ( CONSTELLATION->spacecraft[ii][1L + iProc * OPTIONS->nb_ensemble_min_per_proc].et - twrite ) < 0 ) && ( ( twrite - CONSTELLATION->spacecraft[ii][1L + iProc * OPTIONS->nb_ensemble_min_per_proc].et ) > CONSTELLATION->spacecraft[ii][1L + iProc * OPTIONS->nb_ensemble_min_per_proc].INTEGRATOR.dt_pos_neg ) ){

	  // ENSEMBLES SPACECRAFT (BUT NOT FOR GPS SATELLITES)
	  if (ii < OPTIONS->n_satellites - OPTIONS-> nb_gps){
	    if ( OPTIONS->nb_ensembles_min > 0 ) {
	      for (eee = 1L + iProc * OPTIONS->nb_ensemble_min_per_proc; eee< 1 + iProc * OPTIONS->nb_ensemble_min_per_proc + OPTIONS->nb_ensemble_min_per_proc ; eee++){
	       
		propagate_spacecraft( &CONSTELLATION->spacecraft[ii][eee], PARAMS, starttime, OPTIONS->et_oldest_tle_epoch, &density, GROUND_STATION, OPTIONS, CONSTELLATION, iProc, iDebugLevel, start_ensemble, array_sc );

	      } // end of loop over all ensembles for each iProc<
	    } // enf if ensembles are run
	  } // end if sc is not gps
	} // end while of propagatoin
      } // end of if running ensembles
      } // end of if not using tle_sgp4

      // propagate ONE LAST TIME UNTIL THE SATELLITE EPOCH IS THE SAME AS THE CONSTELLATION EPOCH
      //      if (iProc == 0){
            if (fabs(starttime-CONSTELLATION->spacecraft[ii][0].et) > OPTIONS->dt/1.e6){ //// if sc.et not equal to constellation initial epoch
      if (start_ensemble[ii] == 0){ // if this iProc runs main sc ii
	CONSTELLATION->spacecraft[ii][0].INTEGRATOR.dt = fabs(starttime - CONSTELLATION->spacecraft[ii][0].et);
	CONSTELLATION->spacecraft[ii][0].INTEGRATOR.dt_pos_neg = starttime - CONSTELLATION->spacecraft[ii][0].et;
	propagate_spacecraft( &CONSTELLATION->spacecraft[ii][0], PARAMS, starttime, OPTIONS->et_oldest_tle_epoch, &density, GROUND_STATION, OPTIONS, CONSTELLATION, iProc, iDebugLevel, start_ensemble, array_sc );
      }      //      }
      // // ENSEMBLES SPACECRAFT (BUT NOT FOR GPS SATELLITES)
      if (ii < OPTIONS->n_satellites - OPTIONS-> nb_gps){ // if main sc is not a gps
	if ( OPTIONS->nb_ensembles_min > 0 ) { // if running ensembles
	  for (eee = 1L + iProc * OPTIONS->nb_ensemble_min_per_proc; eee< 1 + iProc * OPTIONS->nb_ensemble_min_per_proc + OPTIONS->nb_ensemble_min_per_proc ; eee++){ // go over all ensemble sc
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.dt = fabs(starttime - CONSTELLATION->spacecraft[ii][eee].et);
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.dt_pos_neg = starttime - CONSTELLATION->spacecraft[ii][eee].et;
	    propagate_spacecraft( &CONSTELLATION->spacecraft[ii][eee], PARAMS, starttime, OPTIONS->et_oldest_tle_epoch, &density, GROUND_STATION, OPTIONS, CONSTELLATION, iProc, iDebugLevel,  start_ensemble, array_sc );
	  } // end of go over all ensemble sc
	} // end of if running ensembles
      } // end of if main sc is not a gps
      if ((strcmp( OPTIONS->type_orbit_initialisation, "tle" ) == 0) || (strcmp(OPTIONS->type_orbit_initialisation, "tle_sgp4" ) == 0 )){
      write_output( CONSTELLATION->spacecraft[ii] , 0, choose_tle_to_initialise_orbit, ii, OPTIONS->n_satellites,OPTIONS->nb_gps, OPTIONS->nb_ensembles_min,  OPTIONS->nb_ensemble_min_per_proc,  iProc,  OPTIONS->nb_ensembles_output_files, OPTIONS->filename_output_ensemble,previous_lat,OPTIONS,PARAMS->EARTH.earth_fixed_frame,1,1, ( CONSTELLATION->et + OPTIONS->dt_pos_neg ), nProcs, iDebugLevel, compute_collisions, start_ensemble, array_sc, CONSTELLATION, PARAMS);
      }
      // END OF PROPAGATE ONE LAST TIME UNTIL THE SATELLITE EPOCH IS THE SAME AS THE CONSTELLATION EPOCH
  } // if sc.et not "equal" to constellation initial epoch
      // GO BACK TO DT AND FORCES CHOSEN IN THE INPUT FILE
      //      if (iProc == 0){
      if (start_ensemble[ii] == 0){ // if this iProc runs main sc ii
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.dt = OPTIONS->dt;
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.dt_pos_neg = OPTIONS->dt_pos_neg;
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.include_drag  = save_include_drag; // the drag is now computed (if the user decided so) because the epoch reached the epoch start time of the constellation
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.include_solar_pressure  = save_include_solar_pressure; // the solar_pressure is now computed (if the user decided so) because the epoch reached the epoch start time of the constellation
            CONSTELLATION->spacecraft[ii][0].INTEGRATOR.include_earth_pressure  = save_include_earth_pressure; // the earth_pressure is now computed (if the user decided so) because the epoch reached the epoch start time of the constellation
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.solar_cell_efficiency  = save_solar_cell_efficiency; // the solar_power is now computed (if the user decided so) because the epoch reached the epoch start time of the constellation
      } // end of if this iProc runs main sc ii
      //      }
      // // ENSEMBLES SPACECRAFT (BUT NOT FOR GPS SATELLITES)
      if (ii < OPTIONS->n_satellites - OPTIONS-> nb_gps){
	if ( OPTIONS->nb_ensembles_min > 0 ) {
	  for (eee = 1L + iProc * OPTIONS->nb_ensemble_min_per_proc; eee< 1 + iProc * OPTIONS->nb_ensemble_min_per_proc + OPTIONS->nb_ensemble_min_per_proc ; eee++){
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.dt = OPTIONS->dt;
	    	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.dt_pos_neg = OPTIONS->dt_pos_neg;
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_drag  = save_include_drag;
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_solar_pressure  = save_include_solar_pressure;
	    	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_earth_pressure  = save_include_earth_pressure;
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.solar_cell_efficiency  = save_solar_cell_efficiency;
	  }
	}
      }
      //      if (iProc == 0){
      if (start_ensemble[ii] == 0){ // if this iProc runs main sc ii
	if ( CONSTELLATION->spacecraft[ii][0].INTEGRATOR.isGPS == 0 ){
	  if ((strcmp( OPTIONS->type_orbit_initialisation, "tle" ) == 0) || (strcmp(OPTIONS->type_orbit_initialisation, "tle_sgp4" ) == 0 ) ){
	fclose(CONSTELLATION->spacecraft[ii][0].fptle);
	  }
	}
      }
      //      }
    } // end of go through all main sc except gps
  } // end of if the orbit initialization of the satellites was done with TLEs, then propagate each satellite until their epoch reaches the same as the constellation's epoch start time

  // Propagate each GPS satellite until their epoch reaches the same as the constellation's epoch start time
  //  if (iProc == 0){
  if (OPTIONS->nb_gps > 0){
if (iProc == 0){
    printf("Propagating the GPS satellites until the epoch start time of the constellation...\n");
 }
    for (ii = OPTIONS->n_satellites - OPTIONS->nb_gps; ii < OPTIONS->n_satellites; ii++){ // go through all GPS sc
      if (start_ensemble[ii] == 0){ // if this iProc runs main sc ii

	time_between_last_gps_epoch_and_constellation_epoch_starttime = -1.0;
	new_dt_for_gps = -1.0;
	while ( ( ( CONSTELLATION->spacecraft[ii][0].et - twrite ) < 0 ) && ( ( twrite - CONSTELLATION->spacecraft[ii][0].et ) > CONSTELLATION->spacecraft[ii][0].INTEGRATOR.dt_pos_neg ) ){
	  propagate_spacecraft( &CONSTELLATION->spacecraft[ii][0], PARAMS, starttime, OPTIONS->et_oldest_tle_epoch, &density, GROUND_STATION, OPTIONS, CONSTELLATION, iProc , iDebugLevel, start_ensemble, array_sc); // don't care about starttime here because this is used for linear interpolation with the density drivers (F10.7, Ap, ...) and the attitude, which we do not care for the GPS satellites (for now!)
	}

	time_between_last_gps_epoch_and_constellation_epoch_starttime =fabs( starttime - CONSTELLATION->spacecraft[ii][0].et);
	new_dt_for_gps = time_between_last_gps_epoch_and_constellation_epoch_starttime;
	CONSTELLATION->spacecraft[ii][0].INTEGRATOR.dt = new_dt_for_gps;
	propagate_spacecraft( &CONSTELLATION->spacecraft[ii][0], PARAMS, starttime, OPTIONS->et_oldest_tle_epoch, &density, GROUND_STATION, OPTIONS, CONSTELLATION, iProc, iDebugLevel, start_ensemble, array_sc );

	write_output( CONSTELLATION->spacecraft[ii] , 0, choose_tle_to_initialise_orbit, ii, OPTIONS->n_satellites,OPTIONS->nb_gps,OPTIONS->nb_ensembles_min,  OPTIONS->nb_ensemble_min_per_proc,  iProc,  OPTIONS->nb_ensembles_output_files, OPTIONS->filename_output_ensemble, previous_lat,OPTIONS,PARAMS->EARTH.earth_fixed_frame,1,1, ( CONSTELLATION->et + OPTIONS->dt_pos_neg ), nProcs, iDebugLevel, compute_collisions, start_ensemble, array_sc, CONSTELLATION, PARAMS);

	CONSTELLATION->spacecraft[ii][0].INTEGRATOR.dt = OPTIONS->dt;
		CONSTELLATION->spacecraft[ii][0].INTEGRATOR.dt_pos_neg = OPTIONS->dt_pos_neg;
      } // end of if this iProc runs main sc ii
    } // end of going through all GPS sc
  } // end of if (OPTIONS->nb_gps > 0)
  //  } // end of iProc == 0
  // end of Propagate each GPS satellite until their epoch reaches the same as the constellation's epoch start time



 if (compute_collisions == 1){  // if computing the probability of collision
    if (nProcs > OPTIONS->nb_ensembles_min){
      print_error_any_iproc(iProc, "If computing the probability of collision, the number of processors used in MPI has to be smaller than (or equal to) the number of ensembles (this bug will be fixed in the future)");
    }


   if ( ( iDebugLevel >= 2 ) ){
     printf("--- (generate_ephemerides) Allocating memory for some of the collision assessment variables. (iProc %d) \n", iProc);
   }

   save_tca = malloc( OPTIONS->nb_time_steps * sizeof( double ) );
   if ( save_tca == NULL ){
     print_error_any_iproc(iProc, "No memory space for save_tca");
   }
   et_time_step_of_save_tca = malloc( OPTIONS->nb_time_steps * sizeof( double ) );
   if ( et_time_step_of_save_tca == NULL ){
     print_error_any_iproc(iProc, "No memory space for et_time_step_of_save_tca");
   }

   save_dca = malloc( OPTIONS->nb_time_steps * sizeof( double ) );
   if ( save_dca == NULL ){
     print_error_any_iproc(iProc, "No memory space for save_dca");
   }
   if (( iDebugLevel >= 2 ) ){
     printf("--- (generate_ephemerides) Done allocating memory for some of the collision assessment variables. (iProc %d)\n", iProc);
   }

 } // end of if computing the probability of collision



 // if (iProc == 0){
   if (compute_collisions == 1){ // if computing the probability of collision 2020-11-13 i don't want to worry about dt_pos_neg (backward propagation) for collision avoidance
     if ( ( iDebugLevel >= 2 ) ){
       printf("--- (generate_ephemerides) Allocating memory for other collision assessment variables (only by main node). (iProc %d)\n", iProc);
     }

     // allocate memory for previous_eci_r[nsatref][3]:
     previous_eci_r = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof( double * ) );
     if (  previous_eci_r == NULL ){
       print_error(iProc, "Not enough memory for previous_eci_r");
     }
     for ( isc_ref = 0; isc_ref < OPTIONS->nb_satellites_not_including_gps ; isc_ref++){
       previous_eci_r[isc_ref] = malloc( 3 * sizeof( double * ) );
       if (  previous_eci_r[isc_ref] == NULL ){
	 print_error(iProc, "Not enough memory for previous_eci_r");
       }
     }

     // allocate memory for previous_eci_v[nsatref][3]:
     previous_eci_v = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof( double * ) );
     if (  previous_eci_v == NULL ){
       print_error(iProc, "Not enough memory for previous_eci_v");
     }
     for ( isc_ref = 0; isc_ref < OPTIONS->nb_satellites_not_including_gps ; isc_ref++){
       previous_eci_v[isc_ref] = malloc( 3 * sizeof( double * ) );
       if (  previous_eci_v[isc_ref] == NULL ){
	 print_error(iProc, "Not enough memory for previous_eci_v");
       }
     }

     // allocate memory for previous_eci_a[nsatref][3]:
     previous_eci_a = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof( double * ) );
     if (  previous_eci_a == NULL ){
       print_error(iProc, "Not enough memory for previous_eci_a");
     }
     for ( isc_ref = 0; isc_ref < OPTIONS->nb_satellites_not_including_gps ; isc_ref++){
       previous_eci_a[isc_ref] = malloc( 3 * sizeof( double * ) );
       if (  previous_eci_a[isc_ref] == NULL ){
	 print_error(iProc, "Not enough memory for previous_eci_a");
       }
     }


     total_collision_each_dt = malloc( OPTIONS->nb_time_steps * sizeof(double));
     for (ccc = 0; ccc < OPTIONS->nb_time_steps ; ccc++){
       total_collision_each_dt[ccc] = 0;
     }
 
     //////////// IF COLLISION ASSESSMENT IS ON, FIRST PROPAGATE ONLY THE UNPERTURBED ORBITS
     // If collision assessment is on, then first propagate all reference sc (unpertubed orbits) to compute times of close approach (and allocate memory for variables)
     save_first_distance_unpertubed_orbit = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof( double *) );
     if (save_first_distance_unpertubed_orbit == NULL){
       print_error_any_iproc(iProc, "No memory for save_first_distance_unpertubed_orbit");
     }
     for (ii = 0; ii < OPTIONS->nb_satellites_not_including_gps; ii++) {
       save_first_distance_unpertubed_orbit[ii] = malloc(3 * sizeof(double));
       if (save_first_distance_unpertubed_orbit[ii] == NULL){
	 print_error_any_iproc(iProc, "No memory for save_first_distance_unpertubed_orbit");
       }
     }
     save_last_distance_unpertubed_orbit = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof( double *) );
     if (save_last_distance_unpertubed_orbit == NULL){
       print_error_any_iproc(iProc, "No memory for save_last_distance_unpertubed_orbit");
     }
     for (ii = 0; ii < OPTIONS->nb_satellites_not_including_gps; ii++) {
       save_last_distance_unpertubed_orbit[ii] = malloc(3 * sizeof(double));
       if (save_last_distance_unpertubed_orbit[ii] == NULL){
	 print_error_any_iproc(iProc, "No memory for save_last_distance_unpertubed_orbit");
       }
     }

     for (ii = 0; ii < OPTIONS->nb_satellites_not_including_gps; ii++) {
       if ( start_ensemble[ii] == 0){ // if this iProc runs main sc ii
       v_copy( save_first_distance_unpertubed_orbit[ii],  CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL);
       }
     }

     if (( iDebugLevel >= 2 ) ){
       printf("--- (generate_ephemerides) Done allocating memory for other collision assessment variables (only by main node). (iProc %d)\n", iProc);
     }

     if (( iDebugLevel >= 2 ) ){
       printf("--- (generate_ephemerides) Propagating all unpertubed orbits to determine times of close approach (only by main node). (iProc %d)\n", iProc);
     }

     OPTIONS->first_run_to_find_tca_before_collision_assessment = 1;

     while ( ( fabs(CONSTELLATION->et - endtime) < 0.01 ) && ( min_altitude_constellation > 100.0) ){ // propagate all unrpertubed orbits to determine times of close approach
       already_propagated_ref_sc = 1;


       // PRINT PROGRESS TO SCREEN
       if (iProc == 0){
              printf("\033[A\33[2K\rPropagating the unpertubed orbits to compute times of close approach... %.0f%%\n", ( CONSTELLATION->et - starttime ) / fabs( endtime - starttime ) *100.0);
       }

       for (ii = 0; ii < OPTIONS->nb_satellites_not_including_gps; ii++) {
	 if ( start_ensemble[ii] == 0){ // if this iProc runs main sc ii
	 previous_lat = CONSTELLATION->spacecraft[ii][0].GEODETIC.latitude*RAD2DEG;
	 v_copy( previous_eci_r[ii], CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL);
	 v_copy( previous_eci_v[ii], CONSTELLATION->spacecraft[ii][0].v_i2cg_INRTL);
	 v_copy( previous_eci_a[ii], CONSTELLATION->spacecraft[ii][0].a_i2cg_INRTL);


	 /* printf("\n"); */
	 /* etprint(CONSTELLATION->et, "t"); */
	 /* pti(ii, "ii"); */
	 /* printf("(%.10f; %.10f; %.10f) (%.10f; %.10f; %.10f)\n", CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL[0]*1000., CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL[1]*1000., CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL[2]*1000., CONSTELLATION->spacecraft[ii][0].v_i2cg_INRTL[0]*1000., CONSTELLATION->spacecraft[ii][0].v_i2cg_INRTL[1]*1000., CONSTELLATION->spacecraft[ii][0].v_i2cg_INRTL[2]*1000.); */
	  propagate_spacecraft( &CONSTELLATION->spacecraft[ii][0], PARAMS, starttime, OPTIONS->et_oldest_tle_epoch, &density, GROUND_STATION, OPTIONS, CONSTELLATION, iProc, iDebugLevel, start_ensemble, array_sc );

	 /* printf("\n"); */
	 /* etprint(CONSTELLATION->et + OPTIONS->dt, "time"); */
	 /* v_print(CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL, "r"); */
	 /* v_print(CONSTELLATION->spacecraft[ii][0].v_i2cg_INRTL, "v"); */
	 
	 // iProc 0 needs to send previous_eci_r/v/a and r/v/a_i2cg_INRT to iProc 1 so that iProc 1 can find if there is a close approach. But this is only if iProc 1 runs main sc 1. If it's iProc 0 that runs main sc 1, then no need to send anything. // !!!!!!!!!! FOR NOW WORKS ONLY IF TWO REFERENCE SATELLIES ONLY
	 if ( (iProc == 0) && ( which_iproc_is_running_main_sc[1] == 1 ) ){
	   for (ccc = 0; ccc < 3; ccc++){

	   MPI_Send(&previous_eci_r[0][ccc], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
	   MPI_Send(&previous_eci_v[0][ccc], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
	   MPI_Send(&previous_eci_a[0][ccc], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
	   MPI_Send(&CONSTELLATION->spacecraft[0][0].r_i2cg_INRTL[ccc], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
	   MPI_Send(&CONSTELLATION->spacecraft[0][0].v_i2cg_INRTL[ccc], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
	   MPI_Send(&CONSTELLATION->spacecraft[0][0].a_i2cg_INRTL[ccc], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
	   }
	 }

	 if (ii == 1){ // !!!!!!!!!! FOR NOW WORKS ONLY IF TWO REFERENCE SATELLIES ONLY
	   // iProc 1 needs to receive previous_eci_r/v/a and r/v/a_i2cg_INRT from iProc 0 so that iProc 1 can find if there is a close approach. But this is only if iProc 1 runs main sc 1. If it's iProc 0 that runs main sc 1, then no need to receive anything
	   if ( (iProc == 1) && ( which_iproc_is_running_main_sc[1] == 1 ) ){
	     for (ccc = 0; ccc < 3; ccc++){
	       MPI_Recv(&previous_eci_r[0][ccc], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	       MPI_Recv(&previous_eci_v[0][ccc], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	       MPI_Recv(&previous_eci_a[0][ccc], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	       MPI_Recv(&CONSTELLATION->spacecraft[0][0].r_i2cg_INRTL[ccc], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	       MPI_Recv(&CONSTELLATION->spacecraft[0][0].v_i2cg_INRTL[ccc], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	       MPI_Recv(&CONSTELLATION->spacecraft[0][0].a_i2cg_INRTL[ccc], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	     }
	   }
	   tca1 = -1; dca1 = 1e6; tca2 = -1; dca2 = 1e6; tca3 = -1; dca3 = 1e6;
	   close_approach_exists = 0;

	   if ( ( CONSTELLATION->et + OPTIONS->dt > starttime + 3 * OPTIONS->dt )  && (CONSTELLATION->et + OPTIONS->dt <  endtime - OPTIONS->dt)  ){ // To use the collision assessment algorithm, the propagation needs to start at least three 3 time steps before TCA and needs to end at least one time step after TCA. Therefore, do not try to find TCA earlier than 3 times steps after initial epoch or later than 1 time step before final epoch. 2020-11-13 i don't want to worry about dt_pos_neg (backward propagation) for collision avoidance
	     //	     	     print_test(); printf("XXXXXXXXXXXXXX\n%d\nXXXXXXXXX\n",iProc);

	     ancas_existence_of_min_using_two_values_of_function_and_two_values_of_its_derivative( CONSTELLATION->et, OPTIONS->dt, previous_eci_r[ii-1], previous_eci_v[ii-1], previous_eci_a[ii-1], CONSTELLATION->spacecraft[ii-1][0].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii-1][0].v_i2cg_INRTL, CONSTELLATION->spacecraft[ii-1][0].a_i2cg_INRTL, previous_eci_r[ii], previous_eci_v[ii], previous_eci_a[ii], CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][0].v_i2cg_INRTL, CONSTELLATION->spacecraft[ii][0].a_i2cg_INRTL, &close_approach_exists, &gamma0, &gamma1, &gamma2, &gamma3);

	   }


	   if ( close_approach_exists == 1 ){

	     ancas(OPTIONS->min_dist_close_approach, CONSTELLATION->et, OPTIONS->dt, previous_eci_r[ii-1], previous_eci_v[ii-1], previous_eci_a[ii-1], CONSTELLATION->spacecraft[ii-1][0].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii-1][0].v_i2cg_INRTL, CONSTELLATION->spacecraft[ii-1][0].a_i2cg_INRTL, previous_eci_r[ii], previous_eci_v[ii], previous_eci_a[ii], CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][0].v_i2cg_INRTL, CONSTELLATION->spacecraft[ii][0].a_i2cg_INRTL, gamma0, gamma1, gamma2, gamma3, &tca1, &dca1, &tca2, &dca2, &tca3, &dca3 );
	   }
	   if ( tca1 > -1 ){
	     save_tca[nb_tca] = tca1;
	     et_time_step_of_save_tca[nb_tca] = CONSTELLATION->et;
	     save_dca[nb_tca] = dca1;
	     //	              etprint(save_tca[nb_tca], "tca1");
	     nb_tca = nb_tca + 1;
	   }
	   if ( tca2 > -1 ){
	     save_tca[nb_tca] = tca2;
	     et_time_step_of_save_tca[nb_tca] = CONSTELLATION->et;
	     save_dca[nb_tca] = dca2;
	     //      	    etprint(save_tca[nb_tca], "tca2");
	     nb_tca = nb_tca + 1;
	   }
	   if ( tca3 > -1 ){
	     save_tca[nb_tca] = tca3;
	     et_time_step_of_save_tca[nb_tca] = CONSTELLATION->et;
	     save_dca[nb_tca] = dca3;
	     //  etprint(save_tca[nb_tca], "tca3");
	     nb_tca = nb_tca + 1;
	   }
	 } // end of if ii == 1

	 // compute min altitude at this time step
	 
/* 	 if (CONSTELLATION->spacecraft[ii][0].GEODETIC.altitude < */
/* 	     min_altitude_constellation) */
/* 	   min_altitude_constellation = */
/* 	     CONSTELLATION->spacecraft[ii][0].GEODETIC.altitude; */

	      

	 // WRITE OUTPUT
	 if ( (CONSTELLATION->spacecraft[ii][0].et - twrite) >= OPTIONS->dt - 0.01) { // write output
	   if ( ( fmod( CONSTELLATION->spacecraft[ii][0].et - starttime, OPTIONS->dt_output ) < OPTIONS->dt / 2. ) || ( fabs( fmod( CONSTELLATION->spacecraft[ii][0].et - starttime, OPTIONS->dt_output ) - OPTIONS->dt_output ) < OPTIONS->dt / 2. ) || ( CONSTELLATION->spacecraft[ii][0].et > endtime - 0.01) )  {


	     write_output( CONSTELLATION->spacecraft[ii] , 0, choose_tle_to_initialise_orbit,ii, OPTIONS->n_satellites,OPTIONS->nb_gps,  OPTIONS->nb_ensembles_min, OPTIONS->nb_ensemble_min_per_proc,  iProc,  OPTIONS->nb_ensembles_output_files, OPTIONS->filename_output_ensemble, previous_lat,OPTIONS,PARAMS->EARTH.earth_fixed_frame,1,0, ( CONSTELLATION->et + OPTIONS->dt_pos_neg) , nProcs, iDebugLevel, compute_collisions, start_ensemble, array_sc, CONSTELLATION, PARAMS);
	     et2utc_c(CONSTELLATION->spacecraft[ii][0].et, "ISOC" ,6 ,255 , times_att);
/* 	     if (write_density == 1){ */
/* 	       	       fprintf(density_file,"%s %e %f %f %e\n", times_att, density, CONSTELLATION->spacecraft[ii][0].INTEGRATOR.Ta, CONSTELLATION->spacecraft[ii][0].INTEGRATOR.cd_tot_norm, CONSTELLATION->spacecraft[ii][0].INTEGRATOR.A_ref_tot/1000000.); // in output A_ref_tot in km^2 */
/* 	     } */

	   }
	   // !!!!!! twrite used to be this below before the parallel constellation version
/* 	   if ( ii == (OPTIONS->n_satellites - 1) ) { */
/* 	     twrite = CONSTELLATION->spacecraft[ii][0].et; */
/* 	   } */
	   // !!!!!! end of twrite used to be this below before the parallel constellation version
	 } // end of write output

	 /* etprint(CONSTELLATION->et + OPTIONS->dt, "t"); */
	 /* pti(ii, "ii"); */
	 /* v_print(CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL, "r"); */
	 /* v_print(CONSTELLATION->spacecraft[ii][0].v_i2cg_INRTL, "v"); */
	
	 } // end of if this iProc runs main sc ii
       } //  END OF LOOP OVER ALL MAIN SATELLITES
       // UPDATE THE CONSTELLATION TIME
       CONSTELLATION->et = CONSTELLATION->et + OPTIONS->dt;
       twrite = CONSTELLATION->et ;
       //     print_test();
     } // end of  propagate all unrpertubed orbits to determine times of close approach a

     if ( ( iDebugLevel >= 2 ) ){
       printf("--- (generate_ephemerides) Done propagating all unpertubed orbits to determine times of close approach (only by main node). (iProc %d)\n", iProc);
     }

     for (ii = 0; ii < OPTIONS->nb_satellites_not_including_gps; ii++) {
       if ( start_ensemble[ii] == 0){ // if this iProc runs main sc ii
       v_copy( save_last_distance_unpertubed_orbit[ii],  CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL);
       }
     }

     //     if nb_tca, save_dca, et_time_step_of_save_tca, save_tca were calculated by iProc 1 then send them to iProc 0 (and iProc 0 needs to receive them). They were calculated by iProc if main sc 1 is run by iProc 1. Otherwise, they were calculated by iProc 0. Also, need to send save_last_distance_unpertubed_orbit[1] and save_first_distance_unpertubed_orbit[1] to iProc 0 (if save_last_distance_unpertubed_orbit[1] and save_first_distance_unpertubed_orbit[1] was calculated by iProc 1, so if main sc 1 was run by iProc 1). Note: iProc > 1 will never get here because main sc 1 (or 0) can never be run by iProc > 1. However, when looking at the close approach between the ensembles, iProc > 1 will be doing stuff.
     if ( which_iproc_is_running_main_sc[1] == 1 ){
       if (iProc == 1){
	 MPI_Send(&nb_tca, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	 for (ccc = 0; ccc < nb_tca; ccc++){
	   MPI_Send(&save_dca[ccc], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	   MPI_Send(&save_tca[ccc], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	   MPI_Send(&et_time_step_of_save_tca[ccc], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	 }
	 for (ccc = 0; ccc < 3; ccc++){
	   MPI_Send(&save_last_distance_unpertubed_orbit[1][ccc], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	   MPI_Send(&save_first_distance_unpertubed_orbit[1][ccc], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	 }
       }
       if (iProc == 0){
	 MPI_Recv(&nb_tca, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	 for (ccc = 0; ccc < nb_tca; ccc++){
	   MPI_Recv(&save_dca[ccc], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	   MPI_Recv(&save_tca[ccc], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	   MPI_Recv(&et_time_step_of_save_tca[ccc], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	 }
	 for (ccc = 0; ccc < 3; ccc++){
	   MPI_Recv(&save_last_distance_unpertubed_orbit[1][ccc], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	   MPI_Recv(&save_first_distance_unpertubed_orbit[1][ccc], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	 }
       }
     }

   
     // // Reinitialize time to propagate ensembles satellites from initial epoch to furthest TCA ( + CONSTELLATION->collision_time_span )
     CONSTELLATION->et = starttime;
     twrite = CONSTELLATION->et;
     if (iProc == 0){ // only need one proc to do the job here
     // We look at the distance at initial epoch and the distance at final epoch (unless a time span of a tca includes them, in which case we don't look at the distance at initial epoch and the distance at final epoch). If these distances are smaller than OPTIONS->min_dist_close_approach, then we consider them as a TCA. Indeed, it is not because the 2 unperturbed have not reached a closest approach in the propagation that they should be ignored in the probability of collision.
     if ( nb_tca <= 0){ // if no close approach between initial and final epoch
       // !!!!!!!!!!!!!!!!!! delete next line and uncomment below ti
       //	if ( distance_between_two_sc( save_first_distance_unpertubed_orbit[0], save_first_distance_unpertubed_orbit[1] ) < 1e6){ // if the distance at initial epoch is smaller than OPTIONS->min_dist_close_approach then consider this at a tca (that starts at inital epoch + period / 2)
       if ( distance_between_two_sc( save_first_distance_unpertubed_orbit[0], save_first_distance_unpertubed_orbit[1] ) < OPTIONS->min_dist_close_approach){ // if the distance at initial epoch is smaller than OPTIONS->min_dist_close_approach then consider this at a tca (that starts at inital epoch + period / 2)
	 save_tca[nb_tca] = starttime + CONSTELLATION->collision_time_span / 2.;
	 et_time_step_of_save_tca[nb_tca] = starttime + ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt );
	 save_dca[nb_tca] = distance_between_two_sc( save_first_distance_unpertubed_orbit[0], save_first_distance_unpertubed_orbit[1] ); // !!!!!!!!! this is incorrect. What would be correct would be to take the distance at inital epoch + period / 2 but we did not save this distance. Since save_dca is never used after, we don't really care, except in the collision output file where dca is reported
	 nb_tca = nb_tca + 1;
	  
       }
       // !!!!!!!!!!!!!!!!!! delete next line and uncomment below it
       //if (distance_between_two_sc( save_last_distance_unpertubed_orbit[0], save_last_distance_unpertubed_orbit[1] ) < 1e6) {  // if the distance at final epoch is smaller than OPTIONS->min_dist_close_approach then consider this at a tca (that starts at final epoch - period / 2)
       if (distance_between_two_sc( save_last_distance_unpertubed_orbit[0], save_last_distance_unpertubed_orbit[1] ) < OPTIONS->min_dist_close_approach) {  // if the distance at final epoch is smaller than OPTIONS->min_dist_close_approach then consider this at a tca (that starts at final epoch - period / 2)
	 save_tca[nb_tca] = endtime - CONSTELLATION->collision_time_span / 2.;
	 et_time_step_of_save_tca[nb_tca] = endtime - ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt );
	 save_dca[nb_tca] = distance_between_two_sc( save_last_distance_unpertubed_orbit[0], save_last_distance_unpertubed_orbit[1] ); // !!!!!!!!! this is incorrect. What would be correct would be to take the distance at final epoch - period / 2 but we did not save this distance. Since save_dca is never used after, we don't really care
	 nb_tca = nb_tca + 1;
       }
     } // end if no close approach between initial and final epoch
     else{ // if there is at a least one close approach between initial and final epoch
       if ( et_time_step_of_save_tca[0] - ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt ) > starttime ){ // if the span of the first tca starts after the initial epoch
	 if ( distance_between_two_sc( save_first_distance_unpertubed_orbit[0], save_first_distance_unpertubed_orbit[1] ) < OPTIONS->min_dist_close_approach){ // if the distance at initial epoch is smaller than OPTIONS->min_dist_close_approach then consider this at a tca (that starts at inital epoch + period / 2)
	   // need to move all previously calculated tca by one index
	   for ( ccc = nb_tca; ccc > 0; ccc-- ){
	     save_tca[ccc] = save_tca[ccc-1];
	     et_time_step_of_save_tca[ccc] = et_time_step_of_save_tca[ccc-1];
	     save_dca[ccc] = save_dca[ccc-1];
	   }
	   save_tca[0] = starttime + CONSTELLATION->collision_time_span / 2.;
	   et_time_step_of_save_tca[0] = starttime + ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt );
	   save_dca[0] = distance_between_two_sc( save_first_distance_unpertubed_orbit[0], save_first_distance_unpertubed_orbit[1] ) < OPTIONS->min_dist_close_approach; // !!!!!!!!! this is incorrect. What would be correct would be to take the distance at inital epoch + period / 2 but we did not save this distance. Since save_dca is never used after, we don't really care, except in the collision output file where dca is reported
	   nb_tca = nb_tca + 1;
	  
	 }
	  
       } // end of if the span of the first tca starts after the initial epoch

       if ( et_time_step_of_save_tca[nb_tca - 1] + ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt ) < endtime ){ // if the span of the last tca ends before the final epoch

	 if (distance_between_two_sc( save_last_distance_unpertubed_orbit[0], save_last_distance_unpertubed_orbit[1] ) < OPTIONS->min_dist_close_approach) {  // if the distance at final epoch is smaller than OPTIONS->min_dist_close_approach then consider this at a tca (that starts at final epoch - period / 2)
	   save_tca[nb_tca] = endtime - CONSTELLATION->collision_time_span / 2.;
	   et_time_step_of_save_tca[nb_tca] = endtime - ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt );
	   save_dca[nb_tca] = distance_between_two_sc( save_last_distance_unpertubed_orbit[0], save_last_distance_unpertubed_orbit[1] ) < OPTIONS->min_dist_close_approach; // !!!!!!!!! this is incorrect. What would be correct would be to take the distance at final epoch - period / 2 but we did not save this distance. Since save_dca is never used after, we don't really care
	   nb_tca = nb_tca + 1;
	 }

       } // end of if the span of the last tca ends before the final epoch

     } // end of if there is at a least one close approach between initial and final epoch

     // Update compute_collisions: if there is no close approach between the unpertubed orbits, then don't compute collisions
     if ( nb_tca  <= 0 ){
       compute_collisions = 0;
       compute_collisions_was_on_but_no_tca_found = 1;
     }

   if (nProcs > 1){
     for (ccc = 1; ccc < nProcs_that_are_gonna_run_ensembles; ccc++){

       MPI_Send(&nb_tca, 1, MPI_INT, ccc, 0, MPI_COMM_WORLD);
       MPI_Send(&compute_collisions, 1, MPI_INT, ccc, 0, MPI_COMM_WORLD);
       MPI_Send(&already_propagated_ref_sc, 1, MPI_INT, ccc, 0, MPI_COMM_WORLD);
     }
   }
     } // end of if iProc == 0
 else{
   if (nProcs > 1){
if (iProc < nProcs_that_are_gonna_run_ensembles){
     MPI_Recv(&nb_tca, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     MPI_Recv(&compute_collisions, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     MPI_Recv(&already_propagated_ref_sc, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     //            printf("Iproc %d has value %d\n", iProc, nb_tca);
 }
   }
 }

   } // end of if compute_collision
 if (compute_collisions == 1){
   done_with_tca = malloc( nb_tca * sizeof(int) );
   if ( done_with_tca == NULL ){
     print_error(iProc, "Not enough memory for done_with_tca");
   }

   // allocate meory for r/v/a eci for all sc at each time step of the time spanning each TCA
   if ( ( iDebugLevel >= 2 ) ){
     printf("--- (generate_ephemerides) Allocating big variables for collision assessment. (iProc %d)\n", iProc);
   }

   if (( iDebugLevel >= 3 ) ){
     printf("--- (generate_ephemerides) Allocating save_x_i2cg_INRTL, save_y_i2cg_INRTL, save_z_i2cg_INRTL, save_vx_i2cg_INRTL, save_vy_i2cg_INRTL, save_vz_i2cg_INRTL, save_ax_i2cg_INRTL, save_ay_i2cg_INRTL, save_az_i2cg_INRTL for collision assessment. (iProc %d)\n", iProc);
   }

   allocate_memory_r_a_v_in_span(&save_x_i2cg_INRTL, &save_y_i2cg_INRTL, &save_z_i2cg_INRTL, &save_vx_i2cg_INRTL, &save_vy_i2cg_INRTL, &save_vz_i2cg_INRTL, &save_ax_i2cg_INRTL, &save_ay_i2cg_INRTL, &save_az_i2cg_INRTL,  nb_tca, OPTIONS->nb_satellites_not_including_gps,  total_ensemble_final_with_ref, nb_time_steps_in_tca_time_span, iProc, iDebugLevel);

   if ( ( iDebugLevel >= 3 ) ){
     printf("--- (generate_ephemerides) Done allocating save_x_i2cg_INRTL, save_y_i2cg_INRTL, save_z_i2cg_INRTL, save_vx_i2cg_INRTL, save_vy_i2cg_INRTL, save_vz_i2cg_INRTL, save_ax_i2cg_INRTL, save_ay_i2cg_INRTL, save_az_i2cg_INRTL for collision assessment.(iProc %d)\n", iProc);
   }
    
   // allocate memory for nb_coll_per_step_per_iproc_in_tca nref * nb_tca * nProcs * ndtspan
   nb_coll_per_step_per_iproc_in_tca = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof(int ***) );
   if ( nb_coll_per_step_per_iproc_in_tca == NULL ){
     print_error_any_iproc(iProc, "Not enough memory for nb_coll_per_step_per_iproc_in_tca");
   }
   for (ii = 0; ii < OPTIONS->nb_satellites_not_including_gps; ii++){ // all reference sc
     nb_coll_per_step_per_iproc_in_tca[ii] = malloc( nb_tca * sizeof( int ** ) );
     for (iiitca = 0; iiitca < nb_tca; iiitca++){ // all TCA
       nb_coll_per_step_per_iproc_in_tca[ii][iiitca] = malloc( nProcs * sizeof( int *) );
       for (ccc = 0; ccc < nProcs; ccc++){ // all procs
	 nb_coll_per_step_per_iproc_in_tca[ii][iiitca][ccc] = malloc( nb_time_steps_in_tca_time_span * sizeof(int) );
	 if ( nb_coll_per_step_per_iproc_in_tca[ii][iiitca][ccc] == NULL ){
	   print_error_any_iproc(iProc, "Not enough memory for nb_coll_per_step_per_iproc_in_tca");
	 }
	 for (ppp = 0; ppp < nb_time_steps_in_tca_time_span; ppp++){
	   nb_coll_per_step_per_iproc_in_tca[ii][iiitca][ccc][ppp] = 0;
	 }
       } // end of all procs
       if ( nb_coll_per_step_per_iproc_in_tca[ii][iiitca] == NULL ){
	 print_error_any_iproc(iProc, "Not enough memory for nb_coll_per_step_per_iproc_in_tca");
       }
     } // end of all TCA

     if ( nb_coll_per_step_per_iproc_in_tca[ii] == NULL ){
       print_error_any_iproc(iProc, "Not enough memory for nb_coll_per_step_per_iproc_in_tca");
     }
   }// end of all reference sc
    // end of allocate memory for nb_coll_per_step_per_iproc_in_tca

   if (( iDebugLevel >= 2 ) ){
     printf("--- (generate_ephemerides) End of allocating big variables for collision assessment.(iProc %d)\n", iProc);
   }


   if (iProc == 0){
     // allocate memory for nb_coll_per_step_in_TCA[ii][iiitca][ppp]
     nb_coll_per_step_in_TCA = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof(int **) );
     if ( nb_coll_per_step_in_TCA == NULL ){
       print_error_any_iproc(iProc, "Not enough memory for nb_coll_per_step_in_TCA");
     }
     for (ii = 0; ii < OPTIONS->nb_satellites_not_including_gps; ii++){ // all reference sc
       nb_coll_per_step_in_TCA[ii] = malloc( nb_tca * sizeof( int * ) );
       for (iiitca = 0; iiitca < nb_tca; iiitca++){ // all TCA
	 nb_coll_per_step_in_TCA[ii][iiitca] = malloc( nb_time_steps_in_tca_time_span * sizeof( int ) );
	 if ( nb_coll_per_step_in_TCA[ii][iiitca] == NULL ){
	   print_error_any_iproc(iProc, "Not enough memory for nb_coll_per_step_in_TCA");
	 }
	 for (ppp = 0; ppp < nb_time_steps_in_tca_time_span; ppp++){
	   nb_coll_per_step_in_TCA[ii][iiitca][ppp] = 0;
	 }
       } // end of all TCA
       if ( nb_coll_per_step_in_TCA[ii] == NULL ){
	 print_error_any_iproc(iProc, "Not enough memory for nb_coll_per_step_in_TCA");
       }
     }// end of all reference sc
      // end of allocate memory for nb_coll_per_step_in_TCA[ii][iiitca][ppp]


      // allocate memory for total_nb_collisions[ii][iiitca]
     total_nb_collisions = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof(int *) );
     if ( total_nb_collisions == NULL ){
       print_error_any_iproc(iProc, "Not enough memory for total_nb_collisions");
     }
     for (ii = 0; ii < OPTIONS->nb_satellites_not_including_gps; ii++){ // all reference sc
       total_nb_collisions[ii] = malloc( nb_tca * sizeof( int ) );
       if ( total_nb_collisions[ii] == NULL ){
	 print_error_any_iproc(iProc, "Not enough memory for total_nb_collisions");
       }
       for (iiitca = 0; iiitca < nb_tca; iiitca++){
	 total_nb_collisions[ii][iiitca] = 0;
       }
     }// end of all reference sc
      // end of allocate memory for total_nb_collisions[ii][iiitca]

      // Compute max tca (furthest time of close approach). This will be the end epoch to propagate ensembles until
     max_tca = et_time_step_of_save_tca[0];

     dca_at_max_tca = save_dca[0];
   
     for (ccc = 0; ccc < nb_tca; ccc++){
       for (ppp = 1; ppp < nProcs_that_are_gonna_run_ensembles; ppp++){
	 if (nProcs > 1){
	   MPI_Send(&save_tca[ccc], 1, MPI_DOUBLE, ppp, 0, MPI_COMM_WORLD);
	   MPI_Send(&et_time_step_of_save_tca[ccc], 1, MPI_DOUBLE, ppp, 0, MPI_COMM_WORLD);
	 }
       }
       printf("\n");
       etprint(save_tca[ccc], "tca");
       ptd(save_dca[ccc], "dca");
       printf("tca %d out of %d\n", ccc+1, nb_tca);

       /* ptd( save_dca[ccc], "dca"); */
       if ( et_time_step_of_save_tca[ccc] > max_tca ){
	 max_tca =  et_time_step_of_save_tca[ccc];
	 dca_at_max_tca = save_dca[ccc];
       }
     }
     if (nProcs > 1){
       for (ppp = 1; ppp < nProcs_that_are_gonna_run_ensembles; ppp++){
	 MPI_Send(&max_tca, 1, MPI_DOUBLE, ppp, 0, MPI_COMM_WORLD);
       }
     }
   } // end of if iProc == 0
   else{
     //  pti(iProc, "iProc");
     for (ccc = 0; ccc < nb_tca; ccc++){
       if (nProcs > 1){
if (iProc < nProcs_that_are_gonna_run_ensembles){
	 MPI_Recv(&save_tca[ccc], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	 MPI_Recv(&et_time_step_of_save_tca[ccc], 1, MPI_DOUBLE, 0,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 }
       }
     }
     if (nProcs > 1){
if (iProc < nProcs_that_are_gonna_run_ensembles){
       MPI_Recv(&max_tca, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 }
     }
   }
 } // end of if collision assessment is on, then first propagate all reference sc (unpertubed orbits) to compute times of close approach (and allocate memory for variables)
  //////////// END IF COLLISION ASSESSMENT IS ON, FIRST PROPAGATE ONLY THE UNPERTURBED ORBITS



  //     exitall();
  OPTIONS->first_run_to_find_tca_before_collision_assessment = 0;
  if (iProc == 0){
    printf("\n");
  }

  if ( ( iDebugLevel >= 1 ) ){
    printf("--- (generate_ephemerides) Propagating all spacecraft (reference and ensembles spacecraft). (iProc %d)\n", iProc);
  }

  if (compute_collisions == 1){ // some reinitialization of paramters
    if ((strcmp(OPTIONS->test_omniweb_or_external_file, "swpc_mod") == 0) && (OPTIONS->swpc_need_predictions)){
      for (ii = 0; ii < OPTIONS->n_satellites; ii++) { // go through all reference satellites (and their associated ensembles deeper in the loop)
	  CONSTELLATION->aaa_mod[ii] =  OPTIONS->aaa_mod[ii];
	for (eee = 1L + iProc * OPTIONS->nb_ensemble_min_per_proc; eee< 1 + iProc * OPTIONS->nb_ensemble_min_per_proc + OPTIONS->nb_ensemble_min_per_proc ; eee++){ // go through all ensembles of reference sc ii
	  //	    printf("%d %d\n", iProc, eee);

	  CONSTELLATION->sum_sigma_for_f107_average[ii][eee] = 0;
	}
      }
    }
  }
  //  printf("XXXXXXXXXX %d %d\n",CONSTELLATION->aaa_mod[0], CONSTELLATION->aaa_mod[1]);exitf();
  /* for (ii = 0; ii < OPTIONS->n_satellites; ii++) { */
  /*   if (iProc == 0){ */
  /*   fclose( CONSTELLATION->spacecraft[ii][0].fp ); */
  /*   fclose( CONSTELLATION->spacecraft[ii][0].fpecef ); */
  /*   fclose( CONSTELLATION->spacecraft[ii][0].fpout ); */
  /*   }} */
  /*************************************************** UNCOMMENT BELOW UNTIL THIS END OF SCRIPT ************************/
  //////////////////// PROPAGATION OF ALL SC (UNPERTURBED (except if collision assessment is on) AND ENSEMBLES)

  double sma_ini = CONSTELLATION->spacecraft[0][0].OE.sma;
  // !!!!!!!!!! REMOVE LINE BELOW AND UNCOMMENT NEXT ONE
  // while ( ( fabs( CONSTELLATION->et - endtime ) > 0.01 ) && ( min_altitude_constellation > 200.0) && ( CONSTELLATION->spacecraft[0][0].OE.sma > sma_ini - 20 )){ 
    while ( ( fabs( CONSTELLATION->et - endtime ) > OPTIONS->dt / 2. ) && ( min_altitude_constellation > 100.0) ){ // propagate the constellation by dt //  while ( ( CONSTELLATION->et < endtime -0.01 ) && ( min_altitude_constellation > 200.0) ){ // propagate the constellation by dt
    /* printf("\n"); */
    /* pti(iProc, "iProc"); */
      //                           etprint(CONSTELLATION->et, "t");


    // Time to stop the propgation
    min_end_time = endtime;
    //    etprint(min_end_time, "end");
    if ( ( compute_collisions == 1 ) ){//2020-11-13 i don't want to worry about dt_pos_neg (backward propagation) for collision avoidance
      if ( max_tca + ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt ) < endtime ){
	min_end_time =  max_tca + ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt );
      }
    }

    // Print progress to screen
    if (iProc == 0){
                                                                                     print_progress( min_end_time, CONSTELLATION->et , starttime, iProc, OPTIONS->nb_gps )  ;
    }

    if ( ( compute_collisions == 1 )  ){       // start collision assessment when the secondary sc time enters the span of time around TCA. !!!!!!!!! CONSTELLATION->et  is temporary and assumes the reference sc have the same epoch start. If it's not the case then this needs to be changed (below and at other locations in the code) //2020-11-13 i don't want to worry about dt_pos_neg (backward propagation) for collision avoidance
      itca = -1; // if itca is different from -1 then it means that the time is in the time spanning TCA (unpertubed)
      time_step_of_tca = -1;
      for ( ccc = start_itca; ccc < nb_tca; ccc++ ){
    	if ( ( CONSTELLATION->et + OPTIONS->dt >= et_time_step_of_save_tca[ccc] - ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt ) )  ){
    	  itca = ccc;
	  et_current_tca = et_time_step_of_save_tca[ccc];
    	  et_start_of_span = et_time_step_of_save_tca[ccc] - ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt );
    	  et_end_of_span = et_time_step_of_save_tca[ccc] + ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt );
	  //	  etprint(et_end_of_span, "et_end_of_span");
    	  time_step_of_tca = (int)(nb_time_steps_in_tca_time_span / 2);
    	  done_with_tca[itca] = 0;
	  //	  print_test();
    	}
      }
    } // end of  start collision assessment when the secondary sc time enters the span of time around TCA. !!!!!!!!! CONSTELLATION->et  is temporary and assumes the reference sc have the same epoch start. If it's not the case then this needs to be changed (below and at other locations in the code)

    for (ii = 0; ii < OPTIONS->n_satellites; ii++) { // go through all reference satellites (and their associated ensembles deeper in the loop)
      //      if (iProc == 0){
      if ( start_ensemble[ii] == 0){ // if this iProc runs main sc ii
	if ( already_propagated_ref_sc == 0 ){ // i unpertubed orbits (eee = 0) would have already been propagated so don't propagate them again here
	  // main spacecraft only (eee = 0 -> "[ii][0]")
	  previous_lat = CONSTELLATION->spacecraft[ii][0].GEODETIC.latitude*RAD2DEG;

	  propagate_spacecraft( &CONSTELLATION->spacecraft[ii][0], PARAMS, starttime, OPTIONS->et_oldest_tle_epoch, &density, GROUND_STATION, OPTIONS, CONSTELLATION, iProc, iDebugLevel, start_ensemble, array_sc );

	  //	  printf("%f\n", CONSTELLATION->spacecraft[ii][0].INTEGRATOR.surface[0].Cd );
	  if (CONSTELLATION->spacecraft[ii][0].GEODETIC.altitude <  min_altitude_constellation) {
	    min_altitude_constellation = CONSTELLATION->spacecraft[ii][0].GEODETIC.altitude;
	  }
	} // end of if collision assessment is on (compute_collisions = 1), then the unpertubed orbits (eee = 0) have already been propagated so don't propagate them again her
      } // end of if this iProc runs main sc ii
      //} // end of iproc == 0
      // ensembles spacecraft (but not for gps satellites)



      if (ii < OPTIONS->n_satellites - OPTIONS-> nb_gps){ // if the reference satellite is not a GPS, we propagate the associated ensembles
	if ( OPTIONS->nb_ensembles_min > 0 ) { // if ensembles are computed
	  if ( array_sc[1] > 0 )  { // if this iProc runs ensembles (otherwise array_sc[1] = - 1)
	  for (eee = 1L + iProc * OPTIONS->nb_ensemble_min_per_proc; eee< 1 + iProc * OPTIONS->nb_ensemble_min_per_proc + OPTIONS->nb_ensemble_min_per_proc ; eee++){ // go through all ensembles of reference sc ii
	    //	    printf("%d %d\n", iProc, eee);


	    propagate_spacecraft( &CONSTELLATION->spacecraft[ii][eee], PARAMS, starttime, OPTIONS->et_oldest_tle_epoch, &density, GROUND_STATION, OPTIONS, CONSTELLATION, iProc  , iDebugLevel, start_ensemble, array_sc );

	    // min altitude
	    /* 	    if (CONSTELLATION->spacecraft[ii][eee].GEODETIC.altitude < min_altitude_constellation){ */
	    /* 	      min_altitude_constellation =  CONSTELLATION->spacecraft[ii][eee].GEODETIC.altitude; */
	    /* 	    } */

	    ispan_constellation = (int)( fabs( CONSTELLATION->et + OPTIONS->dt_pos_neg - et_start_of_span ) / OPTIONS->dt );
	    /* if (iProc == 1){ */
	    /* etprint(CONSTELLATION->et + OPTIONS->dt, "out" ); */
	    /* printf("%d %d\n", ispan_constellation,nb_time_steps_in_tca_time_span-1); */
	    /* //      etprint(et_start_of_span, "start span"); */
	    /* //    	      etprint(et_end_of_span, "end span");  */
	    /* } */
 	    if ( ( compute_collisions == 1 )  && ( itca != -1 ) && ( ispan_constellation < nb_time_steps_in_tca_time_span ) ){ // if the time is in the time span around TCA (!!!!!!!! this temporary and assumes the reference sc have the same epoch start. If it's not the case then this needs to be changed (below and at other locations in the code))
 	      //CONSTELLATION->spacecraft[ii][eee].ispan = (int)( ( CONSTELLATION->et + OPTIONS->dt - et_start_of_span ) / OPTIONS->dt ); // CONSTELLATION->et + OPTIONS->dt corresponds to the new time after propagation (recall that CONSTELLATION->et us pidated only at the end of the loop so that's why we need to add + OPTIONS->dt here)
	      /* if (iProc == 1){ */
	      /* etprint(CONSTELLATION->et + OPTIONS->dt, "in" ); */
	      /* pti(ispan_constellation, "ispan_constellation"); */
	      /* } */

 	      //	      printf("%f\n", save_x_i2cg_INRTL[0][0][0][0]);
 	      // // Save positions of all ensembles
 	      save_x_i2cg_INRTL[ii][eee][ispan_constellation] = CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[0];
 	      save_y_i2cg_INRTL[ii][eee][ispan_constellation] = CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[1];
 	      save_z_i2cg_INRTL[ii][eee][ispan_constellation] = CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[2];

 	      save_vx_i2cg_INRTL[ii][eee][ispan_constellation] = CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[0];
 	      save_vy_i2cg_INRTL[ii][eee][ispan_constellation] = CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[1];
 	      save_vz_i2cg_INRTL[ii][eee][ispan_constellation] = CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[2];

 	      save_ax_i2cg_INRTL[ii][eee][ispan_constellation] = CONSTELLATION->spacecraft[ii][eee].a_i2cg_INRTL[0];
 	      save_ay_i2cg_INRTL[ii][eee][ispan_constellation] = CONSTELLATION->spacecraft[ii][eee].a_i2cg_INRTL[1];
 	      save_az_i2cg_INRTL[ii][eee][ispan_constellation] = CONSTELLATION->spacecraft[ii][eee].a_i2cg_INRTL[2];


	      if ( ii == 0 ){//  send the position/velocity/acceleration of a primary sc computed by a iProc to all other iProc  !!!!!!!!!! FOR NOW WORKS ONLY IF TWO REFERENCE SATELLIES ONLY  ("TWO" because see next if condition). Note that main sc 0 is always run by iProc 0 and never by iProc > 0.

		if (nProcs > 1){

		  for (ccc = 0; ccc < nProcs_that_are_gonna_run_ensembles; ccc++){ // send the position/velocity/acceleration computed by a iProc to all other iProc
		    if ( ccc != iProc ){
		      MPI_Send(&save_x_i2cg_INRTL[ii][eee][ispan_constellation], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD);
		      MPI_Send(&save_y_i2cg_INRTL[ii][eee][ispan_constellation], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD);
		      MPI_Send(&save_z_i2cg_INRTL[ii][eee][ispan_constellation], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD);

		      MPI_Send(&save_vx_i2cg_INRTL[ii][eee][ispan_constellation], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD);
		      MPI_Send(&save_vy_i2cg_INRTL[ii][eee][ispan_constellation], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD);
		      MPI_Send(&save_vz_i2cg_INRTL[ii][eee][ispan_constellation], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD);

		      MPI_Send(&save_ax_i2cg_INRTL[ii][eee][ispan_constellation], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD);
		      MPI_Send(&save_ay_i2cg_INRTL[ii][eee][ispan_constellation], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD);
		      MPI_Send(&save_az_i2cg_INRTL[ii][eee][ispan_constellation], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD);

		    }
		  } // end of send the position/velocity/acceleration of a primary computed by a iProc to all other iProc
		} // end of if nProcs > 1
		//			MPI_Barrier(MPI_COMM_WORLD);
		for (ccc = 0; ccc < nProcs; ccc++){ 	  // receive the positions by iProc sent by all other iProc
		  if ( ccc != iProc ){

		    eee_all_other_iproc = ccc * OPTIONS->nb_ensemble_min_per_proc +  eee - iProc * OPTIONS->nb_ensemble_min_per_proc;
		    if (nProcs > 1){
		      if (iProc < nProcs_that_are_gonna_run_ensembles){
		      MPI_Recv(&save_x_i2cg_INRTL[ii][eee_all_other_iproc][ispan_constellation], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		      //            printf("save_x_i2cg_INRTL[%d][%d][%d][%d]: %e (iProc: %d)\n", itca, ii, eee_all_other_iproc, ispan_constellation, save_x_i2cg_INRTL[ii][eee_all_other_iproc][ispan_constellation], iProc);
		      MPI_Recv(&save_y_i2cg_INRTL[ii][eee_all_other_iproc][ispan_constellation], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		      MPI_Recv(&save_z_i2cg_INRTL[ii][eee_all_other_iproc][ispan_constellation], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		      MPI_Recv(&save_vx_i2cg_INRTL[ii][eee_all_other_iproc][ispan_constellation], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		      MPI_Recv(&save_vy_i2cg_INRTL[ii][eee_all_other_iproc][ispan_constellation], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		      MPI_Recv(&save_vz_i2cg_INRTL[ii][eee_all_other_iproc][ispan_constellation], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		      MPI_Recv(&save_ax_i2cg_INRTL[ii][eee_all_other_iproc][ispan_constellation], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		      MPI_Recv(&save_ay_i2cg_INRTL[ii][eee_all_other_iproc][ispan_constellation], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		      MPI_Recv(&save_az_i2cg_INRTL[ii][eee_all_other_iproc][ispan_constellation], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		      }
		    }
		  }
		} // end of receive the positions by iProc sent by all other iProc
	      }  // end of if ii == 0
 	    }	   //  	     end of if the time is in the time span around TCA
	    

	    //	    	    test_print_iproc( iProc, "C");

	    if (ii == 1){ // only get in this loop if secondary sc
	      if ( ( compute_collisions == 1 )  && ( itca != -1 ) && ( ispan_constellation == nb_time_steps_in_tca_time_span - 1 ) ){// eee just got out of the if condition "if ( ( compute_collisions == 1 )  && ( itca != -1 ) && ( ispan_constellation < nb_time_steps_in_tca_time_span ) )". This means that the secondary spacecraft eee waas propagated all they way to TCA + span/2 and is now ready to be assessed for collisions with all eee primary (called eee_prim in this if condition)


		if ( ( iDebugLevel >= 5 ) ){
		  printf("----- (generate_ephemerides) Just got out of time spanning TCA %d and now computing collisions for this time span (iProc %d)\n", itca+1, iProc);
		}

		for (eee_prim = 1; eee_prim< 1 + nProcs * OPTIONS->nb_ensemble_min_per_proc ; eee_prim++){ // go over all primary ensemble sc (eee_prim) to look for a collision with secondary ensemble sc eee
		  /* //!!!!!!!!!! */
		  /* int istep = 0; */
		  /* distance_temp = 1e6; */
		  /* while ( (istep < nb_time_steps_in_tca_time_span) && (distance_temp > OPTIONS->min_dist_collision) ){ */
		  /*     distance_temp = sqrt( (save_x_i2cg_INRTL[1][eee][istep] - save_x_i2cg_INRTL[0][eee_prim][istep])*(save_x_i2cg_INRTL[1][eee][istep] - save_x_i2cg_INRTL[0][eee_prim][istep]) + */
		  /* 			  (save_y_i2cg_INRTL[1][eee][istep] - save_y_i2cg_INRTL[0][eee_prim][istep])*(save_y_i2cg_INRTL[1][eee][istep] - save_y_i2cg_INRTL[0][eee_prim][istep]) + */
		  /* 			  (save_z_i2cg_INRTL[1][eee][istep] - save_z_i2cg_INRTL[0][eee_prim][istep])*(save_z_i2cg_INRTL[1][eee][istep] - save_z_i2cg_INRTL[0][eee_prim][istep]) ); */

		  /*   if ( distance_temp < OPTIONS->min_dist_collision){ */
		  /*     ncol[iProc] = ncol[iProc] + 1; */
		  /*     //	    printf("ncol[%d]: %d (final step: %d)\n",iProc, ncol[iProc], istep); */
		  /*   } */
		  /*   istep = istep + 1; */
		  /* } */
		  /* //!!!!!!!!!!!! */
		  //	printf("eee: %d | eee_prim: %d (iProc :%d)\n", eee, eee_prim, iProc);
		  compute_collision_between_one_secondary_and_all_primary(save_x_i2cg_INRTL[1][eee], save_y_i2cg_INRTL[1][eee], save_z_i2cg_INRTL[1][eee], save_vx_i2cg_INRTL[1][eee], save_vy_i2cg_INRTL[1][eee], save_vz_i2cg_INRTL[1][eee], save_ax_i2cg_INRTL[1][eee], save_ay_i2cg_INRTL[1][eee], save_az_i2cg_INRTL[1][eee],
									  save_x_i2cg_INRTL[0][eee_prim], save_y_i2cg_INRTL[0][eee_prim], save_z_i2cg_INRTL[0][eee_prim], save_vx_i2cg_INRTL[0][eee_prim], save_vy_i2cg_INRTL[0][eee_prim], save_vz_i2cg_INRTL[0][eee_prim], save_ax_i2cg_INRTL[0][eee_prim], save_ay_i2cg_INRTL[0][eee_prim], save_az_i2cg_INRTL[0][eee_prim],
									  OPTIONS, iProc,  nb_coll_per_step_per_iproc_in_tca[1][itca][iProc], et_time_step_of_save_tca, nb_time_steps_in_tca_time_span, itca, eee_prim, eee, tca_file, dca_file, sample_file, OPTIONS->write_collision_files, &eee_prim_that_collide );// !!!!!!!!!! FOR NOW WORKS ONLY IF TWO REFERENCE SATELLIES ONLY
		  //		  printf("out %d\n", eee_prim_that_collide);

		} // end of go over all primary ensemble sc (eee_prim) to look for a collision with secondary ensemble sc eee
		if ( ( iDebugLevel >= 5 ) ){
		  printf("----- (generate_ephemerides) DONE just got out of time spanning TCA %d and now computing collisions for this time span (iProc %d)\n", itca+1, iProc);
		}

	      }
	    } // end of if ii == 1

	  } // end of go through all ensembles of reference sc ii
	  } // end of if this iProc runs ensembles (otherwise array_sc[1] = - 1)
	  //	  test_print_iproc( iProc, "D");
	} // end of if ensembles are computed
      } // end of if reference satellite is not GPS, we propagate the associated ensembles
      //      test_print_iproc( iProc, "E");
      // WRITE OUTPUT

      /* if ( ( iDebugLevel >= 4 ) ){ */
      /*   printf("----- (generate_ephemerides) Writing the output (iProc %d).\n",  iProc); */
      /* } */

      if ( ( already_propagated_ref_sc == 1 )   ){ //2020-11-13 i don't want to worry about dt_pos_neg (backward propagation) for ensemble propagation
	if ( output_only_at_tca != 1 ){
	  if ( (array_sc[start_ensemble[ii]] >= 0) && ((CONSTELLATION->spacecraft[ii][1 + iProc * OPTIONS->nb_ensemble_min_per_proc].et - twrite) >= OPTIONS->dt - 0.01) ){// && ( itca != -1 ) ) { // (array_sc[start_ensemble[ii]] >= 0): if this iProc runs main sc ii (start_ensemble[ii] = 0, array_sc[0] = 0) or this iProc does not run main sc ii but runs ensembles for main sc ii (start_ensemble[ii] = 1, array_sc[1] > 0) (if none of these two, so if this iProc does not run main sc ii and that this iProc does not run ensembles for main sc ii, then we don't want this iProc to write anything for main sc ii (start_ensemble[ii] = 1 and array_sc[1] = -1).

	    if ( ( fmod( CONSTELLATION->spacecraft[ii][1 + iProc * OPTIONS->nb_ensemble_min_per_proc].et - starttime, OPTIONS->dt_output ) < OPTIONS->dt / 2.) || ( fabs( fmod( CONSTELLATION->spacecraft[ii][1 + iProc * OPTIONS->nb_ensemble_min_per_proc].et - starttime, OPTIONS->dt_output ) - OPTIONS->dt_output ) < OPTIONS->dt / 2. ) || ( CONSTELLATION->spacecraft[ii][1 + iProc * OPTIONS->nb_ensemble_min_per_proc].et > min_end_time - 0.01) )  {
	      write_output( CONSTELLATION->spacecraft[ii] , 0, choose_tle_to_initialise_orbit,ii, OPTIONS->n_satellites,OPTIONS->nb_gps,  OPTIONS->nb_ensembles_min, OPTIONS->nb_ensemble_min_per_proc,  iProc,  OPTIONS->nb_ensembles_output_files, OPTIONS->filename_output_ensemble, previous_lat,OPTIONS,PARAMS->EARTH.earth_fixed_frame,0,1, ( CONSTELLATION->et + OPTIONS->dt_pos_neg ), nProcs, iDebugLevel, compute_collisions, start_ensemble, array_sc, CONSTELLATION, PARAMS);
	    }
	  }
	  if ( ( ii == (OPTIONS->n_satellites - 1) ) && (compute_collisions == 0) ) {
	    twrite = CONSTELLATION->et + OPTIONS->dt;
	  }
	  if ( ( ii == (OPTIONS->n_satellites - 1) ) && (compute_collisions == 1) ) {
	    twrite = CONSTELLATION->et + OPTIONS->dt;
	  }

	  /* etprint(et_current_tca, "tca"); */
	  /* etprint(CONSTELLATION->spacecraft[ii][1 + iProc * OPTIONS->nb_ensemble_min_per_proc].et, "sc"); */
	}
	// to output only at tca of unperturbed orbit !!!!!!!!!!!!!!!!!
	else{
	  if (itca != -1){
	    // if ( fabs(CONSTELLATION->spacecraft[ii][1 + iProc * OPTIONS->nb_ensemble_min_per_proc].et - et_current_tca) <= 0.01) {// && ( itca != -1 ) ) { // write output
	      if ( ( CONSTELLATION->spacecraft[ii][1 + iProc * OPTIONS->nb_ensemble_min_per_proc].et >= et_current_tca - OPTIONS->dt ) && ( CONSTELLATION->spacecraft[ii][1 + iProc * OPTIONS->nb_ensemble_min_per_proc].et <= et_current_tca + OPTIONS->dt  ) ){ // !!!!!!!!! REMOVE THIS IF AND UNCOMMENT THE ONE RIGHT ABOVE                                     \
	      //		      print_test();
	      write_output( CONSTELLATION->spacecraft[ii] , 0, choose_tle_to_initialise_orbit,ii, OPTIONS->n_satellites,OPTIONS->nb_gps,  OPTIONS->nb_ensembles_min, OPTIONS->nb_ensemble_min_per_proc,  iProc,  OPTIONS->nb_ensembles_output_files, OPTIONS->filename_output_ensemble, previous_lat,OPTIONS,PARAMS->EARTH.earth_fixed_frame,0,1, ( CONSTELLATION->et + OPTIONS->dt_pos_neg ), nProcs, iDebugLevel, compute_collisions, start_ensemble, array_sc, CONSTELLATION, PARAMS);
	    }
	    // end of to output only at tca of unperturbed orbit !!!!!!!!!!!
	  }
	}
    
      } // end of if ( ( already_propagated_ref_sc == 1 )
      else{
	//  old (when no parallel programming for constellation)      if ( ((CONSTELLATION->spacecraft[ii][start_ensemble + iProc * OPTIONS->nb_ensemble_min_per_proc].et - twrite) >= OPTIONS->dt - 0.01) ){// && ( itca != -1 ) ) { // write output

/* 	  if (iProc > 1){ */
/* 	    printf("iProc %d has array_sc[start_ensemble[%d]] = %d \n", iProc, ii, array_sc[start_ensemble[ii]] ); */
/* 	  } */

	if ( (array_sc[start_ensemble[ii]] >= 0) && (fabs(CONSTELLATION->spacecraft[ii][array_sc[start_ensemble[ii]]].et - twrite) >= OPTIONS->dt - 0.01) ){// && ( itca != -1 ) ) { // (array_sc[start_ensemble[ii]] >= 0): if this iProc runs main sc ii (start_ensemble[ii] = 0, array_sc[0] = 0) or this iProc does not run main sc ii but runs ensembles for main sc ii (start_ensemble[ii] = 1, array_sc[1] > 0) (if none of these two, so if this iProc does not run main sc ii and that this iProc does not run ensembles for main sc ii, then we don't want this iProc to write anything for main sc ii (start_ensemble[ii] = 1 and array_sc[1] = -1). array_sc[start_ensemble[ii]] in the second condition (after '&&') represents main sc ii if this iProc runs main sc ii (because array_sc[start_ensemble[ii]] = 0), or it represents the first ensemble run by this iProc if this iProc does not run main sc ii
	  //	    test_print("A")
	  //    etprint(CONSTELLATION->et, "1");
	  /* 	    etprint(starttime, "start"); */
	  /* 	    etprint(CONSTELLATION->spacecraft[ii][array_sc[start_ensemble[ii]]].et, "et"); */
	  //	    printf("%e %e %e %f %d || %f\n", CONSTELLATION->spacecraft[ii][array_sc[start_ensemble[ii]]].et-starttime, OPTIONS->dt_output, starttime, fmod( CONSTELLATION->spacecraft[ii][array_sc[start_ensemble[ii]]].et - starttime, OPTIONS->dt_output ) , ( fmod( CONSTELLATION->spacecraft[ii][array_sc[start_ensemble[ii]]].et - starttime, OPTIONS->dt_output ) < 0.01), fabs( fmod( CONSTELLATION->spacecraft[ii][array_sc[start_ensemble[ii]]].et - starttime, OPTIONS->dt_output ) - OPTIONS->dt_output ) );
	  
	  //	  printf("%f\n", CONSTELLATION->spacecraft[ii][array_sc[start_ensemble[ii]]].et);
	  if (OPTIONS->dt_pos_neg >= 0){
	    if ( ( fmod( CONSTELLATION->spacecraft[ii][array_sc[start_ensemble[ii]]].et - starttime, OPTIONS->dt_output ) < OPTIONS->dt / 2.) || ( fabs( fmod( CONSTELLATION->spacecraft[ii][array_sc[start_ensemble[ii]]].et - starttime, OPTIONS->dt_output ) - OPTIONS->dt_output ) < OPTIONS->dt / 2.) || ( CONSTELLATION->spacecraft[ii][array_sc[start_ensemble[ii]]].et > min_end_time - 0.01) )  {
	    //    test_print("B");
/* 	  if (iProc > 4){ */
/* 	    printf("%d\n", iProc); */
/* 	  } */

	      write_output( CONSTELLATION->spacecraft[ii] , 0, choose_tle_to_initialise_orbit,ii, OPTIONS->n_satellites,OPTIONS->nb_gps,  OPTIONS->nb_ensembles_min, OPTIONS->nb_ensemble_min_per_proc,  iProc,  OPTIONS->nb_ensembles_output_files, OPTIONS->filename_output_ensemble, previous_lat,OPTIONS,PARAMS->EARTH.earth_fixed_frame,1,1, ( CONSTELLATION->et + OPTIONS->dt_pos_neg ), nProcs, iDebugLevel, compute_collisions, start_ensemble, array_sc, CONSTELLATION, PARAMS);
      	    et2utc_c(CONSTELLATION->spacecraft[ii][0].et, "ISOC" ,6 ,255 , times_att);
/* 	     if (write_density == 1){ */
/* 	       fprintf(density_file,"%s %e %f %f %e\n", times_att, density, CONSTELLATION->spacecraft[ii][0].INTEGRATOR.Ta, CONSTELLATION->spacecraft[ii][0].INTEGRATOR.cd_tot_norm, CONSTELLATION->spacecraft[ii][0].INTEGRATOR.A_ref_tot/1000000.); // in output A_ref_tot in km^2 */
/* 	       //	       printf("%e %e\n", density, CONSTELLATION->spacecraft[ii][0].density_here); */
/* 	     } */
/* char  time_temp_str[256]; */
/*   double time_temp; */
/*   strcpy(time_temp_str, "2017-05-04T00:00:40.999261"); */
/*   str2et_c(time_temp_str, &time_temp); */
/*   if (CONSTELLATION->spacecraft[ii][0].et >= time_temp){ */
/*     MPI_Finalize(); exit(0); */
/*   } */

      	  }
	  }
	  else{
	    if ( ( fmod(fabs( CONSTELLATION->spacecraft[ii][array_sc[start_ensemble[ii]]].et - starttime), OPTIONS->dt_output ) < OPTIONS->dt / 2.) || ( fabs( fmod( CONSTELLATION->spacecraft[ii][array_sc[start_ensemble[ii]]].et - starttime, OPTIONS->dt_output ) - OPTIONS->dt_output ) < OPTIONS->dt / 2.) || ( CONSTELLATION->spacecraft[ii][array_sc[start_ensemble[ii]]].et < min_end_time + 0.01) )  {
	    //    test_print("B");
/* 	  if (iProc > 4){ */
/* 	    printf("%d\n", iProc); */
/* 	  } */

	      write_output( CONSTELLATION->spacecraft[ii] , 0, choose_tle_to_initialise_orbit,ii, OPTIONS->n_satellites,OPTIONS->nb_gps,  OPTIONS->nb_ensembles_min, OPTIONS->nb_ensemble_min_per_proc,  iProc,  OPTIONS->nb_ensembles_output_files, OPTIONS->filename_output_ensemble, previous_lat,OPTIONS,PARAMS->EARTH.earth_fixed_frame,1,1, ( CONSTELLATION->et + OPTIONS->dt_pos_neg ), nProcs, iDebugLevel, compute_collisions, start_ensemble, array_sc, CONSTELLATION, PARAMS);
      	    et2utc_c(CONSTELLATION->spacecraft[ii][0].et, "ISOC" ,6 ,255 , times_att);
/* 	     if (write_density == 1){ */
/* 	       fprintf(density_file,"%s %e %f %f %e\n", times_att, density, CONSTELLATION->spacecraft[ii][0].INTEGRATOR.Ta, CONSTELLATION->spacecraft[ii][0].INTEGRATOR.cd_tot_norm, CONSTELLATION->spacecraft[ii][0].INTEGRATOR.A_ref_tot/1000000.); // in output A_ref_tot in km^2 */
/* 	       //	       printf("%e %e\n", density, CONSTELLATION->spacecraft[ii][0].density_here); */
/* 	     } */
/* char  time_temp_str[256]; */
/*   double time_temp; */
/*   strcpy(time_temp_str, "2017-05-04T00:00:40.999261"); */
/*   str2et_c(time_temp_str, &time_temp); */
/*   if (CONSTELLATION->spacecraft[ii][0].et >= time_temp){ */
/*     MPI_Finalize(); exit(0); */
/*   } */

      	  }

	  }
	}

      	
      	if ( ( ii == (OPTIONS->n_satellites - 1) ) && (compute_collisions == 0) ) {
      	  twrite = CONSTELLATION->et + OPTIONS->dt_pos_neg;
      	}
      	if ( ( ii == (OPTIONS->n_satellites - 1) ) && (compute_collisions == 1) ) {
      	  twrite = CONSTELLATION->et + OPTIONS->dt_pos_neg;
      	}

      } // end of else (so end of if already_propagated_ref_sc != 1)

      /* if ( ( iDebugLevel >= 4 ) ){ */
      /*   printf("----- (generate_ephemerides) Done writing the output (iProc %d).\n",  iProc); */
      /* } */

      
      //          test_print_iproc( iProc, "F");
    } // go through all reference satellites (and their associated ensembles deeper in the loop)  (END OF LOOP OVER ALL MAIN SATELLITES)

    if ( ( compute_collisions == 1 )){
      if ( ( CONSTELLATION->et + OPTIONS->dt ) > max_tca + ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt )   ){// if collision assessment is on, then the ensembles are propagated until the time of closest approach (+ a certain time, ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt ) )
      //        test_print_iproc( iProc, "B");
      break;
      }
    }
    /// !!!!!!!!!!!!!!!!!!!!!! COMMENT THIS BLOCK!!!!!!!!!!!!!!!!!!!!!
    /* if (( CONSTELLATION->et + OPTIONS->dt ) > remove_this_var) { */
    /*   break; */
    /* } */
    /// !!!!!!!!!!!!!!!!!!!!!! END OF COMMENT THIS BLOCK!!!!!!!!!!!!!!!!!!!!!
    // UPDATE THE CONSTELLATION TIME
    
    CONSTELLATION->et = CONSTELLATION->et + OPTIONS->dt_pos_neg;
  } // end of propagate the constellation by dt  (END OF WHILE PROPAGATION OF ALL SATELLITES AND ENSEMBLES)
  //////////////////// END OF PROPAGATION



  if (( iDebugLevel >= 1 ) ){
    printf("--- (generate_ephemerides) Done propagating all spacecraft (reference and ensembles spacecraft). (iProc %d)\n", iProc);
  }

  if (iProc == 0){
    printf("\n");
  }

  //  exitall();
     /* test_print_iproc(iProc, "A"); */

  
  // //  COMPUTE COLLISION FOR EACH TCA ///////////////////////////////////
    /* !!!!!!!! */
    /* int *ncol; */
    /* ncol = malloc(nProcs * sizeof(int)); */
    /* ncol[iProc] = 0; */
    /* double distance_temp; */
    /* !!!!!!! */
  if ( ( compute_collisions == 1 ) ){ // if collision assessment is on



    /* for (iiitca = start_itca; iiitca < nb_tca; iiitca++){ // for each tca, collision assessment */
    
    /*   for (eee_sec = 1L + iProc * OPTIONS->nb_ensemble_min_per_proc; eee_sec< 1 + iProc * OPTIONS->nb_ensemble_min_per_proc + OPTIONS->nb_ensemble_min_per_proc ; eee_sec++){ // go through all ensembles of secondary sc */
      
    /*     //                 print_progress_collision(eee_sec, iProc, OPTIONS->nb_ensemble_min_per_proc, nb_tca); */
    /*     for (eee_prim = 1; eee_prim< 1 + nProcs * OPTIONS->nb_ensemble_min_per_proc ; eee_prim++){ // go over all primary ensemble sc (eee_prim) to look for a collision with secondary ensemble sc eee */
    /* 	/\* //!!!!!!!!!! *\/ */
    /* 	/\* int istep = 0; *\/ */
    /* 	/\* distance_temp = 1e6; *\/ */
    /* 	/\* while ( (istep < nb_time_steps_in_tca_time_span) && (distance_temp > OPTIONS->min_dist_collision) ){ *\/ */
    /* 	/\*     distance_temp = sqrt( (save_x_i2cg_INRTL[iiitca][1][eee_sec][istep] - save_x_i2cg_INRTL[iiitca][0][eee_prim][istep])*(save_x_i2cg_INRTL[iiitca][1][eee_sec][istep] - save_x_i2cg_INRTL[iiitca][0][eee_prim][istep]) + *\/ */
    /* 	/\* 			  (save_y_i2cg_INRTL[iiitca][1][eee_sec][istep] - save_y_i2cg_INRTL[iiitca][0][eee_prim][istep])*(save_y_i2cg_INRTL[iiitca][1][eee_sec][istep] - save_y_i2cg_INRTL[iiitca][0][eee_prim][istep]) + *\/ */
    /* 	/\* 			  (save_z_i2cg_INRTL[iiitca][1][eee_sec][istep] - save_z_i2cg_INRTL[iiitca][0][eee_prim][istep])*(save_z_i2cg_INRTL[iiitca][1][eee_sec][istep] - save_z_i2cg_INRTL[iiitca][0][eee_prim][istep]) ); *\/ */

    /* 	/\*   if ( distance_temp < OPTIONS->min_dist_collision){ *\/ */
    /* 	/\*     ncol[iProc] = ncol[iProc] + 1; *\/ */
    /* 	/\*     //	    printf("ncol[%d]: %d (final step: %d)\n",iProc, ncol[iProc], istep); *\/ */
    /* 	/\*   } *\/ */
    /* 	/\*   istep = istep + 1; *\/ */
    /* 	/\* } *\/ */
    /* 	/\* //!!!!!!!!!!!! *\/ */
    /* 	//	printf("eee_sec: %d | eee_prim: %d (iProc :%d)\n", eee_sec, eee_prim, iProc); */
    /*     compute_collision_between_one_secondary_and_all_primary(save_x_i2cg_INRTL[iiitca][1][eee_sec], save_y_i2cg_INRTL[iiitca][1][eee_sec], save_z_i2cg_INRTL[iiitca][1][eee_sec], save_vx_i2cg_INRTL[iiitca][1][eee_sec], save_vy_i2cg_INRTL[iiitca][1][eee_sec], save_vz_i2cg_INRTL[iiitca][1][eee_sec], save_ax_i2cg_INRTL[iiitca][1][eee_sec], save_ay_i2cg_INRTL[iiitca][1][eee_sec], save_az_i2cg_INRTL[iiitca][1][eee_sec], */
    /* 							      save_x_i2cg_INRTL[iiitca][0][eee_prim], save_y_i2cg_INRTL[iiitca][0][eee_prim], save_z_i2cg_INRTL[iiitca][0][eee_prim], save_vx_i2cg_INRTL[iiitca][0][eee_prim], save_vy_i2cg_INRTL[iiitca][0][eee_prim], save_vz_i2cg_INRTL[iiitca][0][eee_prim], save_ax_i2cg_INRTL[iiitca][0][eee_prim], save_ay_i2cg_INRTL[iiitca][0][eee_prim], save_az_i2cg_INRTL[iiitca][0][eee_prim], */
    /* 							      OPTIONS, iProc,  nb_coll_per_step_per_iproc_in_tca[1][iiitca][iProc], et_time_step_of_save_tca, nb_time_steps_in_tca_time_span, iiitca, eee_prim, eee_sec );// !!!!!!!!!! FOR NOW WORKS ONLY IF TWO REFERENCE SATELLIES ONLY */

    /*   } // end of go over all primary ensemble sc (eee_prim) to look for a collision with secondary ensemble sc eee */
    /*   } // end of go through all ensembles of secondary sc */
    /* }// end of for each tca, collision assessment */

    //////////////////////// END OF COMPUTE COLLISION FOR EACH TCA ////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////

    //test_print_iproc(iProc, "V");
    /* MPI_Finalize(); */
    /* exit(0); */
    // !!!!!!!!!!!
    /* MPI_Barrier(MPI_COMM_WORLD); */

    /* if ( iProc !=0 ){ */
    /*   MPI_Send(&ncol[iProc], 1, MPI_INT, 0, 0, MPI_COMM_WORLD); */
    /* } */
    /* if (iProc == 0){ */
    /*   int ncol_total = ncol[0]; */
    /*     for (ccc = 1; ccc < nProcs; ccc++){ // send the position/velocity/acceleration computed by a iProc to all other iProc */
    /* 	MPI_Recv(&ncol[ccc], 1, MPI_INT, ccc,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); */
    /* 	ncol_total = ncol_total + ncol[ccc]; */
    /*     } */
    /*     printf("<<<<<<<<< TOTAL COLLISIONS BUT ONLY AT STEPS: %d >>>>>>>>>\n", ncol_total); */
    /* } */

    // !!!!!!!!!!!
    // COLLISION

    if ( ( iDebugLevel >= 4 ) ){
      printf("----- (generate_ephemerides) Sending nb_coll_per_step_per_iproc_in_tca to main node... (iProc %d)\n", iProc);
    }

    ii = 1; // !!!!!!!!! FOR NOW WORKS ONLY IF TWO REFERENCE SATELLIES ONLY
    for ( iiitca = start_itca; iiitca < nb_tca; iiitca++ ){
      for ( ppp = 0; ppp < nb_time_steps_in_tca_time_span; ppp++){
	if (iProc != 0){ 	      // send step_coll to main node
	  if (nProcs > 1){
	    MPI_Send(&nb_coll_per_step_per_iproc_in_tca[ii][iiitca][iProc][ppp], 1, MPI_INT, 0, 2 * ii *   + iiitca * ( nb_time_steps_in_tca_time_span ) *  ( nProcs * OPTIONS->nb_ensemble_min_per_proc + 1 + 1 ) + ppp * ( 1 + nProcs * OPTIONS->nb_ensemble_min_per_proc ) + iProc , MPI_COMM_WORLD);
	  }
	} // end of send step_coll to main node
      }
    }
    if ( ( iDebugLevel >= 4 ) ){
      printf("----- (generate_ephemerides) Done sending nb_coll_per_step_per_iproc_in_tca to main node. (iProc %d)\n", iProc);
    }

    if ( iProc == 0 ){
      if ( ( iDebugLevel >= 4 ) ){
	printf("----- (generate_ephemerides) Main node receiving nb_coll_per_step_per_iproc_in_tca... (iProc %d)\n", iProc);
      }

      nb_tca_without_collisions = 0;
      nb_tca_with_collisions = nb_tca - start_itca;
      for ( iiitca = start_itca; iiitca < nb_tca; iiitca++ ){
	for ( ppp = 1; ppp < nb_time_steps_in_tca_time_span-1-1; ppp++){ // "ppp = 1" and "ppp < nb_time_steps_in_tca_time_span-1" and not "ppp = 0" and "ppp < nb_time_steps_in_tca_time_span": see explanation below (comments starting with "ABOUT THE TIME SPAN").
	
	  for (ccc = 0; ccc < nProcs; ccc++){
	    if (ccc > 0){
	      if (nProcs > 1){
		MPI_Recv(&nb_coll_per_step_per_iproc_in_tca[ii][iiitca][ccc][ppp], 1, MPI_DOUBLE, ccc,2 * ii *   + iiitca * ( nb_time_steps_in_tca_time_span ) *  ( nProcs * OPTIONS->nb_ensemble_min_per_proc + 1 + 1 ) + ppp * ( 1 + nProcs * OPTIONS->nb_ensemble_min_per_proc ) + ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      }
	    }
	    nb_coll_per_step_in_TCA[ii][iiitca][ppp] = nb_coll_per_step_in_TCA[ii][iiitca][ppp] + nb_coll_per_step_per_iproc_in_tca[ii][iiitca][ccc][ppp];


	    /* 	  if ( ppp == 348 ){ */
	    /* if (nb_coll_per_step_per_iproc_in_tca[ii][iiitca][ccc][ppp] != 0){ */
	    /*   printf("nb_coll_per_step_in_TCA[%d][%d][%d]: %d | nb_coll_per_step_per_iproc_in_tca[%d][%d][%d][%d]: %d\n",  ii, iiitca, ppp, nb_coll_per_step_in_TCA[ii][iiitca][ppp], ii, iiitca, ccc, ppp, nb_coll_per_step_per_iproc_in_tca[ii][iiitca][ccc][ppp]); */
	    /* } */
	    /* } */

	    total_nb_collisions[ii][iiitca]  = total_nb_collisions[ii][iiitca]  + nb_coll_per_step_per_iproc_in_tca[ii][iiitca][ccc][ppp];
	  }
	}
	if (total_nb_collisions[ii][iiitca] == 0){
	  nb_tca_without_collisions = nb_tca_without_collisions + 1;
	}
      }
      nb_tca_with_collisions = nb_tca - start_itca - nb_tca_without_collisions;

      if ( ( iDebugLevel >= 4 ) ){
	printf("----- (generate_ephemerides) Main node done receiving nb_coll_per_step_per_iproc_in_tca. (iProc %d)\n", iProc);
      }

    } // end of if iproc == 0

  } // end of if collision assessment is on
   // END OF  COLLISION


  /* // PROPAGATE ONE LAST TIME UNTIL THE SATELLITE EPOCH IS THE SAME AS THE CONSTELLATION EPOCH */
  /*    for (ii = 0; ii < OPTIONS->n_satellites; ii++) { */
  /*      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.dt = min_end_time - CONSTELLATION->spacecraft[ii][0].et; */
  /*      //	et2utc_c(CONSTELLATION->spacecraft[ii][0].et, "ISOC" ,0 ,255 , times_att); */
  /* 	//		printf("before one last time: <%s>\n", times_att); */

  /*      propagate_spacecraft( &CONSTELLATION->spacecraft[ii][0], PARAMS, starttime, &density, GROUND_STATION ); */

  /*      // // ENSEMBLES SPACECRAFT (BUT NOT FOR GPS SATELLITES)  */
  /*      if (ii < OPTIONS->n_satellites - OPTIONS-> nb_gps){ */
  /* 	if ( OPTIONS->nb_ensembles_min > 0 ) { */
  /* 	  for (eee = 1L + iProc * OPTIONS->nb_ensemble_min_per_proc; eee< 1 + iProc * OPTIONS->nb_ensemble_min_per_proc + OPTIONS->nb_ensemble_min_per_proc ; eee++){ */
  /* 	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.dt = min_end_time - CONSTELLATION->spacecraft[ii][eee].et; */
  /* 	    propagate_spacecraft( &CONSTELLATION->spacecraft[ii][eee], PARAMS, starttime, &density, GROUND_STATION ); */
  /* 	  }  */
  /* 	}  */
  /*      } */

  /*      write_output( CONSTELLATION->spacecraft[ii] , 0, choose_tle_to_initialise_orbit, ii, OPTIONS->n_satellites,OPTIONS->nb_gps, OPTIONS->nb_ensembles_min,  OPTIONS->nb_ensemble_min_per_proc,  iProc,  OPTIONS->nb_ensembles_output_files, OPTIONS->filename_output_ensemble,previous_lat,OPTIONS);// name_file_ensembles_attitude, OPTIONS,PARAMS->EARTH.earth_fixed_frame); */
  /*    } */
  /* // END OF PROPAGATE ONE LAST TIME UNTIL THE SATELLITE EPOCH IS THE SAME AS THE CONSTELLATION EPOCH */


  for (ii = 0; ii < OPTIONS->n_satellites; ii++) {
    //if (iProc == 0){
    if ( start_ensemble[ii] == 0){ // if this iProc runs main sc ii
      fclose( CONSTELLATION->spacecraft[ii][0].fpecef ); // this was open by all iProc
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  int find_file_name;
  char *next;
  if (iProc == 0){
    if (OPTIONS->nb_gps > 0){
      char ch;
CONSTELLATION->file_CYGNSS_constellation_for_specular = fopen(CONSTELLATION->filename_CYGNSS_constellation_for_specular, "w+");
      for (ii = 0;ii< OPTIONS->n_satellites - OPTIONS->nb_gps;ii++){
      //      printf("<%s>\n", CONSTELLATION->spacecraft[ii][0].filenameecef);
      CONSTELLATION->spacecraft[ii][0].fpecef = fopen( CONSTELLATION->spacecraft[ii][0].filenameecef, "r");   // (ECEF format) (CBV)
    

/* 	//rewind(CONSTELLATION->spacecraft[ii][0].fpecef); */
    
	while( ( ch = fgetc(CONSTELLATION->spacecraft[ii][0].fpecef) ) != EOF )
	  fputc(ch,CONSTELLATION->file_CYGNSS_constellation_for_specular);
	fprintf(CONSTELLATION->file_CYGNSS_constellation_for_specular,"\n");
	fprintf(CONSTELLATION->file_CYGNSS_constellation_for_specular,"\n");
    fclose( CONSTELLATION->spacecraft[ii][0].fpecef ); // this was open by all iProc
      }

         CONSTELLATION->file_GPS_constellation_for_specular = fopen(CONSTELLATION->filename_GPS_constellation_for_specular, "w+");

      for (ii = OPTIONS->n_satellites - OPTIONS->nb_gps;ii< OPTIONS->n_satellites;ii++){
          CONSTELLATION->spacecraft[ii][0].fpecef = fopen( CONSTELLATION->spacecraft[ii][0].filenameecef, "r");   // (ECEF format) (CBV)
	//rewind(CONSTELLATION->spacecraft[ii][0].fpecef);
    
	while( ( ch = fgetc(CONSTELLATION->spacecraft[ii][0].fpecef) ) != EOF )
	  fputc(ch,CONSTELLATION->file_GPS_constellation_for_specular);
	fprintf(CONSTELLATION->file_GPS_constellation_for_specular,"\n");
	fprintf(CONSTELLATION->file_GPS_constellation_for_specular,"\n");
	fclose( CONSTELLATION->spacecraft[ii][0].fpecef ); // this was open by all iProc
      }
      fclose(CONSTELLATION->file_CYGNSS_constellation_for_specular);
      fclose(CONSTELLATION->file_GPS_constellation_for_specular);

    }
  } // end of if iProc == 0
  // COLLISION
  MPI_Barrier(MPI_COMM_WORLD);
  if ( iProc == 0 ){
    if ( (compute_collisions == 1) ){ // if collision assessment is on
      ii = 1; // !!!!!!!!! FOR NOW WORKS ONLY IF TWO REFERENCE SATELLIES ONLY

      file_collision = fopen(CONSTELLATION->filename_collision, "w+");

      fprintf(file_collision, "Total number of ensembles: %d.\nTime step of integration DT: %.2f seconds.\nNumber of close approaches: %d (+ %d without any collision).\n", nProcs * OPTIONS->nb_ensemble_min_per_proc, OPTIONS->dt, nb_tca_with_collisions, nb_tca_without_collisions);
      for ( iiitca = start_itca; iiitca < nb_tca; iiitca++ ){
	fprintf(file_collision, "Cumulative probability of collision during the time spanning TCA %d: %f.\n",iiitca + 1, (double)total_nb_collisions[ii][iiitca] / ( nProcs * OPTIONS->nb_ensemble_min_per_proc * nProcs * OPTIONS->nb_ensemble_min_per_proc ) );
      }
      for ( iiitca = start_itca; iiitca < nb_tca; iiitca++ ){
	//    pti(total_nb_collisions[ii][iiitca], "Total number of collisions");
	printf("CA %d - Total cumulative probability of collision: %.5f\n", iiitca + 1,(double)total_nb_collisions[ii][iiitca] / ( nProcs * OPTIONS->nb_ensemble_min_per_proc * nProcs * OPTIONS->nb_ensemble_min_per_proc )) ;
	  fprintf(file_collision, "\n#Detailed results for close approach %d\n", iiitca+1);
	  et2utc_c(save_tca[iiitca], "ISOC", 3, 255, time_itca);
	  fprintf(file_collision, "TCA: %s.\n", time_itca);
	  fprintf(file_collision, "DCA: %.1f m.\n", save_dca[iiitca] * 1000.); // use dca with caution (see comment "this is incorrect. What would be correct would be to take the distance at inital epoch + period / 2 but we did not save this distance. Since save_dca is never used after, we don't really care, except in the collision output file where dca is reported")
	  fprintf(file_collision, "Primary object sample near conjunction point: %d\n", eee_prim_that_collide);// this is one of the samples that is near the conjnuction point (see function compute_collision_between_one_secondary_and_all_primary for the definition of "near")
	if ( total_nb_collisions[ii][iiitca] > 0 ){
	  fprintf(file_collision, "Cumulative Probability of collision during the time spanning TCA: %f.\n", (double)total_nb_collisions[ii][iiitca] / ( nProcs * OPTIONS->nb_ensemble_min_per_proc * nProcs * OPTIONS->nb_ensemble_min_per_proc ) );
	  fprintf(file_collision, "Total number of collisions: %d.\n", total_nb_collisions[ii][iiitca] );
	  et2utc_c(et_time_step_of_save_tca[iiitca] - ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt) + OPTIONS->dt , "ISOC", 3, 255, time_itca); // "+ OPTIONS->dt": see explanation below (comments starting with "ABOUT THE TIME SPAN")
	  fprintf(file_collision, "Time spanning TCA %d: %s - ", iiitca + 1, time_itca);
	  et2utc_c(et_time_step_of_save_tca[iiitca] + ((int)(nb_time_steps_in_tca_time_span / 2)* OPTIONS->dt )- OPTIONS->dt, "ISOC", 3, 255, time_itca); // "- OPTIONS->dt": see explanation below (comments starting with "ABOUT THE TIME SPAN")
	  fprintf(file_collision, "%s (duration of span: %.3f seconds (%d time steps)).\n", time_itca, CONSTELLATION->collision_time_span - 2 * OPTIONS->dt , nb_time_steps_in_tca_time_span - 2 ); // "- 2 * OPTIONS->dt" and "- 2": see explanation below (comments starting with "ABOUT THE TIME SPAN").
	  fprintf(file_collision, "##REAL_TIME NB_COLLISIONS_PER_DT\n" );
	  for ( ppp = 1; ppp < nb_time_steps_in_tca_time_span-1-1; ppp++){// "ppp = 1" and "ppp < nb_time_steps_in_tca_time_span-1-1" and not "ppp = 0" and "ppp < nb_time_steps_in_tca_time_span-1": see explanation below (comments starting with "ABOUT THE TIME SPAN"). The reason why we remove -1 (in addition to removing another -1 (reason explained in "ABOUT THE TIME SPAN")) is that the last time step number does not have any collision in the time step following it because it is the time step number that is at the border of the interval. So better explained: if tn is the last time step number of the time span, by definition there is no collision in [tn, tn+1] because we stopped the algorithm at the interval [tn-1, tn]. The number of collisions in [tn-1, tn] corresponds to the time step number tn-1. This is why we don't look at the time step number tn.
	    et_step_collision = et_time_step_of_save_tca[iiitca] - ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt ) + ppp * OPTIONS->dt;
	    et2utc_c(et_step_collision, "ISOC", 3, 255, time_itca);
	    et2utc_c(et_step_collision+ OPTIONS->dt, "ISOC", 3, 255, time_itca2);
	    fprintf(file_collision, "%s -> %s %d\n", time_itca, time_itca2, nb_coll_per_step_in_TCA[ii][iiitca][ppp]);
	  }
	  fprintf(file_collision, "#End of results for close approach %d\n", iiitca+1);
	}
      }

      fclose(file_collision);
    } // end of if collision assessment is on
    else if (compute_collisions_was_on_but_no_tca_found == 1){ // the user chose to compute collision but to TCA was found between the initial epoch and the final epoch
      file_collision = fopen(CONSTELLATION->filename_collision, "w+");
      fprintf(file_collision, "No close approach (< %.2f meters) for the unperturbed orbits was found between the initial epoch (%s) and the final epoch (%s).\n", OPTIONS->min_dist_close_approach*1000, OPTIONS->initial_epoch, OPTIONS->final_epoch);
      fclose(file_collision);
    }

  } // end of if iProc == 0
/* // ABOUT THE TIME SPAN */
/*   The algorithm is unlikely to find a min distance in the last time step of the span (and in the first time step of the span). */
/*     Reason: imagine there is a min distance in the last time step dt. We note tn-1 the start of this time step (dn-1 the distance at tn-1), and tn the end of the time step (dn the distance at tn). The algorithm finds the min distance only if dn > dn-1. If dn has still not reached a value > dn-1 after the min in dt, then the min will me missed by the algorithm. To find it, the algorithm would need to look at dn+1. dn+1 is for sure > dn (because the min was between tn-1 and tn) so the algorithm will know that there is min somewhere between tn-1 and tn+1, and will therefore call ANCAS to find it. */
/*   -> that’s ok, there is nothing we can do about it (there actually is but then the algorithm gets complicated and it’s not worth it). A good solution though is to tell SpOCK that the span is until tn, knowing that the collisions between tn-1 and tn will be missed, and write the report as if the span ends at tn-1. If for some reason the user really wants to look at collisions between tn-1 and tn, then the user just has to tell SpOCK that the span ends at tn+1. Therefore, SpOCK will miss the collisions between tn and tn+1 (the “new” last step), but not the ones between tn-1 and tn.							 */
/*   The same reasoning applies to the first time step. As collisions between t0 and t1 will be missed, then write the report starting at t1. If the user wants to compute the collisions between t0 and t1, the user just needs to tell SpOCK to start the span at t-1. */
/*   This implies that the span will no matter what end no later than one time step before the propagation ends, and will start one time step after the propagation starts (actually 3 time steps but this is for another reason…). So if for instance the user wants to find collisions until 12:10:00 and the time step is 60 s, then the end epoch (2nd line of section #TIME in the main input file) must be at least 12:11:00 (it can also be 12:12:00, 12:13:00, …!). Same for the start epoch: is the user wants to compute collisions starting at 09:33:00 and the time step is 60 s then the start epoch must be 09:30:00 (as we said, it’s 3 time steps for the initial epoch, one of the 3 is for the reason we have mentioned in this slide, the 2 others are for another reason). */
/* -> message to get: if the time spanning the TCA is delta_t, collisions are actually computed only for delta_t minus two time steps (the first and the last time steps are not taken into account in the collision assessment). */
/* // end of ABOUT THE TIME SPAN */

  // end of COLLISION
  
  // Close files needed

  for (ii = 0; ii < OPTIONS->n_satellites; ii++) {
    if (iProc == 0){
      //    fclose( CONSTELLATION->spacecraft[ii][0].fpecef ); // this was open by all iProc
    }
    if ( start_ensemble[ii] == 0){ // if this iProc runs main sc ii
      if ( CONSTELLATION->spacecraft[ii][0].INTEGRATOR.isGPS == 0 ){
	fclose( CONSTELLATION->spacecraft[ii][0].fp );
	fclose( CONSTELLATION->spacecraft[ii][0].fpout );
	fclose( CONSTELLATION->spacecraft[ii][0].fprho );
	fclose( CONSTELLATION->spacecraft[ii][0].fpatt );
	}
  //  fclose( CONSTELLATION->spacecraft[ii][0].fpecef ); // this was open by all iProc

      // comment line below if you did not open this file at the beginning of this function
      if (ii < OPTIONS->n_satellites - OPTIONS->nb_gps){
	fclose(CONSTELLATION->spacecraft[ii][0].INTEGRATOR.file_given_output);
      }
      if (GROUND_STATION->nb_ground_stations > 0){
	if (ii < OPTIONS->n_satellites - OPTIONS->nb_gps){ //do not compute the ground station coverage for gps
	  for (iground = 0; iground < GROUND_STATION->nb_ground_stations; iground++){
	    fclose(CONSTELLATION->spacecraft[ii][0].fp_coverage_ground_station[iground]);
	  }
	}
      }
      if (ii < OPTIONS->n_satellites - OPTIONS-> nb_gps){
	if (CONSTELLATION->spacecraft[ii][0].INTEGRATOR.solar_cell_efficiency != -1){
	  fclose( CONSTELLATION->spacecraft[ii][0].fpower );
	  fclose( CONSTELLATION->spacecraft[ii][0].fpeclipse );
	}
      }
    } // end of if this iProc runs main sc ii
    //} // end of if iproc == 0

    if (ii < OPTIONS->n_satellites - OPTIONS-> nb_gps){
      if ( OPTIONS->nb_ensembles_min > 0 ){
	if ( array_sc[1] > 0 )  { // if this iProc runs ensembles (otherwise array_sc[1] = - 1)
	for (fff = 0; fff < OPTIONS->nb_ensembles_output_files ; fff ++){
	  if ( (strcmp(OPTIONS->filename_output_ensemble[fff], "tca" )!= 0 ) && (strcmp(OPTIONS->filename_output_ensemble[fff], "dca" )!= 0 ) && (strcmp(OPTIONS->filename_output_ensemble[fff], "sample" )!= 0 )){
	    fclose(CONSTELLATION->spacecraft[ii][0].fpiproc[fff]);
	  }
	  else if (ii == 0){
	    if (strcmp(OPTIONS->filename_output_ensemble[fff], "tca" ) == 0){
	      fclose( tca_file );
	    }
	    if (strcmp(OPTIONS->filename_output_ensemble[fff], "dca" ) == 0){
	      fclose( dca_file );
	    }
/* 	    if (strcmp(OPTIONS->filename_output_ensemble[fff], "sample" ) == 0){ */
/* 	      fclose( sample_file ); */
/* 	    } */

	  }
	}
	} // end of if this iProc runs ensembles (otherwise array_sc[1] = - 1)
      }
      /* if ( OPTIONS->nb_ensembles_attitude > 0 ){ */
      /* 	if (ii == 0){ // for now, ensemble on attitude if there is only one satellite in the constellation */
      /* 	for (fff = 0; fff < n_files_ensembles_attitude ; fff ++){ */
      /* 	  fclose(CONSTELLATION->spacecraft[ii][0].fpiproc_attitude[iProc][fff]); */
      /* 	} */
      /* 	} */
      /* } */
    }
    free( CONSTELLATION->spacecraft[ii]);
  }




  // computes if the mean spacecraft flies over the storm
  /* if (OPTIONS->nb_storm > 0){ */
  /* for (ii = 0; ii < OPTIONS->n_satellites - OPTIONS-> nb_gps; ii++) { */
  /*   fly_storm( &CONSTELLATION->spacecraft[ii][0], OPTIONS, PARAMS, STORM, ii); */
  /* } */
  /* } */

  /* fclose(fp_temp_att1); */
  /* fclose(fp_temp_att_us); */
  /* fclose(fp_temp_att2); */

  //  fclose(density_file);




  if (( iDebugLevel >= 1 ) ){
    printf("-- (generate_ephemerides) Just got out from generate_ephemerides.(iProc %d)\n", iProc);
  }


   
 return 0;

}


/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           write_output
//  Purpose:        Writes output for 2 files per spacecraft.
//                      File scXout.txt contains various keplerian, inertial, and ecef parameters
//                      File scXout_LLA.txt is Dr. Aaron Ridley's requested output for gitm
//                      See header files for more info
//  Assumptions:    None
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | J. Getchius   | 05/20/2015    |   ---     | Initial Implementation
//
/////////////////////////////////////////////////////////////////////////////////////////

int write_output(   SPACECRAFT_T    *SC,
                    int             init_flag,
                    int             choose_tle_to_initialise_orbit,
		    int             sc_index,
                    int             n_sc,
		    int             n_gps,
		    int             nb_ensembles_min,
		    //		    int             nb_ensembles_attitude,
		    int             nb_ensemble_min_per_proc,
		    //		    int             nb_ensemble_attitude_per_proc,
		    int             iProc,
		    int             nb_ensembles_output_files,
		    char            name_file_ensembles[30][1000],// !!! the number of columns has to correpond to n_files_ensembles
		    // char            name_file_ensembles_attitude[256][1], // !!! the number of columns has to correpond to n_files_ensembles_attitude
		    double          previous_lat,
		    OPTIONS_T       *OPTIONS,
		    char earth_fixed_frame[100], 
		    int write_reference_sc,
		    int write_ensembles,
		    double constellation_et,
		    int nProcs,
		    int iDebugLevel,
		    int compute_collisions,
		    int *start_ensemble,
		    int *array_sc,
		    CONSTELLATION_T *CONSTELLATION,
		    PARAMS_T *PARAMS
		    )
{

  // if set to 0 then don't write in the given files. if set to 1 the write in given files

  int write_lla = 1;
  int write_rho = 1;
    int write_attitude = 1;
  int write_tle = 1;
  int write_ecef = 1;
  int write_all = 1;
  if ( SC[0].INTEGRATOR.isGPS == 1 ){
    write_lla = 0;
    write_rho = 0;
    write_attitude = 0;
    write_tle = 0;
    write_all = 0;
  }


  /* char str_nProcs[100]; */
  /* strcpy(str_nProcs, ""); */
  /* sprintf(str_nProcs, "%d", nProcs); */

  if ( (iProc == 0) && ( iDebugLevel >= 2 ) ){
    printf("-- (write_output) Just got in write_output.(iProc %d)\n", iProc);
  }
  double et_intial_epoch;
  str2et_c(OPTIONS->initial_epoch, &et_intial_epoch);

  // Declarations
  SpiceDouble       xform[6][6];
  double estate[6], jstate[6];

  double total_power;
  int hr, mn, sc;
  char time_local[256], ampm[256];
  double lon_an_dn = 0 ;

  double x[3];
  double lt;
  int eee, fff, sss;
  char    times[256], times_with_milliseconds[256];
  /*   char    mm[256]; */
  /*   char    yy[256]; */
  int     ii;
  /*   int     rv; */
  char text[256];
  char text_time_with_milliseconds[256];
    char text_without_millisecond[256];

  int nProcs_that_are_gonna_run_ensembles;
  nProcs_that_are_gonna_run_ensembles = nProcs;
  if ( nProcs > OPTIONS->nb_ensembles_min ){
    nProcs_that_are_gonna_run_ensembles = OPTIONS->nb_ensembles_min;
  }



  // TLE epoch of the GPS or the satellite (if the orbit initialisation of the satellite was done from a TLE)
  //  if (iProc == 0){
  if (start_ensemble[sc_index] == 0){ // if this iProc runs main sc sc_index
  et2utc_c(SC[0].et, "ISOC" , 6, 255 , times);
  sscanf(times, "%4[^\n]",text);
  strncat(text, "/",1);
  strncat(text, &times[5],1);
  strncat(text, &times[6],1);
  strncat(text, "/",1);
  strncat(text, &times[8],1);
  strncat(text, &times[9],1);
  strncat(text, " ",1);
  for (ii = 0; ii<15; ii++){
    strncat(text, &times[11+ii],1);
  }

    et2utc_c(SC[0].et, "ISOC" , 0, 255 , times);
  sscanf(times, "%4[^\n]",text_without_millisecond);
  strncat(text_without_millisecond, "/",1);
  strncat(text_without_millisecond, &times[5],1);
  strncat(text_without_millisecond, &times[6],1);
  strncat(text_without_millisecond, "/",1);
  strncat(text_without_millisecond, &times[8],1);
  strncat(text_without_millisecond, &times[9],1);
  strncat(text_without_millisecond, " ",1);
  for (ii = 0; ii<9; ii++){
    strncat(text_without_millisecond, &times[11+ii],1);
  }



  et2utc_c(SC[0].et, "ISOC" ,6 ,255 , times_with_milliseconds);
  sscanf(times_with_milliseconds, "%4[^\n]",text_time_with_milliseconds);
  strncat(text_time_with_milliseconds, "/",1);
  strncat(text_time_with_milliseconds, &times_with_milliseconds[5],1);
  strncat(text_time_with_milliseconds, &times_with_milliseconds[6],1);
  strncat(text_time_with_milliseconds, "/",1);
  strncat(text_time_with_milliseconds, &times_with_milliseconds[8],1);
  strncat(text_time_with_milliseconds, &times_with_milliseconds[9],1);
  strncat(text_time_with_milliseconds, " ",1);
  for (ii = 0; ii<15; ii++){
    strncat(text_time_with_milliseconds, &times_with_milliseconds[11+ii],1);
  }
  } // end of  if this iProc runs main sc sc_index
  //  } // end of if iproc == 0

  if (init_flag == 1) {


    //if (iProc == 0){
    if (start_ensemble[sc_index] == 0){ // if this iProc runs main sc sc_index

    // Initialization header
    if (write_reference_sc){

      if (write_lla == 1){

      fprintf(SC[0].fpout, "// ---------------------------------------------------------------------------------- \n");

      fprintf(SC[0].fpout, "// \n");
      fprintf(SC[0].fpout, "// Please note this file is auto generated by the Spacecraft Orbital Characterization Kit (SpOCK)\n");
      fprintf(SC[0].fpout, "// Trajectory Spacecraft %s in a constellation of %d spacecraft and %d GPS satellites \n", SC[0].name_sat, n_sc-n_gps, n_gps);
      if ( ( sc_index >= (n_sc - n_gps) ) || ( choose_tle_to_initialise_orbit == 1 ) )
	fprintf(SC[0].fpout ,"// Epoch of last TLE for %s: %s \n",  SC[0].name_sat, text);
      fprintf(SC[0].fpout, "// \n");
      fprintf(SC[0].fpout, "// Trajectory is specified with:  TIME LONGITUDE(DEG) LATITUDE(DEG) ALTITUDE(KM) \n");
      fprintf(SC[0].fpout, "// \n");
      fprintf(SC[0].fpout, "// Version control of SpOCK is under Joel Getchius's Mac and Charles Bussy-Virat's Mac \n");
      fprintf(SC[0].fpout, "// \n");
      fprintf(SC[0].fpout, "// ---------------------------------------------------------------------------------- \n");
      fprintf(SC[0].fpout, "#START \n");
      }

      if (write_rho == 1){
      fprintf(SC[0].fprho, "// ---------------------------------------------------------------------------------- \n");
      fprintf(SC[0].fprho, "// \n");
      fprintf(SC[0].fprho, "// Please note this file is auto generated by the Spacecraft Orbital Characterization Kit (SpOCK) \n");
      fprintf(SC[0].fprho, "// Trajectory Spacecraft %s in a constellation of %d spacecraft and %d GPS satellites \n", SC[0].name_sat, n_sc-n_gps, n_gps);
      if (sc_index >= (n_sc - n_gps))
	fprintf(SC[0].fprho ,"// Epoch of last TLE for %s: %s \n",  SC[0].name_sat, text);
      fprintf(SC[0].fprho, "// \n");
      fprintf(SC[0].fprho, "// Version control of SpOCK is under Joel Getchius's Mac and Charles Bussy-Virat's Mac \n");
      fprintf(SC[0].fprho, "// \n");
      fprintf(SC[0].fprho, "// Atmosphere info: \n");
      fprintf(SC[0].fprho, "// TIME DENSITY(kg/km^3) TEMP(K) CD_TOT_NORM A_REF_TOT(km^2) BC*mass(m^2)");

      /* for (fff = 0; fff < OPTIONS->nb_storm; fff++){ */
      /*   fprintf(SC[0].fprho, "SEE %s ", STORM->storm_name[fff]); */
      /* } */
   
      fprintf(SC[0].fprho,"\n");

      fprintf(SC[0].fprho, "// ---------------------------------------------------------------------------------- \n");
      }



      

            if (write_attitude == 1){
      fprintf(SC[0].fpatt, "// ---------------------------------------------------------------------------------- \n");



      fprintf(SC[0].fpatt, "// \n");
      fprintf(SC[0].fpatt, "// Please note this file is auto generated by the Spacecraft Orbital Characterization Kit (SpOCK)\n");
      fprintf(SC[0].fpatt, "// Attitude Spacecraft %s in a constellation of %d spacecraft and %d GPS satellites \n", SC[0].name_sat, n_sc-n_gps, n_gps);
      if ( ( sc_index >= (n_sc - n_gps) ) || ( choose_tle_to_initialise_orbit == 1 ) )
	fprintf(SC[0].fpatt ,"// Epoch of last TLE for %s: %s \n",  SC[0].name_sat, text);
      fprintf(SC[0].fpatt, "// \n");
      fprintf(SC[0].fpatt, "// Attitude is specified with:  TIME PITCH(DEG) ROLL(DEG) YAW(DEG) ORDER_PITCH ORDER_ROLL ORDER_YAW \n");
      fprintf(SC[0].fpatt, "// \n");
      fprintf(SC[0].fpatt, "// Version control of SpOCK is under Joel Getchius's Mac and Charles Bussy-Virat's Mac \n");
      fprintf(SC[0].fpatt, "// \n");
      fprintf(SC[0].fpatt, "// ---------------------------------------------------------------------------------- \n");
      fprintf(SC[0].fpatt, "#START \n");
      }

      // Initialization header for TLE output file
      if (sc_index < n_sc - n_gps){
      if ( choose_tle_to_initialise_orbit == 1 ) {
	if ( write_tle == 1 ){
	fprintf(SC[0].fptle, "// ---------------------------------------------------------------------------------- \n");
	fprintf(SC[0].fptle, "// \n");
	fprintf(SC[0].fptle, "// Please note this file is auto generated by the Spacecraft Orbital Characterization Kit (SpOCK) \n");
	fprintf(SC[0].fptle, "// Trajectory Spacecraft %s in a constellation of %d spacecraft and %d GPS satellites \n", SC[0].name_sat, n_sc-n_gps, n_gps);
	fprintf(SC[0].fptle ,"// This file shows the position of satellite %s from its TLE epoch start (%s) to the constellation epoch start (%s)\n",  SC[0].name_sat, text_time_with_milliseconds, OPTIONS->initial_epoch);
	fprintf(SC[0].fptle, "// \n");
	fprintf(SC[0].fptle, "// Trajectory is specified with: TIME  \t\tr_i2cg_INRTL(3)(KM) \t\t     v_i2cg_INRTL(3)(KM/S) \t      LONG(DEG)    LAT(DEG)     ALT(KM)       SMA(KM)     INC(DEG)      ECC       TRUE ANO(DEG)  RAAN(DEG) ARG PERIG(DEG) RIGHT ASC(DEG) LOCAL TIME(DEG)  a_i2cg_INRTL(3)(KM/S^2) r_ecef2cg_ECEF(3)(KM) v_ecef2cg_ECEF(3)(KM)\n");
	fprintf(SC[0].fptle, "// \n");
	fprintf(SC[0].fptle, "// Version control of SpOCK is under Joel Getchius's Mac and Charles Bussy-Virat's Mac \n");
	fprintf(SC[0].fptle, "// \n");
	fprintf(SC[0].fptle, "// ---------------------------------------------------------------------------------- \n");
	fprintf(SC[0].fptle, "#START \n");
	}
      }
      }

      // Initialization header
      if (write_ecef == 1){
      fprintf(SC[0].fpecef, "// ---------------------------------------------------------------------------------- \n");
      fprintf(SC[0].fpecef, "// \n");
      fprintf(SC[0].fpecef, "// Please note this file is auto generated by the Spacecraft Orbital Characterization Kit (SpOCK) \n");
      fprintf(SC[0].fpecef, "// Trajectory Spacecraft %s in a constellation of %d spacecraft and %d GPS satellites (iProc %d)\n", SC[0].name_sat, n_sc-n_gps, n_gps, iProc);
      if (sc_index >= (n_sc - n_gps))
	fprintf(SC[0].fpecef ,"// Epoch of last TLE for %s: %s \n",  SC[0].name_sat, text);
      fprintf(SC[0].fpecef, "// \n");
      fprintf(SC[0].fpecef, "// Trajectory is specified with: \n");
      fprintf(SC[0].fpecef, "//\tTIME  \t\t      r_ecef2cg_ECEF(3)(KM) \t\t  v_ecef2cg_ECEF(3)(KM/S)\n");
      fprintf(SC[0].fpecef, "// \n");
      fprintf(SC[0].fpecef, "// Version control of SpOCK is under Joel Getchius's Mac and Charles Bussy-Virat's Mac \n");
      fprintf(SC[0].fpecef, "// ---------------------------------------------------------------------------------- \n");
      fprintf(SC[0].fpecef, "#START \n");
    
      }

      // Initialization header
      if (write_all == 1){
      fprintf(SC[0].fp, "// ---------------------------------------------------------------------------------- \n");
      fprintf(SC[0].fp, "// \n");
      fprintf(SC[0].fp, "// Please note this file is auto generated by the Spacecraft Orbital Characterization Kit (SpOCK) \n");
      fprintf(SC[0].fp, "// Trajectory Spacecraft %s in a constellation of %d spacecraft and %d GPS satellites \n", SC[0].name_sat, n_sc-n_gps, n_gps);
      if (sc_index >= (n_sc - n_gps))
	fprintf(SC[0].fp ,"// Epoch of last TLE for %s: %s \n",  SC[0].name_sat, text);
      fprintf(SC[0].fp, "// \n");
      fprintf(SC[0].fp, "// Version control of SpOCK is under Joel Getchius's Mac and Charles Bussy-Virat's Mac \n");
      fprintf(SC[0].fp, "// \n");
      fprintf(SC[0].fp, "// Trajectory is specified with: \n");
      fprintf(SC[0].fp, "//    ET  \t\tr_i2cg_INRTL(3)(KM) \t\t     v_i2cg_INRTL(3)(KM/S) \t      LONG(DEG)    LAT(DEG)     ALT(KM)       SMA(KM)     INC(DEG)      ECC       TRUE ANO(DEG)  RAAN(DEG) ARG PERIG(DEG) RIGHT ASC(DEG) LOCAL TIME(DEG)  a_i2cg_INRTL(3)(KM/S^2) a_i2cg_LVLH(3)(KM/S^2) a_i2cg_LVLH_gravity(3)(KM/S^2) a_i2cg_LVLH_drag(3)(KM/S^2) a_i2cg_INRTL_drag(3)(KM/S^2) a_i2cg_INRTL_gravity(3)(KM/S^2) beta_angle(deg) PHASE_ANGLE(DEG) ARG_PERIG_AVER(DEG) SMA_AVER(KM) ECC_AVER SOLAR_ZENITH(DEG)");

      /* for (fff = 0; fff < OPTIONS->nb_storm; fff++){ */
      /*   fprintf(SC[0].fp, "SEE %s ", STORM->storm_name[fff]); */
      /* } */
   
      fprintf(SC[0].fp,"\n");

      fprintf(SC[0].fp, "// ---------------------------------------------------------------------------------- \n");
      }
    }


    if (sc_index < n_sc - n_gps){
      if (write_reference_sc == 1){
	// Initialization header for power file
	if (SC[0].INTEGRATOR.solar_cell_efficiency != -1){
	  fprintf(SC[0].fpower, "// ---------------------------------------------------------------------------------- \n");
	  fprintf(SC[0].fpower, "// \n");
	  fprintf(SC[0].fpower, "// Please note this file is auto generated by the Spacecraft Orbital Characterization Kit (SpOCK) \n");
	  fprintf(SC[0].fpower, "// Spacecraft power %s in a constellation of %i spacecraft and %d GPS satellites\n", SC[0].name_sat, n_sc - n_gps, n_gps);
	  fprintf(SC[0].fpower, "// Number of surfaces: %f \n", SC[0].INTEGRATOR.nb_surfaces);
	  fprintf(SC[0].fpower, "// \n");
	  fprintf(SC[0].fpower, "// Trajectory is specified with:  TIME POWER(W) LIGHT/SHADOW SUN_ELEVATION\n");
	  fprintf(SC[0].fpower, "// Sun elevation angle counted from body xy plane (positive if Sun in the direction of +z body) (-999 if satellite is not in the light)\n");
	  fprintf(SC[0].fpower, "// \n");
	  fprintf(SC[0].fpower, "// Version control of SpOCK is under Joel Getchius's Mac and Charles Bussy-Virat's Mac \n");
	  fprintf(SC[0].fpower, "// \n");
	  fprintf(SC[0].fpower, "// ---------------------------------------------------------------------------------- \n");
	  fprintf(SC[0].fpower, "#START \n");
	}
      }
    }
    } // end of  if this iProc runs main sc sc_index
    //    } // end if iproc == 0

    if (sc_index < n_sc - n_gps){
      // Initialization header
      // ensembles in COE
      //      if ( write_ensembles == 1 ){
	if ( nb_ensemble_min_per_proc > 0 ) {
	if ( array_sc[1] > 0 )  { // if this iProc runs ensembles (otherwise array_sc[1] = - 1)
	  for (fff = 0; fff < nb_ensembles_output_files ; fff ++){
	    if ( (strcmp(OPTIONS->filename_output_ensemble[fff], "tca" )!= 0 ) && (strcmp(OPTIONS->filename_output_ensemble[fff], "dca" )!= 0 ) && (strcmp(OPTIONS->filename_output_ensemble[fff], "sample" )!= 0 )){
/* 	    if (iProc == 1){ */
/* 	      printf("iProc %d opening file <%s> for sc %d\n", iProc, SC[0].filenameiproc[fff], sc_index);  */
/*  	    } */

	    fprintf(SC[0].fpiproc[fff], "// ---------------------------------------------------------------------------------- \n");
	    fprintf(SC[0].fpiproc[fff], "// \n");
	    fprintf(SC[0].fpiproc[fff], "// Please note this file is auto generated by the Spacecraft Orbital Characterization Kit (SpOCK) \n");
	    ///   fprintf(SC[0].fpiproc[fff], "// Ensembles for Spacecraft %s in a constellation of %d spacecraft and %d GPS satellites \n", SC[0].name_sat, n_sc-n_gps, n_gps);
   fprintf(SC[0].fpiproc[fff], "// Ensembles in a constellation of %d spacecraft and %d GPS satellites \n",  n_sc-n_gps, n_gps);

	    fprintf(SC[0].fpiproc[fff], "// \n");
	    fprintf(SC[0].fpiproc[fff], "// Trajectory is specified with:  time %s (%d ensembles - %d processors) \n", name_file_ensembles[fff], nb_ensembles_min, nProcs_that_are_gonna_run_ensembles);
	    fprintf(SC[0].fpiproc[fff], "// \n");
	    fprintf(SC[0].fpiproc[fff], "// Version control of SpOCK is under Joel Getchius's Mac and Charles Bussy-Virat's Mac \n");
	    fprintf(SC[0].fpiproc[fff], "// \n");
	    fprintf(SC[0].fpiproc[fff], "// ---------------------------------------------------------------------------------- \n");
	    fprintf(SC[0].fpiproc[fff], "#START");
	  }
	  }
	} // end of if this iProc runs ensembles (otherwise array_sc[1] = - 1)

	}
	//      }
      // ensemble in attitude
      /* if ( nb_ensemble_attitude_per_proc > 0 ) {
      /* 	for (fff = 0; fff < ( sizeof(name_file_ensembles_attitude[256])/sizeof(name_file_ensembles_attitude[256][0]) ) ; fff ++){  */
      /* 	    fprintf(SC[0].fpiproc_attitude[iProc][fff], "// ---------------------------------------------------------------------------------- \n"); */
      /* 	    fprintf(SC[0].fpiproc_attitude[iProc][fff], "// \n"); */
      /* 	    fprintf(SC[0].fpiproc_attitude[iProc][fff], "// Please note this file is auto generated by the Spacecraft Orbital Characterization Kit (SpOCK) \n"); */
      /* 	    fprintf(SC[0].fpiproc_attitude[iProc][fff], "// Ensembles for Spacecraft %s in a constellation of %d spacecraft and %d GPS satellites \n", SC[0].name_sat, n_sc-n_gps, n_gps);     */
      /* 	    fprintf(SC[0].fpiproc_attitude[iProc][fff], "// \n"); */
      /* 	    fprintf(SC[0].fpiproc_attitude[iProc][fff], "// Power is specified with:  time %s (%d ensembles)\n", name_file_ensembles_attitude[fff], nb_ensembles_attitude); */
      /* 	    fprintf(SC[0].fpiproc_attitude[iProc][fff], "// \n"); */
      /* 	    fprintf(SC[0].fpiproc_attitude[iProc][fff], "// Version control of SpOCK is under Joel Getchius's Mac and Charles Bussy-Virat's Mac \n"); */
      /* 	    fprintf(SC[0].fpiproc_attitude[iProc][fff], "// \n"); */
      /* 	    fprintf(SC[0].fpiproc_attitude[iProc][fff], "// ---------------------------------------------------------------------------------- \n"); */
      /* 	    fprintf(SC[0].fpiproc_attitude[iProc][fff], "#START"); */
      /* 	} */
      /*       } */


    }

    if ( ( sc_index >= n_sc - n_gps) || ( choose_tle_to_initialise_orbit == 1 ) ){

      return 0;
    }
    
  }


  /* Write the results in the TLE format file */
  //  if ( iProc == 0){
  if (start_ensemble[sc_index] == 0){ // if this iProc runs main sc sc_index
  if (write_reference_sc == 1){
    if ( init_flag == 2  ){ // init_flag = 2 means that we write the positions of the satellite from its TLE epoch start to the constellation epoch start
      if (write_tle==1){
    fprintf(SC[0].fptle, "%s", text_time_with_milliseconds );
    fprintf(SC[0].fptle, " " );

    fprintf(SC[0].fptle, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
    	    SC[0].r_i2cg_INRTL[0],
    	    SC[0].r_i2cg_INRTL[1],
    	    SC[0].r_i2cg_INRTL[2],
    	    SC[0].v_i2cg_INRTL[0],
    	    SC[0].v_i2cg_INRTL[1],
    	    SC[0].v_i2cg_INRTL[2],
    	    (SC[0].GEODETIC.longitude * RAD2DEG),
    	    (SC[0].GEODETIC.latitude  * RAD2DEG),
    	    SC[0].GEODETIC.altitude,
    	    SC[0].OE.sma,
    	    (SC[0].OE.inclination * RAD2DEG),
    	    (SC[0].OE.eccentricity),
    	    (SC[0].OE.f * RAD2DEG),
    	    (SC[0].OE.long_an*RAD2DEG),
    	    (SC[0].OE.w*RAD2DEG),
    	    SC[0].OE.ra*RAD2DEG,
    	    SC[0].r_ecef2cg_ECEF[0],
    	    SC[0].r_ecef2cg_ECEF[1],
    	    SC[0].r_ecef2cg_ECEF[2],
    	    SC[0].v_ecef2cg_ECEF[0],
    	    SC[0].v_ecef2cg_ECEF[1],
    	    SC[0].v_ecef2cg_ECEF[2]);
      }

      /* fprintf(SC[0].fptle, "%s", text ); */
      /* fprintf(SC[0].fptle, " " ); */
      /* fprintf(SC[0].fptle, "%f %f %f %6.3f %6.3f %6.3f \n", */
      /* 	      SC[0].r_i2cg_INRTL[0], */
      /* 	      SC[0].r_i2cg_INRTL[1], */
      /* 	      SC[0].r_i2cg_INRTL[2], */
      /* 	      SC[0].GEODETIC.longitude * RAD2DEG, */
      /* 	      SC[0].GEODETIC.latitude  * RAD2DEG, */
      /* 	      SC[0].GEODETIC.altitude); */
      return 0;
    } // end of init_flag = 2 means that we write the positions of the satellite from its TLE epoch start to the constellation epoch start
  } // end of if (write_reference_sc == 1)

  /* Write the results in the detailed format file */
  if (write_reference_sc == 1){

    if (write_all == 1){
    fprintf(SC[0].fp, "%s", text );
    fprintf(SC[0].fp, " " );
    fprintf(SC[0].fp, "%.10f %.10f %.10f %.10f %.10f %.10f %f %f %f %f %f %f %f %f %f %f",
	    SC[0].r_i2cg_INRTL[0],
	    SC[0].r_i2cg_INRTL[1],
	    SC[0].r_i2cg_INRTL[2],
	    SC[0].v_i2cg_INRTL[0],
	    SC[0].v_i2cg_INRTL[1],
	    SC[0].v_i2cg_INRTL[2],
	    (SC[0].GEODETIC.longitude * RAD2DEG),
	    (SC[0].GEODETIC.latitude  * RAD2DEG),
	    SC[0].GEODETIC.altitude,
	    SC[0].OE.sma,
	    (SC[0].OE.inclination * RAD2DEG),
	    (SC[0].OE.eccentricity),
	    (SC[0].OE.f * RAD2DEG),
	    (SC[0].OE.long_an*RAD2DEG),
	    (SC[0].OE.w*RAD2DEG),
	    SC[0].OE.ra*RAD2DEG);
    /* for (fff = 0; fff < OPTIONS->nb_storm; fff++){ */
    /*   fprintf(SC[0].fp," %d", SC[0].see_storm[fff]); */
    /* } */


    // !!! OUTPUT THE LOCAL TIME OF THE SATELLITE (IN DEGREES). UNCOMMENT IF YOU WANT TO OUTPUT IT!
    et2lst_c ( SC[0].et,  399,  SC[0].GEODETIC.longitude, "PLANETOCENTRIC", 51, 51,
	       &hr, &mn,  &sc,  time_local, ampm             );

    lon_an_dn = fmod( ( hr + mn / 60.0 + sc / 3600.0 ) * 15 , 360);
    fprintf(SC[0].fp, " %f", lon_an_dn);
    // !!! OUTPUT THE LOCAL TIME OF THE SATELLITE (IN DEGREES). UNCOMMENT IF YOU WANT TO OUTPUT IT!
    fprintf(SC[0].fp, " %11.10f %11.10f %11.10f", SC[0].a_i2cg_INRTL[0], SC[0].a_i2cg_INRTL[1], SC[0].a_i2cg_INRTL[2]);
    fprintf(SC[0].fp, " %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e", SC[0].a_i2cg_LVLH[0], SC[0].a_i2cg_LVLH[1], SC[0].a_i2cg_LVLH[2], SC[0].a_i2cg_LVLH_gravity[0], SC[0].a_i2cg_LVLH_gravity[1], SC[0].a_i2cg_LVLH_gravity[2], SC[0].a_i2cg_LVLH_drag[0], SC[0].a_i2cg_LVLH_drag[1], SC[0].a_i2cg_LVLH_drag[2], SC[0].a_i2cg_INRTL_drag[0], SC[0].a_i2cg_INRTL_drag[1], SC[0].a_i2cg_INRTL_drag[2],  SC[0].a_i2cg_INRTL_gravity[0], SC[0].a_i2cg_INRTL_gravity[1], SC[0].a_i2cg_INRTL_gravity[2]);

    fprintf(SC[0].fp, " %f", SC[0].INTEGRATOR.beta_angle*RAD2DEG);


    fprintf(SC[0].fp, " %f", SC[0].OE.an_to_sc*RAD2DEG);


    char orbnumber[10];
    strcpy(orbnumber, "");
    //if (SC->OE.ave_increm == 1){ // means one full orbit has been travelled

    if  ( ( (SC[0].et - SC[0].et_last_orbit ) < OPTIONS->dt_output ) ){
    fprintf(SC[0].fp, " %f", SC[0].OE.w_ave*RAD2DEG);
    fprintf(SC[0].fp, " %f", SC[0].OE.sma_ave);
    fprintf(SC[0].fp, " %f", SC[0].OE.ecc_ave);
    }
    else{
    fprintf(SC[0].fp, " 9999.999999"); // w_ave
    fprintf(SC[0].fp, " 9999.999999"); // sma_ave
    fprintf(SC[0].fp, " 9999.999999"); // ecc_ave

    }
    fprintf(SC[0].fp, " %f", SC[0].OE.zenith*RAD2DEG);
    if  ( ( (SC[0].et - SC[0].et_last_orbit ) < OPTIONS->dt_output ) ){
      fprintf(SC[0].fp, " ORB");
      strcpy(orbnumber, "");
      sprintf(orbnumber, "%d", SC[0].orbit_number);
      fprintf(SC[0].fp, "%s", orbnumber);
    }

    if ( ( SC[0].GEODETIC.latitude  * RAD2DEG < 0.15 ) && ( SC[0].GEODETIC.latitude  * RAD2DEG > -0.15 ) ){
      spkez_c(10, SC[0].et, "J2000", "NONE", 399, x, &lt); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration.

      et2lst_c ( SC[0].et,  399,  SC[0].GEODETIC.longitude, "PLANETOCENTRIC", 51, 51,
		 &hr, &mn,  &sc,  time_local, ampm             );

      lon_an_dn = fmod( ( hr + mn / 60.0 + sc / 3600.0 ) * 15 , 360);

      if ( SC[0].GEODETIC.latitude  * RAD2DEG - previous_lat > 0 ){
	fprintf(SC[0].fp, " %s %s", time_local, "AN");
      }
      else{
	fprintf(SC[0].fp, " %5.2f %s", lon_an_dn, "DN");
	et2utc_c(SC[0].et, "ISOC" ,0 ,255 , times);
      }
    }


    fprintf(SC[0].fp, "\n");

    }
    /* Write the results in the ECEF file */
    /* char text_ecef[256];  */

    /* et2utc_c(SC[0].et, "ISOC" ,0 ,255 , times);     */
    /* sscanf(times, "%4[^\n]",text_ecef); */
    /* strncat(text_ecef, "/",1); */
    /* strncat(text_ecef, &times[5],1); */
    /* strncat(text_ecef, &times[6],1); */
    /* strncat(text_ecef, "/",1); */
    /* strncat(text_ecef, &times[8],1); */
    /* strncat(text_ecef, &times[9],1); */
    /* strncat(text_ecef, " ",1); */

    /* for (ii = 0; ii<12; ii++){ */
    /*   strncat(text_ecef, &times[11+ii],1); */
    /* } */

    //  double T_J2000_to_ECEF[3][3];
    double geodetic[3];
    if (init_flag == 1){
      estate[0] = SC[0].r_i2cg_INRTL[0];estate[1] = SC[0].r_i2cg_INRTL[1];estate[2] = SC[0].r_i2cg_INRTL[2];
      estate[3] = SC[0].v_i2cg_INRTL[0];estate[4] = SC[0].v_i2cg_INRTL[1];estate[5] = SC[0].v_i2cg_INRTL[2];
      sxform_c (  "J2000", earth_fixed_frame,  SC[0].et,    xform  );
      mxvg_c   (  xform,       estate,   6,  6, jstate );
      SC[0].r_ecef2cg_ECEF[0] = jstate[0]; SC[0].r_ecef2cg_ECEF[1] = jstate[1]; SC[0].r_ecef2cg_ECEF[2] = jstate[2];
      SC[0].v_ecef2cg_ECEF[0] = jstate[3]; SC[0].v_ecef2cg_ECEF[1] = jstate[4]; SC[0].v_ecef2cg_ECEF[2] = jstate[5];

/*    // Update Geodetic State */
/* 	eci2lla(SC[0].r_i2cg_INRTL , SC[0].et, geodetic ); */
	
/* 	SC[0].GEODETIC.altitude = geodetic[2]; */
/* 	SC[0].GEODETIC.latitude = geodetic[0]; */
/* 	SC[0].GEODETIC.longitude = geodetic[1];  */

/* 	// update the planet fixed state		 */
/*   double earth_flattening    = 1/298.257223560; */
/*   double earth_radius        = 6378.137; */

/* 	  geodetic_to_geocentric(earth_flattening,             */
/* 				 SC[0].GEODETIC.altitude, */
/* 				 SC[0].GEODETIC.latitude, */
/* 				 SC[0].GEODETIC.longitude, */
/* 				 earth_radius,        */
/* 				 SC[0].r_ecef2cg_ECEF) ; */
    }
    if (write_ecef == 1){
    fprintf(SC[0].fpecef, "%s", text_without_millisecond);  //  fprintf(SC[0].fpecef, "%s", text_ecef );
    fprintf(SC[0].fpecef, " " );
    fprintf(SC[0].fpecef, "%.10f %.10f %.10f %.10f %.10f %.10f\n",
	    SC[0].r_ecef2cg_ECEF[0],
	    SC[0].r_ecef2cg_ECEF[1],
	    SC[0].r_ecef2cg_ECEF[2],
	    SC[0].v_ecef2cg_ECEF[0],
	    SC[0].v_ecef2cg_ECEF[1],
	    SC[0].v_ecef2cg_ECEF[2]);
    }

    /* write the results in the density file */
    if (write_rho == 1){
      fprintf(SC[0].fprho,"%s %e %f %f %e %e\n", text, SC[0].density_here, SC[0].INTEGRATOR.Ta, SC[0].INTEGRATOR.cd_tot_norm, SC[0].INTEGRATOR.A_ref_tot/1000000., SC[0].INTEGRATOR.sum_cd_a_cos * 1000000); // in output A_ref_tot in km^2 */
    }

    /* Write the results in the LLA format file */
    if(  write_lla ==1 ){
    fprintf(SC[0].fpout, "%s", text );
    fprintf(SC[0].fpout, " " );

    // !!! DO NOT UNCOMMENT IF USED WITH STORMS !!!!!
    // !!!!!!!!!!!! CONVERSION BELOW TO ERASE !!!!!!!! DO NOT UNCOMMENT IF USED WITH STORMS !!!!!!!
    /* if ( SC[0].GEODETIC.longitude * RAD2DEG > 180.0 ){ */
    /*   SC[0].GEODETIC.longitude = SC[0].GEODETIC.longitude - 2 * M_PI; */
    /* } */
    // !!!!!!!!!!!! END OF CONVERSION BELOW TO ERASE

    fprintf(SC[0].fpout, "%6.3f %6.3f %6.3f \n",
	    SC[0].GEODETIC.longitude * RAD2DEG,
	    SC[0].GEODETIC.latitude  * RAD2DEG,
	    SC[0].GEODETIC.altitude);
    }


        if(  write_attitude ==1 ){
    fprintf(SC[0].fpatt, "%s", text );
    fprintf(SC[0].fpatt, " " );
    fprintf(SC[0].fpatt, "%6.3f %6.3f %6.3f %d %d %d \n",
	    SC[0].INTEGRATOR.attitude.pitch_current,
	    SC[0].INTEGRATOR.attitude.roll_current,
	    SC[0].INTEGRATOR.attitude.yaw_current,
	    SC[0].INTEGRATOR.attitude.order_pitch_current,
	    SC[0].INTEGRATOR.attitude.order_roll_current,
	    SC[0].INTEGRATOR.attitude.order_yaw_current	    );
    }
  }
  } // end of if this iProc runs main sc sc_index
  //  } // end of if iproc == 0

  if (sc_index < n_sc - n_gps){ // if main sc is not a GPS

    if ((init_flag != 1 ) || (compute_collisions != 1)){
      if (write_ensembles == 1){
	char text[256];
	et2utc_c(constellation_et, "ISOC" ,0 ,255 , times);

	sscanf(times, "%4[^\n]",text);
	strncat(text, "/",1);
	strncat(text, &times[5],1);
	strncat(text, &times[6],1);
	strncat(text, "/",1);
	strncat(text, &times[8],1);
	strncat(text, &times[9],1);
	strncat(text, " ",1);
	for (ii = 0; ii<8; ii++){
	  strncat(text, &times[11+ii],1);
	}



	/* Write the results in the ensemble file */
	// ensemble in COE
	if ( nb_ensembles_min > 0 ){
	  if (( iDebugLevel >= 3 ) ){
	    printf("---- (write_output) Writing results of ensembles.(iProc %d)\n", iProc);
	  }
	  if ( array_sc[1] > 0 )  { // if this iProc runs ensembles (otherwise array_sc[1] = - 1)
	    for (fff = 0; fff < nb_ensembles_output_files ; fff ++){
	      if ( (strcmp(OPTIONS->filename_output_ensemble[fff], "tca" )!= 0 ) && (strcmp(OPTIONS->filename_output_ensemble[fff], "dca" )!= 0 ) && (strcmp(OPTIONS->filename_output_ensemble[fff], "sample" )!= 0 ) ){
		//		if ( ( strcmp( name_file_ensembles[fff], "cd" ) != 0 ) && ( strcmp( name_file_ensembles[fff], "srp" ) != 0 )){// || ( ( strcmp( name_file_ensembles[fff], "cd" ) == 0 ) && ( constellation_et < et_intial_epoch + OPTIONS->dt_output/2. ) ) ) {
		fprintf(SC[0].fpiproc[fff], "\n%s", text );
		fprintf(SC[0].fpiproc[fff], " " );
		//		}
		if ( strcmp( name_file_ensembles[fff], "x_eci" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].r_i2cg_INRTL[0]);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "y_eci" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].r_i2cg_INRTL[1]);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "z_eci" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].r_i2cg_INRTL[2]);
		    //	    	    printf("%.9f\n", SC[eee + iProc * nb_ensemble_min_per_proc].r_i2cg_INRTL[2]);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "vx_eci" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].v_i2cg_INRTL[0]);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "vy_eci" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].v_i2cg_INRTL[1]);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "vz_eci" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].v_i2cg_INRTL[2]);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "latitude" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].GEODETIC.latitude * RAD2DEG);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "longitude" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].GEODETIC.longitude * RAD2DEG);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "altitude" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].GEODETIC.altitude);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "power" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    total_power = 0.0;
		    for (sss = 0; sss < SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.nb_surfaces; sss++){ //SC[0].INTEGRATOR.nb_surfaces
		      total_power = total_power + SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.surface[sss].power_per_surface;
		    }
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    total_power);
		  }
		}
		// // !!!!!!!!!!!! OUTPUT THE ATTITUDE FOR THE ENSEMBLES ASSUMES THAT EITHER ONE OF THESE WAS CHOSEN IN THE INPUT FILE: SOLAR POWER, SOLAR PRESSURE, DRAG, COVERAGE GROUND STATION
		else if ( strcmp( name_file_ensembles[fff], "pitch" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.attitude.pitch_current);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "roll" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.attitude.roll_current);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "yaw" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.attitude.yaw_current);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "sma" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].OE.sma);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "inclination" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].OE.inclination*RAD2DEG);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "eccentricity" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].OE.eccentricity);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "true_anomaly" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].OE.f*RAD2DEG);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "RAAN" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].OE.long_an*RAD2DEG);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "argument_perigee" ) == 0 ){
		  for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		    fprintf(SC[0].fpiproc[fff], "%.9f ",
			    SC[eee + iProc * nb_ensemble_min_per_proc].OE.w*RAD2DEG);
		  }
		}
		else if ( strcmp( name_file_ensembles[fff],"rho" ) == 0 ){
		  if (OPTIONS->include_drag == 1){
		    //	    print_test();
		    for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		      fprintf(SC[0].fpiproc[fff], "%.9f ",
			      SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.density[SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.index_in_driver_interpolated]);
		    }
		  }
		}

		else if ( strcmp( name_file_ensembles[fff], "f107" ) == 0 ){
		  if (OPTIONS->include_drag == 1){
		    //	    print_test();
		    for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		      if ( strcmp(OPTIONS->format_density_driver, "static") == 0 ){
			fprintf(SC[0].fpiproc[fff], "%.9f ",
				SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.f107_static);
		      }
		      else{
			fprintf(SC[0].fpiproc[fff], "%.9f ",
				SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.f107[SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.index_in_driver_interpolated]);

		      }
		    }
		    //	    print_test();
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "f107a" ) == 0 ){
		  if (OPTIONS->include_drag == 1){
		    for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		      if ( strcmp(OPTIONS->format_density_driver, "static") == 0 ){
			fprintf(SC[0].fpiproc[fff], "%.9f ",
				SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.f107A_static);
		      }
		      else{
			fprintf(SC[0].fpiproc[fff], "%.9f ",
				SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.f107A[SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.index_in_driver_interpolated]);
		      }
		    }
		  }
		}
		else if ( strcmp( name_file_ensembles[fff], "ap" ) == 0 ){
		  if (OPTIONS->include_drag == 1){
		    for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		      if ( strcmp(OPTIONS->format_density_driver, "static") == 0 ){
			fprintf(SC[0].fpiproc[fff], "%.9f ",
				SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.Ap_static);
		      }
		      else{
			fprintf(SC[0].fpiproc[fff], "%.9f ",
				SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.Ap[SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.index_in_driver_interpolated]);
		      }
		    }
		  }
		}


		else if ( strcmp( name_file_ensembles[fff], "cd" ) == 0 ){
		  if (OPTIONS->include_drag == 1){
		    //		    if ( constellation_et < et_intial_epoch + OPTIONS->dt_output/2. ){
/* 		    if (iProc == 2){ */
/* 		       printf("<%s>\n", times); */
/* 		    } */

//		    fprintf(SC[0].fpiproc[fff], "\n");
		      if (strcmp(OPTIONS->type_orbit_initialisation, "collision_vcm" ) == 0 ){ // if collision with vcm then the bc is input rom the VCM (or CDM). here not a cd but a bc.
		    
			//	if (SC[0].already_output_cd_ensemble == 0){
		    for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		      fprintf(SC[0].fpiproc[fff], "%.9e ",
			      SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.bc_vcm* 1000. * 1000); // in m2/kg  

		    }
		    //			}
			SC[0].already_output_cd_ensemble = 1;

		      }
		      else{

		       
		    for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		      fprintf(SC[0].fpiproc[fff], "%.9e ",
			      SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.surface[0].Cd); // only the first surface

		    }
		      }
		      //		  }
		}
		}

		else if ( strcmp( name_file_ensembles[fff], "srp" ) == 0 ){
		  if (OPTIONS->include_solar_pressure == 1){
		    //		    if ( constellation_et < et_intial_epoch + OPTIONS->dt_output/2. ){
/* 		    if (iProc == 2){ */
/* 		       printf("<%s>\n", times); */
/* 		    } */
//		    fprintf(SC[0].fpiproc[fff], "\n");

		      if (strcmp(OPTIONS->type_orbit_initialisation, "collision_vcm" ) == 0 ){ // if collision with vcm then the bc is input rom the VCM (or CDM). here not a cd but a bc.
			//	if (SC[0].already_output_srp_ensemble == 0){
		    for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		      fprintf(SC[0].fpiproc[fff], "%.9e ",
			      SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.srp_vcm * 1000. * 1000); // in m2/kg

		    }
		    //			}
			SC[0].already_output_srp_ensemble = 1;

		      }
		      else{
		    for (eee = 1L; eee< nb_ensemble_min_per_proc + 1; eee++){
		      fprintf(SC[0].fpiproc[fff], "%.9e ",
			      SC[eee + iProc * nb_ensemble_min_per_proc].INTEGRATOR.surface[0].solar_radiation_coefficient); // only the first surface

		    }
		      }
		      //		  }
		}
		}

	      } // end of if ensemble to plot is not dca or tca
	    } // end of going through all ensemble files

	  }   // end of if this iProc runs ensembles (otherwise array_sc[1] = - 1)
	}

	if (( iDebugLevel >= 3 ) ){
	  printf("---- (write_output) Done results of ensembles.(iProc %d)\n", iProc);
	}

      }
    }
  
    // ensembles in attitude
    /* if ( nb_ensembles_attitude > 0 ){ */
    /*   for (fff = 0; fff < ( sizeof(name_file_ensembles_attitude[256])/sizeof(name_file_ensembles_attitude[256][0]) ) ; fff ++){ */
    /* 	fprintf(SC[0].fpiproc_attitude[iProc][fff], "\n%s", text ); */
    /* 	fprintf(SC[0].fpiproc_attitude[iProc][fff], " " ); */
    /* 	if ( strcmp( name_file_ensembles_attitude[fff], "power" ) == 0 ){ */
    /* 	  for (eee = 1L; eee< nb_ensemble_attitude_per_proc + 1; eee++){ */
    /* 	    total_power = 0.0; */
    /* 	    for (sss = 0; sss < SC[eee + iProc * nb_ensemble_attitude_per_proc].INTEGRATOR.nb_surfaces; sss++){ //SC[0].INTEGRATOR.nb_surfaces */
    /* 	      total_power = total_power + SC[eee + iProc * nb_ensemble_attitude_per_proc].INTEGRATOR.surface[sss].power_per_surface; */
    /* 	    } */
    /* 	    fprintf(SC[0].fpiproc_attitude[iProc][fff], "%f ", */
    /* 		    total_power); */
    /* 	  } */
    /* 	} */
    /*   } */
    /* } */
    //    if (iProc == 0){
    if (start_ensemble[sc_index] == 0){ // if this iProc runs main sc sc_index
      if (write_reference_sc == 1){
	/* Write the results on the power in the power file */
	if (SC[0].INTEGRATOR.solar_cell_efficiency != -1){
	  fprintf(SC[0].fpower, "%s", text );
	  fprintf(SC[0].fpower, " " );
	  for (sss = 0; sss < SC[0].INTEGRATOR.nb_surfaces; sss++){ //SC[0].INTEGRATOR.nb_surfaces
	    fprintf(SC[0].fpower,"%5.2f ",
		    SC[0].INTEGRATOR.surface[sss].power_per_surface);
	  }
      
	  fprintf(SC[0].fpower,"%s ",
		  SC[0].INTEGRATOR.shadow);

	  fprintf(SC[0].fpower,"%f ",
		  SC[0].INTEGRATOR.sun_elevation*RAD2DEG);

	  fprintf(SC[0].fpower, "\n");

	  // sc in elcipse of Moon if sc in shadow of Moon and sc not in shadow of Earth
	  if ( ( strcmp(SC[0].INTEGRATOR.shadow, "light") == 0 ) && ( strcmp(SC[0].INTEGRATOR.shadow_moon, "light") != 0 ) ){
	    fprintf(SC[0].fpeclipse,"%s\n", text);	    
	  }

/* 	  if ( strcmp(SC[0].INTEGRATOR.shadow_moon, "light") != 0 ) { */
/* 	    printf("%d: %s\n",sc_index, text ); */
/* 	  } */

	}
      }
    } // end of if this iProc runs main sc sc_index
    //        } // end of iproc == 0
  }  // end of if main sc is not a GPS






  if ( (iProc == 0) && ( iDebugLevel >= 2 ) ){
    printf("-- (write_output) Just got out of write_output.(iProc %d)\n", iProc);
  }



  return 0;
    
}

/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           interpolate_position
//  Purpose:        Interpolate the positions of a satellite
//  Assumptions:    None
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 10/04/2015    |   ---     | Initial implementation
//      | C. Bussy-Virat| 10/14/2015    |   ---     | Modifications of inputs
//
/////////////////////////////////////////////////////////////////////////////////////////
int interpolate_position(   char filenamella[300],
			    char filenameecef[300],
                            double          dt_before,
                            double          dt_after, 
			    int             sat_number,
			    OPTIONS_T *OPTIONS,
			    PARAMS_T *PARAMS)

{

  /* Declarations */
  double heading; // heading
  FILE *fplla = NULL, *fpecef = NULL;
  double earth_flattening    = 1/298.257223560;
  double earth_radius        = 6378.137;
  char name_sat[300];
  char filename_position_interpolated[300];
  FILE *fp_position_interpolated;
  char text_name_interpolated[300];
  double ecef_x_previous_step = 0, ecef_y_previous_step = 0, ecef_z_previous_step = 0;
    double ecef_vx_previous_step = 0, ecef_vy_previous_step = 0, ecef_vz_previous_step = 0; // heading
  double et_to_interpolate;
  char *line_ecef = NULL;
  size_t len_ecef = 0;
  char *line_lla = NULL;
  size_t len_lla = 0;
  ssize_t read;
  double lla_to_interpolate[3];
  double ecef_to_interpolate[3];
    double ecef_v_to_interpolate[3]; // heading
  char time_to_interpolate[300];
  double lla_interpolated[3];
  double ecef_interpolated[3];
    double ecef_v_interpolated[3];// heading 
  double time_for_interpolation=0;
  double save_previous_et=0;
  int line_number = 0;
  char times[300];
  char *next;
  char text_output[300], text[300];
  int find_file_name;
  int found_eoh; 

  strcpy(name_sat,OPTIONS->filename_output[sat_number]);
  next = &name_sat[0];
  strcpy(text_output," ");
  find_file_name =  (int)(strchr(next, '.') - next);
  strncat(text_output, next, find_file_name);
  //  printf("Interpolating the positions of%s...\n", text_output);

  /* Algorithm */
  // Read satellite position file
  fplla = fopen(filenamella, "r");
  if (fplla == NULL){
    printf("***! Could not find the file %s. The program will stop. !***\n", filenamella); MPI_Finalize();exit(0);
  }


  fpecef= fopen(filenameecef, "r");

  if (fpecef == NULL){
    printf("***! Could not find the file %s. The program will stop. !***\n", filenameecef); MPI_Finalize();exit(0);
  }



  // Skip the header ECEF
  found_eoh = 0;
  while ( found_eoh == 0 ) {
    getline(&line_ecef, &len_ecef, fpecef);
    sscanf(line_ecef, "%s", text);
    if (  strcmp( "#START", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  // Skip the header LLA
  found_eoh = 0;
  while ( found_eoh == 0 ) {
    getline(&line_lla, &len_lla, fplla);
    sscanf(line_lla, "%s", text);
    if (  strcmp( "#START", text  ) == 0 )  {
      found_eoh = 1;
    }
  }

  strcpy(text_name_interpolated, OPTIONS->dir_output_run_name_sat_name[sat_number]);
  strcat(text_name_interpolated, "/");
  strcat(text_name_interpolated,"interpolated_position_LLA_ECEF_");
  strcat(text_name_interpolated, name_sat);
  strcpy(filename_position_interpolated, text_name_interpolated);
  fp_position_interpolated = fopen(filename_position_interpolated, "w+");

  while ( (read = getline(&line_lla, &len_lla,  fplla)) != -1 ) { // ecef and lla files have the same number of time steps
    getline(&line_ecef, &len_ecef, fpecef);
    sscanf(line_lla, "%19[^\n] %lf %lf %lf", time_to_interpolate, &lla_to_interpolate[0], &lla_to_interpolate[1], &lla_to_interpolate[2]);
    //    sscanf(line_ecef, "%19[^\n] %lf %lf %lf", time_to_interpolate, &ecef_to_interpolate[0], &ecef_to_interpolate[1], &ecef_to_interpolate[2]);
        sscanf(line_ecef, "%19[^\n] %lf %lf %lf %lf %lf %lf", time_to_interpolate, &ecef_to_interpolate[0], &ecef_to_interpolate[1], &ecef_to_interpolate[2], &ecef_v_to_interpolate[0], &ecef_v_to_interpolate[1], &ecef_v_to_interpolate[2]); // heading
    str2et_c(time_to_interpolate, &et_to_interpolate);
    

    if ( line_number > 0 ){
      time_for_interpolation = save_previous_et + dt_after;
      while (time_for_interpolation - 0.00001 <= et_to_interpolate){

	ecef_interpolated[0] = ecef_x_previous_step + ( ( time_for_interpolation - save_previous_et ) / ( et_to_interpolate - save_previous_et )  ) * ( ecef_to_interpolate[0] - ecef_x_previous_step ) ;
	ecef_interpolated[1] = ecef_y_previous_step + ( ( time_for_interpolation - save_previous_et ) / ( et_to_interpolate - save_previous_et ) ) * ( ecef_to_interpolate[1] - ecef_y_previous_step ) ;
	ecef_interpolated[2] = ecef_z_previous_step + ( ( time_for_interpolation - save_previous_et ) / ( et_to_interpolate - save_previous_et ) ) * ( ecef_to_interpolate[2] - ecef_z_previous_step ) ;

	// heading
	ecef_v_interpolated[0] = ecef_vx_previous_step + ( ( time_for_interpolation - save_previous_et ) / ( et_to_interpolate - save_previous_et )  ) * ( ecef_v_to_interpolate[0] - ecef_vx_previous_step ) ;
	ecef_v_interpolated[1] = ecef_vy_previous_step + ( ( time_for_interpolation - save_previous_et ) / ( et_to_interpolate - save_previous_et ) ) * ( ecef_v_to_interpolate[1] - ecef_vy_previous_step ) ;
	ecef_v_interpolated[2] = ecef_vz_previous_step + ( ( time_for_interpolation - save_previous_et ) / ( et_to_interpolate - save_previous_et ) ) * ( ecef_v_to_interpolate[2] - ecef_vz_previous_step ) ;

	// end of heading

	
	geocentric_to_geodetic(
			       ecef_interpolated,
			       &earth_radius,
			       &earth_flattening,
			       &lla_interpolated[2], &lla_interpolated[1], &lla_interpolated[0]);

		// calculate heading of satellite
	compute_heading(&heading, ecef_v_interpolated, lla_interpolated[0], lla_interpolated[1], PARAMS->EARTH.flattening); // heading


	lla_interpolated[2] = lla_interpolated[2];
	lla_interpolated[1] = lla_interpolated[1] *RAD2DEG;
	lla_interpolated[0] = lla_interpolated[0] * RAD2DEG;

	// write in the extrapolated file the interpolated positions
	et2utc_c(time_for_interpolation, "ISOC" ,0 ,255 , times);
		fprintf(fp_position_interpolated, "%s %f %f %f %f %f %f %f \n", times, lla_interpolated[0], lla_interpolated[1], lla_interpolated[2], ecef_interpolated[0], ecef_interpolated[1], ecef_interpolated[2], heading*RAD2DEG); // heading
		//	fprintf(fp_position_interpolated, "%s %f %f %f %f %f %f \n", times, lla_interpolated[0], lla_interpolated[1], lla_interpolated[2], ecef_interpolated[0], ecef_interpolated[1], ecef_interpolated[2]);
	time_for_interpolation = time_for_interpolation + dt_after;

      }
       
      save_previous_et = et_to_interpolate;
      ecef_x_previous_step = ecef_to_interpolate[0];
      ecef_y_previous_step = ecef_to_interpolate[1];
      ecef_z_previous_step = ecef_to_interpolate[2];
      ecef_vx_previous_step = ecef_v_to_interpolate[0];// heading 
      ecef_vy_previous_step = ecef_v_to_interpolate[1];// heading 
      ecef_vz_previous_step = ecef_v_to_interpolate[2];// heading 

    }
    else{
      save_previous_et = et_to_interpolate;
      ecef_x_previous_step = ecef_to_interpolate[0];
      ecef_y_previous_step = ecef_to_interpolate[1];
      ecef_z_previous_step = ecef_to_interpolate[2];
      ecef_vx_previous_step = ecef_v_to_interpolate[0];// heading 
      ecef_vy_previous_step = ecef_v_to_interpolate[1];// heading 
      ecef_vz_previous_step = ecef_v_to_interpolate[2];// heading 


      geocentric_to_geodetic(
			     ecef_to_interpolate,
			     &earth_radius,
			     &earth_flattening,
			     &lla_to_interpolate[2], &lla_to_interpolate[1], &lla_to_interpolate[0]);

	// calculate heading of satellite
	compute_heading(&heading, ecef_v_to_interpolate, lla_to_interpolate[0], lla_to_interpolate[1], PARAMS->EARTH.flattening); // heading

      
      lla_to_interpolate[2] = lla_to_interpolate[2] ;
      lla_to_interpolate[1] = lla_to_interpolate[1] *RAD2DEG;
      lla_to_interpolate[0] = lla_to_interpolate[0] * RAD2DEG;


      et2utc_c(et_to_interpolate, "ISOC" ,0 ,255 , times);     // just to get the same time format as the rest of the file
      fprintf(fp_position_interpolated, "%s %f %f %f %f %f %f %f \n", times, lla_to_interpolate[0], lla_to_interpolate[1], lla_to_interpolate[2],  ecef_to_interpolate[0], ecef_to_interpolate[1], ecef_to_interpolate[2], heading*RAD2DEG); // heading
      
      //      fprintf(fp_position_interpolated, "%s %f %f %f %f %f %f \n", times, lla_to_interpolate[0], lla_to_interpolate[1], lla_to_interpolate[2],  ecef_to_interpolate[0], ecef_to_interpolate[1], ecef_to_interpolate[2]);

    }
    line_number = line_number + 1;
    
  }

  fclose(fplla);
  fclose(fpecef);
  fclose(fp_position_interpolated);


  return 0;

}


int ancas(double min_dist_close_approach, double et_start, double dt_interval, double r1_start[3], double v1_start[3], double a1_start[3], double r1_end[3],double  v1_end[3], double  a1_end[3],double r2_start[3], double v2_start[3], double a2_start[3], double r2_end[3],double  v2_end[3], double  a2_end[3], double gamma0, double gamma1, double gamma2, double gamma3, double *tca1, double *dca1, double *tca2, double *dca2, double *tca3, double *dca3){ // source: Alfvano 1994

  double gamma0_normalized,  gamma1_normalized, gamma2_normalized;
  double root1 = -1e6, root2 = -1e6, root3 = -1e6 ;
  double tca2_before = *tca2, dca2_before = *dca2;
  double tca3_before = *tca3, dca3_before = *dca3;

  // Normalize the coefficients so that the coefficient of order 3 is 1
  gamma0_normalized = gamma0 / gamma3;
  gamma1_normalized = gamma1 / gamma3;
  gamma2_normalized = gamma2 / gamma3;



  //	  printf("I found a min or max distance ||  r2_end[0] = %e\n", r2_end[0]);



  gsl_poly_solve_cubic( gamma2_normalized, gamma1_normalized, gamma0_normalized, &root1, &root2, &root3 );
  //  printf("%e %e %e || r2_end[0] = %e\n",root1, root2, root3, r2_end[0]);
  // test if it's a minimum and not a maximum
  double C_dot_of_f_dot_at_root;
  double found_closest_approach_root1 = 0, found_closest_approach_root2 = 0, found_closest_approach_root3 = 0;
  if (( root1 >=0 ) && ( root1 <= 1 ) ){
    C_dot_of_f_dot_at_root = 3 * gamma3 * root1*root1 + 2*gamma2 * root1 + gamma1 ;
    if ( C_dot_of_f_dot_at_root > 0 ){
      //  ptd(root1, "root1");
      // Find the associated range of closest approach (noted rmin_root1)
      double rmin_root1, q_rdi_root1, q_rdj_root1, q_rdk_root1;
      //      double rmin_root1_fake, q_rdi_root1_fake, q_rdj_root1_fake, q_rdk_root1_fake;
      double alpha0, alpha1, alpha2, alpha3, alpha4, alpha5;
      double dt_interval_square = dt_interval*dt_interval;
      double root1_to_the_square = root1 * root1;
      double root1_to_the_third = root1 * root1 * root1;
      double root1_to_the_fourth = root1 * root1 * root1 * root1 ;
      double root1_to_the_fifth = root1 * root1 * root1 * root1 * root1;
      // // I component
      double frdi_start, frdi_dot_start, frdi_dot_dot_start;
      double frdi_end, frdi_dot_end, frdi_dot_dot_end;
      frdi_start = r2_start[0] - r1_start[0];
      frdi_dot_start = v2_start[0] - v1_start[0];
      frdi_dot_dot_start = a2_start[0] - a1_start[0];
      frdi_end = r2_end[0] - r1_end[0];
      frdi_dot_end = v2_end[0] - v1_end[0];
      frdi_dot_dot_end = a2_end[0] - a1_end[0];

      alpha0 = frdi_start;
      alpha1 = frdi_dot_start * dt_interval;
      alpha2 = 0.5 * frdi_dot_dot_start * dt_interval_square;
      alpha3 = -10 * frdi_start - 6 * frdi_dot_start * dt_interval -1.5 * frdi_dot_dot_start * dt_interval_square 
	+ 10 * frdi_end - 4 * frdi_dot_end * dt_interval + 0.5 * frdi_dot_dot_end * dt_interval_square;
      alpha4 = 15 * frdi_start + 8 * frdi_dot_start * dt_interval + 1.5 * frdi_dot_dot_start * dt_interval_square
	-15 * frdi_end + 7 * frdi_dot_end * dt_interval - frdi_dot_dot_end * dt_interval_square;
      alpha5 = -6 * frdi_start - 3 * frdi_dot_start * dt_interval - 0.5 * frdi_dot_dot_start * dt_interval_square
	+6 * frdi_end - 3 * frdi_dot_end * dt_interval + 0.5 * frdi_dot_dot_end * dt_interval_square;
	  
      q_rdi_root1 = alpha5 * root1_to_the_fifth + alpha4 * root1_to_the_fourth + alpha3 * root1_to_the_third + alpha2 * root1_to_the_square + alpha1 * root1 + alpha0;
      //      q_rdi_root1_fake = alpha5 * root1_to_the_fifth + alpha4 * root1_to_the_fourth + alpha3 * root1_to_the_third + alpha2 * root1_to_the_square + alpha0;

      // // J component
      double frdj_start, frdj_dot_start, frdj_dot_dot_start;
      double frdj_end, frdj_dot_end, frdj_dot_dot_end;
      frdj_start = r2_start[1] - r1_start[1];
      frdj_dot_start = v2_start[1] - v1_start[1];
      frdj_dot_dot_start = a2_start[1] - a1_start[1];
      frdj_end = r2_end[1] - r1_end[1];
      frdj_dot_end = v2_end[1] - v1_end[1];
      frdj_dot_dot_end = a2_end[1] - a1_end[1];

      alpha0 = frdj_start;
      alpha1 = frdj_dot_start * dt_interval;
      alpha2 = 0.5 * frdj_dot_dot_start * dt_interval_square;
      alpha3 = -10 * frdj_start - 6 * frdj_dot_start * dt_interval -1.5 * frdj_dot_dot_start * dt_interval_square 
	+ 10 * frdj_end - 4 * frdj_dot_end * dt_interval + 0.5 * frdj_dot_dot_end * dt_interval_square;
      alpha4 = 15 * frdj_start + 8 * frdj_dot_start * dt_interval + 1.5 * frdj_dot_dot_start * dt_interval_square
	-15 * frdj_end + 7 * frdj_dot_end * dt_interval - frdj_dot_dot_end * dt_interval_square;
      alpha5 = -6 * frdj_start - 3 * frdj_dot_start * dt_interval - 0.5 * frdj_dot_dot_start * dt_interval_square
	+6 * frdj_end - 3 * frdj_dot_end * dt_interval + 0.5 * frdj_dot_dot_end * dt_interval_square;
	  
      q_rdj_root1 = alpha5 * root1_to_the_fifth + alpha4 * root1_to_the_fourth + alpha3 * root1_to_the_third + alpha2 * root1_to_the_square  + alpha1 * root1 + alpha0;
      //      q_rdj_root1_fake = alpha5 * root1_to_the_fifth + alpha4 * root1_to_the_fourth + alpha3 * root1_to_the_third + alpha2 * root1_to_the_square  +  alpha0;

      // // K component
      double frdk_start, frdk_dot_start, frdk_dot_dot_start;
      double frdk_end, frdk_dot_end, frdk_dot_dot_end;
      frdk_start = r2_start[2] - r1_start[2];
      frdk_dot_start = v2_start[2] - v1_start[2];
      frdk_dot_dot_start = a2_start[2] - a1_start[2];
      frdk_end = r2_end[2] - r1_end[2];
      frdk_dot_end = v2_end[2] - v1_end[2];
      frdk_dot_dot_end = a2_end[2] - a1_end[2];

      alpha0 = frdk_start;
      alpha1 = frdk_dot_start * dt_interval;
      alpha2 = 0.5 * frdk_dot_dot_start * dt_interval_square;
      alpha3 = -10 * frdk_start - 6 * frdk_dot_start * dt_interval -1.5 * frdk_dot_dot_start * dt_interval_square 
	+ 10 * frdk_end - 4 * frdk_dot_end * dt_interval + 0.5 * frdk_dot_dot_end * dt_interval_square;
      alpha4 = 15 * frdk_start + 8 * frdk_dot_start * dt_interval + 1.5 * frdk_dot_dot_start * dt_interval_square
	-15 * frdk_end + 7 * frdk_dot_end * dt_interval - frdk_dot_dot_end * dt_interval_square;
      alpha5 = -6 * frdk_start - 3 * frdk_dot_start * dt_interval - 0.5 * frdk_dot_dot_start * dt_interval_square
	+6 * frdk_end - 3 * frdk_dot_end * dt_interval + 0.5 * frdk_dot_dot_end * dt_interval_square;
	  
      q_rdk_root1 = alpha5 * root1_to_the_fifth + alpha4 * root1_to_the_fourth + alpha3 * root1_to_the_third + alpha2 * root1_to_the_square + alpha1 * root1 + alpha0;
      //    q_rdk_root1_fake = alpha5 * root1_to_the_fifth + alpha4 * root1_to_the_fourth + alpha3 * root1_to_the_third + alpha2 * root1_to_the_square + alpha0;
      // // Compute rmin
      rmin_root1 = sqrt( q_rdi_root1 * q_rdi_root1 + q_rdj_root1 * q_rdj_root1 + q_rdk_root1 * q_rdk_root1 );
      //      rmin_root1_fake = sqrt( q_rdi_root1_fake * q_rdi_root1_fake + q_rdj_root1_fake * q_rdj_root1_fake + q_rdk_root1_fake * q_rdk_root1_fake );

      if ( ( rmin_root1 <= min_dist_close_approach ) && ( rmin_root1 >= 0 )){

	/* ptd(rmin_root1, "rmin_root1"); */
	/* ptd(distance_between_two_sc(r1_start, r2_start), "d_start"); */
	/* ptd(distance_between_two_sc(r1_end, r2_end), "d_end"); */

	//	exitall();
	//	ptd(rmin_root1_fake, "rmin_root1_fake");
	*dca1 = rmin_root1;
	// Find the associated time of closest approach (noted t_root1)
	*tca1 = et_start+ root1 * dt_interval;
	
	/* char t_root1_str[300]; */
	/* et2utc_c( t_root1, "ISOC", 3, 256, t_root1_str ); */
	/* printf("%s: %e (root1 | %s - %s)\n", t_root1_str, rmin_root1, et_start_str, et_end_str); */
	/* printf("%e %e %e\n", root1, root2, root3); */
	found_closest_approach_root1 = 1;
      }
    }
  }
  if ( ( root2 >=0 ) && ( root2 <= 1 ) ){
    C_dot_of_f_dot_at_root = 3 * gamma3 * root2*root2 + 2*gamma2 * root2 + gamma1 ;
    if ( C_dot_of_f_dot_at_root > 0 ){
      //          ptd(root2, "root2");
      // Find the associated range of closest approach (noted rmin_root2)
      double rmin_root2, q_rdi_root2, q_rdj_root2, q_rdk_root2;
      double alpha0, alpha1, alpha2, alpha3, alpha4, alpha5;
      double dt_interval_square = dt_interval*dt_interval;
      double root2_to_the_square = root2 * root2;
      double root2_to_the_third = root2 * root2 * root2;
      double root2_to_the_fourth = root2 * root2 * root2 * root2 ;
      double root2_to_the_fifth = root2 * root2 * root2 * root2 * root2;
      // // I component
      double frdi_start, frdi_dot_start, frdi_dot_dot_start;
      double frdi_end, frdi_dot_end, frdi_dot_dot_end;
      frdi_start = r2_start[0] - r1_start[0];
      frdi_dot_start = v2_start[0] - v1_start[0];
      frdi_dot_dot_start = a2_start[0] - a1_start[0];
      frdi_end = r2_end[0] - r1_end[0];
      frdi_dot_end = v2_end[0] - v1_end[0];
      frdi_dot_dot_end = a2_end[0] - a1_end[0];

      alpha0 = frdi_start;
      alpha1 = frdi_dot_start * dt_interval;
      alpha2 = 0.5 * frdi_dot_dot_start * dt_interval_square;
      alpha3 = -10 * frdi_start - 6 * frdi_dot_start * dt_interval -1.5 * frdi_dot_dot_start * dt_interval_square 
	+ 10 * frdi_end - 4 * frdi_dot_end * dt_interval + 0.5 * frdi_dot_dot_end * dt_interval_square;
      alpha4 = 15 * frdi_start + 8 * frdi_dot_start * dt_interval + 1.5 * frdi_dot_dot_start * dt_interval_square
	-15 * frdi_end + 7 * frdi_dot_end * dt_interval - frdi_dot_dot_end * dt_interval_square;
      alpha5 = -6 * frdi_start - 3 * frdi_dot_start * dt_interval - 0.5 * frdi_dot_dot_start * dt_interval_square
	+6 * frdi_end - 3 * frdi_dot_end * dt_interval + 0.5 * frdi_dot_dot_end * dt_interval_square;
	  
      q_rdi_root2 = alpha5 * root2_to_the_fifth + alpha4 * root2_to_the_fourth + alpha3 * root2_to_the_third + alpha2 * root2_to_the_square + alpha1 * root2 + alpha0;

      // // J component
      double frdj_start, frdj_dot_start, frdj_dot_dot_start;
      double frdj_end, frdj_dot_end, frdj_dot_dot_end;
      frdj_start = r2_start[1] - r1_start[1];
      frdj_dot_start = v2_start[1] - v1_start[1];
      frdj_dot_dot_start = a2_start[1] - a1_start[1];
      frdj_end = r2_end[1] - r1_end[1];
      frdj_dot_end = v2_end[1] - v1_end[1];
      frdj_dot_dot_end = a2_end[1] - a1_end[1];

      alpha0 = frdj_start;
      alpha1 = frdj_dot_start * dt_interval;
      alpha2 = 0.5 * frdj_dot_dot_start * dt_interval_square;
      alpha3 = -10 * frdj_start - 6 * frdj_dot_start * dt_interval -1.5 * frdj_dot_dot_start * dt_interval_square 
	+ 10 * frdj_end - 4 * frdj_dot_end * dt_interval + 0.5 * frdj_dot_dot_end * dt_interval_square;
      alpha4 = 15 * frdj_start + 8 * frdj_dot_start * dt_interval + 1.5 * frdj_dot_dot_start * dt_interval_square
	-15 * frdj_end + 7 * frdj_dot_end * dt_interval - frdj_dot_dot_end * dt_interval_square;
      alpha5 = -6 * frdj_start - 3 * frdj_dot_start * dt_interval - 0.5 * frdj_dot_dot_start * dt_interval_square
	+6 * frdj_end - 3 * frdj_dot_end * dt_interval + 0.5 * frdj_dot_dot_end * dt_interval_square;
	  
      q_rdj_root2 = alpha5 * root2_to_the_fifth + alpha4 * root2_to_the_fourth + alpha3 * root2_to_the_third + alpha2 * root2_to_the_square + alpha1 * root2 + alpha0;

      // // K component
      double frdk_start, frdk_dot_start, frdk_dot_dot_start;
      double frdk_end, frdk_dot_end, frdk_dot_dot_end;
      frdk_start = r2_start[2] - r1_start[2];
      frdk_dot_start = v2_start[2] - v1_start[2];
      frdk_dot_dot_start = a2_start[2] - a1_start[2];
      frdk_end = r2_end[2] - r1_end[2];
      frdk_dot_end = v2_end[2] - v1_end[2];
      frdk_dot_dot_end = a2_end[2] - a1_end[2];

      alpha0 = frdk_start;
      alpha1 = frdk_dot_start * dt_interval;
      alpha2 = 0.5 * frdk_dot_dot_start * dt_interval_square;
      alpha3 = -10 * frdk_start - 6 * frdk_dot_start * dt_interval -1.5 * frdk_dot_dot_start * dt_interval_square 
	+ 10 * frdk_end - 4 * frdk_dot_end * dt_interval + 0.5 * frdk_dot_dot_end * dt_interval_square;
      alpha4 = 15 * frdk_start + 8 * frdk_dot_start * dt_interval + 1.5 * frdk_dot_dot_start * dt_interval_square
	-15 * frdk_end + 7 * frdk_dot_end * dt_interval - frdk_dot_dot_end * dt_interval_square;
      alpha5 = -6 * frdk_start - 3 * frdk_dot_start * dt_interval - 0.5 * frdk_dot_dot_start * dt_interval_square
	+6 * frdk_end - 3 * frdk_dot_end * dt_interval + 0.5 * frdk_dot_dot_end * dt_interval_square;
	  
      q_rdk_root2 = alpha5 * root2_to_the_fifth + alpha4 * root2_to_the_fourth + alpha3 * root2_to_the_third + alpha2 * root2_to_the_square + alpha1 * root2 + alpha0;

      // // Compute rmin
      rmin_root2 = sqrt( q_rdi_root2 * q_rdi_root2 + q_rdj_root2 * q_rdj_root2 + q_rdk_root2 * q_rdk_root2 );

      if ( ( rmin_root2 <= min_dist_close_approach ) && ( rmin_root2 >= 0 ) ){

	*dca2 = rmin_root2;
	// Find the associated time of closest approach (noted t_root2)
	*tca2 = et_start+ root2 * dt_interval;
	
	//	char t_root2_str[300];
	//	et2utc_c( t_root2, "ISOC", 3, 256, t_root2_str );
	//	printf("%s: %e (root2 | %s - %s)\n", t_root2_str, rmin_root2, et_start_str, et_end_str);
	///	printf("%e %e %e\n", root1, root2, root3);
	found_closest_approach_root2 = 1;
      }
    }

  }
  if ( ( root3 >=0 ) && ( root3 <= 1 ) ){
    C_dot_of_f_dot_at_root = 3 * gamma3 * root3*root3 + 2*gamma2 * root3 + gamma1 ;
    if ( C_dot_of_f_dot_at_root > 0 ){
      //          ptd(root3, "root3");
      // Find the associated range of closest approach (noted rmin_root3)
      double rmin_root3, q_rdi_root3, q_rdj_root3, q_rdk_root3;
      double alpha0, alpha1, alpha2, alpha3, alpha4, alpha5;
      double dt_interval_square = dt_interval*dt_interval;
      double root3_to_the_square = root3 * root3;
      double root3_to_the_third = root3 * root3 * root3;
      double root3_to_the_fourth = root3 * root3 * root3 * root3 ;
      double root3_to_the_fifth = root3 * root3 * root3 * root3 * root3;
      // // I component
      double frdi_start, frdi_dot_start, frdi_dot_dot_start;
      double frdi_end, frdi_dot_end, frdi_dot_dot_end;
      frdi_start = r2_start[0] - r1_start[0];
      frdi_dot_start = v2_start[0] - v1_start[0];
      frdi_dot_dot_start = a2_start[0] - a1_start[0];
      frdi_end = r2_end[0] - r1_end[0];
      frdi_dot_end = v2_end[0] - v1_end[0];
      frdi_dot_dot_end = a2_end[0] - a1_end[0];

      alpha0 = frdi_start;
      alpha1 = frdi_dot_start * dt_interval;
      alpha2 = 0.5 * frdi_dot_dot_start * dt_interval_square;
      alpha3 = -10 * frdi_start - 6 * frdi_dot_start * dt_interval -1.5 * frdi_dot_dot_start * dt_interval_square 
	+ 10 * frdi_end - 4 * frdi_dot_end * dt_interval + 0.5 * frdi_dot_dot_end * dt_interval_square;
      alpha4 = 15 * frdi_start + 8 * frdi_dot_start * dt_interval + 1.5 * frdi_dot_dot_start * dt_interval_square
	-15 * frdi_end + 7 * frdi_dot_end * dt_interval - frdi_dot_dot_end * dt_interval_square;
      alpha5 = -6 * frdi_start - 3 * frdi_dot_start * dt_interval - 0.5 * frdi_dot_dot_start * dt_interval_square
	+6 * frdi_end - 3 * frdi_dot_end * dt_interval + 0.5 * frdi_dot_dot_end * dt_interval_square;
	  
      q_rdi_root3 = alpha5 * root3_to_the_fifth + alpha4 * root3_to_the_fourth + alpha3 * root3_to_the_third + alpha2 * root3_to_the_square + alpha1 * root3 + alpha0;

      // // J component
      double frdj_start, frdj_dot_start, frdj_dot_dot_start;
      double frdj_end, frdj_dot_end, frdj_dot_dot_end;
      frdj_start = r2_start[1] - r1_start[1];
      frdj_dot_start = v2_start[1] - v1_start[1];
      frdj_dot_dot_start = a2_start[1] - a1_start[1];
      frdj_end = r2_end[1] - r1_end[1];
      frdj_dot_end = v2_end[1] - v1_end[1];
      frdj_dot_dot_end = a2_end[1] - a1_end[1];

      alpha0 = frdj_start;
      alpha1 = frdj_dot_start * dt_interval;
      alpha2 = 0.5 * frdj_dot_dot_start * dt_interval_square;
      alpha3 = -10 * frdj_start - 6 * frdj_dot_start * dt_interval -1.5 * frdj_dot_dot_start * dt_interval_square 
	+ 10 * frdj_end - 4 * frdj_dot_end * dt_interval + 0.5 * frdj_dot_dot_end * dt_interval_square;
      alpha4 = 15 * frdj_start + 8 * frdj_dot_start * dt_interval + 1.5 * frdj_dot_dot_start * dt_interval_square
	-15 * frdj_end + 7 * frdj_dot_end * dt_interval - frdj_dot_dot_end * dt_interval_square;
      alpha5 = -6 * frdj_start - 3 * frdj_dot_start * dt_interval - 0.5 * frdj_dot_dot_start * dt_interval_square
	+6 * frdj_end - 3 * frdj_dot_end * dt_interval + 0.5 * frdj_dot_dot_end * dt_interval_square;
	  
      q_rdj_root3 = alpha5 * root3_to_the_fifth + alpha4 * root3_to_the_fourth + alpha3 * root3_to_the_third + alpha2 * root3_to_the_square + alpha1 * root3 + alpha0;

      // // K component
      double frdk_start, frdk_dot_start, frdk_dot_dot_start;
      double frdk_end, frdk_dot_end, frdk_dot_dot_end;
      frdk_start = r2_start[2] - r1_start[2];
      frdk_dot_start = v2_start[2] - v1_start[2];
      frdk_dot_dot_start = a2_start[2] - a1_start[2];
      frdk_end = r2_end[2] - r1_end[2];
      frdk_dot_end = v2_end[2] - v1_end[2];
      frdk_dot_dot_end = a2_end[2] - a1_end[2];

      alpha0 = frdk_start;
      alpha1 = frdk_dot_start * dt_interval;
      alpha2 = 0.5 * frdk_dot_dot_start * dt_interval_square;
      alpha3 = -10 * frdk_start - 6 * frdk_dot_start * dt_interval -1.5 * frdk_dot_dot_start * dt_interval_square 
	+ 10 * frdk_end - 4 * frdk_dot_end * dt_interval + 0.5 * frdk_dot_dot_end * dt_interval_square;
      alpha4 = 15 * frdk_start + 8 * frdk_dot_start * dt_interval + 1.5 * frdk_dot_dot_start * dt_interval_square
	-15 * frdk_end + 7 * frdk_dot_end * dt_interval - frdk_dot_dot_end * dt_interval_square;
      alpha5 = -6 * frdk_start - 3 * frdk_dot_start * dt_interval - 0.5 * frdk_dot_dot_start * dt_interval_square
	+6 * frdk_end - 3 * frdk_dot_end * dt_interval + 0.5 * frdk_dot_dot_end * dt_interval_square;
	  
      q_rdk_root3 = alpha5 * root3_to_the_fifth + alpha4 * root3_to_the_fourth + alpha3 * root3_to_the_third + alpha2 * root3_to_the_square + alpha1 * root3 + alpha0;

      // // Compute rmin
      rmin_root3 = sqrt( q_rdi_root3 * q_rdi_root3 + q_rdj_root3 * q_rdj_root3 + q_rdk_root3 * q_rdk_root3 );
      if ( ( rmin_root3 <= min_dist_close_approach ) && ( rmin_root3 >= 0 ) ){

	// Find the associated time of closest approach (noted t_root3)
	*dca3 = rmin_root3;
	*tca3 = et_start + root3 * dt_interval;
	
	/* char t_root3_str[300]; */
	/* et2utc_c( t_root3, "ISOC", 3, 256, t_root3_str ); */
	/* printf("%s: %e (root3 | %s - %s)\n", t_root3_str, rmin_root3, et_start_str, et_end_str); */
	found_closest_approach_root3 = 1;
      }
    }
  }


  if ( ( found_closest_approach_root1 == 1 ) || ( found_closest_approach_root2 == 1 ) || ( found_closest_approach_root3 == 1 ) ){

    if ( fabs( root1 - root2 ) < 0.00001 ){ // accuracy if 0.00001 second 
      *tca2 = tca2_before;
      *dca2 = dca2_before;
    }
    if ( fabs( root3 - root1 ) < 0.00001 ){ // accuracy if 0.00001 second 
      *tca3 = tca3_before;
      *dca3 = dca3_before;
    }
    if ( ( *tca2 != tca2_before  ) && (fabs( root2 - root3 ) < 0.00001 ) ){ // accuracy if 0.00001 second
      *tca3 = tca3_before;
      *dca3 = dca3_before      ;
    }
    /* etprint(et_start, "time"); */
    /* printf("%f %f %f\n", root1, root2, root3); */
  }
  
  

  //  printf("OK\n");


  return 0;
}

int close_approach_ensemble( double *eci_x_primary_sc_in_span, double *eci_y_primary_sc_in_span, double *eci_z_primary_sc_in_span, double  *eci_vx_primary_sc_in_span, double *eci_vy_primary_sc_in_span, double *eci_vz_primary_sc_in_span, double *eci_ax_primary_sc_in_span, double *eci_ay_primary_sc_in_span, double *eci_az_primary_sc_in_span,  
double *eci_x_secondary_sc_in_span, double *eci_y_secondary_sc_in_span, double *eci_z_secondary_sc_in_span, double  *eci_vx_secondary_sc_in_span, double *eci_vy_secondary_sc_in_span, double *eci_vz_secondary_sc_in_span, double *eci_ax_secondary_sc_in_span, double *eci_ay_secondary_sc_in_span, double *eci_az_secondary_sc_in_span,  
			     int time_step_of_tca,double *gamma0, double *gamma1, double *gamma2, double *gamma3,  int *min_exists, double r1_start[3], double v1_start[3], double a1_start[3], double r1_end[3],double  v1_end[3], double  a1_end[3],double r2_start[3], double v2_start[3], double a2_start[3], double r2_end[3],double  v2_end[3], double a2_end[3], int *time_step_start_interval, double *min_distance_in_time_spanning_tca, double *direction_distance, int initial_epoch_time_step_in_span, int final_epoch_time_step_in_span, int iProc, double et_time_step_of_save_tca, OPTIONS_T *OPTIONS ){ 
    // close_approach_ensemble determines if there is a close approach (defined as a minimum of the distance) between the primary ensemble sc eee_prim and the secondary ensemble sc eee in the interval of time spanning TCA (unperturbed). If there is, it returns the position, velocity, and acceleration of of sc eee_prim and sc eee at both the beginning and end of interval of time 3 * OPTIONS->dt in which the close approach between sc eee_prim and sc eee is expected to be. It also returns the coefficients gamma that are then used to calculate the min distance using the function ancas (if there is a min distanc)

  // ASSUMPTIONS:
  // - FOR NOW WORKS ONLY IF TWO REFERENCE SATELLIES ONLY  
  // - THE 2 REFERENCE SC (UNPERTURBED ORBITS) HAVE THE SAME EPOCH START.

  //  NOTES:
  // - in arguments of close_approach_ensemble, start means the oldest time of 4 four latest points, and end means the most recent time of 4 four latest points (see Alfano 2009)

  /* Declarations */    
  int nb_time_steps_in_tca_time_span = time_step_of_tca * 2 + 1; // recall that nb_time_steps_in_tca_time_span is odd
    double et_start_of_span = et_time_step_of_save_tca - ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt );

    char time_et_time_step_of_save_tca[256], time_end_epoch[256], time_start_epoch[256];
  double eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary;

  double p1, p2, p3, p4;
  double D;
  double tau1, tau1_square, tau1_cube;
  double tau2, tau2_square, tau2_cube;
  double dist_at_step_plus_n_minus_one_dt, dist_at_step_plus_n_dt;
  double dist_at_step, dist_at_step_minus_one_dt, dist_at_step_minus_two_dt;
  double eci_r_secondary_minus_eci_r_primary[3], eci_v_secondary_minus_eci_v_primary[3];

  /* Algorithm */ // see Alfvano 2009 end of section II (Satellite Conjunction Monte Carlo Analysis) (same notations)

  // the initial epoch needs to start at least three 3 time steps (because of the algorithm) before TCA. So check that it's the case and if it's not the return an error message 
  //  pti(final_epoch_time_step_in_span, "final_epoch_time_step_in_span");
  
  
  if (initial_epoch_time_step_in_span > time_step_of_tca - 2){
    et2utc_c(et_time_step_of_save_tca, "ISOC", 3, 255, time_et_time_step_of_save_tca);
    et2utc_c(et_start_of_span + initial_epoch_time_step_in_span * OPTIONS->dt , "ISOC", 3, 255, time_start_epoch);
    printf("***! The propagation needs to start at least two time steps before TCA. Here the propagation starts at %s, but the first TCA is %s. The program will stop. (iProc: %d) !***\n", time_start_epoch , time_et_time_step_of_save_tca, iProc); MPI_Finalize(); exit(0);
  }

  if ( final_epoch_time_step_in_span < time_step_of_tca ){
    et2utc_c(et_time_step_of_save_tca, "ISOC", 3, 255, time_et_time_step_of_save_tca);
    et2utc_c(et_start_of_span + final_epoch_time_step_in_span * OPTIONS->dt , "ISOC", 3, 255, time_end_epoch);
    printf("***! The propation can not end before TCA (it can end at TCA). Here the propagation ends at %s, but the last TCA is %s. The program will stop. (iProc: %d) !***\n", time_end_epoch , time_et_time_step_of_save_tca, iProc); MPI_Finalize(); exit(0);
  }

  // D(TCA)  
  eci_r_secondary_minus_eci_r_primary[0] = eci_x_secondary_sc_in_span[time_step_of_tca] - eci_x_primary_sc_in_span[time_step_of_tca];
  eci_r_secondary_minus_eci_r_primary[1] = eci_y_secondary_sc_in_span[time_step_of_tca] - eci_y_primary_sc_in_span[time_step_of_tca];
  eci_r_secondary_minus_eci_r_primary[2] = eci_z_secondary_sc_in_span[time_step_of_tca] - eci_z_primary_sc_in_span[time_step_of_tca];
  v_mag(&dist_at_step, eci_r_secondary_minus_eci_r_primary );
  // D(TCA-dt)
  eci_r_secondary_minus_eci_r_primary[0] = eci_x_secondary_sc_in_span[time_step_of_tca-1] - eci_x_primary_sc_in_span[time_step_of_tca-1];
  eci_r_secondary_minus_eci_r_primary[1] = eci_y_secondary_sc_in_span[time_step_of_tca-1] - eci_y_primary_sc_in_span[time_step_of_tca-1];
  eci_r_secondary_minus_eci_r_primary[2] = eci_z_secondary_sc_in_span[time_step_of_tca-1] - eci_z_primary_sc_in_span[time_step_of_tca-1];
  v_mag(&dist_at_step_minus_one_dt, eci_r_secondary_minus_eci_r_primary );
  // D(TCA-2*dt)
  eci_r_secondary_minus_eci_r_primary[0] = eci_x_secondary_sc_in_span[time_step_of_tca-2] - eci_x_primary_sc_in_span[time_step_of_tca-2];
  eci_r_secondary_minus_eci_r_primary[1] = eci_y_secondary_sc_in_span[time_step_of_tca-2] - eci_y_primary_sc_in_span[time_step_of_tca-2];
  eci_r_secondary_minus_eci_r_primary[2] = eci_z_secondary_sc_in_span[time_step_of_tca-2] - eci_z_primary_sc_in_span[time_step_of_tca-2];
  v_mag(&dist_at_step_minus_two_dt, eci_r_secondary_minus_eci_r_primary );

  /* ptd(dist_at_step, "dist_at_step"); */
  /* ptd(dist_at_step_minus_one_dt,"dist_at_step_minus_one_dt"); */
  /* ptd(dist_at_step_minus_two_dt,"dist_at_step_minus_two_dt"); */
  // Determine direction
  *direction_distance = dist_at_step_minus_two_dt - dist_at_step_minus_one_dt;
  
  // Move backward if *direction_distance < 0, move forward if *direction_distance > 0
  int istep;
  int found_greater_distance = 0;
  double dist_at_step_minus_n_dt, dist_at_step_minus_n_minus_one_dt;
  dist_at_step_minus_n_minus_one_dt = dist_at_step_minus_two_dt;
  dist_at_step_minus_n_dt = dist_at_step_minus_two_dt; // this will be overwritten in the loop below. But in the case where time_step_of_tca-3 <  initial_epoch_time_step_in_span (so WHERE TIME_STEP_OF_TCA = INITIAL_EPOCH_TIME_STEP_IN_SPAN + 2 because time_step_of_tca can't be < initial_epoch_time_step_in_span + 2 because of if condiction at beginning of this function), we don't get in the loop so dist_at_step_minus_n_dt has to be equal to the distance at time_step_of_tca - 2, which is the distance at initial epoch. This is a border limit case. 
  if ( *direction_distance <= 0 ){
    istep = 3;
    while ( ( istep < time_step_of_tca ) && ( found_greater_distance == 0 ) && ( time_step_of_tca-istep >=  initial_epoch_time_step_in_span ) ){ // go bakward in time until either finding a greater distance or reaching the beginning of the time span (ie istep = time_step_of_tca)
      eci_r_secondary_minus_eci_r_primary[0] = eci_x_secondary_sc_in_span[time_step_of_tca-istep] - eci_x_primary_sc_in_span[time_step_of_tca-istep];
      eci_r_secondary_minus_eci_r_primary[1] = eci_y_secondary_sc_in_span[time_step_of_tca-istep] - eci_y_primary_sc_in_span[time_step_of_tca-istep];
      eci_r_secondary_minus_eci_r_primary[2] = eci_z_secondary_sc_in_span[time_step_of_tca-istep] - eci_z_primary_sc_in_span[time_step_of_tca-istep];
      v_mag(&dist_at_step_minus_n_dt, eci_r_secondary_minus_eci_r_primary );
      if ( dist_at_step_minus_n_dt >= dist_at_step_minus_n_minus_one_dt ){
	found_greater_distance = 1;
      }
      /* ptd(dist_at_step_minus_n_dt, "dist_at_step_minus_n_dt"); */
    
      dist_at_step_minus_n_minus_one_dt = dist_at_step_minus_n_dt;
      istep = istep + 1;
    }
    istep = istep - 1;
  }

  else{
    istep = 1;
    dist_at_step_plus_n_minus_one_dt = dist_at_step;
    dist_at_step_plus_n_dt = dist_at_step; // this will be overwritten in the loop below. But in the case where time_step_of_tca+1 >  final_epoch_time_step_in_span (so WHERE TIME_STEP_OF_TCA + 1 = FINAL_EPOCH_TIME_STEP_IN_SPAN + 1 (so TIME_STEP_OF_TCA = FINAL_EPOCH_TIME_STEP_IN_SPAN) because time_step_of_tca + 1 can't be > final_epoch_time_step_in_span + 1 (because time_step_of_tca can't be > final_epoch_time_step_in_span) because of if condiction at beginning of this function), we don't get in the loop so dist_at_step_plus_n_dt has to be equal to the distance at time_step_of_tca, which is the distance at final epoch. This is a border limit case. 
    while ( ( istep < time_step_of_tca ) && ( found_greater_distance == 0 ) && (time_step_of_tca+istep <= final_epoch_time_step_in_span) ){ // go foward in time until either finding a greater distance or reaching the end of the time span (ie istep = time_step_of_tca)
      eci_r_secondary_minus_eci_r_primary[0] = eci_x_secondary_sc_in_span[time_step_of_tca+istep] - eci_x_primary_sc_in_span[time_step_of_tca+istep];
      eci_r_secondary_minus_eci_r_primary[1] = eci_y_secondary_sc_in_span[time_step_of_tca+istep] - eci_y_primary_sc_in_span[time_step_of_tca+istep];
      eci_r_secondary_minus_eci_r_primary[2] = eci_z_secondary_sc_in_span[time_step_of_tca+istep] - eci_z_primary_sc_in_span[time_step_of_tca+istep];

      v_mag(&dist_at_step_plus_n_dt, eci_r_secondary_minus_eci_r_primary );
      if ( dist_at_step_plus_n_dt >= dist_at_step_plus_n_minus_one_dt ){
	found_greater_distance = 1;
      }
      //           ptd(dist_at_step_plus_n_dt, "dist_at_step_plus_n_dt positve");
      dist_at_step_plus_n_minus_one_dt = dist_at_step_plus_n_dt;
      istep = istep + 1;
    }
    istep = istep - 1;
    //    pti(time_step_of_tca, "time_step_of_tca");
      /*  if (istep == time_step_of_tca  -1){ */
      /* printf("%d %d %d\n", found_greater_distance, istep, time_step_of_tca * 2); */
      /*   } */
  }

  /////////////////////////////////////////////////////////////
  ///////////////// INITIALIZE ANCAS //////////////////////////
  // By initialioze ANCS, we mean calculate the gamma coefficients (see Alfano 1994) that are then used in the function ancas (we call the cunfction ancas only if we know in this function here that there is a minimum (see test at the end of this function))
  // If found_greater_distance = 1, then use the last four consecutive distances in ANCAS algorithm to determine minimum distance
  if ( found_greater_distance == 1 ){ // if a greater distance have been found, then compute minimum distance with ANCAS
    // The reason we don't use the function ancas is because we use 4 values of the function to determine the min, while the function ancas uses 2 values of the function and 2 values of the derivative of the function (we actually could do that too but we want to use a similar algorithm to Alfano 2009 to then compare the results)

    // Four last points
    if ( *direction_distance <= 0 ){
      eci_r_secondary_minus_eci_r_primary[0] = eci_x_secondary_sc_in_span[time_step_of_tca-istep] - eci_x_primary_sc_in_span[time_step_of_tca-istep];
      eci_r_secondary_minus_eci_r_primary[1] = eci_y_secondary_sc_in_span[time_step_of_tca-istep] - eci_y_primary_sc_in_span[time_step_of_tca-istep];
      eci_r_secondary_minus_eci_r_primary[2] = eci_z_secondary_sc_in_span[time_step_of_tca-istep] - eci_z_primary_sc_in_span[time_step_of_tca-istep];

      eci_v_secondary_minus_eci_v_primary[0] = eci_vx_secondary_sc_in_span[time_step_of_tca-istep] - eci_vx_primary_sc_in_span[time_step_of_tca-istep];
      eci_v_secondary_minus_eci_v_primary[1] = eci_vy_secondary_sc_in_span[time_step_of_tca-istep] - eci_vy_primary_sc_in_span[time_step_of_tca-istep];
      eci_v_secondary_minus_eci_v_primary[2] = eci_vz_secondary_sc_in_span[time_step_of_tca-istep] - eci_vz_primary_sc_in_span[time_step_of_tca-istep];

      v_dot( &eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary, eci_r_secondary_minus_eci_r_primary, eci_v_secondary_minus_eci_v_primary );
      p1 = 2 * eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary;

      // p2
      eci_r_secondary_minus_eci_r_primary[0] = eci_x_secondary_sc_in_span[time_step_of_tca-istep + 1] - eci_x_primary_sc_in_span[time_step_of_tca-istep + 1];
      eci_r_secondary_minus_eci_r_primary[1] = eci_y_secondary_sc_in_span[time_step_of_tca-istep + 1] - eci_y_primary_sc_in_span[time_step_of_tca-istep + 1];
      eci_r_secondary_minus_eci_r_primary[2] = eci_z_secondary_sc_in_span[time_step_of_tca-istep + 1] - eci_z_primary_sc_in_span[time_step_of_tca-istep + 1];

      eci_v_secondary_minus_eci_v_primary[0] = eci_vx_secondary_sc_in_span[time_step_of_tca-istep + 1] - eci_vx_primary_sc_in_span[time_step_of_tca-istep + 1];
      eci_v_secondary_minus_eci_v_primary[1] = eci_vy_secondary_sc_in_span[time_step_of_tca-istep + 1] - eci_vy_primary_sc_in_span[time_step_of_tca-istep + 1];
      eci_v_secondary_minus_eci_v_primary[2] = eci_vz_secondary_sc_in_span[time_step_of_tca-istep + 1] - eci_vz_primary_sc_in_span[time_step_of_tca-istep + 1];

      v_dot( &eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary, eci_r_secondary_minus_eci_r_primary, eci_v_secondary_minus_eci_v_primary );
      p2 = 2 * eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary;
 
      // p3
      eci_r_secondary_minus_eci_r_primary[0] = eci_x_secondary_sc_in_span[time_step_of_tca-istep + 2] - eci_x_primary_sc_in_span[time_step_of_tca-istep + 2];
      eci_r_secondary_minus_eci_r_primary[1] = eci_y_secondary_sc_in_span[time_step_of_tca-istep + 2] - eci_y_primary_sc_in_span[time_step_of_tca-istep + 2];
      eci_r_secondary_minus_eci_r_primary[2] = eci_z_secondary_sc_in_span[time_step_of_tca-istep + 2] - eci_z_primary_sc_in_span[time_step_of_tca-istep + 2];

      eci_v_secondary_minus_eci_v_primary[0] = eci_vx_secondary_sc_in_span[time_step_of_tca-istep + 2] - eci_vx_primary_sc_in_span[time_step_of_tca-istep + 2];
      eci_v_secondary_minus_eci_v_primary[1] = eci_vy_secondary_sc_in_span[time_step_of_tca-istep + 2] - eci_vy_primary_sc_in_span[time_step_of_tca-istep + 2];
      eci_v_secondary_minus_eci_v_primary[2] = eci_vz_secondary_sc_in_span[time_step_of_tca-istep + 2] - eci_vz_primary_sc_in_span[time_step_of_tca-istep + 2];
      v_dot( &eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary, eci_r_secondary_minus_eci_r_primary, eci_v_secondary_minus_eci_v_primary );
      p3 = 2 * eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary;

      // p4
      eci_r_secondary_minus_eci_r_primary[0] = eci_x_secondary_sc_in_span[time_step_of_tca-istep + 3] - eci_x_primary_sc_in_span[time_step_of_tca-istep + 3];
      eci_r_secondary_minus_eci_r_primary[1] = eci_y_secondary_sc_in_span[time_step_of_tca-istep + 3] - eci_y_primary_sc_in_span[time_step_of_tca-istep + 3];
      eci_r_secondary_minus_eci_r_primary[2] = eci_z_secondary_sc_in_span[time_step_of_tca-istep + 3] - eci_z_primary_sc_in_span[time_step_of_tca-istep + 3];

      eci_v_secondary_minus_eci_v_primary[0] = eci_vx_secondary_sc_in_span[time_step_of_tca-istep + 3] - eci_vx_primary_sc_in_span[time_step_of_tca-istep + 3];
      eci_v_secondary_minus_eci_v_primary[1] = eci_vy_secondary_sc_in_span[time_step_of_tca-istep + 3] - eci_vy_primary_sc_in_span[time_step_of_tca-istep + 3];
      eci_v_secondary_minus_eci_v_primary[2] = eci_vz_secondary_sc_in_span[time_step_of_tca-istep + 3] - eci_vz_primary_sc_in_span[time_step_of_tca-istep + 3];
      v_dot( &eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary, eci_r_secondary_minus_eci_r_primary, eci_v_secondary_minus_eci_v_primary );
      p4 = 2 * eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary;

    }
    else{

      // p4
      eci_r_secondary_minus_eci_r_primary[0] = eci_x_secondary_sc_in_span[time_step_of_tca+istep] - eci_x_primary_sc_in_span[time_step_of_tca+istep];
      eci_r_secondary_minus_eci_r_primary[1] = eci_y_secondary_sc_in_span[time_step_of_tca+istep] - eci_y_primary_sc_in_span[time_step_of_tca+istep];
      eci_r_secondary_minus_eci_r_primary[2] = eci_z_secondary_sc_in_span[time_step_of_tca+istep] - eci_z_primary_sc_in_span[time_step_of_tca+istep];

      eci_v_secondary_minus_eci_v_primary[0] = eci_vx_secondary_sc_in_span[time_step_of_tca+istep] - eci_vx_primary_sc_in_span[time_step_of_tca+istep];
      eci_v_secondary_minus_eci_v_primary[1] = eci_vy_secondary_sc_in_span[time_step_of_tca+istep] - eci_vy_primary_sc_in_span[time_step_of_tca+istep];
      eci_v_secondary_minus_eci_v_primary[2] = eci_vz_secondary_sc_in_span[time_step_of_tca+istep] - eci_vz_primary_sc_in_span[time_step_of_tca+istep];
      v_dot( &eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary, eci_r_secondary_minus_eci_r_primary, eci_v_secondary_minus_eci_v_primary );
      p4 = 2 * eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary;

      // p3
      eci_r_secondary_minus_eci_r_primary[0] = eci_x_secondary_sc_in_span[time_step_of_tca+istep - 1] - eci_x_primary_sc_in_span[time_step_of_tca+istep - 1];
      eci_r_secondary_minus_eci_r_primary[1] = eci_y_secondary_sc_in_span[time_step_of_tca+istep - 1] - eci_y_primary_sc_in_span[time_step_of_tca+istep - 1];
      eci_r_secondary_minus_eci_r_primary[2] = eci_z_secondary_sc_in_span[time_step_of_tca+istep - 1] - eci_z_primary_sc_in_span[time_step_of_tca+istep - 1];

      eci_v_secondary_minus_eci_v_primary[0] = eci_vx_secondary_sc_in_span[time_step_of_tca+istep - 1] - eci_vx_primary_sc_in_span[time_step_of_tca+istep - 1];
      eci_v_secondary_minus_eci_v_primary[1] = eci_vy_secondary_sc_in_span[time_step_of_tca+istep - 1] - eci_vy_primary_sc_in_span[time_step_of_tca+istep - 1];
      eci_v_secondary_minus_eci_v_primary[2] = eci_vz_secondary_sc_in_span[time_step_of_tca+istep - 1] - eci_vz_primary_sc_in_span[time_step_of_tca+istep - 1];
      v_dot( &eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary, eci_r_secondary_minus_eci_r_primary, eci_v_secondary_minus_eci_v_primary );
      p3 = 2 * eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary;

      // p2
      eci_r_secondary_minus_eci_r_primary[0] = eci_x_secondary_sc_in_span[time_step_of_tca+istep - 2] - eci_x_primary_sc_in_span[time_step_of_tca+istep - 2];
      eci_r_secondary_minus_eci_r_primary[1] = eci_y_secondary_sc_in_span[time_step_of_tca+istep - 2] - eci_y_primary_sc_in_span[time_step_of_tca+istep - 2];
      eci_r_secondary_minus_eci_r_primary[2] = eci_z_secondary_sc_in_span[time_step_of_tca+istep - 2] - eci_z_primary_sc_in_span[time_step_of_tca+istep - 2];

      eci_v_secondary_minus_eci_v_primary[0] = eci_vx_secondary_sc_in_span[time_step_of_tca+istep - 2] - eci_vx_primary_sc_in_span[time_step_of_tca+istep - 2];
      eci_v_secondary_minus_eci_v_primary[1] = eci_vy_secondary_sc_in_span[time_step_of_tca+istep - 2] - eci_vy_primary_sc_in_span[time_step_of_tca+istep - 2];
      eci_v_secondary_minus_eci_v_primary[2] = eci_vz_secondary_sc_in_span[time_step_of_tca+istep - 2] - eci_vz_primary_sc_in_span[time_step_of_tca+istep - 2];

      v_dot( &eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary, eci_r_secondary_minus_eci_r_primary, eci_v_secondary_minus_eci_v_primary );
      p2 = 2 * eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary;

      // p1
      eci_r_secondary_minus_eci_r_primary[0] = eci_x_secondary_sc_in_span[time_step_of_tca+istep - 3] - eci_x_primary_sc_in_span[time_step_of_tca+istep - 3];
      eci_r_secondary_minus_eci_r_primary[1] = eci_y_secondary_sc_in_span[time_step_of_tca+istep - 3] - eci_y_primary_sc_in_span[time_step_of_tca+istep - 3];
      eci_r_secondary_minus_eci_r_primary[2] = eci_z_secondary_sc_in_span[time_step_of_tca+istep - 3] - eci_z_primary_sc_in_span[time_step_of_tca+istep - 3];

      eci_v_secondary_minus_eci_v_primary[0] = eci_vx_secondary_sc_in_span[time_step_of_tca+istep - 3] - eci_vx_primary_sc_in_span[time_step_of_tca+istep - 3];
      eci_v_secondary_minus_eci_v_primary[1] = eci_vy_secondary_sc_in_span[time_step_of_tca+istep - 3] - eci_vy_primary_sc_in_span[time_step_of_tca+istep - 3];
      eci_v_secondary_minus_eci_v_primary[2] = eci_vz_secondary_sc_in_span[time_step_of_tca+istep - 3] - eci_vz_primary_sc_in_span[time_step_of_tca+istep - 3];
      v_dot( &eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary, eci_r_secondary_minus_eci_r_primary, eci_v_secondary_minus_eci_v_primary );
      p1 = 2 * eci_r_secondary_minus_eci_r_primary_dot_eci_v_secondary_minus_eci_v_primary;
    }

    // times
    tau1 = 1./3;
    tau2 = 2./3;
    tau1_cube = tau1 * tau1 * tau1;
    tau2_cube = tau2 * tau2 * tau2;
    tau1_square = tau1 * tau1;
    tau2_square = tau2 * tau2;
    
    // D
    D = tau1_cube * tau2_square + tau1_square * tau2 + tau1 * tau2_cube - tau1 * tau2_square - tau1_cube * tau2 - tau1_square * tau2_cube;

    // ANCAS's equations
    *gamma0 = p1;
    *gamma1 = ( ( tau2_cube - tau2_square ) * ( p2 - p1 ) + ( tau1_square - tau1_cube ) * ( p3 - p1 ) + ( tau1_cube * tau2_square - tau1_square * tau2_cube) * ( p4 - p1 ) ) /  D;
    *gamma2 = ( ( tau2 - tau2_cube ) * ( p2 - p1 ) + ( tau1_cube - tau1 ) * ( p3 - p1 ) + ( tau1 * tau2_cube - tau1_cube * tau2 ) * ( p4 - p1 ) ) /  D;
    *gamma3 = ( ( tau2_square - tau2 ) * ( p2 - p1 ) + ( tau1 - tau1_square ) * ( p3 - p1 ) + ( tau1_square * tau2 - tau1 * tau2_square ) * ( p4 - p1 ) ) /  D;


    // Check if there is a root, because if there's not then it means there is min distance between the 2 sc so no need to compute gsl_poly_solve_cubic in ancas function   
    int real_root_exists = 1;
    double min1;
    double max1;
    if ( *gamma0 > 0 ){
      min1 = *gamma1;
      if (*gamma1 + *gamma2 < min1){
	min1 = *gamma1 + *gamma2 ;
      }
      if (*gamma1 + *gamma2 + *gamma3 < min1){
	min1 = *gamma1 + *gamma2 + *gamma3 ;
      }

      if (min1 > -*gamma0){
	real_root_exists = 0;
      }
    }
    else{
      max1 = *gamma1;
      if (*gamma1 + *gamma2 > max1){
	max1 = *gamma1 + *gamma2 ;
      }
      if(*gamma1 + *gamma2 + *gamma3 > max1){
	max1 = *gamma1 + *gamma2 + *gamma3;
      }
      if (max1 < -*gamma0){
	real_root_exists = 0;
      } 
    }




    /*   ////////////////// !!!!!!!!!!!!! REMOVE BLOCK BELOW */
    /* // Normalize the coefficients so that the coefficient of order 3 is 1 */
    /*   double gamma0_normalized, gamma1_normalized, gamma2_normalized; */
    /* gamma0_normalized = *gamma0 / *gamma3; */
    /* gamma1_normalized = *gamma1 / *gamma3; */
    /* gamma2_normalized = *gamma2 / *gamma3; */
    /* double root1, root2, root3 ; */
    /* gsl_poly_solve_cubic( gamma2_normalized, gamma1_normalized, gamma0_normalized, &root1, &root2, &root3 ); */
    /* printf("%e %e %e\n", root1, root2, root3); */

    /*   ////////////////// !!!!!!!!!!!!! end of REMOVE BLOCK BELOW */

    if ( real_root_exists == 1 ){
      *min_exists = 1;

      // Record position, velocity, and acceration at P1 and P4 of each of the 2 sc. These will be used in ancas to calculate the min distance
      if ( *direction_distance <= 0 ){
	*time_step_start_interval =  time_step_of_tca-istep;
	
	r1_start[0] = eci_x_primary_sc_in_span[time_step_of_tca-istep];	r1_start[1] = eci_y_primary_sc_in_span[time_step_of_tca-istep];	r1_start[2] = eci_z_primary_sc_in_span[time_step_of_tca-istep];
	v1_start[0] = eci_vx_primary_sc_in_span[time_step_of_tca-istep];	v1_start[1] = eci_vy_primary_sc_in_span[time_step_of_tca-istep];	v1_start[2] = eci_vz_primary_sc_in_span[time_step_of_tca-istep];
	a1_start[0] = eci_ax_primary_sc_in_span[time_step_of_tca-istep];	a1_start[1] = eci_ay_primary_sc_in_span[time_step_of_tca-istep];	a1_start[2] = eci_az_primary_sc_in_span[time_step_of_tca-istep];

	r2_start[0] = eci_x_secondary_sc_in_span[time_step_of_tca-istep];	r2_start[1] = eci_y_secondary_sc_in_span[time_step_of_tca-istep];	r2_start[2] = eci_z_secondary_sc_in_span[time_step_of_tca-istep];
	v2_start[0] = eci_vx_secondary_sc_in_span[time_step_of_tca-istep];	v2_start[1] = eci_vy_secondary_sc_in_span[time_step_of_tca-istep];	v2_start[2] = eci_vz_secondary_sc_in_span[time_step_of_tca-istep];
	a2_start[0] = eci_ax_secondary_sc_in_span[time_step_of_tca-istep];	a2_start[1] = eci_ay_secondary_sc_in_span[time_step_of_tca-istep];	a2_start[2] = eci_az_secondary_sc_in_span[time_step_of_tca-istep];


	r1_end[0] = eci_x_primary_sc_in_span[time_step_of_tca-istep+3];	r1_end[1] = eci_y_primary_sc_in_span[time_step_of_tca-istep+3];	r1_end[2] = eci_z_primary_sc_in_span[time_step_of_tca-istep+3];
	v1_end[0] = eci_vx_primary_sc_in_span[time_step_of_tca-istep+3];	v1_end[1] = eci_vy_primary_sc_in_span[time_step_of_tca-istep+3];	v1_end[2] = eci_vz_primary_sc_in_span[time_step_of_tca-istep+3];
	a1_end[0] = eci_ax_primary_sc_in_span[time_step_of_tca-istep+3];	a1_end[1] = eci_ay_primary_sc_in_span[time_step_of_tca-istep+3];	a1_end[2] = eci_az_primary_sc_in_span[time_step_of_tca-istep+3];

	r2_end[0] = eci_x_secondary_sc_in_span[time_step_of_tca-istep+3];	r2_end[1] = eci_y_secondary_sc_in_span[time_step_of_tca-istep+3];	r2_end[2] = eci_z_secondary_sc_in_span[time_step_of_tca-istep+3];
	v2_end[0] = eci_vx_secondary_sc_in_span[time_step_of_tca-istep+3];	v2_end[1] = eci_vy_secondary_sc_in_span[time_step_of_tca-istep+3];	v2_end[2] = eci_vz_secondary_sc_in_span[time_step_of_tca-istep+3];
	a2_end[0] = eci_ax_secondary_sc_in_span[time_step_of_tca-istep+3];	a2_end[1] = eci_ay_secondary_sc_in_span[time_step_of_tca-istep+3];	a2_end[2] = eci_az_secondary_sc_in_span[time_step_of_tca-istep+3 ];
      }

      else{
	*time_step_start_interval =  time_step_of_tca+istep-3;


	
	r1_start[0] = eci_x_primary_sc_in_span[time_step_of_tca+istep-3];	r1_start[1] = eci_y_primary_sc_in_span[time_step_of_tca+istep-3];	r1_start[2] = eci_z_primary_sc_in_span[time_step_of_tca+istep-3];
	v1_start[0] = eci_vx_primary_sc_in_span[time_step_of_tca+istep-3];	v1_start[1] = eci_vy_primary_sc_in_span[time_step_of_tca+istep-3];	v1_start[2] = eci_vz_primary_sc_in_span[time_step_of_tca+istep-3];
	a1_start[0] = eci_ax_primary_sc_in_span[time_step_of_tca+istep-3];	a1_start[1] = eci_ay_primary_sc_in_span[time_step_of_tca+istep-3];	a1_start[2] = eci_az_primary_sc_in_span[time_step_of_tca+istep-3];

	r2_start[0] = eci_x_secondary_sc_in_span[time_step_of_tca+istep-3];	r2_start[1] = eci_y_secondary_sc_in_span[time_step_of_tca+istep-3];	r2_start[2] = eci_z_secondary_sc_in_span[time_step_of_tca+istep-3];
	v2_start[0] = eci_vx_secondary_sc_in_span[time_step_of_tca+istep-3];	v2_start[1] = eci_vy_secondary_sc_in_span[time_step_of_tca+istep-3];	v2_start[2] = eci_vz_secondary_sc_in_span[time_step_of_tca+istep-3];
	a2_start[0] = eci_ax_secondary_sc_in_span[time_step_of_tca+istep-3];	a2_start[1] = eci_ay_secondary_sc_in_span[time_step_of_tca+istep-3];	a2_start[2] = eci_az_secondary_sc_in_span[time_step_of_tca+istep-3];


	r1_end[0] = eci_x_primary_sc_in_span[time_step_of_tca+istep];	r1_end[1] = eci_y_primary_sc_in_span[time_step_of_tca+istep];	r1_end[2] = eci_z_primary_sc_in_span[time_step_of_tca+istep];
	v1_end[0] = eci_vx_primary_sc_in_span[time_step_of_tca+istep];	v1_end[1] = eci_vy_primary_sc_in_span[time_step_of_tca+istep];	v1_end[2] = eci_vz_primary_sc_in_span[time_step_of_tca+istep];
	a1_end[0] = eci_ax_primary_sc_in_span[time_step_of_tca+istep];	a1_end[1] = eci_ay_primary_sc_in_span[time_step_of_tca+istep];	a1_end[2] = eci_az_primary_sc_in_span[time_step_of_tca+istep];

	r2_end[0] = eci_x_secondary_sc_in_span[time_step_of_tca+istep];	r2_end[1] = eci_y_secondary_sc_in_span[time_step_of_tca+istep];	r2_end[2] = eci_z_secondary_sc_in_span[time_step_of_tca+istep];
	v2_end[0] = eci_vx_secondary_sc_in_span[time_step_of_tca+istep];	v2_end[1] = eci_vy_secondary_sc_in_span[time_step_of_tca+istep];	v2_end[2] = eci_vz_secondary_sc_in_span[time_step_of_tca+istep];
	a2_end[0] = eci_ax_secondary_sc_in_span[time_step_of_tca+istep];	a2_end[1] = eci_ay_secondary_sc_in_span[time_step_of_tca+istep];	a2_end[2] = eci_az_secondary_sc_in_span[time_step_of_tca-istep+3 ];

	/* if (iProc == 2){ */
	/* if (istep == time_step_of_tca  -1){ */

	/*   printf("XXXXXXX\nXXXXXXXX\nI found a min or max distance ||  r2_end[0] = %e\nXXXXXXX\nXXXXXXXX\n\n", r2_end[0]); */
	/* } */
	/* } */

      }

    


    }
    else{
      *min_exists = 0;
    }
     
  } // end of if a greater distance have been found, then compute minimum distance with ANCAS

  else{
    *min_exists = 0;
    // ptd( eci_x_secondary_sc_in_span[time_step_of_tca-istep], " eci_x_secondary_sc_in_span[time_step_of_tca-istep]");

    if ( *direction_distance <= 0 ){
    *min_distance_in_time_spanning_tca = dist_at_step_minus_n_dt;// if no min in distance has been found (in other words if the derivative of the distance with respect to time does not cancel), then save the last distance dist_at_step_minus_n_dt. This last distance dist_at_step_minus_n_dt corresponds to the distance at the beginning (if *direction_distance < 0) or end (if *direction_distance > 0) of the time spanning TCA, which is the minimum distance between over this entire time spanning TCA
    }
    else{
      //         ptd(dist_at_step_plus_n_dt, "dist_at_step_plus_n_dt");
      *min_distance_in_time_spanning_tca = dist_at_step_plus_n_dt;
    }

  }

  ///////////////////////// END OF INITIALIZE ANCAS ///////////
  /////////////////////////////////////////////////////////////

  return 0;
}


int ancas_existence_of_min_using_two_values_of_function_and_two_values_of_its_derivative( double et_start, double dt_interval, double r1_start[3], double v1_start[3], double a1_start[3], double r1_end[3],double  v1_end[3], double  a1_end[3],double r2_start[3], double v2_start[3], double a2_start[3], double r2_end[3],double  v2_end[3], double  a2_end[3], int *min_exists, double *gamma0, double *gamma1, double *gamma2, double *gamma3){ // source: Alfvano 1994)

  double fd_start, fd_dot_start, fd_dot_dot_start;
  double fd_dot_start_temp, fd_dot_dot_start_temp, fd_dot_dot_start_temp2;
  double fd_dot_end_temp, fd_dot_dot_end_temp, fd_dot_dot_end_temp2;
  double fd_end, fd_dot_end, fd_dot_dot_end;
  double rd_start[3], rd_end[3];
  double rd_dot_start[3], rd_dot_end[3];
  double rd_dot_dot_start[3], rd_dot_dot_end[3];

  char et_start_str[256];
  et2utc_c(et_start, "ISOC", 3, 255, et_start_str);
  char et_end_str[256];
  et2utc_c(et_start + dt_interval, "ISOC", 3, 255, et_end_str);

  // Set up variables
  // // rd
/* 	     print_test(); printf("XXXXXXXXXXXXXX\n\nXXXXXXXXX\n"); */
/* 	     v_print(r2_start, "r2_start"); */
/* 	     v_print(r1_start, "r1_start"); */
  v_sub( rd_start, r2_start, r1_start);
/* 	     print_test(); printf("XXXXXXXXXXXXXX\n\nXXXXXXXXX\n"); */
  v_sub( rd_end, r2_end, r1_end);
  // // rd_dot
  v_sub( rd_dot_start, v2_start, v1_start);
  v_sub( rd_dot_end, v2_end, v1_end);
  // // rd_dot_dot
  v_sub( rd_dot_dot_start, a2_start, a1_start);
  v_sub( rd_dot_dot_end, a2_end, a1_end);
  // // f
  v_dot(&fd_start, rd_start, rd_start);
  v_dot(&fd_end, rd_end, rd_end);
  // // fd_dot
  v_dot(&fd_dot_start_temp, rd_dot_start, rd_start);
  fd_dot_start = fd_dot_start_temp * 2;
  v_dot(&fd_dot_end_temp, rd_dot_end, rd_end);
  fd_dot_end = fd_dot_end_temp * 2;
  // // fd_dot_dot
  v_dot( &fd_dot_dot_start_temp, rd_dot_dot_start, rd_start );
  v_dot( &fd_dot_dot_start_temp2, rd_dot_start, rd_dot_start );
  fd_dot_dot_start = 2 * ( fd_dot_dot_start_temp + fd_dot_dot_start_temp2 );
  v_dot( &fd_dot_dot_end_temp, rd_dot_dot_end, rd_end );
  v_dot( &fd_dot_dot_end_temp2, rd_dot_end, rd_dot_end );
  fd_dot_dot_end = 2 * ( fd_dot_dot_end_temp + fd_dot_dot_end_temp2 );

  // Calculate coefficients of order 3 polynomial
  *gamma0 = fd_dot_start;
  *gamma1 = fd_dot_dot_start * dt_interval;
  *gamma2 = -3 * fd_dot_start - 2 * fd_dot_dot_start * dt_interval + 3 * fd_dot_end - fd_dot_dot_end * dt_interval;
  *gamma3 = 2 * fd_dot_start + fd_dot_dot_start * dt_interval - 2 * fd_dot_end + fd_dot_dot_end * dt_interval;

  // Check if there is a root, because if there's not then it means there is min distance between the 2 sc so no need to compute gsl_poly_solve_cubic in ancas function
  int real_root_exists = 1;
  double min1;
  double max1;
  if ( *gamma0 > 0 ){
    min1 = *gamma1;
    if (*gamma1 + *gamma2 < min1){
      min1 = *gamma1 + *gamma2 ;
    }
    if (*gamma1 + *gamma2 + *gamma3 < min1){
      min1 = *gamma1 + *gamma2 + *gamma3 ;
    }

    if (min1 > -*gamma0){
      real_root_exists = 0;
    }
  }
  else{
    max1 = *gamma1;
    if (*gamma1 + *gamma2 > max1){
      max1 = *gamma1 + *gamma2 ;
    }
    if(*gamma1 + *gamma2 + *gamma3 > max1){
      max1 = *gamma1 + *gamma2 + *gamma3;
    }
    if (max1 < -*gamma0){
      real_root_exists = 0;
    } 
  }

  if ( real_root_exists == 1 ){
    *min_exists = 1;
  }
  else{
    *min_exists = 0;
  }

  return 0;
}


int print_progress(double min_end_time, double et , double starttime, int iProc, int nb_gps){

  if (iProc == 0){ // print progress


    //               etprint(et, "time");
                printf("\033[A\33[2K\rSpOCK is propagating the spacecraft... %.0f%%\n", ( et - starttime ) / ( min_end_time - starttime ) *100.0);
	    //      printf("Propagating the constellation... %.0f%%\r", ( et - starttime ) / ( min_end_time - starttime ) *100.0);

	    
            //    etprint(et, "");
/*       frac1 = (et- starttime ) / ( min_end_time - starttime )*100.0; */
/*       frac2 = (et- starttime -dt) / ( min_end_time - starttime )*100.0; */

/*       if ( (int) frac1 > (int) frac2) { */
/*       	printf("Propagating the constellation... %.0f%%\r", */
/*       	       frac1); */
/*       } */

    
  } // end of print progress


  return 0;
}

int print_progress_epoch_sc_to_epoch_constellation(double min_end_time, double et , double starttime, int iProc, int nb_gps){

  printf("\033[A\33[2K\rPropagating the satellites until the epoch start time of the constellation... %.0f%%\n", ( et - starttime ) / ( min_end_time - starttime ) *100.0);
    
  return 0;
}

int print_progress_kalman(double min_end_time, double et , double starttime, int iProc, int nb_gps){

  if (iProc == 0){ 
    printf("\033[A\33[2K\rRunning the Kalman filter... %.0f%%\n", ( et - starttime ) / ( min_end_time - starttime ) *100.0);
  } 

  return 0;
}

int print_progress_collision(int eee_sec, int iProc, int nb_ensemble_min_per_proc, int nb_tca){

  double progress;
  progress = (double)(eee_sec - nb_ensemble_min_per_proc * iProc ) / nb_ensemble_min_per_proc / nb_tca;
  printf("Collision assessment by node %d... %.0f%%\r", iProc, progress * 100.0);
  return 0;
}


/* int send_r_v_a_in_tca_span( double ****save_r_i2cg_INRTL, double ****save_v_i2cg_INRTL, double ****save_a_i2cg_INRTL, int iProc, int nProcs, int ii, int eee, int nb_ensemble_min_per_proc, CONSTELLATION_T *CONSTELLATION ){ // this functions shares the position, velocity, and acceleration of all ensembles eee (associated with reference satellite ii) between all nodes for the entire time spanning TCA */

/*   int ccc; */
/*   int eee_all_other_iproc; */
/*   for (ccc = 0; ccc < nProcs; ccc++){ // send the position computed by a iProc to all other iProc */
/*     if ( ccc != iProc ){ */
/*       MPI_Send(&save_r_i2cg_INRTL[ii][eee][CONSTELLATION->spacecraft[ii][eee].ispan][0], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD); */
/*       MPI_Send(&save_r_i2cg_INRTL[ii][eee][CONSTELLATION->spacecraft[ii][eee].ispan][1], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD); */
/*       MPI_Send(&save_r_i2cg_INRTL[ii][eee][CONSTELLATION->spacecraft[ii][eee].ispan][2], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD); */

/*       MPI_Send(&save_v_i2cg_INRTL[ii][eee][CONSTELLATION->spacecraft[ii][eee].ispan][0], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD); */
/*       MPI_Send(&save_v_i2cg_INRTL[ii][eee][CONSTELLATION->spacecraft[ii][eee].ispan][1], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD); */
/*       MPI_Send(&save_v_i2cg_INRTL[ii][eee][CONSTELLATION->spacecraft[ii][eee].ispan][2], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD); */

/*       MPI_Send(&save_a_i2cg_INRTL[ii][eee][CONSTELLATION->spacecraft[ii][eee].ispan][0], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD); */
/*       MPI_Send(&save_a_i2cg_INRTL[ii][eee][CONSTELLATION->spacecraft[ii][eee].ispan][1], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD); */
/*       MPI_Send(&save_a_i2cg_INRTL[ii][eee][CONSTELLATION->spacecraft[ii][eee].ispan][2], 1, MPI_DOUBLE, ccc, iProc, MPI_COMM_WORLD); */
/*       //      printf("Iproc %d has value %f\n", iProc, et_initial); */
/*     } */
/*   } // end of send the position of a primary computed by a iProc to all other iProc */



/*   return 0; */
/* } */



int receive_r_v_a_in_tca_span( double ****save_r_i2cg_INRTL, double ****save_v_i2cg_INRTL, double ****save_a_i2cg_INRTL, int iProc, int nProcs, int ii, int eee, int nb_ensemble_min_per_proc, CONSTELLATION_T *CONSTELLATION, int nProcs_that_are_gonna_run_ensembles ){ // this functions shares the position, velocity, and acceleration of all ensembles eee (associated with reference satellite ii) between all nodes for the entire time spanning TCA

  int ccc;
  int eee_all_other_iproc;


  MPI_Barrier(MPI_COMM_WORLD);
  for (ccc = 0; ccc < nProcs; ccc++){ 	  // receive the positions by iProc sent by all other iProc
    if ( ccc != iProc ){
      eee_all_other_iproc = ccc * nb_ensemble_min_per_proc +  eee - iProc * nb_ensemble_min_per_proc;
      if (nProcs > 1){
	if (iProc < nProcs_that_are_gonna_run_ensembles){
      MPI_Recv(&save_r_i2cg_INRTL[ii][eee_all_other_iproc][CONSTELLATION->spacecraft[ii][eee_all_other_iproc].ispan][0], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&save_r_i2cg_INRTL[ii][eee_all_other_iproc][CONSTELLATION->spacecraft[ii][eee_all_other_iproc].ispan][1], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&save_r_i2cg_INRTL[ii][eee_all_other_iproc][CONSTELLATION->spacecraft[ii][eee_all_other_iproc].ispan][2], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      MPI_Recv(&save_v_i2cg_INRTL[ii][eee_all_other_iproc][CONSTELLATION->spacecraft[ii][eee_all_other_iproc].ispan][0], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&save_v_i2cg_INRTL[ii][eee_all_other_iproc][CONSTELLATION->spacecraft[ii][eee_all_other_iproc].ispan][1], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&save_v_i2cg_INRTL[ii][eee_all_other_iproc][CONSTELLATION->spacecraft[ii][eee_all_other_iproc].ispan][2], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      MPI_Recv(&save_a_i2cg_INRTL[ii][eee_all_other_iproc][CONSTELLATION->spacecraft[ii][eee_all_other_iproc].ispan][0], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&save_a_i2cg_INRTL[ii][eee_all_other_iproc][CONSTELLATION->spacecraft[ii][eee_all_other_iproc].ispan][1], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&save_a_i2cg_INRTL[ii][eee_all_other_iproc][CONSTELLATION->spacecraft[ii][eee_all_other_iproc].ispan][2], 1, MPI_DOUBLE, ccc, ccc, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
      }
      /* 	    if (iProc == 2){ */
      /* etprint(CONSTELLATION->et, ""); */
      /* printf("iProc %d received values from iProc %d (index %d for eee = %d) of x = %f \n", iProc, ccc, eee_all_other_iproc, eee, CONSTELLATION->spacecraft[ii][eee_all_other_iproc].r_i2cg_INRTL[0]); */
      /* 	    } */
    }
  } // end of receive the positions by iProc sent by all other iProc

  return 0;
}


double distance_between_two_sc( double eci_r_sc1[3], double eci_r_sc2[3] ){

  double diffce[3];
  v_sub(diffce, eci_r_sc1, eci_r_sc2);
  double distance;
  v_mag( &distance, diffce );

  return distance;
}

int allocate_memory_r_a_v_in_span(double ****save_x_i2cg_INRTL, double ****save_y_i2cg_INRTL, double ****save_z_i2cg_INRTL, double ****save_vx_i2cg_INRTL, double ****save_vy_i2cg_INRTL, double ****save_vz_i2cg_INRTL, double ****save_ax_i2cg_INRTL, double ****save_ay_i2cg_INRTL, double ****save_az_i2cg_INRTL, int nb_tca, int nb_satellites_not_including_gps, int total_ensemble_final_with_ref, int nb_time_steps_in_tca_time_span, int iProc, int iDebugLevel){
  //  int iiitca;
  int ii;
  int eee;
  if ( iDebugLevel >= 3 ){
    printf("---- (generate_ephemerides)(allocate_memory_r_a_v_in_span) iProc %d just got in allocate_memory_r_a_v_in_span.\n", iProc);
  }


    // allocate memory for (*save_x_i2cg_INRTL): x eci for all sc at each time step of the time spanning TCA
    // // (*save_x_i2cg_INRTL)[OPTIONS->nb_satellites_not_including_gps][total_ensemble_final_with_ref][nb_time_steps_in_tca_time_span]
    (*save_x_i2cg_INRTL) = malloc( nb_tca * sizeof( double ** ) );
    if ( (*save_x_i2cg_INRTL) == NULL ){
      print_error_any_iproc(iProc, "Not enough memory for (*save_x_i2cg_INRTL)");
    }
      for (ii = 0; ii < nb_satellites_not_including_gps; ii++){ // all ref sc
	(*save_x_i2cg_INRTL)[ii] = malloc( total_ensemble_final_with_ref * sizeof( double * ) );
	if ( (*save_x_i2cg_INRTL)[ii] == NULL ){
	  print_error_any_iproc(iProc, "Not enough memory for (*save_x_i2cg_INRTL)");
	}
	for (eee = 0; eee < total_ensemble_final_with_ref; eee++){ // all ensembles
	  (*save_x_i2cg_INRTL)[ii][eee] = malloc( nb_time_steps_in_tca_time_span * sizeof( double ) ) ;
	  if ( (*save_x_i2cg_INRTL)[ii][eee] == NULL ){
	    print_error_any_iproc(iProc, "Not enough memory for (*save_x_i2cg_INRTL)");
	  }
	} // all ensembles
      }// all ref sc


    // allocate memory for (*save_y_i2cg_INRTL): x eci for all sc at each time step of the time spanning TCA
    // // (*save_y_i2cg_INRTL)[OPTIONS->nb_satellites_not_including_gps][total_ensemble_final_with_ref][nb_time_steps_in_tca_time_span]
    (*save_y_i2cg_INRTL) = malloc( nb_tca * sizeof( double ** ) );
    if ( (*save_y_i2cg_INRTL) == NULL ){
      print_error_any_iproc(iProc, "Not enough memory for (*save_y_i2cg_INRTL)");
    }
      for (ii = 0; ii < nb_satellites_not_including_gps; ii++){ // all ref sc
	(*save_y_i2cg_INRTL)[ii] = malloc( total_ensemble_final_with_ref * sizeof( double * ) );
	if ( (*save_y_i2cg_INRTL)[ii] == NULL ){
	  print_error_any_iproc(iProc, "Not enough memory for (*save_y_i2cg_INRTL)");
	}
	for (eee = 0; eee < total_ensemble_final_with_ref; eee++){ // all ensembles
	  (*save_y_i2cg_INRTL)[ii][eee] = malloc( nb_time_steps_in_tca_time_span * sizeof( double ) ) ;
	  if ( (*save_y_i2cg_INRTL)[ii][eee] == NULL ){
	    print_error_any_iproc(iProc, "Not enough memory for (*save_y_i2cg_INRTL)");
	  }
	} // all ensembles
      }// all ref sc

    // allocate memory for (*save_z_i2cg_INRTL): x eci for all sc at each time step of the time spanning TCA
    // // (*save_z_i2cg_INRTL)[OPTIONS->nb_satellites_not_including_gps][total_ensemble_final_with_ref][nb_time_steps_in_tca_time_span]
    (*save_z_i2cg_INRTL) = malloc( nb_tca * sizeof( double ** ) );
    if ( (*save_z_i2cg_INRTL) == NULL ){
      print_error_any_iproc(iProc, "Not enough memory for (*save_z_i2cg_INRTL)");
    }
      for (ii = 0; ii < nb_satellites_not_including_gps; ii++){ // all ref sc
	(*save_z_i2cg_INRTL)[ii] = malloc( total_ensemble_final_with_ref * sizeof( double * ) );
	if ( (*save_z_i2cg_INRTL)[ii] == NULL ){
	  print_error_any_iproc(iProc, "Not enough memory for (*save_z_i2cg_INRTL)");
	}
	for (eee = 0; eee < total_ensemble_final_with_ref; eee++){ // all ensembles
	  (*save_z_i2cg_INRTL)[ii][eee] = malloc( nb_time_steps_in_tca_time_span * sizeof( double ) ) ;
	  if ( (*save_z_i2cg_INRTL)[ii][eee] == NULL ){
	    print_error_any_iproc(iProc, "Not enough memory for (*save_z_i2cg_INRTL)");
	  }
	} // all ensembles
      }// all ref sc


    // allocate memory for (*save_vx_i2cg_INRTL): x eci for all sc at each time step of the time spanning TCA
    // // (*save_vx_i2cg_INRTL)[OPTIONS->nb_satellites_not_including_gps][total_ensemble_final_with_ref][nb_time_steps_in_tca_time_span]
    (*save_vx_i2cg_INRTL) = malloc( nb_tca * sizeof( double ** ) );
    if ( (*save_vx_i2cg_INRTL) == NULL ){
      print_error_any_iproc(iProc, "Not enough memory for (*save_vx_i2cg_INRTL)");
    }
      for (ii = 0; ii < nb_satellites_not_including_gps; ii++){ // all ref sc
	(*save_vx_i2cg_INRTL)[ii] = malloc( total_ensemble_final_with_ref * sizeof( double * ) );
	if ( (*save_vx_i2cg_INRTL)[ii] == NULL ){
	  print_error_any_iproc(iProc, "Not enough memory for (*save_vx_i2cg_INRTL)");
	}
	for (eee = 0; eee < total_ensemble_final_with_ref; eee++){ // all ensembles
	  (*save_vx_i2cg_INRTL)[ii][eee] = malloc( nb_time_steps_in_tca_time_span * sizeof( double ) ) ;
	  if ( (*save_vx_i2cg_INRTL)[ii][eee] == NULL ){
	    print_error_any_iproc(iProc, "Not enough memory for (*save_vx_i2cg_INRTL)");
	  }
	} // all ensembles
      }// all ref sc


    // allocate memory for (*save_vy_i2cg_INRTL): x eci for all sc at each time step of the time spanning TCA
    // // (*save_vy_i2cg_INRTL)[OPTIONS->nb_satellites_not_including_gps][total_ensemble_final_with_ref][nb_time_steps_in_tca_time_span]
    (*save_vy_i2cg_INRTL) = malloc( nb_tca * sizeof( double ** ) );
    if ( (*save_vy_i2cg_INRTL) == NULL ){
      print_error_any_iproc(iProc, "Not enough memory for (*save_vy_i2cg_INRTL)");
    }
      for (ii = 0; ii < nb_satellites_not_including_gps; ii++){ // all ref sc
	(*save_vy_i2cg_INRTL)[ii] = malloc( total_ensemble_final_with_ref * sizeof( double * ) );
	if ( (*save_vy_i2cg_INRTL)[ii] == NULL ){
	  print_error_any_iproc(iProc, "Not enough memory for (*save_vy_i2cg_INRTL)");
	}
	for (eee = 0; eee < total_ensemble_final_with_ref; eee++){ // all ensembles
	  (*save_vy_i2cg_INRTL)[ii][eee] = malloc( nb_time_steps_in_tca_time_span * sizeof( double ) ) ;
	  if ( (*save_vy_i2cg_INRTL)[ii][eee] == NULL ){
	    print_error_any_iproc(iProc, "Not enough memory for (*save_vy_i2cg_INRTL)");
	  }
	} // all ensembles
      }// all ref sc

    // allocate memory for (*save_vz_i2cg_INRTL): x eci for all sc at each time step of the time spanning TCA
    // // (*save_vz_i2cg_INRTL)[OPTIONS->nb_satellites_not_including_gps][total_ensemble_final_with_ref][nb_time_steps_in_tca_time_span]
    (*save_vz_i2cg_INRTL) = malloc( nb_tca * sizeof( double ** ) );
    if ( (*save_vz_i2cg_INRTL) == NULL ){
      print_error_any_iproc(iProc, "Not enough memory for (*save_vz_i2cg_INRTL)");
    }
      for (ii = 0; ii < nb_satellites_not_including_gps; ii++){ // all ref sc
	(*save_vz_i2cg_INRTL)[ii] = malloc( total_ensemble_final_with_ref * sizeof( double * ) );
	if ( (*save_vz_i2cg_INRTL)[ii] == NULL ){
	  print_error_any_iproc(iProc, "Not enough memory for (*save_vz_i2cg_INRTL)");
	}
	for (eee = 0; eee < total_ensemble_final_with_ref; eee++){ // all ensembles
	  (*save_vz_i2cg_INRTL)[ii][eee] = malloc( nb_time_steps_in_tca_time_span * sizeof( double ) ) ;
	  if ( (*save_vz_i2cg_INRTL)[ii][eee] == NULL ){
	    print_error_any_iproc(iProc, "Not enough memory for (*save_vz_i2cg_INRTL)");
	  }
	} // all ensembles
      }// all ref sc

    // allocate memory for (*save_ax_i2cg_INRTL): x eci for all sc at each time step of the time spanning TCA
    // // (*save_ax_i2cg_INRTL)[OPTIONS->nb_satellites_not_including_gps][total_ensemble_final_with_ref][nb_time_steps_in_tca_time_span]
    (*save_ax_i2cg_INRTL) = malloc( nb_tca * sizeof( double ** ) );
    if ( (*save_ax_i2cg_INRTL) == NULL ){
      print_error_any_iproc(iProc, "Not enough memory for (*save_ax_i2cg_INRTL)");
    }
      for (ii = 0; ii < nb_satellites_not_including_gps; ii++){ // all ref sc
	(*save_ax_i2cg_INRTL)[ii] = malloc( total_ensemble_final_with_ref * sizeof( double * ) );
	if ( (*save_ax_i2cg_INRTL)[ii] == NULL ){
	  print_error_any_iproc(iProc, "Not enough memory for (*save_ax_i2cg_INRTL)");
	}
	for (eee = 0; eee < total_ensemble_final_with_ref; eee++){ // all ensembles
	  (*save_ax_i2cg_INRTL)[ii][eee] = malloc( nb_time_steps_in_tca_time_span * sizeof( double ) ) ;
	  if ( (*save_ax_i2cg_INRTL)[ii][eee] == NULL ){
	    print_error_any_iproc(iProc, "Not enough memory for (*save_ax_i2cg_INRTL)");
	  }
	} // all ensembles
      }// all ref sc


    // allocate memory for (*save_ay_i2cg_INRTL): x eci for all sc at each time step of the time spanning TCA
    // // (*save_ay_i2cg_INRTL)[OPTIONS->nb_satellites_not_including_gps][total_ensemble_final_with_ref][nb_time_steps_in_tca_time_span]
    (*save_ay_i2cg_INRTL) = malloc( nb_tca * sizeof( double ** ) );
    if ( (*save_ay_i2cg_INRTL) == NULL ){
      print_error_any_iproc(iProc, "Not enough memory for (*save_ay_i2cg_INRTL)");
    }
      for (ii = 0; ii < nb_satellites_not_including_gps; ii++){ // all ref sc
	(*save_ay_i2cg_INRTL)[ii] = malloc( total_ensemble_final_with_ref * sizeof( double * ) );
	if ( (*save_ay_i2cg_INRTL)[ii] == NULL ){
	  print_error_any_iproc(iProc, "Not enough memory for (*save_ay_i2cg_INRTL)");
	}
	for (eee = 0; eee < total_ensemble_final_with_ref; eee++){ // all ensembles
	  (*save_ay_i2cg_INRTL)[ii][eee] = malloc( nb_time_steps_in_tca_time_span * sizeof( double ) ) ;
	  if ( (*save_ay_i2cg_INRTL)[ii][eee] == NULL ){
	    print_error_any_iproc(iProc, "Not enough memory for (*save_ay_i2cg_INRTL)");
	  }
	} // all ensembles
      }// all ref sc

    // allocate memory for (*save_az_i2cg_INRTL): x eci for all sc at each time step of the time spanning TCA
    // // (*save_az_i2cg_INRTL)[OPTIONS->nb_satellites_not_including_gps][total_ensemble_final_with_ref][nb_time_steps_in_tca_time_span]
    (*save_az_i2cg_INRTL) = malloc( nb_tca * sizeof( double ** ) );
    if ( (*save_az_i2cg_INRTL) == NULL ){
      print_error_any_iproc(iProc, "Not enough memory for (*save_az_i2cg_INRTL)");
    }
      for (ii = 0; ii < nb_satellites_not_including_gps; ii++){ // all ref sc
	(*save_az_i2cg_INRTL)[ii] = malloc( total_ensemble_final_with_ref * sizeof( double * ) );
	if ( (*save_az_i2cg_INRTL)[ii] == NULL ){
	  print_error_any_iproc(iProc, "Not enough memory for (*save_az_i2cg_INRTL)");
	}
	for (eee = 0; eee < total_ensemble_final_with_ref; eee++){ // all ensembles
	  (*save_az_i2cg_INRTL)[ii][eee] = malloc( nb_time_steps_in_tca_time_span * sizeof( double ) ) ;
	  if ( (*save_az_i2cg_INRTL)[ii][eee] == NULL ){
	    print_error_any_iproc(iProc, "Not enough memory for (*save_az_i2cg_INRTL)");
	  }
	} // all ensembles
      }// all ref sc




  if ( iDebugLevel >= 3 ){
    printf("---- (generate_ephemerides)(allocate_memory_r_a_v_in_span) iProc %d just got out from allocate_memory_r_a_v_in_span.\n", iProc);
  }


  return 0;
}




// for a given TCA, compute_collision_between_one_secondary_and_all_primary computes the collision between a secondary sc and all primary sc during the time spanning the TCA
int compute_collision_between_one_secondary_and_all_primary(double *save_x_i2cg_INRTL_sec, double *save_y_i2cg_INRTL_sec, double *save_z_i2cg_INRTL_sec, double *save_vx_i2cg_INRTL_sec, double *save_vy_i2cg_INRTL_sec, double *save_vz_i2cg_INRTL_sec,double *save_ax_i2cg_INRTL_sec, double *save_ay_i2cg_INRTL_sec, double *save_az_i2cg_INRTL_sec, double *save_x_i2cg_INRTL_prim, double *save_y_i2cg_INRTL_prim, double *save_z_i2cg_INRTL_prim, double *save_vx_i2cg_INRTL_prim, double *save_vy_i2cg_INRTL_prim, double *save_vz_i2cg_INRTL_prim,double *save_ax_i2cg_INRTL_prim, double *save_ay_i2cg_INRTL_prim, double *save_az_i2cg_INRTL_prim,OPTIONS_T *OPTIONS, int iProc,  int *nb_coll_per_step_per_iproc_in_tca, double *et_time_step_of_save_tca, int nb_time_steps_in_tca_time_span, int iiitca, int eee_prim, int eee_sec, FILE *tca_file, FILE *dca_file, FILE *sample_file, int write_collision_files, int *eee_prim_that_collide){ 


  /* Declarations */
  eee_prim = eee_prim;  eee_sec = eee_sec; // don't worry about it (to avoid warnings at compilations (we don't want to get rid of these variables because they can be useful when debugging))
  //  char time_tca_ensemble[256];
  double direction_distance = 0;
  double min_distance_in_time_spanning_tca = 1e6;
  double tca_ensemble_1, dca_ensemble_1, tca_ensemble_2, dca_ensemble_2, tca_ensemble_3, dca_ensemble_3;
  int time_step_start_interval;
  double r_primary_start[3], v_primary_start[3], a_primary_start[3], r_primary_end[3], v_primary_end[3],  a_primary_end[3],r_secondary_start[3], v_secondary_start[3], a_secondary_start[3], r_secondary_end[3], v_secondary_end[3], a_secondary_end[3];
  int close_approach_exists;
  double gamma0,  gamma1, gamma2, gamma3;
  double et_start_of_span; 
     /* double et_end_of_span; */
     /* et_end_of_span = et_time_step_of_save_tca[iiitca] + ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt ); */
  et_start_of_span = et_time_step_of_save_tca[iiitca] - ((int)(nb_time_steps_in_tca_time_span / 2) * OPTIONS->dt );
  double et_initial_epoch, et_final_epoch;
  str2et_c(OPTIONS->initial_epoch, &et_initial_epoch);
  str2et_c(OPTIONS->final_epoch, &et_final_epoch);
  int initial_epoch_time_step_in_span = (int)(nearbyint((et_initial_epoch - et_start_of_span) / OPTIONS->dt));// technically, the cast here is not necessary. Indeed, by definition et_start_of_span is at a step 
  int final_epoch_time_step_in_span = (int)(nearbyint((et_final_epoch - et_start_of_span) / OPTIONS->dt)); // technically, the cast here is not necessary. Indeed, by definition et_start_of_span is at a step and so is et_final_epoch (even if the user chose a dt such that the final epoch does not fall in a multiple of dt + inital epoch (-> see in load_option))
  /* etprint(et_start_of_span + final_epoch_time_step_in_span * OPTIONS->dt, "final_epoch_time_step_in_span"); */
  /* etprint(et_start_of_span +  nb_time_steps_in_tca_time_span * OPTIONS->dt , "et_end_of_span" ); */
  /* etprint(et_start_of_span + initial_epoch_time_step_in_span * OPTIONS->dt, "initial_epoch_time_step_in_span")  ; */
  /* etprint(et_start_of_span,"et_start_of_span"); */
  //  pti(final_epoch_time_step_in_span,"final_epoch_time_step_in_span");
  int min_time_step_final_epoch_and_time_step_end_span = nb_time_steps_in_tca_time_span-1;
  if ( final_epoch_time_step_in_span < min_time_step_final_epoch_and_time_step_end_span ){
    min_time_step_final_epoch_and_time_step_end_span = final_epoch_time_step_in_span;
  }
  int max_time_step_initial_epoch_and_time_step_start_span = 0;
  if ( initial_epoch_time_step_in_span > max_time_step_initial_epoch_and_time_step_start_span ){
    max_time_step_initial_epoch_and_time_step_start_span = initial_epoch_time_step_in_span;
  }
  //    etprint(max_time_step_initial_epoch_and_time_step_start_span * OPTIONS->dt+ et_start_of_span, "max_time_step_initial_epoch_and_time_step_start_span");
  int time_step_of_tca =  (int)(nb_time_steps_in_tca_time_span / 2);
  int step_coll;
  /* Algorithm */
  //printf("secondary (%e %e %e) | primary (%e %e %e) || (iProc: %d)\n", save_x_i2cg_INRTL_sec[time_step_of_tca], save_y_i2cg_INRTL_sec[time_step_of_tca], save_z_i2cg_INRTL_sec[time_step_of_tca], save_x_i2cg_INRTL_prim[time_step_of_tca],  save_y_i2cg_INRTL_prim[time_step_of_tca],  save_z_i2cg_INRTL_prim[time_step_of_tca], iProc );

    tca_ensemble_1 = -1; dca_ensemble_1 = 1e6; tca_ensemble_2 = -1; dca_ensemble_2 = 1e6; tca_ensemble_3 = -1;  dca_ensemble_3 = 1e6;
    time_step_start_interval = -1;
    // close_approach_ensemble determines if there is a close approach (defined as a minimum of the distance) between the primary ensemble sc eee_prim and the secondary ensemble sc eee in the interval of time spanning TCA (unperturbed). If there is, it returns the position, velocity, and acceleration of of sc eee_prim and sc eee at both the beginning and end of interval of time 3 * OPTIONS->dt in which the close approach between sc eee_prim and sc eee is expected to be. It also returns the coefficients gamma that are then used to calculate the min distance using the function ancas (if there is a min distanc)
    //            printf("iProc: %d | eee_prim: %d | eee_sec: %d\n", iProc, eee_prim, eee_sec);
    close_approach_ensemble(  save_x_i2cg_INRTL_prim, save_y_i2cg_INRTL_prim, save_z_i2cg_INRTL_prim, save_vx_i2cg_INRTL_prim, save_vy_i2cg_INRTL_prim, save_vz_i2cg_INRTL_prim, save_ax_i2cg_INRTL_prim, save_ay_i2cg_INRTL_prim, save_az_i2cg_INRTL_prim,
			      save_x_i2cg_INRTL_sec, save_y_i2cg_INRTL_sec, save_z_i2cg_INRTL_sec, save_vx_i2cg_INRTL_sec, save_vy_i2cg_INRTL_sec, save_vz_i2cg_INRTL_sec, save_ax_i2cg_INRTL_sec, save_ay_i2cg_INRTL_sec, save_az_i2cg_INRTL_sec,
			      time_step_of_tca,   &gamma0,  &gamma1,  &gamma2,  &gamma3, &close_approach_exists, r_primary_start, v_primary_start, a_primary_start, r_primary_end, v_primary_end,  a_primary_end,r_secondary_start, v_secondary_start, a_secondary_start, r_secondary_end, v_secondary_end, a_secondary_end, &time_step_start_interval, &min_distance_in_time_spanning_tca, &direction_distance, initial_epoch_time_step_in_span, final_epoch_time_step_in_span, iProc, et_time_step_of_save_tca[iiitca], OPTIONS); // here start means the oldest time of 4 four latest points, and end means the most recent time of 4 four latest points (see Alfano 2009) // time_step_start_interval represents the time step in the time spanning TCA (unperturbed) for which the interval of four points starts (interval for which a min distance between the primary ensemble sc (represented by eee_prim) and the secondary ensemble sc (represented by eee) is found (in other words, this interal is the one represnted by P1, P2, P3, and P4 in Alfano 2009)). For instance, if time_step_start_interval = 3, then it is at the third time step of the interval of time spanning TCA (unperturbed) (recall: the interval starts at TCA - span/2) that the interval of four points (P1, P2, P3, and P4) starts. So the close approach between the primary and the secondary sc is expected to be between the third and the seventh time step of the time spanning TCA (unpertubed). In this discussion, there is 2 intervals. The first one is the intervl spanning TCA (unpertubed): varies from TCA - span/2 to TCA + span/2. The second interval we are talking about is the interval in which the min distance between the primary ensemble sc and the secondary ensemble sc is expected to be (it is represented by the four points P1, P2, P3, and P4 in Alfano 2009
    /* ptd(distance_between_two_sc( save_r_i2cg_INRTL[0][time_step_of_tca], save_r_i2cg_INRTL[1][eee][time_step_of_tca]), "D" ); */
    /* etprint(save_tca[iiitca], "TCA unpertubed"); */
    /* etprint(et_time_step_of_save_tca[iiitca], "Time step of TCA untpertubed"); */
    /* etprint(et_start_of_span , "span start"); */
    /* etprint(et_end_of_span , "span end"); */
    /* printf("%f = %d\n", nb_time_steps_in_tca_time_span / 2., time_step_of_tca); */
    /* etprint(et_start_of_span + time_step_of_tca * OPTIONS->dt, "check TCA"); */

    /* pti(close_approach_exists, "1"); */
    //		  test_print_iproc(iProc, "Q");
    /* printf("close_approach_exists: %d (iProc: %d)\n", close_approach_exists, iProc); */
    if ( close_approach_exists == 1 ){ // if a close approach in the time spanning TCA (of the unpertubed orbits) has been found betwee eee_prim and eee (which represents a secondary sc) then compute ANCAS to find TCA and DCA of this close approach
      tca_ensemble_1 = -1; dca_ensemble_1 = 1e6; tca_ensemble_2 = -1; dca_ensemble_2 = 1e6; tca_ensemble_3 = -1; dca_ensemble_3 = 1e6;
      //		    	    test_print_iproc( iProc, "X");
      //            printf("iProc: %d | eee_prim: %d | eee_sec: %d\n", iProc, eee_prim, eee_sec);
      ancas(OPTIONS->min_dist_close_approach, et_start_of_span + time_step_start_interval * OPTIONS->dt, 3 * OPTIONS->dt, r_primary_start, v_primary_start, a_primary_start, r_primary_end, v_primary_end,  a_primary_end,r_secondary_start, v_secondary_start, a_secondary_start, r_secondary_end, v_secondary_end, a_secondary_end, gamma0, gamma1, gamma2, gamma3, &tca_ensemble_1, &dca_ensemble_1, &tca_ensemble_2, &dca_ensemble_2, &tca_ensemble_3, &dca_ensemble_3 ); // 3 * OPTIONS->dt necause the interval over which the min distance is calculated of time covers 4 points (it is represented by the four points P1, P2, P3, and P4 in Alfano 2009)
      /* ptd(distance_between_two_sc( save_r_i2cg_INRTL[0][time_step_of_tca], save_r_i2cg_INRTL[1][eee][time_step_of_tca]), "dstance between primary and secondary ensemble" ); */
      /* etprint(et_start_of_span + time_step_start_interval * OPTIONS->dt, "start interval Pi"); */

      /* if ( ( tca_ensemble_1 > et_end_of_span - 2 * OPTIONS-> dt ) || ( tca_ensemble_2 > et_end_of_span - 2 * OPTIONS-> dt ) || ( tca_ensemble_3 > et_end_of_span - 2 * OPTIONS-> dt ) ){ */

      /* 	etprint(tca_ensemble_1, "tca_ensemble_1"); */
      /* 	etprint(tca_ensemble_2, "tca_ensemble_2"); */
      /* 	etprint(tca_ensemble_3, "tca_ensemble_3"); */
      /* } */ 

/*       if ( (dca_ensemble_1 <= OPTIONS->min_dist_collision * 2)){// || (dca_ensemble_2 <= OPTIONS->min_dist_collision * 10) || (dca_ensemble_3 <= OPTIONS->min_dist_collision * 10) ){ */
/* 	*eee_prim_that_collide = eee_prim; */
/* 	//	printf("collide %d\n", *eee_prim_that_collide); */
/* 	if (write_collision_files == 1){ */
/* 	  //	  printf("sample file: <%s>\n", sample_file); */
/* 	  fprintf(sample_file, "[%d %d %f]",  eee_prim, eee_sec, dca_ensemble_1); //	  	  	fprintf(sample_file, "%d ",  eee_prim); */
/* 	} */
/*       } */

      /* test_print_iproc( iProc, "Y"); */
      /* printf("%e %e %e %e (iProc = %d)\n", dca_ensemble_1, dca_ensemble_2, dca_ensemble_3, OPTIONS->min_dist_collision, iProc); */
      if (  dca_ensemble_1 <= OPTIONS->min_dist_collision )  {
    	step_coll = (int)( (tca_ensemble_1 - et_start_of_span) / OPTIONS->dt ); // time step in time spanning TCA where the collision between secondary sc ewee and primary sc eee_prim occurs
    	nb_coll_per_step_per_iproc_in_tca[step_coll] = nb_coll_per_step_per_iproc_in_tca[step_coll] + 1;
	if (write_collision_files == 1){
	fprintf(tca_file, "%f ", tca_ensemble_1 - et_initial_epoch );
	fprintf(dca_file, "%f ", dca_ensemble_1 );

	}
	/* if (iProc == 1){ */
	/*   if (step_coll == 348){ */
	/* /\* printf("HHHHHHHHHHHHHHHHHHHH\n"); *\/ */
	/* /\* et2utc_c(et_start_of_span + time_step_start_interval * OPTIONS->dt, "ISOC", 3, 255, time_tca_ensemble); *\/ */
	/* /\* printf("Primary sc #%d collide with secondary sc #%d between %s and ", eee_prim, eee_sec, time_tca_ensemble); *\/ */
	/* /\* et2utc_c(et_start_of_span + time_step_start_interval * OPTIONS->dt + 3 * OPTIONS->dt, "ISOC", 3, 255, time_tca_ensemble); *\/ */
	/* /\* printf("%s\n", time_tca_ensemble); *\/ */
	/* /\* et2utc_c(tca_ensemble_1, "ISOC", 3, 255, time_tca_ensemble); *\/ */
	/* /\* printf("TCA: %s\nDCA: %e km\n",time_tca_ensemble, dca_ensemble_1 ); *\/ */
	/* /\* printf("START\n"); *\/ */
	/* /\* v_print(r_primary_start,"r_primary_start"); *\/ */
	/* /\* v_print(v_primary_start,"v_primary_start"); *\/ */
	/* /\* v_print(a_primary_start,"a_primary_start"); *\/ */
	/* /\* v_print(r_secondary_start,"r_secondary_start"); *\/ */


	/* /\* v_print(r_primary_end,"r_primary_end"); *\/ */
	/* /\* v_print(r_secondary_end,"r_secondary_end"); *\/ */
	/* /\* et2utc_c(et_start_of_span, "ISOC", 3, 255, time_tca_ensemble); *\/ */
	/*     /\* et2utc_c(tca_ensemble_1, "ISOC", 3, 255, time_tca_ensemble); *\/ */
	/*     /\* 	printf("(d1) %s - nb_coll_per_step_per_iproc[%d]: %d <prim: %d x: %e | sec: %d> (iProc: %d)\n", time_tca_ensemble, step_coll,  nb_coll_per_step_per_iproc_in_tca[step_coll], eee_prim, save_x_i2cg_INRTL_prim[step_coll] , eee_sec, iProc); *\/ */
	/*   } */
	/* 	} */

      }
	    	
      if (  dca_ensemble_2 <= OPTIONS->min_dist_collision )  {
    	step_coll = (int)( (tca_ensemble_2 - et_start_of_span) / OPTIONS->dt ); // time step in time spanning TCA where the collision between secondary sc ewee and primary sc eee_prim occurs
	//    	printf("[iiitca = %d][iProc = %d][step_coll = %d]\n",  iiitca, iProc, step_coll);
    	nb_coll_per_step_per_iproc_in_tca[step_coll] = nb_coll_per_step_per_iproc_in_tca[step_coll] + 1;
	if (write_collision_files == 1){
	fprintf(tca_file, "%f ", tca_ensemble_2 - et_initial_epoch );
	fprintf(dca_file, "%f ", dca_ensemble_2  );
	}
	/* if (iProc == 1){ */
	/*   	if (step_coll == 348){ */
	/* et2utc_c(tca_ensemble_2, "ISOC", 3, 255, time_tca_ensemble); */
    	/* printf("(d2) %s - nb_coll_per_step_per_iproc[%d]: %d <prim: %d x: %e | sec: %d> (iProc: %d)\n", time_tca_ensemble, step_coll,  nb_coll_per_step_per_iproc_in_tca[step_coll], eee_prim, save_x_i2cg_INRTL_prim[step_coll] , eee_sec, iProc); */

	/*   } */
	/* 	} */

      }
      if (  dca_ensemble_3 <= OPTIONS->min_dist_collision )  {
    	step_coll = (int)( (tca_ensemble_3 - et_start_of_span) / OPTIONS->dt ); // time step in time spanning TCA where the collision between secondary sc ewee and primary sc eee_prim occurs
	//    	printf("[iiitca = %d][iProc = %d][step_coll = %d]\n",  iiitca, iProc, step_coll);
    	nb_coll_per_step_per_iproc_in_tca[step_coll] = nb_coll_per_step_per_iproc_in_tca[step_coll] + 1;
	if (write_collision_files == 1){
	fprintf(tca_file, "%f ", tca_ensemble_3 - et_initial_epoch );
	fprintf(dca_file, "%f ", dca_ensemble_3);
	}
	/* if (iProc == 1){ */
	/*   	if (step_coll == 348){ */
	/* et2utc_c(tca_ensemble_3, "ISOC", 3, 255, time_tca_ensemble); */
    	/* printf("(d3) %s - nb_coll_per_step_per_iproc[%d]: %d <prim: %d x: %e | sec: %d> (iProc: %d)\n", time_tca_ensemble, step_coll,  nb_coll_per_step_per_iproc_in_tca[step_coll], eee_prim, save_x_i2cg_INRTL_prim[step_coll] , eee_sec, iProc); */

	/* 	  } */
	/* } */
      }

      //		    	    test_print_iproc( iProc, "Z");

    } // end of if a close approach in the time spanning TCA (of the unpertubed orbits) has been found betwee eee_prim and eee (which represents a secondary sc) then compute ANCAS to find TCA and DCA of this close approach

    /* else if (min_distance_in_time_spanning_tca < OPTIONS->min_dist_collision ){ // if the the distance between the 2 sc has not reached a minimum during the time spanning the TCA then look at the distances at the borders of the span (the one at the beginning of the span and the one at the end of the span). If one (or two) of these distances is < min_dist_collision, then count this case as a collision (only one collision, even if both distances are < min_dist_collision) */
      
    /*   //      printf("iProc: %d | eee_prim: %d | eee_sec: %d || %e %e\n", iProc, eee_prim, eee_sec, min_distance_in_time_spanning_tca, direction_distance); */
    /*   if (direction_distance >= 0){ // this will add colllision at the end of the span */
    /*    	  ptd(min_distance_in_time_spanning_tca, "min_distance_in_time_spanning_tca"); */
    /* 	  nb_coll_per_step_per_iproc_in_tca[min_time_step_final_epoch_and_time_step_end_span] = nb_coll_per_step_per_iproc_in_tca[min_time_step_final_epoch_and_time_step_end_span] + 1; */
    /* 	} */
    /*   else{ // this will add collisions at the start of the span */
    /* 	  //	  ptd(min_distance_in_time_spanning_tca, "min_distance_in_time_spanning_tca"); */
    /* 	      	  nb_coll_per_step_per_iproc_in_tca[ max_time_step_initial_epoch_and_time_step_start_span ] = nb_coll_per_step_per_iproc_in_tca[ max_time_step_initial_epoch_and_time_step_start_span] + 1; */
    /* 	} */
    /*   } */
      


  return 0;
}


							      // computes the heading of a sc based on its ECEF velocity. 0 is North, 90 is East, 180 is South, 270 is West
							      int compute_heading(double *heading, double v_ecef[3], double lon, double lat, double earth_flattening){ 
								// lon, lat, heading in radians
  
								double local_north_in_enu[3];
								double local_east_in_enu[3];
								double T_enu_to_ecef[3][3];
								double T_ecef_to_enu[3][3];
								double v_enu[3];
  local_north_in_enu[0] = 0; local_north_in_enu[1] = 1; local_north_in_enu[2] = 0;
  local_east_in_enu[0] = 1; local_east_in_enu[1] = 0; local_east_in_enu[2] = 0;
      // // convert velocity from ECEF to ENU coordinates of the sub satellite location
      compute_T_enu_to_ecef( T_enu_to_ecef, lat, lon, earth_flattening);
      m_trans(T_ecef_to_enu, T_enu_to_ecef);
      m_x_v(v_enu, T_ecef_to_enu, v_ecef);
      // // Dot product between the local north and the vector ground station to sc in ENU reference system
      double unit_v_enu[3];
      v_norm(unit_v_enu, v_enu);
      double unit_v_enu_dot_local_north_in_enu;
      v_dot(&unit_v_enu_dot_local_north_in_enu, unit_v_enu, local_north_in_enu);

      double unit_v_enu_dot_local_east_in_enu;
      v_dot(&unit_v_enu_dot_local_east_in_enu, unit_v_enu, local_east_in_enu);
      if ( unit_v_enu_dot_local_east_in_enu >= 0 ){
	*heading = acos( unit_v_enu_dot_local_north_in_enu );
      }
      else{
	*heading = 2*M_PI - acos( unit_v_enu_dot_local_north_in_enu );
      }

  return 0;
}


