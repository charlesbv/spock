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

#include <sys/time.h>
#include <mpi.h>
#include "propagator.h"
#include "options.h"
#include <time.h>
#include <limits.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "f2c.h"
#include <gsl/gsl_statistics.h>
#include "gsl/gsl_sf_legendre.h"

/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           initialize_constellation
//  Purpose:        Initializes a constellation based on the period of the 1st spacecraft
//  Assumptions:    None.
//  References      Various
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | J. Getchius   | 05/20/2015    |   ---     | Initial Implementation
//      | C. Bussy-Virat| 07/23/2015    |   ---     | Initialization of ALL the orbital elements from the input file; Replace altitude by apogee altitude in the input file
//
/////////////////////////////////////////////////////////////////////////////////////////
int initialize_constellation(CONSTELLATION_T *CONSTELLATION, OPTIONS_T       *OPTIONS, PARAMS_T        *PARAMS, GROUND_STATION_T *GROUND_STATION, int iDebugLevel, int iProc, int nProcs )

{
 


  double **r_alongtrack_all, **r_crosstrack_all, **r_radial_all, **bc_pert_all, **srp_pert_all;
  double **v_alongtrack_all, **v_crosstrack_all, **v_radial_all;
  // Declarations	
  double geodetic[3];
  int bbb;
  int m_time_span;

  /* gsl_rng * r_gaussian_generator; */
  /* long seed; */
  /* r_gaussian_generator = gsl_rng_alloc(gsl_rng_mt19937); */

  //  int start_ensemble;
  double period;
  int ccc;
  double collision_equinoctial_sigma_in_diag_basis[6];
  double collision_equinoctial_sigma_in_equinoctial_basis[6];

  double collision_eci_sigma_in_diag_basis[6];
  double collision_eci_sigma_in_eci_basis[6];
  
  /* double collision_x_eci_sigma, collision_y_eci_sigma, collision_z_eci_sigma, collision_vx_eci_sigma, collision_vy_eci_sigma, collision_vz_eci_sigma; */
  /* double collision_x_eci_sigma_in_diag_basis, collision_y_eci_sigma_in_diag_basis, collision_z_eci_sigma_in_diag_basis, collision_vx_eci_sigma_in_diag_basis, collision_vy_eci_sigma_in_diag_basis, collision_vz_eci_sigma_in_diag_basis; */
  /* int aaa_count; */
  /* int swpc_forecast_day; */
  /* double nb_index_in_81_days; */

  //  int aaa_sigma;

  int write_attitude_file = 0;
  char time_attitude[256];
  FILE *file_attitude[N_SATS];
  char filename_attitude[N_SATS][256];
  double et_initial_epoch, et_final_epoch;
  int iground;
  SpiceDouble       xform[6][6];
  double estate[6], jstate[6];

  int hhh;
  char *next;
  int find_file_name;
 
  double altitude_perigee;
  int index_altitude_right_below_perigee_save = 1000000000;

  int iHeaderLength_gitm;
  FILE *gitm_file;
  int i,j,k;
  int ggg;
  double random_pitch_angular_velocity, random_roll_angular_velocity, random_yaw_angular_velocity;
  int time_before_reset;
  double inclination_sigma = 0;
  double w_sigma = 0;
  double long_an_sigma = 0;
  double f_sigma = 0;
  double eccentricity_sigma = 0;
  double sma_sigma = 0;
  int eee;
  double et;
  /* double T_J2000_to_ECEF[3][3]; */
  /* double T_ECEF_to_J2000[3][3]; */
  int ii,sss,nnn,aaa;
  double radius_perigee;
  FILE *gps_tle_file;

  //  Begin Calcs
  // Initialize Time
  str2et_c(OPTIONS->initial_epoch, &et); // Convert a string representing an epoch to a double precision value representing the number of TDB seconds past the J2000 epoch corresponding to the input epoch (in /Users/cbv/cspice/src/cspice/str2et_c.c) (CBV)
  OPTIONS->et_initial_epoch = et; // should have done that in load_options.c
  str2et_c(OPTIONS->initial_epoch, &et_initial_epoch);
  str2et_c(OPTIONS->final_epoch, &et_final_epoch);

  r_alongtrack_all = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof(double *) );
  r_crosstrack_all = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof(double *) );
  r_radial_all = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof(double *) );

  v_alongtrack_all = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof(double *) );
  v_crosstrack_all = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof(double *) );
  v_radial_all = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof(double *) );

  bc_pert_all = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof(double *) );
 srp_pert_all = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof(double *) );


  CONSTELLATION->et = et;
  CONSTELLATION->aaa_sigma = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof(int) );
  CONSTELLATION->aaa_mod = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof(int) );
  CONSTELLATION->sum_sigma_for_f107_average = malloc( ( OPTIONS->nb_satellites_not_including_gps ) * sizeof( double * ) );
  if (CONSTELLATION->sum_sigma_for_f107_average == NULL){
    print_error_any_iproc(iProc, "Not enough memory for sum_sigma_for_f107_average");
  }
  if (OPTIONS->nb_ensembles_density) { // the user chose to run ensembles on the density using data from SWPC 
    if (OPTIONS->swpc_need_predictions){ // if future predictions of F10.7 and Ap (not only past observations)
      CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time = malloc( nProcs * sizeof( double *));
      CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time_sorted = malloc( nProcs * sizeof( double *));
      CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time = malloc( nProcs * sizeof( double *));
      CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time_sorted = malloc( nProcs * sizeof( double *));

    }
  }

  CONSTELLATION->spacecraft = malloc( OPTIONS->n_satellites * sizeof(SPACECRAFT_T*) ); 

  if ( CONSTELLATION->spacecraft == NULL){
    printf("***! Could not allow memory to CONSTELLATION->spacecraft. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }




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
  }
    
  //    if (iProc == 0){
  /*     for (iProcf = 0; iProcf < nProcs; iProcf++){ */
  /*       printf("%d - %d\n", iStart_save[iProcf], iEnd_save[iProcf] - 1) ; */
  /*     } */
  //}

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
    //    printf("iProc %d | start_ensemble[%d] %d\n", iProc, ii, start_ensemble[ii]);
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
    //   printf("iProc %d: array_sc[%d] = %d\n", iProc, ielt, array_sc[ielt]);
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



  /*************************************************************************************/
  /*************************************************************************************/
  /*************************************************************************************/
  /*************************************************************************************/
  /*************************************************************************************/
  /**************************** SPACECRAFTS OTHER THAN GPS *****************************/
  /*************************************************************************************/
  /*************************************************************************************/
  /*************************************************************************************/
  /*************************************************************************************/
  /************************************************************************************/

  if (iDebugLevel >= 1){
    if (iProc == 0) printf("-- Number of spacecraft: %d\n", OPTIONS->n_satellites - OPTIONS->nb_gps);
  }

  for (ii = 0; ii < OPTIONS->n_satellites - OPTIONS->nb_gps; ii++){ // go through all main SC other than GPS

    CONSTELLATION->aaa_sigma[ii] =  OPTIONS->aaa_sigma;
        CONSTELLATION->aaa_mod[ii] =  OPTIONS->aaa_mod[ii];
   
    //  printf("ooooooooooooooooo %d %d\n",CONSTELLATION->aaa_mod[0], CONSTELLATION->aaa_mod[1]);
    //    printf(" CONSTELLATION->aaa_mod[%d] %d\n", ii,  CONSTELLATION->aaa_mod[ii]);
    CONSTELLATION->spacecraft[ii] = malloc( ( OPTIONS->nb_ensembles_min + 1 ) * sizeof(SPACECRAFT_T) ); 
r_alongtrack_all[ii] = malloc( ( OPTIONS->nb_ensembles_min ) * sizeof(double) ); 
r_crosstrack_all[ii] = malloc( ( OPTIONS->nb_ensembles_min) * sizeof(double) ); 
r_radial_all[ii] = malloc( ( OPTIONS->nb_ensembles_min) * sizeof(double) ); 

v_alongtrack_all[ii] = malloc( ( OPTIONS->nb_ensembles_min ) * sizeof(double) ); 
v_crosstrack_all[ii] = malloc( ( OPTIONS->nb_ensembles_min) * sizeof(double) ); 
v_radial_all[ii] = malloc( ( OPTIONS->nb_ensembles_min) * sizeof(double) ); 

 bc_pert_all[ii] = malloc( ( OPTIONS->nb_ensembles_min) * sizeof(double) ); 
 srp_pert_all[ii] = malloc( ( OPTIONS->nb_ensembles_min) * sizeof(double) ); 


    if ( CONSTELLATION->spacecraft[ii] == NULL){
      printf("***! Could not allow memory to CONSTELLATION->spacecraft[ii]. The program will stop. !***\n"); MPI_Finalize(); exit(0);
    }

    CONSTELLATION->sum_sigma_for_f107_average[ii] = malloc( ( OPTIONS->nb_ensemble_min_per_proc * nProcs + 1 ) * sizeof( double ) );
    if (CONSTELLATION->sum_sigma_for_f107_average[ii] == NULL){
      print_error_any_iproc(iProc, "Not enough memory for sum_sigma_for_f107_average[ii]");
    }


    // Initialize the ensemble parameters
    srand (time (NULL));
    // Initialize the ensemble parameters on COE
    if (( OPTIONS->nb_ensembles > 0 ) && ( strcmp(OPTIONS->type_orbit_initialisation, "oe" ) == 0 )){ // if we run ensembles on the orbital elements
      sma_sigma = OPTIONS->apogee_alt_sigma[ii] / ( 1 + OPTIONS->eccentricity[ii] );
      inclination_sigma    = OPTIONS->inclination_sigma[ii] * DEG2RAD;
      w_sigma               = OPTIONS->w_sigma[ii] * DEG2RAD;                                // argument of perigee (CBV)
      long_an_sigma        = OPTIONS->long_an_sigma[ii] * DEG2RAD;                          // RAAN (CBV)
      f_sigma               = OPTIONS->f_sigma[ii] * DEG2RAD;                                // true anomaly (CBV)
      eccentricity_sigma    = OPTIONS->eccentricity_sigma[ii];
    }

    /*     if (iProc == 0){ */
    /*       start_ensemble = 0; */
    /*   } */
    /*     else{ */
    /*       start_ensemble = 1; */
    /*     } */

    //      for (eee = start_ensemble + iProc * OPTIONS->nb_ensemble_min_per_proc; eee< 1 + iProc * OPTIONS->nb_ensemble_min_per_proc + OPTIONS->nb_ensemble_min_per_proc ; eee++){ // go through all sc (including ensembels, if you run ensembles) other than GPS
    for ( ielt = 0; ielt < OPTIONS->nb_ensemble_min_per_proc + 1; ielt ++ ){ // go through all sc (including ensembels, if you run ensembles) other than GPS. Each iProc runs OPTIONS->nb_ensemble_min_per_proc ensemble sc. But not all iProc run a main sc. If the iProc runs or not a main sc (and which main sc) is decided in the variable start_ensemble[ii]: if start_ensemble[ii] = 0 for this iProc, then this iProc runs main sc ii. If start_ensemble[ii] = 1 for this iProc, then this iProc does NOT run main sc ii. 
      eee = array_sc[ielt]; // eee represents the sc that is run: eee = 0 represents a main sc. eee > 0 represents an ensemble sc
      if ( (  start_ensemble[ii] == 0 ) || ( ( start_ensemble[ii] != 0 ) && ( eee > 0 ) ) ){ // if main sc ii is run by this iProc OR  ( if this main sc is not run by this iProc and eee corresponds to an ensemble sc (so this iProc does not run main sc ii but runs an ensemble corresponding to main sc ii))
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.isGPS = 0;
	/*             if (iProc == 0){ */
	/*       printf("iProc %d runs main sc %d and ensemble sc %d\n", iProc, ii, eee); */
	/*             } */
	/* if (iProc == 1){ */
	/* printf("iProc: %d | eee = %d | from %d to %d | ii = %d\n", iProc, eee, start_ensemble + iProc * OPTIONS->nb_ensemble_min_per_proc, iProc * OPTIONS->nb_ensemble_min_per_proc + OPTIONS->nb_ensemble_min_per_proc, ii); */
	/* } */
	/*************************************************************************************/
	/*************************************************************************************/
	/*************************************************************************************/
	/*************************************************************************************/
	/********************* INITIALIZE: et, r_i2cg_INRTL, v_i2cg_INRTL, OE ****************/
	/*************************************************************************************/
	/*************************************************************************************/
	/*************************************************************************************/
	/*************************************************************************************/

	/*************************************************************************************/
	/*************************************************************************************/
	/*************************************************************************************/
	/************** INITIALIZE: et, r_i2cg_INRTL, v_i2cg_INRTL, OE FOR MAIN SC ***********/
	/*************************************************************************************/
	/*************************************************************************************/
	/*************************************************************************************/

	//      if (eee == 0){ // eee = 0 represents the main spacecraft (eee > 0 is a sc from an ensemble)
	if ( (start_ensemble[ii] == 0) && (eee == 0)){ // if this iProc runs main sc ii and that eee corresponds to the main sc. eee = 0 represents the main spacecraft (eee > 0 is a sc from an ensemble). start_ensemble[ii] = 0 if main sc is run by this iProc

	  /*************************************************************************************/
	  /*************************************************************************************/
	  /********** INITIALIZE: et, r_i2cg_INRTL, v_i2cg_INRTL, OE FOR MAIN SC BY TLE ************/
	  /*************************************************************************************/
	  /*************************************************************************************/

	  if ( (strcmp( OPTIONS->type_orbit_initialisation, "tle" ) == 0 ) || (strcmp(OPTIONS->type_orbit_initialisation, "tle_sgp4" ) == 0 )){
	    if (iDebugLevel >= 1){
	      if (iProc == 0) printf("-- (initialize_constellation) Initializing et, r_i2cg_INRTL, v_i2cg_INRTL, orbital elements from TLE.\n");
	    }

	    /*** Convert the TLEs into inertial state (postion, velocity) ***/
	    // Initialize TLE parameters
	    SpiceInt frstyr = 1961; // do not really care about this line (as long as the spacecrafts are not flying after 2060)
	    //	    SpiceDouble geophs[8];
	    int pp;
	    int lineln=0;
	    size_t len = 0;
	    char *line_temp = NULL;
	    //	    SpiceDouble elems[10];
	    SpiceDouble state[6];
	    SpiceDouble epoch_sat; // epoch of TLE of the sat
	    char sat_tle_file_temp[256];
	    FILE *sat_tle_file;
	    ORBITAL_ELEMENTS_T OE_temp;

	    /* Set up the geophysical quantities.  At last check these were the values used by Space Command and SGP4 */
	    /* PARAMS->geophs[ 0 ] =    1.082616e-3;   // J2 */
	    /* PARAMS->geophs[ 1 ] =   -2.53881e-6;    // J3 */
	    /* PARAMS->geophs[ 2 ] =   -1.65597e-6;    // J4 */
	    /* PARAMS->geophs[ 3 ] =    7.43669161e-2; // KE */
	    /* PARAMS->geophs[ 4 ] =    120.0;         // QO */
	    /* PARAMS->geophs[ 5 ] =    78.0;          // SO */
	    /* PARAMS->geophs[ 6 ] =    6378.135;      // ER */
	    /* PARAMS->geophs[ 7 ] =    1.0;           // AE */

	    /* Read in the next two lines from the text file that contains the TLEs. */
	  //newstructure
/* 	  strcpy(sat_tle_file_temp, OPTIONS->dir_input_tle); */
/* 	  strcat(sat_tle_file_temp,"/"); */
	  strcpy(sat_tle_file_temp,"");
	  //newstructure
	    if ( OPTIONS->one_file_with_tle_of_all_sc == 0) { // one tle per file
	      strcat(sat_tle_file_temp,OPTIONS->tle_initialisation[ii]);
	      sat_tle_file = fopen(sat_tle_file_temp,"r");
	    }
	    else{ // all tles in one file
	      strcat(sat_tle_file_temp,OPTIONS->tle_initialisation[0]);

	      sat_tle_file = fopen(sat_tle_file_temp,"r");
	      for (bbb = 0; bbb < ii; bbb++){ // skip the TLEs of the sc before the current sc
		getline( &line_temp, &len, sat_tle_file ); 
		getline( &line_temp, &len, sat_tle_file );
	      }
	    }
	  
	    // First line
	    getline( &line_temp, &len, sat_tle_file );
	    lineln = strlen(line_temp)-1;
	    SpiceChar line[2][lineln];
	    for (pp = 0; pp<lineln-1 ; pp++)
	      line[0][pp] = line_temp[pp];
	    line[0][ lineln-1 ] = '\0';
	    // Second line
	    getline( &line_temp, &len, sat_tle_file );
	    for (pp = 0; pp<lineln-1 ; pp++)
	      line[1][pp] = line_temp[pp];
	    line[1][ lineln-1 ] = '\0';

	    fclose(sat_tle_file);
	    // Convert the elements of the TLE into "elems" and "epoch" that can be then read by the SPICE routine ev2lin_ to convert into the inertial state
	    //	    zzgetelm(frstyr, line, &epoch_sat,elems);
	      	    getelm_c( frstyr, lineln, line, &epoch_sat, CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.elems );

	    CONSTELLATION->spacecraft[ii][0].INTEGRATOR.bstar = CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.elems[2];


	    // Now propagate the state using ev2lin_ to the epoch of interest
	    extern /* Subroutine */ int ev2lin_(SpiceDouble *, SpiceDouble *,
						SpiceDouble *, SpiceDouble *);

	    ev2lin_( &epoch_sat, PARAMS->geophs, CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.elems, state );
	    CONSTELLATION->spacecraft[ii][eee].et = epoch_sat;
	    CONSTELLATION->spacecraft[ii][eee].et_sc_initial = epoch_sat;
	    

			      

	    static doublereal  precm[36];
	        static doublereal invprc[36]	/* was [6][6] */;
		    extern /* Subroutine */ int zzteme_(doublereal *, doublereal *);
		        extern /* Subroutine */ int invstm_(doublereal *, doublereal *);
			    extern /* Subroutine */ int  mxvg_(
	    doublereal *, doublereal *, integer *, integer *, doublereal *);

	        zzteme_(&CONSTELLATION->spacecraft[ii][eee].et, precm);

/*     ...now convert STATE to J2000. Invert the state transformation */
/*     operator (important to correctly do this). */

    invstm_(precm, invprc);
    static integer c__6 = 6;
        static doublereal tmpsta[6];
	    mxvg_(invprc, state, &c__6, &c__6, tmpsta);
    moved_(tmpsta, &c__6, state);
	    for (pp = 0; pp<3; pp++){
	      //	      	      	      printf("%f ", state[pp]);
	      CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[pp] = state[pp];
	      CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[pp] = state[pp+3];
	    }
	    	    /* printf("\n"); */
	    	    /* for (pp = 0; pp<3; pp++){ */
	      	    /*   	      printf("%f ", state[pp+3]); */
		    /* } */
	    /* printf("\n"); */
	    /* MPI_Finalize();exit(0); */
	    
	    // Initialize the Keplerian Elements
	    cart2kep( &OE_temp, CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].et,  PARAMS->EARTH.GRAVITY.mu );

	    CONSTELLATION->spacecraft[ii][eee].OE.sma          =OE_temp.sma;
	    CONSTELLATION->spacecraft[ii][eee].OE.eccentricity =OE_temp.eccentricity;
	    CONSTELLATION->spacecraft[ii][eee].OE.inclination  =OE_temp.inclination;
	    CONSTELLATION->spacecraft[ii][eee].OE.long_an      =OE_temp.long_an;
	    CONSTELLATION->spacecraft[ii][eee].OE.w            =OE_temp.w;
	    CONSTELLATION->spacecraft[ii][eee].OE.f            =OE_temp.f;
	    CONSTELLATION->spacecraft[ii][eee].OE.tp           =OE_temp.tp;
	    CONSTELLATION->spacecraft[ii][eee].OE.E            =OE_temp.E;
	    CONSTELLATION->spacecraft[ii][eee].OE.ra            =OE_temp.ra;
	    CONSTELLATION->spacecraft[ii][eee].OE.initial_an_to_sc = OE_temp.an_to_sc;
	    CONSTELLATION->spacecraft[ii][eee].OE.w_ave_temp = OE_temp.w;
	    CONSTELLATION->spacecraft[ii][eee].OE.w_ave = 9999.9999/RAD2DEG;
	    CONSTELLATION->spacecraft[ii][eee].OE.sma_ave_temp = OE_temp.sma;
	    CONSTELLATION->spacecraft[ii][eee].OE.sma_ave = 9999.9999;
	    CONSTELLATION->spacecraft[ii][eee].OE.ecc_ave_temp = OE_temp.eccentricity;
	    CONSTELLATION->spacecraft[ii][eee].OE.ecc_ave = 9999.9999;

	    CONSTELLATION->spacecraft[ii][eee].OE.ave_increm = 1;
	    CONSTELLATION->spacecraft[ii][eee].et_last_orbit = CONSTELLATION->spacecraft[ii][eee].et;
	    CONSTELLATION->spacecraft[ii][eee].orbit_number = 0;
    
	    /* printf("ecc = %e | %e\n", CONSTELLATION->spacecraft[ii][eee].OE.eccentricity, elems[5]); */
	    /* exit(0); */

	    radius_perigee = CONSTELLATION->spacecraft[ii][eee].OE.sma * ( 1 - CONSTELLATION->spacecraft[ii][eee].OE.eccentricity );
	    if (radius_perigee < PARAMS->EARTH.radius){
	      printf("***! The orbit of satellite %d intersects the Earth (altitude of perigee = %f km). The program will stop. !***\n",ii, radius_perigee-PARAMS->EARTH.radius);
	      MPI_Finalize();
	      exit(0);
	    }
	    if (iDebugLevel >= 1){
	      if (iProc == 0) printf("-- (initialize_constellation) Done initializing et, r_i2cg_INRTL, v_i2cg_INRTL, orbital elements from TLE.\n");
	    }

	  }

	  /*************************************************************************************/
	  /*************************************************************************************/
	  /********** INITIALIZE: et, r_i2cg_INRTL, v_i2cg_INRTL, OE FOR MAIN SC BY COE ********/
	  /*************************************************************************************/
	  /*************************************************************************************/

	  if ( strcmp( OPTIONS->type_orbit_initialisation, "oe" ) == 0 ){
	    if (iDebugLevel >= 1){
	      if (iProc == 0) printf("-- (initialize_constellation) Initializing et, r_i2cg_INRTL, v_i2cg_INRTL, orbital elements from orbital elements.\n");
	    }

	    CONSTELLATION->spacecraft[ii][eee].et = et;
	    CONSTELLATION->spacecraft[ii][eee].et_sc_initial = et;
	    // Initialize the Keplerian Elements
	    CONSTELLATION->spacecraft[ii][eee].OE.sma            = ( PARAMS->EARTH.radius + OPTIONS->apogee_alt[ii] ) / ( 1 + OPTIONS->eccentricity[ii] ); // semi-major axis (CBV)

	    //    printf("sma = %f\n", CONSTELLATION->spacecraft[ii][eee].OE.sma);
	    radius_perigee = CONSTELLATION->spacecraft[ii][eee].OE.sma * ( 1 - OPTIONS->eccentricity[ii] );
	    if (radius_perigee < PARAMS->EARTH.radius){
	      printf("***! The orbit of satellite %d intersects the Earth (altitude of perigee = %f km). The program will stop. !***\n",ii, radius_perigee-PARAMS->EARTH.radius);
	      MPI_Finalize();
	      exit(0);
	    }
    
	    CONSTELLATION->spacecraft[ii][eee].OE.inclination    = OPTIONS->inclination[ii] * DEG2RAD;
	    CONSTELLATION->spacecraft[ii][eee].OE.w              = OPTIONS->w[ii] * DEG2RAD;                                // argument of perigee (CBV)
	    CONSTELLATION->spacecraft[ii][eee].OE.long_an        = OPTIONS->long_an[ii] * DEG2RAD;                          // RAAN (CBV)
	    CONSTELLATION->spacecraft[ii][eee].OE.f              = OPTIONS->f[ii] * DEG2RAD;                                // true anomaly (CBV)
	    CONSTELLATION->spacecraft[ii][eee].OE.eccentricity   = OPTIONS->eccentricity[ii];
	    CONSTELLATION->spacecraft[ii][eee].OE.initial_an_to_sc =  fmod(CONSTELLATION->spacecraft[ii][eee].OE.w + CONSTELLATION->spacecraft[ii][eee].OE.f, 2*M_PI);

	    // Initialize the inertial state
	    kep2cart(   CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL,
			CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL,
			&PARAMS->EARTH.GRAVITY.mu,
			&CONSTELLATION->spacecraft[ii][eee].OE); // Computes the ECI coordinates based on the Keplerian inputs (orbital elements and mu) (in propagator.c) (CBV)

	    // !!!!!!!!!!!!!!!!!!!!!! REMOVE TWO LINES BELOW!!!!!!!!!  

	    /* //	  for molniya */
	    /* CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[0] = -1529.8942870000; CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[1] =  -2672.8773570000; CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[2] =  -6150.1153400000; */
	    /* CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[0] =  8.7175180000; CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[1] =  -4.9897090000; CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[2] =  0.0000000000; */

	    /* 	  // for geo */
	    /* CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[0] = 36607.3548590000;  CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[1] =  -20921.7217480000; CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[2] = -0.0000000000; */
	    /* CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[0] = 1.5256360000;  CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[1] = 2.6694510000;  CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[2] =  0.0000000000; */

	    // for iss
	    /* CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[0] = -4467.3083710453;  CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[1] =  -5053.5032161387; CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[2] =  -427.6792780851; */
	    /* CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[0] = 3.8279470371;  CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[1] = -2.8757493806;  CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[2] =  -6.0045254246; */

	    // !!!!!!!!!!!!!!!!!!!!!! END OF REMOVE TWO LINES BELOW!!!!!!!!!

	    // Right ascension
	    CONSTELLATION->spacecraft[ii][eee].OE.ra =  atan2(CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[1], CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[0]);
	    if ( CONSTELLATION->spacecraft[ii][eee].OE.ra < 0){
	      CONSTELLATION->spacecraft[ii][eee].OE.ra = 2*M_PI + CONSTELLATION->spacecraft[ii][eee].OE.ra ;
	    }



	    // !!!!!!!!!!!!!!!!!!!!!!!!!! THE BLOCK BELOW IS TO INITIALIZE SATELLLITE 2 WITH THE CORRECT SPACING WITH RESPECT TO SATELLITE 1. IT IS JUST FOR A TRY, AND SHOULD BE REMOVED AND IMPLEMENTED IN THE CODE ITSELF. SO IF IT'S STILL HERE AFTER JUNE 10TH, 2016 THEN REMOVE IT!
	    /* if ( strcmp( OPTIONS->type_orbit_initialisation, "oe" ) == 0 ){ */
	    /*   if ( ii > 0){ */
	    /*     v_copy( CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL, CONSTELLATION->spacecraft[0][0].r_i2cg_INRTL); */
	    /*     v_copy( CONSTELLATION->spacecraft[ii][0].v_i2cg_INRTL, CONSTELLATION->spacecraft[0][0].v_i2cg_INRTL); */

	    /*     CONSTELLATION->spacecraft[ii][0].OE.sma            = ( PARAMS->EARTH.radius + OPTIONS->apogee_alt[0] ) / ( 1 + OPTIONS->eccentricity[0] ); // semi-major axis (CBV) */

	    /*     //    printf("sma = %f\n", CONSTELLATION->spacecraft[ii][0].OE.sma); */
	    /*     radius_perigee = CONSTELLATION->spacecraft[ii][0].OE.sma * ( 1 - OPTIONS->eccentricity[0] ); */
	    /*     if (radius_perigee < PARAMS->EARTH.radius){ */
	    /*       printf("The orbit of satellite %d intersects the Earth (altitude of perigee = %f km). The program will stop. !***\n",ii, radius_perigee-PARAMS->EARTH.radius); */
	    /*       exit(0); */
	    /*     } */
    
	    /*     CONSTELLATION->spacecraft[ii][0].OE.inclination    = OPTIONS->inclination[0] * DEG2RAD; */
	    /*     CONSTELLATION->spacecraft[ii][0].OE.w              = OPTIONS->w[0] * DEG2RAD;                                // argument of perigee (CBV) */
	    /*     CONSTELLATION->spacecraft[ii][0].OE.long_an        = OPTIONS->long_an[0] * DEG2RAD;                          // RAAN (CBV) */
	    /*     CONSTELLATION->spacecraft[ii][0].OE.f              = OPTIONS->f[0] * DEG2RAD;                                // true anomaly (CBV) */
	    /*     CONSTELLATION->spacecraft[ii][0].OE.eccentricity   = OPTIONS->eccentricity[0]; */
	    /*     CONSTELLATION->spacecraft[ii][0].OE.ra =  atan2(CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL[1], CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL[0]); */
	    /*     if ( CONSTELLATION->spacecraft[ii][0].OE.ra < 0){ */
	    /*       CONSTELLATION->spacecraft[ii][0].OE.ra = 2*M_PI + CONSTELLATION->spacecraft[ii][0].OE.ra ; */
	    /*     } */

	    /*   } */
	    /* } */
	    // !!!!!!!!!!!!!!!!!!!!!!!!!! END OF THE BLOCK BELOW IS TO INITIALIZE SATELLLITE 2 WITH THE CORRECT SPACING WITH RESPECT TO SATELLITE 1. IT IS JUST FOR A TRY, AND SHOULD BE REMOVED AND IMPLEMENTED IN THE CODE ITSELF. SO IF IT'S STILL HERE AFTER JUNE 10TH, 2016 THEN REMOVE IT!

	    if (iDebugLevel >= 1){
	      if (iProc == 0) printf("-- (initialize_constellation) Done initializing et, r_i2cg_INRTL, v_i2cg_INRTL, orbital elements from orbital elements.\n");
	    }

	  } // end of initiliazing et, r_i2cg_INRTL, v_i2cg_INRTL, OE FOR MAIN SC BY COE

	  /*************************************************************************************/
	  /*************************************************************************************/
	  /** INITIALIZE: et, r_i2cg_INRTL, v_i2cg_INRTL, r_ecef2cg_ECEF, v_ecef2cg_ECEF, OE FOR MAIN SC BY ECEF STATE (position and velocity in ECEF) **/
	  /*************************************************************************************/
	  /*************************************************************************************/

	  if ( strcmp( OPTIONS->type_orbit_initialisation, "state_ecef" ) == 0 ){
	    if (iDebugLevel >= 1){
	      if (iProc == 0) printf("-- (initialize_constellation) Initializing et, r_i2cg_INRTL, v_i2cg_INRTL, r_ecef2cg_ECEF, v_ecef2cg_ECEF, orbital elements from ECEF state (position and velocity).\n");
	    }

	    CONSTELLATION->spacecraft[ii][eee].et = et;
	    CONSTELLATION->spacecraft[ii][eee].et_sc_initial = et;

	    // ECEF position and velocity
	    CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF[0] = OPTIONS->x_ecef[ii];
	    CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF[1] = OPTIONS->y_ecef[ii];
	    CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF[2] = OPTIONS->z_ecef[ii];
	    CONSTELLATION->spacecraft[ii][eee].v_ecef2cg_ECEF[0] = OPTIONS->vx_ecef[ii];
	    CONSTELLATION->spacecraft[ii][eee].v_ecef2cg_ECEF[1] = OPTIONS->vy_ecef[ii];
	    CONSTELLATION->spacecraft[ii][eee].v_ecef2cg_ECEF[2] = OPTIONS->vz_ecef[ii];

	    // ECI position and velocity
	    estate[0] = CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF[0];estate[1] = CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF[1];estate[2] = CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF[2];
	    estate[3] = CONSTELLATION->spacecraft[ii][eee].v_ecef2cg_ECEF[0];estate[4] = CONSTELLATION->spacecraft[ii][eee].v_ecef2cg_ECEF[1];estate[5] = CONSTELLATION->spacecraft[ii][eee].v_ecef2cg_ECEF[2];
	    sxform_c (  PARAMS->EARTH.earth_fixed_frame,  "J2000", CONSTELLATION->spacecraft[ii][eee].et,    xform  ); 
	    mxvg_c   (  xform,       estate,   6,  6, jstate ); 
	    CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[0] = jstate[0]; CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[1] = jstate[1]; CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[2] = jstate[2];
	    CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[0] = jstate[3]; CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[1] = jstate[4]; CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[2] = jstate[5];
	    /* pxform_c( PARAMS->EARTH.earth_fixed_frame, "J2000", CONSTELLATION->spacecraft[ii][eee].et, T_ECEF_to_J2000); */
	    /* m_x_v(CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, T_ECEF_to_J2000, CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF); */
	    /* m_x_v(CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL, T_ECEF_to_J2000, CONSTELLATION->spacecraft[ii][eee].v_ecef2cg_ECEF); */

	    // Orbital elements
	    cart2kep( &CONSTELLATION->spacecraft[ii][eee].OE, CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].et,  PARAMS->EARTH.GRAVITY.mu);

	    CONSTELLATION->spacecraft[ii][eee].OE.initial_an_to_sc = CONSTELLATION->spacecraft[ii][eee].OE.an_to_sc;
	    CONSTELLATION->spacecraft[ii][eee].OE.w_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.w;
	    CONSTELLATION->spacecraft[ii][eee].OE.w_ave = 9999.999999/RAD2DEG;
	    CONSTELLATION->spacecraft[ii][eee].OE.sma_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.sma;
	    CONSTELLATION->spacecraft[ii][eee].OE.sma_ave = 9999.999999;
	    CONSTELLATION->spacecraft[ii][eee].OE.ecc_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.eccentricity;
	    CONSTELLATION->spacecraft[ii][eee].OE.ecc_ave = 9999.999999;

	    CONSTELLATION->spacecraft[ii][eee].OE.ave_increm = 1;
	    CONSTELLATION->spacecraft[ii][eee].et_last_orbit = CONSTELLATION->spacecraft[ii][eee].et;
	    CONSTELLATION->spacecraft[ii][eee].orbit_number = 0;


	    radius_perigee = CONSTELLATION->spacecraft[ii][eee].OE.sma * ( 1 - CONSTELLATION->spacecraft[ii][eee].OE.eccentricity  );
	    if (radius_perigee < PARAMS->EARTH.radius){
	      printf("***! The orbit of satellite %d intersects the Earth (altitude of perigee = %f km). The program will stop. !***\n",ii, radius_perigee-PARAMS->EARTH.radius);
	      MPI_Finalize();
	      exit(0);
	    }

	    if (iDebugLevel >= 1){
	      if (iProc == 0) printf("-- (initialize_constellation) Done initializing et, r_i2cg_INRTL, v_i2cg_INRTL, r_ecef2cg_ECEF, v_ecef2cg_ECEF, orbital elements from ECEF state (position and velocity).\n");
	    }

	  } // end of initializing et, r_i2cg_INRTL, v_i2cg_INRTL, r_ecef2cg_ECEF, v_ecef2cg_ECEF, OE FOR MAIN SC BY ECEF STATE


	  /*************************************************************************************/
	  /*************************************************************************************/
	  /** INITIALIZE: et, r_i2cg_INRTL, v_i2cg_INRTL, OE FOR MAIN SC BY ECI STATE (position and velocity in ECI) **/
	  /*************************************************************************************/
	  /*************************************************************************************/

	  if ( strcmp( OPTIONS->type_orbit_initialisation, "state_eci" ) == 0 ){
	    if (iDebugLevel >= 1){
	      if (iProc == 0) printf("-- (initialize_constellation) Initializing et, r_i2cg_INRTL, v_i2cg_INRTL, orbital elements from ECI state (position and velocity).\n");
	    }

	    CONSTELLATION->spacecraft[ii][eee].et = et;
	    CONSTELLATION->spacecraft[ii][eee].et_sc_initial = et;

	    // ECI position and velocity
	    CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[0] = OPTIONS->x_eci[ii];
	    CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[1] = OPTIONS->y_eci[ii];
	    CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[2] = OPTIONS->z_eci[ii];
	    CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[0] = OPTIONS->vx_eci[ii];
	    CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[1] = OPTIONS->vy_eci[ii];
	    CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[2] = OPTIONS->vz_eci[ii];

	    // Orbital elements
	    cart2kep( &CONSTELLATION->spacecraft[ii][eee].OE, CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].et,  PARAMS->EARTH.GRAVITY.mu);

	    CONSTELLATION->spacecraft[ii][eee].OE.initial_an_to_sc = CONSTELLATION->spacecraft[ii][eee].OE.an_to_sc;
	    
	    CONSTELLATION->spacecraft[ii][eee].OE.w_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.w;
	    CONSTELLATION->spacecraft[ii][eee].OE.w_ave = 9999.999999/RAD2DEG;
	    CONSTELLATION->spacecraft[ii][eee].OE.sma_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.sma;
	    CONSTELLATION->spacecraft[ii][eee].OE.sma_ave = 9999.999999;
	    CONSTELLATION->spacecraft[ii][eee].OE.ecc_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.eccentricity;
	    CONSTELLATION->spacecraft[ii][eee].OE.ecc_ave = 9999.999999;

	    CONSTELLATION->spacecraft[ii][eee].OE.ave_increm = 1;
	    CONSTELLATION->spacecraft[ii][eee].et_last_orbit = CONSTELLATION->spacecraft[ii][eee].et;
	    CONSTELLATION->spacecraft[ii][eee].orbit_number = 0;


	    radius_perigee = CONSTELLATION->spacecraft[ii][eee].OE.sma * ( 1 - CONSTELLATION->spacecraft[ii][eee].OE.eccentricity  );
	    if (radius_perigee < PARAMS->EARTH.radius){
	      printf("***! The orbit of satellite %d intersects the Earth (altitude of perigee = %f km). The program will stop. !***\n",ii, radius_perigee-PARAMS->EARTH.radius);
	      MPI_Finalize();
	      exit(0);
	    }

	    if (iDebugLevel >= 1){
	      if (iProc == 0) printf("-- (initialize_constellation) Done initializing et, r_i2cg_INRTL, v_i2cg_INRTL, r_ecef2cg_ECEF, v_ecef2cg_ECEF from ECI state (position and velocity).\n");
	    }

	  } // end of initializing et, r_i2cg_INRTL, v_i2cg_INRTL, r_ecef2cg_ECEF, v_ecef2cg_ECEF, OE FOR MAIN SC BY ECI STATE


	  /*************************************************************************************/
	  /*************************************************************************************/
	  /** INITIALIZE: et, r_i2cg_INRTL, v_i2cg_INRTL, OE FOR MAIN SC BY ECI STATE FROM COLLISION INPUT FILE (position and velocity in ECI) **/
	  /*************************************************************************************/
	  /*************************************************************************************/

	  if ( (strcmp( OPTIONS->type_orbit_initialisation, "collision" ) == 0 ) || (strcmp( OPTIONS->type_orbit_initialisation, "collision_vcm" ) == 0 )){
	    if (iDebugLevel >= 1){
	      if (iProc == 0) printf("-- (initialize_constellation) Initializing et, r_i2cg_INRTL, v_i2cg_INRTL, orbital elements from ECI state reead from the input collision file (position and velocity).\n");  }


	    if (strcmp( OPTIONS->type_orbit_initialisation, "collision" ) == 0 ){
	    CONSTELLATION->spacecraft[ii][eee].et_sc_initial = et;
	    CONSTELLATION->spacecraft[ii][eee].et = et;
	    }
	    else{
	    CONSTELLATION->spacecraft[ii][eee].et = OPTIONS->et_vcm[ii];
	      CONSTELLATION->spacecraft[ii][eee].et_sc_initial = OPTIONS->et_vcm[ii]; // the eoch time of the two VCMs are different
/* 	      if ( (OPTIONS->et_vcm[ii] > OPTIONS->swpc_et_first_prediction) && (OPTIONS->swpc_need_predictions == 1)){// overwrite aaa_mod written before */
/* previous_index( &OPTIONS->aaa_mod, OPTIONS->et_mod_f107_ap, OPTIONS->et_interpo[0] - OPTIONS->swpc_et_first_prediction, ( OPTIONS->nb_time_steps + (int)(( 2 * 24 * 3600. ) / OPTIONS->dt)));  */
/* 	      CONSTELLATION->aaa_mod[ii] =  OPTIONS->aaa_mod + OPTIONS->et_vcm[ii] - CONSTE; */
/* 	      } */
/* 	      else{ */
/* 		CONSTELLATION->aaa_mod[ii] =  0; */
/* 		" */
	    }

	    // ECI position and velocity
	    CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[0] = OPTIONS->x_eci[ii];
	    CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[1] = OPTIONS->y_eci[ii];
	    CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[2] = OPTIONS->z_eci[ii];
	    CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[0] = OPTIONS->vx_eci[ii];
	    CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[1] = OPTIONS->vy_eci[ii];
	    CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[2] = OPTIONS->vz_eci[ii];

	    // Orbital elements
	    cart2kep( &CONSTELLATION->spacecraft[ii][eee].OE, CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].et,  PARAMS->EARTH.GRAVITY.mu);

	    CONSTELLATION->spacecraft[ii][eee].OE.initial_an_to_sc = CONSTELLATION->spacecraft[ii][eee].OE.an_to_sc;
	    
	    CONSTELLATION->spacecraft[ii][eee].OE.w_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.w;
	    CONSTELLATION->spacecraft[ii][eee].OE.w_ave = 9999.999999/RAD2DEG;
	    CONSTELLATION->spacecraft[ii][eee].OE.sma_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.sma;
	    CONSTELLATION->spacecraft[ii][eee].OE.sma_ave = 9999.999999;
	    CONSTELLATION->spacecraft[ii][eee].OE.ecc_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.eccentricity;
	    CONSTELLATION->spacecraft[ii][eee].OE.ecc_ave = 9999.999999;

	    CONSTELLATION->spacecraft[ii][eee].OE.ave_increm = 1;
	    CONSTELLATION->spacecraft[ii][eee].et_last_orbit = CONSTELLATION->spacecraft[ii][eee].et;
	    CONSTELLATION->spacecraft[ii][eee].orbit_number = 0;


	    radius_perigee = CONSTELLATION->spacecraft[ii][eee].OE.sma * ( 1 - CONSTELLATION->spacecraft[ii][eee].OE.eccentricity  );
	    if (radius_perigee < PARAMS->EARTH.radius){
	      printf("***! The orbit of satellite %d intersects the Earth (altitude of perigee = %f km). The program will stop. !***\n",ii, radius_perigee-PARAMS->EARTH.radius);
	      MPI_Finalize();
	      exit(0);
	    }

	    if (iDebugLevel >= 1){
	      if (iProc == 0) printf("-- (initialize_constellation) Initializing et, r_i2cg_INRTL, v_i2cg_INRTL, orbital elements from ECI state reead from the input collision file (position and velocity).\n");
	    }

	  } // end of initializing et, r_i2cg_INRTL, v_i2cg_INRTL, r_ecef2cg_ECEF, v_ecef2cg_ECEF, OE FOR MAIN SC BY ECI STATE FROM COLLISION INPUT FILE


	  /*************************************************************************************/
	  /*************************************************************************************/
	  /** INITIALIZE: et, r_i2cg_INRTL, v_i2cg_INRTL, OE FOR MAIN SC BY DEPLOYMENT **/
	  /*************************************************************************************/
	  /*************************************************************************************/

	  if ( strcmp( OPTIONS->type_orbit_initialisation, "deployment" ) == 0 ){
	    if (iDebugLevel >= 1){
	      if (iProc == 0) printf("-- (initialize_constellation) Initializing et, r_i2cg_INRTL, v_i2cg_INRTL, orbital elements from a deployment module.\n");
	    }

	    CONSTELLATION->spacecraft[ii][eee].et = et;
	    CONSTELLATION->spacecraft[ii][eee].et_sc_initial = et;
	  
	    // Convert the deployment module orbital elements in ECI coordinates
	    ORBITAL_ELEMENTS_T deployment_module_orbital_elements;
	    deployment_module_orbital_elements.sma            = ( PARAMS->EARTH.radius + OPTIONS->deployment_module_orbital_elements[ii][0] ) / ( 1 + OPTIONS->deployment_module_orbital_elements[ii][5] ); 
	    radius_perigee = deployment_module_orbital_elements.sma * ( 1 - OPTIONS->deployment_module_orbital_elements[ii][5] );
	    if (radius_perigee < PARAMS->EARTH.radius){
	      printf("***! The orbit of satellite %d intersects the Earth (altitude of perigee = %f km). The program will stop. !***\n",ii, radius_perigee-PARAMS->EARTH.radius);
	      MPI_Finalize();
	      exit(0);
	    }    
	    deployment_module_orbital_elements.inclination    = OPTIONS->deployment_module_orbital_elements[ii][1] * DEG2RAD;
	    deployment_module_orbital_elements.w              = OPTIONS->deployment_module_orbital_elements[ii][2] * DEG2RAD; 
	    deployment_module_orbital_elements.long_an        = OPTIONS->deployment_module_orbital_elements[ii][3] * DEG2RAD; 
	    deployment_module_orbital_elements.f              = OPTIONS->deployment_module_orbital_elements[ii][4] * DEG2RAD; 
	    deployment_module_orbital_elements.eccentricity   = OPTIONS->deployment_module_orbital_elements[ii][5];
	    //	  printf("%f %f %f %f %f %f %f\n", OPTIONS->deployment_module_orbital_elements[ii][0], deployment_module_orbital_elements.sma,deployment_module_orbital_elements.inclination*RAD2DEG,deployment_module_orbital_elements.w,deployment_module_orbital_elements.long_an,deployment_module_orbital_elements.f,deployment_module_orbital_elements.eccentricity);
	    double eci_position_deployment_module[3]; double eci_velocity_deployment_module[3];
	    kep2cart(eci_position_deployment_module, eci_velocity_deployment_module, &PARAMS->EARTH.GRAVITY.mu, &deployment_module_orbital_elements);

	    // Convert the deployment velocity of the satellite from the deployment module LVLH coordinate system to ECI velocity
	    double velocity_satellte_deployed_lvlh[3]; double velocity_satellte_deployed_eci[3];
	    velocity_satellte_deployed_lvlh[0] = OPTIONS->deployment_speed[ii] * cos(OPTIONS->deployment_angle[ii]*DEG2RAD) / 1000.; // / 1000. because the input is in m/s but we want km/s
	    velocity_satellte_deployed_lvlh[1] = -OPTIONS->deployment_speed[ii] * sin(OPTIONS->deployment_angle[ii]*DEG2RAD) / 1000.; // "-" because LVLH _Y points to the left and we are considering that an angle of for example 30 degrees is counted from the in-track direction to the right (clockwise). If the sc is ejected to the left with an angle of for example 90 degrees, then deployment_angle is equal to 270 degrees
	    velocity_satellte_deployed_lvlh[2] = 0; // for now we assume the satellites are deployed in the LVLH_X/Y plane
	    double T_inrtl_2_lvlh_for_deployment[3][3]; 	  double T_lvlh_2_inrtl_for_deployment[3][3];
	    compute_T_inrtl_2_lvlh(T_inrtl_2_lvlh_for_deployment, eci_position_deployment_module, eci_velocity_deployment_module);
	    m_trans(T_lvlh_2_inrtl_for_deployment, T_inrtl_2_lvlh_for_deployment);
	    m_x_v(velocity_satellte_deployed_eci, T_lvlh_2_inrtl_for_deployment, velocity_satellte_deployed_lvlh);
	  
	    // Add deployment velocity of the satellite to the velocity of the deployment module (in ECI)
	    v_add(CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL, eci_velocity_deployment_module, velocity_satellte_deployed_eci);

	    // ECI position of the satellite is the same as the ECI position of the deployment module
	    v_copy(CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, eci_position_deployment_module);

	    // Right ascension
	    CONSTELLATION->spacecraft[ii][eee].OE.ra =  atan2(CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[1], CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[0]);
	    if ( CONSTELLATION->spacecraft[ii][eee].OE.ra < 0){
	      CONSTELLATION->spacecraft[ii][eee].OE.ra = 2*M_PI + CONSTELLATION->spacecraft[ii][eee].OE.ra ;
	    }

	    // Orbital elements
	    cart2kep( &CONSTELLATION->spacecraft[ii][eee].OE, CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].et,  PARAMS->EARTH.GRAVITY.mu);

	    CONSTELLATION->spacecraft[ii][eee].OE.initial_an_to_sc = CONSTELLATION->spacecraft[ii][eee].OE.an_to_sc;
	    
	    CONSTELLATION->spacecraft[ii][eee].OE.w_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.w;
	    CONSTELLATION->spacecraft[ii][eee].OE.w_ave = 9999.999999/RAD2DEG;
	    CONSTELLATION->spacecraft[ii][eee].OE.sma_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.sma;
	    CONSTELLATION->spacecraft[ii][eee].OE.sma_ave = 9999.999999;
	    CONSTELLATION->spacecraft[ii][eee].OE.ecc_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.eccentricity;
	    CONSTELLATION->spacecraft[ii][eee].OE.ecc_ave = 9999.999999;

	    CONSTELLATION->spacecraft[ii][eee].OE.ave_increm = 1;
	    CONSTELLATION->spacecraft[ii][eee].et_last_orbit = CONSTELLATION->spacecraft[ii][eee].et;
	    CONSTELLATION->spacecraft[ii][eee].orbit_number = 0;


	    radius_perigee = CONSTELLATION->spacecraft[ii][eee].OE.sma * ( 1 - CONSTELLATION->spacecraft[ii][eee].OE.eccentricity  );
	    if (radius_perigee < PARAMS->EARTH.radius){
	      printf("***! The orbit of satellite %d intersects the Earth (altitude of perigee = %f km). The program will stop. !***\n",ii, radius_perigee-PARAMS->EARTH.radius);
	      MPI_Finalize();
	      exit(0);
	    }

	    if (iDebugLevel >= 1){
	      if (iProc == 0) printf("-- (initialize_constellation) Done initializing et, r_i2cg_INRTL, v_i2cg_INRTL, orbital elements from a deployment module.\n");
	    }

	  } // end of initializing et, r_i2cg_INRTL, v_i2cg_INRTL, OE FOR MAIN SC BY DEPLOYMENT



	
	  for (ccc = 0; ccc < nProcs; ccc++){ // the iProc that runs main sc ii sends to the other iProc (that do not run main sc ii) the orbital elements of main sc ii
	    if (ccc != iProc){
	      MPI_Send(&CONSTELLATION->spacecraft[ii][0].OE.sma, 1, MPI_DOUBLE, ccc, 0, MPI_COMM_WORLD);
	      MPI_Send(&CONSTELLATION->spacecraft[ii][0].OE.eccentricity , 1, MPI_DOUBLE, ccc, 0, MPI_COMM_WORLD);
	      MPI_Send(&CONSTELLATION->spacecraft[ii][0].OE.inclination , 1, MPI_DOUBLE, ccc, 0, MPI_COMM_WORLD);
	      MPI_Send(&CONSTELLATION->spacecraft[ii][0].OE.w , 1, MPI_DOUBLE, ccc, 0, MPI_COMM_WORLD);
	      MPI_Send(&CONSTELLATION->spacecraft[ii][0].OE.long_an , 1, MPI_DOUBLE, ccc, 0, MPI_COMM_WORLD);
	      MPI_Send(&CONSTELLATION->spacecraft[ii][0].OE.f , 1, MPI_DOUBLE, ccc, 0, MPI_COMM_WORLD);
	      //   printf("iProc %d sends to iProc %d for ii = %d\n", iProc, ccc, ii);
	
	    }
	  }

	  /* 	if (nProcs > 1){ */
	  /* 	  for (ccc = 0; ccc < nProcs; ccc++){ */
	  /* 	    MPI_Send(&CONSTELLATION->spacecraft[ii][0].OE.sma, 1, MPI_DOUBLE, ccc, 0, MPI_COMM_WORLD); */
	  /* 	    MPI_Send(&CONSTELLATION->spacecraft[ii][0].OE.eccentricity , 1, MPI_DOUBLE, ccc, 0, MPI_COMM_WORLD); */
	  /* 	    MPI_Send(&CONSTELLATION->spacecraft[ii][0].OE.inclination , 1, MPI_DOUBLE, ccc, 0, MPI_COMM_WORLD); */
	  /* 	    MPI_Send(&CONSTELLATION->spacecraft[ii][0].OE.w , 1, MPI_DOUBLE, ccc, 0, MPI_COMM_WORLD); */
	  /* 	    MPI_Send(&CONSTELLATION->spacecraft[ii][0].OE.long_an , 1, MPI_DOUBLE, ccc, 0, MPI_COMM_WORLD); */
	  /* 	    MPI_Send(&CONSTELLATION->spacecraft[ii][0].OE.f , 1, MPI_DOUBLE, ccc, 0, MPI_COMM_WORLD); */
	  /* 	    //	  printf("Sent to %d (out of %d)\n", ccc, nProcs); */
	  /* 	  } */
	  /* 	} */



	  /* 	if (iProc == 0) printf("BEFORE\n"); */
	  /*     MPI_Barrier(MPI_COMM_WORLD); */
	  /*     if (iProc == 0) printf("RIG after\n"); */
	  /* 	  MPI_Bcast(&CONSTELLATION->spacecraft[ii][0].OE.sma, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD); */
	  /* 	  MPI_Bcast(&CONSTELLATION->spacecraft[ii][0].OE.eccentricity , 1, MPI_DOUBLE,  0, MPI_COMM_WORLD); */
	  /* 	  MPI_Bcast(&CONSTELLATION->spacecraft[ii][0].OE.inclination , 1, MPI_DOUBLE,  0, MPI_COMM_WORLD); */
	  /* 	  MPI_Bcast(&CONSTELLATION->spacecraft[ii][0].OE.w , 1, MPI_DOUBLE,  0, MPI_COMM_WORLD); */
	  /* 	  MPI_Bcast(&CONSTELLATION->spacecraft[ii][0].OE.long_an , 1, MPI_DOUBLE,  0, MPI_COMM_WORLD); */
	  /* 	  MPI_Bcast(&CONSTELLATION->spacecraft[ii][0].OE.f , 1, MPI_DOUBLE,  0, MPI_COMM_WORLD); */
	  /* 	if (iProc == 0) printf("AFTER\n"); */
	  if (iDebugLevel >= 1){
	    if (iProc == 0) printf("-- (initialize_constellation) Done initiliazing et, r_i2cg_INRTL, v_i2cg_INRTL, OE FOR MAIN SC by iProc 0 (voluntarily prnting all iProc here). (iProc %d)\n", iProc); 
	  }
	}// end of initiliazing et, r_i2cg_INRTL, v_i2cg_INRTL, OE FOR MAIN SC BY COE and TLE

	/*************************************************************************************/
	/*************************************************************************************/
	/*************************************************************************************/
	/************** INITIALIZE: et, r_i2cg_INRTL, v_i2cg_INRTL, OE, number_of_collisions FOR ENSEMBLE SC ***********/
	/*************************************************************************************/
	/*************************************************************************************/
	/*************************************************************************************/

	else if ( ( start_ensemble[ii] != 0 ) && ( eee == array_sc[start_ensemble[ii]] ) ){// if this iProc is not runnning the main sc, then receives the orbital elements of main sc ii, which was run (and sent) by iProc which_iproc_is_running_main_sc[ii] (eee = array_sc[start_ensemble[ii]] is so that you receive only once)

	  MPI_Recv(&CONSTELLATION->spacecraft[ii][0].OE.sma, 1, MPI_DOUBLE, which_iproc_is_running_main_sc[ii], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(&CONSTELLATION->spacecraft[ii][0].OE.eccentricity, 1, MPI_DOUBLE, which_iproc_is_running_main_sc[ii], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(&CONSTELLATION->spacecraft[ii][0].OE.inclination, 1, MPI_DOUBLE, which_iproc_is_running_main_sc[ii], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(&CONSTELLATION->spacecraft[ii][0].OE.w, 1, MPI_DOUBLE, which_iproc_is_running_main_sc[ii], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(&CONSTELLATION->spacecraft[ii][0].OE.long_an, 1, MPI_DOUBLE, which_iproc_is_running_main_sc[ii], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(&CONSTELLATION->spacecraft[ii][0].OE.f, 1, MPI_DOUBLE, which_iproc_is_running_main_sc[ii], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


	  //	printf("iProc %d receives from iProc %d for ii = %d\n", iProc, which_iproc_is_running_main_sc[ii], ii);
	}

	if (eee > 0){ // here eee > 0 so we are initializing et, r_i2cg_INRTL, v_i2cg_INRTL, OE, number_of_collisions for the sc from an ensemble
	  /*       	if (nProcs > 1){ */
	  /* 	  if ( iProc > 0){ */
	  /* 	    if (eee == start_ensemble + iProc * OPTIONS->nb_ensemble_min_per_proc){ // only receives this one time */
	  /*       MPI_Recv(&CONSTELLATION->spacecraft[ii][0].OE.sma, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); */
	  /*       MPI_Recv(&CONSTELLATION->spacecraft[ii][0].OE.eccentricity, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); */
	  /*       MPI_Recv(&CONSTELLATION->spacecraft[ii][0].OE.inclination, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); */
	  /*       MPI_Recv(&CONSTELLATION->spacecraft[ii][0].OE.w, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); */
	  /*       MPI_Recv(&CONSTELLATION->spacecraft[ii][0].OE.long_an, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); */
	  /*       MPI_Recv(&CONSTELLATION->spacecraft[ii][0].OE.f, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); */
	  /* 	    } */
	  /* 	  } */
	  /*       	} */
	  if (iDebugLevel >= 1){
	    if (iProc == 0) printf("-- (initialize_constellation) Initializing et, r_i2cg_INRTL, v_i2cg_INRTL for ensemble spacecraft.\n");
	  }

	  if (strcmp( OPTIONS->type_orbit_initialisation, "collision_vcm" ) == 0 ){
	      CONSTELLATION->spacecraft[ii][eee].et_sc_initial = OPTIONS->et_vcm[ii]; // the eoch time of the two VCMs are different
	  CONSTELLATION->spacecraft[ii][eee].et = OPTIONS->et_vcm[ii];
	  }
	  else{
	    CONSTELLATION->spacecraft[ii][eee].et_sc_initial = et;
	  CONSTELLATION->spacecraft[ii][eee].et = CONSTELLATION->et;
	  }
	  // COLLISION
	  CONSTELLATION->spacecraft[ii][eee].number_of_collisions = 0;
	  // end of COLLISION

	  /*****************************************************************************************************************************/
	  /********************************************************* ENSEMBLES ON COE **************************************************/
	  /*****************************************************************************************************************************/
	  if (( OPTIONS->nb_ensembles > 0 ) && ( strcmp(OPTIONS->type_orbit_initialisation, "oe" ) == 0 )){ // initialize COE if we run ensembles on the COE

	    // Initiate the COE to each ensemble as the same COE as for the spacecraft
	    CONSTELLATION->spacecraft[ii][eee].OE.sma            = ( PARAMS->EARTH.radius + OPTIONS->apogee_alt[ii] ) / ( 1 + OPTIONS->eccentricity[ii] ); // semi-major axis (CBV)

	    CONSTELLATION->spacecraft[ii][eee].OE.inclination    = OPTIONS->inclination[ii] * DEG2RAD;
	    CONSTELLATION->spacecraft[ii][eee].OE.w              = OPTIONS->w[ii] * DEG2RAD;
	    CONSTELLATION->spacecraft[ii][eee].OE.long_an        = OPTIONS->long_an[ii] * DEG2RAD;
	    CONSTELLATION->spacecraft[ii][eee].OE.f              = OPTIONS->f[ii] * DEG2RAD;
	    CONSTELLATION->spacecraft[ii][eee].OE.eccentricity   = OPTIONS->eccentricity[ii];
	    // For each COE to ensemble, replace it by a random normal number
	    // eccentricity
	    if ( OPTIONS->coe_to_ensemble[ii][0] == 1 ) {
	      CONSTELLATION->spacecraft[ii][eee].OE.eccentricity = randn( OPTIONS->eccentricity[ii], eccentricity_sigma );
	      // printf("eccentricity = %f\n", CONSTELLATION->spacecraft[ii][eee].OE.eccentricity);
	    }
	    // true anomaly
	    if ( OPTIONS->coe_to_ensemble[ii][1] == 1 ) {
	      CONSTELLATION->spacecraft[ii][eee].OE.f = randn( OPTIONS->f[ii] * DEG2RAD, f_sigma );
	      // printf("true ano = %f\n", CONSTELLATION->spacecraft[ii][eee].OE.f * RAD2DEG);
	    }
	    // RAAN
	    if ( OPTIONS->coe_to_ensemble[ii][2] == 1 ){
	      CONSTELLATION->spacecraft[ii][eee].OE.long_an = randn( OPTIONS->long_an[ii] * DEG2RAD, long_an_sigma );
	      // printf("RANN = %f\n", CONSTELLATION->spacecraft[ii][eee].OE.long_an * RAD2DEG);
	    }
	    // argument of perigee
	    if ( OPTIONS->coe_to_ensemble[ii][3] == 1 ) {
	      CONSTELLATION->spacecraft[ii][eee].OE.w = randn( OPTIONS->w[ii] * DEG2RAD, w_sigma );
	      // printf("argument of perigee = %f\n", CONSTELLATION->spacecraft[ii][eee].OE.w * RAD2DEG);
	    }
	    // inclination
	    if ( OPTIONS->coe_to_ensemble[ii][4] == 1 ) {
	      CONSTELLATION->spacecraft[ii][eee].OE.inclination = randn( OPTIONS->inclination[ii] * DEG2RAD, inclination_sigma );
	      //	  printf("inclination = %f\n",CONSTELLATION->spacecraft[ii][eee].OE.inclination * RAD2DEG);
	    }
	    // semi-major axis
	    if ( OPTIONS->coe_to_ensemble[ii][5] == 1 ) {
	      CONSTELLATION->spacecraft[ii][eee].OE.sma = randn( ( PARAMS->EARTH.radius + OPTIONS->apogee_alt[ii] ) / ( 1 + OPTIONS->eccentricity[ii] ), sma_sigma );
	      //	  printf("sma = %f\n", CONSTELLATION->spacecraft[ii][eee].OE.sma );
	    }
	  }

	  else if (( OPTIONS->nb_ensembles > 0 ) && ( strcmp(OPTIONS->type_orbit_initialisation, "state_eci" ) == 0 )){ // initialize state eci and OE if we run ensembles on state eci


	    CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[0] = randn( OPTIONS->x_eci[ii], OPTIONS->x_eci_sigma[ii] );
	    CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[1] = randn( OPTIONS->y_eci[ii], OPTIONS->y_eci_sigma[ii] );
	    CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[2] = randn( OPTIONS->z_eci[ii], OPTIONS->z_eci_sigma[ii] );

	    CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[0] = randn( OPTIONS->vx_eci[ii], OPTIONS->vx_eci_sigma[ii] );
	    CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[1] = randn( OPTIONS->vy_eci[ii], OPTIONS->vy_eci_sigma[ii] );
	    CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[2] = randn( OPTIONS->vz_eci[ii], OPTIONS->vz_eci_sigma[ii] );

	    // Orbital elements
	    cart2kep( &CONSTELLATION->spacecraft[ii][eee].OE, CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].et,  PARAMS->EARTH.GRAVITY.mu);

	    CONSTELLATION->spacecraft[ii][eee].OE.initial_an_to_sc = CONSTELLATION->spacecraft[ii][eee].OE.an_to_sc;
	   
	    CONSTELLATION->spacecraft[ii][eee].OE.w_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.w;
	    CONSTELLATION->spacecraft[ii][eee].OE.w_ave = 9999.999999/RAD2DEG;
	    CONSTELLATION->spacecraft[ii][eee].OE.sma_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.sma;
	    CONSTELLATION->spacecraft[ii][eee].OE.sma_ave = 9999.999999;
	    CONSTELLATION->spacecraft[ii][eee].OE.ecc_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.eccentricity;
	    CONSTELLATION->spacecraft[ii][eee].OE.ecc_ave = 9999.999999;

	    CONSTELLATION->spacecraft[ii][eee].OE.ave_increm = 1;
	    CONSTELLATION->spacecraft[ii][eee].et_last_orbit = CONSTELLATION->spacecraft[ii][eee].et;
	    CONSTELLATION->spacecraft[ii][eee].orbit_number = 0;


	    radius_perigee = CONSTELLATION->spacecraft[ii][eee].OE.sma * ( 1 - CONSTELLATION->spacecraft[ii][eee].OE.eccentricity  );
	    if (radius_perigee < PARAMS->EARTH.radius){
	      printf("***! The orbit of satellite %d intersects the Earth (altitude of perigee = %f km). The program will stop. !***\n",ii, radius_perigee-PARAMS->EARTH.radius);
	      MPI_Finalize();
	      exit(0);
	    }

	  }

	  else if (( OPTIONS->nb_ensembles > 0 ) && (( strcmp(OPTIONS->type_orbit_initialisation, "collision" ) == 0 ) || ( strcmp(OPTIONS->type_orbit_initialisation, "collision_vcm" ) == 0 ) )){ // initialize state eci and OE if we run ensembles on state eci from a collision input file
	    if ( strcmp(OPTIONS->type_orbit_initialisation, "collision_vcm" ) == 0 ){
	    // GOOD REFERENCE: http://www.prepacom.net/HEC2/math/cours/Changement%20de%20bases.pdf
	    // generate random uncertainties from teh eigenvalues. These uncertainties corresponds to the the basis of the principal axes of the ellipsoid

	    /* // !!!!!!! DELETE BLOCK BELOW AND UNCOMMENT THE ONE BELOW IT */
	    /* struct timeval t1; */
	    /* gettimeofday(&t1, NULL); */
	    /* seed = t1.tv_usec * t1.tv_sec;// time (NULL) + iProc; //\* getpid(); */
	    /* gsl_rng_set (r_gaussian_generator, seed);                  // set seed */
	    /* collision_equinoctial_sigma_in_diag_basis[0] = gsl_ran_gaussian(r_gaussian_generator, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][0]  ) ); */
	    /* collision_equinoctial_sigma_in_diag_basis[1] = gsl_ran_gaussian(r_gaussian_generator, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][1]  ) ); */
	    /* collision_equinoctial_sigma_in_diag_basis[2] = gsl_ran_gaussian(r_gaussian_generator, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][2]  ) ); */
	    /* collision_equinoctial_sigma_in_diag_basis[3] = gsl_ran_gaussian(r_gaussian_generator, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][3]  ) ); */
	    /* collision_equinoctial_sigma_in_diag_basis[4] = gsl_ran_gaussian(r_gaussian_generator, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][4]  ) ); */
	    /* collision_equinoctial_sigma_in_diag_basis[5] = gsl_ran_gaussian(r_gaussian_generator, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][5]  ) ); */

	    collision_equinoctial_sigma_in_diag_basis[0] = randn( 0, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][0]  ) );
	    collision_equinoctial_sigma_in_diag_basis[1] = randn( 0, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][1]  ) );
	    collision_equinoctial_sigma_in_diag_basis[2] = randn( 0, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][2]  ) );
	    collision_equinoctial_sigma_in_diag_basis[3] = randn( 0, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][3]  ) );
	    collision_equinoctial_sigma_in_diag_basis[4] = randn( 0, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][4]  ) );
	    collision_equinoctial_sigma_in_diag_basis[5] = randn( 0, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][5]  ) );
	    /* collision_equinoctial_sigma_in_diag_basis[6] = randn( 0, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][6]  ) ); */
	    /* collision_equinoctial_sigma_in_diag_basis[7] = randn( 0, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][7]  ) ); */
	    /* collision_equinoctial_sigma_in_diag_basis[8] = randn( 0, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][8]  ) ); */

	    
	    // convert thes e unertainties rfom the ellipsoid basis to ECI p
	    m_x_v6( collision_equinoctial_sigma_in_equinoctial_basis, OPTIONS->rotation_matrix_for_diagonalization[ii], collision_equinoctial_sigma_in_diag_basis ); /// !!!!! maybe invert of rotation_matrix_for_diagonalization?



/* 	    //	  printf("\n"); */
/* 	    for (ccc = 0; ccc < 6; ccc++){ */
/* 	      collision_equinoctial_sigma_in_equinoctial_basis[ccc] = collision_equinoctial_sigma_in_equinoctial_basis[ccc] / 1000.; // / 1000. to convert the eigenvalues from m to km (or from m/s to km/s) */
/* 	      //	    printf("%15.10e %15.10e %d %d %d\n", collision_equinoctial_sigma_in_equinoctial_basis[ccc], OPTIONS->eigenvalue_covariance_matrix[ii][ccc] / 1000. ,eee, ccc, ii); */
/* 	    } */


	    double af, ag, lequin, nequin, chi, psi;
	    double af_mean, ag_mean, lequin_mean, nequin_mean, chi_mean, psi_mean;
	    // calculate the mean equinotication elemnts (ie the eq elts of the mean psoiton and velcoity)
	    double mu = 398600.4418; // km^3/s^2 (ideally, should read from propagator.c -> load_params but the call to load_params is later)
	    double fr;
	    // fr is 1 for prograde orbits, -1 for retrograde orbits. !!!! this is a convention, but there are others that assume fr to be 1 all the time. make sure that the VCM was generated using the same convention
	    ORBITAL_ELEMENTS_T oe_temp;
	    double rvec[3], vvec[3];
	    rvec[0] = OPTIONS->x_eci[ii]; rvec[1] = OPTIONS->y_eci[ii]; rvec[2] = OPTIONS->z_eci[ii];
	    vvec[0] = OPTIONS->vx_eci[ii]; vvec[1] = OPTIONS->vy_eci[ii]; vvec[2] = OPTIONS->vz_eci[ii];

	    cart2kep(&oe_temp, rvec, vvec, OPTIONS->et_vcm[ii] , mu);


	    if (oe_temp.inclination <= M_PI/2.){
	      fr = 1;
	    }
	    else{
	      fr = -1;
	    }

	    cart_to_equin( &af_mean, &ag_mean, &lequin_mean, &nequin_mean, &chi_mean, &psi_mean,  mu,  fr, rvec, vvec); // r and v in km km/s

	    af = af_mean + collision_equinoctial_sigma_in_equinoctial_basis[0];
	    ag = ag_mean + collision_equinoctial_sigma_in_equinoctial_basis[1];
	    lequin = lequin_mean + collision_equinoctial_sigma_in_equinoctial_basis[2];
	    nequin = nequin_mean + collision_equinoctial_sigma_in_equinoctial_basis[3] * nequin_mean; // elemts of the covaraicen amtrix in the VCM were adimentionless, ie divied by their mean value. so need to multiply here
	    chi = chi_mean + collision_equinoctial_sigma_in_equinoctial_basis[4];
	    psi = psi_mean + collision_equinoctial_sigma_in_equinoctial_basis[5];

	    // bc and srp for the main sc will be calcucalted later
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.bc_vcm = OPTIONS->bc_cdm[ii] +  randn( 0, (OPTIONS->bc_cdm_std[ii]));
	      //OPTIONS->bc_vcm[ii] +  randn( 0, sqrt(OPTIONS->covariance_matrix_equinoctial[ii][6][6])) *  OPTIONS->bc_vcm[ii];  // elemts of the covaraicen amtrix in the VCM were adimentionless, ie divied by their mean value. so need to multiply here
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.srp_vcm = OPTIONS->srp_cdm[ii] +  randn( 0, (OPTIONS->srp_cdm_std[ii]));

	      //OPTIONS->srp_vcm[ii] + randn( 0, sqrt(OPTIONS->covariance_matrix_equinoctial[ii][8][8])) *  OPTIONS->srp_vcm[ii];  // elemts of the covaraicen amtrix in the VCM were adimentionless, ie divied by their mean value. so need to multiply here ([7] is for bdot)
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.bc_vcm = CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.bc_vcm  / 1000. / 1000.; // OPTIONS->bc_vcm in m2/kg. convert in km2/kg 
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.srp_vcm = CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.srp_vcm  / 1000. / 1000.; // OPTIONS->srp_vcm in m2/kg. convert in km2/kg 
	  

	    double rvec_pert[3], vvec_pert[3];
	    equin_to_cart(rvec_pert, vvec_pert, af, ag, lequin, nequin, chi, psi, mu, fr);

	    CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[0] = rvec_pert[0];
	    CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[1] = rvec_pert[1];
	    CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[2] = rvec_pert[2];

	    CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[0] = vvec_pert[0];
	    CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[1] = vvec_pert[1];
	    CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[2] = vvec_pert[2];

/* 	    v_print(CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, "r_pert"); */
/* 	    v_print(rvec, "r_mean"); */
	    double dr_pert[3], dr_pert_lvlh[3];
	    double dv_pert[3], dv_pert_lvlh[3];
	    v_sub(dr_pert, CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, rvec);
	    v_sub(dv_pert, CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL, vvec);
      double T_inrtl_2_lvlh_pert[3][3];
      compute_T_inrtl_2_lvlh(T_inrtl_2_lvlh_pert, rvec, vvec);
      m_x_v(dr_pert_lvlh, T_inrtl_2_lvlh_pert, dr_pert);
        m_x_v(dv_pert_lvlh, T_inrtl_2_lvlh_pert, dv_pert);
/*     v_print(dr_pert_lvlh, "dr_lvlh"); */
/*     v_print(dv_pert_lvlh, "dv_lvlh"); */
    r_alongtrack_all[ii][eee] = dr_pert_lvlh[0]; r_crosstrack_all[ii][eee] = dr_pert_lvlh[1]; r_radial_all[ii][eee] = dr_pert_lvlh[2];
        v_alongtrack_all[ii][eee] = dv_pert_lvlh[0]; v_crosstrack_all[ii][eee] = dv_pert_lvlh[1]; v_radial_all[ii][eee] = dv_pert_lvlh[2];
	bc_pert_all[ii][eee] = CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.bc_vcm;
	srp_pert_all[ii][eee] = CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.srp_vcm;
	//	bc_pert_all[ii][eee] = 
	    /* if ( ( ii == 0 )  && ( eee == start_ensemble + iProc * OPTIONS->nb_ensemble_min_per_proc ) ){  */
	    /*   printf("delta_vx_eci: %.10e - delta_vx_diag: %.10e - eee: %d - iProc: %d - vxsigma: %e - total vx_eci: %.10e\n", collision_equinoctial_sigma_in_equinoctial_basis[3], collision_equinoctial_sigma_in_diag_basis[3]/1000., eee, iProc, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][3] )/1000., CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[0] ); */
	    /*   //	    m_print6(OPTIONS->inverse_rotation_matrix_for_diagonalization[ii], "OPTIONS->inverse_rotation_matrix_for_diagonalization[ii]"); */
	    /*   } */

	    // This will be used to save position of sc in the time span of the TCA if collisions assessment is on
	    CONSTELLATION->spacecraft[ii][eee].ispan = 0;

	  
	  
	    // Orbital elements
	    cart2kep( &CONSTELLATION->spacecraft[ii][eee].OE, CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].et,  PARAMS->EARTH.GRAVITY.mu);
	    CONSTELLATION->spacecraft[ii][eee].OE.initial_an_to_sc = CONSTELLATION->spacecraft[ii][eee].OE.an_to_sc;
	    
	    CONSTELLATION->spacecraft[ii][eee].OE.w_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.w;
	    CONSTELLATION->spacecraft[ii][eee].OE.w_ave = 9999.999999/RAD2DEG;
	    CONSTELLATION->spacecraft[ii][eee].OE.sma_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.sma;
	    CONSTELLATION->spacecraft[ii][eee].OE.sma_ave = 9999.999999;
	    CONSTELLATION->spacecraft[ii][eee].OE.ecc_ave_temp = CONSTELLATION->spacecraft[ii][eee].OE.eccentricity;
	    CONSTELLATION->spacecraft[ii][eee].OE.ecc_ave = 9999.999999;

	    CONSTELLATION->spacecraft[ii][eee].OE.ave_increm = 1;
	    CONSTELLATION->spacecraft[ii][eee].et_last_orbit = CONSTELLATION->spacecraft[ii][eee].et;
	    CONSTELLATION->spacecraft[ii][eee].orbit_number = 0;


	    radius_perigee = CONSTELLATION->spacecraft[ii][eee].OE.sma * ( 1 - CONSTELLATION->spacecraft[ii][eee].OE.eccentricity  );
	    if (radius_perigee < PARAMS->EARTH.radius){
	      printf("***! The orbit of satellite %d intersects the Earth (altitude of perigee = %f km). The program will stop. !***\n",ii, radius_perigee-PARAMS->EARTH.radius);
	      MPI_Finalize();
	      exit(0);
	    }
	  }
	    else if ( strcmp(OPTIONS->type_orbit_initialisation, "collision" ) == 0 ){
	    // GOOD REFERENCE: http://www.prepacom.net/HEC2/math/cours/Changement%20de%20bases.pdf
	    // generate random uncertainties from teh eigenvalues. These uncertainties corresponds to the the basis of the principal axes of the ellipsoid

	    /* // !!!!!!! DELETE BLOCK BELOW AND UNCOMMENT THE ONE BELOW IT */
	    /* struct timeval t1; */
	    /* gettimeofday(&t1, NULL); */
	    /* seed = t1.tv_usec * t1.tv_sec;// time (NULL) + iProc; //\* getpid(); */
	    /* gsl_rng_set (r_gaussian_generator, seed);                  // set seed */
	    /* collision_eci_sigma_in_diag_basis[0] = gsl_ran_gaussian(r_gaussian_generator, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][0]  ) ); */
	    /* collision_eci_sigma_in_diag_basis[1] = gsl_ran_gaussian(r_gaussian_generator, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][1]  ) ); */
	    /* collision_eci_sigma_in_diag_basis[2] = gsl_ran_gaussian(r_gaussian_generator, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][2]  ) ); */
	    /* collision_eci_sigma_in_diag_basis[3] = gsl_ran_gaussian(r_gaussian_generator, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][3]  ) ); */
	    /* collision_eci_sigma_in_diag_basis[4] = gsl_ran_gaussian(r_gaussian_generator, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][4]  ) ); */
	    /* collision_eci_sigma_in_diag_basis[5] = gsl_ran_gaussian(r_gaussian_generator, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][5]  ) ); */

	    collision_eci_sigma_in_diag_basis[0] = randn( 0, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][0]  ) );
	    collision_eci_sigma_in_diag_basis[1] = randn( 0, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][1]  ) );
	    collision_eci_sigma_in_diag_basis[2] = randn( 0, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][2]  ) );
	    collision_eci_sigma_in_diag_basis[3] = randn( 0, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][3]  ) );
	    collision_eci_sigma_in_diag_basis[4] = randn( 0, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][4]  ) );
	    collision_eci_sigma_in_diag_basis[5] = randn( 0, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][5]  ) );



	    // convert thes e unertainties rfom the ellipsoid basis to ECI p
	    m_x_v6( collision_eci_sigma_in_eci_basis, OPTIONS->rotation_matrix_for_diagonalization[ii], collision_eci_sigma_in_diag_basis );



	    //	  printf("\n");
	    for (ccc = 0; ccc < 6; ccc++){
	      collision_eci_sigma_in_eci_basis[ccc] = collision_eci_sigma_in_eci_basis[ccc] / 1000.; // / 1000. to convert the eigenvalues from m to km (or from m/s to km/s)
	      //	    printf("%15.10e %15.10e %d %d %d\n", collision_eci_sigma_in_eci_basis[ccc], OPTIONS->eigenvalue_covariance_matrix[ii][ccc] / 1000. ,eee, ccc, ii);
	    }


	    CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[0] = OPTIONS->x_eci[ii] + collision_eci_sigma_in_eci_basis[0];
	    CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[1] = OPTIONS->y_eci[ii] + collision_eci_sigma_in_eci_basis[1];
	    CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[2] = OPTIONS->z_eci[ii] + collision_eci_sigma_in_eci_basis[2]; 

	    CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[0] = OPTIONS->vx_eci[ii] + collision_eci_sigma_in_eci_basis[3];
	    CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[1] = OPTIONS->vy_eci[ii] + collision_eci_sigma_in_eci_basis[4];
	    CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[2] = OPTIONS->vz_eci[ii] + collision_eci_sigma_in_eci_basis[5];

	    /* if ( ( ii == 0 )  && ( eee == start_ensemble + iProc * OPTIONS->nb_ensemble_min_per_proc ) ){  */
	    /*   printf("delta_vx_eci: %.10e - delta_vx_diag: %.10e - eee: %d - iProc: %d - vxsigma: %e - total vx_eci: %.10e\n", collision_eci_sigma_in_eci_basis[3], collision_eci_sigma_in_diag_basis[3]/1000., eee, iProc, sqrt( OPTIONS->eigenvalue_covariance_matrix[ii][3] )/1000., CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[0] ); */
	    /*   //	    m_print6(OPTIONS->inverse_rotation_matrix_for_diagonalization[ii], "OPTIONS->inverse_rotation_matrix_for_diagonalization[ii]"); */
	    /*   } */

	    // This will be used to save position of sc in the time span of the TCA if collisions assessment is on
	    CONSTELLATION->spacecraft[ii][eee].ispan = 0;

	  
	  
	    // Orbital elements
	    cart2kep( &CONSTELLATION->spacecraft[ii][eee].OE, CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].et,  PARAMS->EARTH.GRAVITY.mu);


	    radius_perigee = CONSTELLATION->spacecraft[ii][eee].OE.sma * ( 1 - CONSTELLATION->spacecraft[ii][eee].OE.eccentricity  );
	    if (radius_perigee < PARAMS->EARTH.radius){
	      printf("***! The orbit of satellite %d intersects the Earth (altitude of perigee = %f km). The program will stop. !***\n",ii, radius_perigee-PARAMS->EARTH.radius);
	      MPI_Finalize();
	      exit(0);
	    }

	    }
	    
	  } // end of initialize state eci and OE if we run ensembles on state eci from a collision input file

	  else { // initialize COE if we do not run ensembles on the orbital elements or on state eci (including collision case)

	    CONSTELLATION->spacecraft[ii][eee].OE.sma            = CONSTELLATION->spacecraft[ii][0].OE.sma; // semi-major axis (CBV)

	    //    printf("sma = %f\n", CONSTELLATION->spacecraft[ii][0].OE.sma);
	    radius_perigee = CONSTELLATION->spacecraft[ii][eee].OE.sma * ( 1 - CONSTELLATION->spacecraft[ii][0].OE.eccentricity );
	    if (radius_perigee < PARAMS->EARTH.radius){
	      printf("***! The orbit of satellite %d intersects the Earth (altitude of perigee = %f km). The program will stop. !***\n",ii, radius_perigee-PARAMS->EARTH.radius);
	      MPI_Finalize();
	      exit(0);
	    }
	    CONSTELLATION->spacecraft[ii][eee].OE.inclination    = CONSTELLATION->spacecraft[ii][0].OE.inclination ;

	    CONSTELLATION->spacecraft[ii][eee].OE.w              = CONSTELLATION->spacecraft[ii][0].OE.w ;                                // argument of perigee (CBV)
	    CONSTELLATION->spacecraft[ii][eee].OE.long_an        = CONSTELLATION->spacecraft[ii][0].OE.long_an ;                          // RAAN (CBV)
	    CONSTELLATION->spacecraft[ii][eee].OE.f              = CONSTELLATION->spacecraft[ii][0].OE.f ;                                // true anomaly (CBV)

	    CONSTELLATION->spacecraft[ii][eee].OE.eccentricity   = CONSTELLATION->spacecraft[ii][0].OE.eccentricity;

	  } // end of initialize COE if we do not run ensembles on the orbital elements or on state eci (including collision case)
	  if ( ( OPTIONS->nb_ensembles <= 0 ) || ( ( strcmp(OPTIONS->type_orbit_initialisation, "state_eci" ) != 0 ) && ( strcmp(OPTIONS->type_orbit_initialisation, "collision" ) != 0 )  && ( strcmp(OPTIONS->type_orbit_initialisation, "collision_vcm" ) != 0 ) ) ){ // intialize eci r/v if we do not run ensembles on state eci or on collision

	    // Initialize the inertial state
	    kep2cart(   CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL,
			CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL,
			&PARAMS->EARTH.GRAVITY.mu,
			&CONSTELLATION->spacecraft[ii][eee].OE); // Computes the ECI coordinates based on the Keplerian inputs (orbital elements and mu) (in propagator.c) (CBV)
	
	    // Right ascension
	    CONSTELLATION->spacecraft[ii][eee].OE.ra =  atan2(CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[1], CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[0]);
	    if ( CONSTELLATION->spacecraft[ii][eee].OE.ra < 0){
	      CONSTELLATION->spacecraft[ii][eee].OE.ra = 2*M_PI + CONSTELLATION->spacecraft[ii][eee].OE.ra ;
	    }

	    if (iDebugLevel >= 1){
	      if (iProc == 0) printf("-- (initialize_constellation) Done initializing et, r_i2cg_INRTL, v_i2cg_INRTL for ensemble spacecraft.\n");
	    }
	  }  // end of intialize eci r/v if we do not run ensembles on state eci or on collision
	} // end of the initialization of et, r_i2cg_INRTL, v_i2cg_INRTL, OE, number_of_collisions for the sc from an ensemble
      
	/*************************************************************************************/
	/*************************************************************************************/
	/************* COE AND TLE INITIALIZE for main SC and for ensemble SC:
- r_ecef2cg_ECEF, v_ecef2cg_ECEF
- geodetic
- filename, filenameecef, filenameout, filenamepower for main SC ONLY
- name_sat
- INTEGRATOR: dt, nb_surfaces, mass, solar_cell_efficiency, degree, order, include_drag/solar_pressure/earth_pressure/moon/sun, Ap, Ap_hist, f107, f107A, density, initialize_geo_with_bstar, sc_main_nb, sc_ensemble_nb
	**************/
	/*************************************************************************************/
	/*************************************************************************************/
	if (iDebugLevel >= 1){
	  if (iProc == 0) printf("-- (initialize_constellation) Initializing r_ecef2cg_ECEF, v_ecef2cg_ECEF, geodetic, output file names, name_sat, dt, nb_surfaces, mass, solar_cell_efficiency, degree, order, include_drag/solar_pressure/earth_pressure/moon/sun, Ap, Ap_hist, f107, f107A, density, initialize_geo_with_bstar for all spacecraft.\n");
	}
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.sc_main_nb = ii;
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.sc_ensemble_nb = eee;

	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.initialize_geo_with_bstar = OPTIONS->initialize_geo_with_bstar;

	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.coll_vcm = OPTIONS->coll_vcm;

	// !!!! REMOVE BLOCK BELOW
	//	OPTIONS->bc_vcm_std[ii] = 0.23 * OPTIONS->bc_vcm[ii];
	//OPTIONS->srp_vcm_std[ii] = 0;//0.23 * OPTIONS->srp_vcm[ii];
	// !!!! END OF REMOVE BLOCK BELOW
	if (eee == 0){ // bc_vcm and srp_vcm were calculated prevously (eigen value block). here only the main sc
	  ///	 CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.bc_vcm = OPTIONS->bc_vcm[ii] / 1000. / 1000.; // OPTIONS->bc_vcm in m2/kg. convert in km2/kg
	  //	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.srp_vcm = OPTIONS->srp_vcm[ii] / 1000. / 1000.; // OPTIONS->srp_vcm in m2/kg. convert in km2/kg
	 CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.bc_vcm = OPTIONS->bc_cdm[ii] / 1000. / 1000.; // OPTIONS->bc_vcm in m2/kg. convert in km2/kg
		  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.srp_vcm = OPTIONS->srp_cdm[ii] / 1000. / 1000.; // OPTIONS->srp_vcm in m2/kg. convert in km2/kg
		  CONSTELLATION->spacecraft[ii][eee].already_output_cd_ensemble = 0;
		  CONSTELLATION->spacecraft[ii][eee].already_output_srp_ensemble = 0;

	  //	  printf("MAIN %d %d %e (%e) %e (%e)\n", ii, eee, CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.bc_vcm * 1e6, OPTIONS->bc_vcm_std[ii], CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.srp_vcm  * 1e6, OPTIONS->srp_vcm_std[ii]);
	}
	
/* 	else{ */

/* 	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.bc_vcm = fabs(randn(OPTIONS->bc_vcm[ii], OPTIONS->bc_vcm_std[ii]));  */

/* 	//	  printf("bc sc[%d][%d] %e %e\n", ii, eee, CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.bc_vcm, OPTIONS->bc_vcm[ii]); */
/* 	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.bc_vcm = CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.bc_vcm  / 1000. / 1000.; // OPTIONS->bc_vcm in m2/kg. convert in km2/kg  */

/* 	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.srp_vcm = fabs(randn(OPTIONS->srp_vcm[ii], OPTIONS->srp_vcm_std[ii]));  */
/* 	  //printf("srp sc[%d][%d] %e %e\n", ii, eee, CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.srp_vcm, OPTIONS->srp_vcm[ii]); */
/* 	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.srp_vcm = CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.srp_vcm  / 1000. / 1000.; // OPTIONS->srp_vcm in m2/kg. convert in km2/kg  */
/* 	} */


	  // uncomment block below to ingore uncertainties in Bc or SRP
	 /* CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.bc_vcm = OPTIONS->bc_vcm[ii] / 1000. / 1000.; // OPTIONS->bc_vcm in m2/kg. convert in km2/kg  */
	
	 /*  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.srp_vcm = OPTIONS->srp_vcm[ii] / 1000. / 1000.; // OPTIONS->srp_vcm in m2/kg. convert in km2/kg  */
	 /*  // end of uncomment block below to ingore uncertainties in Bc or SRP */

	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.et_vcm = OPTIONS->et_vcm[ii];


/* 	// Initialize the geodetic state- */

/* 	if (iDebugLevel >= 2){ */
/* 	  if (iProc == 0) printf("--- (initialize_constellation) Initializing geodetic for all spacecraft.\n"); */
/* 	} */


/* 	eci2lla(CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL , CONSTELLATION->spacecraft[ii][eee].et, geodetic ); */

	
/* 	CONSTELLATION->spacecraft[ii][eee].GEODETIC.altitude = geodetic[2]; */
/* 	CONSTELLATION->spacecraft[ii][eee].GEODETIC.latitude = geodetic[0]; */
/* 	CONSTELLATION->spacecraft[ii][eee].GEODETIC.longitude = geodetic[1];  */
/* 	if (iDebugLevel >= 2){ */
/* 	  if (iProc == 0) printf("--- (initialize_constellation) Done initializing geodetic for all spacecraft.\n"); */
/* 	} */

	// Initialize the planet fixed state		
	if (iDebugLevel >= 2){
	  if (iProc == 0) printf("--- (initialize_constellation) Initializing r_ecef2cg_ECEF and v_ecef2cg_ECEF for all spacecraft.\n");
	}
	if ( strcmp( OPTIONS->type_orbit_initialisation, "state_ecef" ) != 0 ){ // already initialized if state_ecef chosen by user
/* 	  geodetic_to_geocentric(PARAMS->EARTH.flattening,             */
/* 				 CONSTELLATION->spacecraft[ii][eee].GEODETIC.altitude, */
/* 				 CONSTELLATION->spacecraft[ii][eee].GEODETIC.latitude, */
/* 				 CONSTELLATION->spacecraft[ii][eee].GEODETIC.longitude, */
/* 				 PARAMS->EARTH.radius,        */
/* 				 CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF) ; */
/*   v_print(CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, "ECI"); */
/*   v_print(CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF, "ECEF"); */
/*  printf("(%f, %f) # %f\n", CONSTELLATION->spacecraft[ii][eee].GEODETIC.latitude*RAD2DEG, CONSTELLATION->spacecraft[ii][eee].GEODETIC.longitude*RAD2DEG, CONSTELLATION->spacecraft[ii][eee].GEODETIC.altitude); */

	  estate[0] = CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[0];estate[1] = CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[1];estate[2] = CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[2];
	  estate[3] = CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[0];estate[4] = CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[1];estate[5] = CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[2];
	  sxform_c (  "J2000", PARAMS->EARTH.earth_fixed_frame,  CONSTELLATION->spacecraft[ii][eee].et,    xform  );
	  mxvg_c   (  xform,       estate,   6,  6, jstate );
	  CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF[0] = jstate[0]; CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF[1] = jstate[1]; CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF[2] = jstate[2];
	  CONSTELLATION->spacecraft[ii][eee].v_ecef2cg_ECEF[0] = jstate[3]; CONSTELLATION->spacecraft[ii][eee].v_ecef2cg_ECEF[1] = jstate[4]; CONSTELLATION->spacecraft[ii][eee].v_ecef2cg_ECEF[2] = jstate[5];
 /*  v_print(CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, "ECI"); */
 /*  v_print(CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF, "ECEF"); */
 /*  v_print(CONSTELLATION->spacecraft[ii][eee].v_ecef2cg_ECEF, "velocity ECEF"); */
 /* printf("(%f, %f) # %f\n", CONSTELLATION->spacecraft[ii][eee].GEODETIC.latitude*RAD2DEG, CONSTELLATION->spacecraft[ii][eee].GEODETIC.longitude*RAD2DEG, CONSTELLATION->spacecraft[ii][eee].GEODETIC.altitude); */


	}


/*  double T_J2000_to_ECEF[3][3]; */
/*  printf("CONSTELLATION->spacecraft[ii][eee].et = %f\n",CONSTELLATION->spacecraft[ii][eee].et); */
/* 	pxform_c( "J2000", PARAMS->EARTH.earth_fixed_frame, CONSTELLATION->spacecraft[ii][eee].et, T_J2000_to_ECEF); // Return the matrix (here T_J2000_to_ECEF) that transforms position vectors from one specified frame (here J2000) to another (here ITRF93) at a specified epoch  (in /Users/cbv/cspice/src/cspice/pxform_c.c) (CBV) */

/* 	  	m_x_v(CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF, T_J2000_to_ECEF, CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL); // multiply a matrix by a vector (in prop_math.c). So here we convert ECI (J2000) to ECEF coordinates (CBV) */
/* 	v_print(CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF,"uu"); */

	if (iDebugLevel >= 2){
	  if (iProc == 0) printf("--- (initialize_constellation) Done initializing r_ecef2cg_ECEF and v_ecef2cg_ECEF for all spacecraft.\n");
	}


	geocentric_to_geodetic(
			       CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF,
			       &PARAMS->EARTH.radius,
			       &PARAMS->EARTH.flattening,
			       &CONSTELLATION->spacecraft[ii][eee].GEODETIC.altitude,
			       &CONSTELLATION->spacecraft[ii][eee].GEODETIC.latitude,
			       &CONSTELLATION->spacecraft[ii][eee].GEODETIC.longitude ); // Computes lat/long/altitude based on ECEF coordinates and planet fixed state (planetary semimajor axis, flattening parameter) (in propagator.c) (CBV)



	// !!!!!!!!!!!!!!! ERASE
	/*       double x[6]; */
	/* double lt; */
	/*     spkez_c(10, CONSTELLATION->spacecraft[ii][eee].et, "J2000", "NONE", 399, x, &lt); //   Return the state (position and velocity) of a target bodyrelative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration. */
	/*     double sun_norm[3]; */
	/*     double sun[3]; */
	/*     sun[0] = x[0]; sun[1] = x[1]; sun[2] = x[2]; */
	/*     v_norm(sun_norm, sun); */
	/*     v_scale(sun_norm,sun_norm,6.878137e+03); */
	/*     v_print(sun_norm,"sun_norm"); */
	/*     v_print(CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL,"CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL"); */
	/*     v_print(CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF,"CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF"); */
	/*     printf("lat = %f| lon = %f| alt = %f\n", CONSTELLATION->spacecraft[ii][eee].GEODETIC.latitude*RAD2DEG, CONSTELLATION->spacecraft[ii][eee].GEODETIC.longitude*RAD2DEG,CONSTELLATION->spacecraft[ii][eee].GEODETIC.altitude); */
	/*     v_print(sun,"sun"); */
	/*     exit(0); */
	// !!!!!!!!!!!!!!! END OF ERASE
 
	// !!!!!!!!!!!!!!!!!!!!!!!!!! THE BLOCK BELOW IS TO INITIALIZE SATELLLITE 2 WITH THE CORRECT SPACING WITH RESPECT TO SATELLITE 1. IT IS JUST FOR A TRY, AND SHOULD BE REMOVED AND IMPLEMENTED IN THE CODE ITSELF. SO IF IT'S STILL HERE AFTER JUNE 10TH, 2016 THEN REMOVE IT!
	/* 	if ( strcmp( OPTIONS->type_orbit_initialisation, "oe" ) == 0 ){ */
	/* if (ii > 0){ */
	/* 	    v_copy( CONSTELLATION->spacecraft[ii][eee].r_ecef2cg_ECEF, CONSTELLATION->spacecraft[0][eee].r_ecef2cg_ECEF); */
	/* 	    v_copy( CONSTELLATION->spacecraft[ii][eee].v_ecef2cg_ECEF, CONSTELLATION->spacecraft[0][eee].v_ecef2cg_ECEF); */
	/* 	    CONSTELLATION->spacecraft[ii][eee].GEODETIC.altitude = CONSTELLATION->spacecraft[0][eee].GEODETIC.altitude ; */
	/* 	    CONSTELLATION->spacecraft[ii][eee].GEODETIC.latitude = CONSTELLATION->spacecraft[0][eee].GEODETIC.latitude ; */
	/* 	    CONSTELLATION->spacecraft[ii][eee].GEODETIC.longitude = CONSTELLATION->spacecraft[0][eee].GEODETIC.longitude ; */
	/* } */
	/* 	} */
	// !!!!!!!!!!!!!!!!!!!!!!!!!! END OF THE BLOCK BELOW IS TO INITIALIZE SATELLLITE 2 WITH THE CORRECT SPACING WITH RESPECT TO SATELLITE 1. IT IS JUST FOR A TRY, AND SHOULD BE REMOVED AND IMPLEMENTED IN THE CODE ITSELF. SO IF IT'S STILL HERE AFTER JUNE 10TH, 2016 THEN REMOVE IT!
	// Output file for each spacecraft
	if (iDebugLevel >= 2){
	  if (iProc == 0) printf("--- (initialize_constellation) Initializing the output file names and name_sat.\n");
	}

	if (eee == 0){ // filenames are only for the main sc 
	  // COLLISION
	  strcpy(CONSTELLATION->filename_collision, OPTIONS->dir_output_run_name);
	  strcat( CONSTELLATION->filename_collision, "/");
	  strcat( CONSTELLATION->filename_collision,OPTIONS->filename_collision);
	  // end of COLLISION

	  strcpy(CONSTELLATION->spacecraft[ii][0].filename, OPTIONS->dir_output_run_name_sat_name[ii]);
	  strcat(CONSTELLATION->spacecraft[ii][0].filename, "/");
	  strcat(CONSTELLATION->spacecraft[ii][0].filename, OPTIONS->filename_output[ii]);

	  if (( strcmp( OPTIONS->type_orbit_initialisation, "tle" ) == 0) || (strcmp(OPTIONS->type_orbit_initialisation, "tle_sgp4" ) == 0 )){
	    strcpy(CONSTELLATION->spacecraft[ii][0].filenametle, OPTIONS->dir_output_run_name_sat_name[ii]);
	    strcat(CONSTELLATION->spacecraft[ii][0].filenametle, "/");
	    strcat(CONSTELLATION->spacecraft[ii][0].filenametle, "TLE_");
	    strcat(CONSTELLATION->spacecraft[ii][0].filenametle, OPTIONS->filename_output[ii]);
	  }

	  // Now CONSTELLATION->spacecraft[ii][0].filenameecef is moved to generate_ephemerides because we want iProc 0 to know CONSTELLATION->spacecraft[ii][0].filenameecef even for the main sc ii that it does not run (this is because iProc 0 will gather all ECEF files at the end of the propagation)
	  /* 	strcpy(CONSTELLATION->spacecraft[ii][0].filenameecef, OPTIONS->dir_output_run_name_sat_name[ii]); */
	  /* 	strcat(CONSTELLATION->spacecraft[ii][0].filenameecef, "/"); */
	  /* 	strcat(CONSTELLATION->spacecraft[ii][0].filenameecef, "ECEF_"); */
	  /* 	strcat(CONSTELLATION->spacecraft[ii][0].filenameecef, OPTIONS->filename_output[ii]); */

	  strcpy(CONSTELLATION->spacecraft[ii][0].filenameout, OPTIONS->dir_output_run_name_sat_name[ii]);
	  strcat(CONSTELLATION->spacecraft[ii][0].filenameout, "/");
	  strcat(CONSTELLATION->spacecraft[ii][0].filenameout, "LLA_");
	  strcat(CONSTELLATION->spacecraft[ii][0].filenameout,OPTIONS->filename_output[ii]);


	  strcpy(CONSTELLATION->spacecraft[ii][0].filenamerho, OPTIONS->dir_output_run_name_sat_name[ii]);
	  strcat(CONSTELLATION->spacecraft[ii][0].filenamerho, "/");
	  strcat(CONSTELLATION->spacecraft[ii][0].filenamerho, "density_");
	  strcat(CONSTELLATION->spacecraft[ii][0].filenamerho,OPTIONS->filename_output[ii]);

	  strcpy(CONSTELLATION->spacecraft[ii][0].filenameatt, OPTIONS->dir_output_run_name_sat_name[ii]);
	  strcat(CONSTELLATION->spacecraft[ii][0].filenameatt, "/");
	  strcat(CONSTELLATION->spacecraft[ii][0].filenameatt, "attitude_");
	  strcat(CONSTELLATION->spacecraft[ii][0].filenameatt,OPTIONS->filename_output[ii]);


	  strcpy(CONSTELLATION->spacecraft[ii][0].filenamekalman, OPTIONS->dir_output_run_name_sat_name[ii]);
	  strcat(CONSTELLATION->spacecraft[ii][0].filenamekalman, "/");
	  strcat(CONSTELLATION->spacecraft[ii][0].filenamekalman, "kalman_");
	  strcat(CONSTELLATION->spacecraft[ii][0].filenamekalman,OPTIONS->filename_output[ii]);


	  strcpy(CONSTELLATION->spacecraft[ii][0].filenamekalman_meas, OPTIONS->dir_output_run_name_sat_name[ii]);
	  strcat(CONSTELLATION->spacecraft[ii][0].filenamekalman_meas, "/");
	  strcat(CONSTELLATION->spacecraft[ii][0].filenamekalman_meas, "meas_converted_kalman_");
	  strcat(CONSTELLATION->spacecraft[ii][0].filenamekalman_meas,OPTIONS->filename_output[ii]);

	  strcpy(CONSTELLATION->spacecraft[ii][0].filename_kalman_init, OPTIONS->filename_kalman_init);

	  strcpy(CONSTELLATION->spacecraft[ii][0].INTEGRATOR.filename_given_output, OPTIONS->dir_output_run_name_sat_name[ii]);
	  strcat(CONSTELLATION->spacecraft[ii][0].INTEGRATOR.filename_given_output, "/");
	  strcat(CONSTELLATION->spacecraft[ii][0].INTEGRATOR.filename_given_output, "given_output_");
	  strcat(CONSTELLATION->spacecraft[ii][0].INTEGRATOR.filename_given_output,OPTIONS->filename_output[ii]);


	  if (OPTIONS->nb_ground_stations > 0){
	    for ( iground = 0; iground < OPTIONS->nb_ground_stations; iground ++){
	      strcpy(CONSTELLATION->spacecraft[ii][0].filename_coverage_ground_station[iground], OPTIONS->dir_output_run_name_sat_name_coverage[ii]);
	      strcat(CONSTELLATION->spacecraft[ii][0].filename_coverage_ground_station[iground], "/");
	      strcat(CONSTELLATION->spacecraft[ii][0].filename_coverage_ground_station[iground], OPTIONS->name_ground_station[iground]);
	      strcat(CONSTELLATION->spacecraft[ii][0].filename_coverage_ground_station[iground], "_by_");
	      strcat(CONSTELLATION->spacecraft[ii][0].filename_coverage_ground_station[iground],OPTIONS->filename_output[ii]);
	    }
	  }


	  if (OPTIONS->solar_cell_efficiency != -1){
	    strcpy(CONSTELLATION->spacecraft[ii][0].filenamepower, OPTIONS->dir_output_run_name_sat_name[ii]);
	    strcat(CONSTELLATION->spacecraft[ii][0].filenamepower, "/");
	    strcat(CONSTELLATION->spacecraft[ii][0].filenamepower, "power_");
	    strcat(CONSTELLATION->spacecraft[ii][0].filenamepower,OPTIONS->filename_output[ii]);

	    strcpy(CONSTELLATION->spacecraft[ii][0].filenameeclipse, OPTIONS->dir_output_run_name_sat_name[ii]);
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameeclipse, "/");
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameeclipse, "eclipse_");
	    strcat(CONSTELLATION->spacecraft[ii][0].filenameeclipse,OPTIONS->filename_output[ii]);

	  }


	} // end of eee = 0

	//	    if (eee == start_ensemble + iProc * OPTIONS->nb_ensemble_min_per_proc){ // only receives it one time 
	if (eee == array_sc[start_ensemble[ii]]){ // if eee is the first sc run by iProc for this given main sc ii
	  if (write_attitude_file == 1){
	    strcpy(filename_attitude[ii], OPTIONS->dir_output_run_name_sat_name[ii]);
	    strcat(filename_attitude[ii], "/");
	    strcat(filename_attitude[ii], "attitude_");
	    strcat(filename_attitude[ii],OPTIONS->filename_output[ii]);
	    file_attitude[ii] = NULL;
	    file_attitude[ii] = fopen(filename_attitude[ii], "w+");
	    if (file_attitude[ii] == NULL){
	      printf("***! Could not open the output file for the attitude called %s. The program will stop. !***\n", filename_attitude[ii]); MPI_Finalize(); exit(0);
	    }
	    fprintf(file_attitude[ii], "This file shows the attitude of each ensemble spacecraft of the main spacecraft %s. One block per ensemble spacecraft. First line is the pitch, roll, and yaw angular velocities, as well as the order of rotation (pitch, roll, yaw). Recall that this order is the same for all ensemble spacecraft. Following lines show time vs pitch, roll, yaw (rotation of the body reference system with respect to the LVLH reference system). All angles are given in degrees.\n", OPTIONS->name_sat[ii]);
	  
	  }
	}

	// Name of each satellites
	strcpy(CONSTELLATION->spacecraft[ii][eee].name_sat, "");
	strcpy(CONSTELLATION->spacecraft[ii][eee].name_sat, OPTIONS->name_sat[ii]);
	if (iDebugLevel >= 2){
	  if (iProc == 0) printf("--- (initialize_constellation) Done initializing the output file names and name_sat.\n");
	}

        // Compute the Orbital period of the first state
	if ( ( ii  == 0 ) && ( eee == 0 ) ){ // !!!!!!!!!! FOR NOW WORKS ONLY IF TWO REFERENCE SATELLIES ONLY 
	  period = pow( CONSTELLATION->spacecraft[ii][0].OE.sma, 3.0);
	  period = period / PARAMS->EARTH.GRAVITY.mu;
	  period = 2.0 * M_PI * sqrt( period );
	  CONSTELLATION->collision_time_span =  period / 10. ;    // the time span is half of the orbit of the first sc (considered to be the primary sc)
	  //            ptd(CONSTELLATION->collision_time_span , "old");
	  // Actually, we re-evaluate the time span so it is an EVEN multiple of OPTIONS->dt (the closest possible to period/2)
	  if ( fmod(CONSTELLATION->collision_time_span, OPTIONS->dt) == 0 ){
	    if ( fmod( CONSTELLATION->collision_time_span / OPTIONS->dt, 2 )  == 1  ){ // if time span is an odd multiple of dt then remove dt from it to make it even
	      CONSTELLATION->collision_time_span = CONSTELLATION->collision_time_span - OPTIONS->dt;
	    }
	  }
	  else{
	    m_time_span = (int)( CONSTELLATION->collision_time_span / OPTIONS->dt );
	    if ( fmod( m_time_span, 2 ) == 1 ){
	      m_time_span = m_time_span + 1;
	      CONSTELLATION->collision_time_span = m_time_span * OPTIONS->dt;
	    }
	    else{
	      CONSTELLATION->collision_time_span = m_time_span * OPTIONS->dt;
	    }  
	  }
	  CONSTELLATION->collision_time_span = CONSTELLATION->collision_time_span + 2 * OPTIONS->dt; // for the reason why we add 2 * OPTIONS->dt, see comment "ABOUT THE TIME SPAN" at the end of generate_ephemerides
	  /* ptd(CONSTELLATION->collision_time_span , "new"); */
	  /* exitall(); */
	  // end of actually, we re-evaluate the time span so it is an EVEN multiple of OPTIONS->dt (the closest possible to period/2)
	  if (nProcs > 1){
	    for (ccc = 1; ccc < nProcs; ccc++){ // main sc 0 is for sure run by iProc 0, not matter how many main sc/ensemble sc/iProc there are
	      MPI_Send(&CONSTELLATION->collision_time_span, 1, MPI_DOUBLE, ccc, 0, MPI_COMM_WORLD);
	    }
	  }
	}
      	if (nProcs > 1){
	  if ( iProc > 0){
	    if ((ii == 0) && (eee == array_sc[start_ensemble[ii]])){ // only receives it one time  (array_sc[start_ensemble[ii]] is the first sc run by this iProc. Since main sc 0 is never run by iProc > 0, array_sc[start_ensemble[ii]] corresponds to the first ensemble sc run by this iProc)
	      MPI_Recv(&CONSTELLATION->collision_time_span, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	  }
	}

	if (iDebugLevel >= 2){
	  if (iProc == 0) printf("--- (initialize_constellation) Initializing dt, nb_surfaces, mass, solar_cell_efficiency, degree, order, include_drag/solar_pressure/earth_pressure/moon/sun, Ap, Ap_hist, f107, f107A, density, initialize_geo_with_bstar for all spacecraft.\n");
	}

	// Integrator
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.dt            = OPTIONS->dt;
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.dt_pos_neg            = OPTIONS->dt_pos_neg;
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.nb_surfaces   = OPTIONS->n_surfaces;   // Number of surfaces on the SC
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.nb_surfaces_eff   = OPTIONS->n_surfaces_eff;   // Number of surfaces on the SC
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.mass          = OPTIONS->mass;         // Mass of spacecraft
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.solar_cell_efficiency = OPTIONS->solar_cell_efficiency; // Solar cell efficiency
	strcpy( CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.opengl_filename_solar_power, OPTIONS->opengl_filename_solar_power ); // Solar cell efficiency
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.opengl = OPTIONS->opengl;
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.opengl_power = OPTIONS->opengl_power;
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.degree        = OPTIONS->degree;       // Gravity degree
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.order         = OPTIONS->order;        // Gravity order
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_drag  = OPTIONS->include_drag; // include drag (CBV)
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.thrust  = OPTIONS->thrust;
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_solar_pressure  = OPTIONS->include_solar_pressure; // include drag (CBV)
		CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_earth_pressure  = OPTIONS->include_earth_pressure; // include drag (CBV)
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_sun   = OPTIONS->include_sun;  // include Sun perturbations (CBV)
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_moon  = OPTIONS->include_moon; // include Moon perturbations (CBV)
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.density_mod = OPTIONS->density_mod; // // the desnity given by msis is multiplied by density_mod + density_amp * sin(2*pi*t/T + density_phase*T) where T is the orbital period  
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.density_mod_amp = OPTIONS->density_mod_amp; 
	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.density_mod_phase = OPTIONS->density_mod_phase; 


	if ( CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_drag == 1 ){ // if the user wants to use drag
	  strcpy(CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.format_density_driver, OPTIONS->format_density_driver);

	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.density = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );

	  if (  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.density == NULL ){
	    printf("***! Could not allow memory space for  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.density \n. The program will stop. !***\n");
	    MPI_Finalize();
	    exit(0);
	  }


	  if ( ( strcmp(OPTIONS->format_density_driver, "density_file") != 0 ) && ( strcmp(OPTIONS->format_density_driver, "gitm") != 0 ) ){ // the user chooses f107 and Ap for the density
	    if ( strcmp(OPTIONS->format_density_driver, "static") == 0 ){ // if the user chooses a constant f107 and Ap for the density  

	      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.Ap_static            = OPTIONS->Ap_static;           // magnetic index(daily)
	      for (hhh = 0; hhh < 7; hhh++){
		CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.Ap_hist_static[hhh]            = OPTIONS->Ap_hist_static[hhh];           // magnetic index(historical)
	      }
	      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.f107_static          = OPTIONS->f107_static;         // Daily average of F10.7 flux
	      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.f107A_static         = OPTIONS->f107A_static;        // 81 day average of F10.7 flux

	    }
	    else{  // if the user chooses a time-varying f107 and Ap for the density  


	      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.Ap = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) ); // "* 2.0" because of the Runge Kunta order 4 method
	      if (OPTIONS->use_ap_hist == 1){
		CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.Ap_hist = malloc(  7 * sizeof(double *) ); // historical ap

		for (hhh = 0; hhh < 7; hhh++){
		  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.Ap_hist[hhh] = malloc(OPTIONS->nb_time_steps * 2 * sizeof(double ));
		}
	      }

	      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.f107 = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
	      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.f107A = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );

	      if (  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.Ap == NULL ){
		printf("***! Could not allow memory space for  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.Ap \n. The program will stop. !***\n");
		MPI_Finalize();
		exit(0);
	      }
	      if (OPTIONS->use_ap_hist == 1){
		if (  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.Ap_hist == NULL ){
		  printf("***! Could not allow memory space for  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.Ap_hist \n. The program will stop. !***\n");
		  MPI_Finalize();
		  exit(0);
		}
	      }
	      if (  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.f107 == NULL ){
		printf("***! Could not allow memory space for  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.f107 \n. The program will stop. !***\n");
		MPI_Finalize();
		exit(0);
	      }
	      if (  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.f107A == NULL ){
		printf("***! Could not allow memory space for  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.f107A \n. The program will stop. !***\n");
		MPI_Finalize();
		exit(0);
	      }
	      /* for (aaa = 0; aaa< OPTIONS->nb_time_steps * 2; aaa++){  // "* 2.0" because of the Runge Kunta order 4 method */
	      /*   if ( OPTIONS->et_interpo[aaa] <= et_final_epoch + 0.01 ) {  */
	      /* 	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.Ap[aaa]            = OPTIONS->Ap[aaa];           // magnetic index(daily) */
	      /* 	if (OPTIONS->use_ap_hist == 1){ */
	      /* 	  for (hhh = 0; hhh < 7; hhh++){ */
	      /* 	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.Ap_hist[hhh][aaa]            = OPTIONS->Ap_hist[hhh][aaa];           // magnetic index(historical) */
	      /* 	  } */
	      /* 	} */

	      /* 	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.f107[aaa]          = OPTIONS->f107[aaa];         // Daily average of F10.7 flux */

	      /*   } */
	      /* } */
	      if ((strcmp(OPTIONS->test_omniweb_or_external_file, "swpc_mod") == 0) && (OPTIONS->swpc_need_predictions)){
		CONSTELLATION->sum_sigma_for_f107_average[ii][eee] = 0;
	      }

	      if (OPTIONS->nb_ensembles_density) { // the user chose to run ensembles on the density using data from SWPC 
		if (OPTIONS->swpc_need_predictions){ // if future predictions of F10.7 and Ap (not only past observations)
		  CONSTELLATION->sum_sigma_for_f107_average[ii][eee] = 0;

		  if (eee == 1 + iProc * OPTIONS->nb_ensemble_min_per_proc){ // allocate memory for arrays only once per iProc  
		    CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time[iProc] = malloc( OPTIONS->nb_ensemble_min_per_proc * sizeof( double ) ); 
		    if (CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time[iProc] == NULL){
		      printf("***! (propagate_spacecraft)(compute_drag) (generate_ensemble_f107_ap) There is not enough memory for CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time[iProc]. The program will stop. !***\n"); MPI_Finalize();exit(0);
		    }
		    CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time_sorted[iProc] = malloc( OPTIONS->nb_ensemble_min_per_proc *  sizeof( double ) ); 
		    if (CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time_sorted[iProc] == NULL){
		      printf("***! (propagate_spacecraft)(compute_drag) (generate_ensemble_f107_ap) There is not enough memory for CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time_sorted[iProc]. The program will stop. !***\n"); MPI_Finalize();exit(0);
		    }
		    CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time[iProc] = malloc( OPTIONS->nb_ensemble_min_per_proc * sizeof( double ) ); 
		    if (CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time[iProc] == NULL){
		      printf("***! (propagate_spacecraft)(compute_drag) (generate_ensemble_ap_ap) There is not enough memory for CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time[iProc]. The program will stop. !***\n"); MPI_Finalize();exit(0);
		    }
		    CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time_sorted[iProc] = malloc( OPTIONS->nb_ensemble_min_per_proc *  sizeof( double ) ); 
		    if (CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time_sorted[iProc] == NULL){
		      printf("***! (propagate_spacecraft)(compute_drag) (generate_ensemble_ap_ap) There is not enough memory for CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time_sorted[iProc]. The program will stop. !***\n"); MPI_Finalize();exit(0);
		    }
		  } // end of allocate memory for arrays only once per iProc  
		} // end of if future predictions of F10.7 and Ap (not only past observations)
	      } // end of the user chose to run ensembles on the density using data from SWPC

	      /* if ( (OPTIONS->nb_ensembles_density) && (OPTIONS->swpc_need_predictions) && ( et_initial_epoch >= OPTIONS->swpc_et_first_prediction ) ) {  */
	      /* 	    // if the user chose to run ensembles on the density using data from SWPC  */
	      /* 	  	    // if future predictions of F10.7 and Ap (not only past observations) (past values (value before the current time at running) are perfectly known (result from observations, not predictions)) */
	      /* 	    // only overwrite predictions of ensembles (not observations). OPTIONS->et_interpo[aaa] corresponds to the first say of predictions */
	      /*   start_ensemble_bis = 1; */
	      /*   nb_index_in_81_days =  81. * 24 * 3600 / (OPTIONS->dt/2.) + 1 ; */
	      /* 	    if ( eee == start_ensemble_bis + iProc * OPTIONS->nb_ensemble_min_per_proc){ // initialize array only once per iProc */
	      /* 	    // // Generate nb_ensembles_density normal random values */
	      /* 	    for ( eee_bis = 0; eee_bis < OPTIONS->nb_ensemble_min_per_proc; eee_bis++){ */
	      /* 	      CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time[iProc][eee_bis]   = randn( OPTIONS->f107[0], OPTIONS->sigma_f107[CONSTELLATION->aaa_sigma]); */
	      /* 	      CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time[iProc][eee_bis]   = randn(  OPTIONS->Ap[0], OPTIONS->sigma_ap[CONSTELLATION->aaa_sigma]); */

	      /* 	      /\* if (eee_bis == 0){ *\/ */
	      /* 	      /\* etprint(OPTIONS->et_interpo[aaa], "time"); *\/ */
	      /* 	      /\* printf("Ap[%d]: %f | sigma_ap[%d]: %f \n",aaa, OPTIONS->Ap[aaa],CONSTELLATION->aaa_sigma, OPTIONS->sigma_ap[CONSTELLATION->aaa_sigma]); *\/ */
	      /* 	      /\* } *\/ */
	      /* 	    } */

	      /* 	    // // Order values in ascending order */
	      /* 	    sort_asc_order(CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time_sorted[iProc],  CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time[iProc], OPTIONS->nb_ensemble_min_per_proc); */
	      /* 	    sort_asc_order(CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time_sorted[iProc],  CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time[iProc], OPTIONS->nb_ensemble_min_per_proc); */


	

	      /* 	      // Initialization of swpc_et_first_prediction for the calculation of F10.7A for an ensemble */
	      /* 		  if ( et_initial_epoch > OPTIONS->swpc_et_first_prediction){ // if the propagation starts in the future. Otherwise, sum_sigma_for_f107_average = 0 for all ensemble sc */

	      /* 		    if (  CONSTELLATION->sum_sigma_for_f107_average[ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_main_nb][ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb] == 0 )  {// for the first time that variable et gets bigger than swpc_et_first_prediction (0.01 for numerical reasons)		     */
	      /* 		   // initialize sum_sigma_for_f107_average as the sum of all sigma on F10.7 from first prediction to inital epoch. There is no mathematical logic in here. The exact solution would be to sum all the deviations between f107_ensemble and f107_refce from first prediciton until intial epoch, and sum them up. This is KIND OF similar. The reason I don't do the correct approach is that I did not calculate f107_ensemble for times before initial epoch */
	      /* 		      for (aaa_a = 0; aaa_a < CONSTELLATION->aaa_sigma; aaa_a++){ */

	      /* 		      CONSTELLATION->sum_sigma_for_f107_average[ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_main_nb][ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb] = CONSTELLATION->sum_sigma_for_f107_average[ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_main_nb][ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb] + OPTIONS->sigma_f107[aaa_a];		     */
	      /* 		      //		      ptd( OPTIONS->sigma_f107[aaa_a], "s"); */
	      /* 		      } */

	      /* 		      //		      printf("eee_bis: %d | index: %d | iProc: %d | sum: %f\n",  CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb, index_in_driver_interpolated, iProc, CONSTELLATION->sum_sigma_for_f107_average[ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_main_nb][ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb]);exit(0); */
		   
	      /* 		    } // end of for the first time that variable et gets bigger than swpc_et_first_prediction */
	      /* 		  }// end of if the propagation starts in the future   */
	      /* 	      // End of initialization of swpc_et_first_prediction for the calculation of F10.7A for an ensemble */

	      /* 	  } // initialize array only once per iProc */
	      /* 	    else if ( CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb == 0){ // also need to calculate f107 and ap for reference sc (no perturbation so directly equal to the values in the prediction files) */
	      /* 	     CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.Ap[0]            = OPTIONS->Ap[0];           // magnetic index(daily) */
	      /* 	    if (OPTIONS->use_ap_hist == 1){ */
	      /* 	      for (hhh = 0; hhh < 7; hhh++){ */
	      /* 		 CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.Ap_hist[hhh][0]            = OPTIONS->Ap_hist[hhh][0];           // magnetic index(historical) */
	      /* 	      } */
	      /* 	    } */

	      /* 	     CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.f107[0]          = OPTIONS->f107[0];         // Daily average of F10.7 flux */
	      /* 	      /\* print_test(); *\/ */
	      /* 	      /\* printf(" CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.f107[%d]: %f\n", 0,  CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.f107[0]); *\/ */

	      /* 	     CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.f107A[0]         = OPTIONS->f107A[0];        // 81 day average of F10.7 flux */


	      /* 	    } // end of also need to calculate f107 and ap for reference sc (no perturbation so directly equal to the values in the prediction files) */
	      /* 	    if ( CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb != 0) { // don't overwrite previously written values of f107 and Ap for reference sc */
	      /* 	     CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.Ap[0]            = CONSTELLATION->ensemble_array_per_iproc_ap_at_given_time_sorted[iProc][ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb-1-iProc*OPTIONS->nb_ensemble_min_per_proc]; */
	      /* 	     CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.f107[0]          = CONSTELLATION->ensemble_array_per_iproc_f107_at_given_time_sorted[iProc][ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb-1-iProc*OPTIONS->nb_ensemble_min_per_proc];  */
	
	      /* 	    CONSTELLATION->sum_sigma_for_f107_average[ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_main_nb][ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb] = CONSTELLATION->sum_sigma_for_f107_average[ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_main_nb][ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb] +  (  CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.f107[0]  - OPTIONS->f107[0] ); */
	      /* 	     CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.f107A[0]         = OPTIONS->f107A[0] + CONSTELLATION->sum_sigma_for_f107_average[ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_main_nb][ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb] / nb_index_in_81_days; // derivation of F10.7A considering uncerainties in F10.7 */
	      /* 	    //	    printf("eee_bis: %d | index: %d | iProc: %d | sum: %f | 81: %f | opeion %f\n",  CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb, index_in_driver_interpolated, iProc, CONSTELLATION->sum_sigma_for_f107_average[ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_main_nb][ CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb], nb_index_in_81_days, OPTIONS->dt); */
	      /* 	    //	    print_test(); */
	      /* 	  /\* if (iProc == 0){ *\/ */
	      /* 	  /\*   //	    if (( CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb == start_ensemble_bis + iProc * OPTIONS->nb_ensemble_min_per_proc + 1)){ *\/ */
	      /* 	  /\*   if (index_in_driver_interpolated == 4){ *\/ */
	      /* 	  /\*     printf("eee_bis: %d | f107:  %f | index: %d | iProc: %d\n",  CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb,  CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.f107[index_in_driver_interpolated], index_in_driver_interpolated, iProc); *\/ */
	      /* 	  /\* } *\/ */
	      /* 	  /\* } *\/ */
	      /* 	    } // end of don't overwrite previously written values of f107 and Ap for reference sc */
	      /* 	  } // end of: */
	      /* 	    // if the user chose to run ensembles on the density using data from SWPC AND */
	      /* 	  	    // if future predictions of F10.7 and Ap (not only past observations) (past values (value before the current time at running) are perfectly known (result from observations, not predictions)) */
	      /* 	  // only overwrite predictions of ensembles (not observations). OPTIONS->et_interpo[aaa] corresponds to the first say of predictions */
	      /* 	  else{ // if we don't run ensembles on F10.7/Ap or that we run ensemble on F10.7/Ap but that the initial epoch is before the first prediction (so there is no uncertainty in F10.7/Ap because the initial epooch corresponds to an observation, not a prediction) */
	      /* 	    //	    print_test(); */
	      /* 	     CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.Ap[0]            = OPTIONS->Ap[0];           // magnetic index(daily) */
	      /* 	    if (OPTIONS->use_ap_hist == 1){ */
	      /* 	      for (hhh = 0; hhh < 7; hhh++){ */
	      /* 		 CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.Ap_hist[hhh][0]            = OPTIONS->Ap_hist[hhh][0];           // magnetic index(historical) */
	      /* 	      } */
	      /* 	    } */

	      /* 	     CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.f107[0]          = OPTIONS->f107[0];         // Daily average of F10.7 flux */
	      /* 	     CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.f107A[0]         = OPTIONS->f107A[0];        // 81 day average of F10.7 flux */

	      /* 	    //  	    printf("%d %f %d %f\n",  CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.sc_ensemble_nb,  CONSTELLATION->spacecraft[ii][eee_bis].INTEGRATOR.f107[index_in_driver_interpolated], index_in_driver_interpolated, OPTIONS->f107[index_in_driver_interpolated]); */
	      /* 	  } // end of if we don't run ensembles on F10.7/Ap or that we run ensemble on F10.7/Ap but that the initial epoch is before the first prediction (so there is no uncertainty in F10.7/Ap because the inital epoch corresponds to an observation, not a prediction) */



	    } // end of if the user chooses a time-varying f107 and Ap for the density 

	  } // end of the user chooses f107 and Ap and Ap_hist for the density

	  else if ( strcmp(OPTIONS->format_density_driver, "density_file") == 0 ){ // the user chooses to directly input the density from a file
	    for (aaa = 0; aaa< OPTIONS->nb_time_steps * 2; aaa++){  // "* 2.0" because of the Runge Kunta order 4 method
	      if (OPTIONS->et_interpo[aaa] <= et_final_epoch + 0.01 ) {
		CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.density[aaa]            = OPTIONS->density[aaa] * 1e9; // Convert from kg/m^3 to kg/km^3
	      }
	    }

	  } // end of  the user chooses to directly input the density from a file

	  else if ( strcmp(OPTIONS->format_density_driver, "gitm") == 0 ){ // the user chooses GITM
	    if ( (ii ==0) && (eee == 0) ) { // here we initialize ATMOSPHERE if gitm is chosen so it doesn't have to be done for all satellites, just the first one (since ATMOSPHERE is the same for all satellites/ensembles)
	      if (OPTIONS->nb_gitm_file < 1){
		printf("***! No GITM file has been found. The program will stop. !***\n");
		MPI_Finalize();
		exit(0);
	      }
	      else{
		for (ggg = 0; ggg < OPTIONS->nb_gitm_file; ggg++){

		  PARAMS->ATMOSPHERE.nb_gitm_file = OPTIONS->nb_gitm_file;
		  strcpy(PARAMS->ATMOSPHERE.array_gitm_file[ggg], OPTIONS->array_gitm_file[ggg]);
		  PARAMS->ATMOSPHERE.array_gitm_date[ggg][0] = OPTIONS->array_gitm_date[ggg][0] ;
		  PARAMS->ATMOSPHERE.array_gitm_date[ggg][1]= OPTIONS->array_gitm_date[ggg][1] ;
		  PARAMS->ATMOSPHERE.is_first_step_of_run = 1;
		}
		// Read the first data file to get the number of lon/lat/alt and allocate memory accordingly

		gitm_file = fopen(PARAMS->ATMOSPHERE.array_gitm_file[0],"rb");

		// NLONS, NLATS, NALTS
		fseek(gitm_file, 20, SEEK_SET	 );
		fread(&PARAMS->ATMOSPHERE.nLons_gitm,sizeof(PARAMS->ATMOSPHERE.nLons_gitm),1,gitm_file);
		fread(&PARAMS->ATMOSPHERE.nLats_gitm,sizeof(PARAMS->ATMOSPHERE.nLats_gitm),1,gitm_file);
		fread(&PARAMS->ATMOSPHERE.nAlts_gitm,sizeof(PARAMS->ATMOSPHERE.nAlts_gitm),1,gitm_file);
		fseek(gitm_file , 4, SEEK_CUR	 );
 
		// Allocate for the GITM file which date is right BEFORE the epoch of the sc

		PARAMS->ATMOSPHERE.longitude_gitm = malloc(  PARAMS->ATMOSPHERE.nLons_gitm * sizeof(double **) );
		for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){
		  PARAMS->ATMOSPHERE.longitude_gitm[i] = malloc(PARAMS->ATMOSPHERE.nLats_gitm * sizeof(double *));
		  for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){
		    PARAMS->ATMOSPHERE.longitude_gitm[i][j] = malloc(PARAMS->ATMOSPHERE.nAlts_gitm * sizeof(double));
		  }
		}
		PARAMS->ATMOSPHERE.latitude_gitm = malloc(  PARAMS->ATMOSPHERE.nLons_gitm * sizeof(double **) );
		for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){
		  PARAMS->ATMOSPHERE.latitude_gitm[i] = malloc(PARAMS->ATMOSPHERE.nLats_gitm * sizeof(double *));
		  for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){
		    PARAMS->ATMOSPHERE.latitude_gitm[i][j] = malloc(PARAMS->ATMOSPHERE.nAlts_gitm * sizeof(double));
		  }
		}
		PARAMS->ATMOSPHERE.altitude_gitm = malloc(  PARAMS->ATMOSPHERE.nLons_gitm * sizeof(double **) );
		for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){
		  PARAMS->ATMOSPHERE.altitude_gitm[i] = malloc(PARAMS->ATMOSPHERE.nLats_gitm * sizeof(double *));
		  for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){
		    PARAMS->ATMOSPHERE.altitude_gitm[i][j] = malloc(PARAMS->ATMOSPHERE.nAlts_gitm * sizeof(double));
		  }
		}


		PARAMS->ATMOSPHERE.density_gitm_right_before = malloc(  PARAMS->ATMOSPHERE.nLons_gitm * sizeof(double **) );
		for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){
		  PARAMS->ATMOSPHERE.density_gitm_right_before[i] = malloc(PARAMS->ATMOSPHERE.nLats_gitm * sizeof(double *));
		  for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){
		    PARAMS->ATMOSPHERE.density_gitm_right_before[i][j] = malloc(PARAMS->ATMOSPHERE.nAlts_gitm * sizeof(double));
		  }
		}

		if (  PARAMS->ATMOSPHERE.longitude_gitm == NULL ){
		  printf("***! Could not allow memory for PARAMS->ATMOSPHERE.longitude_gitm. The program will stop. !***\n");
		  MPI_Finalize();
		  exit(0);
		}
		if (  PARAMS->ATMOSPHERE.latitude_gitm == NULL ){
		  printf("***! Could not allow memory for PARAMS->ATMOSPHERE.latitude_gitm. The program will stop. !***\n");
		  MPI_Finalize();
		  exit(0);
		}
		if (  PARAMS->ATMOSPHERE.altitude_gitm == NULL ){
		  printf("***! Could not allow memory for PARAMS->ATMOSPHERE.altitude_gitm. The program will stop. !***\n");
		  MPI_Finalize();
		  exit(0);
		}

		if (  PARAMS->ATMOSPHERE.density_gitm_right_before == NULL ){
		  printf("***! Could not allow memory for PARAMS->ATMOSPHERE.density_gitm_right_before. The program will stop. !***\n");
		  MPI_Finalize();
		  exit(0);
		}
		// Allocate for the GITM file which date is right AFTER the epoch of the sc


		PARAMS->ATMOSPHERE.density_gitm_right_after = malloc(  PARAMS->ATMOSPHERE.nLons_gitm * sizeof(double **) );


		PARAMS->ATMOSPHERE.density_gitm_right_after = malloc(  PARAMS->ATMOSPHERE.nLons_gitm * sizeof(double **) );
		for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){
		  PARAMS->ATMOSPHERE.density_gitm_right_after[i] = malloc(PARAMS->ATMOSPHERE.nLats_gitm * sizeof(double *));
		  for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){
		    PARAMS->ATMOSPHERE.density_gitm_right_after[i][j] = malloc(PARAMS->ATMOSPHERE.nAlts_gitm * sizeof(double));
		  }
		}


		if (  PARAMS->ATMOSPHERE.density_gitm_right_after == NULL ){
		  printf("***! Could not allow memory for PARAMS->ATMOSPHERE.density_gitm_right_after. The program will stop. !***\n");
		  MPI_Finalize();
		  exit(0);
		}

		// Read the lon/lat/alt array (same for all files in the propagation) and nVars
		// NVARS
		fseek(gitm_file, 4, SEEK_CUR	 );
		fread(&PARAMS->ATMOSPHERE.nVars_gitm,sizeof(PARAMS->ATMOSPHERE.nVars_gitm),1,gitm_file);
		fseek(gitm_file, 4, SEEK_CUR	 );
		iHeaderLength_gitm = 8L + 4+4 +	3*4 + 4+4 + 4 + 4+4 + PARAMS->ATMOSPHERE.nVars_gitm*40 + PARAMS->ATMOSPHERE.nVars_gitm*(4+4) +  7*4 + 4+4; 
		fseek(gitm_file, iHeaderLength_gitm, SEEK_SET);
		// READ THE LONGITUDE
		fseek(gitm_file, 4, SEEK_CUR );
		for (k = 0; k < PARAMS->ATMOSPHERE.nAlts_gitm; k++){
		  for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){
		    for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){
		      fread(&PARAMS->ATMOSPHERE.longitude_gitm[i][j][k], sizeof(PARAMS->ATMOSPHERE.longitude_gitm[i][j][k]), 1, gitm_file);
		    }
		  }
		}
		fseek(gitm_file, 4, SEEK_CUR );
		//printf("\n %f\n", PARAMS->ATMOSPHERE.longitude_gitm[11][40][21]);
     

		// READ THE LATITUDE
		fseek(gitm_file, 4, SEEK_CUR );
		for (k = 0; k < PARAMS->ATMOSPHERE.nAlts_gitm; k++){
		  for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){
		    for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){
		      fread(&PARAMS->ATMOSPHERE.latitude_gitm[i][j][k], sizeof(PARAMS->ATMOSPHERE.latitude_gitm[i][j][k]), 1, gitm_file);
		    }
		  }
		}
		fseek(gitm_file, 4, SEEK_CUR );
		//     printf("\n %f\n", PARAMS->ATMOSPHERE.latitude_gitm[11][40][21]);
     

		// READ THE ALTITUDE
		fseek(gitm_file, 4, SEEK_CUR );
		for (k = 0; k < PARAMS->ATMOSPHERE.nAlts_gitm; k++){
		  for (j = 0; j < PARAMS->ATMOSPHERE.nLats_gitm; j++){
		    for (i = 0; i < PARAMS->ATMOSPHERE.nLons_gitm; i++){
		      fread(&PARAMS->ATMOSPHERE.altitude_gitm[i][j][k], sizeof(PARAMS->ATMOSPHERE.altitude_gitm[i][j][k]), 1, gitm_file);
		      PARAMS->ATMOSPHERE.altitude_gitm[i][j][k] = PARAMS->ATMOSPHERE.altitude_gitm[i][j][k] / 1000.0; // conversion meters to kilometers
		    }
		  }
		}
   
		//     printf("\n %f\n", PARAMS->ATMOSPHERE.altitude_gitm[11][40][21]);

		fclose(gitm_file);
	      }
	      //  printf("%s | %d\n", PARAMS->ATMOSPHERE.array_gitm_file[ggg], ggg);

	    }


	    if (eee == 0){
	      // find the altitude that is right below the perigee altitude. Note: if the duration of the run is long (so that the sc looses a good amount of altitude, so like 6 months) then index_altitude_right_below_perigee is set to 0 so that we go over all the altitudes and not only the ones that start below the perigee calcualted at the initialization.
	      //if many main satellites with different altitude of the perigee, then take the min of these perigee to calculate index_altitude_right_below_perigee
	      int found_altitude_right_below_perigee = 0;
	      if ( et_final_epoch - et_initial_epoch < 6 * 31 * 24 * 3600.0 ){ // if the duration of the run is long (so that the sc looses a good amount of altitude, so like 6 months) then index_altitude_right_below_perigee is set to 0 so that we go over all the altitudes and not only the ones that start below the perigee calcualted at the initialization.
		altitude_perigee = radius_perigee - PARAMS->EARTH.radius;
		k = 0;
		while ( (k < PARAMS->ATMOSPHERE.nAlts_gitm ) && (found_altitude_right_below_perigee == 0) ){
		  if (PARAMS->ATMOSPHERE.altitude_gitm[0][0][k] >= altitude_perigee){
		    if ( (k-1) < index_altitude_right_below_perigee_save ){
		      index_altitude_right_below_perigee_save = k-1;
		    } 
		    found_altitude_right_below_perigee = 1;
		  }
		  k = k+1;
		}
		if (ii == OPTIONS->n_satellites - OPTIONS->nb_gps - 1){
		  PARAMS->ATMOSPHERE.index_altitude_right_below_perigee = index_altitude_right_below_perigee_save;

		}
	      }
	      else{
		PARAMS->ATMOSPHERE.index_altitude_right_below_perigee = 0;
	      }
	    }
	
	  } // end of the user chooses GITM

	} // end of if the user wants to use drag
	if (iDebugLevel >= 2){
	  if (iProc == 0) printf("--- (initialize_constellation) Done initializing dt, nb_surfaces, mass, solar_cell_efficiency, degree, order, include_drag/solar_pressure/earth_pressure/moon/sun, Ap, Ap_hist, f107, f107A, density, initialize_geo_with_bstar for all spacecraft.\n");
	}


	if (iDebugLevel >= 1){
	  if (iProc == 0) printf("-- (initialize_constellation) Done initializing r_ecef2cg_ECEF, geodetic, output file names, name_sat, dt, nb_surfaces, mass, solar_cell_efficiency, degree, order, include_drag/solar_pressure/earth_pressure/moon/sun, Ap, Ap_hist, f107, f107A, density, initialize_geo_with_bstar for all spacecraft.\n");
	}

	/*************************************************************************************/
	/*************************************************************************************/
	/************* END OF COE AND TLE INITIALIZE for main SC and for ensemble SC:
- r_ecef2cg_ECEF
- geodetic
- filename, filenameecef, filenameout, filenamepower for main SC ONLY
- name_sat
- INTEGRATOR: dt, nb_surfaces, mass, solar_cell_efficiency, degree, order, include_drag/solar_pressure/earth_pressure/moon/sun, Ap, Ap_hist, f107, f107A, density, initialize_geo_with_bstar, sc_main_nb, sc_ensemble_nb
	**************/
	/*************************************************************************************/
	/*************************************************************************************/



	/*************************************************************************************/
	/*************************************************************************************/
	/*********************** COE AND TLE INITIALIZE: attitude ****************************/
	/*************************************************************************************/
	/*************************************************************************************/
	if (iDebugLevel >= 1){
	  if (iProc == 0) printf("-- (initialize_constellation) Initializing the attitude.\n");
	}

	if ( ( strcmp(OPTIONS->attitude_profile, "ensemble_angular_velocity") != 0 ) && ( strcmp(OPTIONS->attitude_profile, "ensemble_initial_attitude") != 0 ) ) { // if we do not run ensembles on the initial angular velocity

	  // Give memory to attitude variables
	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.file_is_quaternion = OPTIONS->file_is_quaternion;
		if ( OPTIONS->file_is_quaternion == 0){
	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) ); // "* 2.0" because of the Runge Kunta order 4 method

	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_pitch = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) ); 
	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_roll = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_yaw = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
		}
		else{
	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.quaternion = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double*) ); // "* 2.0" because of the Runge Kunta order 4 method
	  for (ccc = 0; ccc < OPTIONS->nb_time_steps * 2; ccc++){
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.quaternion[ccc] = malloc( 4* sizeof(double) );
	  }


		}
		 CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.index_in_attitude_interpolated_first = (int)( CONSTELLATION->spacecraft[ii][eee].et - OPTIONS->et_oldest_tle_epoch ) / OPTIONS->dt *2.;
		  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.index_in_driver_interpolated_first = (int)( CONSTELLATION->spacecraft[ii][eee].et - OPTIONS->et_oldest_tle_epoch ) / OPTIONS->dt * 2.;

	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.index_in_attitude_interpolated = (int)( CONSTELLATION->spacecraft[ii][eee].et - OPTIONS->et_oldest_tle_epoch ) / OPTIONS->dt *2.;
	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.index_in_driver_interpolated = (int)( CONSTELLATION->spacecraft[ii][eee].et - OPTIONS->et_oldest_tle_epoch ) / OPTIONS->dt * 2.;
	  //	  printf("%d %d | %d %d | %d %d\n", ii, eee, CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.index_in_attitude_interpolated, CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.index_in_driver_interpolated, CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.index_in_attitude_interpolated_first , CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.index_in_driver_interpolated_first );
		if ( OPTIONS->file_is_quaternion == 1){
	  if (  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.quaternion == NULL ){
	    printf("***! Could not allow memory space for  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.quaternion \n. The program will stop. !***\n");
	    MPI_Finalize();
	    exit(0);
	  }
		}
		else{
	  if (  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch == NULL ){
	    printf("***! Could not allow memory space for  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch \n. The program will stop. !***\n");
	    MPI_Finalize();
	    exit(0);
	  }
	  if (  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll == NULL ){
	    printf("***! Could not allow memory space for  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll \n. The program will stop. !***\n");
	    MPI_Finalize();
	    exit(0);
	  }
	  if (  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw == NULL ){
	    printf("***! Could not allow memory space for  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw \n. The program will stop. !***\n");
	    MPI_Finalize();
	    exit(0);
	  }
	  if (  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_pitch == NULL ){
	    printf("***! Could not allow memory space for  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_pitch \n. The program will stop. !***\n");
	    MPI_Finalize();
	    exit(0);
	  }
	  if (  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_roll == NULL ){
	    printf("***! Could not allow memory space for  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_roll \n. The program will stop. !***\n");
	    MPI_Finalize();
	    exit(0);
	  }
	  if (  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_yaw == NULL ){
	    printf("***! Could not allow memory space for  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_yaw \n. The program will stop. !***\n");
	    MPI_Finalize();
	    exit(0);
	  }
		}
	  strcpy(CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.attitude_profile, OPTIONS->attitude_profile); // nadir, sun_pointed, ...
	  /*************************************************************************************/
	  /******************* COE AND TLE INITIALIZE: attitude for main SC ********************/
	  /*************************************************************************************/
	  if (eee == 0){ // eee = 0 represents the main spacecraft (eee > 0 is a sc from an ensemble)

	    for (aaa = 0; aaa< OPTIONS->nb_time_steps * 2; aaa++){  // "* 2.0" because of the Runge Kunta order 4 method
	      if (OPTIONS->et_interpo[aaa] <= et_final_epoch + 0.01) {
		if ( OPTIONS->file_is_quaternion == 0){

		CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch[aaa] = OPTIONS->pitch[aaa];
		CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll[aaa]  = OPTIONS->roll[aaa];
		CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw[aaa]   = OPTIONS->yaw[aaa];
		CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_pitch[aaa] = OPTIONS->order_pitch[aaa];
		CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_roll[aaa]  = OPTIONS->order_roll[aaa];
		CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_yaw[aaa]   = OPTIONS->order_yaw[aaa];
		  }
		else{

		CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.quaternion[aaa][0] = OPTIONS->quaternion[aaa][0];

		CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.quaternion[aaa][1] = OPTIONS->quaternion[aaa][1];
		CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.quaternion[aaa][2] = OPTIONS->quaternion[aaa][2];
		CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.quaternion[aaa][3] = OPTIONS->quaternion[aaa][3];

		}
	      }
	    }

	  } // end of initializing the attitude for the main SC

	  /*************************************************************************************/
	  /******************* COE AND TLE INITIALIZE: attitude for ensemble SC ****************/
	  /*************************************************************************************/

	  /*****************************************************************************************************************************/
	  /**************************************************** ENSEMBLES ON ATTITUDE **************************************************/
	  /*****************************************************************************************************************************/
	  else{ // here eee > 0 so we are initializing the attitude for the sc from an ensemble

	    if (OPTIONS->nb_ensembles_attitude > 0){ // initialize the attitude if we run ensembles on the attitude. So it's not with ensemble_angular_velocity and not with ensemble_initial_attitude (these would have had to be put in section #ATTITUDE). But it corresponds to a drift of the sc from a reference attitude (specified in section #ATTITUDE (nadir, sun_pointed, ...)) for with a given time (set by time_before_reset) with a random angular velocity.

	      //	    if (iProc == 0){
	      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.pitch_sigma_ensemble = OPTIONS->pitch_sigma_ensemble;
	      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.roll_sigma_ensemble = OPTIONS->roll_sigma_ensemble;
	      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.yaw_sigma_ensemble = OPTIONS->yaw_sigma_ensemble;
	      //	    }

	      //for (aaa = 0; aaa< OPTIONS->nb_time_steps * 2; aaa++){  // "* 2.0" because of the Runge Kunta order 4 method
	      //	  last_aaa_previous_loop = 0;
	      aaa = 0;
	      if (write_attitude_file == 1){
		fprintf(file_attitude[ii], "#START_ATTITUDE_ENSEMBLE for ensemble %d\n", eee);
	      }


	      while ( aaa < OPTIONS->nb_time_steps * 2 ){
		//		printf("%d %d %d || %f %f\n", iProc, aaa, OPTIONS->nb_time_steps * 2, OPTIONS->et_interpo[aaa], et_final_epoch);
		if (OPTIONS->et_interpo[aaa] <= et_final_epoch + 0.01 ) {

		  time_before_reset = 0;
		  // every OPTIONS->attitude_reset_delay seconds, the attitude is set equal to the attitude of the main spacecraft (nadir or sun pointed or from an atttitude file)
		  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch[aaa] =  OPTIONS->pitch[aaa];
		  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll[aaa] = OPTIONS->roll[aaa];
		  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw[aaa] = OPTIONS->yaw[aaa];
		  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_pitch[aaa] = OPTIONS->order_pitch[aaa];
		  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_roll[aaa]  = OPTIONS->order_roll[aaa];
		  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_yaw[aaa]   = OPTIONS->order_yaw[aaa];
		  // during OPTIONS->attitude_reset_delay seconds, the satellite's attitude changes with a constant random angular velocity
		  random_pitch_angular_velocity = randn( 0.0, OPTIONS->pitch_sigma_angular_velocity_ensemble);
		  random_roll_angular_velocity = randn( 0.0, OPTIONS->roll_sigma_angular_velocity_ensemble);
		  random_yaw_angular_velocity = randn( 0.0, OPTIONS->yaw_sigma_angular_velocity_ensemble);
		  if (write_attitude_file == 1){
		    if (aaa==0){
		      fprintf(file_attitude[ii], "%e %e %e %d %d %d\n", random_pitch_angular_velocity, random_roll_angular_velocity, random_yaw_angular_velocity, CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_pitch[aaa], CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_roll[aaa],CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_yaw[aaa] );		      
		    }
		    et2utc_c(OPTIONS->et_oldest_tle_epoch+aaa*OPTIONS->dt/2., "ISOC" ,0 ,255 , time_attitude);
		    fprintf(file_attitude[ii], "%s %e %e %e\n", time_attitude,CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch[aaa], CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll[aaa], CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw[aaa]);
		  }
		  //		    if (iProc == 0){
/* 		      if ( (ii == 1) && (eee == array_sc[start_ensemble[ii]]) ){ */
/* 			printf("iProc %d | aaa = %d out of %d | time reset %f out of %f\n", iProc, aaa, OPTIONS->nb_time_steps * 2, time_before_reset * OPTIONS->dt / 2.0 , OPTIONS->attitude_reset_delay); */
/* 		      } */
		      //}


		  while ( ( time_before_reset * OPTIONS->dt / 2.0 < ( OPTIONS->attitude_reset_delay - OPTIONS->dt / 2.0 ) ) && ( aaa < OPTIONS->nb_time_steps * 2 ) ){
		    time_before_reset = time_before_reset + 1 ;
		    aaa = aaa + 1;

		    /* 	if (iProc == 2){ */
/* 				  printf("iProc %d | aaa = %d out of %d | time reset %f out of %f\n", iProc, aaa, OPTIONS->nb_time_steps * 2, time_before_reset * OPTIONS->dt / 2.0 , OPTIONS->attitude_reset_delay); */
/* 				  	  } */

		    if (aaa < OPTIONS->nb_time_steps * 2){
		      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch[aaa] = random_pitch_angular_velocity * time_before_reset * OPTIONS->dt / 2.0 + OPTIONS->pitch[aaa]; // the angular velocity is distributed as gaussian with a 0 rad/s mean and a standard deviation OPTIONS->pitch/roll/yaw_sigma_angular_velocity_ensemble chosen by the user. The attitude is calculated from this random angular velocity around the attitude of the main spacecraft (nadir or from an attitude file (FOR NOW NOT POSSIBLE TO RUN ENSEMBLES IF ATTITUDE IS SUN_POINTED)): OPTIONS->pitch/roll/yaw[aaa]
		      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll[aaa] = random_roll_angular_velocity * time_before_reset * OPTIONS->dt / 2.0 + OPTIONS->roll[aaa];
		      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw[aaa] = random_yaw_angular_velocity * time_before_reset * OPTIONS->dt / 2.0 + OPTIONS->yaw[aaa];
		      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_pitch[aaa] = OPTIONS->order_pitch[aaa];
		      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_roll[aaa]  = OPTIONS->order_roll[aaa];
		      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_yaw[aaa]   = OPTIONS->order_yaw[aaa];
/* 		    if (iProc == 0){ */
/* 		      if ( (ii == 0) && (eee == 1) ){ */
/* 			printf("IN iProc %d | aaa = %d out of %d | time reset %f out of %f | pitch %f\n", iProc, aaa, OPTIONS->nb_time_steps * 2, time_before_reset * OPTIONS->dt / 2.0 , OPTIONS->attitude_reset_delay, CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch[aaa]); */
/* 		      } */
/* 		    } */

		      if (write_attitude_file == 1){
			if (fmod(aaa, 2) == 0){
			  et2utc_c(OPTIONS->et_oldest_tle_epoch+aaa*OPTIONS->dt/2., "ISOC" ,0 ,255 , time_attitude);
			  fprintf(file_attitude[ii], "%s %e %e %e\n", time_attitude,CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch[aaa], CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll[aaa], CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw[aaa] );
			}
		      }
		  
		    }
		  } // end of while loop for the attitude with a constant random angular velocity
		  aaa = aaa + 1;
		  //last_aaa_previous_loop = aaa;
		}
	      } // end of initializing the attitude for the ensemble spacecrafts
	      if (write_attitude_file == 1){
		fprintf(file_attitude[ii], "#END_ATTITUDE_ENSEMBLE\n\n");
	      }

	    } // end of initialize the attitude if we run ensembles on the attitude
	  
	    else{ // initialize the attitude if we do not run ensembles on the attitude
	      // WE DON'T INITALIZE THE ATTIUDE HERE ANYMORE BUT IN THE FUNCTION SET_ATTITUDE CALLED IN PROPAGATE_SPACECRAFT. This is to avoid going through all time steps for all ensembles. Note that for all other cases, ie if we run any kind of ensembles on the attitude, the initialization of the attitude is still set in initialize_constellation (and not in propagate_spacecraft)
	      /* for (aaa = 0; aaa< OPTIONS->nb_time_steps * 2; aaa++){ */
	      /*   if (OPTIONS->et_interpo[aaa] <= et_final_epoch + 0.01 ) { */
	      /* 	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch[aaa] = OPTIONS->pitch[aaa]; */
	      /* 	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll[aaa]  = OPTIONS->roll[aaa]; */
	      /* 	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw[aaa]   = OPTIONS->yaw[aaa]; */
	      /* 	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_pitch[aaa] = OPTIONS->order_pitch[aaa]; */
	      /* 	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_roll[aaa]  = OPTIONS->order_roll[aaa]; */
	      /* 	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.order_yaw[aaa]   = OPTIONS->order_yaw[aaa]; */
	      /*   } */
	      /* } */
	      //end of  WE DON'T INITALIZE THE ATTIUDE HERE ANYMORE BUT IN THE FUNCTION SET_ATTITUDE CALLED IN PROPAGATE_SPACECRAFT
	    } // end of initialize the attitude if we do not run ensembles on the attitude
	  } // end of initializinf the attitude for the ensemble spacecraft

	} // end of if we do not run ensembles on the initial angular velocity

	/*********************************************************************************************************************************/
	/******************************************* ENSEMBLES ON INITIAL ANGULAR VELOCITY ***********************************************/
	/*********************************************************************************************************************************/
	else if ( strcmp(OPTIONS->attitude_profile, "ensemble_angular_velocity") == 0 ){  // if we run ensembles on the initial angular velocity. The attitude is defined by its initial angle (pitch; roll; yaw), the mean angular velocity (mean ang velo pitch; mean ang velo roll; mean ang velo yaw), and the standard deviation on the angular velo (sigma ang velo pitch; sigma ang velo roll; sigma ang velo yaw)
	  //	if ( (eee == 0) && (ii == 0)){ // initialize attitude for main spacecraft if  we run ensembles on the initial angular velocity // !!!!!!!!!!!! for now works only if there one satellite in the constellation
	  if ( eee == 0){ // initialize attitude for main spacecraft if  we run ensembles on the initial angular velocity 
	    strcpy(CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.attitude_profile, OPTIONS->attitude_profile); // here attitude_profile = "ensemble_angular_velocity"
	    // pitch
	    CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.pitch_ini_ensemble = OPTIONS->pitch_ini_ensemble;
	    CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.pitch_angular_velocity_ensemble = OPTIONS->pitch_mean_angular_velocity_ensemble;
	    CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.pitch_sigma_angular_velocity_ensemble = OPTIONS->pitch_sigma_angular_velocity_ensemble / 2.0;
	    // roll
	    CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.roll_ini_ensemble = OPTIONS->roll_ini_ensemble;
	    CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.roll_angular_velocity_ensemble = OPTIONS->roll_mean_angular_velocity_ensemble;
	    CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.roll_sigma_angular_velocity_ensemble = OPTIONS->roll_sigma_angular_velocity_ensemble / 2.0;
	    // yaw
	    CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.yaw_ini_ensemble = OPTIONS->yaw_ini_ensemble;
	    CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.yaw_angular_velocity_ensemble = OPTIONS->yaw_mean_angular_velocity_ensemble;
	    CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.yaw_sigma_angular_velocity_ensemble = OPTIONS->yaw_sigma_angular_velocity_ensemble / 2.0;

	  } // end of initialize attitude for main spacecraft if  we run ensembles on the initial angular velocity
	  //	else if (ii == 0){ // initialize attitude for ensemble spacecraft if  we run ensembles on the initial angular velocity
	  else{ // initialize attitude for ensemble spacecraft if  we run ensembles on the initial angular velocity
	    strcpy(CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.attitude_profile, OPTIONS->attitude_profile); // here attitude_profile = "ensemble_angular_velocity"
	    // calculate the random angular velocity = normal distribution around the mean angular velocity OPTIONS->pitch/roll/yaw_mean_angular_velocity_ensemble with a standard deviation OPTIONS->pitch/roll/yaw_sigma_angular_velocity_ensemble

	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch_angular_velocity_ensemble = randn( OPTIONS->pitch_mean_angular_velocity_ensemble, OPTIONS->pitch_sigma_angular_velocity_ensemble / 2.0);
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll_angular_velocity_ensemble = randn( OPTIONS->roll_mean_angular_velocity_ensemble, OPTIONS->roll_sigma_angular_velocity_ensemble / 2.0);
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw_angular_velocity_ensemble = randn( OPTIONS->yaw_mean_angular_velocity_ensemble, OPTIONS->yaw_sigma_angular_velocity_ensemble / 2.0);
	     
	    //	     fprintf(fp_temp, "%f \n", CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch_angular_velocity_ensemble);
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch_ini_ensemble = OPTIONS->pitch_ini_ensemble;
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll_ini_ensemble = OPTIONS->roll_ini_ensemble;
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw_ini_ensemble = OPTIONS->yaw_ini_ensemble;


	    /* 	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch_angular_velocity_ensemble = randn( CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.pitch_angular_velocity_ensemble, CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.pitch_sigma_angular_velocity_ensemble); */
	    /* 	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll_angular_velocity_ensemble = randn( CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.roll_angular_velocity_ensemble, CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.roll_sigma_angular_velocity_ensemble); */
	    /* 	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw_angular_velocity_ensemble = randn( CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.yaw_angular_velocity_ensemble, CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.yaw_sigma_angular_velocity_ensemble); */
	     
	    /* 	  //	     fprintf(fp_temp, "%f \n", CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch_angular_velocity_ensemble); */
	    /* 	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch_ini_ensemble = CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.pitch_ini_ensemble; */
	    /* 	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll_ini_ensemble = CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.roll_ini_ensemble; */
	    /* 	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw_ini_ensemble = CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.yaw_ini_ensemble; */

	  } // end of initialize attitude for ensemble spacecraft if  we run ensembles on the initial angular velocity
	} // end of if we run ensembles on the initial angular velocity
	else if ( strcmp(OPTIONS->attitude_profile, "ensemble_initial_attitude") == 0 ){ //if we run ensembles on the initial attitude. The ensembles all have different initial attitude but same constant angular velocities
	  //	if ( (eee == 0) && (ii == 0)){ // initialize attitude for main spacecraft if we run ensembles on the initial attitude // !!!!!!!!!!!! for now works only if there one satellite in the constellation
	  if ( eee == 0 ){ // initialize attitude for main spacecraft if we run ensembles on the initial attitude 
	    strcpy(CONSTELLATION->spacecraft[ii][0].INTEGRATOR.attitude.attitude_profile, OPTIONS->attitude_profile); // here attitude_profile = "ensemble_initial_attitude"
	    // the reference satellite has the same attitude as the mean attitude chosen by the user (first line of ###ENSEMBLES_INITIAL_ATTITUDE)

	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch_for_attitude_ensemble = OPTIONS->pitch_mean_ensemble;
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll_for_attitude_ensemble = OPTIONS->roll_mean_ensemble;
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw_for_attitude_ensemble = OPTIONS->yaw_mean_ensemble;

	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch_angular_velocity_constant = OPTIONS->pitch_angular_velocity_constant;
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll_angular_velocity_constant = OPTIONS->roll_angular_velocity_constant;
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw_angular_velocity_constant = OPTIONS->yaw_angular_velocity_constant;
	  
	  } // end of initialize attitude for main spacecraft if  we run ensembles on the initial angular velocity
	  //	else if (ii == 0){ // initialize attitude for ensemble spacecraft if  we run ensembles on the initial angular velocity
	  else{ // initialize attitude for ensemble spacecraft if  we run ensembles on the initial angular velocity

	    strcpy(CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.attitude_profile, OPTIONS->attitude_profile); // here attitude_profile = "ensemble_angular_velocity"
	    // calculate the random angular velocity = normal distribution around the mean angular velocity OPTIONS->pitch/roll/yaw_mean_angular_velocity_ensemble with a standard deviation OPTIONS->pitch/roll/yaw_sigma_angular_velocity_ensemble

	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch_for_attitude_ensemble = randn( OPTIONS->pitch_mean_ensemble, OPTIONS->pitch_sigma_for_ensemble_initial_attitude / 2.);
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll_for_attitude_ensemble = randn( OPTIONS->roll_mean_ensemble, OPTIONS->roll_sigma_for_ensemble_initial_attitude / 2.);
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw_for_attitude_ensemble = randn( OPTIONS->yaw_mean_ensemble, OPTIONS->yaw_sigma_for_ensemble_initial_attitude / 2.);

	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.pitch_angular_velocity_constant = OPTIONS->pitch_angular_velocity_constant;
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.roll_angular_velocity_constant = OPTIONS->roll_angular_velocity_constant;
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude.yaw_angular_velocity_constant = OPTIONS->yaw_angular_velocity_constant;

	  } // end of initialize attitude for ensemble spacecraft if  we run ensembles on the initial angular velocity
      
	} // end of if we run ensembles on the initial attitude
	if (iDebugLevel >= 1){
	  if (iProc == 0) printf("-- (initialize_constellation) Done initializing the attitude.\n");
	}


	/*************************************************************************************/
	/*************************************************************************************/
	/********************** END OF COE AND TLE INITIALIZE: attitude **********************/
	/*************************************************************************************/
	/*************************************************************************************/
    
	/*************************************************************************************/
	/*************************************************************************************/
	/*********************** COE AND TLE INITIALIZE: surface *****************************/
	/*************************************************************************************/
	/*************************************************************************************/

	if (iDebugLevel >= 1){
	  if (iProc == 0) printf("-- (initialize_constellation) Initializing the geometry.\n");
	}

	for (sss = 0; sss < OPTIONS->n_surfaces; sss++){

	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].area                  = OPTIONS->surface[sss].area / (1e10);

	  // !!!!!!!!!!!!!!!!!! TO ERASE
	  /* if (ii == 3 || ii == 7){ */
	  /*   	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[8].area                  = 5000.0 / (1e10); */
	  /* } */
	  // !!!!!!!!!!!!!!!!!! END OF TO ERASE

	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].area_solar_panel      = OPTIONS->surface[sss].area_solar_panel / (1e10);
	  /*************************************************************************************/
	  /*********************** COE AND TLE INITIALIZE: ensembles on Cd *****************************/
	  /*************************************************************************************/
	  if (OPTIONS->nb_ensembles_cd > 0){ // if we run ensembles on Cd
	    if (eee == 0){
	      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].Cd = OPTIONS->surface[sss].Cd * OPTIONS->cd_modification[ii];
	    }
	    else{
	      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].Cd = randn( OPTIONS->surface[sss].Cd * OPTIONS->cd_modification[ii], OPTIONS->surface[sss].Cd_sigma);
	      //	    	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].Cd = OPTIONS->surface[sss].Cd * 2 * eee;  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ERASE THIS LINE AND UNCOMMENT THE PREVIOUS ONE!!!
	    }
	
	  } // end of if we run ensembles on Cd
	  else{ // if we do not run ensembles on Cd
	    // new cd
	    if (OPTIONS->new_cd == 1){
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].acco_coeff  = 
	      OPTIONS->surface[sss].acco_coeff;        
	    }
	    else{

	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].Cd  =
	      OPTIONS->surface[sss].Cd * OPTIONS->cd_modification[ii];           // Coefficient of drag
	    }
	    // end of new cd
	    // !!!!!!!!!!!!!!!!!!!!!!!! ERASE THIS BLOCK BELOW
	    //if ( ( ii == 2 ) || ( ii == 7 ) ){
	    //CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].Cd                    = OPTIONS->surface[sss].Cd*6;
	    //}
	    // !!!!!!!!!!!!!!!!!!!!!!!! END OF ERASE THIS BLOCK BELOW

	  }// end of if we do not run ensembles on Cd

	  // !!!!!!!!!!! THIS BLOCK IS USING THE EQUATION FROM STK (http://www.agi.com/resources/help/online/stk/10.1/index.html?page=source%2Fhpop%2Fhpop-05.htm). UNCOMMENT THE BLOCK BELOW THAT USES VALLADO AND COMMENT THIS STK BLOCK IF YOU WANT TO USE VALLADO'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)
	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].solar_radiation_coefficient  = OPTIONS->surface[sss].solar_radiation_coefficient;
	  // !!!!!!!!!!! END OF THIS BLOCK IS USING THE EQUATION FROM STK (http://www.agi.com/resources/help/online/stk/10.1/index.html?page=source%2Fhpop%2Fhpop-05.htm). UNCOMMENT THE BLOCK BELOW THAT USES VALLADO AND COMMENT THIS STK BLOCK IF YOU WANT TO USE VALLADO'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)
	  // !!!!!!!!!! THIS BLOCK USES VALLADO'S EQUATIONS. COMMENT IT AND UNCOMMENT THE BLOCK ABOVE THAT USES STK IF YOU WANT TO USE STK'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)
	  /* CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].specular_reflectivity = OPTIONS->surface[sss].specular_reflectivity; */
	  /* CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].diffuse_reflectivity  = OPTIONS->surface[sss].diffuse_reflectivity; */
	  // !!!!!!!!!! END OF THIS BLOCK USES VALLADO'S EQUATIONS. COMMENT IT AND UNCOMMENT THE BLOCK ABOVE THAT USES STK IF YOU WANT TO USE STK'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)

	  strcpy(CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].name_of_surface, OPTIONS->surface[sss].name_of_surface); //

	  for (nnn = 0; nnn < 3; nnn++){
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].normal[nnn]         = OPTIONS->surface[sss].normal[nnn];           // normal vector in the SC reference system
	  }
	}
	for (sss = 0; sss < OPTIONS->n_surfaces_eff; sss++){

	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface_eff[sss].area                  = OPTIONS->surface_eff[sss].area / (1e10);

	  // !!!!!!!!!!!!!!!!!! TO ERASE
	  /* if (ii == 3 || ii == 7){ */
	  /*   	CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[8].area                  = 5000.0 / (1e10); */
	  /* } */
	  // !!!!!!!!!!!!!!!!!! END OF TO ERASE


	  /*************************************************************************************/
	  /*********************** COE AND TLE INITIALIZE: ensembles on Cd *****************************/
	  /*************************************************************************************/
	  if (OPTIONS->nb_ensembles_cd > 0){ // if we run ensembles on Cd
	    if (eee == 0){
	      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface_eff[sss].Cd = OPTIONS->surface[0].Cd * OPTIONS->cd_modification[ii]; // !!!!!!!!! assumes that for all effective have the same cd
	    }
	    else{
	      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface_eff[sss].Cd = randn( OPTIONS->surface[0].Cd * OPTIONS->cd_modification[ii], OPTIONS->surface[0].Cd_sigma); // !!!!!!!!! assumes that for all effective have the same cd
	      //	    	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].Cd = OPTIONS->surface[sss].Cd * 2 * eee;  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ERASE THIS LINE AND UNCOMMENT THE PREVIOUS ONE!!!
	    }
	
	  } // end of if we run ensembles on Cd
	  else{ // if we do not run ensembles on Cd
	    // new cd
	    if (OPTIONS->new_cd == 1){
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface_eff[sss].acco_coeff  = 
	      OPTIONS->surface[0].acco_coeff;        // !!!!!!!!! assumes that for all effective have the same acco_coeff
	    }
	    else{

	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface_eff[sss].Cd  =
	      OPTIONS->surface[0].Cd * OPTIONS->cd_modification[ii];           // Coefficient of drag
	    }
	    // end of new cd
	    // !!!!!!!!!!!!!!!!!!!!!!!! ERASE THIS BLOCK BELOW
	    //if ( ( ii == 2 ) || ( ii == 7 ) ){
	    //CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].Cd                    = OPTIONS->surface[sss].Cd*6;
	    //}
	    // !!!!!!!!!!!!!!!!!!!!!!!! END OF ERASE THIS BLOCK BELOW

	  }// end of if we do not run ensembles on Cd

	  // !!!!!!!!!!! THIS BLOCK IS USING THE EQUATION FROM STK (http://www.agi.com/resources/help/online/stk/10.1/index.html?page=source%2Fhpop%2Fhpop-05.htm). UNCOMMENT THE BLOCK BELOW THAT USES VALLADO AND COMMENT THIS STK BLOCK IF YOU WANT TO USE VALLADO'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)
	  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface_eff[sss].solar_radiation_coefficient  = OPTIONS->surface[0].solar_radiation_coefficient; // !!!!!!!!! assumes that for all effective have the same solar_radiation_coefficient
	  // !!!!!!!!!!! END OF THIS BLOCK IS USING THE EQUATION FROM STK (http://www.agi.com/resources/help/online/stk/10.1/index.html?page=source%2Fhpop%2Fhpop-05.htm). UNCOMMENT THE BLOCK BELOW THAT USES VALLADO AND COMMENT THIS STK BLOCK IF YOU WANT TO USE VALLADO'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)
	  // !!!!!!!!!! THIS BLOCK USES VALLADO'S EQUATIONS. COMMENT IT AND UNCOMMENT THE BLOCK ABOVE THAT USES STK IF YOU WANT TO USE STK'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)
	  /* CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].specular_reflectivity = OPTIONS->surface[sss].specular_reflectivity; */
	  /* CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface[sss].diffuse_reflectivity  = OPTIONS->surface[sss].diffuse_reflectivity; */
	  // !!!!!!!!!! END OF THIS BLOCK USES VALLADO'S EQUATIONS. COMMENT IT AND UNCOMMENT THE BLOCK ABOVE THAT USES STK IF YOU WANT TO USE STK'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)

	  strcpy(CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface_eff[sss].name_of_surface, "Area effective - ignore name"); //

	  for (nnn = 0; nnn < 3; nnn++){
	    CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface_eff[sss].normal[nnn]         = OPTIONS->surface_eff[sss].normal[nnn];           // normal vector in the SC reference system
	  }
	} // end of surface effective
	
  if (OPTIONS->opengl == 1){
  
    CONSTELLATION->area_attitude_opengl_phi0 = OPTIONS->area_attitude_opengl_phi0;
  CONSTELLATION->area_attitude_opengl_theta0 = OPTIONS->area_attitude_opengl_theta0;
  CONSTELLATION->area_attitude_opengl_dtheta = OPTIONS->area_attitude_opengl_dtheta;
  CONSTELLATION->area_attitude_opengl_dphi = OPTIONS->area_attitude_opengl_dphi;

  int nb_phi, nb_theta;
  nb_theta = (int)((180 - CONSTELLATION->area_attitude_opengl_theta0) / CONSTELLATION->area_attitude_opengl_dtheta ) + 1;// 180 is included (because different from theta = 0)    
  nb_phi = (int)((360 - CONSTELLATION->area_attitude_opengl_phi0) / CONSTELLATION->area_attitude_opengl_dphi ) ; // 360 is not included (because same as phi = 0)


  CONSTELLATION->area_attitude_opengl = NULL;
  CONSTELLATION->area_attitude_opengl = malloc(nb_theta * sizeof(double **));
  if ( CONSTELLATION->area_attitude_opengl == NULL){
    printf("***! Could not allow memory to CONSTELLATION->area_attitude_opengl. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }

  CONSTELLATION->which_face_theta_phi = NULL;
  CONSTELLATION->which_face_theta_phi = malloc(nb_theta * sizeof(int **));
  if ( CONSTELLATION->which_face_theta_phi == NULL){
    printf("***! Could not allow memory to CONSTELLATION->which_face_theta_phi. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }


  CONSTELLATION->area_attitude_opengl_total = NULL;
  CONSTELLATION->area_attitude_opengl_total = malloc(nb_theta * sizeof(double *));
  if ( CONSTELLATION->area_attitude_opengl_total == NULL){
    printf("***! Could not allow memory to CONSTELLATION->area_attitude_opengl_total. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }

  if (( OPTIONS->solar_cell_efficiency != -1) && (OPTIONS->opengl_power == 1)){
    CONSTELLATION->area_solar_panel_attitude_opengl = NULL;
  CONSTELLATION->area_solar_panel_attitude_opengl = malloc(nb_theta * sizeof(double *));
  if ( CONSTELLATION->area_solar_panel_attitude_opengl == NULL){
    printf("***! Could not allow memory to CONSTELLATION->area_solar_panel_attitude_opengl. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }
  }


  CONSTELLATION->nb_faces = OPTIONS->nb_faces;
  CONSTELLATION->nb_faces_theta_phi = NULL;
  CONSTELLATION->nb_faces_theta_phi = malloc(nb_theta * sizeof(double *));
  if ( CONSTELLATION->nb_faces_theta_phi == NULL){
    printf("***! Could not allow memory to CONSTELLATION->nb_faces_theta_phi. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }

  CONSTELLATION->normal_face = NULL;
  CONSTELLATION->normal_face = malloc(OPTIONS->nb_faces * sizeof(double *));
  if ( CONSTELLATION->normal_face == NULL){
    printf("***! Could not allow memory to CONSTELLATION->normal_face. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }

  int itheta, iphi, iface, iface_here;
  for (iface = 0; iface < OPTIONS->nb_faces; iface++){
    CONSTELLATION->normal_face[iface] = malloc(3 * sizeof(double));
    if ( CONSTELLATION->normal_face[iface] == NULL){
      printf("***! Could not allow memory to CONSTELLATION->normal_face[iface]. The program will stop. !***\n"); MPI_Finalize(); exit(0);
    }
    CONSTELLATION->normal_face[iface][0] = OPTIONS->normal_face[iface][0];
    CONSTELLATION->normal_face[iface][1] = OPTIONS->normal_face[iface][1];
    CONSTELLATION->normal_face[iface][2] = OPTIONS->normal_face[iface][2];
  }
  
  for (itheta = 0; itheta < nb_theta; itheta++){
    CONSTELLATION->area_attitude_opengl[itheta] = malloc(nb_phi * sizeof(double *));
  if ( CONSTELLATION->area_attitude_opengl[itheta] == NULL){
    printf("***! Could not allow memory to CONSTELLATION->area_attitude_opengl[itheta]. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }
    CONSTELLATION->which_face_theta_phi[itheta] = malloc(nb_phi * sizeof(int *));
  if ( CONSTELLATION->which_face_theta_phi[itheta] == NULL){
    printf("***! Could not allow memory to CONSTELLATION->which_face_theta_phi[itheta]. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }


    CONSTELLATION->nb_faces_theta_phi[itheta] = malloc(nb_phi * sizeof(double));
  if ( CONSTELLATION->nb_faces_theta_phi[itheta] == NULL){
    printf("***! Could not allow memory to CONSTELLATION->nb_faces_theta_phi[itheta]. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }



    CONSTELLATION->area_attitude_opengl_total[itheta] = malloc(nb_phi * sizeof(double));
  if ( CONSTELLATION->area_attitude_opengl_total[itheta] == NULL){
    printf("***! Could not allow memory to CONSTELLATION->area_attitude_opengl_total[itheta]. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }

  if (( OPTIONS->solar_cell_efficiency != -1) && (OPTIONS->opengl_power == 1)){
    CONSTELLATION->area_solar_panel_attitude_opengl[itheta] = malloc(nb_phi * sizeof(double));
  if ( CONSTELLATION->area_solar_panel_attitude_opengl[itheta] == NULL){
    printf("***! Could not allow memory to CONSTELLATION->area_solar_panel_attitude_opengl[itheta]. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }
  }

    for (iphi = 0; iphi < nb_phi; iphi++){
      CONSTELLATION->nb_faces_theta_phi[itheta][iphi]  = OPTIONS->nb_faces_theta_phi[itheta][iphi];
      CONSTELLATION->area_attitude_opengl[itheta][iphi] = malloc(OPTIONS->nb_faces * sizeof(double)); 
      CONSTELLATION->which_face_theta_phi[itheta][iphi] = malloc(OPTIONS->nb_faces * sizeof(int)); 
  CONSTELLATION->area_attitude_opengl_total[itheta][iphi] = OPTIONS->area_attitude_opengl_total[itheta][iphi] / (1e10);

  if (( OPTIONS->solar_cell_efficiency != -1) && (OPTIONS->opengl_power == 1)){
    CONSTELLATION->area_solar_panel_attitude_opengl[itheta][iphi] = OPTIONS->area_solar_panel_attitude_opengl[itheta][iphi]/ (1e10);
  }

      for (iface = 0; iface < OPTIONS->nb_faces_theta_phi[itheta][iphi]; iface++){
	iface_here = OPTIONS->which_face_theta_phi[itheta][iphi][iface];
	CONSTELLATION->area_attitude_opengl[itheta][iphi][iface_here] = OPTIONS->area_attitude_opengl[itheta][iphi][iface_here] / (1e10);  // OPTIONS->area_attitude_opengl is in cm^2. Convert it here to km^2
	CONSTELLATION->which_face_theta_phi[itheta][iphi][iface] = OPTIONS->which_face_theta_phi[itheta][iphi][iface] ;
      }
    }
  }

  


/*   for (itheta = 0; itheta < nb_theta; itheta++){ */
/*     for (iphi = 0; iphi < nb_phi; iphi++){ */
/*       printf("# %f %f %d\n", itheta * CONSTELLATION->area_attitude_opengl_dtheta + CONSTELLATION->area_attitude_opengl_theta0, iphi * CONSTELLATION->area_attitude_opengl_dphi + CONSTELLATION->area_attitude_opengl_phi0, CONSTELLATION->nb_faces_theta_phi[itheta][iphi]); */
/*       for (iface = 0; iface < CONSTELLATION->nb_faces_theta_phi[itheta][iphi]; iface++){ */
/* 	iface_here = CONSTELLATION->which_face_theta_phi[itheta][iphi][iface]; */
/* 	  printf("%d %f\n", iface_here, CONSTELLATION->area_attitude_opengl[itheta][iphi][iface_here]*1e10); */
/*       } */
/*       printf("%f\n", CONSTELLATION->area_attitude_opengl_total[itheta][iphi]*1e10); */
/*     } */
/*   } */
//  exitf();  



/*   for (itheta = 0; itheta < nb_theta; itheta++){ */
/*     for (iphi = 0; iphi < nb_phi; iphi++){ */
/*             printf("ddd # %f %f %f\n", itheta * OPTIONS->area_attitude_opengl_dtheta + OPTIONS->area_attitude_opengl_theta0, iphi * OPTIONS->area_attitude_opengl_dphi + OPTIONS->area_attitude_opengl_phi0, OPTIONS->area_solar_panel_attitude_opengl[itheta][iphi]); */
/*     } */
/*   } */
/*   exitf(); */

  }

	if (iDebugLevel >= 1){
	  if (iProc == 0) printf("-- (initialize_constellation) Done initializing the geometry.\n");
	}

	/*************************************************************************************/
	/*************************************************************************************/
	/********************* END OF COE AND TLE INITIALIZE: surface ************************/
	/*************************************************************************************/
	/*************************************************************************************/
	//      printf("iProc: %d | eee = %d | from %d to %d | ii = %d\n", iProc, eee, start_ensemble_bis + iProc * OPTIONS->nb_ensemble_min_per_proc, iProc * OPTIONS->nb_ensemble_min_per_proc + OPTIONS->nb_ensemble_min_per_proc, ii);


// the conversion of the acceleration from inrtl to lvlh is done in propagate_spacecraft but if init_flag = 1 the function propagate_spacecraft has not been called yet so the conversion inrtl to lvlh has not been done. Actually, the acceleration in the inertial frame has not been calculated yet either because it's calculated in propagate_spacecraft (with compute_dxdt)
      double v_dummy[3];
      double test_density;

      if (eee == 0){

	CONSTELLATION->spacecraft[ii][0].INTEGRATOR.file_given_output = fopen( CONSTELLATION->spacecraft[ii][0].INTEGRATOR.filename_given_output, "w+" );
      }
      CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.write_given_output = 1;
      //      printf("iProc %d start_ensemble[%d (%d)]: %d | dxdt\n", iProc, ii, eee, start_ensemble[ii]);

  // the only case where the attitute needs to be set is  if there is no ensemble at all run on the attitude
  if ( ( strcmp(OPTIONS->attitude_profile, "ensemble_angular_velocity") != 0 ) && ( strcmp(OPTIONS->attitude_profile, "ensemble_initial_attitude") != 0 ) ) { // if we do not run ensembles on the initial angular velocity

    if (OPTIONS->nb_ensembles_attitude <= 0){ // if we dont' run ensembles of any kind on the attitude

      if ( ( CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.solar_cell_efficiency != -1) || (GROUND_STATION->nb_ground_stations > 0) || ( CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_solar_pressure == 1) || ( CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_earth_pressure == 1)  || ( CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.include_drag == 1) ){ // these are the case where SpOCK uses the attitude

	if ( CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.isGPS == 0){ // no attitude to set up for GPS because we dont compute the drag, solar rariatin pressu,re, power, and gorund station coverage (attitude is needed for coverage because we calculate azimuth and levation angles form the spaceecraft reference system)
	  
	  set_attitude( CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.attitude,  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.index_in_attitude_interpolated, OPTIONS,  CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.file_is_quaternion);

	}
      }
    }
  }
  // end of if the only case where the attitute needs to be set is the only case where it has not been set initially in initialize_constellation, which is if there is no ensemble at all run on the attitude

      compute_dxdt( v_dummy, CONSTELLATION->spacecraft[ii][eee].a_i2cg_INRTL, &CONSTELLATION->spacecraft[ii][eee].et, CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL, PARAMS, &CONSTELLATION->spacecraft[ii][eee].INTEGRATOR, et_initial_epoch, OPTIONS->et_oldest_tle_epoch, &test_density, CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.index_in_attitude_interpolated, CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.index_in_driver_interpolated, CONSTELLATION, OPTIONS, iProc, iDebugLevel, &CONSTELLATION->spacecraft[ii][eee]);

      double T_inrtl_2_lvlh[3][3];
  compute_T_inrtl_2_lvlh(T_inrtl_2_lvlh, CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL);
  m_x_v(CONSTELLATION->spacecraft[ii][eee].a_i2cg_LVLH, T_inrtl_2_lvlh, CONSTELLATION->spacecraft[ii][eee].a_i2cg_INRTL);

	
      } // end of if main sc ii is run by this iProc OR  ( if this main sc is not run by this iProc and eee corresponds to an ensemble sc

 /* printf("%d\n", CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.nb_surfaces_eff); */
 /*  for (sss = 0; sss < CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.nb_surfaces_eff; sss++){ */
 /*    printf("normal[%d]: (%f, %f, %f) | total area: %f\n", sss, CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface_eff[sss].normal[0], CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface_eff[sss].normal[1], CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface_eff[sss].normal[2], CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.surface_eff[sss].area*1e10); */
 /*  } */
 /*  exitf(); */
 
      
    } // end go through all sc (including ensembels, if you run ensembles) other than GPS

    if (write_attitude_file == 1){
      fclose(file_attitude[ii]);
    }
    

/* // the conversion of the acceleration from inrtl to lvlh is done in propagate_spacecraft but if init_flag = 1 the function propagate_spacecraft has not been called yet so the conversion inrtl to lvlh has not been done. Actually, the acceleration in the inertial frame has not been calculated yet either because it's calculated in propagate_spacecraft (with compute_dxdt) */
/*       double v_dummy[3]; */
/*       double test_density; */

/*       if (eee == 0){ */

/* 	CONSTELLATION->spacecraft[ii][0].INTEGRATOR.file_given_output = fopen( CONSTELLATION->spacecraft[ii][0].INTEGRATOR.filename_given_output, "w+" ); */
/*       } */
/*       CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.write_given_output = 1; */
/*       printf("iProc %d start_ensemble[%d (%d)]: %d | dxdt\n", iProc, ii, eee, start_ensemble[ii]); */
/*       compute_dxdt( v_dummy, CONSTELLATION->spacecraft[ii][eee].a_i2cg_INRTL, &CONSTELLATION->spacecraft[ii][eee].et, CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL, PARAMS, &CONSTELLATION->spacecraft[ii][eee].INTEGRATOR, et_initial_epoch, OPTIONS->et_oldest_tle_epoch, &test_density, CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.index_in_attitude_interpolated, CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.index_in_driver_interpolated, CONSTELLATION, OPTIONS, iProc, iDebugLevel, &CONSTELLATION->spacecraft[ii][eee]); */

/*       double T_inrtl_2_lvlh[3][3]; */
/*   compute_T_inrtl_2_lvlh(T_inrtl_2_lvlh, CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL); */
/*   m_x_v(CONSTELLATION->spacecraft[ii][eee].a_i2cg_LVLH, T_inrtl_2_lvlh, CONSTELLATION->spacecraft[ii][eee].a_i2cg_INRTL); */

    //  !!!!!!!!!!!!! below can be computed and printed only with  one iproc
    if ((nProcs == 1) && (strcmp( OPTIONS->type_orbit_initialisation, "collision_vcm" ) == 0 )){
    double std_r_alongtrack, std_r_crosstrack, std_r_radial, std_v_alongtrack, std_v_crosstrack, std_v_radial, std_bc, std_srp;
    double mean_r_alongtrack, mean_r_crosstrack, mean_r_radial, mean_v_alongtrack, mean_v_crosstrack, mean_v_radial, mean_bc, mean_srp;
	     std_r_alongtrack =  gsl_stats_sd(r_alongtrack_all[ii], 1, ( OPTIONS->nb_ensembles_min  ) );
	     std_r_crosstrack =  gsl_stats_sd(r_crosstrack_all[ii], 1, ( OPTIONS->nb_ensembles_min  ) );
	     std_r_radial =  gsl_stats_sd(r_radial_all[ii], 1, ( OPTIONS->nb_ensembles_min  ) );

	     std_v_alongtrack =  gsl_stats_sd(v_alongtrack_all[ii], 1, ( OPTIONS->nb_ensembles_min  ) );
	     std_v_crosstrack =  gsl_stats_sd(v_crosstrack_all[ii], 1, ( OPTIONS->nb_ensembles_min  ) );
	     std_v_radial =  gsl_stats_sd(v_radial_all[ii], 1, ( OPTIONS->nb_ensembles_min  ) );

	     std_bc = gsl_stats_sd(bc_pert_all[ii], 1, ( OPTIONS->nb_ensembles_min  ) );
	     std_srp = gsl_stats_sd(srp_pert_all[ii], 1, ( OPTIONS->nb_ensembles_min  ) );

	     mean_r_alongtrack =  gsl_stats_mean(r_alongtrack_all[ii], 1, ( OPTIONS->nb_ensembles_min  ) );
	     mean_r_crosstrack =  gsl_stats_mean(r_crosstrack_all[ii], 1, ( OPTIONS->nb_ensembles_min  ) );
	     mean_r_radial =  gsl_stats_mean(r_radial_all[ii], 1, ( OPTIONS->nb_ensembles_min  ) );

	     mean_v_alongtrack =  gsl_stats_mean(v_alongtrack_all[ii], 1, ( OPTIONS->nb_ensembles_min  ) );
	     mean_v_crosstrack =  gsl_stats_mean(v_crosstrack_all[ii], 1, ( OPTIONS->nb_ensembles_min  ) );
	     mean_v_radial =  gsl_stats_mean(v_radial_all[ii], 1, ( OPTIONS->nb_ensembles_min  ) );

	     mean_bc = gsl_stats_mean(bc_pert_all[ii], 1, ( OPTIONS->nb_ensembles_min  ) );
	     mean_srp = gsl_stats_mean(srp_pert_all[ii], 1, ( OPTIONS->nb_ensembles_min  ) );

	     printf("\nSC %d by iProc %d\nSTD:\n",  ii, iProc);
	     printf("pos: %.4f (radial) %.4f (along) %.4f (cross) km\n", std_r_radial, std_r_alongtrack, std_r_crosstrack);
	     printf("vel: %.4e (radial) %.4e (along) %.4e (cross) km/s\n", std_v_radial, std_v_alongtrack, std_v_crosstrack );
	     //	     printf("bc: %.4e (%.4e) | srp %.4e (%.4e) m2/kg\n", std_bc * 1e6, OPTIONS->bc_vcm[ii] * sqrt( OPTIONS->covariance_matrix_equinoctial[ii][6][6] ) , std_srp * 1e6 ,  OPTIONS->srp_vcm[ii] * sqrt( OPTIONS->covariance_matrix_equinoctial[ii][8][8] ) );  //i noticed that if the normalized variance in the vcm for bc or srp if very small (ie the std is very small compared to the mean vvalue) then my randm algoirthm can't compute distibutino to reprouce this very small std. it will run fine but the std will tend to be larger (by a factor 10 for sitnace). it deosnt really matter for wha ti want to do because if the std on bc is super small, its exact value doesnt matter too much (even by a factor 10)

	     printf("bc: %.4e (%.4e) | srp %.4e (%.4e) m2/kg\n", std_bc*std_bc * 1e12,  OPTIONS->bc_cdm_std[ii]*OPTIONS->bc_cdm_std[ii]  , std_srp*std_bc * 1e12 ,  OPTIONS->srp_cdm_std[ii] * OPTIONS->srp_cdm_std[ii] );  

	     printf("\n MEAN:\n");
	     printf("pos: %.4e (radial) %.4e (along) %.4e (cross) m\n", mean_r_radial*1000, mean_r_alongtrack*1000, mean_r_crosstrack*1000);
	     printf("vel: %.4e (radial) %.4e (along) %.4e (cross) m/s\n", mean_v_radial*1000, mean_v_alongtrack*1000, mean_v_crosstrack*1000 );
	     //	     printf("bc: %.4e (%.4e) | srp %.4e (%.4e) m2/kg\n", mean_bc * 1e6, OPTIONS->bc_vcm[ii] , mean_srp * 1e6 ,  OPTIONS->srp_vcm[ii]  );
	     printf("bc: %.4e (%.4e) | srp %.4e (%.4e) m2/kg\n", mean_bc * 1e6, OPTIONS->bc_cdm[ii] , mean_srp * 1e6 ,  OPTIONS->srp_cdm[ii]  );
    }



    
  } // end of go through all main SC other than GPS


  
 
  /*************************************************************************************/
  /*************************************************************************************/
  /*************************************************************************************/
  /*************************************************************************************/
  /*************************************************************************************/
  /*************************************** GPS *****************************************/
  /*************************************************************************************/
  /*************************************************************************************/
  /*************************************************************************************/
  /*************************************************************************************/
  /*************************************************************************************/


  if (iDebugLevel >= 1){
    if (iProc == 0) printf("-- (initialize_constellation) Initializing the orbit for all GPS (%d GPS).\n", OPTIONS->nb_gps);
  }


  if (OPTIONS->nb_gps > 0){
    gps_tle_file = fopen(OPTIONS->tle_constellation_gps_filename,"r");
  }
  

  for (ii = OPTIONS->n_satellites - OPTIONS->nb_gps ; ii < OPTIONS->n_satellites; ii++){ // go over all GPS sc

    CONSTELLATION->spacecraft[ii] = malloc( sizeof(SPACECRAFT_T) ); // no ensemble for GPS satellites
    if ( start_ensemble[ii] == 0){ // if this iProc runs main sc ii  

      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.isGPS = 1;

      /*** Convert the TLEs into inertial state (postion, velocity) ***/
      SpiceInt frstyr = 1961; // do not really care about this line (as long as the spacecrafts are not flying after 2060)
      //      SpiceDouble geophs[8];
      int pp;
      int lineln=0;
      size_t len = 0;
      char *line_temp = NULL;
      //      SpiceDouble elems[10];
      SpiceDouble state[6];      


      if ( ii == iStart_save[iProc] ){
	for (bbb = 0; bbb < ii - (OPTIONS->n_satellites - OPTIONS->nb_gps); bbb++){ // for the first time the iProc opens the TLE file, it needs to skip the TLEs of the GPS that it is not running
	  getline( &line_temp, &len, gps_tle_file ); 
	  getline( &line_temp, &len, gps_tle_file );
	  getline( &line_temp, &len, gps_tle_file );
	}
      }

      /* /\* Set up the geophysical quantities.  At last check these were the values used by Space Command and SGP4 *\/ */
      /* PARAMS->geophs[ 0 ] =    1.082616e-3;   // J2 */
      /* PARAMS->geophs[ 1 ] =   -2.53881e-6;    // J3 */
      /* PARAMS->geophs[ 2 ] =   -1.65597e-6;    // J4 */
      /* PARAMS->geophs[ 3 ] =    7.43669161e-2; // KE */
      /* PARAMS->geophs[ 4 ] =    120.0;         // QO */
      /* PARAMS->geophs[ 5 ] =    78.0;          // SO */
      /* PARAMS->geophs[ 6 ] =    6378.135;      // ER */
      /* PARAMS->geophs[ 7 ] =    1.0;           // AE */

      /* Read in the next two lines from the text file that contains the TLEs. */


      // Skip the name of the GPS satellite that is being read
      getline( &line_temp, &len, gps_tle_file );
      //      printf("GPS %d by iProc %d:\n<%s>\n",ii-(OPTIONS->n_satellites - OPTIONS->nb_gps),iProc , line_temp);
      // First line
      getline( &line_temp, &len, gps_tle_file );
      lineln = strlen(line_temp)-1;
      SpiceChar line[2][lineln];
      for (pp = 0; pp<lineln-1 ; pp++)
	line[0][pp] = line_temp[pp];
      line[0][ lineln-1 ] = '\0';
      // Second line
      getline( &line_temp, &len, gps_tle_file );
      for (pp = 0; pp<lineln-1 ; pp++)
	line[1][pp] = line_temp[pp];
      line[1][ lineln-1 ] = '\0';

      // Convert the elements of the TLE into "elems" and "epoch" that can be then read by the SPICE routine ev2lin_ to convert into the inertial state
      getelm_c( frstyr, lineln, line, &OPTIONS->epoch_gps[ii-(OPTIONS->n_satellites - OPTIONS->nb_gps )], CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.elems );
      // Now propagate the state using ev2lin_ to the epoch of interest
      extern /* Subroutine */ int ev2lin_(SpiceDouble *, SpiceDouble *,
					  SpiceDouble *, SpiceDouble *);

      ev2lin_( &OPTIONS->epoch_gps[ii-(OPTIONS->n_satellites - OPTIONS->nb_gps )], PARAMS->geophs, CONSTELLATION->spacecraft[ii][eee].INTEGRATOR.elems, state );

            CONSTELLATION->spacecraft[ii][0].et = OPTIONS->epoch_gps[ii-(OPTIONS->n_satellites - OPTIONS->nb_gps )];
      CONSTELLATION->spacecraft[ii][eee].et_sc_initial = OPTIONS->epoch_gps[ii-(OPTIONS->n_satellites - OPTIONS->nb_gps )];

	    static doublereal  precm[36];
	        static doublereal invprc[36]	/* was [6][6] */;
		    extern /* Subroutine */ int zzteme_(doublereal *, doublereal *);
		        extern /* Subroutine */ int invstm_(doublereal *, doublereal *);
			    extern /* Subroutine */ int  mxvg_(
	    doublereal *, doublereal *, integer *, integer *, doublereal *);

	        zzteme_(&CONSTELLATION->spacecraft[ii][0].et, precm);

/*     ...now convert STATE to J2000. Invert the state transformation */
/*     operator (important to correctly do this). */

    invstm_(precm, invprc);
    static integer c__6 = 6;
        static doublereal tmpsta[6];
	    mxvg_(invprc, state, &c__6, &c__6, tmpsta);
    moved_(tmpsta, &c__6, state);


      for (pp = 0; pp<3; pp++){
	CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL[pp] = state[pp];
	CONSTELLATION->spacecraft[ii][0].v_i2cg_INRTL[pp] = state[pp+3];
      }

      // Initialize the Keplerian Elements
      ORBITAL_ELEMENTS_T OE_temp;
      cart2kep( &OE_temp, CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][0].v_i2cg_INRTL, CONSTELLATION->spacecraft[ii][0].et,  PARAMS->EARTH.GRAVITY.mu);

	    CONSTELLATION->spacecraft[ii][eee].OE.initial_an_to_sc = OE_temp.an_to_sc;
	    CONSTELLATION->spacecraft[ii][eee].OE.w_ave_temp = OE_temp.w;
	    CONSTELLATION->spacecraft[ii][eee].OE.w_ave = 9999.999999/RAD2DEG;
	    CONSTELLATION->spacecraft[ii][eee].OE.sma_ave_temp = OE_temp.sma;
	    CONSTELLATION->spacecraft[ii][eee].OE.sma_ave = 9999.999999;
	    CONSTELLATION->spacecraft[ii][eee].OE.ecc_ave_temp = OE_temp.eccentricity;
	    CONSTELLATION->spacecraft[ii][eee].OE.ecc_ave = 9999.999999;

	    CONSTELLATION->spacecraft[ii][eee].OE.ave_increm = 1;
	    CONSTELLATION->spacecraft[ii][eee].et_last_orbit = CONSTELLATION->spacecraft[ii][eee].et;
	    CONSTELLATION->spacecraft[ii][eee].orbit_number = 0;


      CONSTELLATION->spacecraft[ii][0].OE.sma          =OE_temp.sma;
      CONSTELLATION->spacecraft[ii][0].OE.eccentricity =OE_temp.eccentricity;
      CONSTELLATION->spacecraft[ii][0].OE.inclination  =OE_temp.inclination;
      CONSTELLATION->spacecraft[ii][0].OE.long_an      =OE_temp.long_an;
      CONSTELLATION->spacecraft[ii][0].OE.w            =OE_temp.w;
      CONSTELLATION->spacecraft[ii][0].OE.f            =OE_temp.f;
      CONSTELLATION->spacecraft[ii][0].OE.tp           =OE_temp.tp;
      CONSTELLATION->spacecraft[ii][0].OE.E            =OE_temp.E;
      CONSTELLATION->spacecraft[ii][0].OE.ra            =OE_temp.ra;

      radius_perigee = CONSTELLATION->spacecraft[ii][0].OE.sma * ( 1 - CONSTELLATION->spacecraft[ii][0].OE.eccentricity );
      if (radius_perigee < PARAMS->EARTH.radius){
	printf("***! The orbit of satellite %d intersects the Earth (altitude of perigee = %f km). The program will stop. !***\n",ii, radius_perigee-PARAMS->EARTH.radius);
	MPI_Finalize();
	exit(0);
      }

      //// !!!!!! WEIRD: In this block, I test if going back to cart then to kep gives the same OE. It does. BUT, the cart is different (meaning the next line finds different cart than the one calculated a bit above but still results in the same OE). The reason for that is still not clear and if has to be found in the future! Future tests will validate that these results are correct, though.
      // Initialize the inertial state
      /*       kep2cart(   CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][0].v_i2cg_INRTL, */
      /* 		  &PARAMS->EARTH.GRAVITY.mu, */
      /* 		  &CONSTELLATION->spacecraft[ii][0].OE); // Computes the ECI coordinates based on the Keplerian inputs (orbital elements and mu) (in propagator.c) (CBV) */


      /*       ORBITAL_ELEMENTS_T OE_temp2; */
      /*       cart2kep( &OE_temp2, CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL, CONSTELLATION->spacecraft[ii][0].v_i2cg_INRTL, CONSTELLATION->spacecraft[ii][0].et ,  PARAMS->EARTH.GRAVITY.mu); */

      /*       printf("\nsma = %f | %f\n",CONSTELLATION->spacecraft[ii][0].OE.sma, OE_temp2.sma); */
      /*       printf("inclination = %f | %f\n",CONSTELLATION->spacecraft[ii][0].OE.inclination*RAD2DEG, OE_temp2.inclination*RAD2DEG); */
      /*       printf("eccentricity = %f | %f\n", CONSTELLATION->spacecraft[ii][0].OE.eccentricity, OE_temp2.eccentricity); */
      /*       printf("long_an = %f | %f\n", CONSTELLATION->spacecraft[ii][0].OE.long_an, OE_temp2.long_an); */
      /*       printf("w = %f | %f\n", CONSTELLATION->spacecraft[ii][0].OE.w*RAD2DEG, OE_temp2.w*RAD2DEG); */
      /*       printf("tp = %f | %f\n", CONSTELLATION->spacecraft[ii][0].OE.tp, OE_temp2.tp); */
      /*       printf("E = %f | %f\n", CONSTELLATION->spacecraft[ii][0].OE.E, OE_temp2.E); */


/* 	eci2lla(CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL , CONSTELLATION->spacecraft[ii][eee].et, geodetic ); */
	
/* 	CONSTELLATION->spacecraft[ii][0].GEODETIC.altitude = geodetic[2]; */
/* 	CONSTELLATION->spacecraft[ii][0].GEODETIC.latitude = geodetic[0]; */
/* 	CONSTELLATION->spacecraft[ii][0].GEODETIC.longitude = geodetic[1];  */

/* 	// Initialize the planet fixed state		 */

/* 	  geodetic_to_geocentric(PARAMS->EARTH.flattening,             */
/* 				 CONSTELLATION->spacecraft[ii][0].GEODETIC.altitude, */
/* 				 CONSTELLATION->spacecraft[ii][0].GEODETIC.latitude, */
/* 				 CONSTELLATION->spacecraft[ii][0].GEODETIC.longitude, */
/* 				 PARAMS->EARTH.radius,        */
/* 				 CONSTELLATION->spacecraft[ii][0].r_ecef2cg_ECEF) ; */


      // Initialize the planet fixed state
      estate[0] = CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL[0];estate[1] = CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL[1];estate[2] = CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL[2];
      estate[3] = CONSTELLATION->spacecraft[ii][0].v_i2cg_INRTL[0];estate[4] = CONSTELLATION->spacecraft[ii][0].v_i2cg_INRTL[1];estate[5] = CONSTELLATION->spacecraft[ii][0].v_i2cg_INRTL[2];
      sxform_c (  "J2000", PARAMS->EARTH.earth_fixed_frame,  CONSTELLATION->spacecraft[ii][0].et,    xform  );
      mxvg_c   (  xform,       estate,   6,  6, jstate );
      CONSTELLATION->spacecraft[ii][0].r_ecef2cg_ECEF[0] = jstate[0]; CONSTELLATION->spacecraft[ii][0].r_ecef2cg_ECEF[1] = jstate[1]; CONSTELLATION->spacecraft[ii][0].r_ecef2cg_ECEF[2] = jstate[2];
      CONSTELLATION->spacecraft[ii][0].v_ecef2cg_ECEF[0] = jstate[3]; CONSTELLATION->spacecraft[ii][0].v_ecef2cg_ECEF[1] = jstate[4]; CONSTELLATION->spacecraft[ii][0].v_ecef2cg_ECEF[2] = jstate[5];
      /* pxform_c( "J2000", PARAMS->EARTH.earth_fixed_frame, CONSTELLATION->spacecraft[ii][0].et, T_J2000_to_ECEF); // Return the matrix (here T_J2000_to_ECEF) that transforms position vectors from one specified frame (here J2000) to another (here ITRF93) at a specified epoch  (in /Users/cbv/cspice/src/cspice/pxform_c.c) (CBV) */
      /* m_x_v(CONSTELLATION->spacecraft[ii][0].r_ecef2cg_ECEF, T_J2000_to_ECEF, CONSTELLATION->spacecraft[ii][0].r_i2cg_INRTL); // multiply a matrix by a vector (in prop_math.c). So here we convert ECI (J2000) to ECEF coordinates (CBV) */
        
      // initialize the geodetic state-
      geocentric_to_geodetic(
			     CONSTELLATION->spacecraft[ii][0].r_ecef2cg_ECEF,
			     &PARAMS->EARTH.radius,
			     &PARAMS->EARTH.flattening,
			     &CONSTELLATION->spacecraft[ii][0].GEODETIC.altitude,
			     &CONSTELLATION->spacecraft[ii][0].GEODETIC.latitude,
			     &CONSTELLATION->spacecraft[ii][0].GEODETIC.longitude ); // Computes lat/long/altitude based on ECEF coordinates and planet fixed state (planetary semimajor axis, flattening parameter) (in propagator.c) (CBV)
        

      // Output file for each spacecraft
      strcpy(CONSTELLATION->spacecraft[ii][0].filename, OPTIONS->dir_output_run_name_sat_name[ii]);
      strcat(CONSTELLATION->spacecraft[ii][0].filename, "/");
      strcat(CONSTELLATION->spacecraft[ii][0].filename, OPTIONS->filename_output[ii]);
    
      // Now CONSTELLATION->spacecraft[ii][0].filenameecef is moved to generate_ephemerides because we want iProc 0 to know CONSTELLATION->spacecraft[ii][0].filenameecef even for the main sc ii that it does not run (this is because iProc 0 will gather all ECEF files at the end of the propagation)
      /*     strcpy(CONSTELLATION->spacecraft[ii][0].filenameecef, OPTIONS->dir_output_run_name_sat_name[ii]); */
      /*     strcat(CONSTELLATION->spacecraft[ii][0].filenameecef, "/"); */
      /*     strcat(CONSTELLATION->spacecraft[ii][0].filenameecef, "ECEF_"); */
      /*     strcat(CONSTELLATION->spacecraft[ii][0].filenameecef, OPTIONS->filename_output[ii]); */

      strcpy(CONSTELLATION->spacecraft[ii][0].filenamerho, OPTIONS->dir_output_run_name_sat_name[ii]);
      strcat(CONSTELLATION->spacecraft[ii][0].filenamerho, "/");
      strcat(CONSTELLATION->spacecraft[ii][0].filenamerho, "density_");
      strcat(CONSTELLATION->spacecraft[ii][0].filenamerho,OPTIONS->filename_output[ii]);


      strcpy(CONSTELLATION->spacecraft[ii][0].filenameout, OPTIONS->dir_output_run_name_sat_name[ii]);
      strcat(CONSTELLATION->spacecraft[ii][0].filenameout, "/");
      strcat(CONSTELLATION->spacecraft[ii][0].filenameout, "LLA_");
      strcat(CONSTELLATION->spacecraft[ii][0].filenameout,OPTIONS->filename_output[ii]);

            strcpy(CONSTELLATION->spacecraft[ii][0].filenameatt, OPTIONS->dir_output_run_name_sat_name[ii]);
      strcat(CONSTELLATION->spacecraft[ii][0].filenameatt, "/");
      strcat(CONSTELLATION->spacecraft[ii][0].filenameatt, "LLA_");
      strcat(CONSTELLATION->spacecraft[ii][0].filenameatt,OPTIONS->filename_output[ii]);

      strcpy(CONSTELLATION->spacecraft[ii][0].filenamekalman, OPTIONS->dir_output_run_name_sat_name[ii]);
      strcat(CONSTELLATION->spacecraft[ii][0].filenamekalman, "/");
      strcat(CONSTELLATION->spacecraft[ii][0].filenamekalman, "kalman_");
      strcat(CONSTELLATION->spacecraft[ii][0].filenamekalman,OPTIONS->filename_output[ii]);

      strcpy(CONSTELLATION->spacecraft[ii][0].filenamekalman_meas, OPTIONS->dir_output_run_name_sat_name[ii]);
      strcat(CONSTELLATION->spacecraft[ii][0].filenamekalman_meas, "/");
      strcat(CONSTELLATION->spacecraft[ii][0].filenamekalman_meas, "meas_converted_kalman_");
      strcat(CONSTELLATION->spacecraft[ii][0].filenamekalman_meas,OPTIONS->filename_output[ii]);
      
      strcpy(CONSTELLATION->spacecraft[ii][0].filename_kalman_init, OPTIONS->filename_kalman_init);
      // I COMMENTED BELOW EBCAUSE I DON'T SEE ANY REASON TO COMPUTE THE GROUND STATION COVERAGE FOR GPS 
      /*     if (OPTIONS->nb_ground_stations > 0){ */
      /*       for ( iground = 0; iground < OPTIONS->nb_ground_stations; iground ++){ */
      /* 	strcpy(CONSTELLATION->spacecraft[ii][0].filename_coverage_ground_station[iground], OPTIONS->dir_output_run_name_sat_name_coverage[ii]); */
      /* 	strcat(CONSTELLATION->spacecraft[ii][0].filename_coverage_ground_station[iground], "/"); */
      /* 	strcat(CONSTELLATION->spacecraft[ii][0].filename_coverage_ground_station[iground], OPTIONS->name_ground_station[iground]); */
      /* 	strcat(CONSTELLATION->spacecraft[ii][0].filename_coverage_ground_station[iground], "_by_"); */
      /* 	strcat(CONSTELLATION->spacecraft[ii][0].filename_coverage_ground_station[iground],OPTIONS->filename_output[ii]); */
      /*       } */
      /*     } */


      // cbv COMMENTED THIS BLOCK ON APRIL 22 218 BECASESUE SOLAR POWER IS NOT COMPUTED FOR GPS
/*       if (CONSTELLATION->spacecraft[ii][0].INTEGRATOR.solar_cell_efficiency != -1){ */
/* 	strcpy(CONSTELLATION->spacecraft[ii][0].filenamepower, OPTIONS->dir_output_run_name_sat_name[ii]); */
/* 	strcat(CONSTELLATION->spacecraft[ii][0].filenamepower, "/"); */
/* 	strcat(CONSTELLATION->spacecraft[ii][0].filenamepower, "power_"); */
/* 	strcat(CONSTELLATION->spacecraft[ii][0].filenamepower,OPTIONS->filename_output[ii]); */

/* 	strcpy(CONSTELLATION->spacecraft[ii][0].filenameeclipse, OPTIONS->dir_output_run_name_sat_name[ii]); */
/* 	strcat(CONSTELLATION->spacecraft[ii][0].filenameeclipse, "/"); */
/* 	strcat(CONSTELLATION->spacecraft[ii][0].filenameeclipse, "eclipse_"); */
/* 	strcat(CONSTELLATION->spacecraft[ii][0].filenameeclipse,OPTIONS->filename_output[ii]); */
      
/*       } */
      // end of cbv COMMENTED THIS BLOCK ON APRIL 22 218 BECASESUE SOLAR POWER IS NOT COMPUTED FOR GPS

      // Name of each satellite
      strcpy(CONSTELLATION->spacecraft[ii][0].name_sat, "");
      strcpy(CONSTELLATION->spacecraft[ii][0].name_sat, OPTIONS->gps_file_name[ii - (OPTIONS->n_satellites - OPTIONS->nb_gps)]);

    
      /*   // Compute the Orbital period of the first state */
      /*       period = pow( CONSTELLATION->spacecraft[ii][0].OE.sma, 3.0); */
      /*       period = period / PARAMS->EARTH.GRAVITY.mu; */
      /*       period = 2.0 * M_PI * sqrt( period ); */
    
      // Integrator
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.sc_main_nb = ii;
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.sc_ensemble_nb = 0;

      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.dt            = OPTIONS->dt;
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.dt_pos_neg            = OPTIONS->dt_pos_neg;
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.nb_surfaces   = 1;
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.nb_surfaces_eff   = 1;
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.mass          = 1630.0; // Wikipedia...
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.solar_cell_efficiency = -1; // do not care for now
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.degree        = OPTIONS->degree;       // Gravity degree
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.order         = OPTIONS->order;        // Gravity order
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.include_drag  = 0;
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.include_solar_pressure  = 0; // do not care for now
            CONSTELLATION->spacecraft[ii][0].INTEGRATOR.include_earth_pressure  = 0; // do not care for now
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.include_sun   = OPTIONS->include_sun;  // include Sun perturbations (CBV)
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.include_moon  = OPTIONS->include_moon; // include Moon perturbations (CBV)
      
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.index_in_attitude_interpolated = 0; // do not care for now because the attitude is for the power, solar radiation pressure, earth_pressure, or drag. For GPS satellite, we do not compute these (the main reason is that the TLE do not give the geometry). However, we need to set it to 0 (like a regular spacecraft) because to compute coverage, the variable index_in_attitude_interpolated is used
      CONSTELLATION->spacecraft[ii][0].INTEGRATOR.index_in_driver_interpolated = 0; // do not care for now because the attitude is for the power, solar radiation pressure, or drag. For GPS satellite, we do not compute these (the main reason is that the TLE do not give the geometry)

	CONSTELLATION->spacecraft[ii][0].INTEGRATOR.density_mod = OPTIONS->density_mod; // // the desnity given by msis is multiplied by density_mod + density_amp * sin(2*pi*t/T + density_phase*T) where T is the orbital period  
	CONSTELLATION->spacecraft[ii][0].INTEGRATOR.density_mod_amp = OPTIONS->density_mod_amp; 
	CONSTELLATION->spacecraft[ii][0].INTEGRATOR.density_mod_phase = OPTIONS->density_mod_phase; 

    } // end of if this iProc runs main sc ii  
  } // end of go over all GPS sc

  if (OPTIONS->nb_gps > 0){
    fclose(gps_tle_file);
  }
  if (iDebugLevel >= 1){
    if (iProc == 0) printf("-- (initialize_constellation) Done initializing the orbit for all GPS (%d GPS).\n", OPTIONS->nb_gps);
  }
  
  /* For the specular points computation */
  if (OPTIONS->nb_gps > 0){

    strcpy(CONSTELLATION->filename_CYGNSS_constellation_for_specular, OPTIONS->dir_output_run_name);
    strcat(CONSTELLATION->filename_CYGNSS_constellation_for_specular, "/CONSTELLATION_CYGNSS_for_run_");
    next = OPTIONS->filename_output[0];
    find_file_name =  (int)(strchr(next, '.') - next);
    strncat(CONSTELLATION->filename_CYGNSS_constellation_for_specular, next, find_file_name);
    strcat(CONSTELLATION->filename_CYGNSS_constellation_for_specular, ".txt");

    strcpy(CONSTELLATION->filename_GPS_constellation_for_specular, OPTIONS->dir_output_run_name);
    strcat(CONSTELLATION->filename_GPS_constellation_for_specular, "/CONSTELLATION_GPS_for_run_");
    next = OPTIONS->filename_output[0];
    find_file_name =  (int)(strchr(next, '.') - next);
    strncat(CONSTELLATION->filename_GPS_constellation_for_specular, next, find_file_name);
    strcat(CONSTELLATION->filename_GPS_constellation_for_specular, ".txt");

  }
  /* end of for the specular points computation */



  if (OPTIONS->nb_ground_stations > 0){
    for (iground = 0; iground < OPTIONS->nb_ground_stations; iground++){
      strcpy( GROUND_STATION->name_ground_station[iground], OPTIONS->name_ground_station[iground]);
      GROUND_STATION->latitude_ground_station[iground] = OPTIONS->latitude_ground_station[iground] * DEG2RAD;
      GROUND_STATION->longitude_ground_station[iground] = OPTIONS->longitude_ground_station[iground] * DEG2RAD;
      GROUND_STATION->altitude_ground_station[iground] = OPTIONS->altitude_ground_station[iground];
      GROUND_STATION->min_elevation_angle_ground_station[iground] = OPTIONS->min_elevation_angle_ground_station[iground] * DEG2RAD;
      // Convert lla position of ground station into ECEF
      geodetic_to_geocentric( PARAMS->EARTH.flattening, GROUND_STATION->altitude_ground_station[iground]/1000., GROUND_STATION->latitude_ground_station[iground], GROUND_STATION->longitude_ground_station[iground], PARAMS->EARTH.radius, GROUND_STATION->ecef_ground_station[iground]);

    }
    GROUND_STATION->nb_ground_stations = OPTIONS->nb_ground_stations;
  }

  /* if ( ( strcmp(OPTIONS->format_density_driver, "density_file") != 0 ) && ( strcmp(OPTIONS->format_density_driver, "gitm") != 0 ) ){ */
  /*   if ( strcmp(OPTIONS->format_density_driver, "dynamic") == 0 ){  */

  /*     /\* free(OPTIONS->Ap); *\/ */
  /*     /\* //	    if (OPTIONS->use_ap_hist == 1){ *\/ */
  /*     /\* free(OPTIONS->Ap_hist); *\/ */
  /*     /\* //	    } *\/ */
  /*     /\* free(OPTIONS->f107); *\/ */
  /*     /\* free(OPTIONS->f107A); *\/ */

  /*   } */
  /* } */
  /* else{ */
  /*   free(OPTIONS->density); */
  /* } */


  /* free(OPTIONS->et_interpo); */
  //    printf("uudpqw90jdwqdpwuuuuuu %f\n", CONSTELLATION->spacecraft[0][0].INTEGRATOR.attitude.pitch[0]);

  return 0;


}



double randu ( double mu, double sigma) // source: http://phoxis.org/2013/05/04/plot-histogram-in-terminal/
{
  double U1;
  struct timeval t1;
  gettimeofday(&t1, NULL);
  srand(t1.tv_usec * t1.tv_sec);

  U1 = -1 + ((double) rand () / RAND_MAX) * 2;
 
  return (mu + sigma * (double) U1);
}



double randn ( double mu, double sigma) // source: http://phoxis.org/2013/05/04/plot-histogram-in-terminal/
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
  struct timeval t1;
  gettimeofday(&t1, NULL);
  srand(t1.tv_usec * t1.tv_sec);
  //	  srand(time(NULL) + iProc); 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}


double randn_iproc (int iProc, double mu, double sigma) // source: http://phoxis.org/2013/05/04/plot-histogram-in-terminal/
// note: the only difference with rand is that the seed here depends on iProc (the seed should actually also depends on iProcs in rand because it depends on the time with microsecond accuracy so each iProcs calls it at a different time)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
  struct timeval t1;
  gettimeofday(&t1, NULL);
  srand(t1.tv_usec * t1.tv_sec + iProc);
  //	  srand(time(NULL) + iProc); 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}


int sign(double x){
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}



/*  LocalWords:  min
 */
