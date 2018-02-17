////////////////////////////////////////////////////////////////////////////////////////
//
//  Michigan Orbit Analyst Tool
//
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | J. Getchius   | 05/09/2015    |   ---     | Initial Implementation
//
//////////////////////////////////////////////////////////////////////////////////////////

// Includes
#include "options.h"
#include "propagator.h" 
#include "moat_prototype.h"
#include "kalman_9state_new_with_tau.h" // kalman


int nProcs;
int iProc;

int main(int argc, char * argv[]) {


  
    // Declarations

    OPTIONS_T        OPTIONS;
    PARAMS_T         PARAMS;
    GROUND_STATION_T GROUND_STATION;
    //    STORM_T         STORM;
    char             filename_input[256];
    char filename_input_no_path[256];
    int ierr;
    CONSTELLATION_T *CONSTELLATION = malloc(sizeof(CONSTELLATION_T));

    ierr = MPI_Init(&argc, &argv);
     
    /* find out MY process ID, and how many processes were started. */
      
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &iProc);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nProcs);


    
/*     char test_wget[1000]; */
/*     strcpy(test_wget, "wget http://services.swpc.noaa.gov/text/45-day-ap-forecast.txt"); */
/*     system(test_wget); */
    
//    printf("enwofjoiewjfoiewjfew\neowifheiowjfewewj\nfeowfjhweoinfeiw\n");
    strcpy( filename_input, "./input/main_input/");
    strcat(filename_input, argv[1]);
    strcpy(filename_input_no_path, "");
    strcat(filename_input_no_path, argv[1]);



    // Debug Level
    int iDebugLevel = -1;
    if (argc > 2){
    iDebugLevel = atoi(argv[2]);
    }
    //  Load Options
    if (iDebugLevel >= 0){
      if (iProc == 0) printf("\n- Loading options... (iProc %d)\n", iProc);
    }
    //    OPTIONS.new_cd = 1; // set it to 1 to iinitalize the drag coeff with the accomoation coefficient in the geometry file (and not the drag coefficient). Does not work if running ensembles on Cd

    load_options( &OPTIONS, filename_input, iProc, nProcs, iDebugLevel, filename_input_no_path);
    /*      print_options(&OPTIONS); */
/*     int pp; */
/*     for (pp = 0; pp < OPTIONS.nb_time_steps*2 -1; pp++){ */
/*       etprint(OPTIONS.et_interpo[pp], ""); */
/*     } */
    //           MPI_Finalize();exit(0);
    // Set up SPICE
    if (iDebugLevel >= 0){
     if (iProc == 0) printf("\n- Loading spice... (iProc %d)\n", iProc);
    }
    if ((iProc == 0) && (iDebugLevel >= 1)){
      printf("-- (main) Spice eop: %s\n", OPTIONS.eop_file);
      printf("-- (main) Spice planet_ephem_: %s\n", OPTIONS.planet_ephem_file);
      printf("-- (main) Spice earth_binary: %s\n", OPTIONS.earth_binary_pck);
    }
    furnsh_c(OPTIONS.eop_file); // EOP parameters file
    furnsh_c(OPTIONS.planet_ephem_file); // Load planetary ephemerides
    furnsh_c(OPTIONS.earth_binary_pck);

    // Load Params
    if (iDebugLevel >= 0){
      if (iProc == 0) printf("\n- Loading parameters... (iProc %d)\n", iProc);
    }


    //newstructure
    //load_params( &PARAMS, main_directory_location, iDebugLevel, OPTIONS.earth_fixed_frame, OPTIONS.use_ap_hist, iProc );
    load_params( &PARAMS,  iDebugLevel, OPTIONS.earth_fixed_frame, OPTIONS.use_ap_hist, iProc, OPTIONS.path_to_spice );
    //newstructure


/*   // Update ECEF state */
/*   estate[0] = ; estate[1] = ; estate[2] = ; */
/*   estate[3] = ; estate[4] = ; estate[5] = ; */
/*   sxform_c (  "J2000", PARAMS->EARTH.earth_fixed_frame,  ,    xform  );  */

/*   mxvg_c   (  xform,       estate,   6,  6, jstate );  */
/*   SC->r_ecef2cg_ECEF[0] = jstate[0]; SC->r_ecef2cg_ECEF[1] = jstate[1]; SC->r_ecef2cg_ECEF[2] = jstate[2]; */
/*   SC->v_ecef2cg_ECEF[0] = jstate[3]; SC->v_ecef2cg_ECEF[1] = jstate[4]; SC->v_ecef2cg_ECEF[2] = jstate[5]; */



/*     MPI_Finalize();exit(0); */





/* // !!!!!!!!!! TEMPORARY - quaternion -> ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/q2m_c.html */
    



/*     SpiceDouble  q[4]; */
/*       SpiceDouble       r[3][3]; */
/*       // if q given as SPICE representation (ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/q2m_c.html) */
/*       q[0] = -1.45417517048816e-001; */
/*       q[1] = 8.73771916300652e-001; */
/*       q[2] =  -4.61631150177178e-001; */
/*       q[3] =  4.76766736018681e-002; */
/*       q[0] = sqrt(2)/2; */
/*       q[1] = 0; */
/*       q[2] = 0; */
/*       q[3] = -sqrt(2)/2; */


/* /\* /\\*       // if q given as engineering representation then convert to SPICE representation (ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/q2m_c.html) *\\/ *\/ */
/* /\*       q[0] = q[3]; *\/ */
/* /\*       q[1] = -q[0]; *\/ */
/* /\*       q[2] = -q[1]; *\/ */
/* /\*       q[3] = -q[2]; *\/ */


/* /\*       q2m_c(q,r); *\/ */
/* /\*       SpiceDouble tau, delta, alpha; *\/ */
/* /\*       m2eul_c ( r, 1, 2, 3, &tau, &delta, &alpha ); *\/ */


/* /\*       int i,j; *\/ */
/* /\*       for (i = 0; i < 3; i++){ *\/ */
/* /\* 	for (j = 0; j<3; j++){ *\/ */
/* /\* 	  printf("%f ", r[i][j]); *\/ */
/* /\* 	} *\/ */
/* /\* 	printf("\n"); *\/ */
/* /\*       } *\/ */
/* /\*       printf("tau %f\ndelta %f\nalpha %f\n", tau*180./M_PI, delta*180./M_PI, alpha*180./M_PI); *\/ */

/* /\*            MPI_Finalize();exit(0); *\/ */
/* // !!!!!!!!!! end of TEMPORARY - quaternion */





    //  Build Constellation
    if (iDebugLevel >= 0){
     
      if (iProc == 0) printf("\n- Initializing Constellation... (iProc: %d)\n", iProc);
    }


    initialize_constellation( CONSTELLATION, &OPTIONS, &PARAMS, &GROUND_STATION, iDebugLevel, iProc, nProcs);

    //    printf("\nSSSSSSSSSSSS %d\n",OPTIONS.nb_time_steps);

     if ( OPTIONS.use_kalman == 1 ){  // if 0 then it uses the classical propagation like in the previous veresions. If set to 1 then it uses Kalman Filter, which means observations are needed as inputs
     // Kalman filter
     MEAS_T MEAS;
     KALMAN_T KF;

     char filename_meas[1000];
  char *line = NULL;
  size_t len = 0;

     KF.sc = CONSTELLATION->spacecraft[0][0]; // !!!!!!! all sc

     //     strcpy(filename_meas, "/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock/spock/true_10s/true_10s1/noise_true_10s1.txt");// "netcdf/output/cyg07.ddmi.s20170609-000000-e20170609-235959.l1.power-brcs.sand004.txt")  ;
     KF.fp_kalman_init = fopen(CONSTELLATION->spacecraft[0][0].filename_kalman_init, "r");

       getline(&line, &len, KF.fp_kalman_init);
       strtok(line, "\n");  strtok(line, "\r");

       strcpy(filename_meas, "");
       sscanf(line, "%s", filename_meas);
       getline(&line, &len, KF.fp_kalman_init);

       //          strcpy(filename_meas, "/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock/spock/true/true1/noise_true1.txt");// "netcdf/output/cyg07.ddmi.s20170609-000000-e20170609-235959.l1.power-brcs.sand004.txt")  ;
  	  //	  strcpy(KF.filename_kalman_init, CONSTELLATION->spacecraft[0][0].filename_kalman_init) ;// !!!!!! all sc
     MEAS.fp_meas = fopen(filename_meas, "r");

     str2et_c(OPTIONS.initial_epoch, &OPTIONS.et_initial_epoch);

  str2et_c(OPTIONS.final_epoch, &OPTIONS.et_final_epoch);
  double min_end_time;
  min_end_time = OPTIONS.et_final_epoch;
  //  int tt=0;
  while ((feof(MEAS.fp_meas) == 0) && (KF.sc.et  <= OPTIONS.et_final_epoch)){// !!!!!! all sc
       // Time to stop the propgation
       
       //    etprint(min_end_time, "end");


       // Print progress to screen
       if (iProc == 0){
	 	 print_progress_kalman( min_end_time, KF.sc.et , OPTIONS.et_initial_epoch, iProc, OPTIONS.nb_gps )  ;
       }

       kalman_filt( &MEAS, &KF, &PARAMS, &OPTIONS, &GROUND_STATION, CONSTELLATION, iProc, iDebugLevel, nProcs);


      
      KF.sc.et_next_time_step = OPTIONS.et_initial_epoch + (int)( ( MEAS.et - OPTIONS.et_initial_epoch ) / OPTIONS.dt ) * OPTIONS.dt + OPTIONS.dt;
/*         etprint(KF.sc.et_next_time_step, "KF.sc.et_next_time_step"); */
/* 	etprint(MEAS.et, "meas"); */

//      etprint(KF.sc.et, "after meas");
      //         MPI_Finalize();exit(0);

       //       tt = tt+1;
/*        if (tt == 2){ */
/*   MPI_Finalize();exit(0); */
/*        } */

  if (MEAS.et > OPTIONS.et_final_epoch){
    break;
  }  

/*        etprint(MEAS.et, "eee"); */
/*        etprint(KF.sc.et, "KF.sc.et"); */
/*        etprint(OPTIONS.et_final_epoch, "OPTIONS.et_final_epoch"); */
     }
     fclose(KF.fp_kalman_init);
     fclose(KF.fp_kalman);
     fclose(KF.fp_kalman_meas);
     fclose(MEAS.fp_meas);
         	   if (iProc == 0) {
    	printf("\n- Done running the Kalman filter.\n");
    	        }
/* 		   // Also propagate without kalman filter (to compare kalman vs no kalman) */
/*     if (iDebugLevel >= 0){ */
/*       if (iProc == 0){ */
/*       printf("\n- Generating the ephemerides... (iProc: %d)\n", iProc); */
/*     } */
/*     } */
/*        generate_ephemerides( CONSTELLATION, &OPTIONS, &PARAMS, &GROUND_STATION, iProc, nProcs, iDebugLevel); */

/*     // Notify Exit */
/*     	   if (iProc == 0) { */
/*     	printf("\n- Done propagating the spacecraft.\n"); */
/*     	        } */

     }
     else{
       // Generate the Ephemerides
    if (iDebugLevel >= 0){
      if (iProc == 0){
      printf("\n- Generating the ephemerides... (iProc: %d)\n", iProc);
    }
    }

       generate_ephemerides( CONSTELLATION, &OPTIONS, &PARAMS, &GROUND_STATION, iProc, nProcs, iDebugLevel);

    // Notify Exit
    	   if (iProc == 0) {
    	printf("Done propagating the spacecraft.\n");
    	        }
     }


	    MPI_Barrier(MPI_COMM_WORLD);
    free(CONSTELLATION);
    
        ierr = MPI_Finalize();

    return 0;
}
