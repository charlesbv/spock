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

#include <mpi.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <unistd.h>
#include "options.h"
#include "propagator.h"
#include <string.h>
#include <stdlib.h>
 #include <stdio.h>
 #include "SpiceUsr.h"
 #include <unistd.h>
 #include <sys/types.h>
 #include <pwd.h>
 #include "time.h"
 #include "dirent.h"
 #include "gsl/gsl_math.h"
 #include "gsl/gsl_eigen.h"
 #include "gsl/gsl_linalg.h"
//spock_eq6_20170722
 /////////////////////////////////////////////////////////////////////////////////////////
 //
 //  Name:           load_options
 //  Purpose:        Loads params based input file
 //  Assumptions:    Spice files have to be in directory_to_cspice/data/ (where directory_to_cspice is indicated in the Makefile at SPICE_DIR)
 //  References      Various
 //
 //  Change Log:
 //      |   Developer   |       Date    |   SCR     |   Notes
 //      | --------------|---------------|-----------|-------------------------------
 //      | J. Getchius   | 05/20/2015    |   ---     | Initial Implementation
 //      | C. Bussy-Virat| 07/22/2015    |   ---     | Add orbital elements in the input file and reorganize it
 //      | C. Bussy-Virat| 07/23/2015    |   ---     | Replace altitude by apogee of altitude in the input file
 //      | C. Bussy-Virat| 07/29/2015    |   ---     | Added number of surfaces and normal to the surface
 //      | C. Bussy-Virat| 07/30/2015    |   ---     | Read from different input files (basic.d, orbit.d, spacecraft.d, attitude.d, force.d)
 //      | C. Bussy-Virat| 07/31/2015    |   ---     | For each surface of the sc, read the area, Cd, and normal vector
 //      | C. Bussy-Virat| 08/04/2015    |   ---     | The end time is read from the input file (not the number of days anymore)
 //
 /////////////////////////////////////////////////////////////////////////////////////////

 int load_options( OPTIONS_T *OPTIONS,
		   char      filename[1000], int iProc, int nProcs,  int iDebugLevel, char filename_input_raw[256]) {

   //newstructure
   char text_spice[256];
   //newstructure

   int do_not_download_file_swpc_or_wget = 0; // set this variable to 1 to prevent from downloading F10.7 and Ap files from SWPC or Omniweb
   /* Declarations */
	     double et_most_recent_tle_epoch;
	     char most_recent_tle_epoch[256];

	     // Initialize TLE parameters
	     SpiceInt frstyr = 1961; // do not really care about this line (as long as the spacecrafts are not flying after 2060)
	     SpiceDouble geophs[8];
	     int pp;
	     int lineln=0;
	     size_t len = 0;
	     char *line_temp = NULL;
	     SpiceDouble elems[10];

	     SpiceDouble epoch_sat; // epoch of TLE of the sat
	     char sat_tle_file_temp[256];
	     FILE *sat_tle_file = NULL;
	     ORBITAL_ELEMENTS_T OE_temp;

	     /* Set up the geophysical quantities.  At last check these were the values used by Space Command and SGP4 */
	     geophs[ 0 ] =    1.082616e-3;   // J2
	     geophs[ 1 ] =   -2.53881e-6;    // J3
	     geophs[ 2 ] =   -1.65597e-6;    // J4
	     geophs[ 3 ] =    7.43669161e-2; // KE
	     geophs[ 4 ] =    120.0;         // QO
	     geophs[ 5 ] =    78.0;          // SO
	     geophs[ 6 ] =    6378.135;      // ER
	     geophs[ 7 ] =    1.0;           // AE


   double   epoch_sat_previous;
	   int bbb;
	   int ii;
    OPTIONS->write_collision_files = 0;
   char filename_without_extension[1000];
   double et_end, et_start;
	 int ccc;
     double et_text_output_dir_now;
   char text_output_dir_now[256];
   double et_initial=0;
	     char filename_f107_ap_pred[1000], wget_filename_f107_ap_pred[1000];
	   char **list_filename_obs_swpc_f107, **list_filename_obs_swpc_ap;
	   int nb_file_obs_swpc = 0, nb_file_obs_swpc_ap = 0;
	   int ooo;
	     int quarter_initial_swpc, quarter_final_swpc, quarter_current, quarter_initial_swpc_ap, quarter_initial_epoch;
   char date_file_obs_initial_swpc[256];
   char date_file_obs_final_swpc[256];

   char initial_filename_f107_obs[256];
   char final_filename_f107_obs[256];
   int nb_storm_temp = 0;
   int find_path_cspice = 0;
   char homedir_str[1000];

   struct passwd *pw = getpwuid(getuid());
   const char *homedir = pw->pw_dir;
   char path_to_spice[256];
   char leap_sec_file_path[256];
   char eop_file_path[256];
   char planet_ephem_file_path[256];
   char earth_binary_pck_path[256];

   char dir_output_run_name_temp[1000], path_output_run_name_temp[100];
   char sat_nb_str[6];
   double ang_velo[6];//  double ang_velo[3];
   char str_no_gs[256];
   double et_correct_end_date_omniweb_minus_40_5_days;
   //  char test_omniweb_or_external_file[256];
   FILE *ap_wget_omniweb;
   FILE *f107_wget_omniweb;
   char filename_f107_to_calculate_f107_average[256];
   FILE *file_f107_to_calculate_f107_average;
   double et_initial_epoch, et_final_epoch;
   char initial_epoch_wget_temp[12], final_epoch_wget_temp[12];
   //newstructure
   //  char initial_epoch_wget[8], final_epoch_wget[8];
   char initial_epoch_wget[10], final_epoch_wget[10];
   //newstructure
   char initial_epoch_wget_minus_57hours_temp[12];
   //newstructure
   //  char initial_epoch_wget_minus_57hours[8];
   char initial_epoch_wget_minus_57hours[10];
   //newstructure
   char initial_epoch_wget_minus40_5days_temp[12], final_epoch_wget_plus40_5days_temp[12];
   //newstructure
   //  char initial_epoch_wget_minus40_5days[8], final_epoch_wget_plus40_5days[8];
   char initial_epoch_wget_minus40_5days[10], final_epoch_wget_plus40_5days[10];
   //newstructure
   char *next_wget_initial;
   int find_wget_initial;
   char *next_wget_final;
   int find_wget_final;
   char str_wget[1000];
   char test_wget_omniweb_error[10];
   //newstructure
   //  char correct_end_date_omniweb_temp[8];
   char correct_end_date_omniweb_temp[10];
   //newstructure
   char correct_end_date_omniweb[256];
   char correct_end_date_omniweb_minus_40_5_days[256];
   double et_correct_end_date_omniweb;
   int ierr;
   char name_sat_temp[256];
   int hhh;
   int have_found_older_time_than_me, count_have_found_older_time_than_me;
   int nnn;
   int i_gitm_file = 0;
   //  char            density_file[256];
   char density_file_temp[256];
   int more_ensembles_output, inc_ensembles_output = 0, inc_ensembles_output_final = 0;
   char text_location[256];
   //  int nb_satellites_not_including_gps;
   FILE *fp;
   char *line = NULL;

   int sss;

   int found_eoh = 0;
   char text[256], text_temp[256];
   char earth_binary_pck_path_temp[256];
   char filename_f107[256];
   char filename_f107A[256];
   char filename_ap[256];
   char *next, *new_next;
   int find_file_name, find_comments, find_file_name_sum = 0;
   char text_output[256];
   char text_surface[256];
   char answer_coe[10];

 //newstructure
   int err;
   struct stat s;
       OPTIONS->which_sc_oldest_tle_epoch = 0;
 /*   // Make sure all input directories exist. If they don't, then create them. */
 /*   if ((iProc == 0) & ( iDebugLevel >= 1 ) ){ */
 /*     printf("-- (load_options) Making sure all input directories exist. If they don't, then creating them... (iProc %d)\n", iProc); */
 /*   } */

 /*   struct stat s; */
 /*   int err; */
 /*   // // INPUT */
 /*   strcpy(OPTIONS->dir_input, ""); */
 /*   strcat(OPTIONS->dir_input, main_directory_location); */
 /*   strcat(OPTIONS->dir_input,"input"); */
 /*   err = stat(OPTIONS->dir_input, &s); */
 /*   if (err !=-1){ */
 /*     if(S_ISDIR(s.st_mode)) { // the dir exists */
 /*     }  */
 /*     else { // exists but it is not a dir */
 /*       strcat(OPTIONS->dir_input,"_dir"); */
 /*       mode_t process_mask = umask(0); */
 /*       mkdir(OPTIONS->dir_input, S_IRWXU | S_IRWXG | S_IRWXO); */
 /*       umask(process_mask); */
 /*     } */
 /*   } */
 /*   else{// the dir does not exist */
 /*     mode_t process_mask = umask(0); */
 /*     mkdir(OPTIONS->dir_input, S_IRWXU | S_IRWXG | S_IRWXO); */
 /*     umask(process_mask); */
 /*   } */

 /*   // // INPUT_GEOMETRY */
 /*   strcpy(OPTIONS->dir_input_geometry, ""); */
 /*   strcat(OPTIONS->dir_input_geometry, OPTIONS->dir_input); */
 /*   strcat(OPTIONS->dir_input_geometry,"/geometry"); */
 /*   err = stat(OPTIONS->dir_input_geometry, &s); */
 /*   if (err !=-1){ */
 /*     if(S_ISDIR(s.st_mode)) { // the dir exists */
 /*     }  */
 /*     else { // exists but it is not a dir */
 /*       strcat(OPTIONS->dir_input_geometry,"_dir"); */
 /*       mode_t process_mask = umask(0); */
 /*       mkdir(OPTIONS->dir_input_geometry, S_IRWXU | S_IRWXG | S_IRWXO); */
 /*       umask(process_mask); */
 /*     } */
 /*   } */
 /*   else{// the dir does not exist */
 /*     mode_t process_mask = umask(0); */
 /*     mkdir(OPTIONS->dir_input_geometry, S_IRWXU | S_IRWXG | S_IRWXO); */
 /*     umask(process_mask); */
 /*   } */

 /*   // // INPUT_COLLISION */
 /*   strcpy(OPTIONS->dir_input_collision, ""); */
 /*   strcat(OPTIONS->dir_input_collision, OPTIONS->dir_input); */
 /*   strcat(OPTIONS->dir_input_collision,"/collision"); */
 /*   err = stat(OPTIONS->dir_input_collision, &s); */
 /*   if (err !=-1){ */
 /*     if(S_ISDIR(s.st_mode)) { // the dir exists */
 /*     }  */
 /*     else { // exists but it is not a dir */
 /*       strcat(OPTIONS->dir_input_collision,"_dir"); */
 /*       mode_t process_mask = umask(0); */
 /*       mkdir(OPTIONS->dir_input_collision, S_IRWXU | S_IRWXG | S_IRWXO); */
 /*       umask(process_mask); */
 /*     } */
 /*   } */
 /*   else{// the dir does not exist */
 /*     mode_t process_mask = umask(0); */
 /*     mkdir(OPTIONS->dir_input_collision, S_IRWXU | S_IRWXG | S_IRWXO); */
 /*     umask(process_mask); */
 /*   } */

 /*   // // INPUT_COVERAGE */
 /*   strcpy(OPTIONS->dir_input_coverage, ""); */
 /*   strcat(OPTIONS->dir_input_coverage, OPTIONS->dir_input); */
 /*   strcat(OPTIONS->dir_input_coverage,"/coverage"); */
 /*   err = stat(OPTIONS->dir_input_coverage, &s); */
 /*   if (err !=-1){ */
 /*     if(S_ISDIR(s.st_mode)) { // the dir exists */
 /*     }  */
 /*     else { // exists but it is not a dir */
 /*       strcat(OPTIONS->dir_input_coverage,"_dir"); */
/*       mode_t process_mask = umask(0); */
/*       mkdir(OPTIONS->dir_input_coverage, S_IRWXU | S_IRWXG | S_IRWXO); */
/*       umask(process_mask); */
/*     } */
/*   } */
/*   else{// the dir does not exist */
/*     mode_t process_mask = umask(0); */
/*     mkdir(OPTIONS->dir_input_coverage, S_IRWXU | S_IRWXG | S_IRWXO); */
/*     umask(process_mask); */
/*   } */

/*   // // INPUT_COVERAGE_STORM */
/*   strcpy(OPTIONS->dir_input_coverage_storm, ""); */
/*   strcat(OPTIONS->dir_input_coverage_storm, OPTIONS->dir_input_coverage); */
/*   strcat(OPTIONS->dir_input_coverage_storm,"/storm"); */
/*   err = stat(OPTIONS->dir_input_coverage_storm, &s); */
/*   if (err !=-1){ */
/*     if(S_ISDIR(s.st_mode)) { // the dir exists */
/*     }  */
/*     else { // exists but it is not a dir */
/*       strcat(OPTIONS->dir_input_coverage_storm,"_dir"); */
/*       mode_t process_mask = umask(0); */
/*       mkdir(OPTIONS->dir_input_coverage_storm, S_IRWXU | S_IRWXG | S_IRWXO); */
/*       umask(process_mask); */
/*     } */
/*   } */
/*   else{// the dir does not exist */
/*     mode_t process_mask = umask(0); */
/*     mkdir(OPTIONS->dir_input_coverage_storm, S_IRWXU | S_IRWXG | S_IRWXO); */
/*     umask(process_mask); */
/*   } */



/*   // // INPUT_COVERAGE_GROUND_STATION */
/*   strcpy(OPTIONS->dir_input_coverage_ground_station, ""); */
/*   strcat(OPTIONS->dir_input_coverage_ground_station, OPTIONS->dir_input_coverage); */
/*   strcat(OPTIONS->dir_input_coverage_ground_station,"/ground_station"); */
/*   err = stat(OPTIONS->dir_input_coverage_ground_station, &s); */
/*   if (err !=-1){ */
/*     if(S_ISDIR(s.st_mode)) { // the dir exists */
/*     }  */
/*     else { // exists but it is not a dir */
/*       strcat(OPTIONS->dir_input_coverage_ground_station,"_dir"); */
/*       mode_t process_mask = umask(0); */
/*       mkdir(OPTIONS->dir_input_coverage_ground_station, S_IRWXU | S_IRWXG | S_IRWXO); */
/*       umask(process_mask); */
/*     } */
/*   } */
/*   else{// the dir does not exist */
/*     mode_t process_mask = umask(0); */
/*     mkdir(OPTIONS->dir_input_coverage_ground_station, S_IRWXU | S_IRWXG | S_IRWXO); */
/*     umask(process_mask); */
/*   } */




/*   // // INPUT_DENSITY */
/*   strcpy(OPTIONS->dir_input_density, ""); */
/*   strcat(OPTIONS->dir_input_density, OPTIONS->dir_input); */
/*   strcat(OPTIONS->dir_input_density,"/density"); */
/*   err = stat(OPTIONS->dir_input_density, &s); */
/*   if (err !=-1){ */
/*     if(S_ISDIR(s.st_mode)) { // the dir exists */
/*     }  */
/*     else { // exists but it is not a dir */
/*       strcat(OPTIONS->dir_input_density,"_dir"); */
/*       mode_t process_mask = umask(0); */
/*       mkdir(OPTIONS->dir_input_density, S_IRWXU | S_IRWXG | S_IRWXO); */
/*       umask(process_mask); */
/*     } */
/*   } */
/*   else{// the dir does not exist */
/*     mode_t process_mask = umask(0); */
/*     mkdir(OPTIONS->dir_input_density, S_IRWXU | S_IRWXG | S_IRWXO); */
/*     umask(process_mask); */
/*   } */


/*   // // INPUT_DENSITY_MSIS */
/*   strcpy(OPTIONS->dir_input_density_msis, ""); */
/*   strcat(OPTIONS->dir_input_density_msis, OPTIONS->dir_input_density); */
/*   strcat(OPTIONS->dir_input_density_msis,"/density_NRLMSIS00e"); */
/*   err = stat(OPTIONS->dir_input_density_msis, &s); */
/*   if (err !=-1){ */
/*     if(S_ISDIR(s.st_mode)) { // the dir exists */
/*     }  */
/*     else { // exists but it is not a dir */
/*       strcat(OPTIONS->dir_input_density_msis,"_dir"); */
/*       mode_t process_mask = umask(0); */
/*       mkdir(OPTIONS->dir_input_density_msis, S_IRWXU | S_IRWXG | S_IRWXO); */
/*       umask(process_mask); */
/*     } */
/*   } */
/*   else{// the dir does not exist */
/*     mode_t process_mask = umask(0); */
/*     mkdir(OPTIONS->dir_input_density_msis, S_IRWXU | S_IRWXG | S_IRWXO); */
/*     umask(process_mask); */
/*   } */


/*   // // INPUT_DENSITY_GITM */
/*   strcpy(OPTIONS->dir_input_density_gitm, ""); */
/*   strcat(OPTIONS->dir_input_density_gitm, OPTIONS->dir_input_density); */
/*   strcat(OPTIONS->dir_input_density_gitm,"/density_gitm"); */
/*   err = stat(OPTIONS->dir_input_density_gitm, &s); */
/*   if (err !=-1){ */
/*     if(S_ISDIR(s.st_mode)) { // the dir exists */
/*     }  */
/*     else { // exists but it is not a dir */
/*       strcat(OPTIONS->dir_input_density_gitm,"_dir"); */
/*       mode_t process_mask = umask(0); */
/*       mkdir(OPTIONS->dir_input_density_gitm, S_IRWXU | S_IRWXG | S_IRWXO); */
/*       umask(process_mask); */
/*     } */
/*   } */
/*   else{// the dir does not exist */
/*     mode_t process_mask = umask(0); */
/*     mkdir(OPTIONS->dir_input_density_gitm, S_IRWXU | S_IRWXG | S_IRWXO); */
/*     umask(process_mask); */
/*   } */

/*   // // INPUT_DENSITY_DENSITY_FILES */
/*   strcpy(OPTIONS->dir_input_density_density_files, ""); */
/*   strcat(OPTIONS->dir_input_density_density_files, OPTIONS->dir_input_density); */
/*   strcat(OPTIONS->dir_input_density_density_files,"/density_density_files"); */
/*   err = stat(OPTIONS->dir_input_density_density_files, &s); */
/*   if (err !=-1){ */
/*     if(S_ISDIR(s.st_mode)) { // the dir exists */
/*     }  */
/*     else { // exists but it is not a dir */
/*       strcat(OPTIONS->dir_input_density_density_files,"_dir"); */
/*       mode_t process_mask = umask(0); */
/*       mkdir(OPTIONS->dir_input_density_density_files, S_IRWXU | S_IRWXG | S_IRWXO); */
/*       umask(process_mask); */
/*     } */
/*   } */
/*   else{// the dir does not exist */
/*     mode_t process_mask = umask(0); */
/*     mkdir(OPTIONS->dir_input_density_density_files, S_IRWXU | S_IRWXG | S_IRWXO); */
/*     umask(process_mask); */
/*   } */


/*   // // INPUT_ATTITUDE */
/*   strcpy(OPTIONS->dir_input_attitude, ""); */
/*   strcat(OPTIONS->dir_input_attitude, OPTIONS->dir_input); */
/*   strcat(OPTIONS->dir_input_attitude,"/attitude"); */
/*   err = stat(OPTIONS->dir_input_attitude, &s); */
/*   if (err !=-1){ */
/*     if(S_ISDIR(s.st_mode)) { // the dir exists */
/*     }  */
/*     else { // exists but it is not a dir */
/*       strcat(OPTIONS->dir_input_attitude,"_dir"); */
/*       mode_t process_mask = umask(0); */
/*       mkdir(OPTIONS->dir_input_attitude, S_IRWXU | S_IRWXG | S_IRWXO); */
/*       umask(process_mask); */
/*     } */
/*   } */
/*   else{// the dir does not exist */
/*     mode_t process_mask = umask(0); */
/*     mkdir(OPTIONS->dir_input_attitude, S_IRWXU | S_IRWXG | S_IRWXO); */
/*     umask(process_mask); */
/*   } */


/*   // // INPUT_TLE */
/*   strcpy(OPTIONS->dir_input_tle, ""); */
/*   strcat(OPTIONS->dir_input_tle, OPTIONS->dir_input); */
/*   strcat(OPTIONS->dir_input_tle,"/tle"); */
/*   err = stat(OPTIONS->dir_input_tle, &s); */
/*   if (err !=-1){ */
/*     if(S_ISDIR(s.st_mode)) { // the dir exists */
/*     }  */
/*     else { // exists but it is not a dir */
/*       strcat(OPTIONS->dir_input_tle,"_dir"); */
/*       mode_t process_mask = umask(0); */
/*       mkdir(OPTIONS->dir_input_tle, S_IRWXU | S_IRWXG | S_IRWXO); */
/*       umask(process_mask); */
/*     } */
/*   } */
/*   else{// the dir does not exist */
/*     mode_t process_mask = umask(0); */
/*     mkdir(OPTIONS->dir_input_tle, S_IRWXU | S_IRWXG | S_IRWXO); */
/*     umask(process_mask); */
/*   } */


/*   // // INPUT_GPS_TLE */
/*   strcpy(OPTIONS->dir_input_tle_gps_tle, ""); */
/*   strcat(OPTIONS->dir_input_tle_gps_tle, OPTIONS->dir_input_tle); */
/*   strcat(OPTIONS->dir_input_tle_gps_tle,"/constellation_gps_tle"); */
/*   err = stat(OPTIONS->dir_input_tle_gps_tle, &s); */
/*   if (err !=-1){ */
/*     if(S_ISDIR(s.st_mode)) { // the dir exists */
/*     }  */
/*     else { // exists but it is not a dir */
/*       strcat(OPTIONS->dir_input_tle_gps_tle,"_dir"); */
/*       mode_t process_mask = umask(0); */
/*       mkdir(OPTIONS->dir_input_tle_gps_tle, S_IRWXU | S_IRWXG | S_IRWXO); */
/*       umask(process_mask); */
/*     } */
/*   } */
/*   else{// the dir does not exist */
/*     mode_t process_mask = umask(0); */
/*     mkdir(OPTIONS->dir_input_tle_gps_tle , S_IRWXU | S_IRWXG | S_IRWXO); */
/*     umask(process_mask); */
/*   } */
/*   if ((iProc == 0) & ( iDebugLevel >= 1 ) ){ */
/*     printf("-- (load_options) The input directories have been correctly created.\n"); */
/*   } */
/*   // end of make sure all input directories exist. If they don't, then create them.  */


//  fp = fopen(filename, "r"); //newstructure
       char filename_input_no_path[1000]; // filename_input_raw is the name of the input file as indicated by the argument of the SpOCK call
       // filename_input_no_path is filename_input_raw without the path to the file 
       strcpy(filename_input_no_path, "");    
       if (strrchr(&filename_input_raw[0], '/')== NULL){ // not '/' character
	 // in the name indicated at 1st line of section #OUTPUT
	 strcpy(filename_input_no_path, filename_input_raw);
       }
       else{ // '/' character in the name indicated at 1st line of section #OUTPUT
	 strncat(filename_input_no_path,strrchr(&filename_input_raw[0], '/')+1, (int)(strrchr(&filename_input_raw[0], '/') + strlen(filename_input_raw) ) -1);
       }

       fp = fopen(filename_input_raw, "r"); //newstructure  

  if (fp == NULL){
    printf("***! The main input file was not found. The program will stop. !***\n");
    ierr =  MPI_Finalize();
    exit(0);
  }
  if ((iProc == 0) & ( iDebugLevel >= 1 ) ){
    if (fp != NULL){
      printf("-- (load_options) The main input file was correctly open.\n");
    }
  }


//newstructure
// Read main input file to figure out the path to spice (where leap second file, and other files are). You need to read the leap second file before reading section TIME. But then you need to read the section TIME before loading the other SPICE kernels
  /* SPICE */
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section #SPICE.\n");
  }

  found_eoh = 0;
  rewind(fp);
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#SPICE", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  if (feof(fp)){
    printf("***! No section #SPICE found in %s. The program will stop. !***\n", filename);
    ierr =  MPI_Finalize();
    exit(0);
  }
  // Path of SPICE
  getline(&line, &len, fp);
  sscanf(line, "%s", OPTIONS->path_to_spice);// ASSUMPTION: the path of spice specified in the main inoput file in the sectioN #SPICE must be the same as the path of SPICE installation in the Makefile, with /data at the end. So if in the Makefile the SPICE directory is hello/hi/ then the path of spice in the main input file in section #SPICE must be hello/hi/data/
  strcat(OPTIONS->path_to_spice, "/"); // don't put the '/' at the end of the path in section #SPICE

  // Does the user wants to use a particular earth_binary_pck file (yes if line below is not empty or does not start with a #)
  getline(&line, &len, fp);
  strcpy(text_spice, "");
  sscanf(line, "%s", text_spice);
  if ( (strlen(text_spice) == 0 ) || ( text_spice[0] == '#' ) ){ // if the user left the line (below the path to SPICE) empty of if it's the start of another section ('#' character) then it means the user wants to use the default earth_binary_pck file
    strcpy(text_spice, "default");
  }

  // // Now load the leap second file. ASSUMPTION: it has to be in directory_to_cspice/data/ (where directory_to_cspice is indicated in the Makefile at SPICE_DIR)
  strcpy(leap_sec_file_path,OPTIONS->path_to_spice);
  strcat(leap_sec_file_path, "naif0012.tls");
  strcpy(OPTIONS->leap_sec_file, leap_sec_file_path );
  furnsh_c( OPTIONS->leap_sec_file);
  // end of sectio SPICE until we get back to it after reading section #TIME
//newstructure

//newstructure                                                                                
/*   // Before reading the time, the leap second file from SPICE needs to be downloaded */
/*   // // open Makefile to read the path to spice */
/*   FILE *makefile_prop = NULL; */
/*   makefile_prop = fopen("../Makefile", "r"); */
/*   if (makefile_prop == NULL){ */
/*     printf("***! The Makefile could not be open. You need to run the propagator from a subfolder of the Propagator (ex: Propagator/subfolder_where_you_make_your_runs). The program will stop.!***\n"); MPI_Finalize(); exit(0); */
/*   } */
/*   int found_spice_dir = 0; */
/*   while ( found_spice_dir == 0 && !feof(makefile_prop)) { */
/*     getline(&line, &len, makefile_prop); */
/*     sscanf(line, "%s", text); */
/*     if (  strcmp( "SPICE_DIR", text  ) == 0 )  { */
/*       found_spice_dir = 1; */
/*     } */
/*   } */
/*   if (feof(makefile_prop)){ */
/*     printf("***! Could not find the path to the spice directory in Makefile (SPICE_DIR). The program will stop.!***\n"); MPI_Finalize(); exit(0); */
/*   } */
/*   RemoveSpaces(line); */
/*   strtok(line, "\n");strtok(line, "\r"); */
/*   next = &line[0]; */
/*   find_path_cspice = (int)(strchr(next, '=') - next); */
/*   strcpy(path_to_spice, ""); */
/*   strncat(path_to_spice, strchr(next, '=') + 1, (int)(strlen(line))-1); */
/*   strcpy(homedir_str, homedir); */
/*   if (strstr(path_to_spice, "${HOME}") != NULL){ */
/*     strcpy(path_to_spice, str_replace(path_to_spice, "${HOME}", homedir_str)); */
/*   } */
/*   strcat(path_to_spice, "/data/"); */
/*   // // Now load the leap second file. ASSUMPTION: it has to be in directory_to_cspice/data/ (where directory_to_cspice is indicated in the Makefile at SPICE_DIR) */
/*   strcpy(leap_sec_file_path,path_to_spice); */
/*   strcat(leap_sec_file_path, "naif0011.tls"); */
/*   strcpy(OPTIONS->leap_sec_file, leap_sec_file_path ); */
/*   furnsh_c( OPTIONS->leap_sec_file); */
//newstructure                                                                                 

  /* TIME */
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section #TIME.\n");
  }

  found_eoh = 0;
  rewind(fp);
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#TIME", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  if (feof(fp)){
    printf("***! No section #TIME found in %s. The program will stop. !***\n", filename);
    ierr =  MPI_Finalize();
    exit(0);

  }

  // Date of initial epoch
  getline(&line, &len, fp);
  sscanf(line, "%s", text);
  if ( strcmp(text, "now") == 0){
    if (iProc == 0){
    time_t t = time(NULL);
    struct tm tm = *gmtime(&t);    
    int year_now = tm.tm_year + 1900;
    int month_now = tm.tm_mon + 1;
    int day_now = tm.tm_mday;
    int hour_now = tm.tm_hour;
    int minute_now = tm.tm_min;
    int second_now = tm.tm_sec;

    char year_now_str[15];
    sprintf(year_now_str, "%d", year_now);
    char month_now_str[15];
    sprintf(month_now_str, "%d", month_now);
    char day_now_str[15];
    sprintf(day_now_str, "%d", day_now);
    char hour_now_str[15];
    sprintf(hour_now_str, "%d", hour_now);
    char minute_now_str[15];
    sprintf(minute_now_str, "%d", minute_now);
    char second_now_str[15];
    sprintf(second_now_str, "%d", second_now);
    
    char date_now[256];
    strcpy(date_now,year_now_str);
    strcat(date_now, "-");
    strcat(date_now, month_now_str);
    strcat(date_now, "-");
    strcat(date_now, day_now_str);
    strcat(date_now, " ");
    strcat(date_now, hour_now_str);
    strcat(date_now, ":");
    strcat(date_now, minute_now_str);
    strcat(date_now, ":");
    strcat(date_now, second_now_str);
  
    strcpy(OPTIONS->initial_epoch, date_now );
	str2et_c(OPTIONS->initial_epoch, &et_initial); // don't worry about it

	if (nProcs > 1){
    	for (ccc = 1; ccc < nProcs; ccc++){
    MPI_Send(&et_initial, 1, MPI_DOUBLE, ccc, 0, MPI_COMM_WORLD);
	}
	}
    }
    else{
      if (nProcs > 1){
      MPI_Recv(&et_initial, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
      //      printf("Iproc %d has value %f\n", iProc, et_initial);
      }
    }


	et2utc_c(et_initial, "ISOC", 6, 255, OPTIONS->initial_epoch); // don't worry about it

    // Date of final epoch (= now + number of hours chosen in the input file input.d)
    rewind(fp);
    found_eoh=0;
    while ( found_eoh == 0 && !feof(fp)) {
      getline(&line, &len, fp);
      sscanf(line, "%s", text);
      if (  strcmp( "#TIME", text  ) == 0 )  {
	found_eoh = 1;
      }
    }
    if (feof(fp)){
      printf("***! No section #TIME found in %s. The program will stop. !***\n", filename);
      ierr =  MPI_Finalize();
      exit(0);

    }

    getline(&line, &len, fp);
    double nb_hours_simu;
    char text_useless[5];
    sscanf(line, "%4[^\n] %lf", text_useless,&nb_hours_simu);
    double et_now;
    str2et_c(OPTIONS->initial_epoch, &et_now); 
    et_end = et_now + nb_hours_simu * 3600.;
    char times_end[300];
    et2utc_c(et_end, "ISOC" ,6 ,255 , times_end);    
    strcpy(OPTIONS->final_epoch, times_end );
    getline(&line, &len, fp);
  }
  else{
    strcpy(OPTIONS->initial_epoch, text );
    // Date of final epoch
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    strcpy(OPTIONS->final_epoch, text );

  }
  // Time step of simulation
  getline(&line, &len, fp);
  RemoveSpaces(line);
  sscanf( line, "%lf", &OPTIONS->dt);

  // Re-evaluate OPTIONS->final epoch and et_end so that final epoch is a multiple of dt + inital epoch (if originally it is not, then take the closest multiple just before the orginal value of final epoch (ex: if initial epoch is at 12:00:00, dt = 20s, and final epoch as written in main input file is 12:00:50 then re-evaluate final epoch to be at 12:00:40))
  
  str2et_c(OPTIONS->final_epoch, &et_end);
  str2et_c(OPTIONS->initial_epoch, &et_start);
  OPTIONS->dt_pos_neg = OPTIONS->dt;
  if (et_start > et_end){
    OPTIONS->dt_pos_neg = -OPTIONS->dt;
  }
  //  if ( fabs(OPTIONS->dt -  fmod( et_end - et_start, OPTIONS->dt ) ) > 0.01 ){ // 0.01 for numerical reasons
  if ( fabs(fmod( et_end - et_start, OPTIONS->dt ) ) > 1e-6 ){ // 0.01 for numerical reasons
    /* printf("X = %f\n", fabs(OPTIONS->dt - fmod( et_end - et_start, OPTIONS->dt ) )); */
    /* printf("old: <%s>\n",OPTIONS->final_epoch); */
    et_end = (int) ( (et_end - et_start) / OPTIONS->dt ) * OPTIONS->dt + et_start;
    et2utc_c(et_end, "ISOC", 6, 300, OPTIONS->final_epoch);

    if (iProc == 0){
      //      printf("***! The final epoch has been reset to %s (so it is equal to initial epoch + a multiple of dt). !***\n", OPTIONS->final_epoch);
      }
    //    printf("new: <%s>\n",OPTIONS->final_epoch);
}
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) End of section #TIME.\n");
  }
  //    etprint(et_end, "new again");


//newstructure
// A bit before, we read main input file to figure out the path to spice (where leap second file, and other files are). Now we figure out the other SPICE kernels to load
  /* SPICE */
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section #SPICE.\n");
  }

  found_eoh = 0;
  rewind(fp);
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#SPICE", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  if (feof(fp)){
    printf("***! No section #SPICE found in %s. The program will stop. !***\n", filename);
    ierr =  MPI_Finalize();
    exit(0);
  }

  if ( strcmp(text_spice, "default") == 0 ){ // if 'default' chosen in section #SPICE. So that means that either the line below the path was empty or that the user put 'default' in this line
    strcpy(earth_binary_pck_path,OPTIONS->path_to_spice);
    strcpy(earth_binary_pck_path_temp, "earth_000101_210404_210111.bpc");// !!!!! used to be earth_000101_191026_190804.bpc
    strcat(earth_binary_pck_path, "earth_000101_210404_210111.bpc");// !!!!! used to be "earth_000101_210404_210111.bpc");
    // bpc to use for runs before 2000: earth_720101_070527_NameModifiedByCbvToAddUnderscore.bpc
  }
  else { // if the user wants to use a particular earth_binary_pck file, then this file needs to be put in path_to_spice  (with the other spice files). IMPORTANT: the name of the file must have the format string1_string2_YYMMDD_string3.bpc where YYMMDD is the date of the last available information on the frame (so we don't know anything for the rotating frame any time after this date). Example: earth_000101_210404_210111.bpc (after 160830 we can't use ITRF93)
    strcpy(earth_binary_pck_path,OPTIONS->path_to_spice);
    strcpy(earth_binary_pck_path_temp, text_spice);
    strcat(earth_binary_pck_path, text_spice);
  }

  strcpy(eop_file_path,OPTIONS->path_to_spice);
  strcat(eop_file_path, "pck00010.tpc");
  strcpy(planet_ephem_file_path,OPTIONS->path_to_spice);
  strcat(planet_ephem_file_path, "de432s.bsp");
  strcpy(OPTIONS->eop_file, eop_file_path);
  strcpy(OPTIONS->planet_ephem_file, planet_ephem_file_path );
  strcpy(OPTIONS->earth_binary_pck, earth_binary_pck_path );
  // Figure out if the ITRF93 Earth body-fixed rotating frame for conversion ECI to ECEF is available for the epoch of the propagation. If it is not, then use the less accurate IAU_EARTH rotating frame
  // // from the name of the earth_binary_pck file, figure out the date after which we can't use ITRF93 (160830 in the example a few lines above)
  next = &earth_binary_pck_path_temp[0];
  char earth_binary_pck_path_temp_new[256];
  strcpy(earth_binary_pck_path_temp_new, "");
  strncat(earth_binary_pck_path_temp_new, strchr(next, '_')+1,  (int)(strrchr(next, '_') - next) - (int)(strchr(next, '_') - next) -1);
  next = &earth_binary_pck_path_temp_new[0];
  char earth_binary_pck_path_temp_new_really[256];
  strcpy(earth_binary_pck_path_temp_new_really, "");
  strncat(earth_binary_pck_path_temp_new_really, strchr(next, '_')+1, (int)(strlen(earth_binary_pck_path_temp_new)-1));
  // // now we know this date, convert it into the correct format (so we can compare it to the end epoch time of the constellation). ASSUMPTION: have to be after 2000
  char earth_binary_pck_date[100];
  strcpy(earth_binary_pck_date, "20");
  strncat(earth_binary_pck_date, &earth_binary_pck_path_temp_new_really[0], 2);
  strcat(earth_binary_pck_date, "-");
  strncat(earth_binary_pck_date, &earth_binary_pck_path_temp_new_really[0]+2, 2);
  strcat(earth_binary_pck_date, "-");
  strncat(earth_binary_pck_date, &earth_binary_pck_path_temp_new_really[0]+4, 2);
  strcat(earth_binary_pck_date, "T00:00");
  double et_earth_binary_pck_date;
  str2et_c(earth_binary_pck_date, &et_earth_binary_pck_date);
  // // Compare this date to the end epoch of the constellation
  str2et_c(OPTIONS->final_epoch, &et_final_epoch);
  int accurate_earth_fixed_frame_available = 0;


  if (et_final_epoch + 24 * 3600. > et_earth_binary_pck_date){ // if the propagations ends later than the date of the pck file then ITRF93 can't be used (with one day margin)
    accurate_earth_fixed_frame_available = 0;
  }
  else{
    accurate_earth_fixed_frame_available = 1;
  }
  strcpy(OPTIONS->earth_fixed_frame, "");
  if (accurate_earth_fixed_frame_available == 1){
    strcpy(OPTIONS->earth_fixed_frame, "ITRF93");
  }
  else{
    strcpy(OPTIONS->earth_fixed_frame, "IAU_EARTH");
  }
  if ((iProc == 0) & ( iDebugLevel >= 3 ) ){
    printf("---- (load_options) Earth body-fixed rotating frame for conversion ECI to ECEF: %s\n", OPTIONS->earth_fixed_frame);
  }

  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) End of section #SPICE\n");
  }

//newstructure

//newstructure                                                                                
/*   /\* SPICEFILES *\/ //  ASSUMPTION: spice files have to be in directory_to_cspice/data/ (where directory_to_cspice is indicated in the Makefile at SPICE_DIR) */
/*   if ((iProc == 0) & ( iDebugLevel >= 2 )){ */
/*     printf("--- (load_options) Beginning of section #SPICEFILES.\n"); // the leap second file was loaded at the very beginning of load_options */
/*   } */
/*   rewind(fp); */
/*   found_eoh = 0; */
/*   while ( found_eoh == 0 && !feof(fp)) { */
/*     getline(&line, &len, fp); */
/*     sscanf(line, "%s", text); */
/*     if (  strcmp( "#SPICEFILES", text  ) == 0 )  { */
/*       found_eoh = 1; */
/*     } */
/*   } */
/*   if (feof(fp)){ // if there is no section #SPICEFILES in the main input file, then the default files are loaded */
/*     strcpy(text, "default"); */
/*   } */
/*   else{ // if there is a section #SPICEFILES in the main input file then there 2 choices: 'default' or the name of the earth_binary_pck_path file to load */
/*     getline(&line,&len,fp); */
/*     sscanf(line,"%s", text); */
/*   } */

/*   if ( strcmp(text, "default") == 0 ){ // if 'default' chosen in section #SPICEFILES */
/*     strcpy(earth_binary_pck_path,path_to_spice); */
/*     strcpy(earth_binary_pck_path_temp, "earth_000101_210404_210111.bpc"); */
/*     strcat(earth_binary_pck_path, "earth_000101_210404_210111.bpc"); */
/*   }  */
/*   else { // if the user wants to use a particular earth_binary_pck file, then this file needs to be put in the cspice/data/ directory (with the other spice files). IMPORTANT: the name of the file must have the format string1_string2_YYMMDD_string3.bpc where YYMMDD is the date of the last available information on the frame (so we don't know anything for the rotating frame any time after this date). Example: earth_000101_191026_190804  .bpc (after 160830 we can't use ITRF93) */
/*     strcpy(earth_binary_pck_path,path_to_spice); */
/*     strcpy(earth_binary_pck_path_temp, text); */
/*     strcat(earth_binary_pck_path, text); */
/*   } */
/*   strcpy(eop_file_path,path_to_spice); */
/*   strcat(eop_file_path, "pck00010.tpc"); */
/*   strcpy(planet_ephem_file_path,path_to_spice); */
/*   strcat(planet_ephem_file_path, "de430.bsp"); */
/*   strcpy(OPTIONS->eop_file, eop_file_path); */
/*   strcpy(OPTIONS->planet_ephem_file, planet_ephem_file_path ); */
/*   strcpy(OPTIONS->earth_binary_pck, earth_binary_pck_path ); */
/*   // Figure out if the ITRF93 Earth body-fixed rotating frame for conversion ECI to ECEF is available for the epoch of the propagation. If it is not, then use the less accurate IAU_EARTH rotating frame  */
/*   // // from the name of the earth_binary_pck file, figure out the date after which we can't use ITRF93 (160830 in the example a few lines above) */
/*   next = &earth_binary_pck_path_temp[0]; */
/*   char earth_binary_pck_path_temp_new[256]; */
/*   strcpy(earth_binary_pck_path_temp_new, ""); */
/*   strncat(earth_binary_pck_path_temp_new, strchr(next, '_')+1,  (int)(strrchr(next, '_') - next) - (int)(strchr(next, '_') - next) -1); */
/*   next = &earth_binary_pck_path_temp_new[0]; */
/*   char earth_binary_pck_path_temp_new_really[256]; */
/*   strcpy(earth_binary_pck_path_temp_new_really, ""); */
/*   strncat(earth_binary_pck_path_temp_new_really, strchr(next, '_')+1, (int)(strlen(earth_binary_pck_path_temp_new)-1)); */
/*   // // now we know this date, convert it into the correct format (so we can compare it to the end epoch time of the constellation). ASSUMPTION: have to be after 2000 */
/*   char earth_binary_pck_date[100]; */
/*   strcpy(earth_binary_pck_date, "20"); */
/*   strncat(earth_binary_pck_date, &earth_binary_pck_path_temp_new_really[0], 2); */
/*   strcat(earth_binary_pck_date, "-"); */
/*   strncat(earth_binary_pck_date, &earth_binary_pck_path_temp_new_really[0]+2, 2); */
/*   strcat(earth_binary_pck_date, "-"); */
/*   strncat(earth_binary_pck_date, &earth_binary_pck_path_temp_new_really[0]+4, 2); */
/*   strcat(earth_binary_pck_date, "T00:00"); */
/*   double et_earth_binary_pck_date; */
/*   str2et_c(earth_binary_pck_date, &et_earth_binary_pck_date);  */
/*   // // Compare this date to the end epoch of the constellation */
/*   str2et_c(OPTIONS->final_epoch, &et_final_epoch); */
/*   int accurate_earth_fixed_frame_available = 0; */
/*   if (et_final_epoch + 24 * 3600. > et_earth_binary_pck_date){ // if the propagations ends later than the date of the pck file then ITRF93 can't be used (with one day margin) */
/*     accurate_earth_fixed_frame_available = 0; */
/*   } */
/*   else{ */
/*     accurate_earth_fixed_frame_available = 1; */
/*   } */
/*   strcpy(OPTIONS->earth_fixed_frame, ""); */
/*   if (accurate_earth_fixed_frame_available == 1){ */
/*     strcpy(OPTIONS->earth_fixed_frame, "ITRF93"); */
/*   } */
/*   else{ */
/*     strcpy(OPTIONS->earth_fixed_frame, "IAU_EARTH"); */
/*   } */
/*   if ((iProc == 0) & ( iDebugLevel >= 3 ) ){ */
/*     printf("---- (load_options) Earth body-fixed rotating frame for conversion ECI to ECEF: %s\n", OPTIONS->earth_fixed_frame); */
/*   } */

/*   if ((iProc == 0) & ( iDebugLevel >= 2 )){ */
/*     printf("--- (load_options) End of section #SPICEFILES.\n"); */
/*   } */
/* //newstructure */


  //  printf("<%s>\n", OPTIONS->earth_fixed_frame);







  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section #DENSITY_MOD.\n");
  }
  rewind(fp);
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#DENSITY_MOD", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  if (feof(fp)){
    OPTIONS->density_mod = 1;
   OPTIONS->density_mod_amp = 0;
   OPTIONS->density_mod_phase = 0;

  }

  else{
    int nbb;
  getline(&line,&len,fp);
   OPTIONS->density_mod = 1;
   OPTIONS->density_mod_amp = 0;
   OPTIONS->density_mod_phase = 0;
  sscanf(line,"%lf %lf %lf",&OPTIONS->density_mod, &OPTIONS->density_mod_amp, &OPTIONS->density_mod_phase); // the desnity given by msis is multiplied by density_mod + density_amp * sin(2*pi*t/T + density_phase*T) where T is the orbital period  
  //  printf("%f %f %f\n", OPTIONS->density_mod, OPTIONS->density_mod_amp, OPTIONS->density_mod_phase); exitf();
  }

  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) End of section #DENSITY_MOD.\n");
  }


    if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section #THRUST.\n");
  }
  rewind(fp);
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#THRUST", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  if (feof(fp)){
    OPTIONS->thrust = 0;
  }

  else{
        OPTIONS->thrust = 1;

	strcpy(OPTIONS->thrust_filename, "");
  getline(&line,&len,fp);
  sscanf(line,"%s", OPTIONS->thrust_filename);
  read_thrust(OPTIONS);

  //  etprint(OPTIONS->et_thrust_start, "start");
  //etprint(OPTIONS->et_thrust_stop, "stop");  
    //       printf("<%f %f %f>\n", OPTIONS->thrust_accel_lvlh[0], OPTIONS->thrust_accel_lvlh[1], OPTIONS->thrust_accel_lvlh[2]);
      
  }

  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) End of section #THRUST.\n");
  }



  /* ORBIT just to determine if collision_vcm optoin is on, nin which case don't look at geometry file. ORBIT sseciton will be read again in details later*/
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section #ORBIT.\n");
  }
  rewind(fp);
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#ORBIT", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  if (feof(fp)){
    printf("***! No section #ORBIT found in %s. The program will stop. !***\n", filename);
    ierr =  MPI_Finalize();
    exit(0);

  }

  getline(&line,&len,fp);
  char type_orbit_initialisation_temp[300];
  sscanf(line,"%s",type_orbit_initialisation_temp);
  

  /* SPACECRAFTS */
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section #SPACECRAFT.\n");
  }
  rewind(fp);
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#SPACECRAFT", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  if (feof(fp)){
    printf("***! No section #SPACECRAFT found in %s. The program will stop. !***\n", filename);
    ierr =  MPI_Finalize();
    exit(0);

  }

  // Number of spacecrafts
  getline(&line,&len,fp);
  sscanf(line,"%d",&OPTIONS->n_satellites);
  // Propagate the GPS satellites or not
  getline(&line,&len,fp);
  sscanf(line,"%s",text_location);
  //newstructure
/*   strcpy(OPTIONS->tle_constellation_gps_filename, OPTIONS->dir_input_tle_gps_tle); */
/*   strcat(OPTIONS->tle_constellation_gps_filename, "/"); */
  strcpy(OPTIONS->tle_constellation_gps_filename,  "");
  //newstructure
  strcat(OPTIONS->tle_constellation_gps_filename, text_location);

  if ( strcmp(text_location, "0") != 0 ){
    read_gps_tle(OPTIONS->tle_constellation_gps_filename, &OPTIONS->nb_gps, OPTIONS->gps_file_name); 
    OPTIONS->n_satellites = OPTIONS->n_satellites + OPTIONS->nb_gps;
  }
  else{
    OPTIONS->nb_gps = 0;
    
  }
  OPTIONS->nb_satellites_not_including_gps = OPTIONS->n_satellites - OPTIONS->nb_gps;

  // Mass of spacecraft
  getline(&line,&len,fp);
  sscanf(line,"%lf", &OPTIONS->mass);
  // Efficiency of the solar cells
  getline(&line, &len, fp);
  strcpy(OPTIONS->opengl_filename_solar_power, "");
  sscanf(line, "%lf %s", &OPTIONS->solar_cell_efficiency, OPTIONS->opengl_filename_solar_power);
  if (strcmp(OPTIONS->opengl_filename_solar_power, "") == 0){
    OPTIONS->opengl_power = 0;
  }
  else{
    OPTIONS->opengl_power = 1;
  }

  //printf("<%s> %f\n", OPTIONS->opengl_filename_solar_power, OPTIONS->solar_cell_efficiency);
  //exitf();
  // Filename of the surface's properties
  getline(&line,&len,fp);
  strtok(line, "\n");  strtok(line, "\r"); 
  next = &line[0];
  strcpy(OPTIONS->filename_surface, "");
  strcpy(OPTIONS->filename_area_attitude_opengl, "");
  int nb_spaces_geo ;
  char *line_copy_geo ;
  line_copy_geo = malloc(sizeof(char) * strlen(line)+1); 
  strcpy(line_copy_geo, line);

  for (nb_spaces_geo=0; line_copy_geo[nb_spaces_geo]; line_copy_geo[nb_spaces_geo]=='.' ? nb_spaces_geo++ : *line_copy_geo++);
  if (nb_spaces_geo > 1){
    sscanf(line, "%s %s", OPTIONS->filename_surface, OPTIONS->filename_area_attitude_opengl);
    OPTIONS->opengl = 1;
  }
  else{
    sscanf(line, "%s", OPTIONS->filename_surface);
    OPTIONS->opengl = 0;
  }
  //newstructure
/*   strcpy(text_location, OPTIONS->dir_input_geometry); */
/*   strcat(text_location, "/"); */
/*   strcpy(OPTIONS->filename_surface, ""); */


  char text_ball_coeff[256];
  //newstructure  
/*   strcpy(text_ball_coeff, OPTIONS->dir_input_geometry); */
/*   strcat(text_ball_coeff, "/"); */
  strcpy(text_ball_coeff, "");
  //newstructure  
  strcat(text_ball_coeff, "ballistic_coefficient");

  if ( (strcmp(OPTIONS->filename_surface, text_ball_coeff) != 0) ){

    if (strcmp(type_orbit_initialisation_temp, "collision_vcm") != 0){
    load_surface( OPTIONS, OPTIONS->filename_surface, nProcs );
    OPTIONS->initialize_geo_with_bstar = 0; 
    }
    else{
      OPTIONS->n_surfaces = 1;
      OPTIONS->n_surfaces_eff = 1;
    }

  }
  else{
    OPTIONS->initialize_geo_with_bstar = 1;
    OPTIONS->n_surfaces = 1;
        OPTIONS->n_surfaces_eff = 1;
  }

  if ( (strcmp(OPTIONS->filename_surface, text_ball_coeff) == 0) && (OPTIONS->solar_cell_efficiency != -1) ){
    printf("***! The solar power can't be computed because you did not indicate a geometry file for the satellite(s) (see section '#SPACECRAFT').!***\n");
    OPTIONS->solar_cell_efficiency  = -1;
  }
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) End of section #SPACECRAFT.\n");
  }

  /* OUTPUT */
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section #OUTPUT.\n");
  }
  rewind(fp);
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#OUTPUT", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  if (feof(fp)){
    printf("***! No section #OUTPUT found in %s. The program will stop. !***\n", filename);
    ierr =  MPI_Finalize();
    exit(0);
  }

  // Names of output files
  getline(&line,&len,fp);
        RemoveSpaces(line);  strtok(line, "\n");  strtok(line, "\r"); 
  sscanf(line, "%s", text);
  if (text[(int)(strlen(text))-1] == '/'){ // if the user wrote a '/' at the very end then remove it
    strcpy(text_temp, text);
    strcpy(text,"");
      strncat(text, &text_temp[0], (int)(strlen(text_temp))-1);
  }

  next = &text[0];

  int find_tild=0;
    find_tild = (int)(strchr(next, '~') - next);
  if (text[0] == '~'){
    strcpy(text_temp, homedir);
    strcat(text_temp, "/");
    strncat(text_temp,next + find_tild+2, next + (int)(strlen(text)));
    strcpy(text,text_temp);
  next = &text[0];
    }

    int find_dir=0;
    find_dir = (int)(strrchr(next, '/') - next);
    if ((find_dir < 0) || (find_dir > 10000)){// no / in the name so the user wants the otuput folder in the directory where he runs the command
    strcpy(text_temp, text);
    strcpy(text,"./");
      strncat(text, &text_temp[0], (int)(strlen(text_temp)));
    find_dir = (int)(strrchr(next, '/') - next);
    //      find_dir = -1;

    }

    strcpy(path_output_run_name_temp, "");
    strncat(path_output_run_name_temp, next, find_dir +1);
    //strncat(text_surface, next, find_file_name);
    strcpy(dir_output_run_name_temp, "");
    strncat(dir_output_run_name_temp, next + find_dir +1, next + (int)(strlen(text)));

    /* printf("<%s>\n",path_output_run_name_temp); */

 
    if (strcmp(dir_output_run_name_temp, "out") == 0 ){
      strcpy(filename_without_extension, "");
      strncat(filename_without_extension, &filename_input_no_path[0], (int)(strrchr(&filename_input_no_path[0], '.') - &filename_input_no_path[0]));
      strcpy(dir_output_run_name_temp, filename_without_extension);

    }
/*     printf("<%s>\n",dir_output_run_name_temp); */
/*     MPI_Finalize();exit(0); */

    if ( strcmp(dir_output_run_name_temp, "now") == 0){
  
  
    if ( iProc == 0 ){

    time_t t_dir_now = time(NULL);
    struct tm tm_dir_now = *gmtime(&t_dir_now);    
    int year_dir_now = tm_dir_now.tm_year + 1900;
    int month_dir_now = tm_dir_now.tm_mon + 1;
    int day_dir_now = tm_dir_now.tm_mday;
    int hour_dir_now = tm_dir_now.tm_hour;
    int minute_dir_now = tm_dir_now.tm_min;
    int second_dir_now = tm_dir_now.tm_sec;

    char year_dir_now_str[15];
    sprintf(year_dir_now_str, "%d", year_dir_now);
    char month_dir_now_str[15];
    sprintf(month_dir_now_str, "%d", month_dir_now);
    char day_dir_now_str[15];
    sprintf(day_dir_now_str, "%d", day_dir_now);
    char hour_dir_now_str[15];
    sprintf(hour_dir_now_str, "%d", hour_dir_now);
    char minute_dir_now_str[15];
    sprintf(minute_dir_now_str, "%d", minute_dir_now);
    char second_dir_now_str[15];
    sprintf(second_dir_now_str, "%d", second_dir_now);
    
    char date_dir_now[256];
    strcpy(date_dir_now,year_dir_now_str);
    strcat(date_dir_now, "-");
    strcat(date_dir_now, month_dir_now_str);
    strcat(date_dir_now, "-");
    strcat(date_dir_now, day_dir_now_str);
    strcat(date_dir_now, " ");
    strcat(date_dir_now, hour_dir_now_str);
    strcat(date_dir_now, ":");
    strcat(date_dir_now, minute_dir_now_str);
    strcat(date_dir_now, ":");
    strcat(date_dir_now, second_dir_now_str);

	str2et_c(date_dir_now, &et_text_output_dir_now); // don't worry about it

	if (nProcs > 1){
    	for (ccc = 1; ccc < nProcs; ccc++){
    MPI_Send(&et_text_output_dir_now, 1, MPI_DOUBLE, ccc, 0, MPI_COMM_WORLD);
	}
	}
    }
    else{
      if (nProcs > 1){
      MPI_Recv(&et_text_output_dir_now, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
    	       MPI_STATUS_IGNORE);
      //      printf("Iproc %d has value %f\n", iProc, et_text_output_dir_now);
      }
    }

	et2utc_c(et_text_output_dir_now, "ISOC", 0, 255, text_output_dir_now); // don't worry about it

	strcpy(dir_output_run_name_temp, "");
	strncat(dir_output_run_name_temp, &text_output_dir_now[5], 2);
	strncat(dir_output_run_name_temp, &text_output_dir_now[8], 2);
	strncat(dir_output_run_name_temp, &text_output_dir_now[2], 2);
	strcat(dir_output_run_name_temp, "T");
	strncat(dir_output_run_name_temp, &text_output_dir_now[11], 2);
	strncat(dir_output_run_name_temp, &text_output_dir_now[14], 2);
	strncat(dir_output_run_name_temp, &text_output_dir_now[17], 2);
    
    MPI_Barrier(MPI_COMM_WORLD);


  }
/*             printf("<%s>\n",path_output_run_name_temp); */
/*         printf("<%s>\n",dir_output_run_name_temp); */
/*         MPI_Finalize();exit(0); */
  for (sss = 0; sss<OPTIONS->n_satellites; sss++){
    strcpy(OPTIONS->name_sat[sss],"");
    if ( (OPTIONS->nb_gps > 0 ) &&  sss >= OPTIONS->nb_satellites_not_including_gps ) {
      strcpy( OPTIONS->filename_output[sss], OPTIONS->gps_file_name[sss-OPTIONS->nb_satellites_not_including_gps]);
      strcpy(OPTIONS->name_sat[sss],OPTIONS->filename_output[sss]);
      strcat( OPTIONS->filename_output[sss],".txt" );
    }
    else{
      strcpy(OPTIONS->filename_output[sss], dir_output_run_name_temp);
      sprintf(sat_nb_str, "%d", sss+1);
      strcat(OPTIONS->filename_output[sss], sat_nb_str);
      strcpy(OPTIONS->name_sat[sss],OPTIONS->filename_output[sss]);
      strcat(OPTIONS->filename_output[sss], ".txt");

    }
  }

  // COLLISION
  strcpy(OPTIONS->filename_collision, dir_output_run_name_temp);
  strcat(OPTIONS->filename_collision, "_collision.txt");
  // end of COLLISION  

  // time step of the output
  getline(&line, &len, fp);
  sscanf( line, "%lf", &OPTIONS->dt_output);
  if (OPTIONS->dt_output < OPTIONS->dt){
    if (iProc == 0){
    printf("***! The time step of the output (%.3f seconds) cannot be bigger that the integration time step (%.3f seconds). Therefore, the time step of the output is set to %.3f seconds. !***\n", OPTIONS->dt_output, OPTIONS->dt, OPTIONS->dt);
    }
    OPTIONS->dt_output = OPTIONS->dt;
  }
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) End of section #OUTPUT.\n");
  }



  /* CD Modification */
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section ##CDMOD.\n");
  }

  for (sss = 0; sss<OPTIONS->nb_satellites_not_including_gps; sss++){
    OPTIONS->cd_modification[sss] = 1.0;
  }
  rewind(fp);
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "##CDMOD", text  ) == 0 )  {
      found_eoh = 1;
    }
  }

  if (found_eoh) {
    for (sss = 0; sss<OPTIONS->nb_satellites_not_including_gps; sss++){
      getline(&line,&len,fp);
      sscanf(line,"%lf",&OPTIONS->cd_modification[sss]);
      // printf("sat: %d; mod : %lf\n",sss,OPTIONS->cd_modification[sss]);
    }
  }
  rewind(fp);
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) End of section ##CDMOD.\n");
  }

  /* ORBIT */
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section #ORBIT.\n");
  }
  rewind(fp);
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#ORBIT", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  if (feof(fp)){
    printf("***! No section #ORBIT found in %s. The program will stop. !***\n", filename);
    ierr =  MPI_Finalize();
    exit(0);

  }


  OPTIONS->aaa_mod = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof(int) ); // used if swpc_mod is on
  int isc;
  for (isc  = 0; isc < OPTIONS->nb_satellites_not_including_gps; isc++){ // this will be calculated again if swpc_mod is on
	      OPTIONS->aaa_mod[isc] = 0;
	}

  getline(&line,&len,fp);
  sscanf(line,"%s",OPTIONS->type_orbit_initialisation);

  if ( strcmp(OPTIONS->type_orbit_initialisation, "oe" ) == 0 ){
    if ( strcmp(OPTIONS->filename_surface, text_ball_coeff) == 0 ){
      printf("***! You need to intialize the orbit with a TLE because you chose 'ballistic_coefficient' instead of a geometry file (see section '#SPACECRAFT'). The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    // COE of each satellite
    for (sss = 0; sss<OPTIONS->nb_satellites_not_including_gps; sss++){
      strcpy(name_sat_temp,"");
      getline(&line, &len, fp);
      sscanf(line, "%s", name_sat_temp);
      sscanf(line, "%lf %lf %lf %lf %lf %lf", &OPTIONS->apogee_alt[sss], &OPTIONS->inclination[sss], &OPTIONS->w[sss], &OPTIONS->long_an[sss], &OPTIONS->f[sss], &OPTIONS->eccentricity[sss]);

      if (name_sat_temp[0] == '\0'){
	printf("Are you sure you correctly entered orbital elements for each of the %d spacecrafts? yes/no \n", OPTIONS->nb_satellites_not_including_gps);
	scanf("%s",answer_coe);
	while ( ( strcmp(answer_coe,"no") != 0 ) && ( strcmp(answer_coe,"yes")  != 0 ) ) {
	  printf("Please type 'yes' or 'no': \n");
	  scanf("%s",answer_coe);
	}
	if (strcmp(answer_coe,"no") == 0){
	  printf("Enter orbital elements for each spacecraft and run the program again! \n");
	  ierr =  MPI_Finalize();
	  exit(0);
	}
	if (strcmp(answer_coe,"yes") == 0){
	  sss = OPTIONS->nb_satellites_not_including_gps;
	}
      }
    }
  }


  if ( strcmp(OPTIONS->type_orbit_initialisation, "state_ecef" ) == 0 ){
    if ( strcmp(OPTIONS->filename_surface, text_ball_coeff) == 0 ){
      printf("***! You need to intialize the orbit with a TLE because you chose 'ballistic_coefficient' instead of a geometry file (see section '#SPACECRAFT'). The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    // COE of each satellite
    for (sss = 0; sss<OPTIONS->nb_satellites_not_including_gps; sss++){
      strcpy(name_sat_temp,"");
      getline(&line, &len, fp);
      sscanf(line, "%s", name_sat_temp);
      sscanf(line, "(%lf; %lf; %lf) (%lf; %lf; %lf)", &OPTIONS->x_ecef[sss], &OPTIONS->y_ecef[sss], &OPTIONS->z_ecef[sss], &OPTIONS->vx_ecef[sss], &OPTIONS->vy_ecef[sss], &OPTIONS->vz_ecef[sss]);


      if (name_sat_temp[0] == '\0'){
	printf("Are you sure you correctly entered orbital elements for each of the %d spacecrafts? yes/no \n", OPTIONS->nb_satellites_not_including_gps);
	scanf("%s",answer_coe);
	while ( ( strcmp(answer_coe,"no") != 0 ) && ( strcmp(answer_coe,"yes")  != 0 ) ) {
	  printf("Please type 'yes' or 'no': \n");
	  scanf("%s",answer_coe);
	}
	if (strcmp(answer_coe,"no") == 0){
	  printf("Enter orbital elements for each spacecraft and run the program again! \n");
	  ierr =  MPI_Finalize();
	  exit(0);
	}
	if (strcmp(answer_coe,"yes") == 0){
	  sss = OPTIONS->nb_satellites_not_including_gps;
	}
      }
    }
  }
  if ( strcmp(OPTIONS->type_orbit_initialisation, "state_eci" ) == 0 ){
    if ( strcmp(OPTIONS->filename_surface, text_ball_coeff) == 0 ){
      printf("***! You need to intialize the orbit with a TLE because you chose 'ballistic_coefficient' instead of a geometry file (see section '#SPACECRAFT'). The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    // COE of each satellite
    for (sss = 0; sss<OPTIONS->nb_satellites_not_including_gps; sss++){
      strcpy(name_sat_temp,"");
      getline(&line, &len, fp);
      sscanf(line, "%s", name_sat_temp);
      RemoveSpaces(line);  strtok(line, "\n");  strtok(line, "\r"); 
      sscanf(line, "(%lf;%lf;%lf)(%lf;%lf;%lf)", &OPTIONS->x_eci[sss], &OPTIONS->y_eci[sss], &OPTIONS->z_eci[sss], &OPTIONS->vx_eci[sss], &OPTIONS->vy_eci[sss], &OPTIONS->vz_eci[sss]);

      if (name_sat_temp[0] == '\0'){
	printf("Are you sure you correctly entered orbital elements for each of the %d spacecrafts? yes/no \n", OPTIONS->nb_satellites_not_including_gps);
	scanf("%s",answer_coe);
	while ( ( strcmp(answer_coe,"no") != 0 ) && ( strcmp(answer_coe,"yes")  != 0 ) ) {
	  printf("Please type 'yes' or 'no': \n");
	  scanf("%s",answer_coe);
	}
	if (strcmp(answer_coe,"no") == 0){
	  printf("Enter orbital elements for each spacecraft and run the program again! \n");
	  ierr =  MPI_Finalize();
	  exit(0);
	}
	if (strcmp(answer_coe,"yes") == 0){
	  sss = OPTIONS->nb_satellites_not_including_gps;
	}
      }
    }
  }

  if ( strcmp(OPTIONS->type_orbit_initialisation, "deployment" ) == 0 ){
    if ( strcmp(OPTIONS->filename_surface, text_ball_coeff) == 0 ){
      printf("***! You need to intialize the orbit with a TLE because you chose 'ballistic_coefficient' instead of a geometry file (see section '#SPACECRAFT'). The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    // COE of each satellite
    for (sss = 0; sss<OPTIONS->nb_satellites_not_including_gps; sss++){
      strcpy(name_sat_temp,"");
      sscanf(line, "%s", name_sat_temp);
      getline(&line, &len, fp);
      sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf", &OPTIONS->deployment_module_orbital_elements[sss][0], &OPTIONS->deployment_module_orbital_elements[sss][1], &OPTIONS->deployment_module_orbital_elements[sss][2], &OPTIONS->deployment_module_orbital_elements[sss][3], &OPTIONS->deployment_module_orbital_elements[sss][4], &OPTIONS->deployment_module_orbital_elements[sss][5], &OPTIONS->deployment_speed[sss], &OPTIONS->deployment_angle[sss]);

      if (name_sat_temp[0] == '\0'){
	printf("Are you sure you correctly entered orbital elements for each of the %d spacecrafts? yes/no \n", OPTIONS->nb_satellites_not_including_gps);
	scanf("%s",answer_coe);
	while ( ( strcmp(answer_coe,"no") != 0 ) && ( strcmp(answer_coe,"yes")  != 0 ) ) {
	  printf("Please type 'yes' or 'no': \n");
	  scanf("%s",answer_coe);
	}
	if (strcmp(answer_coe,"no") == 0){
	  printf("Enter orbital elements for each spacecraft and run the program again! \n");
	  ierr =  MPI_Finalize();
	  exit(0);
	}
	if (strcmp(answer_coe,"yes") == 0){
	  sss = OPTIONS->nb_satellites_not_including_gps;
	}
      }
    }
  }


  OPTIONS->use_sgp4 = 0;
  if ( (strcmp(OPTIONS->type_orbit_initialisation, "tle" ) == 0 ) || (strcmp(OPTIONS->type_orbit_initialisation, "tle_sgp4" ) == 0 )){
    if (strcmp(OPTIONS->type_orbit_initialisation, "tle_sgp4" ) == 0 ){
      OPTIONS->use_sgp4 = 1;
    }
    //    if (iProc == 0){
    sss = 0;
    while (sss<OPTIONS->nb_satellites_not_including_gps){
        
      getline(&line,&len, fp);
	sscanf(line, "%s", OPTIONS->tle_initialisation[sss]);

      if ( ( OPTIONS->tle_initialisation[sss][0] != '\0' ) && ( OPTIONS->tle_initialisation[sss][0] != '#' ) ){
	OPTIONS->one_file_with_tle_of_all_sc = 0;
      }
      else{ // if the user put only one file name, then this file needs to include the TLEs for all sc
	OPTIONS->one_file_with_tle_of_all_sc = 1;
		sss = OPTIONS->nb_satellites_not_including_gps;
	/* printf("Are you sure you correctly entered TLE file names for each of the %d spacecraft? yes/no \n", OPTIONS->nb_satellites_not_including_gps); */
	/* scanf("%s",answer_coe); */
	/* while ( ( strcmp(answer_coe,"no") != 0 ) && ( strcmp(answer_coe,"yes")  != 0 ) ) { */
	/*   printf("Please type 'yes' or 'no': \n"); */
	/*   scanf("%s",answer_coe); */
	/* } */
	/* if (strcmp(answer_coe,"no") == 0){ */
	/*   printf("Enter TLE file names for each spacecraft and run the program again! \n"); */
	/*   ierr =  MPI_Finalize(); */
	/*   exit(0); */
	/* } */
	/* if (strcmp(answer_coe,"yes") == 0){ */
	/*   sss = OPTIONS->nb_satellites_not_including_gps; */
	/* } */
      }
      sss = sss + 1;
    }
    for (ii = 0; ii < OPTIONS->nb_satellites_not_including_gps; ii++){
      if (ii > 0){
	epoch_sat_previous = epoch_sat;
      }

    /*** Convert the TLEs into inertial state (postion, velocity) ***/
	    /* Read in the next two lines from the text file that contains the TLEs. */
	  //newstructure
/* 	  strcpy(sat_tle_file_temp, OPTIONS->dir_input_tle); */
/* 	  strcat(sat_tle_file_temp,"/"); */
	  strcpy(sat_tle_file_temp,"");

	  //newstructure
	    if ( OPTIONS->one_file_with_tle_of_all_sc == 0) { // one tle per file
	      strcat(sat_tle_file_temp,OPTIONS->tle_initialisation[ii]);
	      sat_tle_file = fopen(sat_tle_file_temp,"r");
	      if (sat_tle_file == NULL){
		printf("!***!\nThe TLE file:\n%s\ncould not be found. The program will stop.\n!***!\n", sat_tle_file_temp); MPI_Finalize(); exit(0);
	      }
	    }
	    else{ // all tles in one file
	      strcat(sat_tle_file_temp,OPTIONS->tle_initialisation[0]);

	      sat_tle_file = fopen(sat_tle_file_temp,"r");
	      if (sat_tle_file == NULL){
		printf("!***!\nThe TLE file:\n%s\ncould not be found. The program will stop.\n!***!\n", sat_tle_file_temp); MPI_Finalize(); exit(0);
	      }
		      
	      for (bbb = 0; bbb < ii; bbb++){ // skip the TLEs of the sc before the current sc
		getline( &line_temp, &len, sat_tle_file ); 
		getline( &line_temp, &len, sat_tle_file );
	      }
	    }
	  
	    // First line
	    getline( &line_temp, &len, sat_tle_file );
	    lineln = strlen(line_temp)-1;
	    SpiceChar line_tle[2][lineln];
	    for (pp = 0; pp<lineln-1 ; pp++)
	      line_tle[0][pp] = line_temp[pp];
	    line_tle[0][ lineln-1 ] = '\0';
	    // Second line
	    getline( &line_temp, &len, sat_tle_file );
	    for (pp = 0; pp<lineln-1 ; pp++)
	      line_tle[1][pp] = line_temp[pp];
	    line_tle[1][ lineln-1 ] = '\0';

	    fclose(sat_tle_file);
	    // Convert the elements of the TLE into "elems" and "epoch" that can be then read by the SPICE routine ev2lin_ to convert into the inertial state

	    	    getelm_c( frstyr, lineln, line_tle, &epoch_sat, elems );

		    // compute the oldest TLE -> to download the inputs for the drag model (we don't care about GPS sc since they don't feel any drag)
	    if (ii == 0){
	      epoch_sat_previous = epoch_sat;
	      et_most_recent_tle_epoch = epoch_sat;
	      OPTIONS->et_oldest_tle_epoch = epoch_sat;
	    }
	    else{
	      if ( epoch_sat < OPTIONS->et_oldest_tle_epoch ){
		OPTIONS->et_oldest_tle_epoch = epoch_sat;
				OPTIONS->which_sc_oldest_tle_epoch = ii;
	      }
	      if (  epoch_sat > et_most_recent_tle_epoch){
		et_most_recent_tle_epoch = epoch_sat;
	      }
	    
	    }
	    //	    CONSTELLATION->spacecraft[ii][eee].et = epoch_sat;


	  

    }
    
  et2utc_c(OPTIONS->et_oldest_tle_epoch, "ISOC" ,6 ,255 , OPTIONS->oldest_tle_epoch);
    if (et_start < et_most_recent_tle_epoch){

    et2utc_c(et_most_recent_tle_epoch, "ISOC" ,6 ,255 , most_recent_tle_epoch);
    printf("***! The initial epoch set in section #TIME of the input file needs to be more recent than the TLE epochs. The most recent TLE epoch is %s so the initial epoch needs to be more recent than this. The program will stop. !***\n", most_recent_tle_epoch); MPI_Finalize();exit(0);
  }

  }
  else{ //if tle are not used then itnialize et_oldest_tle_epoch and oldest_tle_epoch to the initial epoch of the constellation, unless the user is computing collisoins and using a VCM as an input file. In this case et_oldest_tle_epoch is the oldest VMC epoch time among the 2 VCMs. This is done later, after reading in read_vcm the epoch of the VCMs
    strcpy(OPTIONS->oldest_tle_epoch, OPTIONS->initial_epoch);
    double et_initial_epoch_temp;
    str2et_c( OPTIONS->initial_epoch, &et_initial_epoch_temp);
    OPTIONS->et_oldest_tle_epoch = et_initial_epoch_temp;
  }
  //  printf("<%s>\n", OPTIONS->oldest_tle_epoch);


  if ( strcmp(OPTIONS->type_orbit_initialisation, "collision" ) == 0 ){

    char filename_input_collision_temp[256];
    getline(&line, &len, fp);     
    //newstructure
/*     strcpy( OPTIONS->filename_input_collision, OPTIONS->dir_input_collision); */
/*     strcat( OPTIONS->filename_input_collision, "/"); */
    strcpy( OPTIONS->filename_input_collision, "");
    //newstructure
    strcpy(filename_input_collision_temp, "");
    sscanf(line, "%s", filename_input_collision_temp);
    strcat(OPTIONS->filename_input_collision, filename_input_collision_temp);
    OPTIONS->coll_vcm = 0;
    ini_collision( OPTIONS, iProc );
  }

  if ( strcmp(OPTIONS->type_orbit_initialisation, "collision_vcm" ) == 0 ){
    if (OPTIONS->nb_satellites_not_including_gps != 2){
	print_error(iProc, "The number of satellites must be 2 if the option 'collision_vcm' was selected");
    }
    char filename_input_collision_temp[256];
    getline(&line, &len, fp);     
    strcpy( OPTIONS->filename_input_collision_vcm1, "");
    strcpy(filename_input_collision_temp, "");
    sscanf(line, "%s", filename_input_collision_temp);
    strcat(OPTIONS->filename_input_collision_vcm1, filename_input_collision_temp);

    getline(&line, &len, fp);     
    strcpy( OPTIONS->filename_input_collision_vcm2, "");
    strcpy(filename_input_collision_temp, "");
    sscanf(line, "%s", filename_input_collision_temp);
    strcat(OPTIONS->filename_input_collision_vcm2, filename_input_collision_temp);

    getline(&line, &len, fp);     
    strcpy( OPTIONS->filename_input_collision_cdm, "");
    strcpy(filename_input_collision_temp, "");
    sscanf(line, "%s", filename_input_collision_temp);
    strcat(OPTIONS->filename_input_collision_cdm, filename_input_collision_temp);

    getline( &line, &len, fp);
    OPTIONS->nb_ensembles = 0;
    OPTIONS->min_dist_close_approach = 0;
    OPTIONS->min_dist_collision = 0;
    sscanf(line,"%lf %lf %d", &OPTIONS->min_dist_close_approach, &OPTIONS->min_dist_collision, &OPTIONS->nb_ensembles);

    OPTIONS->min_dist_close_approach = OPTIONS->min_dist_close_approach / 1000.;  // !!!!!!!!! input is m but need km
    OPTIONS->min_dist_collision = OPTIONS->min_dist_collision / 1000.;  

    //    printf("%f %f %d\n", OPTIONS->min_dist_close_approach*1000, OPTIONS->min_dist_collision*1000, OPTIONS->nb_ensembles);

    OPTIONS->coll_vcm = 1;

    read_vcm(OPTIONS->filename_input_collision_vcm1, OPTIONS, 0);
    read_vcm(OPTIONS->filename_input_collision_vcm2, OPTIONS, 1);
    read_cdm(OPTIONS->filename_input_collision_cdm, OPTIONS);
    //    exitf();
// et_oldest_tle_epoch has already been calculated (equal to the epcoh start of the main input file) but overwrite it to take into account the VCM epochs -> et_oldest_tle_epoch is the oldest of the two epoch
OPTIONS->et_oldest_tle_epoch = OPTIONS->et_vcm[0];
if  (OPTIONS->et_vcm[1] < OPTIONS->et_vcm[0]){
OPTIONS->et_oldest_tle_epoch =  OPTIONS->et_vcm[1];
  OPTIONS->which_sc_oldest_tle_epoch = 1;
}
  et2utc_c(OPTIONS->et_oldest_tle_epoch, "ISOC" ,6 ,255 , OPTIONS->oldest_tle_epoch);

//etprint(OPTIONS->et_oldest_tle_epoch, "OPTIONS->et_oldest_tle_epoch");exitf();
    ini_collision( OPTIONS, iProc );

  }


  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) End of section #ORBIT.\n");
  }


  /** ENSEMBLES ON COE**/ 
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section ##ENSEMBLES_COE.\n");
  }
  if ( (strcmp(OPTIONS->type_orbit_initialisation, "collision" ) != 0 ) && (strcmp(OPTIONS->type_orbit_initialisation, "collision_vcm" ) != 0 )){
    rewind(fp);
    found_eoh = 0;
    while ( found_eoh == 0 && !feof(fp)) {
      getline(&line, &len, fp);
      sscanf(line, "%s", text);
      if (  strcmp( "##ENSEMBLES_COE", text  ) == 0 )  {
	found_eoh = 1;
      }
    }
    if (feof(fp)){
      OPTIONS->nb_ensembles = 0;
    }
    else {
      // COE of each satellite
      //      char useless[256];
      getline(&line, &len, fp);
      OPTIONS->nb_ensembles = 0;
      sscanf(line,"%d", &OPTIONS->nb_ensembles);
      if ( strcmp(OPTIONS->type_orbit_initialisation, "oe" ) == 0 ){
	for (sss = 0; sss<OPTIONS->nb_satellites_not_including_gps; sss++){
	  getline(&line, &len, fp);
	  sscanf(line, "%lf %lf %lf %lf %lf %lf",  &OPTIONS->apogee_alt_sigma[sss], &OPTIONS->inclination_sigma[sss], &OPTIONS->w_sigma[sss], &OPTIONS->long_an_sigma[sss], &OPTIONS->f_sigma[sss], &OPTIONS->eccentricity_sigma[sss]);
	  OPTIONS->coe_to_ensemble[sss][5] = 0;
	  OPTIONS->coe_to_ensemble[sss][4] = 0;
	  OPTIONS->coe_to_ensemble[sss][3] = 0;
	  OPTIONS->coe_to_ensemble[sss][2] = 0;
	  OPTIONS->coe_to_ensemble[sss][1] = 0;
	  OPTIONS->coe_to_ensemble[sss][0] = 0;
      
	  if (OPTIONS->apogee_alt_sigma[sss] != 0){
	    OPTIONS->coe_to_ensemble[sss][5] = 1;
	  }
	  if (OPTIONS->inclination_sigma[sss] !=0){
	    OPTIONS->coe_to_ensemble[sss][4] = 1;
	  }
	  if (OPTIONS->w_sigma[sss] !=0){
	    OPTIONS->coe_to_ensemble[sss][3] = 1;
	  }
	  if (OPTIONS->long_an_sigma[sss] !=0){
	    OPTIONS->coe_to_ensemble[sss][2] = 1;
	  }
	  if (OPTIONS->f_sigma[sss] !=0){
	    OPTIONS->coe_to_ensemble[sss][1] = 1;
	  }
	  if (OPTIONS->eccentricity_sigma[sss] !=0){
	    OPTIONS->coe_to_ensemble[sss][0] = 1;
	  }
      

	  if (line[0] == '\0'){
	    printf("Are you sure you correctly entered ensembles of orbital elements for each of the %d spacecrafts? yes/no \n", OPTIONS->nb_satellites_not_including_gps);
	    scanf("%s",answer_coe);
	    while ( ( strcmp(answer_coe,"no") != 0 ) && ( strcmp(answer_coe,"yes")  != 0 ) ) {
	      printf("Please type 'yes' or 'no': \n");
	      scanf("%s",answer_coe);
	    }
	    if (strcmp(answer_coe,"no") == 0){
	      printf("Enter ensembles of orbital elements for each spacecraft and run the program again! \n");
	      ierr =  MPI_Finalize();
	      exit(0);
	    }
	    if (strcmp(answer_coe,"yes") == 0){
	      sss = OPTIONS->nb_satellites_not_including_gps;
	    }
	  }  
	}
	if (nProcs > 0){
	  OPTIONS->nb_ensemble_per_proc = (int)(OPTIONS->nb_ensembles / nProcs);
	}
      }

      if ( strcmp(OPTIONS->type_orbit_initialisation, "state_eci" ) == 0 ){
	for (sss = 0; sss<OPTIONS->nb_satellites_not_including_gps; sss++){
	  getline(&line, &len, fp);
	  sscanf(line, "(%lf; %lf; %lf) (%lf; %lf; %lf)", &OPTIONS->x_eci_sigma[sss], &OPTIONS->y_eci_sigma[sss], &OPTIONS->z_eci_sigma[sss], &OPTIONS->vx_eci_sigma[sss], &OPTIONS->vy_eci_sigma[sss], &OPTIONS->vz_eci_sigma[sss]);

	  if (line[0] == '\0'){
	    printf("Are you sure you correctly entered ensembles of state eci for each of the %d spacecrafts? yes/no \n", OPTIONS->nb_satellites_not_including_gps);
	    scanf("%s",answer_coe);
	    while ( ( strcmp(answer_coe,"no") != 0 ) && ( strcmp(answer_coe,"yes")  != 0 ) ) {
	      printf("Please type 'yes' or 'no': \n");
	      scanf("%s",answer_coe);
	    }
	    if (strcmp(answer_coe,"no") == 0){
	      printf("Enter ensembles of state eci for each spacecraft and run the program again! \n");
	      ierr =  MPI_Finalize();
	      exit(0);
	    }
	    if (strcmp(answer_coe,"yes") == 0){
	      sss = OPTIONS->nb_satellites_not_including_gps;
	    }
	  }  
	}

	if (nProcs > 0){
	  OPTIONS->nb_ensemble_per_proc = (int)(OPTIONS->nb_ensembles / nProcs);
	}
      }
    }
  }

  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) End of section ##ENSEMBLES_COE.\n");
  }


  /* KALMAN */
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section #KALMAN.\n");
  }
  rewind(fp);
  found_eoh = 0;
  OPTIONS->use_kalman = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#KALMAN", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  if (!feof(fp)){
    
  getline(&line, &len, fp);
  strcpy(OPTIONS->filename_kalman_init, "");
  sscanf(line,"%s", OPTIONS->filename_kalman_init);
  // the observation file has probably not the same times as the propagation time steps, which means that when using the kalman filter within SpOCK, the time step is not constant. This has to be taken into account for the interpolation of the attitude, Ap, F10.7, ...
  compute_time_interpo(OPTIONS);
  OPTIONS->use_kalman = 1;
  }

  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) End of section #KALMAN.\n");
  }


  /* ATTITUDE */
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section #ATTITUDE.\n");
  }

  rewind(fp);
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#ATTITUDE", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  if (feof(fp)){ // if the user does not include a section #ATTITUDE, then the attitude is set to nadir pointing by default
    strcpy(OPTIONS->attitude_profile, "nadir");
    //nb_time_steps(&OPTIONS->nb_time_steps, OPTIONS->initial_epoch, OPTIONS->final_epoch, OPTIONS->dt);
    if (OPTIONS->use_kalman != 1){ // if KF is used, then the nb of time steps has already been computed
      nb_time_steps(&OPTIONS->nb_time_steps, OPTIONS->et_oldest_tle_epoch, OPTIONS->final_epoch, OPTIONS->dt);// tle
    }
    load_attitude( OPTIONS, OPTIONS->attitude_profile, fp, nProcs, ang_velo, iDebugLevel, iProc );
  }
  else{
    getline(&line,&len,fp);
    sscanf(line,"%s", text);
    strcpy(OPTIONS->attitude_profile, text);
    //    nb_time_steps(&OPTIONS->nb_time_steps, OPTIONS->initial_epoch, OPTIONS->final_epoch, OPTIONS->dt);
    if (OPTIONS->use_kalman != 1){ // if KF is used, then the nb of time steps has already been computed
      nb_time_steps(&OPTIONS->nb_time_steps, OPTIONS->et_oldest_tle_epoch, OPTIONS->final_epoch, OPTIONS->dt);// tle
    }

    // if the user put a vector to represent the angular velocity of the satellite
    if (OPTIONS->attitude_profile[0] == '('){
    RemoveSpaces(line);  strtok(line, "\n");  strtok(line, "\r"); 
      sscanf(line, "(%lf;%lf;%lf)(%lf;%lf;%lf)", &ang_velo[0], &ang_velo[1], &ang_velo[2], &ang_velo[3], &ang_velo[4], &ang_velo[5]); // first 3 is the inital atitude, last 3 is the (constant) angular velocity. in degree and degree per sec
      //      sscanf(line, "(%lf;%lf;%lf)", &ang_velo[0], &ang_velo[1], &ang_velo[2]); // in degree per sec
      strcpy(OPTIONS->attitude_profile, "angular_velocity");
    }

    load_attitude( OPTIONS, OPTIONS->attitude_profile, fp, nProcs, ang_velo, iDebugLevel, iProc );

  }
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section ##ENSEMBLES_ATTITUDE.\n");
  }


  /** ENSEMBLES_ATTITUDE **/
  if ( ( strcmp( OPTIONS->attitude_profile, "ensemble_angular_velocity" ) != 0 ) && ( strcmp( OPTIONS->attitude_profile, "ensemble_initial_attitude" ) != 0) ) { // if you run ensembles on the initial angular velocity then this has been calculated in the section #ATTITUDE
    rewind(fp);
    found_eoh = 0;
    while ( found_eoh == 0 && !feof(fp)) {
      getline(&line, &len, fp);
      sscanf(line, "%s", text);
      if (  strcmp( "##ENSEMBLES_ATTITUDE", text  ) == 0 )  {
	found_eoh = 1;
      }
    }
    if (feof(fp)){
      OPTIONS->nb_ensembles_attitude = 0;
    }
    else{
      getline(&line, &len, fp);
      OPTIONS->nb_ensembles_attitude = 0;
      sscanf(line, "%d", &OPTIONS->nb_ensembles_attitude);
      if ( ( OPTIONS->nb_ensembles_attitude > 0 ) && ( strcmp( OPTIONS->attitude_profile, "ensemble_angular_velocity") != 0 ) ){
	if ((iProc == 0) & ( iDebugLevel >= 2 )){
	  printf("--- (load_options) Beginning of section ###ENSEMBLES_ANGULAR_VELOCITY.\n");
	}
	rewind(fp);
	found_eoh = 0;
	while ( found_eoh == 0 && !feof(fp)) {
	  getline(&line, &len, fp);
	  sscanf(line, "%s", text);
	  if (  strcmp( "###ENSEMBLES_ANGULAR_VELOCITY", text  ) == 0 )  {
	    found_eoh = 1;
	  }
	}
	if (feof(fp)){
	  printf("***! You chose to run ensembles from the section ###ENSEMBLES_ANGULAR_VELOCITY but no such section was found in %s. The program will stop. !***\n", filename); MPI_Finalize(); exit(0);
	}
	getline(&line, &len, fp);
	sscanf(line, "%lf", &OPTIONS->attitude_reset_delay);
	getline(&line, &len, fp);
	sscanf(line, "(%lf; %lf; %lf)", &OPTIONS->pitch_sigma_angular_velocity_ensemble, &OPTIONS->roll_sigma_angular_velocity_ensemble, &OPTIONS->yaw_sigma_angular_velocity_ensemble);
	
	if ((iProc == 0) & ( iDebugLevel >= 2 )){
	  printf("--- (load_options) End of section ###ENSEMBLES_ANGULAR_VELOCITY.\n");
	}

      }
    }
  }
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) End of section ##ENSEMBLES_ATTITUDE.\n");
  }

  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) End of section #ATTITUDE.\n");
  }
    if ((OPTIONS->nb_ensembles_attitude > 0) && (OPTIONS->file_is_quaternion == 1)){
	print_error(iProc, "Ensembles on attitude can't be run if the attitude is set with quaternions");
      }




  /* FORCES */
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section #FORCES.\n");
  }
  rewind(fp);
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#FORCES", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  if (feof(fp)){
    printf("***! No section #FORCES found in %s. The program will stop. !***\n", filename);
    ierr =  MPI_Finalize();
    exit(0);
  }
  getline(&line, &len, fp);
  RemoveSpaces(line);  strtok(line, "\n");  strtok(line, "\r");
  char gravity_map_str[100];
  strcpy(gravity_map_str, "");
  OPTIONS->gravity_map = 0;
  sscanf(line,"%lf %s", &OPTIONS->degree, gravity_map_str);
  if (strcmp(gravity_map_str, "map") == 0){
      OPTIONS->gravity_map = 1;
    }
  OPTIONS->order = OPTIONS->degree; 

  
  getline(&line,&len,fp);
  if(strstr(line, "drag") != NULL) {
    OPTIONS->include_drag = 1;
  }
  else{
    OPTIONS->include_drag = 0;
  }
  if(strstr(line, "solar_pressure") != NULL) {
    OPTIONS->include_solar_pressure = 1;
  }
  else{
    OPTIONS->include_solar_pressure = 0;
  }
  if(strstr(line, "earth_pressure") != NULL) {
    OPTIONS->include_earth_pressure = 1;
  }
  else{
    OPTIONS->include_earth_pressure = 0;
  }


  if(strstr(line, "sun_gravity") != NULL) {
    OPTIONS->include_sun = 1;
  }
  else{
    OPTIONS->include_sun = 0;
  }
  if(strstr(line, "moon_gravity") != NULL) {
    OPTIONS->include_moon = 1;
  }
  else{
    OPTIONS->include_moon = 0;
  }
  //  printf("%f %f %d %d %d %d\n", OPTIONS->order, OPTIONS->degree,OPTIONS->include_drag, OPTIONS->include_solar_pressure, OPTIONS->include_sun, OPTIONS->include_moon);exit(0);
  if ( (strcmp(OPTIONS->filename_surface, text_ball_coeff)  == 0) && ((OPTIONS->include_solar_pressure == 1) || (OPTIONS->include_earth_pressure == 1)) ){
    printf("***! The earth and solar radiation pressure force can't be computed because you did not indicate a geometry file for the satellite(s) (see section '#SPACECRAFT'). !***\n");
    OPTIONS->include_solar_pressure = 0;
    OPTIONS->include_earth_pressure = 0;
  }

  getline(&line,&len,fp);
  if ( OPTIONS->include_drag == 1){
    OPTIONS->nb_ensembles_density = 0;
    OPTIONS->nb_ensemble_density_per_proc = 0;

    sscanf(line, "%s",OPTIONS->format_density_driver);

    if ( strcmp( OPTIONS->format_density_driver, "dynamic"  ) == 0 ){     

      OPTIONS->Ap = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) ); // "* 2.0" because of the Runge Kunta orfer 4 method

      OPTIONS->Ap_hist = malloc(  7 * sizeof(double *) ); // historical ap
      for (hhh = 0; hhh < 7; hhh++){
	OPTIONS->Ap_hist[hhh] = malloc(OPTIONS->nb_time_steps * 2 * sizeof(double ));
      }

      OPTIONS->f107 = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
      OPTIONS->f107A = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );

      if (  OPTIONS->Ap == NULL ){
	printf("***! Could not allow memory space for  OPTIONS->Ap \n. The program will stop. !***\n");
	ierr =  MPI_Finalize();
	exit(0);
      }

      if (  OPTIONS->Ap_hist == NULL ){
	printf("***! Could not allow memory space for  OPTIONS->Ap_hist \n. The program will stop. !***\n");
	ierr =  MPI_Finalize();
	exit(0);
      }

      if (  OPTIONS->f107 == NULL ){
	printf("***! Could not allow memory space for  OPTIONS->f107 \n. The program will stop. !***\n");
	ierr =  MPI_Finalize();
	exit(0);
      }
      if (  OPTIONS->f107A == NULL ){
	printf("***! Could not allow memory space for  OPTIONS->f107A \n. The program will stop. !***\n");
	ierr =  MPI_Finalize();
	exit(0);
      }
      
      getline(&line,&len,fp);
      sscanf(line, "%s",OPTIONS->test_omniweb_or_external_file);

      if (strcmp(OPTIONS->test_omniweb_or_external_file, "omniweb") == 0){ // if the user chose to automatically download F10.7 and Ap from omniweb
      //	str2et_c(OPTIONS->initial_epoch, &et_initial_epoch);
	//	et2utc_c(et_initial_epoch, "ISOC" ,0 ,11 , initial_epoch_wget_temp);
	str2et_c(OPTIONS->final_epoch, &et_final_epoch);
	double et_oldest_tle_epoch_corr = OPTIONS->et_oldest_tle_epoch;
	double et_final_epoch_corr = et_final_epoch;
	if (OPTIONS->et_oldest_tle_epoch > et_final_epoch){ // backward propagation
	  et_oldest_tle_epoch_corr = et_final_epoch;
	  et_final_epoch_corr = OPTIONS->et_oldest_tle_epoch;
	  }
		      
	et2utc_c(et_oldest_tle_epoch_corr, "ISOC" ,0 ,11 , initial_epoch_wget_temp); // tle
	strcpy(initial_epoch_wget, "");
	next_wget_initial = &initial_epoch_wget_temp[0];
	find_wget_initial = (int)(strchr(next_wget_initial, '-') - next_wget_initial);
	strncat(initial_epoch_wget, initial_epoch_wget_temp, 4);
	strncat(initial_epoch_wget, initial_epoch_wget_temp+5, 2);
	strncat(initial_epoch_wget, initial_epoch_wget_temp+8, 2);
	et2utc_c(et_final_epoch_corr, "ISOC" ,0 ,11 , final_epoch_wget_temp);
	strcpy(final_epoch_wget, "");
	next_wget_final = &final_epoch_wget_temp[0];
	find_wget_final = (int)(strchr(next_wget_final, '-') - next_wget_final);
	strncat(final_epoch_wget, final_epoch_wget_temp, 4);
	strncat(final_epoch_wget, final_epoch_wget_temp+5, 2);
	strncat(final_epoch_wget, final_epoch_wget_temp+8, 2);
	strcpy(str_wget, "wget --no-check-certificate  --post-data ");
	strcat(str_wget, "\"activity=retrieve&res=hour&spacecraft=omni2&start_date=");
	strcat(str_wget, initial_epoch_wget);
	strcat(str_wget, "&end_date=");
	strcat(str_wget, final_epoch_wget);
	strcat(str_wget,"&vars=50&scale=Linear&ymin=&ymax=&view=0&charsize=&xstyle=0&ystyle=0&symbol=0&symsize=&linestyle=solid&table=0&imagex=640&imagey=480&color=&back=\"");
	strcat(str_wget, " https://omniweb.gsfc.nasa.gov/cgi/nx1.cgi -O ");
	strcpy(filename_f107, "");
	//newstructure
/* 	strcat(filename_f107, OPTIONS->dir_input_density_msis); */
/* 	strcat(filename_f107, "/"); */
	//newstructure 
	strcat(filename_f107, "f107_" );
	strcat(filename_f107, initial_epoch_wget);
	strcat(filename_f107, "_to_" );
	strcat(filename_f107, final_epoch_wget);
	strcat(filename_f107, ".txt");
	strcat(str_wget, filename_f107);
	strcat(str_wget, " >/dev/null 2>&1");
	if (do_not_download_file_swpc_or_wget !=1 ){
	system(str_wget);
	}
	// get the name of F10.7 now
	strcpy(filename_f107A, "");
	//newstructure 
/* 	strcat(filename_f107A, OPTIONS->dir_input_density_msis); */
/* 	strcat(filename_f107A, "/"); */
	//newstructure 
	strcat(filename_f107A, "f107A_" );
	strcat(filename_f107A, initial_epoch_wget);
	strcat(filename_f107A, "_to_" );
	strcat(filename_f107A, final_epoch_wget);
	strcat(filename_f107A, ".txt");

	// Quick check that the start and end dates are ok
	f107_wget_omniweb = fopen(filename_f107, "r");
	getline(&line, &len, f107_wget_omniweb);
	sscanf(line, "%s", test_wget_omniweb_error);


	if (strcmp(test_wget_omniweb_error, "<H1>") == 0){// if file starts with <H1> then it means that the start date or end date is invalid for omniweb

	  // // go get the correct range for omniweb (end date of omniweb only) 
	  strcpy(correct_end_date_omniweb_temp, "");

	  getline(&line, &len, f107_wget_omniweb);

	  next_wget_final = &line[0];
	  find_wget_final = (int)(strchr(next_wget_final,'-')-next_wget_final);
	  strncat(correct_end_date_omniweb_temp, next_wget_final+find_wget_final+2,8);
	  // // test if start or end date of epoch is more recent than correct_end_date_omniweb_temp
	  // // // convert correct_end_date_omniweb_temp in a correct format to compare to et_initial_epoch and et_final_epoch
	  strcpy(correct_end_date_omniweb, "");
	  strncat(correct_end_date_omniweb, correct_end_date_omniweb_temp, 4);
	  strcat(correct_end_date_omniweb,  "-");
	  strncat(correct_end_date_omniweb, correct_end_date_omniweb_temp+4, 2);
	  strcat(correct_end_date_omniweb,  "-");
	  strncat(correct_end_date_omniweb, correct_end_date_omniweb_temp+6, 2);
	  strcat(correct_end_date_omniweb, "T");
	  strcat(correct_end_date_omniweb, "00:00");
	  str2et_c(correct_end_date_omniweb, &et_correct_end_date_omniweb);
	  // // // compare et_correct_end_date_omniweb to et_initial_epoch and et_final_epoch

	  if  ( ( et_correct_end_date_omniweb < et_oldest_tle_epoch_corr ) || ( et_correct_end_date_omniweb < et_final_epoch_corr ) ){ // if start or end date of epoch is older than correct_end_date_omniweb then the user can not run msis dynamic with omniweb because he/she uses the inital or final epoch of the constellation (including TLEs epoch if TLEs are used to initialize the trajectories) older than the most recent f107 available in omniweb
	    printf("***! Omniweb data (used for the thermospheric density model) is not available for this epoch. You need to choose epoch times more recent than: %s. Or you can use data from the Space Weather Prediction Center (SWPC), which is more up to date. For this, replace 'omniweb' by 'swpc' in the main input file.\nThe program will stop. !***\n", correct_end_date_omniweb);
	    ierr =  MPI_Finalize();
	    exit(0);
	  }
	} 
	fclose(f107_wget_omniweb);
	// Download the omniweb f10.7 that starts 40.5 days before initial_epoch_wget and ends 40.5 days after final_epoch_wget to compute F10.7 81-day average.
	// There should be data for 40.5 days before initial_epoch_wget. But there might not be data for 40.5 days after final_epoch_wget. In that case, upload until the most recent date in omniweb. F10.7 average will be computed a this shorter period (instead of 81 days)
	// // et_initial_epoch minus 40.5 days
	et2utc_c(et_oldest_tle_epoch_corr-40.5*24*3600, "ISOC" ,0 ,11 , initial_epoch_wget_minus40_5days_temp);
	strcpy(initial_epoch_wget_minus40_5days, "");
	next_wget_initial = &initial_epoch_wget_minus40_5days_temp[0];
	find_wget_initial = (int)(strchr(next_wget_initial, '-') - next_wget_initial);
	strncat(initial_epoch_wget_minus40_5days, initial_epoch_wget_minus40_5days_temp, 4);
	strncat(initial_epoch_wget_minus40_5days, initial_epoch_wget_minus40_5days_temp+5, 2);
	strncat(initial_epoch_wget_minus40_5days, initial_epoch_wget_minus40_5days_temp+8, 2);

	et2utc_c(et_final_epoch_corr+41.5*24*3600, "ISOC" ,0 ,11 , final_epoch_wget_plus40_5days_temp); // 41.5 to make sure we cover up to et_final + 40.5 days (because the f107 file stops at 23:00 so problem if et_final _40.5 ends between 23:00 and 23:59)
	strcpy(final_epoch_wget_plus40_5days, "");
	next_wget_final = &final_epoch_wget_plus40_5days_temp[0];
	find_wget_final = (int)(strchr(next_wget_final, '-') - next_wget_final);
	strncat(final_epoch_wget_plus40_5days, final_epoch_wget_plus40_5days_temp, 4);
	strncat(final_epoch_wget_plus40_5days, final_epoch_wget_plus40_5days_temp+5, 2);
	strncat(final_epoch_wget_plus40_5days, final_epoch_wget_plus40_5days_temp+8, 2);
	strcpy(str_wget, "wget --no-check-certificate --post-data ");
	strcat(str_wget, "\"activity=retrieve&res=daily&spacecraft=omni2&start_date=");
	strcat(str_wget, initial_epoch_wget_minus40_5days);
	strcat(str_wget, "&end_date=");
	strcat(str_wget, final_epoch_wget_plus40_5days);
	strcat(str_wget,"&vars=50&scale=Linear&ymin=&ymax=&view=0&charsize=&xstyle=0&ystyle=0&symbol=0&symsize=&linestyle=solid&table=0&imagex=640&imagey=480&color=&back=\"");
	strcat(str_wget, " https://omniweb.gsfc.nasa.gov/cgi/nx1.cgi -O ");
	strcpy(filename_f107_to_calculate_f107_average, "");
	  //newstructure 
/* 	  strcat(filename_f107_to_calculate_f107_average, OPTIONS->dir_input_density_msis); */
/* 	  strcat(filename_f107_to_calculate_f107_average, "/"); */
	  //newstructure 
	strcat(filename_f107_to_calculate_f107_average, "f107_for_81daverage_" );
	strcat(filename_f107_to_calculate_f107_average, initial_epoch_wget_minus40_5days);
	strcat(filename_f107_to_calculate_f107_average, "_to_" );
	strcat(filename_f107_to_calculate_f107_average, final_epoch_wget_plus40_5days);
	strcat(filename_f107_to_calculate_f107_average, ".txt");
	strcat(str_wget, filename_f107_to_calculate_f107_average);
	strcat(str_wget, " >/dev/null 2>&1");
	if (do_not_download_file_swpc_or_wget !=1 ){
	system(str_wget);
	}

	// Quick check that the end date is ok
	file_f107_to_calculate_f107_average = fopen(filename_f107_to_calculate_f107_average, "r");
	getline(&line, &len, file_f107_to_calculate_f107_average);
	sscanf(line, "%s", test_wget_omniweb_error);

	if (strcmp(test_wget_omniweb_error, "<H1>") == 0){// if file starts with <H1> then it means that the start date or end date is invalid for omniweb
	  // // go get the correct range for omniweb (end date only) 
	  strcpy(correct_end_date_omniweb_temp, "");

	  getline(&line, &len, file_f107_to_calculate_f107_average);

	  next_wget_final = &line[0];
	  find_wget_final = (int)(strchr(next_wget_final,'-')-next_wget_final);

	  strncat(correct_end_date_omniweb_temp, next_wget_final+find_wget_final+2,8);
	  // // test if epoch end + 40.5 days is more recent than correct_end_date_omniweb_temp (the start date should be ok since the start date for f10.7 was ok)
	  // // // convert correct_end_date_omniweb_temp in a correct format to compare to epoch end + 40.5 days

	  strcpy(correct_end_date_omniweb, "");
	  strncat(correct_end_date_omniweb, correct_end_date_omniweb_temp, 4);
	  strcat(correct_end_date_omniweb,  "-");
	  strncat(correct_end_date_omniweb, correct_end_date_omniweb_temp+4, 2);
	  strcat(correct_end_date_omniweb,  "-");
	  strncat(correct_end_date_omniweb, correct_end_date_omniweb_temp+6, 2);
	  strcat(correct_end_date_omniweb, "T");
	  strcat(correct_end_date_omniweb, "00:00");
	  str2et_c(correct_end_date_omniweb, &et_correct_end_date_omniweb);
	  et_correct_end_date_omniweb_minus_40_5_days = et_correct_end_date_omniweb-40.5*24*3600;
	  et2utc_c(et_correct_end_date_omniweb_minus_40_5_days, "ISOC" ,0 ,255 , correct_end_date_omniweb_minus_40_5_days);
	  // // // compare et_correct_end_date_omniweb to epoch end + 40.5 days (actually 41.5, cf previous comment (look for "41.5" to find comment))

	  if ( et_correct_end_date_omniweb < et_final_epoch_corr + 41.5*24*3600 ){ // the most recent data available at omniweb is too old to compute the average of F10.7 over 81 days. Therefore, upload until the most recent date in omniweb. F10.7 average will be computed a this shorter period (instead of 81 days)
	    printf("*********************************************************************\n***************************** IMPORTANT *****************************\nThe most recent data available at omniweb is too old to compute the average of F10.7 over 81 days. The most recent data available is on %s.\nTherefore, the 81-day average of F10.7 for times later than %s (%s minus 40.5 days) will be calculated over a shorter period than 81 days: 40.5 days + the number of days from this particular time until %s. For example, for the last time step of the propagation, the value of F10.7 81-day average is calculated over %d days.\n*********************************************************************\n*********************************************************************\n", correct_end_date_omniweb, correct_end_date_omniweb_minus_40_5_days, correct_end_date_omniweb, correct_end_date_omniweb, (int)( 41.5 + (et_correct_end_date_omniweb-et_final_epoch_corr)/3600./24. ));
	  }

	  fclose(file_f107_to_calculate_f107_average);
	  //      remove(filename_f107_to_calculate_f107_average);
	  // // Since the most recent data available at omniweb is too old to compute the average of F10.7 over 81 days, download omniweb data until et_correct_end_date_omniweb
	  strcpy(str_wget, "wget --no-check-certificate --post-data ");
	  strcat(str_wget, "\"activity=retrieve&res=daily&spacecraft=omni2&start_date=");
	  strcat(str_wget, initial_epoch_wget_minus40_5days);
	  strcat(str_wget, "&end_date=");
	  strcat(str_wget, correct_end_date_omniweb_temp);
	  strcat(str_wget,"&vars=50&scale=Linear&ymin=&ymax=&view=0&charsize=&xstyle=0&ystyle=0&symbol=0&symsize=&linestyle=solid&table=0&imagex=640&imagey=480&color=&back=\"");
	  strcat(str_wget, " https://omniweb.gsfc.nasa.gov/cgi/nx1.cgi -O ");
	  strcpy(filename_f107_to_calculate_f107_average, "");
	//newstructure 
/* 	strcat(filename_f107_to_calculate_f107_average, OPTIONS->dir_input_density_msis); */
/* 	strcat(filename_f107_to_calculate_f107_average, "/"); */
	//newstructure 
	  strcat(filename_f107_to_calculate_f107_average, "f107_for_81daverage_" );
	  strcat(filename_f107_to_calculate_f107_average, initial_epoch_wget_minus40_5days);
	  strcat(filename_f107_to_calculate_f107_average, "_to_" );
	  strcat(filename_f107_to_calculate_f107_average, correct_end_date_omniweb_temp);
	  strcat(filename_f107_to_calculate_f107_average, ".txt");
	  strcat(str_wget, filename_f107_to_calculate_f107_average);
	  strcat(str_wget, " >/dev/null 2>&1");

	  if (do_not_download_file_swpc_or_wget !=1 ){
	  system(str_wget);
	  }

	}


	else{
	  strcpy(correct_end_date_omniweb_temp, final_epoch_wget_plus40_5days);
	} 

	// Read this file to calculate F10.7 81-day average. Create the file for the results
	//	calculate_f107_average(filename_f107_to_calculate_f107_average, filename_f107A, initial_epoch_wget, final_epoch_wget,correct_end_date_omniweb_temp, et_initial_epoch,  et_final_epoch);
	//      remove(filename_f107_to_calculate_f107_average);

	// Check the dates for Ap (same as we did for F10.7)
	strcpy(str_wget, "wget --no-check-certificate --post-data ");
	strcat(str_wget, "\"activity=retrieve&res=hour&spacecraft=omni2&start_date=");
	// // remove 57 hours ot inital epoch since ap historical starts 57 hours before current time
	et2utc_c(et_oldest_tle_epoch_corr-57*3600., "ISOC" ,0 ,11 , initial_epoch_wget_minus_57hours_temp);
	strcpy(initial_epoch_wget_minus_57hours, "");
	next_wget_initial = &initial_epoch_wget_minus_57hours_temp[0];
	find_wget_initial = (int)(strchr(next_wget_initial, '-') - next_wget_initial);
	strncat(initial_epoch_wget_minus_57hours, initial_epoch_wget_minus_57hours_temp, 4);
	strncat(initial_epoch_wget_minus_57hours, initial_epoch_wget_minus_57hours_temp+5, 2);
	strncat(initial_epoch_wget_minus_57hours, initial_epoch_wget_minus_57hours_temp+8, 2);
	strcat(str_wget, initial_epoch_wget_minus_57hours);
	strcat(str_wget, "&end_date=");
	strcat(str_wget, final_epoch_wget);
	strcat(str_wget,"&vars=49&scale=Linear&ymin=&ymax=&view=0&charsize=&xstyle=0&ystyle=0&symbol=0&symsize=&linestyle=solid&table=0&imagex=640&imagey=480&color=&back=\"");
	strcat(str_wget, " https://omniweb.gsfc.nasa.gov/cgi/nx1.cgi -O ");
	strcpy(filename_ap, "");
	//newstructure 
/* 	strcat(filename_ap, OPTIONS->dir_input_density_msis); */
/* 	strcat(filename_ap, "/"); */
	//newstructure 
	strcat(filename_ap, "ap_" );
	strcat(filename_ap, initial_epoch_wget_minus_57hours);
	strcat(filename_ap, "_to_" );
	strcat(filename_ap, final_epoch_wget);
	strcat(filename_ap, ".txt");
	strcat(str_wget, filename_ap);
	strcat(str_wget, " >/dev/null 2>&1");

	if (do_not_download_file_swpc_or_wget !=1 ){
	system(str_wget);
	}
	// Quick check that the start and end dates are ok
	ap_wget_omniweb = fopen(filename_ap, "r");
	getline(&line, &len, ap_wget_omniweb);
	sscanf(line, "%s", test_wget_omniweb_error);

	if (strcmp(test_wget_omniweb_error, "<H1>") == 0){// if file starts with <H1> then it means that the start date or end date is invalid for omniweb
	  // // go get the correct range for omniweb (end date only) 
	  strcpy(correct_end_date_omniweb_temp, "");
	  getline(&line, &len, ap_wget_omniweb);
	  next_wget_final = &line[0];
	  find_wget_final = (int)(strchr(next_wget_final,'-')-next_wget_final);
	  strncat(correct_end_date_omniweb_temp, next_wget_final+find_wget_final+2,8);
	  // // test if start or end date of epoch is more recent than correct_end_date_omniweb_temp
	  // // // convert correct_end_date_omniweb_temp in a correct format to compare to et_initial_epoch
	  strcpy(correct_end_date_omniweb, "");
	  strncat(correct_end_date_omniweb, correct_end_date_omniweb_temp, 4);
	  strcat(correct_end_date_omniweb,  "-");
	  strncat(correct_end_date_omniweb, correct_end_date_omniweb_temp+4, 2);
	  strcat(correct_end_date_omniweb,  "-");
	  strncat(correct_end_date_omniweb, correct_end_date_omniweb_temp+6, 2);
	  strcat(correct_end_date_omniweb, "T");
	  strcat(correct_end_date_omniweb, "00:00");
	  str2et_c(correct_end_date_omniweb, &et_correct_end_date_omniweb);
	  // // // compare et_correct_end_date_omniweb to et_initial_epoch
	  if  ( ( et_correct_end_date_omniweb < et_oldest_tle_epoch_corr-57*3600. ) || ( et_correct_end_date_omniweb < et_final_epoch_corr ) ){ // if start or end date of epoch is older than correct_end_date_omniweb then the user can not run msis dynamic with omniweb because he/she uses the inital or final epoch of the constellation older than the most recent ap available in omniweb
	    printf("***! Omniweb data (used for the thermospheric density model) is not available for this epoch. You need to choose epoch times more recent than: %s. Or you can use data from the Space Weather Prediction Center (SWPC), which is more up to date. For this, replace 'omniweb' by 'swpc' in the main input file.\nThe program will stop. !***\n", correct_end_date_omniweb);
	    ierr =  MPI_Finalize();
	    exit(0);
	  }
	} 
	fclose(ap_wget_omniweb);      

	// Linear interpolate F10.7 and Ap
	//      calculate_f107_average(filename_f107_to_calculate_f107_average, "omniweb");
	OPTIONS->use_ap_hist = 0;
	
	lin_interpolate(OPTIONS->f107, OPTIONS->f107A, OPTIONS->Ap, OPTIONS->Ap_hist, OPTIONS->et_interpo, &OPTIONS->use_ap_hist, filename_f107_to_calculate_f107_average, filename_ap, "omniweb", OPTIONS->nb_time_steps * 2, OPTIONS->initial_epoch,et_oldest_tle_epoch_corr, OPTIONS->final_epoch,OPTIONS->dt, 999.9,iDebugLevel, iProc);       // "* 2.0" because of the Runge Kunta orfer 4 method

	getline(&line, &len, fp);
/* 	for (ccc = 0; ccc < OPTIONS->nb_time_steps * 2; ccc++){ */
/* 	  etprint(OPTIONS->et_interpo[ccc], ""); */
/* 	  printf("%f %d\n", OPTIONS->f107[ccc], ccc); */
/* 	} */

      }

      else if (strcmp(OPTIONS->test_omniweb_or_external_file, "swpc") == 0){       // if the user chooses to use SWPC data for Ap and F10.7. THe Ap and F10.7 files are automatically downloaded. If the epoch end or the epoch start of the constellation is in the future, then predictions are downloaded from http://www.swpc.noaa.gov/products/usaf-45-day-ap-and-f107cm-flux-forecast (for Ap and F10.7). If the epoch start of the constellation is before the current date, then observations are downloaded from ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/ (DSD: F10.7, DGD: Ap).
	// This is a different approach from a external file or omniweb, a different function than lin_interpolate is used: lin_interpolate_swpc. It is similar to lin_interpolate except that there are two files for F10.7 and Ap (since there could be a mix of observations and predictions (basically, there's always a mix of observations and predictions except if the epoch end is not in the future, in which case there are only observations))
	    OPTIONS->swpc_final_epoch_more_recent_than_45_days_from_current_day = 0;
	// IMPORTANT: we voluntarily look at the epoch start minus 81 days because we need the data for F10.7 81 days before the epoch start in order to calculate F10.7 81-day average.
	int nb_ensembles_density_temp = 0;
	sscanf(line, "%s %d",OPTIONS->test_omniweb_or_external_file, &nb_ensembles_density_temp); // note: the user does not have to put a number after "swpc". If there is no number, then nb_ensembles_density_temp automatically takes the value 0
	OPTIONS->nb_ensembles_density = nb_ensembles_density_temp;

      if (nProcs > 0){
	OPTIONS->nb_ensemble_density_per_proc = (int)(OPTIONS->nb_ensembles_density / nProcs);
      }


	time_t t_current = time(NULL);
	struct tm tm_current = *gmtime(&t_current);    
	int year_current = tm_current.tm_year + 1900;
	int month_current = tm_current.tm_mon + 1;
	int day_current = tm_current.tm_mday;

	char year_current_str[15];
	sprintf(year_current_str, "%d", year_current);
	char month_current_str[15];
	sprintf(month_current_str, "%d", month_current);
	char day_current_str[15];
	sprintf(day_current_str, "%d", day_current);
    
	char date_current[256];
	strcpy(date_current,year_current_str);
	strcat(date_current, "-");
	strcat(date_current, month_current_str);
	strcat(date_current, "-");
	strcat(date_current, day_current_str);
	strcat(date_current, "T");
	strcat(date_current, "00:00:00.000");
	double et_current_day;	
	str2et_c(date_current, &et_current_day);

	//	str2et_c(OPTIONS->initial_epoch, &et_initial_epoch);
	str2et_c(OPTIONS->final_epoch, &et_final_epoch);
	char initial_epoch_midnight[256], final_epoch_midnight[256];
	strcpy(initial_epoch_midnight, "");  strcpy(final_epoch_midnight, "");
	//	strncat(initial_epoch_midnight, &OPTIONS->initial_epoch[0], 10);
	strncat(initial_epoch_midnight, &OPTIONS->oldest_tle_epoch[0], 10); // tle
	strcat(initial_epoch_midnight, "T00:00:00.000");
	strncat(final_epoch_midnight, &OPTIONS->final_epoch[0], 10);
	strcat(final_epoch_midnight, "T00:00:00.000");
	double et_final_epoch_midnight, et_initial_epoch_midnight;
	str2et_c(final_epoch_midnight, &et_final_epoch_midnight);  str2et_c(initial_epoch_midnight, &et_initial_epoch_midnight);

	double et_initial_epoch_minus_81_days;
	char initial_epoch_minus_81_days[256];
	et_initial_epoch_minus_81_days = OPTIONS->et_oldest_tle_epoch - 81. * 24 * 3600;
	et2utc_c(et_initial_epoch_minus_81_days, "ISOC", 0, 255, initial_epoch_minus_81_days);

	if ( et_initial_epoch_minus_81_days >= et_current_day ){
	  OPTIONS->swpc_need_observations = 0; // we wont use observations because the propagation starts in the future
	}
	else{
	  OPTIONS->swpc_need_observations = 1; // we'll use observations because epoch start is older than today
	}
	// Figure out if predictions are used: predictions are needed if the epoch end or the epoch start are in the future compared to the current day (assumed at midnight)
	if ( et_final_epoch_midnight > ( et_current_day - 24*3600.)){ // -24*3600 by security (like this we are sure that if we don't need predictions it means that the final epoch is at least one day before the current day)
	  OPTIONS->swpc_need_predictions = 1; // we will use predictions because the propagations ends in the future
	}
	else{
	  OPTIONS->swpc_need_predictions = 0; // we won't use predictions because the propagation is in the past
	}

	// Download the observations and predictions files
	// // Observations. 
	if ( OPTIONS->swpc_need_observations == 1){

	  // // Arranged by quarters of year if propagation corresponds to the current year. Otherwise arranged by year
	  // // // Figures out if epoch year (start or end) is same as current year. !!! ASSUMPTION: in section #TIME, unless the user chose "now" at the first line, the epoch start and end times have to start with 4 digits, representing the year (so for instance not just 16-08-20 for the 20th of August 2016, but: 2016-08-20)
	  char year_epoch_start_minus_81_days[10];
	  strcpy(year_epoch_start_minus_81_days, "");
	  strncat(year_epoch_start_minus_81_days, &initial_epoch_minus_81_days[0], 4);
	  int year_epoch_start_minus_81_days_int = atoi(year_epoch_start_minus_81_days);
	  char year_epoch_stop[10];
	  strcpy(year_epoch_stop, "");
	  strncat(year_epoch_stop, &OPTIONS->final_epoch[0], 4);
	  int year_epoch_stop_int = atoi(year_epoch_stop);
	  int epoch_year_same_as_current_year;
	  if ( ( year_epoch_start_minus_81_days_int >= year_current ) || ( year_epoch_stop_int >= year_current ) ){
	    epoch_year_same_as_current_year = 1;
	  }
	  else{
	    epoch_year_same_as_current_year = 0;
	  }

	  strcpy(initial_filename_f107_obs, "");
	  //newstructure 
/* 	  strcat(initial_filename_f107_obs, OPTIONS->dir_input_density_msis); */
/* 	  strcat(initial_filename_f107_obs, "/"); */
	  //newstructure 
	  strcpy(final_filename_f107_obs, "");
	  //newstructure 
/* 	  strcat(final_filename_f107_obs, OPTIONS->dir_input_density_msis); */
/* 	  strcat(final_filename_f107_obs, "/"); */
	  //newstructure 

	  if (epoch_year_same_as_current_year == 1){

	    // // // If year of epoch start or end then figure out which quarter of the year epoch start and end correspond to
	    char first_quarter_start[256];
	    strcpy(first_quarter_start, year_current_str);
	    strcat(first_quarter_start, "-01-01T00:00:00.000"); // first quarter starts on march 31st
	    double et_first_quarter_start;
	    str2et_c( first_quarter_start, &et_first_quarter_start );
	    char first_quarter[256];
	    strcpy(first_quarter, year_current_str);
	    strcat(first_quarter, "-04-01T00:00:00.000"); // first quarter ends on march 31st
	    double et_first_quarter;
	    str2et_c( first_quarter, &et_first_quarter );
	    char second_quarter[256];
	    strcpy(second_quarter, year_current_str);
	    strcat(second_quarter, "-07-01T00:00:00.000"); // second quarter ends on june 30th
	    double et_second_quarter;
	    str2et_c( second_quarter, &et_second_quarter );
	    char third_quarter[256];
	    strcpy(third_quarter, year_current_str);
	    strcat(third_quarter, "-10-01T00:00:00.000"); // third quarter ends on september 30th
	    double et_third_quarter;
	    str2et_c( third_quarter, &et_third_quarter );

	    // // // // which quarter for epoch start - 81 days
	    if ( ( et_initial_epoch_minus_81_days < et_first_quarter ) && ( et_initial_epoch_minus_81_days >= et_first_quarter_start) ){
	      strcpy( date_file_obs_initial_swpc, "Q1" );
	      quarter_initial_swpc = 1;
	    }
	    else if ( ( et_initial_epoch_minus_81_days < et_second_quarter ) && ( et_initial_epoch_minus_81_days >= et_first_quarter ) ) {
	      strcpy( date_file_obs_initial_swpc, "Q2" );
	      quarter_initial_swpc = 2;
	    }
	    else if ( ( et_initial_epoch_minus_81_days < et_third_quarter ) && ( et_initial_epoch_minus_81_days >= et_second_quarter ) ) {
	      strcpy( date_file_obs_initial_swpc, "Q3" );
	      quarter_initial_swpc = 3;
	    }

	    else if ( ( et_initial_epoch_minus_81_days >= et_third_quarter ) ){
	      strcpy( date_file_obs_initial_swpc, "Q4" );
	      quarter_initial_swpc = 4;
	    }

	    else{ // year epoch start is older than current year
	      strcpy( date_file_obs_initial_swpc, year_epoch_start_minus_81_days );
	    }

	    // // // // which quarter for epoch end
	    if ( ( et_final_epoch < et_first_quarter ) && ( et_final_epoch >= et_first_quarter_start) ){
	      strcpy( date_file_obs_final_swpc, "Q1" );
	      quarter_final_swpc = 1;
	    }
	    else if ( ( et_final_epoch < et_second_quarter ) && ( et_final_epoch >= et_first_quarter ) ) {
	      strcpy( date_file_obs_final_swpc, "Q2" );
	      quarter_final_swpc = 2;
	    }
	    else if ( ( et_final_epoch < et_third_quarter ) && ( et_final_epoch >= et_second_quarter ) ) {
	      strcpy( date_file_obs_final_swpc, "Q3" );
	      quarter_final_swpc = 3;
	    }

	    else if ( ( et_final_epoch >= et_third_quarter ) ){
	      strcpy( date_file_obs_final_swpc, "Q4" );
	      quarter_final_swpc = 4;
	    }

	    else{ // year epoch end is older than current year
	      strcpy( date_file_obs_final_swpc, year_epoch_stop );
	   
	    }



	    // // // // which quarter for current date
	    if ( ( et_current_day < et_first_quarter ) && ( et_current_day >= et_first_quarter_start) ){
	      quarter_current = 1;
	    }
	    else if ( ( et_current_day < et_second_quarter ) && ( et_current_day >= et_first_quarter ) ) {
	      quarter_current = 2;
	    }
	    else if ( ( et_current_day < et_third_quarter ) && ( et_current_day >= et_second_quarter ) ) {
	      quarter_current = 3;
	    }

	    else if ( ( et_current_day  >= et_third_quarter ) ){
	      quarter_current = 4;
	    }


	    // // // // which quarter for epoch start 
	    if ( ( OPTIONS->et_oldest_tle_epoch < et_first_quarter ) && (  OPTIONS->et_oldest_tle_epoch  >= et_first_quarter_start) ){
	      quarter_initial_epoch = 1;
	    }
	    else if ( (  OPTIONS->et_oldest_tle_epoch  < et_second_quarter ) && (  OPTIONS->et_oldest_tle_epoch  >= et_first_quarter ) ) {
	      quarter_initial_epoch = 2;
	    }
	    else if ( (  OPTIONS->et_oldest_tle_epoch  < et_third_quarter ) && (  OPTIONS->et_oldest_tle_epoch  >= et_second_quarter ) ) {
	      quarter_initial_epoch = 3;
	    }

	    else if ( (  OPTIONS->et_oldest_tle_epoch   >= et_third_quarter ) ){
	      quarter_initial_epoch = 4;
	    }



	  }
	  else{ // both the epoch start and end correspond to a year previous to the current year !!! ASSUMPTION: in section #TIME, unless the user chose "now" at the first line, the epoch start and end times have to start with 4 digits, representing the year (so for instance not just 16-08-20 for the 20th of August 2016, but: 2016-08-20)
	    char year_epoch_start_minus_81_days[10];
	    strcpy(year_epoch_start_minus_81_days, "");
	    strncat(year_epoch_start_minus_81_days, &initial_epoch_minus_81_days[0], 4);
	    strcpy( date_file_obs_initial_swpc, year_epoch_start_minus_81_days);
	    char year_epoch_stop[10];
	    strcpy(year_epoch_stop, "");
	    strncat(year_epoch_stop, &OPTIONS->final_epoch[0], 4);
	    strcpy( date_file_obs_final_swpc, year_epoch_stop);
	  }

	  // // Download the observation files
	  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
	    printf("----- (load_options) Downloading the SWPC observation files of F10.7 and Ap...\n");
	  }



	  int date_file_obs_initial_swpc_int, date_file_obs_final_swpc_int, date_file_obs_initial_swpc_int_ap;

	  char wget_swpc_quarter_to_download[1000], wget_swpc_quarter_to_download_temp[5];
	  int remove_n_file_for_ap = 0;
	  if ( date_file_obs_final_swpc[0] != 'Q' ){// epoch end year != current year
	    date_file_obs_final_swpc_int = atoi(date_file_obs_final_swpc);
	  }
	  if ( date_file_obs_initial_swpc[0] != 'Q' ){// epoch start year != current year
	    date_file_obs_initial_swpc_int = atoi(date_file_obs_initial_swpc);
	  // // we want to avoid downloading the ap file is older than the epoch start date because we don't need it (contrarily to F10.7 for which we need the data 81 days before the current date)
	  char year_epoch_start[10];
	  strcpy(year_epoch_start, "");
	  strncat(year_epoch_start, &OPTIONS->oldest_tle_epoch[0], 4);
	  int year_epoch_start_int = atoi(year_epoch_start);
	  date_file_obs_initial_swpc_int_ap = year_epoch_start_int;
	  if ( date_file_obs_initial_swpc_int_ap != date_file_obs_initial_swpc_int ){
	    remove_n_file_for_ap = remove_n_file_for_ap +1;
	  }
	  }
	  else{
	  if (quarter_initial_swpc < quarter_initial_epoch){
	    quarter_initial_swpc_ap = quarter_initial_swpc + 1;
	    remove_n_file_for_ap = remove_n_file_for_ap + 1;
	  }
	  else{
	    quarter_initial_swpc_ap = quarter_initial_swpc;
	  }
	  }


	  // calculate the number of observations files
	  // // F10.7
	  if ( date_file_obs_final_swpc[0] == 'Q' ){
	      // quarter_final_swpc should be the min of the current quarter and the quarter of the epoch end (because the observations files are not available for future dates than the current date)
	    if ( quarter_current < quarter_final_swpc ){ 
		quarter_final_swpc = quarter_current;
	      }

	    if ( date_file_obs_initial_swpc[0] == 'Q' ){
	      nb_file_obs_swpc = quarter_final_swpc - quarter_initial_swpc + 1;
	    }
	    else{// epoch start year is before current year so add all files of current year quarters until epoch end with files from all years from epoch start year until current year
	      nb_file_obs_swpc = quarter_final_swpc + ( year_current - date_file_obs_initial_swpc_int ); 
	    }
	  }
	  else{
	    nb_file_obs_swpc = date_file_obs_final_swpc_int - date_file_obs_initial_swpc_int + 1;
	  }

	  // // Ap
	  nb_file_obs_swpc_ap = nb_file_obs_swpc - remove_n_file_for_ap;


	  list_filename_obs_swpc_f107  = malloc( nb_file_obs_swpc * sizeof(char *) );
	  list_filename_obs_swpc_ap  = malloc( nb_file_obs_swpc_ap * sizeof(char *) );

	  for (ooo = 0; ooo<nb_file_obs_swpc; ooo++){
	    list_filename_obs_swpc_f107[ooo] = malloc( 1000 * sizeof(char) );
	  }
	  for (ooo = 0; ooo<nb_file_obs_swpc_ap; ooo++){
	    list_filename_obs_swpc_ap[ooo] = malloc( 1000 * sizeof(char) );
	  }


	  // Download observation files
	  

  int ifile_f107 = 0, ifile_ap = 0;
	    if ( date_file_obs_final_swpc[0] == 'Q' ){// epoch end year = current year
	      
	    if ( date_file_obs_initial_swpc[0] == 'Q' ){// epoch start year = current year 
	      for (ooo = quarter_initial_swpc; ooo < quarter_final_swpc +1; ooo++){
		// F10.7
		strcpy(wget_swpc_quarter_to_download, "wget --no-check-certificate ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/");
		strcpy(wget_swpc_quarter_to_download_temp, "");
		sprintf(wget_swpc_quarter_to_download_temp, "%d", ooo );
		strcat(wget_swpc_quarter_to_download, year_current_str);
		strcat(wget_swpc_quarter_to_download, "Q");
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DSD.txt -O ");
		strcat(wget_swpc_quarter_to_download, initial_filename_f107_obs);
		strcat(wget_swpc_quarter_to_download, year_current_str);
		strcat(wget_swpc_quarter_to_download, "Q");
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DSD");
		strcat(wget_swpc_quarter_to_download, ".txt >/dev/null 2>&1");
		//printf("WGET --no-check-certificate <%s>", wget_swpc_quarter_to_download);
		if (do_not_download_file_swpc_or_wget !=1 ){
		system(wget_swpc_quarter_to_download);
		}
		strcpy(list_filename_obs_swpc_f107[ifile_f107], initial_filename_f107_obs);
		strcat(list_filename_obs_swpc_f107[ifile_f107], year_current_str);
		strcat(list_filename_obs_swpc_f107[ifile_f107], "Q");
		strcat(list_filename_obs_swpc_f107[ifile_f107], wget_swpc_quarter_to_download_temp);
		strcat(list_filename_obs_swpc_f107[ifile_f107], "_DSD.txt");
		ifile_f107 = ifile_f107 + 1;
	      }

		// Ap		
	      for (ooo = quarter_initial_swpc_ap; ooo < quarter_final_swpc +1; ooo++){
		strcpy(wget_swpc_quarter_to_download, "wget --no-check-certificate ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/");
		strcpy(wget_swpc_quarter_to_download_temp, "");
		sprintf(wget_swpc_quarter_to_download_temp, "%d", ooo );
		strcat(wget_swpc_quarter_to_download, year_current_str);
		strcat(wget_swpc_quarter_to_download, "Q");
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DGD.txt -O ");
		strcat(wget_swpc_quarter_to_download, initial_filename_f107_obs);
		strcat(wget_swpc_quarter_to_download, year_current_str);
		strcat(wget_swpc_quarter_to_download, "Q");
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DGD");
		strcat(wget_swpc_quarter_to_download, ".txt >/dev/null 2>&1");		
		//printf("WGET --no-check-certificate <%s>", wget_swpc_quarter_to_download);
		if (do_not_download_file_swpc_or_wget !=1 ){
		system(wget_swpc_quarter_to_download);
		  }

		strcpy(list_filename_obs_swpc_ap[ifile_ap], initial_filename_f107_obs);
		strcat(list_filename_obs_swpc_ap[ifile_ap], year_current_str);
		strcat(list_filename_obs_swpc_ap[ifile_ap], "Q");
		strcat(list_filename_obs_swpc_ap[ifile_ap], wget_swpc_quarter_to_download_temp);
		strcat(list_filename_obs_swpc_ap[ifile_ap], "_DGD.txt");
		ifile_ap = ifile_ap + 1;

	      } 
	    }
	    else{ // if start epoch is from a previous year and that end eopch is the current year then download all files from start epoch year to current year and then all files of the current year up to the epoch end quarter
	      //  download all years from initial epoch year to curren year
	      for (ooo = date_file_obs_initial_swpc_int; ooo < year_current; ooo++){
		// F10.7
		strcpy(wget_swpc_quarter_to_download, "wget --no-check-certificate ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/");
		strcpy(wget_swpc_quarter_to_download_temp, "");
		sprintf(wget_swpc_quarter_to_download_temp, "%d", ooo );
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DSD.txt -O ");
		strcat(wget_swpc_quarter_to_download, initial_filename_f107_obs);
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DSD.txt >/dev/null 2>&1");
		//printf("WGET --no-check-certificate <%s>", wget_swpc_quarter_to_download);
		if (do_not_download_file_swpc_or_wget !=1 ){
		system(wget_swpc_quarter_to_download);
		}

		strcpy(list_filename_obs_swpc_f107[ifile_f107], initial_filename_f107_obs);
		strcat(list_filename_obs_swpc_f107[ifile_f107], wget_swpc_quarter_to_download_temp);
		strcat(list_filename_obs_swpc_f107[ifile_f107], "_DSD.txt");
		ifile_f107 = ifile_f107 + 1;

	      }

		// Ap	  
	      for (ooo = date_file_obs_initial_swpc_int_ap; ooo < year_current; ooo++){
		strcpy(wget_swpc_quarter_to_download, "wget --no-check-certificate ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/");
		strcpy(wget_swpc_quarter_to_download_temp, "");
		sprintf(wget_swpc_quarter_to_download_temp, "%d", ooo );
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DGD.txt -O ");
		strcat(wget_swpc_quarter_to_download, initial_filename_f107_obs);
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DGD.txt >/dev/null 2>&1");
		//printf("WGET --no-check-certificate <%s>", wget_swpc_quarter_to_download);
	      
		if (do_not_download_file_swpc_or_wget !=1 ){
		system(wget_swpc_quarter_to_download);
		}

		strcpy(list_filename_obs_swpc_ap[ifile_ap], initial_filename_f107_obs);
		strcat(list_filename_obs_swpc_ap[ifile_ap], wget_swpc_quarter_to_download_temp);
		strcat(list_filename_obs_swpc_ap[ifile_ap], "_DGD.txt");
		ifile_ap = ifile_ap + 1;


		//	      printf("<%s>\n", wget_swpc_quarter_to_download);

	      }

	      // download all quarter files of current year

	      for (ooo = 1; ooo < quarter_final_swpc +1; ooo++){
		// F10.7

		strcpy(wget_swpc_quarter_to_download, "wget --no-check-certificate ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/");
		strcpy(wget_swpc_quarter_to_download_temp, "");
		sprintf(wget_swpc_quarter_to_download_temp, "%d", ooo );
		strcat(wget_swpc_quarter_to_download, year_current_str);
		strcat(wget_swpc_quarter_to_download, "Q");
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DSD.txt -O ");
		strcat(wget_swpc_quarter_to_download, initial_filename_f107_obs);
		strcat(wget_swpc_quarter_to_download, year_current_str);
		strcat(wget_swpc_quarter_to_download, "Q");
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DSD");
		strcat(wget_swpc_quarter_to_download, ".txt >/dev/null 2>&1");
		//printf("WGET --no-check-certificate <%s>", wget_swpc_quarter_to_download);
		if (do_not_download_file_swpc_or_wget !=1 ){
		system(wget_swpc_quarter_to_download);
		}
		strcpy(list_filename_obs_swpc_f107[ifile_f107], initial_filename_f107_obs);
		strcat(list_filename_obs_swpc_f107[ifile_f107], year_current_str);
		strcat(list_filename_obs_swpc_f107[ifile_f107], "Q");
		strcat(list_filename_obs_swpc_f107[ifile_f107], wget_swpc_quarter_to_download_temp);
		strcat(list_filename_obs_swpc_f107[ifile_f107], "_DSD.txt");
		ifile_f107 = ifile_f107 + 1;

		// Ap
		strcpy(wget_swpc_quarter_to_download, "wget --no-check-certificate ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/");
		strcpy(wget_swpc_quarter_to_download_temp, "");
		sprintf(wget_swpc_quarter_to_download_temp, "%d", ooo );
		strcat(wget_swpc_quarter_to_download, year_current_str);
		strcat(wget_swpc_quarter_to_download, "Q");
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DGD.txt -O ");
		strcat(wget_swpc_quarter_to_download, initial_filename_f107_obs);
		strcat(wget_swpc_quarter_to_download, year_current_str);
		strcat(wget_swpc_quarter_to_download, "Q");
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DGD");
		strcat(wget_swpc_quarter_to_download, ".txt >/dev/null 2>&1");
		//printf("WGET --no-check-certificate <%s>", wget_swpc_quarter_to_download);
		if (do_not_download_file_swpc_or_wget !=1 ){
		system(wget_swpc_quarter_to_download);
		}

		strcpy(list_filename_obs_swpc_ap[ifile_ap], initial_filename_f107_obs);
		strcat(list_filename_obs_swpc_ap[ifile_ap], year_current_str);
		strcat(list_filename_obs_swpc_ap[ifile_ap], "Q");
		strcat(list_filename_obs_swpc_ap[ifile_ap], wget_swpc_quarter_to_download_temp);
		strcat(list_filename_obs_swpc_ap[ifile_ap], "_DGD.txt");
		ifile_ap = ifile_ap + 1;
		//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		//   printf("<%s>\n", wget_swpc_quarter_to_download);
	      }

	    }
	  }
	  else{ // end epoch year is older than current year (so start epoch year too)
	    for (ooo = date_file_obs_initial_swpc_int; ooo < date_file_obs_final_swpc_int + 1; ooo++){
	      // F10.7
	      strcpy(wget_swpc_quarter_to_download, "wget --no-check-certificate ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/");
	      strcpy(wget_swpc_quarter_to_download_temp, "");
	      sprintf(wget_swpc_quarter_to_download_temp, "%d", ooo );
	      strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
	      strcat(wget_swpc_quarter_to_download, "_DSD.txt -O ");
	      strcat(wget_swpc_quarter_to_download, initial_filename_f107_obs);
	      strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
	      strcat(wget_swpc_quarter_to_download, "_DSD.txt >/dev/null 2>&1");
		//printf("WGET --no-check-certificate <%s>", wget_swpc_quarter_to_download);
	      if (do_not_download_file_swpc_or_wget !=1 ){
	     system(wget_swpc_quarter_to_download);
	      }

		strcpy(list_filename_obs_swpc_f107[ifile_f107], initial_filename_f107_obs);
		strcat(list_filename_obs_swpc_f107[ifile_f107], wget_swpc_quarter_to_download_temp);
		strcat(list_filename_obs_swpc_f107[ifile_f107], "_DSD.txt");
		ifile_f107 = ifile_f107 + 1;

	    }
	      // Ap
	    for (ooo = date_file_obs_initial_swpc_int_ap; ooo < date_file_obs_final_swpc_int + 1; ooo++){
	      strcpy(wget_swpc_quarter_to_download, "wget --no-check-certificate ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/");
	      strcpy(wget_swpc_quarter_to_download_temp, "");
	      sprintf(wget_swpc_quarter_to_download_temp, "%d", ooo );
	      strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
	      strcat(wget_swpc_quarter_to_download, "_DGD.txt -O ");
	      strcat(wget_swpc_quarter_to_download, initial_filename_f107_obs);
	      strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
	      strcat(wget_swpc_quarter_to_download, "_DGD.txt >/dev/null 2>&1");
		//printf("WGET --no-check-certificate    <%s>", wget_swpc_quarter_to_download);
	      if (do_not_download_file_swpc_or_wget !=1 ){
	     system(wget_swpc_quarter_to_download); 
	      }

		strcpy(list_filename_obs_swpc_ap[ifile_ap], initial_filename_f107_obs);
		strcat(list_filename_obs_swpc_ap[ifile_ap], wget_swpc_quarter_to_download_temp);
		strcat(list_filename_obs_swpc_ap[ifile_ap], "_DGD.txt");
		ifile_ap = ifile_ap + 1;

	      // printf("<%s>\n", wget_swpc_quarter_to_download);

	    }

	  }

	  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
	    printf("----- (load_options) Done downloading the SWPC observation files of F10.7 and Ap.\n");
	  }

	} // end of OPTIONS->swpc_need_observations == 1

	/* for (ooo = 0; ooo < nb_file_obs_swpc; ooo++){ */
	/*   printf("<%s>\n", list_filename_obs_swpc_f107[ooo]); */
	/* } */
	/* /\* //	printf("\n"); *\/ */
	/* for (ooo = 0; ooo < nb_file_obs_swpc_ap; ooo++){ */
	/*   //	  printf("<%s>\n", list_filename_obs_swpc_ap[ooo]); */
	/* } */

	// // Predictions
	if ( OPTIONS->swpc_need_predictions == 1){
	  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
	    printf("----- (load_options) Downloading the SWPC prediction files of F10.7 and Ap.\n");
	  }

	  // // File to download: http://services.swpc.noaa.gov/text/45-day-ap-forecast.txt (always the same name)
	  double et_current_day_plus_45_day;
	  et_current_day_plus_45_day = et_current_day + 45. * 24 * 3600;  
	  char current_day_plus_45_day[256];
	  et2utc_c(et_current_day_plus_45_day, "ISOC", 0, 255, current_day_plus_45_day);
	  if ( (  OPTIONS->et_oldest_tle_epoch  >=  et_current_day_plus_45_day + 24*3600 ) || (  et_final_epoch >=  et_current_day_plus_45_day + 24*3600  ) ){ // actually plus 46  because we have the prediction on day + 45
	    OPTIONS->swpc_final_epoch_more_recent_than_45_days_from_current_day = 1;
	    printf("***! You chose to compute the drag with predictions of the solar activity from SWPC. However, these predictions are limited to 45 days in the future (%s). Therefore, for times beyond %s, Ap and F10.7 are set to constant values, equal to the average of the values predicted by SWPC over 45 days.\n", current_day_plus_45_day, current_day_plus_45_day,current_day_plus_45_day);// MPI_Finalize(); exit(0);
	  }


	    strcpy( wget_filename_f107_ap_pred, "wget --no-check-certificate http://services.swpc.noaa.gov/text/45-day-ap-forecast.txt -O ");
	    strcpy(filename_f107_ap_pred, "");
	    //	    printf("<%s>\n",  OPTIONS->dir_input_density_msis);
/* 	  //newstructure  */
/* /\* 	    strcat(filename_f107_ap_pred, OPTIONS->dir_input_density_msis); *\/ */
/* /\* 	    strcat(filename_f107_ap_pred, "/"); *\/ */
/* 	  //newstructure  */
 	    strcat(filename_f107_ap_pred, "swpc_predictions_f107_ap.txt" );
	    //	    printf("<%s>\n",  filename_f107_ap_pred);
	    strcat(wget_filename_f107_ap_pred, filename_f107_ap_pred);
	    strcat(wget_filename_f107_ap_pred, " >/dev/null 2>&1");
	    //	    printf("PREF <%s>\n", wget_filename_f107_ap_pred);
	    if (do_not_download_file_swpc_or_wget !=1 ){
	    system(wget_filename_f107_ap_pred);
	    }

	  
	  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
	    printf("----- (load_options) Done downloading the SWPC prediction files of F10.7 and Ap.\n");
	  }

	

	 
	if (( iDebugLevel >= 4 ) ){
	    printf("----- (load_options) Allocating memory for sigma_f107, sigma_f107a, sigma_ap, and et_sigma_f107_ap (iProc %d)\n", iProc);
	  }
 MPI_Barrier( MPI_COMM_WORLD ) ; 
	  if (( iDebugLevel >= 4 ) ){
	    printf("----- (load_options) Went through MPI_Barrier (iProc %d)\n", iProc);
	  }

	OPTIONS->sigma_f107 = malloc( ( OPTIONS->nb_time_steps + (int)(( 2 * 24 * 3600. ) / OPTIONS->dt) ) * 2 * sizeof(double) ); // sigma_f107 does ahave that many elements. This is a maximum.  + (int)(( 2 * 24 * 3600. ) / OPTIONS->dt) because sigma_f107 will be filled with the entire day that includes intial_epoch (24 * 3600 would be enough but we double it to be safe)
      if (  OPTIONS->sigma_f107 == NULL ){
	printf("***! Could not allow memory space for  OPTIONS->sigma_f107 \n. The program will stop. !***\n");
	ierr =  MPI_Finalize();
	exit(0);
      }

      OPTIONS->sigma_ap = malloc( ( OPTIONS->nb_time_steps + (int)(( 2 * 24 * 3600. ) / OPTIONS->dt) ) * 2 * sizeof(double) );
      if (  OPTIONS->sigma_ap == NULL ){
	printf("***! Could not allow memory space for  OPTIONS->sigma_ap \n. The program will stop. !***\n");
	ierr =  MPI_Finalize();
	exit(0);
      }

      OPTIONS->et_sigma_f107_ap = malloc( ( OPTIONS->nb_time_steps + (int)(( 2 * 24 * 3600. ) / OPTIONS->dt) ) * 2 * sizeof(double) );

      if (  OPTIONS->et_sigma_f107_ap == NULL ){
	printf("***! Could not allow memory space for  OPTIONS->et_sigma_f107_ap \n. The program will stop. !***\n");
	ierr =  MPI_Finalize();
	exit(0);
      }
	  if ( ( iDebugLevel >= 4 ) ){
	    printf("----- (load_options) Done allocating memory for sigma_f107, sigma_f107a, sigma_ap, and et_sigma_f107_ap (iProc %d)\n", iProc);
	  }
} // end of if swpc_need_predictions == 1

	lin_interpolate_swpc(OPTIONS->f107, OPTIONS->f107A, OPTIONS->Ap, OPTIONS->et_interpo, OPTIONS->sigma_f107, OPTIONS->sigma_ap, OPTIONS->et_sigma_f107_ap, &OPTIONS->use_ap_hist, list_filename_obs_swpc_f107, list_filename_obs_swpc_ap, filename_f107_ap_pred,  OPTIONS->nb_time_steps * 2, OPTIONS->initial_epoch, OPTIONS->et_oldest_tle_epoch, OPTIONS->final_epoch, OPTIONS->dt, 999.9,iDebugLevel, nb_file_obs_swpc, nb_file_obs_swpc_ap, OPTIONS->swpc_need_observations, OPTIONS->swpc_need_predictions, OPTIONS->nb_ensembles_density, OPTIONS->dir_input_density_msis, &OPTIONS->swpc_et_first_prediction, iProc, OPTIONS->swpc_final_epoch_more_recent_than_45_days_from_current_day); 
	  //	  printf("%f %f %d\n", OPTIONS->f107[2], OPTIONS->Ap[2], iProc);
	  if ( OPTIONS->swpc_need_predictions == 1){
	  generate_ensemble_f107_ap( OPTIONS, iDebugLevel, iProc );
	  }
      // "* 2.0" because of the Runge Kunta orfer 4 method

	free(list_filename_obs_swpc_f107);
	free(list_filename_obs_swpc_ap);
	//`	MPI_Finalize();exit(0);
      }


  else if (strcmp(OPTIONS->test_omniweb_or_external_file, "swpc_mod") == 0){       // if the user chooses to use SWPC data for Ap and F10.7. THe Ap and F10.7 files are automatically downloaded. If the epoch end or the epoch start of the constellation is in the future, then predictions are downloaded from http://www.swpc.noaa.gov/products/usaf-45-day-ap-and-f107cm-flux-forecast (for Ap and F10.7). If the epoch start of the constellation is before the current date, then observations are downloaded from ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/ (DSD: F10.7, DGD: Ap).
	// This is a different approach from a external file or omniweb, a different function than lin_interpolate is used: lin_interpolate_swpc (actually lin_interpolate_swpc_mod here but it works exactly the same as lin_interpolate_swpc except if running ensembles on the density or if modication of the nominal value of f107 and Ap). It is similar to lin_interpolate except that there are two files for F10.7 and Ap (since there could be a mix of observations and predictions (basically, there's always a mix of observations and predictions except if the epoch end is not in the future, in which case there are only observations))
    // the difference with swpc_mod is that here SpOCK reads in the text file (name right after 'swpc_mod' in main input file) that tells it to apply a variation on the value from SWPC. For example, if SWPC gives F10.7 = 100 in 2 days and that in this file it is written '20' (respectively '-10') at the line +2 days then the value of 10.7 used in the propagation is 120 (respectively '90')



	// IMPORTANT: we voluntarily look at the epoch start minus 81 days because we need the data for F10.7 81 days before the epoch start in order to calculate F10.7 81-day average.
    char filename_f107_ap_mod[1000]; // name of the file that will ber ead to appply a delta on f107 and ap
    //newstructure
/*     strcpy(filename_f107_ap_mod, OPTIONS->dir_input_density_msis); */
/*     strcat(filename_f107_ap_mod, "/"); */

    if ( (strcmp(OPTIONS->type_orbit_initialisation, "tle" ) == 0 ) || (strcmp(OPTIONS->type_orbit_initialisation, "tle_sgp4" ) == 0 )){
      print_error_any_iproc(iProc, "You chose the initialize the orbits with a TLE and you chose the option swpc_mod for the modeling of the density. However, there is a bug in SpOCK that needs to be fixed and prevents running it in this configuration (TLE + spwc_mod)"); // all the other lin_interpolate functions don't have this bug (lin_interpolate, lin_interpolate_swpc, lin_interpolate_ap_hist, lin_interpolate_attitude) but I was too lazy to fix lin_interpolate_swpc_mod too. The bug is basically die to the fact that SpOCK use to propagate the sc from the TLE epoch until the initial epoch (set in section #TIME) with NO drag. But then I changed it os that drag is also computed for this time lapse. I had to fix a lot of things and I didn't so it for lin_interpolate_swpc_mod
    }

    strcpy(filename_f107_ap_mod, "");
    //newstructure
    char filename_f107_ap_mod_temp[1000];
    char swpc_pred_test[1000];
	    strcpy(filename_f107_ap_pred, "");
	int nb_spaces ;
	char *line_copy ;
	line_copy = malloc(sizeof(char) * strlen(line));
	strcpy(line_copy, line);
	for (nb_spaces=0; line_copy[nb_spaces]; line_copy[nb_spaces]=='.' ? nb_spaces++ : *line_copy++);

	if (nb_spaces > 1){ // the user input their own spwc f107 ap prediction file
	      sscanf(line, "%s %s %s",OPTIONS->test_omniweb_or_external_file, filename_f107_ap_mod_temp, filename_f107_ap_pred);
	}
	else{ // the user did not input a swpc ap f107 prediciton file so it's the deault one that is used: swpc_predictions_f107_ap.txt
	      sscanf(line, "%s %s",OPTIONS->test_omniweb_or_external_file, filename_f107_ap_mod_temp);
	    strcat(filename_f107_ap_pred, "swpc_predictions_f107_ap.txt" );	  
	}

	strcat(filename_f107_ap_mod, filename_f107_ap_mod_temp);
	//	printf("<%s><%s>\n", filename_f107_ap_mod, filename_f107_ap_pred);


	time_t t_current = time(NULL);
	struct tm tm_current = *gmtime(&t_current);    
	tm_current.tm_mday = tm_current.tm_mday ; // !!!!!!!!! change current time. should be tm_current.tm_mday - 0 so if it is not (different number than '0') then change it back to tm_current.tm_mday - 0 unless you really know what you're doing!
	t_current = mktime(&tm_current);
	int year_current = tm_current.tm_year + 1900;  
	//				year_current = 2017;// !!!!!!!!!!!!!!!!!! remove collision paper 2 and paper 3
	//	year_current = 2016;// !!!!!!!!!!!!!!!!!! remove collision paper 1
	int month_current = tm_current.tm_mon + 1;
	int day_current = tm_current.tm_mday;



	/* time_t t_current = time(NULL); */
	/* struct tm tm_current = *gmtime(&t_current); */
	/* int year_current = tm_current.tm_year + 1900; */
	/* int month_current = tm_current.tm_mon + 1; */
	/* int day_current = tm_current.tm_mday; */

	char year_current_str[15];
	sprintf(year_current_str, "%d", year_current);
	char month_current_str[15];
	sprintf(month_current_str, "%d", month_current);
	char day_current_str[15];
	sprintf(day_current_str, "%d", day_current);
    
	char date_current[256];
	strcpy(date_current,year_current_str);
	strcat(date_current, "-");
	strcat(date_current, month_current_str);
	strcat(date_current, "-");
	strcat(date_current, day_current_str);
	strcat(date_current, "T");
	strcat(date_current, "00:00:00.000");
	//	strcpy(date_current,"2016-11-26T00:00:00");// !!!!!!!!!!! remove collision paper 1
	//		strcpy(date_current,"2017-12-01T00:00:00");// !!!!!!!!!!! remove collision paper 2
	//	strcpy(date_current,"2017-11-01T00:00:00");// !!!!!!!!!!! remove collision paper 3 conjunction fm07 november
	//		strcpy(date_current,"2017-07-22T00:00:00");// !!!!!!!!!!! remove collision paper 3 conjunction fm07 july 22
	//			strcpy(date_current,"2017-07-23T00:00:00");// !!!!!!!!!!! remove collision paper 3 conjunction fm07 july 23

	double et_current_day;
	str2et_c(date_current, &et_current_day);
	//	str2et_c(OPTIONS->initial_epoch, &et_initial_epoch);
	str2et_c(OPTIONS->final_epoch, &et_final_epoch);
	char initial_epoch_midnight[256], final_epoch_midnight[256];
	strcpy(initial_epoch_midnight, "");  strcpy(final_epoch_midnight, "");
	//	strncat(initial_epoch_midnight, &OPTIONS->initial_epoch[0], 10);
	strncat(initial_epoch_midnight, &OPTIONS->oldest_tle_epoch[0], 10); // tle
	strcat(initial_epoch_midnight, "T00:00:00.000");
	strncat(final_epoch_midnight, &OPTIONS->final_epoch[0], 10);
	strcat(final_epoch_midnight, "T00:00:00.000");
	double et_final_epoch_midnight, et_initial_epoch_midnight;
	str2et_c(final_epoch_midnight, &et_final_epoch_midnight);  str2et_c(initial_epoch_midnight, &et_initial_epoch_midnight);

	double et_initial_epoch_minus_81_days;
	char initial_epoch_minus_81_days[256];
	et_initial_epoch_minus_81_days =  OPTIONS->et_oldest_tle_epoch  - 81. * 24 * 3600;
	et2utc_c(et_initial_epoch_minus_81_days, "ISOC", 0, 255, initial_epoch_minus_81_days);

	if ( et_initial_epoch_minus_81_days >= et_current_day ){
	  OPTIONS->swpc_need_observations = 0; // we wont use observations because the propagation starts in the future
	}
	else{
	  OPTIONS->swpc_need_observations = 1; // we'll use observations because epoch start is older than today
	}
	// Figure out if predictions are used: predictions are needed if the epoch end or the epoch start are in the future compared to the current day (assumed at midnight)
	if ( et_final_epoch_midnight > ( et_current_day - 24*3600.)){ // -24*3600 by security (like this we are sure that if we don't need predictions it means that the final epoch is at least one day before the current day)
	  OPTIONS->swpc_need_predictions = 1; // we will use predictions because the propagations ends in the future
	}
	else{
	  OPTIONS->swpc_need_predictions = 0; // we won't use predictions because the propagation is in the past
	}

	// Download the observations and predictions files
	// // Observations.
	if ( OPTIONS->swpc_need_observations == 1){
	  // // Arranged by quarters of year if propagation corresponds to the current year. Otherwise arranged by year
	  // // // Figures out if epoch year (start or end) is same as current year. !!! ASSUMPTION: in section #TIME, unless the user chose "now" at the first line, the epoch start and end times have to start with 4 digits, representing the year (so for instance not just 16-08-20 for the 20th of August 2016, but: 2016-08-20)
	  char year_epoch_start_minus_81_days[10];
	  strcpy(year_epoch_start_minus_81_days, "");
	  strncat(year_epoch_start_minus_81_days, &initial_epoch_minus_81_days[0], 4);
	  int year_epoch_start_minus_81_days_int = atoi(year_epoch_start_minus_81_days);
	  char year_epoch_stop[10];
	  strcpy(year_epoch_stop, "");
	  strncat(year_epoch_stop, &OPTIONS->final_epoch[0], 4);
	  int year_epoch_stop_int = atoi(year_epoch_stop);
	  int epoch_year_same_as_current_year;
	  if ( ( year_epoch_start_minus_81_days_int >= year_current ) || ( year_epoch_stop_int >= year_current ) ){
	    epoch_year_same_as_current_year = 1;
	  }
	  else{
	    epoch_year_same_as_current_year = 0;
	  }

	  strcpy(initial_filename_f107_obs, "");
	  //newstructure
/* 	  strcat(initial_filename_f107_obs, OPTIONS->dir_input_density_msis); */
/* 	  strcat(initial_filename_f107_obs, "/"); */
	  //newstructure
	  strcpy(final_filename_f107_obs, "");
	  //newstructure
/* 	  strcat(final_filename_f107_obs, OPTIONS->dir_input_density_msis); */
/* 	  strcat(final_filename_f107_obs, "/"); */
	  //newstructure

	  if (epoch_year_same_as_current_year == 1){

	    // // // If year of epoch start or end then figure out which quarter of the year epoch start and end correspond to
	    char first_quarter_start[256];
	    strcpy(first_quarter_start, year_current_str);
	    strcat(first_quarter_start, "-01-01T00:00:00.000"); // first quarter starts on march 31st
	    double et_first_quarter_start;
	    str2et_c( first_quarter_start, &et_first_quarter_start );
	    char first_quarter[256];
	    strcpy(first_quarter, year_current_str);
	    strcat(first_quarter, "-04-01T00:00:00.000"); // first quarter ends on march 31st
	    double et_first_quarter;
	    str2et_c( first_quarter, &et_first_quarter );
	    char second_quarter[256];
	    strcpy(second_quarter, year_current_str);
	    strcat(second_quarter, "-07-01T00:00:00.000"); // second quarter ends on june 30th
	    double et_second_quarter;
	    str2et_c( second_quarter, &et_second_quarter );
	    char third_quarter[256];
	    strcpy(third_quarter, year_current_str);
	    strcat(third_quarter, "-10-01T00:00:00.000"); // third quarter ends on september 30th
	    double et_third_quarter;
	    str2et_c( third_quarter, &et_third_quarter );

	    // // // // which quarter for epoch start - 81 days
	    if ( ( et_initial_epoch_minus_81_days < et_first_quarter ) && ( et_initial_epoch_minus_81_days >= et_first_quarter_start) ){
	      strcpy( date_file_obs_initial_swpc, "Q1" );
	      quarter_initial_swpc = 1;
	    }
	    else if ( ( et_initial_epoch_minus_81_days < et_second_quarter ) && ( et_initial_epoch_minus_81_days >= et_first_quarter ) ) {
	      strcpy( date_file_obs_initial_swpc, "Q2" );
	      quarter_initial_swpc = 2;
	    }
	    else if ( ( et_initial_epoch_minus_81_days < et_third_quarter ) && ( et_initial_epoch_minus_81_days >= et_second_quarter ) ) {
	      strcpy( date_file_obs_initial_swpc, "Q3" );
	      quarter_initial_swpc = 3;
	    }

	    else if ( ( et_initial_epoch_minus_81_days >= et_third_quarter ) ){
	      strcpy( date_file_obs_initial_swpc, "Q4" );
	      quarter_initial_swpc = 4;
	    }

	    else{ // year epoch start is older than current year
	      strcpy( date_file_obs_initial_swpc, year_epoch_start_minus_81_days );
	    }

	    // // // // which quarter for epoch end
	    if ( ( et_final_epoch < et_first_quarter ) && ( et_final_epoch >= et_first_quarter_start) ){
	      strcpy( date_file_obs_final_swpc, "Q1" );
	      quarter_final_swpc = 1;
	    }
	    else if ( ( et_final_epoch < et_second_quarter ) && ( et_final_epoch >= et_first_quarter ) ) {
	      strcpy( date_file_obs_final_swpc, "Q2" );
	      quarter_final_swpc = 2;
	    }
	    else if ( ( et_final_epoch < et_third_quarter ) && ( et_final_epoch >= et_second_quarter ) ) {
	      strcpy( date_file_obs_final_swpc, "Q3" );
	      quarter_final_swpc = 3;
	    }

	    else if ( ( et_final_epoch >= et_third_quarter ) ){
	      strcpy( date_file_obs_final_swpc, "Q4" );
	      quarter_final_swpc = 4;
	    }

	    else{ // year epoch end is older than current year
	      strcpy( date_file_obs_final_swpc, year_epoch_stop );
	   
	    }



	    // // // // which quarter for current date
	    if ( ( et_current_day < et_first_quarter ) && ( et_current_day >= et_first_quarter_start) ){
	      quarter_current = 1;
	    }
	    else if ( ( et_current_day < et_second_quarter ) && ( et_current_day >= et_first_quarter ) ) {
	      quarter_current = 2;
	    }
	    else if ( ( et_current_day < et_third_quarter ) && ( et_current_day >= et_second_quarter ) ) {
	      quarter_current = 3;
	    }

	    else if ( ( et_current_day  >= et_third_quarter ) ){
	      quarter_current = 4;
	    }


	    // // // // which quarter for epoch start
	    if ( (  OPTIONS->et_oldest_tle_epoch  < et_first_quarter ) && (  OPTIONS->et_oldest_tle_epoch  >= et_first_quarter_start) ){
	      quarter_initial_epoch = 1;
	    }
	    else if ( (  OPTIONS->et_oldest_tle_epoch  < et_second_quarter ) && ( OPTIONS->et_oldest_tle_epoch >= et_first_quarter ) ) {
	      quarter_initial_epoch = 2;
	    }
	    else if ( ( OPTIONS->et_oldest_tle_epoch  < et_third_quarter ) && ( OPTIONS->et_oldest_tle_epoch  >= et_second_quarter ) ) {
	      quarter_initial_epoch = 3;
	    }

	    else if ( ( OPTIONS->et_oldest_tle_epoch   >= et_third_quarter ) ){
	      quarter_initial_epoch = 4;
	    }



	  }
	  else{ // both the epoch start and end correspond to a year previous to the current year !!! ASSUMPTION: in section #TIME, unless the user chose "now" at the first line, the epoch start and end times have to start with 4 digits, representing the year (so for instance not just 16-08-20 for the 20th of August 2016, but: 2016-08-20)
	    char year_epoch_start_minus_81_days[10];
	    strcpy(year_epoch_start_minus_81_days, "");
	    strncat(year_epoch_start_minus_81_days, &initial_epoch_minus_81_days[0], 4);
	    strcpy( date_file_obs_initial_swpc, year_epoch_start_minus_81_days);
	    char year_epoch_stop[10];
	    strcpy(year_epoch_stop, "");
	    strncat(year_epoch_stop, &OPTIONS->final_epoch[0], 4);
	    strcpy( date_file_obs_final_swpc, year_epoch_stop);
	  }

	  // // Download the observation files
	  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
	    printf("----- (load_options) Downloading the SWPC observation files of F10.7 and Ap...\n");
	  }



	  int date_file_obs_initial_swpc_int, date_file_obs_final_swpc_int, date_file_obs_initial_swpc_int_ap;

	  char wget_swpc_quarter_to_download[1000], wget_swpc_quarter_to_download_temp[5];
	  int remove_n_file_for_ap = 0;
	  if ( date_file_obs_final_swpc[0] != 'Q' ){// epoch end year != current year
	    date_file_obs_final_swpc_int = atoi(date_file_obs_final_swpc);
	  }
	  if ( date_file_obs_initial_swpc[0] != 'Q' ){// epoch start year != current year
	    date_file_obs_initial_swpc_int = atoi(date_file_obs_initial_swpc);
	  // // we want to avoid downloading the ap file is older than the epoch start date because we don't need it (contrarily to F10.7 for which we need the data 81 days before the current date)
	  char year_epoch_start[10];
	  strcpy(year_epoch_start, "");
	  strncat(year_epoch_start, &OPTIONS->initial_epoch[0], 4);
	  int year_epoch_start_int = atoi(year_epoch_start);
	  date_file_obs_initial_swpc_int_ap = year_epoch_start_int;
	  if ( date_file_obs_initial_swpc_int_ap != date_file_obs_initial_swpc_int ){
	    remove_n_file_for_ap = remove_n_file_for_ap +1;
	  }
	  }
	  else{
	  if (quarter_initial_swpc < quarter_initial_epoch){
	    quarter_initial_swpc_ap = quarter_initial_swpc + 1;
	    remove_n_file_for_ap = remove_n_file_for_ap + 1;
	  }
	  else{
	    quarter_initial_swpc_ap = quarter_initial_swpc;
	  }
	  }


	  // calculate the number of observations files
	  // // F10.7
	  if ( date_file_obs_final_swpc[0] == 'Q' ){
	      // quarter_final_swpc should be the min of the current quarter and the quarter of the epoch end (because the observations files are not available for future dates than the current date)
	    if ( quarter_current < quarter_final_swpc ){
		quarter_final_swpc = quarter_current;
	      }

	    if ( date_file_obs_initial_swpc[0] == 'Q' ){
	      nb_file_obs_swpc = quarter_final_swpc - quarter_initial_swpc + 1;
	    }
	    else{// epoch start year is before current year so add all files of current year quarters until epoch end with files from all years from epoch start year until current year
	      nb_file_obs_swpc = quarter_final_swpc + ( year_current - date_file_obs_initial_swpc_int );
	    }
	  }
	  else{
	    nb_file_obs_swpc = date_file_obs_final_swpc_int - date_file_obs_initial_swpc_int + 1;
	  }

	  // // Ap
	  nb_file_obs_swpc_ap = nb_file_obs_swpc - remove_n_file_for_ap;


	  list_filename_obs_swpc_f107  = malloc( nb_file_obs_swpc * sizeof(char *) );
	  list_filename_obs_swpc_ap  = malloc( nb_file_obs_swpc_ap * sizeof(char *) );

	  for (ooo = 0; ooo<nb_file_obs_swpc; ooo++){
	    list_filename_obs_swpc_f107[ooo] = malloc( 1000 * sizeof(char) );
	  }
	  for (ooo = 0; ooo<nb_file_obs_swpc_ap; ooo++){
	    list_filename_obs_swpc_ap[ooo] = malloc( 1000 * sizeof(char) );
	  }


	  // Download observation files
	 
  int ifile_f107 = 0, ifile_ap = 0;
	    if ( date_file_obs_final_swpc[0] == 'Q' ){// epoch end year = current year
	    if ( date_file_obs_initial_swpc[0] == 'Q' ){// epoch start year = current year
	      for (ooo = quarter_initial_swpc; ooo < quarter_final_swpc +1; ooo++){
		// F10.7
		strcpy(wget_swpc_quarter_to_download, "wget --no-check-certificate ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/");
		strcpy(wget_swpc_quarter_to_download_temp, "");
		sprintf(wget_swpc_quarter_to_download_temp, "%d", ooo );
		strcat(wget_swpc_quarter_to_download, year_current_str);
		strcat(wget_swpc_quarter_to_download, "Q");
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DSD.txt -O ");
		strcat(wget_swpc_quarter_to_download, initial_filename_f107_obs);
		strcat(wget_swpc_quarter_to_download, year_current_str);
		strcat(wget_swpc_quarter_to_download, "Q");
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DSD");
		strcat(wget_swpc_quarter_to_download, ".txt >/dev/null 2>&1");
		//printf("WGET --no-check-certificate <%s>", wget_swpc_quarter_to_download);
		if (do_not_download_file_swpc_or_wget !=1 ){
		system(wget_swpc_quarter_to_download);
		}
		strcpy(list_filename_obs_swpc_f107[ifile_f107], initial_filename_f107_obs);
		strcat(list_filename_obs_swpc_f107[ifile_f107], year_current_str);
		strcat(list_filename_obs_swpc_f107[ifile_f107], "Q");
		strcat(list_filename_obs_swpc_f107[ifile_f107], wget_swpc_quarter_to_download_temp);
		strcat(list_filename_obs_swpc_f107[ifile_f107], "_DSD.txt");
		ifile_f107 = ifile_f107 + 1;
	      }

		// Ap
	      for (ooo = quarter_initial_swpc_ap; ooo < quarter_final_swpc +1; ooo++){
		strcpy(wget_swpc_quarter_to_download, "wget --no-check-certificate ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/");
		strcpy(wget_swpc_quarter_to_download_temp, "");
		sprintf(wget_swpc_quarter_to_download_temp, "%d", ooo );
		strcat(wget_swpc_quarter_to_download, year_current_str);
		strcat(wget_swpc_quarter_to_download, "Q");
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DGD.txt -O ");
		strcat(wget_swpc_quarter_to_download, initial_filename_f107_obs);
		strcat(wget_swpc_quarter_to_download, year_current_str);
		strcat(wget_swpc_quarter_to_download, "Q");
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DGD");
		strcat(wget_swpc_quarter_to_download, ".txt >/dev/null 2>&1");
		//printf("WGET --NO-CHECK-CERTIFICATE  <%s>", wget_swpc_quarter_to_download);
		if (do_not_download_file_swpc_or_wget !=1 ){
		system(wget_swpc_quarter_to_download);
		}

		strcpy(list_filename_obs_swpc_ap[ifile_ap], initial_filename_f107_obs);
		strcat(list_filename_obs_swpc_ap[ifile_ap], year_current_str);
		strcat(list_filename_obs_swpc_ap[ifile_ap], "Q");
		strcat(list_filename_obs_swpc_ap[ifile_ap], wget_swpc_quarter_to_download_temp);
		strcat(list_filename_obs_swpc_ap[ifile_ap], "_DGD.txt");
		ifile_ap = ifile_ap + 1;

	      }
	    }
	    else{ // if start epoch is from a previous year and that end eopch is the current year then download all files from start epoch year to current year and then all files of the current year up to the epoch end quarter
	      //  download all years from initial epoch year to curren year
	      for (ooo = date_file_obs_initial_swpc_int; ooo < year_current; ooo++){
		// F10.7
		strcpy(wget_swpc_quarter_to_download, "wget --no-check-certificate  ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/");
		strcpy(wget_swpc_quarter_to_download_temp, "");
		sprintf(wget_swpc_quarter_to_download_temp, "%d", ooo );
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DSD.txt -O ");
		strcat(wget_swpc_quarter_to_download, initial_filename_f107_obs);
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DSD.txt >/dev/null 2>&1");
		//printf("WGET --NO-CHECK-CERTIFICATE  <%s>", wget_swpc_quarter_to_download);
		if (do_not_download_file_swpc_or_wget !=1 ){
		system(wget_swpc_quarter_to_download);
		}
		strcpy(list_filename_obs_swpc_f107[ifile_f107], initial_filename_f107_obs);
		strcat(list_filename_obs_swpc_f107[ifile_f107], wget_swpc_quarter_to_download_temp);
		strcat(list_filename_obs_swpc_f107[ifile_f107], "_DSD.txt");
		ifile_f107 = ifile_f107 + 1;

	      }

		// Ap
	      for (ooo = date_file_obs_initial_swpc_int_ap; ooo < year_current; ooo++){
		strcpy(wget_swpc_quarter_to_download, "wget --no-check-certificate  ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/");
		strcpy(wget_swpc_quarter_to_download_temp, "");
		sprintf(wget_swpc_quarter_to_download_temp, "%d", ooo );
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DGD.txt -O ");
		strcat(wget_swpc_quarter_to_download, initial_filename_f107_obs);
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DGD.txt >/dev/null 2>&1");
		//printf("WGET --NO-CHECK-CERTIFICATE  <%s>", wget_swpc_quarter_to_download);
		if (do_not_download_file_swpc_or_wget !=1 ){
		system(wget_swpc_quarter_to_download);
		}

		strcpy(list_filename_obs_swpc_ap[ifile_ap], initial_filename_f107_obs);
		strcat(list_filename_obs_swpc_ap[ifile_ap], wget_swpc_quarter_to_download_temp);
		strcat(list_filename_obs_swpc_ap[ifile_ap], "_DGD.txt");
		ifile_ap = ifile_ap + 1;


		//	      printf("<%s>\n", wget_swpc_quarter_to_download);

	      }

	      // download all quarter files of current year
	      for (ooo = 1; ooo < quarter_final_swpc +1; ooo++){
		// F10.7
		strcpy(wget_swpc_quarter_to_download, "wget --no-check-certificate  ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/");
		strcpy(wget_swpc_quarter_to_download_temp, "");
		sprintf(wget_swpc_quarter_to_download_temp, "%d", ooo );
		strcat(wget_swpc_quarter_to_download, year_current_str);
		strcat(wget_swpc_quarter_to_download, "Q");
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DSD.txt -O ");
		strcat(wget_swpc_quarter_to_download, initial_filename_f107_obs);
		strcat(wget_swpc_quarter_to_download, year_current_str);
		strcat(wget_swpc_quarter_to_download, "Q");
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DSD");
		strcat(wget_swpc_quarter_to_download, ".txt >/dev/null 2>&1");
		//printf("WGET --NO-CHECK-CERTIFICATE  <%s>", wget_swpc_quarter_to_download);
		if (do_not_download_file_swpc_or_wget !=1 ){
		system(wget_swpc_quarter_to_download);
		}
		strcpy(list_filename_obs_swpc_f107[ifile_f107], initial_filename_f107_obs);
		strcat(list_filename_obs_swpc_f107[ifile_f107], year_current_str);
		strcat(list_filename_obs_swpc_f107[ifile_f107], "Q");
		strcat(list_filename_obs_swpc_f107[ifile_f107], wget_swpc_quarter_to_download_temp);
		strcat(list_filename_obs_swpc_f107[ifile_f107], "_DSD.txt");
		ifile_f107 = ifile_f107 + 1;

		// Ap
		strcpy(wget_swpc_quarter_to_download, "wget --no-check-certificate  ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/");
		strcpy(wget_swpc_quarter_to_download_temp, "");
		sprintf(wget_swpc_quarter_to_download_temp, "%d", ooo );
		strcat(wget_swpc_quarter_to_download, year_current_str);
		strcat(wget_swpc_quarter_to_download, "Q");
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DGD.txt -O ");
		strcat(wget_swpc_quarter_to_download, initial_filename_f107_obs);
		strcat(wget_swpc_quarter_to_download, year_current_str);
		strcat(wget_swpc_quarter_to_download, "Q");
		strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
		strcat(wget_swpc_quarter_to_download, "_DGD");
		strcat(wget_swpc_quarter_to_download, ".txt >/dev/null 2>&1");
		//printf("WGET --NO-CHECK-CERTIFICATE  <%s>", wget_swpc_quarter_to_download);
		if (do_not_download_file_swpc_or_wget !=1 ){
		system(wget_swpc_quarter_to_download);
		}
		strcpy(list_filename_obs_swpc_ap[ifile_ap], initial_filename_f107_obs);
		strcat(list_filename_obs_swpc_ap[ifile_ap], year_current_str);
		strcat(list_filename_obs_swpc_ap[ifile_ap], "Q");
		strcat(list_filename_obs_swpc_ap[ifile_ap], wget_swpc_quarter_to_download_temp);
		strcat(list_filename_obs_swpc_ap[ifile_ap], "_DGD.txt");
		ifile_ap = ifile_ap + 1;

		//   printf("<%s>\n", wget_swpc_quarter_to_download);
	      }

	    }
	  }
	  else{ // end epoch year is older than current year (so start epoch year too)
	    for (ooo = date_file_obs_initial_swpc_int; ooo < date_file_obs_final_swpc_int + 1; ooo++){
	      // F10.7
	      strcpy(wget_swpc_quarter_to_download, "wget --no-check-certificate  ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/");
	      strcpy(wget_swpc_quarter_to_download_temp, "");
	      sprintf(wget_swpc_quarter_to_download_temp, "%d", ooo );
	      strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
	      strcat(wget_swpc_quarter_to_download, "_DSD.txt -O ");
	      strcat(wget_swpc_quarter_to_download, initial_filename_f107_obs);
	      strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
	      strcat(wget_swpc_quarter_to_download, "_DSD.txt >/dev/null 2>&1");
		//printf("WGET --NO-CHECK-CERTIFICATE  <%s>", wget_swpc_quarter_to_download);
	      if (do_not_download_file_swpc_or_wget !=1 ){
	     system(wget_swpc_quarter_to_download);
	      }
		strcpy(list_filename_obs_swpc_f107[ifile_f107], initial_filename_f107_obs);
		strcat(list_filename_obs_swpc_f107[ifile_f107], wget_swpc_quarter_to_download_temp);
		strcat(list_filename_obs_swpc_f107[ifile_f107], "_DSD.txt");
		ifile_f107 = ifile_f107 + 1;

	    }
	      // Ap
	    for (ooo = date_file_obs_initial_swpc_int_ap; ooo < date_file_obs_final_swpc_int + 1; ooo++){
	      strcpy(wget_swpc_quarter_to_download, "wget --no-check-certificate  ftp://ftp.swpc.noaa.gov/pub/indices/old_indices/");
	      strcpy(wget_swpc_quarter_to_download_temp, "");
	      sprintf(wget_swpc_quarter_to_download_temp, "%d", ooo );
	      strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
	      strcat(wget_swpc_quarter_to_download, "_DGD.txt -O ");
	      strcat(wget_swpc_quarter_to_download, initial_filename_f107_obs);
	      strcat(wget_swpc_quarter_to_download, wget_swpc_quarter_to_download_temp);
	      strcat(wget_swpc_quarter_to_download, "_DGD.txt >/dev/null 2>&1");
		//printf("WGET --NO-CHECK-CERTIFICATE     <%s>", wget_swpc_quarter_to_download);
	      if (do_not_download_file_swpc_or_wget !=1 ){
	     system(wget_swpc_quarter_to_download);  
	      }
		strcpy(list_filename_obs_swpc_ap[ifile_ap], initial_filename_f107_obs);
		strcat(list_filename_obs_swpc_ap[ifile_ap], wget_swpc_quarter_to_download_temp);
		strcat(list_filename_obs_swpc_ap[ifile_ap], "_DGD.txt");
		ifile_ap = ifile_ap + 1;

	      // printf("<%s>\n", wget_swpc_quarter_to_download);

	    }

	  }

	  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
	    printf("----- (load_options) Done downloading the SWPC observation files of F10.7 and Ap.\n");
	  }

	}

	/* for (ooo = 0; ooo < nb_file_obs_swpc; ooo++){ */
	/*   printf("<%s>\n", list_filename_obs_swpc_f107[ooo]); */
	/* } */
	/* /\* //	printf("\n"); *\/ */
	/* for (ooo = 0; ooo < nb_file_obs_swpc_ap; ooo++){ */
	/*   //	  printf("<%s>\n", list_filename_obs_swpc_ap[ooo]); */
	/* } */

	// // Predictions
	if ( OPTIONS->swpc_need_predictions == 1){
	  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
	    printf("----- (load_options) Downloading the SWPC prediction files of F10.7 and Ap.\n");
	  }

	  // // File to download: http://services.swpc.noaa.gov/text/45-day-ap-forecast.txt (always the same name)
	  double et_current_day_plus_45_day;
	  et_current_day_plus_45_day = et_current_day + 45. * 24 * 3600;

	/* exit(0); */

	  char current_day_plus_45_day[256];
	  et2utc_c(et_current_day_plus_45_day, "ISOC", 0, 255, current_day_plus_45_day);
	  if ( ( OPTIONS->et_oldest_tle_epoch >=  et_current_day_plus_45_day + 24*3600 ) || (  et_final_epoch >=  et_current_day_plus_45_day + 24*3600  ) ){ // actually plus 46  because we have the prediction on day + 45
	    OPTIONS->swpc_final_epoch_more_recent_than_45_days_from_current_day = 1;
	    printf("***! You chose to compute the drag with predictions of the solar activity from SWPC. However, these predictions are limited to 45 days in the future (%s). Therefore, for times beyond %s, Ap and F10.7 are set to constant values, equal to the average of the values predicted by SWPC over 45 days.\n", current_day_plus_45_day, current_day_plus_45_day,current_day_plus_45_day);// MPI_Finalize(); exit(0);
	  }


	    strcpy( wget_filename_f107_ap_pred, "wget --no-check-certificate  http://services.swpc.noaa.gov/text/45-day-ap-forecast.txt -O ");
/* 	    strcpy(filename_f107_ap_pred, ""); */
/* 	  //newstructure  */
/* /\* 	    strcat(filename_f107_ap_pred, OPTIONS->dir_input_density_msis); *\/ */
/* /\* 	    strcat(filename_f107_ap_pred, "/"); *\/ */
/* 	  //newstructure  */
/* 	    strcat(filename_f107_ap_pred, "swpc_predictions_f107_ap.txt" ); */
	    strcat(wget_filename_f107_ap_pred, filename_f107_ap_pred);
	    strcat(wget_filename_f107_ap_pred, " >/dev/null 2>&1");
	    //	    printf("PREF <%s>\n", wget_filename_f107_ap_pred);
	    if (do_not_download_file_swpc_or_wget !=1 ){
	    system(wget_filename_f107_ap_pred);
	    }
	  
	  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
	    printf("----- (load_options) Done downloading the SWPC prediction files of F10.7 and Ap.\n");
	  }

	

	 
	if (( iDebugLevel >= 4 ) ){
	    printf("----- (load_options) Allocating memory for mod_f107, mod_f107a, mod_ap, and et_mod_f107_ap (iProc %d)\n", iProc);
	  }
 MPI_Barrier( MPI_COMM_WORLD ) ;
	  if (( iDebugLevel >= 4 ) ){
	    printf("----- (load_options) Went through MPI_Barrier (iProc %d)\n", iProc);
	  }

	OPTIONS->mod_f107 = malloc( ( OPTIONS->nb_time_steps + (int)(( 2 * 24 * 3600. ) / OPTIONS->dt) ) * 2 * sizeof(double) ); // mod_f107 does ahave that many elements. This is a maximum.  + (int)(( 2 * 24 * 3600. ) / OPTIONS->dt) because mod_f107 will be filled with the entire day that includes intial_epoch (24 * 3600 would be enough but we double it to be safe)
      if (  OPTIONS->mod_f107 == NULL ){
	printf("***! Could not allow memory space for  OPTIONS->mod_f107 \n. The program will stop. !***\n");
	ierr =  MPI_Finalize();
	exit(0);
      }

      OPTIONS->mod_ap = malloc( ( OPTIONS->nb_time_steps + (int)(( 2 * 24 * 3600. ) / OPTIONS->dt) ) * 2 * sizeof(double) );
      if (  OPTIONS->mod_ap == NULL ){
	printf("***! Could not allow memory space for  OPTIONS->mod_ap \n. The program will stop. !***\n");
	ierr =  MPI_Finalize();
	exit(0);
      }

      OPTIONS->et_mod_f107_ap = malloc( ( OPTIONS->nb_time_steps + (int)(( 2 * 24 * 3600. ) / OPTIONS->dt) ) * 2 * sizeof(double) );

      if (  OPTIONS->et_mod_f107_ap == NULL ){
	printf("***! Could not allow memory space for  OPTIONS->et_mod_f107_ap \n. The program will stop. !***\n");
	ierr =  MPI_Finalize();
	exit(0);
      }
	  if ( ( iDebugLevel >= 4 ) ){
	    printf("----- (load_options) Done allocating memory for mod_f107, mod_f107a, mod_ap, and et_mod_f107_ap (iProc %d)\n", iProc);
	  }
} // end of if swpc_need_predictions == 1

	lin_interpolate_swpc_mod(OPTIONS->f107, OPTIONS->f107A, OPTIONS->Ap, OPTIONS->et_interpo, OPTIONS->mod_f107, OPTIONS->mod_ap, OPTIONS->et_mod_f107_ap, &OPTIONS->use_ap_hist, list_filename_obs_swpc_f107, list_filename_obs_swpc_ap, filename_f107_ap_pred,  OPTIONS->nb_time_steps * 2, OPTIONS->oldest_tle_epoch, OPTIONS->final_epoch, OPTIONS->dt, 999.9,iDebugLevel, nb_file_obs_swpc, nb_file_obs_swpc_ap, OPTIONS->swpc_need_observations, OPTIONS->swpc_need_predictions, &OPTIONS->swpc_et_first_prediction, iProc, filename_f107_ap_mod);
      // "* 2.0" because of the Runge Kunta orfer 4 method
	




	for (isc  = 0; isc < OPTIONS->nb_satellites_not_including_gps; isc++){

	  if (strcmp(type_orbit_initialisation_temp, "collision_vcm") != 0){ // it's onlyif collision_vcm is 1 that the two sc might not start at tthe same time step 
	    if ( (OPTIONS->et_interpo[0] > OPTIONS->swpc_et_first_prediction) && (OPTIONS->swpc_need_predictions == 1)){
	      previous_index( &OPTIONS->aaa_mod[isc], OPTIONS->et_mod_f107_ap, OPTIONS->et_interpo[0] - OPTIONS->swpc_et_first_prediction, ( OPTIONS->nb_time_steps + (int)(( 2 * 24 * 3600. ) / OPTIONS->dt))); 
	    }
	    else{
	      OPTIONS->aaa_mod[isc] = 0;
	    }
	  }
	  else{ // if computing collision with VCM as input then the two sc might not start at the same epoch
	    if ( (OPTIONS->et_vcm[isc] > OPTIONS->swpc_et_first_prediction) && (OPTIONS->swpc_need_predictions == 1)){
	      previous_index( &OPTIONS->aaa_mod[isc], OPTIONS->et_mod_f107_ap, OPTIONS->et_vcm[isc] - OPTIONS->swpc_et_first_prediction, ( OPTIONS->nb_time_steps + (int)(( 2 * 24 * 3600. ) / OPTIONS->dt))); 
	    }
	    else{
	      OPTIONS->aaa_mod[isc] = 0;
	    }

	  }
	}

	free(list_filename_obs_swpc_f107);
	free(list_filename_obs_swpc_ap);
	// MPI_Finalize();exit(0);
  } // end of if "swpc_mod"

      // if the user chose to input external files (that he/she wrote him/herself) for F10.7, Ap, and F10.7A
      else{

	sscanf(line, "%s",filename_f107);
	//newstructure
/* 	strcpy(density_file_temp, OPTIONS->dir_input_density_msis); */
/* 	strcat(density_file_temp, "/"); */
	strcpy(density_file_temp, "");
	//newstructure
	strcat(density_file_temp, filename_f107);
	strcpy(filename_f107, density_file_temp);
	getline(&line,&len,fp);
	sscanf(line, "%s",filename_ap);
	//newstructure
/* 	strcpy(density_file_temp, OPTIONS->dir_input_density_msis); */
/* 	strcat(density_file_temp, "/"); */
	strcpy(density_file_temp, "");
	//newstructure
	strcat(density_file_temp, filename_ap);
	strcpy(filename_ap, density_file_temp);

	OPTIONS->use_ap_hist = 0;
	lin_interpolate(OPTIONS->f107, OPTIONS->f107A, OPTIONS->Ap, OPTIONS->Ap_hist, OPTIONS->et_interpo, &OPTIONS->use_ap_hist, filename_f107, filename_ap, "dynamic_manual", OPTIONS->nb_time_steps * 2, OPTIONS->initial_epoch,OPTIONS->et_oldest_tle_epoch, OPTIONS->final_epoch, OPTIONS->dt, 999.9,iDebugLevel, iProc);       // "* 2.0" because of the Runge Kunta orfer 4 method	// used to be (it all depends on the format of the files the user want to use (see the function lin_interpolate for more detail)):	lin_interpolate(OPTIONS->f107, OPTIONS->f107A, OPTIONS->Ap, OPTIONS->Ap_hist, OPTIONS->et_interpo, &OPTIONS->use_ap_hist, filename_f107, filename_ap, "dynamic_manual", OPTIONS->nb_time_steps * 2, OPTIONS->initial_epoch, OPTIONS->final_epoch, OPTIONS->dt, 999.9,iDebugLevel, iProc);       // "* 2.0" because of the Runge Kunta orfer 4 method


      }

    }
    
    else if (  strcmp(OPTIONS->format_density_driver, "density_file") == 0 ){
      printf("***! Sorry, the option 'density_file' (in section #FORCES) is not available anymore. The program will stop. !***\n"); MPI_Finalize(); exit(0);
      /* getline(&line, &len, fp); */
      /* sscanf(line,"%s", density_file); */
      /* strcpy(density_file_temp, OPTIONS->dir_input_density_density_files); */
      /* strcat(density_file_temp, "/"); */
      /* strcat(density_file_temp, density_file); */
      /* strcpy(density_file, density_file_temp);	 */

      /* // !!!!!!!need to create the density file from python and then go back to the line right below here */
      /* OPTIONS->density = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) ); */

      /* if (  OPTIONS->density == NULL ){ */
      /* 	printf("***! Could not allow memory space for  OPTIONS->density \n. The program will stop. !***\n"); */
      /* 	ierr =  MPI_Finalize(); */
      /* 	exit(0); */
      /* } */
      /* lin_interpolate(OPTIONS->density, OPTIONS->et_interpo, density_file, "raid3", OPTIONS->nb_time_steps * 2, OPTIONS->initial_epoch, OPTIONS->dt, 999.9); */
      /* getline(&line, &len, fp); */
      /* getline(&line, &len, fp); */

    }

    else if ( strcmp(OPTIONS->format_density_driver, "static") == 0 ){
      OPTIONS->use_ap_hist = 0; // use Ap daily in NRLMSIS00e
      getline(&line, &len, fp);
      sscanf(line,"%lf",&OPTIONS->f107_static);

      getline(&line, &len, fp);
      sscanf(line,"%lf",&OPTIONS->f107A_static);
      //      OPTIONS->f107A_static = OPTIONS->f107_static;
      getline(&line, &len, fp);
      sscanf(line,"%lf",&OPTIONS->Ap_static);
      for (hhh = 0; hhh < 7; hhh++){
	sscanf(line,"%lf",&OPTIONS->Ap_hist_static[hhh]);
      }
    }

    // the user chooses GITM to calculate the density at the positions of the sc
    else if ( strcmp(OPTIONS->format_density_driver, "gitm") == 0 ){

      // Read all the file names in the GITM data and store the date in an array
      // IMPORTANT assumptions: the GITM files have to be named with the format as in this example: 3DALL_t150318_044500.bin
	
      DIR *dir;
      struct dirent *ent;
      char extension_gitm_file[256];
      char *next;
      int find_extension;
      char start_filename[256];
      char date_gitm_file[256];
      double et_gitm_file;

      getline(&line, &len, fp);
      sscanf(line,"%s",OPTIONS->gitm_directory);
  

      if ((dir = opendir (OPTIONS->gitm_directory)) != NULL) {

	while ((ent = readdir (dir)) != NULL) {

	  // find the extension of each file. We only care about the file if it is a .bin file
	  strcpy(extension_gitm_file, "");
	  next = &(ent->d_name)[0];
	  find_extension = (int)(strchr(next, '.') - next);
	  strncat(extension_gitm_file, next+find_extension,4);

	  if (strcmp(extension_gitm_file, ".bin") == 0) { // We only care about the file if it is a .bin file
      
	    // We only care about the file if its name starts with "3DALL_t" (see IMPORTANT assumption above)

	    strcpy(start_filename, "");
	    strncat(start_filename, &(ent->d_name)[0], 7);

	    if (strcmp(start_filename, "3DALL_t") == 0) { // We only care about the file if its name starts with "3 */DALL_t" (see IMPORTANT assumption above)

	      // Convert the date format of the GITM file name to the Julian Epoch format and arrange the array in the ascending order of times (older to newer)
	      strcpy(date_gitm_file, "20");
	      strncat(date_gitm_file, &(ent->d_name)[7],2);
	      strcat(date_gitm_file, "-");
	      strncat(date_gitm_file, &(ent->d_name)[9],2);
	      strcat(date_gitm_file, "-");
	      strncat(date_gitm_file, &(ent->d_name)[11],2);
	      strcat(date_gitm_file, "T");
	      strncat(date_gitm_file, &(ent->d_name)[14],2);
	      strcat(date_gitm_file, ":");
	      strncat(date_gitm_file, &(ent->d_name)[16],2);
	      strcat(date_gitm_file, ":");
	      strncat(date_gitm_file, &(ent->d_name)[18],2);
	      str2et_c(date_gitm_file, &et_gitm_file);

	      // arrange in ascending order of time (older to newer)
	      if (i_gitm_file > 0){
		have_found_older_time_than_me = 0;
		count_have_found_older_time_than_me = 0;
		while ( (have_found_older_time_than_me == 0) && (count_have_found_older_time_than_me < i_gitm_file)){
		  if ( ( et_gitm_file < OPTIONS->array_gitm_date[i_gitm_file-1-count_have_found_older_time_than_me][0] ) && (i_gitm_file-1-count_have_found_older_time_than_me > 0) ){
		    have_found_older_time_than_me = 0;
		    count_have_found_older_time_than_me = count_have_found_older_time_than_me+1;
		  }
		  else if ( (i_gitm_file-1-count_have_found_older_time_than_me == 0) && ( et_gitm_file < OPTIONS->array_gitm_date[i_gitm_file-1-count_have_found_older_time_than_me][0] ) ){
		    for (nnn = i_gitm_file; nnn >= 1; nnn--){
		      OPTIONS->array_gitm_date[nnn][0] = OPTIONS->array_gitm_date[nnn-1][0] ;
		      OPTIONS->array_gitm_date[nnn][1] = OPTIONS->array_gitm_date[nnn-1][1] ;
		      strcpy(OPTIONS->array_gitm_file[nnn], OPTIONS->array_gitm_file[nnn-1]); 
		    }
		    OPTIONS->array_gitm_date[0][0] = et_gitm_file;
		    OPTIONS->array_gitm_date[0][1] = i_gitm_file-1-count_have_found_older_time_than_me+1;
		    strcpy(OPTIONS->array_gitm_file[0], OPTIONS->gitm_directory);
		    strcat(OPTIONS->array_gitm_file[0], ent->d_name);

		    have_found_older_time_than_me = 1;		
		  }
		  else{
		    for (nnn = i_gitm_file; nnn > i_gitm_file-1-count_have_found_older_time_than_me+1; nnn--){
		      OPTIONS->array_gitm_date[nnn][0] = OPTIONS->array_gitm_date[nnn-1][0] ;
		      OPTIONS->array_gitm_date[nnn][1] = OPTIONS->array_gitm_date[nnn-1][1] ;
		      strcpy(OPTIONS->array_gitm_file[nnn], OPTIONS->array_gitm_file[nnn-1]); 
		    }
		    OPTIONS->array_gitm_date[i_gitm_file-1-count_have_found_older_time_than_me+1][0] = et_gitm_file;
		    OPTIONS->array_gitm_date[i_gitm_file-1-count_have_found_older_time_than_me+1][1] = i_gitm_file-1-count_have_found_older_time_than_me+1;
		    strcpy(OPTIONS->array_gitm_file[i_gitm_file-1-count_have_found_older_time_than_me+1], OPTIONS->gitm_directory);
		    strcat(OPTIONS->array_gitm_file[i_gitm_file-1-count_have_found_older_time_than_me+1], ent->d_name);

		    have_found_older_time_than_me = 1;
		  }
		}
	      }
	      else{
		OPTIONS->array_gitm_date[0][0] = et_gitm_file; 
		OPTIONS->array_gitm_date[0][1] = i_gitm_file; 

		strcpy(OPTIONS->array_gitm_file[0], OPTIONS->gitm_directory);
		strcat(OPTIONS->array_gitm_file[0], ent->d_name);
	      }
	      // end of arrange in ascending order 
	      i_gitm_file = i_gitm_file + 1;
	    } // end of we only care about the file if its name starts with "3DALL_t" (see IMPORTANT assumption above)  

	  } // end of we only care about the file if it is a .bin file

	}

	OPTIONS->nb_gitm_file = i_gitm_file;
	closedir (dir);
      } else {
	/* could not open directory */
	printf("***! Could not open the GITM data directory. The program will stop. !***\n");
	ierr =  MPI_Finalize();
	exit(0);
      }

      getline(&line, &len, fp);
      getline(&line, &len, fp);

    }

    else {
      printf("***! You need to choose between 'static', 'dynamic', 'density_file', and 'gitm' for the density drivers. The program will stop. !*** \n");
      ierr =  MPI_Finalize();
      exit(0);
    }
  }
  else{
    getline(&line, &len, fp);
    getline(&line, &len, fp);
    getline(&line, &len, fp);

  }

  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) End of section #FORCES.\n");
  }


  /* GROUND_STATIONS */
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section #GROUND_STATIONS.\n");
  }
  rewind(fp);
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#GROUND_STATIONS", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  if (feof(fp)){
    OPTIONS->nb_ground_stations = 0;
  }
  else{
    FILE *file_ground_station=NULL;
    // Get the name of the file that has the ground stations
    getline(&line, &len, fp);
    next = &line[0];
    strcpy(text_output,"");
    find_file_name = (int)(strchr(next, ' ') - next)-1;
    strncat(text_output, next, find_file_name);
    strtok(text_output, "\n");strtok(text_output, "\r");
    //newstructure
/*     strcpy(OPTIONS->filename_ground_station, OPTIONS->dir_input_coverage_ground_station); */
/*     strcat( OPTIONS->filename_ground_station, "/"); */
    strcpy(OPTIONS->filename_ground_station, "");
    //newstructure
    strcat( OPTIONS->filename_ground_station, text_output);
    strtok(OPTIONS->filename_ground_station, "\n");strtok(OPTIONS->filename_ground_station, "\r");  
    //newstructure
/*     strcpy(str_no_gs, OPTIONS->dir_input_coverage_ground_station); */
/*     strcat( str_no_gs, "/"); */
    strcpy(str_no_gs, "");
    //newstructure
    strcat( str_no_gs, "0");
    if (strcmp( OPTIONS->filename_ground_station, str_no_gs ) == 0){
      OPTIONS->nb_ground_stations = 0;
    }
    else{
      // Read this file and loads the ground stations (name, position, min elevation angle)
      file_ground_station = fopen(OPTIONS->filename_ground_station, "r");
      if (file_ground_station == NULL){
	printf("***! Could not open the ground station file %s\n. The program will stop. !***\n", OPTIONS->filename_ground_station);
      }
      found_eoh = 0;
      while ( found_eoh == 0 && !feof(file_ground_station)) {
	getline(&line, &len, file_ground_station);
	sscanf(line, "%s", text);
	if (  strcmp( "#START", text  ) == 0 )  {
	  found_eoh = 1;
	}
      }
      found_eoh = 0;
      OPTIONS->nb_ground_stations = 0;
      while ( found_eoh == 0 && !feof(file_ground_station)) {
	getline(&line, &len, file_ground_station);
	sscanf(line, "%s", text);
	if (  strcmp( "#END", text  ) != 0 )  {
	  sscanf(line, "%s %lf %lf %lf %lf", OPTIONS->name_ground_station[OPTIONS->nb_ground_stations], &OPTIONS->latitude_ground_station[OPTIONS->nb_ground_stations], &OPTIONS->longitude_ground_station[OPTIONS->nb_ground_stations], &OPTIONS->altitude_ground_station[OPTIONS->nb_ground_stations], &OPTIONS->min_elevation_angle_ground_station[OPTIONS->nb_ground_stations]);
	  OPTIONS->nb_ground_stations = OPTIONS->nb_ground_stations + 1;
	}
	else{
	  found_eoh = 1;
	}
      }
      fclose(file_ground_station);
    }
  }
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) End of section #GROUND_STATIONS.\n");
  }

  /* STORMS */
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section #STORMS.\n");
  }
  rewind(fp);
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#STORMS", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  if (feof(fp)){
    OPTIONS->nb_storm = 0;
  }
  else{
    getline(&line,&len,fp);

    RemoveSpaces(line);  strtok(line, "\n");  strtok(line, "\r"); 
    if ( line[0] == '0' ){
      OPTIONS->nb_storm = 0;
    }
    else{
      char * line_copy = malloc(strlen(line) + 1); 
      strcpy(line_copy, line);

      for (nb_storm_temp=0; line_copy[nb_storm_temp]; line_copy[nb_storm_temp]==',' ? nb_storm_temp++ : *line_copy++); 
      nb_storm_temp = nb_storm_temp + 1;
      OPTIONS->nb_storm = nb_storm_temp;
      next = &line[0];
      for (sss = 0; sss<OPTIONS->nb_storm; sss++){
	strcpy(text_output,"");
	if (sss < OPTIONS->nb_storm - 1 ){
	  find_file_name = (int)(strchr(next, ',') - next);
	  strncat(text_output, next, find_file_name);
	  strcpy( OPTIONS->filename_storm[sss], text_output);
	  strtok(OPTIONS->filename_storm[sss], "\n");
	  strtok(OPTIONS->filename_storm[sss], "\r");
	  next = strstr(next, ",") +1;
	  find_file_name_sum = find_file_name_sum + find_file_name;
	}
	else{
	  strncat(text_output, next, (int)(strlen(line)) - find_file_name_sum);
	  strcpy( OPTIONS->filename_storm[sss], text_output);
	  strtok(OPTIONS->filename_storm[sss], "\n");
	  strtok(OPTIONS->filename_storm[sss], "\r");
      
	}

      }
    }

  }

  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) End of section #STORMS.\n");
  }



  // calculates the total number of ensembles
  OPTIONS->nb_ensembles_min = 0;
  int nb_ensemble_array[4] = {OPTIONS->nb_ensembles, OPTIONS->nb_ensembles_attitude, OPTIONS->nb_ensembles_cd, OPTIONS->nb_ensembles_density}; // COE, attitude, Cd, density
  for (sss = 0; (unsigned)sss < sizeof(nb_ensemble_array)/sizeof(nb_ensemble_array[0]); sss++){
    if ( nb_ensemble_array[sss] != 0 ){
      OPTIONS->nb_ensembles_min = nb_ensemble_array[sss];
    }
  }
  for (sss = 0; (unsigned)sss < sizeof(nb_ensemble_array)/sizeof(nb_ensemble_array[0]); sss++){
    if ( ( nb_ensemble_array[sss] > 0 ) && ( nb_ensemble_array[sss] < OPTIONS->nb_ensembles_min ) ){
      OPTIONS->nb_ensembles_min = nb_ensemble_array[sss];
    }
  }

  /* if (OPTIONS->nb_ensembles > 0){ */
  /*   OPTIONS->nb_ensembles_min = OPTIONS->nb_ensembles; */
  /*   if ( (OPTIONS->nb_ensembles_attitude < OPTIONS->nb_ensembles) && (OPTIONS->nb_ensembles_attitude>0)){ */
  /*     OPTIONS->nb_ensembles_min = OPTIONS->nb_ensembles_attitude; */
  /*   } */
  /* } */

  /* else if (OPTIONS->nb_ensembles_attitude > 0){ */
  /*   OPTIONS->nb_ensembles_min = OPTIONS->nb_ensembles_attitude; */
  /* } */

  /* else if (OPTIONS->nb_ensembles_cd > 0){ */
  /*   OPTIONS->nb_ensembles_min = OPTIONS->nb_ensembles_cd; */
  /* } */

  int nProcs_that_are_gonna_run_ensembles;
  nProcs_that_are_gonna_run_ensembles = nProcs;
  if ( ( nProcs > OPTIONS->nb_ensembles_min ) && (  OPTIONS->nb_ensembles_min  > 0 ) ) {
    nProcs_that_are_gonna_run_ensembles = OPTIONS->nb_ensembles_min;
  }

  if (nProcs > 0){
    OPTIONS->nb_ensemble_min_per_proc = (int)(OPTIONS->nb_ensembles_min / nProcs_that_are_gonna_run_ensembles);
  }

  /* OUTPUT_ENSEMBLES */
  if ((iProc == 0) & ( iDebugLevel >= 2 )){
    printf("--- (load_options) Beginning of section #OUTPUT_ENSEMBLES.\n");
  }

  if (OPTIONS->nb_ensembles_min > 0){
    // quick test: the numbers of proc should be < number of ensembles (i'll fix that some day)
/*     if (nProcs > OPTIONS->nb_ensembles_min){ */
      
/*       print_error_any_iproc(iProc, "The number of processors used in MPI has to be smaller than (or equal to) the number of ensembles (this bug will be fixed in the future)"); */
/*     } */
    rewind(fp);
    found_eoh = 0;
    while ( found_eoh == 0 && !feof(fp)) {
      getline(&line, &len, fp);
      sscanf(line, "%s", text);
      if (  strcmp( "#OUTPUT_ENSEMBLES", text  ) == 0 )  {
	found_eoh = 1;
      }
    }
    if (feof(fp)){
      if (iProc == 0){
      printf("***! You chose to run ensembles but no section #OUTPUT_ENSEMBLES was found in %s. !***\n", filename); //MPI_Finalize(); exit(0);
      }
      OPTIONS->nb_ensembles_output_files = 0;
    }
    else{
    // Names of output files
    char filename_output_ensemble_temp[20][1000];
    getline(&line,&len,fp);
    next = &line[0];
    RemoveSpaces(next);
    more_ensembles_output = 1;
    new_next = next;
    while (more_ensembles_output == 1){
      find_comments = (int)(strchr(new_next, '/') - new_next);
      if (find_comments < 0){
	find_comments = 1000000; // if there is no comment then to get in the next if condition we give find_comments a big number
      }
      strcpy(text_output,"");
      find_file_name = (int)(strchr(new_next, ',') - new_next);
      if ( (find_file_name < find_comments ) & (find_file_name > 0)){
	strncat(text_output, new_next, find_file_name);
	strcpy( filename_output_ensemble_temp[inc_ensembles_output], text_output);      
	new_next = new_next + find_file_name+1;
	inc_ensembles_output = inc_ensembles_output + 1;
      }
      else if (find_comments > 0) {
	strncat(text_output, new_next, find_comments);
	strcpy(filename_output_ensemble_temp[inc_ensembles_output], text_output);
	inc_ensembles_output = inc_ensembles_output + 1;
	more_ensembles_output = 0;
      }

      else{
	more_ensembles_output = 0;
      }
    }

    for (sss = 0; sss < inc_ensembles_output; sss++){

      RemoveSpaces(filename_output_ensemble_temp[sss]);
      strtok(filename_output_ensemble_temp[sss], "\n");strtok(filename_output_ensemble_temp[sss], "\r");

      if (strcmp(filename_output_ensemble_temp[sss], "eci_r") == 0){
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "x_eci");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "y_eci");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "z_eci");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
      }
      else if (strcmp(filename_output_ensemble_temp[sss], "eci_v") == 0){
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "vx_eci");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "vy_eci");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "vz_eci");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
      }
      else if (strcmp(filename_output_ensemble_temp[sss], "geodetic") == 0){
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "latitude");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "longitude");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "altitude");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
      }
      else if (strcmp(filename_output_ensemble_temp[sss], "attitude") == 0){
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "pitch");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "roll");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "yaw");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
      }
      else if (strcmp(filename_output_ensemble_temp[sss], "power") == 0){
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "power");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
      }
      else if (strcmp(filename_output_ensemble_temp[sss], "oe") == 0){ // sma, inc, ecc, true_ano, raan, arg_per
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "sma");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "inclination");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "eccentricity");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "true_anomaly");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "RAAN");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "argument_perigee");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;

      }
    
      else if (strcmp(filename_output_ensemble_temp[sss], "density") == 0){ 
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "rho");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "f107");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "f107a");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "ap");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
      }

if ((strcmp(filename_output_ensemble_temp[sss], "collision") == 0) || (strcmp(filename_output_ensemble_temp[sss], "collision_vcm") == 0)){
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "tca");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "dca");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
/* 	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "sample"); */
/* 	inc_ensembles_output_final = inc_ensembles_output_final + 1; */

	OPTIONS->write_collision_files = 1;
      }

      if (strcmp(filename_output_ensemble_temp[sss], "cd") == 0){
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "cd");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;

      }



      if (strcmp(filename_output_ensemble_temp[sss], "srp") == 0){
	strcpy( OPTIONS->filename_output_ensemble[inc_ensembles_output_final], "srp");
	inc_ensembles_output_final = inc_ensembles_output_final + 1;
      }



    }
    OPTIONS->nb_ensembles_output_files = inc_ensembles_output_final;
    }
  }
  
  else{
    OPTIONS->nb_ensembles_output_files = 0;
  }

  if ((iProc == 0) & ( iDebugLevel >= 2 ) ){
    printf("--- (load_options) End of section #OUTPUT_ENSEMBLES.\n");
  }


  free(line);    
  fclose(fp);

  //  make sure all output directories exist. If they don't, then create them. 
  //  OUTPUT

  strcpy(OPTIONS->dir_output, path_output_run_name_temp);

  //newstructure
  //  strcat(OPTIONS->dir_output,"spock");

  //  strcat(OPTIONS->dir_output,"output");
  //newstructure
  err = stat(OPTIONS->dir_output, &s);
  if (err !=-1){
    if(S_ISDIR(s.st_mode)) { // the dir exists
      //      printf("AAA\nAAA\nAAA\n");
    } 
    else { // exists but it is not a dir
      //    printf("BBB\nBBB\nBBB\n");
      strcat(OPTIONS->dir_output,"_dir");
      mode_t process_mask = umask(0);
      mkdir(OPTIONS->dir_output, S_IRWXU | S_IRWXG | S_IRWXO);
      umask(process_mask);
    }
  }
  else{// the dir does not exist
    //  printf("CCC\nCCC\nCCC\n");
    mode_t process_mask = umask(0);
    mkdir(OPTIONS->dir_output, S_IRWXU | S_IRWXG | S_IRWXO);

    umask(process_mask);
  }
  //  MPI_Finalize();exit(0);
  //  OUTPUT_RUN_NAME

  strcpy(OPTIONS->dir_output_run_name_temp, dir_output_run_name_temp);
  strcpy(OPTIONS->dir_output_run_name, OPTIONS->dir_output);
  //  strcat(OPTIONS->dir_output_run_name,"/");
  strcat(OPTIONS->dir_output_run_name,dir_output_run_name_temp);
  err = stat(OPTIONS->dir_output_run_name, &s);
  if (err !=-1){
    if(S_ISDIR(s.st_mode)) { // the dir exists
    } 
    else { // exists but it is not a dir
      strcat(OPTIONS->dir_output_run_name,"_dir");
      mode_t process_mask = umask(0);
      mkdir(OPTIONS->dir_output_run_name, S_IRWXU | S_IRWXG | S_IRWXO);
      umask(process_mask);
    }
  }
  else{// the dir does not exist
    mode_t process_mask = umask(0);
    mkdir(OPTIONS->dir_output_run_name, S_IRWXU | S_IRWXG | S_IRWXO);
    umask(process_mask);
  }

  //  OUTPUT_RUN_NAME_SAT_NAME
  for (sss = 0; sss < OPTIONS->n_satellites ; sss++){
    if (sss < OPTIONS->n_satellites - OPTIONS->nb_gps){
      strcpy(OPTIONS->dir_output_run_name_sat_name[sss], OPTIONS->dir_output_run_name);

      strcat(OPTIONS->dir_output_run_name_sat_name[sss],"/");
      next = OPTIONS->filename_output[sss];
      find_file_name =  (int)(strrchr(next, '.') - next);
      strncat(OPTIONS->dir_output_run_name_sat_name[sss], next, find_file_name);

    }
    else{
      strcpy(OPTIONS->dir_output_run_name_sat_name[sss], OPTIONS->dir_output_run_name);
      strcat(OPTIONS->dir_output_run_name_sat_name[sss], "/constellation_GPS/");
      err = stat(OPTIONS->dir_output_run_name_sat_name[sss], &s);
      if (err !=-1){
	if(S_ISDIR(s.st_mode)) { // the dir exists
	} 
	else { // exists but it is not a dir
	  strcat(OPTIONS->dir_output_run_name_sat_name[sss],"_dir");
	  mode_t process_mask = umask(0);
	  mkdir(OPTIONS->dir_output_run_name_sat_name[sss], S_IRWXU | S_IRWXG | S_IRWXO);
	  umask(process_mask);
	}
      }
      else{// the dir does not exist
	mode_t process_mask = umask(0);
	mkdir(OPTIONS->dir_output_run_name_sat_name[sss], S_IRWXU | S_IRWXG | S_IRWXO);
	umask(process_mask);
      }

      next = OPTIONS->filename_output[sss];
      find_file_name =  (int)(strchr(next, '.') - next);
      strncat(OPTIONS->dir_output_run_name_sat_name[sss], next, find_file_name);
    }

    err = stat(OPTIONS->dir_output_run_name_sat_name[sss], &s);
    if (err !=-1){
      if(S_ISDIR(s.st_mode)) { // the dir exists
      } 
      else { // exists but it is not a dir
	strcat(OPTIONS->dir_output_run_name_sat_name[sss],"_dir");
	mode_t process_mask = umask(0);
	mkdir(OPTIONS->dir_output_run_name_sat_name[sss], S_IRWXU | S_IRWXG | S_IRWXO);
	umask(process_mask);
      }
    }
    else{// the dir does not exist
      mode_t process_mask = umask(0);
      mkdir(OPTIONS->dir_output_run_name_sat_name[sss], S_IRWXU | S_IRWXG | S_IRWXO);
      umask(process_mask);
    }
  }


  //  OUTPUT_RUN_NAME_COLLISION
  if (OPTIONS->write_collision_files == 1){   
    strcpy(OPTIONS->dir_output_run_name_collision, OPTIONS->dir_output_run_name);
    strcat(OPTIONS->dir_output_run_name_collision, "/collision");
    err = stat(OPTIONS->dir_output_run_name_collision, &s);
    if (err !=-1){
      if(S_ISDIR(s.st_mode)) { // the dir exists
      } 
      else { // exists but it is not a dir
	strcat(OPTIONS->dir_output_run_name_collision,"_dir");
	mode_t process_mask = umask(0);
	mkdir(OPTIONS->dir_output_run_name_collision, S_IRWXU | S_IRWXG | S_IRWXO);
	umask(process_mask);
      }
    }
    else{// the dir does not exist
      mode_t process_mask = umask(0);
      mkdir(OPTIONS->dir_output_run_name_collision, S_IRWXU | S_IRWXG | S_IRWXO);
      umask(process_mask);
    }
    
  
  //  OUTPUT_RUN_NAME_COLLISION_tca
    strcpy(OPTIONS->dir_output_run_name_collision_tca, OPTIONS->dir_output_run_name_collision);
    strcat(OPTIONS->dir_output_run_name_collision_tca, "/tca");
    err = stat(OPTIONS->dir_output_run_name_collision_tca, &s);
    if (err !=-1){
      if(S_ISDIR(s.st_mode)) { // the dir exists
      } 
      else { // exists but it is not a dir
	strcat(OPTIONS->dir_output_run_name_collision_tca,"_dir");
	mode_t process_mask = umask(0);
	mkdir(OPTIONS->dir_output_run_name_collision_tca, S_IRWXU | S_IRWXG | S_IRWXO);
	umask(process_mask);
      }
    }
    else{// the dir does not exist
      mode_t process_mask = umask(0);
      mkdir(OPTIONS->dir_output_run_name_collision_tca, S_IRWXU | S_IRWXG | S_IRWXO);
      umask(process_mask);
    }
    
  

   

  //  OUTPUT_RUN_NAME_COLLISION_dca

    strcpy(OPTIONS->dir_output_run_name_collision_dca, OPTIONS->dir_output_run_name_collision);
    strcat(OPTIONS->dir_output_run_name_collision_dca, "/dca");
    err = stat(OPTIONS->dir_output_run_name_collision_dca, &s);
    if (err !=-1){
      if(S_ISDIR(s.st_mode)) { // the dir exists
      } 
      else { // exists but it is not a dir
	strcat(OPTIONS->dir_output_run_name_collision_dca,"_dir");
	mode_t process_mask = umask(0);
	mkdir(OPTIONS->dir_output_run_name_collision_dca, S_IRWXU | S_IRWXG | S_IRWXO);
	umask(process_mask);
      }
    }
    else{// the dir does not exist
      mode_t process_mask = umask(0);
      mkdir(OPTIONS->dir_output_run_name_collision_dca, S_IRWXU | S_IRWXG | S_IRWXO);
      umask(process_mask);
    }
    
  
  
/*   //  OUTPUT_RUN_NAME_COLLISION_sample */

/*     strcpy(OPTIONS->dir_output_run_name_collision_sample, OPTIONS->dir_output_run_name_collision); */
/*     strcat(OPTIONS->dir_output_run_name_collision_sample, "/sample"); */
/*     err = stat(OPTIONS->dir_output_run_name_collision_sample, &s); */
/*     if (err !=-1){ */
/*       if(S_ISDIR(s.st_mode)) { // the dir exists */
/*       }  */
/*       else { // exists but it is not a dir */
/* 	strcat(OPTIONS->dir_output_run_name_collision_sample,"_dir"); */
/* 	mode_t process_mask = umask(0); */
/* 	mkdir(OPTIONS->dir_output_run_name_collision_sample, S_IRWXU | S_IRWXG | S_IRWXO); */
/* 	umask(process_mask); */
/*       } */
/*     } */
/*     else{// the dir does not exist */
/*       mode_t process_mask = umask(0); */
/*       mkdir(OPTIONS->dir_output_run_name_collision_sample, S_IRWXU | S_IRWXG | S_IRWXO); */
/*       umask(process_mask); */
/*     } */
    
  }
  
   
  //  OUTPUT_RUN_NAME_COVERAGE
  if (OPTIONS->nb_storm > 0){   
    strcpy(OPTIONS->dir_output_run_name_coverage, OPTIONS->dir_output_run_name);
    strcat(OPTIONS->dir_output_run_name_coverage, "/coverage");
    err = stat(OPTIONS->dir_output_run_name_coverage, &s);
    if (err !=-1){
      if(S_ISDIR(s.st_mode)) { // the dir exists
      } 
      else { // exists but it is not a dir
	strcat(OPTIONS->dir_output_run_name_coverage,"_dir");
	mode_t process_mask = umask(0);
	mkdir(OPTIONS->dir_output_run_name_coverage, S_IRWXU | S_IRWXG | S_IRWXO);
	umask(process_mask);
      }
    }
    else{// the dir does not exist
      mode_t process_mask = umask(0);
      mkdir(OPTIONS->dir_output_run_name_coverage, S_IRWXU | S_IRWXG | S_IRWXO);
      umask(process_mask);
    }
    
  }



  //  OUTPUT_RUN_NAME_COVERAGE_STORM
  if (OPTIONS->nb_storm > 0){   
    strcpy(OPTIONS->dir_output_run_name_coverage_storm, OPTIONS->dir_output_run_name_coverage);
    strcat(OPTIONS->dir_output_run_name_coverage_storm, "/storm");
    err = stat(OPTIONS->dir_output_run_name_coverage_storm, &s);
    if (err !=-1){
      if(S_ISDIR(s.st_mode)) { // the dir exists
      } 
      else { // exists but it is not a dir
	strcat(OPTIONS->dir_output_run_name_coverage_storm,"_dir");
	mode_t process_mask = umask(0);
	mkdir(OPTIONS->dir_output_run_name_coverage_storm, S_IRWXU | S_IRWXG | S_IRWXO);
	umask(process_mask);
      }
    }
    else{// the dir does not exist
      mode_t process_mask = umask(0);
      mkdir(OPTIONS->dir_output_run_name_coverage_storm, S_IRWXU | S_IRWXG | S_IRWXO);
      umask(process_mask);
    }
    
  }

  //  OUTPUT_RUN_NAME_SAT_NAME_ENSEMBLE

  for (sss = 0; sss < OPTIONS->n_satellites - OPTIONS->nb_gps; sss++){
    strcpy(OPTIONS->dir_output_run_name_sat_name_ensemble, OPTIONS->dir_output_run_name_sat_name[sss]);
    strcat(OPTIONS->dir_output_run_name_sat_name_ensemble,"/ensemble");
      
    err = stat(OPTIONS->dir_output_run_name_sat_name_ensemble, &s);
    if (err !=-1){
      if(S_ISDIR(s.st_mode)) { // the dir exists
      } 
      else { // exists but it is not a dir
	strcat(OPTIONS->dir_output_run_name_sat_name_ensemble,"_dir");
	mode_t process_mask = umask(0);
	mkdir(OPTIONS->dir_output_run_name_sat_name_ensemble, S_IRWXU | S_IRWXG | S_IRWXO);
	umask(process_mask);
      }
    }
    else{// the dir does not exist
      mode_t process_mask = umask(0);
      mkdir(OPTIONS->dir_output_run_name_sat_name_ensemble, S_IRWXU | S_IRWXG | S_IRWXO);
      umask(process_mask);
    }
  } 


  //  OUTPUT_RUN_NAME_SAT_NAME_COVERAGE
  for (sss = 0; sss < OPTIONS->n_satellites; sss++){
    strcpy(OPTIONS->dir_output_run_name_sat_name_coverage[sss], OPTIONS->dir_output_run_name_sat_name[sss]);
    strcat(OPTIONS->dir_output_run_name_sat_name_coverage[sss],"/coverage");
      
    err = stat(OPTIONS->dir_output_run_name_sat_name_coverage[sss], &s);
    if (err !=-1){
      if(S_ISDIR(s.st_mode)) { // the dir exists
      } 
      else { // exists but it is not a dir
	strcat(OPTIONS->dir_output_run_name_sat_name_coverage[sss],"_dir");
	mode_t process_mask = umask(0);
	mkdir(OPTIONS->dir_output_run_name_sat_name_coverage[sss], S_IRWXU | S_IRWXG | S_IRWXO);
	umask(process_mask);
      }
    }
    else{// the dir does not exist
      mode_t process_mask = umask(0);
      mkdir(OPTIONS->dir_output_run_name_sat_name_coverage[sss], S_IRWXU | S_IRWXG | S_IRWXO);
      umask(process_mask);
    }
  } 
  // end of make sure all output directories exist. If they don't, then create them. 


  return 0;
  
 }


/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           load_surface
//  Purpose:        Load the surface file (file name with the properties of each surface: Cd, area (cm^2), solar cell area (cm^2), the normal to the surface, the specular and diffusive reflectivities. Can also be the name of a model)
//  Assumptions:    Same properties for each satellite of the constellation
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 08/11/2015    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////


int load_surface( OPTIONS_T *OPTIONS,
		  char      filename[1000],
		  int       nProcs){

  /* Declarations */

  double cos_surf;
  double cos_limit = 0.9998476951563913;
	double cd_or_acco;
	int sss, is_surf_eff;
int found_new_surface;
  double temp_normal[3];
  FILE *fp = NULL;
  int found_eoh = 0;
  char string_nb_word[256];
  ssize_t read;
  char *line = NULL;
  int surface_counter;
  char str_len[256];
  size_t len = 0;
  char text[256];
  OPTIONS->n_surfaces_eff = 0;

  /* Algorithm */
  // Open geometry file
  strtok(filename, "\n");
  strtok(filename, "\r");
  fp = fopen(filename, "r");
  if (fp == NULL){
    printf("!***!\nThe geometry file:\n%s\ncould not be found. The program will stop.\n!***!\n", filename); MPI_Finalize(); exit(0);
  }
  // Find the number of ensembles on Cd in section #NB_ENSEMBLES_CD. If there is no such section, that's ok (the number of ensembles is set to 0)
  int found_ensemble_cd = 0;
  while ( found_ensemble_cd == 0 && !feof(fp) ){
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#NB_ENSEMBLES_CD", text  ) == 0 )  {
      found_ensemble_cd = 1;
    }
  }
  if (feof(fp)){
    OPTIONS->nb_ensembles_cd = 0;
    OPTIONS->nb_ensemble_cd_per_proc = 0;
  }
  else{
    // Number of ensembles on Cd
    getline(&line, &len, fp);
    OPTIONS->nb_ensembles_cd = 0;
    sscanf(line, "%d", &OPTIONS->nb_ensembles_cd);
    if (nProcs > 0){
      OPTIONS->nb_ensemble_cd_per_proc = (int)(OPTIONS->nb_ensembles_cd / nProcs);
    }
  }

  // Now read the info on the surfaces
  // // Skip header of geometry file
  rewind(fp);
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#ENDOFHEADER", text  ) == 0 )  {
      found_eoh = 1;
    }
  }

  // // Count number of surfaces in the spacecraft   
  if ( found_ensemble_cd == 1 ){
    OPTIONS->n_surfaces = -1; // start at -1 because #NB_ENSEMBLES_CD is not a surface
  }
  else{
    OPTIONS->n_surfaces = 0;
  }
  while ( (read = getline(&line, &len, fp)) != -1 ) {
    if (strchr(&line[0],'#') != NULL)
      OPTIONS->n_surfaces = OPTIONS->n_surfaces + 1 ;
  }
  rewind(fp);
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#ENDOFHEADER", text  ) == 0 )  {
      found_eoh = 1;
    }
  }

  /* if (found_ensemble_cd == 1){ */
  /*   getline(&line, &len, fp);getline(&line, &len, fp); */
  /* } */

  // // Go through all surfaces: normal, area, solar cell area, Cd, solar radiation coefficient
  //  while ( (read = getline(&line, &len, fp)) != -1 ) {
  surface_counter = 0;
  while (surface_counter < OPTIONS->n_surfaces){

      found_new_surface = 0; 
      while ( (found_new_surface == 0) && (read = getline(&line, &len, fp)) != -1 ) {
	if (strchr(&line[0],'#') != NULL)
	  found_new_surface = 1;
      }
      sscanf(line, "%s", text);
      if (  strcmp( "#NB_ENSEMBLES_CD", text  ) != 0 ) { // it could correspond to section #NB_ENSEMBLES_CD, in which case it needs to be skipped
	//  printf("<%s>\n", line);
      strcpy(string_nb_word,"%");
      sprintf(str_len, "%lu", strlen(line)-1);
      strcat(string_nb_word, str_len);
      strcat(string_nb_word, "c");
      sscanf( line, string_nb_word, OPTIONS->surface[surface_counter].name_of_surface);
      getline(&line, &len, fp);
      sscanf( line, "(%lf; %lf; %lf)", &temp_normal[0],  &temp_normal[1], &temp_normal[2] );
      v_norm(OPTIONS->surface[surface_counter].normal, temp_normal);// normalize the normal vector in case the user forgot to do so
      getline(&line, &len, fp);
      sscanf( line, "%lf", &OPTIONS->surface[surface_counter].area );
      // is this surface an effective surface? A surface is effective if its normal is in a  direction that is different from that of all the other effective surfaces
      sss = 0;
      is_surf_eff = 1;
      while (sss < surface_counter){
	v_dot(&cos_surf, OPTIONS->surface[surface_counter].normal, OPTIONS->surface[sss].normal);
	if (cos_surf >= cos_limit){ // this means that the directions of the two sirfaces are different by less than acos(cos_limit) degrees, ie that they basically point in the same direction, so surface_counter is not a surface effective
	  is_surf_eff = 0;
	  break;
	}
	sss = sss + 1;
      }
      OPTIONS->surface[surface_counter].ieff = -1;
      if (is_surf_eff == 1){
      	  v_copy(OPTIONS->surface_eff[OPTIONS->n_surfaces_eff].normal, OPTIONS->surface[surface_counter].normal);
	  OPTIONS->surface_eff[OPTIONS->n_surfaces_eff].area = OPTIONS->surface[surface_counter].area;
	  OPTIONS->surface[surface_counter].ieff = OPTIONS->n_surfaces_eff; // ieff is the index of the first surface that has a particular direction, ie the first surface to be effective in this direction
	  OPTIONS->n_surfaces_eff = OPTIONS->n_surfaces_eff + 1;
      }
      else{
	  OPTIONS->surface_eff[OPTIONS->surface[sss].ieff].area = OPTIONS->surface_eff[OPTIONS->surface[sss].ieff].area + OPTIONS->surface[surface_counter].area;
      }
      
      getline(&line, &len, fp);
      sscanf( line, "%lf", &OPTIONS->surface[surface_counter].area_solar_panel );
      getline(&line, &len, fp);
      if ( OPTIONS->nb_ensembles_cd > 0 ){
	sscanf( line, "%lf, %lf", &OPTIONS->surface[surface_counter].Cd, &OPTIONS->surface[surface_counter].Cd_sigma );
      }
      else{


	sscanf( line, "%lf", &cd_or_acco ); 	
	if (cd_or_acco > 1){
	  OPTIONS->new_cd = 0;
	}
	else{
	  OPTIONS->new_cd = 1;
	}
	// new cd
	if (OPTIONS->new_cd == 1){
	sscanf( line, "%lf", &OPTIONS->surface[surface_counter].acco_coeff ); 	
	}

	else{
	sscanf( line, "%lf", &OPTIONS->surface[surface_counter].Cd ); 	
	}
	// end of new cd
      }

      // !!!!!!!!!!! THIS BLOCK IS USING THE EQUATION FROM STK (http://www.agi.com/resources/help/online/stk/10.1/index.html?page=source%2Fhpop%2Fhpop-05.htm). UNCOMMENT THE BLOCK BELOW THAT USES VALLADO AND COMMENT THIS STK BLOCK IF YOU WANT TO USE VALLADO'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)
      getline(&line, &len, fp);
      sscanf( line, "%lf", &OPTIONS->surface[surface_counter].solar_radiation_coefficient ); 
      // !!!!!!!!!!! END OF THIS BLOCK IS USING THE EQUATION FROM STK (http://www.agi.com/resources/help/online/stk/10.1/index.html?page=source%2Fhpop%2Fhpop-05.htm). UNCOMMENT THE BLOCK BELOW THAT USES VALLADO AND COMMENT THIS STK BLOCK IF YOU WANT TO USE VALLADO'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)
      // !!!!!!!!!! THIS BLOCK USES VALLADO'S EQUATIONS. COMMENT IT AND UNCOMMENT THE BLOCK ABOVE THAT USES STK IF YOU WANT TO USE STK'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)
      /* getline(&line, &len, fp); */
      /* sscanf( line, "%lf", &OPTIONS->surface[surface_counter].specular_reflectivity ); */
      /* getline(&line, &len, fp); */
      /* sscanf( line, "%lf", &OPTIONS->surface[surface_counter].diffuse_reflectivity ); */
      // !!!!!!!!!! END OF THIS BLOCK USES VALLADO'S EQUATIONS. COMMENT IT AND UNCOMMENT THE BLOCK ABOVE THAT USES STK IF YOU WANT TO USE STK'S EQUATIONS. ALSO NEED TO CHANGE initialize_constellation.c AND load_options.c TO READ THE SPECULAR AND DIFFUSE REFLECIVITIES IF YOU WANT TO USE VALLADO'S EQUATIONS (SEE COMMENTS IN THESE CODES)
      //  getline(&line, &len, fp); 
      //      printf("\n*** Surface %s ***\nnormal: %f %f %f\narea: %f\nsolar: %f\n Cd: %f %f\nCs: %f\n", OPTIONS->surface[surface_counter].name_of_surface, OPTIONS->surface[surface_counter].normal[0],OPTIONS->surface[surface_counter].normal[1], OPTIONS->surface[surface_counter].normal[2], OPTIONS->surface[surface_counter].area, OPTIONS->surface[surface_counter].area_solar_panel, OPTIONS->surface[surface_counter].Cd, OPTIONS->surface[surface_counter].Cd_sigma, OPTIONS->surface[surface_counter].solar_radiation_coefficient);
      surface_counter = surface_counter + 1;
      }
    }
  //  pti(OPTIONS->n_surfaces, "N");

    //  }

  /* printf("%d %f\n", OPTIONS->n_surfaces_eff, OPTIONS->n_surfaces); */
  /* for (sss = 0; sss < OPTIONS->n_surfaces_eff; sss++){ */
  /*   printf("normal[%d]: (%f, %f, %f) | total area: %f\n", sss, OPTIONS->surface_eff[sss].normal[0], OPTIONS->surface_eff[sss].normal[1],OPTIONS->surface_eff[sss].normal[2], OPTIONS->surface_eff[sss].area); */
  /* } */
  /* exitf(); */
  
  if (OPTIONS->opengl == 1){
    FILE *fp_opengl;
    fp_opengl = fopen(OPTIONS->filename_area_attitude_opengl, "r");
  
  getline(&line, &len, fp_opengl);
  getline(&line, &len, fp_opengl);
  getline(&line, &len, fp_opengl);
  getline(&line, &len, fp_opengl);
  getline(&line, &len, fp_opengl);
  getline(&line, &len, fp_opengl);
  sscanf(line, "%s %d", text, &OPTIONS->nb_faces);
  getline(&line, &len, fp_opengl);
  sscanf(line, "%s %lf", text, &OPTIONS->area_attitude_opengl_dtheta);
  getline(&line, &len, fp_opengl);
  sscanf(line, "%s %lf", text, &OPTIONS->area_attitude_opengl_dphi);
  getline(&line, &len, fp_opengl);
  sscanf(line, "%s %lf", text, &OPTIONS->area_attitude_opengl_theta0);
  getline(&line, &len, fp_opengl);
  sscanf(line, "%s %lf", text, &OPTIONS->area_attitude_opengl_phi0);

  int found_start_area_attitude = 0;
  while ( found_start_area_attitude == 0 && !feof(fp_opengl) ){
    getline(&line, &len, fp_opengl);
    sscanf(line, "%s", text);
    if (  strcmp( "#START", text  ) == 0 )  {
      found_start_area_attitude = 1;
    }
  }
  if (feof(fp_opengl)){
    printf("***! No data found in %s. The program will stop. !***\n", OPTIONS->filename_area_attitude_opengl);
    MPI_Finalize(); exit(0);
  }

  int nb_phi, nb_theta;
  nb_theta = (int)((180 - OPTIONS->area_attitude_opengl_theta0) / OPTIONS->area_attitude_opengl_dtheta ) + 1; // 180 is included (because different from theta = 0)
  nb_phi = (int)((360 - OPTIONS->area_attitude_opengl_phi0) / OPTIONS->area_attitude_opengl_dphi ); /// 360 is not included (because same as phi = 0)
  OPTIONS->area_attitude_opengl = NULL;
  OPTIONS->area_attitude_opengl = malloc(nb_theta * sizeof(double **));
  OPTIONS->nb_faces_theta_phi = malloc(nb_theta * sizeof(double *));
  if ( OPTIONS->area_attitude_opengl == NULL){
    printf("***! Could not allow memory to OPTIONS->area_attitude_opengl. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }


  OPTIONS->area_attitude_opengl_total = NULL;
  OPTIONS->area_attitude_opengl_total = malloc(nb_theta * sizeof(double *));
  if ( OPTIONS->area_attitude_opengl_total == NULL){
    printf("***! Could not allow memory to OPTIONS->area_attitude_opengl_total. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }


  OPTIONS->which_face_theta_phi = NULL;
  OPTIONS->which_face_theta_phi = malloc(nb_theta * sizeof(int **));
  if ( OPTIONS->which_face_theta_phi == NULL){
    printf("***! Could not allow memory to OPTIONS->which_face_theta_phi. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }


  int itheta, iphi;
  double theta_dummy, phi_dummy;
  int iface, iface_here, iface_dummy;
  
  OPTIONS->normal_face = NULL;
  OPTIONS->normal_face = malloc(OPTIONS->nb_faces * sizeof(double *));
    for (iface = 0; iface < OPTIONS->nb_faces; iface++){
      OPTIONS->normal_face[iface] = malloc(3 * sizeof(double));
      getline(&line, &len, fp_opengl);
      sscanf(line, "%d %lf %lf %lf", &iface_dummy, &OPTIONS->normal_face[iface][0], &OPTIONS->normal_face[iface][1], &OPTIONS->normal_face[iface][2]);

    }
    double dummy;
  getline(&line, &len, fp_opengl); // skip line '#START SECOND BLOCK'
  for (itheta = 0; itheta < nb_theta; itheta++){
    OPTIONS->area_attitude_opengl[itheta] = malloc(nb_phi * sizeof(double *));
    OPTIONS->area_attitude_opengl_total[itheta] = malloc(nb_phi * sizeof(double));
    OPTIONS->which_face_theta_phi[itheta] = malloc(nb_phi * sizeof(int *));
    OPTIONS->nb_faces_theta_phi[itheta] = malloc(nb_phi * sizeof(double));
  if ( OPTIONS->area_attitude_opengl[itheta] == NULL){
    printf("***! Could not allow memory to OPTIONS->area_attitude_opengl[itheta]. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }
  if ( OPTIONS->area_attitude_opengl_total[itheta] == NULL){
    printf("***! Could not allow memory to OPTIONS->area_attitude_opengl_total[itheta]. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }
  if ( OPTIONS->which_face_theta_phi[itheta] == NULL){
    printf("***! Could not allow memory to OPTIONS->which_face_theta_phi[itheta]. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }

    for (iphi = 0; iphi < nb_phi; iphi++){
      OPTIONS->area_attitude_opengl[itheta][iphi] = malloc(OPTIONS->nb_faces * sizeof(double)); // allocate nb_faces even if go onlyup to nb_faces_theta_phi. this is so that whaever theta and phi, the number of faces allocated is the same
      OPTIONS->which_face_theta_phi[itheta][iphi] = malloc(OPTIONS->nb_faces * sizeof(int)); 
      getline(&line, &len, fp_opengl);
      sscanf(line, "# %lf %lf %d", &theta_dummy, &phi_dummy, &OPTIONS->nb_faces_theta_phi[itheta][iphi]);
      for (iface = 0; iface < OPTIONS->nb_faces_theta_phi[itheta][iphi]; iface++){
      getline(&line, &len, fp_opengl);
      sscanf(line, "%d %lf", &iface_here, &dummy);
 OPTIONS->area_attitude_opengl[itheta][iphi][iface_here] = dummy; // OPTIONS->area_attitude_opengl in cm^2. Will be converted to km^2 in initialize_constellation.c
      OPTIONS->which_face_theta_phi[itheta][iphi][iface] = iface_here;
      //printf("which face %d\n", OPTIONS->which_face_theta_phi[itheta][iphi][iface]);
      }
      //      exitf();
      getline(&line, &len, fp_opengl);
      sscanf(line, "%lf", &OPTIONS->area_attitude_opengl_total[itheta][iphi]); // total area projected

    }
  }
  //  exitf();
/*   for (itheta = 0; itheta < nb_theta; itheta++){ */
/*     for (iphi = 0; iphi < nb_phi; iphi++){ */
/*       printf("# %f %f %d\n", itheta * OPTIONS->area_attitude_opengl_dtheta + OPTIONS->area_attitude_opengl_theta0, iphi * OPTIONS->area_attitude_opengl_dphi + OPTIONS->area_attitude_opengl_phi0, OPTIONS->nb_faces_theta_phi[itheta][iphi]); */
/*       for (iface = 0; iface < OPTIONS->nb_faces_theta_phi[itheta][iphi]; iface++){ */
/* 	iface_here = OPTIONS->which_face_theta_phi[itheta][iphi][iface]; */
/* 	  printf("%d %f\n", iface_here, OPTIONS->area_attitude_opengl[itheta][iphi][iface_here]); */
/*       } */
/*       printf("%f\n", OPTIONS->area_attitude_opengl_total[itheta][iphi]); */
/*     } */
/*   } */
    
  fclose(fp_opengl);









  // NOW READ FILE WITH CORNERS OF SOLAR PANELS (IF POWER IS COMPUTED)
  if (( OPTIONS->solar_cell_efficiency != -1) && (OPTIONS->opengl_power == 1)){

    FILE *fp_opengl;
    fp_opengl = fopen(OPTIONS->opengl_filename_solar_power, "r");
  
  int found_start_area_attitude = 0;
  //go to block #START FIRST BLOCK
  while ( found_start_area_attitude == 0 && !feof(fp_opengl) ){
    getline(&line, &len, fp_opengl);
    sscanf(line, "%s", text);
    if (  strcmp( "#START", text  ) == 0 )  {
      found_start_area_attitude = 1;
    }
  }
  if (feof(fp_opengl)){
    printf("***! No data found in %s. The program will stop. !***\n", OPTIONS->filename_area_attitude_opengl);
    MPI_Finalize(); exit(0);
  }
    getline(&line, &len, fp_opengl);
    found_start_area_attitude = 0;
  // don't need to read normals (already read, exccpet those of solar panels but we don't care about them) so go to block #START SECOND BLOCK
  while ( found_start_area_attitude == 0 && !feof(fp_opengl) ){
    getline(&line, &len, fp_opengl);
    sscanf(line, "%s", text);
    if (  strcmp( "#START", text  ) == 0 )  {
      found_start_area_attitude = 1;
    }
  }
  if (feof(fp_opengl)){
    printf("***! No data found in %s. The program will stop. !***\n", OPTIONS->filename_area_attitude_opengl);
    MPI_Finalize(); exit(0);
  }


  OPTIONS->area_solar_panel_attitude_opengl = NULL;
  OPTIONS->area_solar_panel_attitude_opengl = malloc(nb_theta * sizeof(double *));
  if ( OPTIONS->area_solar_panel_attitude_opengl == NULL){
    printf("***! Could not allow memory to OPTIONS->area_solar_panel_attitude_opengl. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }


    double dummy;
  for (itheta = 0; itheta < nb_theta; itheta++){
    OPTIONS->area_solar_panel_attitude_opengl[itheta] = malloc(nb_phi * sizeof(double ));

  if ( OPTIONS->area_solar_panel_attitude_opengl[itheta] == NULL){
    printf("***! Could not allow memory to OPTIONS->area_solar_panel_attitude_opengl[itheta]. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  }

    for (iphi = 0; iphi < nb_phi; iphi++){
      getline(&line, &len, fp_opengl);
      //printf("<%s>\n", line);
      sscanf(line, "# %lf %lf %lf", &theta_dummy, &phi_dummy, &OPTIONS->area_solar_panel_attitude_opengl[itheta][iphi]);

    }
  }
  //  exitf();
/*   for (itheta = 0; itheta < nb_theta; itheta++){ */
/*     for (iphi = 0; iphi < nb_phi; iphi++){ */
/*             printf("# %f %f %f\n", itheta * OPTIONS->area_attitude_opengl_dtheta + OPTIONS->area_attitude_opengl_theta0, iphi * OPTIONS->area_attitude_opengl_dphi + OPTIONS->area_attitude_opengl_phi0, OPTIONS->area_solar_panel_attitude_opengl[itheta][iphi]); */
/*     } */
/*   } */
    
  fclose(fp_opengl);

  }




  }

  //  MPI_Finalize();exit(0);
  return 0;

}



/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           load_attitude
//  Purpose:        Load the attitude file (nadir, sun_pointed or manual (= name of external file))
//  Assumptions:    Same attitude configuration for each satellite of the constellation
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 08/12/2015    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////


int load_attitude( OPTIONS_T *OPTIONS,
		   char      attitude_profile[256], FILE *input_file, int nProcs, double ang_velo[3], int iDebugLevel, int iProc){

  /* Declarations */
  int j;
  int ierr;
  int iia;
  char *line = NULL;
  size_t len = 0;
  int found_eoh = 0;
  char text_location[256];
  char temp_copy[256];
  char text[256];
  double et_initial;

    str2et_c(OPTIONS->initial_epoch, &et_initial);
  //char times[256];                    
  /* Algorithm */
  /* nadir */
  if ((iProc == 0) & ( iDebugLevel >= 3 ) ){
    printf("---- (load_options) (load_attitude) Beginning of function load_attitude.\n");
  }

  if (strcmp(attitude_profile, "nadir") == 0){ 
    if (OPTIONS->use_kalman != 1){
    OPTIONS->et_interpo = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    }
    OPTIONS->pitch = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->roll = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->yaw = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->order_pitch = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->order_roll = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->order_yaw = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );

    if (  OPTIONS->pitch == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->pitch \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  OPTIONS->roll == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->roll \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  OPTIONS->yaw == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->yaw \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }

    if (  OPTIONS->order_pitch == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->order_pitch \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  OPTIONS->order_roll == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->order_roll \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  OPTIONS->order_yaw == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->order_yaw \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }

    if (OPTIONS->use_kalman != 1){
    if ( OPTIONS->et_interpo == NULL ){
      printf("***! Could not allow memory space for OPTIONS->et_interpo  \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    }


  iia = 0;
    if (OPTIONS->use_kalman != 1){
  OPTIONS->et_interpo[iia] = OPTIONS->et_oldest_tle_epoch;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    }
      OPTIONS->pitch[iia] = 0;
      OPTIONS->roll[iia] = 0;
      OPTIONS->yaw[iia] = 0;
      OPTIONS->order_pitch[iia] = 1;
      OPTIONS->order_roll[iia] = 2;
      OPTIONS->order_yaw[iia] = 3;

      while (OPTIONS->et_interpo[iia] + OPTIONS->dt +OPTIONS->dt/1000. <= et_initial){ // +OPTIONS->dt/1000.  for numerical reasons//  while (OPTIONS->et_interpo[iia] < et_initial){
    iia = iia + 1;
    if (OPTIONS->use_kalman != 1){
    OPTIONS->et_interpo[iia] = OPTIONS->et_oldest_tle_epoch + OPTIONS->dt*iia / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    }
    //        etprint(OPTIONS->et_interpo[iia], "tle");
      OPTIONS->pitch[iia] = 0;
      OPTIONS->roll[iia] = 0;
      OPTIONS->yaw[iia] = 0;
      OPTIONS->order_pitch[iia] = 1;
      OPTIONS->order_roll[iia] = 2;
      OPTIONS->order_yaw[iia] = 3;
	
  } // leave this loop when x_after_interpo gets newer than constellation epoch (et_initial)

      //if (iia == 0){// in that case, we didnt get in the previous while loop  but we still need ot increment iia by 1
	iia = iia + 1;

	//}

  // one more time to get to et_initial
      if (fabs(et_initial- OPTIONS->et_oldest_tle_epoch) > OPTIONS->dt/1000.){ // +OPTIONS->dt/1000. for numerical reasons

  OPTIONS->et_interpo[iia] = OPTIONS->et_interpo[iia-1] + (et_initial - OPTIONS->et_interpo[iia-1])/2.;
      OPTIONS->pitch[iia] = 0;
      OPTIONS->roll[iia] = 0;
      OPTIONS->yaw[iia] = 0;
      OPTIONS->order_pitch[iia] = 1;
      OPTIONS->order_roll[iia] = 2;
      OPTIONS->order_yaw[iia] = 3;
	

      //      printf("i = %d\n", iia);
      //  etprint(OPTIONS->et_interpo[iia], "one more");
  iia = iia + 1;
  }
      
      
  if (fabs(et_initial- OPTIONS->et_oldest_tle_epoch) > OPTIONS->dt/1000.){
  j = 0;
  }
  else{ // if et_initial == et_oldest_tle_epoch then don't enter twice et_initial in et_interpo
    j = 1;
  }
  while (iia<OPTIONS->nb_time_steps * 2){
    if (OPTIONS->use_kalman != 1){
    OPTIONS->et_interpo[iia] = et_initial + OPTIONS->dt*j / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    }
    //        etprint(OPTIONS->et_interpo[iia], "const");
      OPTIONS->pitch[iia] = 0;
      OPTIONS->roll[iia] = 0;
      OPTIONS->yaw[iia] = 0;
      OPTIONS->order_pitch[iia] = 1;
      OPTIONS->order_roll[iia] = 2;
      OPTIONS->order_yaw[iia] = 3;

	iia = iia +1;
	j = j+1;
  }


/*     for (iia = 0; iia < OPTIONS->nb_time_steps * 2; iia++){ */
/*       OPTIONS->et_interpo[iia] = et_initial + OPTIONS->dt*iia / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*       OPTIONS->pitch[iia] = 0; */
/*       OPTIONS->roll[iia] = 0; */
/*       OPTIONS->yaw[iia] = 0; */
/*       OPTIONS->order_pitch[iia] = 1; */
/*       OPTIONS->order_roll[iia] = 2; */
/*       OPTIONS->order_yaw[iia] = 3; */
/*     } */

  }

  /* sun_pointed */
  else if (strcmp(attitude_profile, "sun_pointed") == 0){
 
    if (OPTIONS->use_kalman != 1){
    OPTIONS->et_interpo = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    }
    OPTIONS->pitch = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->roll = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->yaw = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->order_pitch = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->order_roll = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->order_yaw = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );

    if (  OPTIONS->pitch == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->pitch \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  OPTIONS->roll == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->roll \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  OPTIONS->yaw == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->yaw \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }

    if (  OPTIONS->order_pitch == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->order_pitch \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  OPTIONS->order_roll == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->order_roll \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  OPTIONS->order_yaw == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->order_yaw \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }

    if (OPTIONS->use_kalman != 1){
    if ( OPTIONS->et_interpo == NULL ){
      printf("***! Could not allow memory space for OPTIONS->et_interpo  \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    }

  iia = 0;
    if (OPTIONS->use_kalman != 1){
  OPTIONS->et_interpo[iia] = OPTIONS->et_oldest_tle_epoch;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    }
      OPTIONS->pitch[iia] = 0;
      OPTIONS->roll[iia] = 0;
      OPTIONS->yaw[iia] = 0;
      OPTIONS->order_pitch[iia] = 1;
      OPTIONS->order_roll[iia] = 2;
      OPTIONS->order_yaw[iia] = 3;

  while (OPTIONS->et_interpo[iia] < et_initial){
    iia = iia + 1;
    if (OPTIONS->use_kalman != 1){
    OPTIONS->et_interpo[iia] = OPTIONS->et_oldest_tle_epoch + OPTIONS->dt*iia / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    }
    //        etprint(OPTIONS->et_interpo[iia], "tle");
      OPTIONS->pitch[iia] = 0;
      OPTIONS->roll[iia] = 0;
      OPTIONS->yaw[iia] = 0;
      OPTIONS->order_pitch[iia] = 1;
      OPTIONS->order_roll[iia] = 2;
      OPTIONS->order_yaw[iia] = 3;
	
  } // leave this loop when x_after_interpo gets newer than constellation epoch (et_initial)
  j = 0;
  while (iia<OPTIONS->nb_time_steps * 2){
    if (OPTIONS->use_kalman != 1){
    OPTIONS->et_interpo[iia] = et_initial + OPTIONS->dt*j / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    }
    //        etprint(OPTIONS->et_interpo[iia], "const");
      OPTIONS->pitch[iia] = 0;
      OPTIONS->roll[iia] = 0;
      OPTIONS->yaw[iia] = 0;
      OPTIONS->order_pitch[iia] = 1;
      OPTIONS->order_roll[iia] = 2;
      OPTIONS->order_yaw[iia] = 3;

	iia = iia +1;
	j = j+1;
  }

      
/*     for (iia = 0; iia < OPTIONS->nb_time_steps * 2; iia++){ */
/*       OPTIONS->et_interpo[iia] = et_initial + OPTIONS->dt*iia / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*       OPTIONS->pitch[iia] = 0; */
/*       OPTIONS->roll[iia] = 0; */
/*       OPTIONS->yaw[iia] = 0; */
/*       OPTIONS->order_pitch[iia] = 1; */
/*       OPTIONS->order_roll[iia] = 2; */
/*       OPTIONS->order_yaw[iia] = 3; */
/*     } */
  }

  else if (strcmp(attitude_profile, "ensemble_angular_velocity") == 0){
    if (OPTIONS->use_kalman != 1){
    OPTIONS->et_interpo = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    }
    rewind(input_file);
    found_eoh = 0;
    while ( found_eoh == 0 && !feof(input_file)) {
      getline(&line, &len, input_file);
      sscanf(line, "%s", text);
      if (  strcmp( "##ENSEMBLES_ATTITUDE", text  ) == 0 )  {
	found_eoh = 1;
      }
    }
    if (feof(input_file)){
      printf("***! You chose ot run ensembles on the attitude ('ensemble_angular_velocity' in section #ATTITUDE) but there is no section ##ENSEMBLES_ATTITUDE in the main input file. The program will stop. !***\n"); MPI_Finalize(); exit(0);
    }
    getline(&line, &len, input_file);
    sscanf(line,"%d", &OPTIONS->nb_ensembles_attitude);
    if (nProcs > 0){
      OPTIONS->nb_ensemble_attitude_per_proc = (int)(OPTIONS->nb_ensembles_attitude / nProcs);
    }

    rewind(input_file);
    found_eoh = 0;
    while ( found_eoh == 0 && !feof(input_file)) {
      getline(&line, &len, input_file);
      sscanf(line, "%s", text);
      if (  strcmp( "###ENSEMBLES_ANGULAR_VELOCITY_DEPLOYMENT", text  ) == 0 )  {
	found_eoh = 1;
      }
    }
    if (feof(input_file)){
      printf("***! You chose ot run ensembles on the attitude ('ensemble_angular_velocity' in section #ATTITUDE) but there is no section ###ENSEMBLES_ANGULAR_VELOCITY_DEPLOYMENT in the main input file. The program will stop. !***\n"); MPI_Finalize(); exit(0);
    }

    getline(&line, &len, input_file);
    sscanf(line, "(%lf; %lf; %lf)", &OPTIONS->pitch_ini_ensemble, &OPTIONS->roll_ini_ensemble, &OPTIONS->yaw_ini_ensemble);
    getline(&line, &len, input_file);
    sscanf(line, "(%lf; %lf; %lf)", &OPTIONS->pitch_mean_angular_velocity_ensemble, &OPTIONS->roll_mean_angular_velocity_ensemble, &OPTIONS->yaw_mean_angular_velocity_ensemble);
    getline(&line, &len, input_file);
    sscanf(line, "(%lf; %lf; %lf)", &OPTIONS->pitch_sigma_angular_velocity_ensemble, &OPTIONS->roll_sigma_angular_velocity_ensemble, &OPTIONS->yaw_sigma_angular_velocity_ensemble);


  iia = 0;
    if (OPTIONS->use_kalman != 1){
  OPTIONS->et_interpo[iia] = OPTIONS->et_oldest_tle_epoch;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    }
  while (OPTIONS->et_interpo[iia] < et_initial){
    iia = iia + 1;
    if (OPTIONS->use_kalman != 1){
    OPTIONS->et_interpo[iia] = OPTIONS->et_oldest_tle_epoch + OPTIONS->dt*iia / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    }
    //        etprint(OPTIONS->et_interpo[iia], "tle");
	
  } // leave this loop when x_after_interpo gets newer than constellation epoch (et_initial)
  j = 0;
  while (iia<OPTIONS->nb_time_steps * 2){
    if (OPTIONS->use_kalman != 1){
    OPTIONS->et_interpo[iia] = et_initial + OPTIONS->dt*j / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    }
    //        etprint(OPTIONS->et_interpo[iia], "const");
	iia = iia +1;
	j = j+1;
  }

/*     for (iia = 0; iia < OPTIONS->nb_time_steps * 2; iia++){ */
/*       OPTIONS->et_interpo[iia] = et_initial + OPTIONS->dt*iia / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*     } */
  }

  else if (strcmp(attitude_profile, "ensemble_initial_attitude") == 0){
    if (OPTIONS->use_kalman != 1){
    OPTIONS->et_interpo = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    }
    rewind(input_file);
    found_eoh = 0;
    while ( found_eoh == 0 && !feof(input_file)) {
      getline(&line, &len, input_file);
      sscanf(line, "%s", text);
      if (  strcmp( "##ENSEMBLES_ATTITUDE", text  ) == 0 )  {
	found_eoh = 1;
      }
    }
    if (feof(input_file)){
      printf("***! You chose ot run ensembles on the attitude ('ensemble_initial_attitude' in section #ATTITUDE) but there is no section ##ENSEMBLES_ATTITUDE in the main input file. The program will stop. !***\n"); MPI_Finalize(); exit(0);
    }

    getline(&line, &len, input_file);
    sscanf(line,"%d", &OPTIONS->nb_ensembles_attitude);
    if (nProcs > 0){
      OPTIONS->nb_ensemble_attitude_per_proc = (int)(OPTIONS->nb_ensembles_attitude / nProcs);
    }
    rewind(input_file);
    found_eoh = 0;
    while ( found_eoh == 0 && !feof(input_file)) {
      getline(&line, &len, input_file);
      sscanf(line, "%s", text);
      if (  strcmp( "###ENSEMBLES_INITIAL_ATTITUDE", text  ) == 0 )  {
	found_eoh = 1;
      }
    }
    if (feof(input_file)){
      printf("***! You chose ot run ensembles on the attitude ('ensemble_initial_attitude' in section #ATTITUDE) but there is no section ###ENSEMBLES_INITIAL_ATTITUDE in the main input file. The program will stop. !***\n"); MPI_Finalize(); exit(0);
    }

    getline(&line, &len, input_file);
    sscanf(line, "(%lf; %lf; %lf)", &OPTIONS->pitch_mean_ensemble, &OPTIONS->roll_mean_ensemble, &OPTIONS->yaw_mean_ensemble);
    getline(&line, &len, input_file);
    sscanf(line, "(%lf; %lf; %lf)", &OPTIONS->pitch_sigma_for_ensemble_initial_attitude, &OPTIONS->roll_sigma_for_ensemble_initial_attitude, &OPTIONS->yaw_sigma_for_ensemble_initial_attitude);
    getline(&line, &len, input_file);
    sscanf(line, "(%lf; %lf; %lf)", &OPTIONS->pitch_angular_velocity_constant, &OPTIONS->roll_angular_velocity_constant, &OPTIONS->yaw_angular_velocity_constant);

  iia = 0;
    if (OPTIONS->use_kalman != 1){
  OPTIONS->et_interpo[iia] = OPTIONS->et_oldest_tle_epoch;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    }
  while (OPTIONS->et_interpo[iia] < et_initial){
    iia = iia + 1;
    if (OPTIONS->use_kalman != 1){
    OPTIONS->et_interpo[iia] = OPTIONS->et_oldest_tle_epoch + OPTIONS->dt*iia / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    }
    //        etprint(OPTIONS->et_interpo[iia], "tle");
	
  } // leave this loop when x_after_interpo gets newer than constellation epoch (et_initial)
  j = 0;
  while (iia<OPTIONS->nb_time_steps * 2){
    if (OPTIONS->use_kalman != 1){
    OPTIONS->et_interpo[iia] = et_initial + OPTIONS->dt*j / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    }
    //        etprint(OPTIONS->et_interpo[iia], "const");
	iia = iia +1;
	j = j+1;
  }

/*     for (iia = 0; iia < OPTIONS->nb_time_steps * 2; iia++){ */
/*       OPTIONS->et_interpo[iia] = et_initial + OPTIONS->dt*iia / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*     } */

  }

  else if (strcmp(attitude_profile, "angular_velocity") == 0){ 
    if (OPTIONS->use_kalman != 1){
      OPTIONS->et_interpo = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    }
    OPTIONS->pitch = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->roll = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->yaw = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->order_pitch = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->order_roll = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->order_yaw = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );

    if (  OPTIONS->pitch == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->pitch \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  OPTIONS->roll == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->roll \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  OPTIONS->yaw == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->yaw \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }

    if (  OPTIONS->order_pitch == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->order_pitch \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  OPTIONS->order_roll == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->order_roll \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  OPTIONS->order_yaw == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->order_yaw \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }

    if (OPTIONS->use_kalman != 1){
    if ( OPTIONS->et_interpo == NULL ){
      printf("***! Could not allow memory space for OPTIONS->et_interpo  \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    }

  iia = 0;
    if (OPTIONS->use_kalman != 1){
  OPTIONS->et_interpo[iia] = OPTIONS->et_oldest_tle_epoch;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    }
      OPTIONS->pitch[iia] = iia * ang_velo[3] * OPTIONS->dt / 2. + ang_velo[0];
      OPTIONS->roll[iia] = iia * ang_velo[4] * OPTIONS->dt / 2. + ang_velo[1];
      OPTIONS->yaw[iia] = iia * ang_velo[5] * OPTIONS->dt / 2. + ang_velo[2] ;
/*       OPTIONS->pitch[iia] = iia * ang_velo[0] * OPTIONS->dt / 2.; */
/*       OPTIONS->roll[iia] = iia * ang_velo[1] * OPTIONS->dt / 2.; */
/*       OPTIONS->yaw[iia] = iia * ang_velo[2] * OPTIONS->dt / 2.; */
      OPTIONS->order_pitch[iia] = 1;
      OPTIONS->order_roll[iia] = 2;
      OPTIONS->order_yaw[iia] = 3;

  while (OPTIONS->et_interpo[iia] < et_initial){
    iia = iia + 1;
    if (OPTIONS->use_kalman != 1){
    OPTIONS->et_interpo[iia] = OPTIONS->et_oldest_tle_epoch + OPTIONS->dt*iia / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    }
      OPTIONS->pitch[iia] = iia * ang_velo[3] * OPTIONS->dt / 2. + ang_velo[0];
      OPTIONS->roll[iia] = iia * ang_velo[4] * OPTIONS->dt / 2. + ang_velo[1];
      OPTIONS->yaw[iia] = iia * ang_velo[5] * OPTIONS->dt / 2. + ang_velo[2] ;
/*       OPTIONS->pitch[iia] = iia * ang_velo[0] * OPTIONS->dt / 2.; */
/*       OPTIONS->roll[iia] = iia * ang_velo[1] * OPTIONS->dt / 2.; */
/*       OPTIONS->yaw[iia] = iia * ang_velo[2] * OPTIONS->dt / 2.; */
      OPTIONS->order_pitch[iia] = 1;
      OPTIONS->order_roll[iia] = 2;
      OPTIONS->order_yaw[iia] = 3;

    //        etprint(OPTIONS->et_interpo[iia], "tle");
	
  } // leave this loop when x_after_interpo gets newer than constellation epoch (et_initial)
  j = 0;
  while (iia<OPTIONS->nb_time_steps * 2){
    if (OPTIONS->use_kalman != 1){
    OPTIONS->et_interpo[iia] = et_initial + OPTIONS->dt*j / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    }
      OPTIONS->pitch[iia] = iia * ang_velo[3] * OPTIONS->dt / 2. + ang_velo[0];
      OPTIONS->roll[iia] = iia * ang_velo[4] * OPTIONS->dt / 2. + ang_velo[1];
      OPTIONS->yaw[iia] = iia * ang_velo[5] * OPTIONS->dt / 2. + ang_velo[2] ;
/*       OPTIONS->pitch[iia] = iia * ang_velo[0] * OPTIONS->dt / 2.; */
/*       OPTIONS->roll[iia] = iia * ang_velo[1] * OPTIONS->dt / 2.; */
/*       OPTIONS->yaw[iia] = iia * ang_velo[2] * OPTIONS->dt / 2.; */
      OPTIONS->order_pitch[iia] = 1;
      OPTIONS->order_roll[iia] = 2;
      OPTIONS->order_yaw[iia] = 3;

    //        etprint(OPTIONS->et_interpo[iia], "const");
	iia = iia +1;
	j = j+1;
  }


/*     for (iia = 0; iia < OPTIONS->nb_time_steps * 2; iia++){ */
/*       OPTIONS->et_interpo[iia] = et_initial + OPTIONS->dt*iia / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*       OPTIONS->pitch[iia] = iia * ang_velo[3] * OPTIONS->dt / 2. + ang_velo[0]; */
/*       OPTIONS->roll[iia] = iia * ang_velo[4] * OPTIONS->dt / 2. + ang_velo[1]; */
/*       OPTIONS->yaw[iia] = iia * ang_velo[5] * OPTIONS->dt / 2. + ang_velo[2] ; */
/* /\*       OPTIONS->pitch[iia] = iia * ang_velo[0] * OPTIONS->dt / 2.; *\/ */
/* /\*       OPTIONS->roll[iia] = iia * ang_velo[1] * OPTIONS->dt / 2.; *\/ */
/* /\*       OPTIONS->yaw[iia] = iia * ang_velo[2] * OPTIONS->dt / 2.; *\/ */
/*       OPTIONS->order_pitch[iia] = 1; */
/*       OPTIONS->order_roll[iia] = 2; */
/*       OPTIONS->order_yaw[iia] = 3; */
/*     } */

  }

  /* if the attitude representation is in angle (pitch/roll) then convert the representation into the cartesian representation */
  else{ 


    strcpy(temp_copy, OPTIONS->attitude_profile);
    //newstructure
/*     strcpy(text_location, OPTIONS->dir_input_attitude); */
/*     strcat(text_location, "/"); */
    strcpy(text_location, "");
    //newstructure
    strcpy(OPTIONS->attitude_profile, text_location);
    strcat(OPTIONS->attitude_profile,temp_copy);
    if (OPTIONS->use_kalman != 1){
    OPTIONS->et_interpo = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    }
    OPTIONS->pitch = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->roll = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->yaw = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->order_pitch = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->order_roll = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->order_yaw = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double) );
    OPTIONS->quaternion = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double*) );
    for (iia = 0; iia < OPTIONS->nb_time_steps * 2; iia++){
    OPTIONS->quaternion[iia] = malloc( 4 * sizeof(double) );
    }
    if (  OPTIONS->pitch == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->pitch \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  OPTIONS->roll == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->roll \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  OPTIONS->yaw == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->yaw \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }

    if (  OPTIONS->order_pitch == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->order_pitch \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  OPTIONS->order_roll == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->order_roll \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  OPTIONS->order_yaw == NULL ){
      printf("***! Could not allow memory space for  OPTIONS->order_yaw \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }

    if (OPTIONS->use_kalman != 1){
    if ( OPTIONS->et_interpo == NULL ){
      printf("***! Could not allow memory space for OPTIONS->et_interpo  \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    }
    lin_interpolate_attitude(OPTIONS->quaternion,OPTIONS->pitch, OPTIONS->roll, OPTIONS->yaw, OPTIONS->order_pitch, OPTIONS->order_roll, OPTIONS->order_yaw, OPTIONS->et_interpo, OPTIONS->attitude_profile, OPTIONS->nb_time_steps * 2, OPTIONS->initial_epoch, OPTIONS->et_oldest_tle_epoch, OPTIONS->dt, iProc, &OPTIONS->file_is_quaternion, OPTIONS->use_kalman);



  }

  if ((iProc == 0) & ( iDebugLevel >= 3 ) ){
    printf("---- (load_options) (load_attitude) End of function load_attitude.\n");
  }


  return 0;

}



/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           lin_interpolate
//  Purpose:        Linear interpolation of the drivers file name
//  Assumptions:    - The F10.7 and Ap files must have the same time step: value of F10.7 (or Ap) every hour, or every day, or every month
//                  - If this time step is 1 day or a month then NRLMSIS00e uses daily Ap 
//                  - If the time step is one hour, then NRLMSIS00e uses historical Ap (please refer to the NRLMSIS00e documentation for further details). In this case, the function lin_interpolate_ap_hist is used. Please refer to the assumptions of this function.
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 08/14/2015    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////

int lin_interpolate(double *f107_after_interpo,
		    double *f107A_after_interpo,
		    double *ap_after_interpo,
		    double **ap_hist_after_interpo,
		    double *x_after_interpo,	  
		    double *use_ap_hist,
		    char f107_filename[256],
		    char ap_filename[256],
		    char src_file[256],
		    int nb_time_steps_simu,
		    char initial_epoch[256],
		    double et_oldest_tle_epoch, 
		    char final_epoch[256], 
		    double dt,
		    double missing_data_value,
		    int iDebugLevel, int iProc){



  /* Declarations */
  int last_date_input_file_older_final_epoch_ap = 0;
  int ggg;
  //  char times[256];
  int  last_date_input_file_older_final_epoch = 0;
  double sum_81_days = 0 ;
  int ddd;
  double time_step_before_interpo;
  double lat_sat, long_sat, alt_sat, loc_time,dens_norm400, dens_norm410, dens_uncert, ss;
  int nb_data_points;
  int x_min_index, save_first_index_in_f107_before_interpo;
  double *f107_before_interpo = NULL;
  double *f107A_before_interpo = NULL;
  double *ap_before_interpo = NULL;
  double *raid3_before_interpo = NULL;
  double *x_f107_before_interpo = NULL;
  double *x_f107A_before_interpo = NULL;
  double *x_ap_before_interpo = NULL;
  double *x_raid3_before_interpo = NULL;
  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  int found_eoh = 0;
  char text[256];
  int line_num, index_to_stop_in_before_interpo=0, index_to_stop_in_before_interpo_for_f107_average =0, index_to_stop_in_before_interpo_ap = 0;
  char yy[256], doy[256], hh[256];
  double driver_temp;
  double x_f107_before_interpo_temp,  x_ap_before_interpo_temp,  x_raid3_before_interpo_temp;
  int i;
  int nb_elements_in_file_f107=0, nb_elements_in_file_ap=0, nb_elements_in_file_raid3=0;
  double et_initial;  double et_final;
  double *a = NULL;
  double *b = NULL;
  double x_min, y_min, x_max, y_max;


  if ((iProc == 0) & ( iDebugLevel >= 3 ) ){
    printf("---- (load_options) (lin_interpolate) Starting to linear interpolate F10.7 and Ap.\n");
  }

  /* Algorithm */
  /* Calculates the x array on which the interpolation is done */

  str2et_c(initial_epoch, &et_initial);  
  str2et_c(final_epoch, &et_final);  

  ///////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////// F10.7 //////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
    printf("----- (load_options) (lin_interpolate) Starting to linear interpolate F10.7.\n");
  }

  fp = fopen(f107_filename, "r");
  if (fp == NULL){
    printf("***! (load_options) (lin_interpolate) The file %s could not be open. The program will stop. ***!\n", f107_filename); MPI_Finalize(); exit(0);
  }

  if ( ( strcmp(src_file, "omniweb") == 0 ) || ( strcmp(src_file, "dynamic_manual") == 0 ) ){

    /* Calculates the x and y arrays to interpolate */
    if ( strcmp(src_file, "dynamic_manual") == 0 ){
      nb_elements_file(&nb_elements_in_file_f107, f107_filename, "YEAR","#ENDOFFILE");
    }
    if ( strcmp(src_file, "omniweb") == 0 ){
      nb_elements_file(&nb_elements_in_file_f107, f107_filename, "YEAR","</pre><hr><HR>");
    }
    x_f107_before_interpo = malloc( nb_elements_in_file_f107 * sizeof(double) );
    f107_before_interpo = malloc( nb_elements_in_file_f107 * sizeof(double) );

    found_eoh = 0;
    while ( found_eoh == 0 && !feof(fp)) {
      getline(&line, &len, fp);
      sscanf(line, "%s", text);
      if (  strcmp( "YEAR", text  ) == 0 )  {
	found_eoh = 1;
      }
    }


    for (line_num = 0; line_num<nb_elements_in_file_f107; line_num++){
      getline(&line, &len, fp);

      sscanf(line, "%s %s %s %lf", yy, doy, hh,&driver_temp);
      strcpy(text,yy);
      strcat(text,"-");
      strcat(text,doy);
      strcat(text,"/");
      strcat(text,hh);
      strcat(text,":");
      strcat(text,"00");
      str2et_c (text, &x_f107_before_interpo_temp );

      x_f107_before_interpo[line_num] = x_f107_before_interpo_temp;
      f107_before_interpo[line_num] = driver_temp;
      if ( line_num > 0 ){
	if ( f107_before_interpo[line_num] > ( missing_data_value - 0.01 ) ) 
	  f107_before_interpo[line_num] = f107_before_interpo[line_num-1];
      }
      if ( x_f107_before_interpo[line_num] > et_final ){
	index_to_stop_in_before_interpo = line_num;
      }
 
      if ( x_f107_before_interpo[line_num] > et_final + 40.5 * 24 * 3600){ // to further calculate F10.7 81 day average, we need to read F10.7 file up to 40.5 days after the end epoch
	index_to_stop_in_before_interpo_for_f107_average = line_num;
	break;
      }

    }

    	index_to_stop_in_before_interpo_for_f107_average = line_num;
	//etprint(x_f107_before_interpo[line_num-1], "x_f107_before_interpo[line_num]");
    time_step_before_interpo = x_f107_before_interpo[1] - x_f107_before_interpo[0];

    if ( time_step_before_interpo > 23 * 3600 ){ // recall the assumption: the time step of the F10.7/Ap file has to be either one hour, or one day, or one month. If it is one day or one month, then we have NRLMSIS00e use daily Ap 
      *use_ap_hist = 0;
    }

    else{ // recall the assumption: the time step of the F10.7/Ap file has to be either one hour, or one day, or one month. If it is one hour, then we have NRLMSIS00e use historical Ap. Remember: having a one hour time step on Ap does not mean the Ap file is not 3-hourly (see comments in lin_interpolate_ap_hist)
      *use_ap_hist = 1;

    }

    if (index_to_stop_in_before_interpo == 0){ // this means that the last date in the F10.7 file is older than et_final. This can happen only when dynamic_manual is chosen 

      last_date_input_file_older_final_epoch = 1;
      index_to_stop_in_before_interpo = nb_elements_in_file_f107 - 1;

    }
  }

  else if ( strcmp(src_file, "raid3") == 0 ){
    /* Calculates the x and y arrays to interpolate */
    nb_elements_file(&nb_elements_in_file_raid3, f107_filename, "Two-digit","#ENDOFFILE");
    x_raid3_before_interpo = malloc( nb_elements_in_file_raid3 * sizeof(double) );
    raid3_before_interpo = malloc( nb_elements_in_file_raid3 * sizeof(double) );

    getline(&line, &len, fp);
    for (line_num = 0; line_num<nb_elements_in_file_raid3; line_num++){
      getline(&line, &len, fp);
      sscanf(line, "%s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", yy, doy, &ss, &lat_sat, &long_sat, &alt_sat, &loc_time, &raid3_before_interpo[line_num], &dens_norm400, &dens_norm410, &dens_uncert, &nb_data_points );
      //      printf("%s |  %s | %f %f %f %f %f %e %e %e %e  %d",yy, doy, ss, lat_sat, long_sat, alt_sat, loc_time, raid3_before_interpo[line_num], dens_norm400, dens_norm410, dens_uncert, nb_data_points );
      strcpy(text,yy);
      strcat(text,"-");
      strcat(text,doy);
      strcat(text,"/");
      strcat(text,"00");
      strcat(text,":");
      strcat(text,"00");
      str2et_c (text, &x_raid3_before_interpo_temp );
      x_raid3_before_interpo_temp = x_raid3_before_interpo_temp + ss;
      x_raid3_before_interpo[line_num] = x_raid3_before_interpo_temp;
    }
  }

  
  //  free(line);
  fclose(fp);

  //char times[256];
  //  etprint(et_initial, "ini");

/*   i = 0; */
/*   x_after_interpo[i] = et_oldest_tle_epoch;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*   while (x_after_interpo[i] < et_initial){ */
/*     i = i + 1; */
/*     x_after_interpo[i] = et_oldest_tle_epoch + dt*i / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*     //        etprint(x_after_interpo[i], "tle"); */
	
/*   } // leave this loop when x_after_interpo gets newer than constellation epoch (et_initial) */
/*   int j = 0; */
/*   while (i<nb_time_steps_simu){ */
/*        x_after_interpo[i] = et_initial + dt*j / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*     //            etprint(x_after_interpo[i], "const"); */
/* 	i = i +1; */
/* 	j = j+1; */
/*   } */

  /* Compute y_after_interpo */
  a = malloc( nb_time_steps_simu * sizeof(double) );
  b = malloc( nb_time_steps_simu * sizeof(double) );


  //  char times[256], times2[256];

  /* Linear interpolate F10.7 */
  
  for (i = 1; i < nb_time_steps_simu-1; i++){
    previous_index(&x_min_index, x_f107_before_interpo, x_after_interpo[i], nb_elements_in_file_f107);
    if (i == 1){
      save_first_index_in_f107_before_interpo = x_min_index;
    }
    /* et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times); */
    /* printf("\nx_after_interpo[%d] = %s \n", i, times); */
    /* et2utc_c(x_f107_before_interpo[x_min_index], "C" ,3 ,255 , times); */
    /* et2utc_c(x_f107_before_interpo[x_min_index+1], "C" ,3 ,255 , times2); */
    /* printf("x_f107_before_interpo[x_min_index] = %s || %s\n",  times,  times2); */
    x_min = x_f107_before_interpo[x_min_index];
    
    if (  x_after_interpo[i] - x_f107_before_interpo[index_to_stop_in_before_interpo]  > 0.01 ){
      f107_after_interpo[i] = f107_before_interpo[index_to_stop_in_before_interpo];
    }
    else {

      y_min = f107_before_interpo[x_min_index];
      x_max = x_f107_before_interpo[x_min_index+1];
      y_max = f107_before_interpo[x_min_index+1];
      a[i- 1] = (y_max - y_min) / (x_max - x_min);
      b[i-1] = y_max - a[i-1]*x_max;

      f107_after_interpo[i] = a[i-1]*x_after_interpo[i] + b[i-1];
    }

      /* // !!!!!!!!!!!!!!!!!!! remove storm block below */
      /* if ( (x_after_interpo[i] - et_initial) / 3600.  >= 6 ){ */
      /* 	f107_after_interpo[i] = f107_after_interpo[i] * 1.5; */
      /* } */
      /* // !!!!!!!!!!!!!!!!!!! end of remove storm block below */

/*  	etprint(x_after_interpo[i], ""); */
/* 	printf("%d %f \n",i,f107_after_interpo[i]); */

  }


/*       for (i = 1; i < nb_time_steps_simu-1; i++){ */
/* 	etprint(x_after_interpo[i], ""); */
/* 	printf("%f %d\n", f107_after_interpo[i], i); */
/*       } */

//      exit(0);

 
  if( fabs( x_after_interpo[0] - x_f107_before_interpo[0] ) > 0.01 ) // if the first time in the driver file is different (by more than 0.01s here) from the first time of propagation then the first interpolated value of the driver is equal to the second one. 
    f107_after_interpo[0] = f107_after_interpo[1];
  else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    f107_after_interpo[0] = f107_before_interpo[0];
  if( fabs( x_after_interpo[nb_time_steps_simu-1] - x_f107_before_interpo[index_to_stop_in_before_interpo-1] ) > 0.01 ) // same comments as the two previous lines
    f107_after_interpo[nb_time_steps_simu-1] = f107_after_interpo[nb_time_steps_simu-2];
  else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    f107_after_interpo[nb_time_steps_simu-1] = f107_before_interpo[index_to_stop_in_before_interpo-1];
  /* Done with linear interpolating F10.7 */
  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
    printf("----- (load_options) (lin_interpolate) Done linear interpolating F10.7.\n");
  }
  //  printf("X %d\n", index_to_stop_in_before_interpo);
  ///////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// DONE WITH F10.7 ////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////


  if ( strcmp(src_file, "raid3") == 0 ){
    return 0;
  }
  else{
    ///////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////// Ap /////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    /* Linear interpolate Ap if use_ap_hist = 0, otherwise call lin_interpolate_ap_hist */

    if (*use_ap_hist == 0){

      if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
	printf("----- (load_options) (lin_interpolate) Starting to linear interpolate Ap daily/monthly.\n");
      }

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////// Ap (if time step of Ap file is not one hour but one day or one month /////////////////////////////
      fp = fopen(ap_filename, "r");
      if (fp == NULL){
	printf("***! (load_options) (lin_interpolate) The file %s could not be open. The program will stop. ***!\n", ap_filename); MPI_Finalize(); exit(0);
      }

      nb_elements_in_file_ap = 0;

      if ( ( strcmp(src_file, "omniweb") == 0 ) || ( strcmp(src_file, "dynamic_manual") == 0 ) ){
	/* Calculates the x and y arrays to interpolate */
	if ( strcmp(src_file, "dynamic_manual") == 0 ){
	  nb_elements_file(&nb_elements_in_file_ap, ap_filename, "YEAR","#ENDOFFILE");
	}

	if ( strcmp(src_file, "omniweb") == 0 ){
	  nb_elements_file(&nb_elements_in_file_ap, ap_filename, "YEAR","</pre><hr><HR>");
	}

	x_ap_before_interpo = malloc( nb_elements_in_file_ap * sizeof(double) );
	ap_before_interpo = malloc( nb_elements_in_file_ap * sizeof(double) );

	found_eoh = 0;
	while ( found_eoh == 0 && !feof(fp)) {
	  getline(&line, &len, fp);
	  sscanf(line, "%s", text);
	  if (  strcmp( "YEAR", text  ) == 0 )  {
	    found_eoh = 1;
	  }
	}

	for (line_num = 0; line_num<nb_elements_in_file_ap; line_num++){
	  getline(&line, &len, fp);
	  sscanf(line, "%s %s %s %lf", yy, doy, hh,&driver_temp);
	  strcpy(text,yy);
	  strcat(text,"-");
	  strcat(text,doy);
	  strcat(text,"/");
	  strcat(text,hh);
	  strcat(text,":");
	  strcat(text,"00");
	  str2et_c (text, &x_ap_before_interpo_temp );
	  x_ap_before_interpo[line_num] = x_ap_before_interpo_temp;
	  ap_before_interpo[line_num] = driver_temp;
	  if ( line_num > 0 ){
	    if ( ap_before_interpo[line_num] > ( missing_data_value - 0.01 ) )
	      ap_before_interpo[line_num] = ap_before_interpo[line_num-1];
	  }
/* 	  etprint(x_ap_before_interpo[line_num], ""); */
/* 	  printf("%f\n", ap_before_interpo[line_num]); */
	  if ( x_ap_before_interpo[line_num] > et_final){
	    index_to_stop_in_before_interpo_ap = line_num;
	    break;
	  }

	}

    if (index_to_stop_in_before_interpo_ap == 0){ // this means that the last date in the Ap file is older than et_final. This can happen only when dynamic_manual is chosen 
      last_date_input_file_older_final_epoch_ap = 1;
      index_to_stop_in_before_interpo_ap = nb_elements_in_file_ap - 1;

    }


   
      }

      //      free(line);

      fclose(fp);

      /* Linear interpolate Ap (if NRLMSIS00e uses daily Ap) */
      for (i = 1; i < nb_time_steps_simu-1; i++){
	//printf("%d %d\n", i, nb_time_steps_simu-1);
	previous_index(&x_min_index, x_ap_before_interpo, x_after_interpo[i], nb_elements_in_file_ap);
	/*     et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times);  */
	/*     printf("\nx_after_interpo[%d] = %s \n", i, times); */
	/*     et2utc_c(x_ap_before_interpo[x_min_index], "C" ,3 ,255 , times);  */
	/*     et2utc_c(x_ap_before_interpo[x_min_index+1], "C" ,3 ,255 , times2);  */
	/*     printf("x_ap_before_interpo[x_min_index] = %s || %s\n",  times,  times2); */
	x_min = x_ap_before_interpo[x_min_index];
	if (  x_after_interpo[i] - x_ap_before_interpo[index_to_stop_in_before_interpo_ap]  > 0.01 ){
	  ap_after_interpo[i] = ap_before_interpo[index_to_stop_in_before_interpo_ap];
	}
	else {

	  y_min = ap_before_interpo[x_min_index];
	  x_max = x_ap_before_interpo[x_min_index+1];
	  y_max = ap_before_interpo[x_min_index+1];
	  a[i- 1] = (y_max - y_min) / (x_max - x_min);
	  b[i-1] = y_max - a[i-1]*x_max;

	  ap_after_interpo[i] = a[i-1]*x_after_interpo[i] + b[i-1];

	}
	if (ap_after_interpo[i] > 990){ // this happens for example in omniweb if the value of ap is unknonw they put 999
	  ap_after_interpo[i] = ap_after_interpo[i-1];
	}
/* 	etprint(x_after_interpo[i], ""); */
/* 	printf("%f %f \n", ap_before_interpo[x_min_index],ap_after_interpo[i]); */
      }

      if( fabs( x_after_interpo[0] - x_ap_before_interpo[0] ) > 0.01 ) // if the first time in the driver file is different (by more than 0.01s here) from the first time of propagation then the first interpolated value of the driver is equal to the second one.
	ap_after_interpo[0] = ap_after_interpo[1];
      else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
	ap_after_interpo[0] = ap_before_interpo[0];
      if( fabs( x_after_interpo[nb_time_steps_simu-1] - x_ap_before_interpo[index_to_stop_in_before_interpo_ap] ) > 0.01 ) // same comments as the two previous lines
	ap_after_interpo[nb_time_steps_simu-1] = ap_after_interpo[nb_time_steps_simu-2];
      else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
	ap_after_interpo[nb_time_steps_simu-1] = ap_before_interpo[index_to_stop_in_before_interpo_ap];

      /* Done with linear interpolating Ap (if NRLMSIS00e uses daily Ap) */
      if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
	printf("----- (load_options) (lin_interpolate) Done linear interpolating Ap daily/monthly.\n");
      }

    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////// Done with Ap (if time step of Ap file is not one hour but one day or one month //////////////////

    else{
      if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
	printf("----- (load_options) (lin_interpolate) Starting to linear interpolate historical Ap.\n");
      }

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////// Ap (if time step of Ap file is one hour and not one day or one month /////////////////////////////

             lin_interpolate_ap_hist(ap_hist_after_interpo,x_after_interpo, ap_filename, src_file, nb_time_steps_simu, initial_epoch,et_oldest_tle_epoch, dt, missing_data_value, iDebugLevel, iProc);     
      for (ddd = 0; ddd<nb_time_steps_simu; ddd++){
	ap_after_interpo[ddd] = ap_hist_after_interpo[0][ddd];

	if (ap_after_interpo[ddd] > 990){ // this happens for example in omniweb if the value of ap is unknonw they put 999
	  ap_after_interpo[ddd] = ap_after_interpo[ddd-1];
	}
      }
      if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
	printf("----- (load_options) (lin_interpolate) Done linear interpolating historical Ap.\n");
      }

    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////// Done with Ap (if time step of Ap file is one hour and not one day or one month ///////////////////


    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// Done with Ap ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////








    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// F10.7 AVERAGE ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
      printf("----- (load_options) (lin_interpolate) Starting to calculate and interpolate F10.7 81 days average.\n");
    }
    f107A_before_interpo = malloc( nb_elements_in_file_f107* sizeof(double)) ; // used to be the following but really we want an upper boud so it's ok to allocate more memory than needed malloc( (index_to_stop_in_before_interpo - save_first_index_in_f107_before_interpo +1 ) * sizeof(double) );
    x_f107A_before_interpo = malloc( nb_elements_in_file_f107* sizeof(double)) ; //malloc( (index_to_stop_in_before_interpo - save_first_index_in_f107_before_interpo +1) * sizeof(double) );
				     //    printf("%d %d %d\n",(index_to_stop_in_before_interpo - save_first_index_in_f107_before_interpo +1 ), index_to_stop_in_before_interpo, save_first_index_in_f107_before_interpo);

    /////////////////////// CALCULATE f107A_before_interpo /////////////////////////////
    /* Initialize sum_81_days as the sum of all F10.7 values from the initial epoch - 40.5 days to the initial epoch + 40.5 days */
    //    char current_time[256];


    int start_index_in_f107_before_interpo_to_calculate_f107_average = save_first_index_in_f107_before_interpo;
    int stop_index_in_f107_before_interpo_to_calculate_f107_average = save_first_index_in_f107_before_interpo + 1;
    int current_index_in_f107_average = 0; //save_first_index_in_f107_before_interpo;
    int current_index_in_f107_before_interpo = save_first_index_in_f107_before_interpo;
    //  // in f107_before_interpo, go back from epoch start to as far in the past as we can until reaching 40.5 days
    while ( ( ( save_first_index_in_f107_before_interpo - start_index_in_f107_before_interpo_to_calculate_f107_average ) * time_step_before_interpo < 40.5 * 24 * 3600. ) && ( start_index_in_f107_before_interpo_to_calculate_f107_average > 0 ) ) {
      sum_81_days = sum_81_days + f107_before_interpo[start_index_in_f107_before_interpo_to_calculate_f107_average];
      start_index_in_f107_before_interpo_to_calculate_f107_average = start_index_in_f107_before_interpo_to_calculate_f107_average - 1;
    }
    //  // in f107_before_interpo, go ahead from epoch start to as far in the future as we can until reaching 40.5 days

    while ( ( ( stop_index_in_f107_before_interpo_to_calculate_f107_average - save_first_index_in_f107_before_interpo ) * time_step_before_interpo < 40.5 * 24 * 3600. ) && ( stop_index_in_f107_before_interpo_to_calculate_f107_average < nb_elements_in_file_f107 - 1 ) ) {
      sum_81_days = sum_81_days + f107_before_interpo[stop_index_in_f107_before_interpo_to_calculate_f107_average];
      stop_index_in_f107_before_interpo_to_calculate_f107_average = stop_index_in_f107_before_interpo_to_calculate_f107_average + 1;
      //      printf("%d %f\n", stop_index_in_f107_before_interpo_to_calculate_f107_average, f107_before_interpo[stop_index_in_f107_before_interpo_to_calculate_f107_average]);

    }
    /* Now that we have initialized sum_81_days, we just have to update it at each time step by removing from it the first element of f107_before_interpo and adding the last element of f107_before_interpo (by first and last I mean in the window time of 81 days) */

    f107A_before_interpo[current_index_in_f107_average] = sum_81_days / ( stop_index_in_f107_before_interpo_to_calculate_f107_average - start_index_in_f107_before_interpo_to_calculate_f107_average - 1 ); // here: current_index_in_f107_average = 0
    x_f107A_before_interpo[current_index_in_f107_average] = x_f107_before_interpo[current_index_in_f107_before_interpo];

    while (current_index_in_f107_average < (index_to_stop_in_before_interpo - save_first_index_in_f107_before_interpo   ) ){
      current_index_in_f107_before_interpo = current_index_in_f107_before_interpo + 1;
      if ( ( start_index_in_f107_before_interpo_to_calculate_f107_average >= 0 ) && ( ( current_index_in_f107_before_interpo - start_index_in_f107_before_interpo_to_calculate_f107_average - 1 ) * time_step_before_interpo >= 40.5 * 24 * 3600. ) ){ // we get in here if there is at least 40.5 days before the current index

    	start_index_in_f107_before_interpo_to_calculate_f107_average = start_index_in_f107_before_interpo_to_calculate_f107_average + 1;
    	sum_81_days = sum_81_days - f107_before_interpo[start_index_in_f107_before_interpo_to_calculate_f107_average];
      }



      if ( stop_index_in_f107_before_interpo_to_calculate_f107_average + 1 < index_to_stop_in_before_interpo_for_f107_average ){ // index_to_stop_in_before_interpo_for_f107_average = index in f107_before_interpo corresponding to end epoch + 40.5 days

	
	stop_index_in_f107_before_interpo_to_calculate_f107_average = stop_index_in_f107_before_interpo_to_calculate_f107_average + 1;
    	sum_81_days = sum_81_days + f107_before_interpo[stop_index_in_f107_before_interpo_to_calculate_f107_average];
      }


      current_index_in_f107_average = current_index_in_f107_average + 1;
      x_f107A_before_interpo[current_index_in_f107_average] = x_f107_before_interpo[current_index_in_f107_before_interpo];
      f107A_before_interpo[current_index_in_f107_average] = sum_81_days / ( stop_index_in_f107_before_interpo_to_calculate_f107_average - start_index_in_f107_before_interpo_to_calculate_f107_average - 1 );
/*       etprint(x_f107A_before_interpo[current_index_in_f107_average], ""); */
/*       printf("%f | %d - %d | %f %f\n", f107A_before_interpo[current_index_in_f107_average] , start_index_in_f107_before_interpo_to_calculate_f107_average, stop_index_in_f107_before_interpo_to_calculate_f107_average, f107_before_interpo[start_index_in_f107_before_interpo_to_calculate_f107_average],f107_before_interpo[stop_index_in_f107_before_interpo_to_calculate_f107_average]); */

/*       etprint(x_f107A_before_interpo[current_index_in_f107_average], ""); */
/*       printf("%f %d\n", f107A_before_interpo[current_index_in_f107_average], current_index_in_f107_average); */


      if (f107A_before_interpo[current_index_in_f107_average] > 300){
	print_error_any_iproc(iProc, "(load_options)(lin_interpolate) A value of F10.7A is greater than 300, something is weird");
	
	MPI_Finalize();exit(0);
      }
      

    }



    if (last_date_input_file_older_final_epoch == 1){ // // this means that the last date in the F10.7 file is older than et_final. This can happen only when dynamic_manual is chosen -> f107A_before_interpo for the index after the last date in the file is constant, equal to the last value calcuilated in the loop previously
       ggg = current_index_in_f107_average;
      while ( x_f107A_before_interpo[ggg]  <= et_final - 3600 ){

	f107A_before_interpo[ggg] = f107A_before_interpo[current_index_in_f107_average];
	ggg = ggg + 1;
      }



    }

    /*   //    char current_time[256]; */
    /*   for (ppp = 0; ppp < (index_to_stop_in_before_interpo - save_first_index_in_f107_before_interpo ); ppp++){ */
    /*   et2utc_c(x_f107A_before_interpo[ppp], "C" ,3 ,255 , current_time); */
    /*   printf("%s: %f\n", current_time, f107A_before_interpo[ppp]); */
    /* } */


    /////////////////////// LINEAR INTERPOLATE f107A_before_interpo /////////////////////////////
    int nb_elements_in_file_f107A = (index_to_stop_in_before_interpo - save_first_index_in_f107_before_interpo ) + 1;
    if (last_date_input_file_older_final_epoch == 1){ // // this means that the last date in the F10.7 file is older than et_final. This can happen only when dynamic_manual is chosen
      nb_elements_in_file_f107A = nb_elements_in_file_f107A + ggg -1 - current_index_in_f107_average ;
    }

    /* Linear interpolate F10.7 average */
    for (i = 1; i < nb_time_steps_simu-1; i++){
      previous_index(&x_min_index, x_f107A_before_interpo, x_after_interpo[i], nb_elements_in_file_f107A);
      x_min = x_f107A_before_interpo[x_min_index];
      //      printf("%f\n", f107A_before_interpo[x_min_index]);
      
      if (  x_after_interpo[i] - x_f107A_before_interpo[nb_elements_in_file_f107A-1]  > 0.01 ){
	f107A_after_interpo[i] = f107A_before_interpo[nb_elements_in_file_f107A-1];
      }
      else {
	y_min = f107A_before_interpo[x_min_index];
	x_max = x_f107A_before_interpo[x_min_index+1];
	y_max = f107A_before_interpo[x_min_index+1];
	a[i- 1] = (y_max - y_min) / (x_max - x_min);
	b[i-1] = y_max - a[i-1]*x_max;
	f107A_after_interpo[i] = a[i-1]*x_after_interpo[i] + b[i-1];
      }
      //	printf("%f %f %f  %d\n",f107A_after_interpo[i], f107A_before_interpo[x_min_index], f107A_before_interpo[2], x_min_index);


    }

    if( fabs( x_after_interpo[0] - x_f107A_before_interpo[0] ) > 0.01 ) // if the first time in the driver file is different (by more than 0.01s here) from the first time of propagation then the first interpolated value of the driver is equal to the second one.
      f107A_after_interpo[0] = f107A_after_interpo[1];
    else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
      f107A_after_interpo[0] = f107A_before_interpo[0];
    if( fabs( x_after_interpo[nb_time_steps_simu-1] - x_f107A_before_interpo[nb_elements_in_file_f107A-1] ) > 0.01 ) // same comments as the two previous lines
      f107A_after_interpo[nb_time_steps_simu-1] = f107A_after_interpo[nb_elements_in_file_f107A-2];
    else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
      f107A_after_interpo[nb_time_steps_simu-1] = f107A_before_interpo[nb_elements_in_file_f107A-1];
    /* Done with linear interpolating F10.7 */


    if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
      printf("----- (load_options) (lin_interpolate) Done calculating and interpolating F10.7 81 days average.\n");
    }
    //    printf("%f %f %f  %d\n",f107A_before_interpo[0], f107A_before_interpo[1], f107A_before_interpo[2], x_min_index);    
     


    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// Done with F10.7 AVERAGE /////////////////////
    ///////////////////////////////////////////////////////////////////////////////////



  }

/*   i = 31892; */
/* etprint(x_after_interpo[i], "");  */
/*   printf("%f\n",f107A_after_interpo[i]); */

/*       for (i = 1; i < nb_time_steps_simu-1; i++){ */
/* 	etprint(x_after_interpo[i], ""); */
/* 	printf("wwq%f %f %f %d\n", f107_after_interpo[i], f107A_after_interpo[i], ap_after_interpo[i], i); */
/*       } */
//     	MPI_Finalize();exit(0);

  free(a);
    free(b);
    free(x_f107_before_interpo);
  //free(x_raid3_before_interpo);

    free(f107_before_interpo);
    if (*use_ap_hist == 0){
      free(ap_before_interpo);
      free(x_ap_before_interpo);  
    }
  free(line);
  if ((iProc == 0) & ( iDebugLevel >= 3 ) ){
    printf("---- (load_options) (lin_interpolate) Done linear interpolating F10.7 and Ap.\n");
  }
  

  return 0;


}

/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           lin_interpolate_swpc
//  Purpose:        Linear interpolation of F10.7 and Ap, as well as derivation and linear interpolation of F10.7A 
//  Assumptions:    - Files are previously (automatically) downloaded from SWPC
//                  - The F10.7 and Ap files must have the same time step: 1 day
//                  - prediciton file has Ap and F10.7 in same file AND 45 days predictions, 9 rows, 5 columns
//                  - the uncertainty file on F10.7 and Ap is assumed to be written before running the propagator. It has to be named ./input/density/density_NRLMSIS00e/sigma_f107_ap.txt and have the correct format 
//                  - the uncertainty file on F10.7 and Ap is assumed to be one day time step, predictions for day+1, day+2, ..., day+45 
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 09/23/2016    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////

int lin_interpolate_swpc(double *f107_after_interpo, // !!!!!!!!! here the 81 d average of 107 is over the 81 days preceding the current time, not centered aound the current time. to use the 81 days centerd around current itme,, the user must not use the "swpc" option but either the "omniweb" or external files -> this will call the functino lin_interpolate, which uses the 81 day average centered around the current time
			 double *f107A_after_interpo,
			 double *ap_after_interpo,
			 double *x_after_interpo,
			 double *sigma_f107_after_interpo, 
			 double *sigma_ap_after_interpo,
			 double *x_sigma_f107_ap_after_interpo,
			 double *use_ap_hist,
			 char **list_filename_obs_swpc_f107,
			 char **list_filename_obs_swpc_ap,
			 char filename_f107_ap_pred[256],
			 int nb_time_steps_simu,
			 char initial_epoch[256],
			 double et_oldest_tle_epoch,
			 char final_epoch[256],
			 double dt,
			 double missing_data_value,
			 int iDebugLevel,
			 int nb_file_obs_swpc_f107,
			 int nb_file_obs_swpc_ap,
			 int swpc_need_observations,
			 int swpc_need_predictions,
			 int nb_ensembles_density,
			 char dir_input_density_msis[1000],
			 double *swpc_et_first_prediction, int iProc,
			 int swpc_final_epoch_more_recent_than_45_days_from_current_day
			 ){

  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
    printf( "----- (load_options) (lin_interpolate_swpc) Just got in lin_interpolate_swpc.\n");
  }


  /* Declarations */
  double mean_f107_45_days = 0;
  double mean_ap_45_days = 0;
  int k;
  double swpc_et_last_observation;
  double et_final_epoch;
  int index_ap_first_pred ;
  int ap_at_least_one_obs = 0; 
  int f107_at_least_one_obs = 0;
  int ifakeobs;
  /* int ooo; */
  /*     char time_x_f107_before_interpo[256];	   char time_x_ap_before_interpo[256]; */
  double et_first_prediction = 0;
  double ap_before_interpo_pred1, ap_before_interpo_pred2, ap_before_interpo_pred3, ap_before_interpo_pred4, ap_before_interpo_pred5;
  double f107_before_interpo_pred1, f107_before_interpo_pred2, f107_before_interpo_pred3, f107_before_interpo_pred4, f107_before_interpo_pred5;
  int nb_day_pred = 0;
  double et_date_pred_temp;
  char date_pred1[256], date_pred2[256], date_pred3[256], date_pred4[256], date_pred5[256], date_pred_temp[256];
  int nrow_pred = 9; // you can't just chage this number, you also need to modify the code accordingly

  int irow;
  int ifile;
  FILE *file_obs=NULL, *file_pred=NULL;
  int  index_nb_f107, index_nb_ap;
  int start_saving_data;
  char text2[256];
  char year_obs[6], month_obs[6], day_obs[6], str_trash[100];
  double float_trash;
  double f107_obs, ap_obs;
  char time_obs[256];
  double et_time_obs;
  char initial_epoch_midnight_minus_81_days[256];
  char initial_epoch_midnight_minus_81_days_temp[256];


  double sum_81_days = 0 ;
  double time_step_before_interpo;
  int x_min_index, save_first_index_in_f107_before_interpo;
  double *f107_before_interpo = NULL;
  double *f107A_before_interpo = NULL;
  double *ap_before_interpo = NULL;
  double *x_f107_before_interpo = NULL;
  double *x_f107A_before_interpo = NULL;
  double *x_ap_before_interpo = NULL;
  double *sigma_f107_before_interpo = NULL;
  double *sigma_ap_before_interpo = NULL;
  double *x_sigma_f107_ap_before_interpo = NULL;
  
  char *line = NULL;
  size_t len = 0;
  int found_eoh = 0;
  char text[256];

  int i;
  int nb_elements_in_file_f107=0, nb_elements_in_file_ap=0;
  double et_initial;  
  double *a = NULL;
  double *b = NULL;
  b = malloc(  ( nb_time_steps_simu + (int)(( 2 * 24 * 3600. ) / dt) )  * sizeof(double) );
  a = malloc( ( nb_time_steps_simu + (int)(( 2 * 24 * 3600. ) / dt) ) * sizeof(double) );


  double x_min, y_min, x_max, y_max;
  // J2000 time of intial and final epoch at midnight
  char initial_epoch_midnight[256], final_epoch_midnight[256];
  str2et_c(initial_epoch, &et_initial); 
  strcpy(initial_epoch_midnight, "");  strcpy(final_epoch_midnight, "");
  //      strncat(initial_epoch_midnight, &initial_epoch[0], 10);
  char oldest_tle_epoch[256];
  et2utc_c(et_oldest_tle_epoch, "ISOC" ,6 ,255 , oldest_tle_epoch);
  strncat(initial_epoch_midnight, &oldest_tle_epoch[0], 10);// tle
  strcat(initial_epoch_midnight, "T00:00:00.000");
  strncat(final_epoch_midnight, &final_epoch[0], 10);
  strcat(final_epoch_midnight, "T00:00:00.000");
  double et_final_epoch_midnight, et_initial_epoch_midnight;
  str2et_c(final_epoch_midnight, &et_final_epoch_midnight);  str2et_c(initial_epoch_midnight, &et_initial_epoch_midnight);

  double et_initial_epoch_midnight_minus_81_days_temp = et_initial_epoch_midnight - 81. * 24 * 3600;
  double et_initial_epoch_midnight_minus_81_days;
  et2utc_c(et_initial_epoch_midnight_minus_81_days_temp, "ISOC", 0, 255, initial_epoch_midnight_minus_81_days_temp);
  strcpy(initial_epoch_midnight_minus_81_days, "");
  strncat(initial_epoch_midnight_minus_81_days, &initial_epoch_midnight_minus_81_days_temp[0], 10);
  strcat(initial_epoch_midnight_minus_81_days, "T00:00:00.000");
  str2et_c(initial_epoch_midnight_minus_81_days, &et_initial_epoch_midnight_minus_81_days);	
  // number of days between epoch start day - 81 days and epoch end day
  int nb_days_initial_minus_81_days_to_final = ceil(( et_final_epoch_midnight - et_initial_epoch_midnight_minus_81_days )) / 3600. / 24 + 1;// ceil is necessary for numerical reasons
  
  int nb_days_initial_to_final = ceil(( et_final_epoch_midnight - et_initial_epoch_midnight )) / 3600. / 24 + 1; // ceil is necessary for numerical reasons

  time_step_before_interpo  = 24 * 3600; // the SWPC files  are assumed to be 1 day time step

  // The number of elements of x_f107_before_interpo and f107_before_interpo is equal to the number of days between epoch start day - 81 days and epoch end day
  nb_elements_in_file_f107 = nb_days_initial_minus_81_days_to_final + 1; // +1 because if only observations (no predictions), we want to have the day right after the peoch day so that we can linear interpolate values during end epoch day
  nb_elements_in_file_ap = nb_days_initial_to_final +1;
  x_f107_before_interpo = malloc( ( nb_elements_in_file_f107 + 4 )* sizeof(double) ); // '+ 4' because in the code i look at the index + 4 for some reason that I dont remember now
  f107_before_interpo = malloc( ( nb_elements_in_file_f107 + 4 ) * sizeof(double) );  // '+ 4' because in the code i look at the index + 4 for some reason that I dont remember now
  x_ap_before_interpo = malloc( ( nb_elements_in_file_ap + 4 ) * sizeof(double) );  // '+ 4' because in the code i look at the index + 4 for some reason that I dont remember now
  ap_before_interpo = malloc( ( nb_elements_in_file_ap + 4 ) * sizeof(double) );  // '+ 4' because in the code i look at the index + 4 for some reason that I dont remember now

    
  // Read observation files
  if (swpc_need_observations == 1){
    if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
      printf("------ (load_options) (lin_interpolate_swpc) Reading observation files\n");
    }

    // // First read the first day of the predction file: this is the day where to stop the observations (if pred start on 09-08 then stop obs at 09-07). We do that because   sometimes the first pred corresponds to the last obs (so basically we ignore the last obs as it is often flawed because it corresponds to the current day)
    // // // Read prediction files
    if (swpc_need_predictions == 1){
      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	printf("------ (load_options) (lin_interpolate_swpc) First reading the first day of the predction file.\n");
      }

      // //  // // open prediction file
      //      printf("<<%s>>\n", filename_f107_ap_pred);
      file_pred = fopen(filename_f107_ap_pred, "r");
      if (file_pred == NULL){
	print_error_any_iproc(iProc, "(load_options)(lin_interpolate_swpc) Could not find filename_f107_ap_pred");
      }
      // // // //Ap block
      // // //  // // skip ap header
      found_eoh = 0;
      while ( ( found_eoh == 0 ) && ( !feof(file_pred) ) ){
	getline(&line, &len, file_pred);
	sscanf(line, "%s", text);
	if (strcmp(text, "45-DAY") == 0){
	  found_eoh = 1;
	}
      }

      // // // // // // read first line of predictions !!! ASSUMPTION: 45 days predictions, 9 rows, 5 columns
      getline(&line, &len, file_pred);      
      //	printf("HH <%s>\n",line);
      sscanf(line, "%s %lf %s %lf %s %lf %s %lf %s %lf", date_pred1, &float_trash, date_pred2, &float_trash, date_pred3, &float_trash, date_pred4, &float_trash, date_pred5, &float_trash);

      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred1[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred1[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred1[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &et_date_pred_temp);	
      et_first_prediction = et_date_pred_temp; 
      fclose(file_pred);
      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	etprint(et_first_prediction, "------ (load_options) (lin_interpolate_swpc) First day of prediction");
      }

      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	printf("------ (load_options) (lin_interpolate_swpc) Done first reading the first day of the predction file.\n");
      }

    }
    else{
      et_first_prediction = -1;
    }
    *swpc_et_first_prediction = et_first_prediction;

    //    // Go over each F10.7 observation file
    index_nb_f107 = 0;
    start_saving_data = 0;

    for (ifile = 0; ifile < nb_file_obs_swpc_f107 ; ifile++){
      // // // open f10.7 observation file
      file_obs = fopen(list_filename_obs_swpc_f107[ifile], "r");
      if (file_obs == NULL){
	printf("***! (load_options)(lin_interpolate_swpc) The F10.7 observation file was not found. The program will stop. !***\n"); MPI_Finalize(); exit(0);
      }
      //    printf("<%s>\n", list_filename_obs_swpc_f107[ifile]);
      // // // skip header
      
      found_eoh = 0;
      while ( found_eoh == 0 && !feof(file_obs)) {
	getline(&line, &len, file_obs);
	sscanf(line, "%s", text);
	strcpy(text2, "");
	strncat(text2, &text[0], 4);
	if (  strcmp( "#---", text2  ) == 0 )  {
	  found_eoh = 1;
	}
      }


      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	printf("------ (load_options) (lin_interpolate_swpc) Reading f10.7 file_obs.\n");
      }

      // // // Read file

      while (!feof(file_obs)){
	getline(&line, &len, file_obs);

	if (feof(file_obs)){
	  /* printf("<%s>\n",list_filename_obs_swpc_f107[ifile]); */
	  /* 	  test_print("P"); */
	  break;
	}

	//	printf("\n<%s>\n", line);
	sscanf(line, "%s %s %s %lf %s", year_obs, month_obs, day_obs, &f107_obs, str_trash);

	if ( strcmp(year_obs, ":Product:") == 0 ){ // the end of the quarter of the year file is reached so skip header of new quarter

	  found_eoh = 0;
	  while ( found_eoh == 0 && !feof(file_obs)) {
	    getline(&line, &len, file_obs);
	    sscanf(line, "%s", text);
	    strcpy(text2, "");
	    strncat(text2, &text[0], 4);
	    if (  strcmp( "#---", text2  ) == 0 )  {
	      found_eoh = 1;
	      getline(&line, &len, file_obs);
	      //    printf("\n<%s>\n", line);
	      sscanf(line, "%s %s %s %lf %s", year_obs, month_obs, day_obs, &f107_obs, str_trash);
	    }
	  }
	}
	//	  printf("XXXXXXXXXXXXXXXXX <%s>", line);
	strcpy(time_obs, year_obs);
	strcat(time_obs, "-");
	strcat(time_obs, month_obs);
	strcat(time_obs, "-");
	strcat(time_obs, day_obs);
	strcat(time_obs, "T00:00:00.000");
	str2et_c(time_obs, &et_time_obs);


	    
	if (fabs(et_time_obs - et_initial_epoch_midnight_minus_81_days) < 0.01){// 0.01 for numerical reasons
	  start_saving_data = 1;
	}

	// // // // start saving values only at the epoch start minus 81 days
	if ( start_saving_data == 1){
	  /* 	    etprint(et_time_obs, "et_time_obs"); */
	  /* 	    etprint(et_initial_epoch_midnight_minus_81_days, "et_initial_epoch_midnight_minus_81_days"); */

	  if (feof(file_obs)){
	    break;
	  }
	 

	  // // // // stop reading file at epoch end or one day before first prediction
	  if ( ifile == nb_file_obs_swpc_f107 - 1 ){
	    if (fabs(et_time_obs - ( et_final_epoch_midnight + 2* 24 * 3600 ) ) < 0.01){// 0.01 for numerical reasons
	      break;
	    }
	    if (swpc_need_predictions == 1){
	      if ( et_time_obs > et_first_prediction  + 0.01 ){
		break;
	      }
	    }
	  }
	  //	    printf("AAAAAAA = %d %d\n", index_nb_f107, nb_elements_in_file_f107);

	  x_f107_before_interpo[index_nb_f107] = et_time_obs;

	  f107_before_interpo[index_nb_f107] = f107_obs;
	  /* 	      etprint(x_f107_before_interpo[index_nb_f107], "" ); */
	  /* 	           printf("%f %d\n", f107_before_interpo[index_nb_f107], index_nb_f107); */
	  if ( index_nb_f107 > 0 ){
	    if ( f107_before_interpo[index_nb_f107] > ( missing_data_value - 0.01 ) )
	      f107_before_interpo[index_nb_f107] = f107_before_interpo[index_nb_f107-1];
	  }
	  /*  	etprint(x_f107_before_interpo[index_nb_f107] ,""); */
	  /* 	printf("%f (%d)\n", f107_before_interpo[index_nb_f107], index_nb_f107); */


	  //	printf("<%s>\n", line);
	  f107_at_least_one_obs = 1;
	  index_nb_f107 = index_nb_f107 + 1;
	}
      }

      fclose(file_obs);
      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	printf("------ (load_options) (lin_interpolate_swpc) Done reading f10.7 file_obs.\n");
      }

    } // end of going over f107 observaation files
    // below was commented because the situtation it covers can not happen since if the final epoch is more recent than the last date in the F10.7 file then we use predictions
    /* 	if (nb_elements_in_file_f107 >= index_nb_f107){ // this happens when the final epoch is more recent than the last date in the F10.7 file. In this case, for all dates more recent than the final date in f107 file, f107 is set to constant equal to the last value in the f107 file */

    /* 	  for (k = index_nb_f107; k < nb_elements_in_file_f107; k++){ */
    /* 	    f107_before_interpo[k] = f107_before_interpo[index_nb_f107-1]; */
    /* 	    x_f107_before_interpo[k] = x_f107_before_interpo[index_nb_f107-1]; */
	      
    /* 	  } */

    /* 	} */
    // end of below was commented because the situtation it covers can not happen since if the final epoch is more recent than the last date in the F10.7 file then we use predictions


    // Ap
    index_nb_ap = 0;
    start_saving_data = 0;
    //    // Go over each ap observation file
    for (ifile = 0; ifile < nb_file_obs_swpc_ap ; ifile++){
      // // // open ap observation file
      file_obs = fopen(list_filename_obs_swpc_ap[ifile], "r");
      if (file_obs == NULL){
	printf("***! (load_options)(lin_interpolate_swpc) The Ap observation file was not found. The program will stop. !***\n"); MPI_Finalize(); exit(0);
      }
      //	  printf("<%s>\n", list_filename_obs_swpc_ap[ifile]);
      // // // skip header
      
      found_eoh = 0;
      while ( found_eoh == 0 && !feof(file_obs)) {
	getline(&line, &len, file_obs);
	sscanf(line, "%s %s", text, text2);
	if (  strcmp( "Date", text2  ) == 0 )  {
	  found_eoh = 1;
	}
      }
      


      // // // Read file
      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	printf("------ (load_options) (lin_interpolate_swpc) Reading Ap file_obs.\n");
      }

      while (!feof(file_obs)){
	getline(&line, &len, file_obs);
	
	if (feof(file_obs)){
	  break;
	}

	//	printf("\n<%s>\n", line);
	sscanf(line, "%s %s", text, text2);
	//	
	  while (text[0] == '#'){ // sometimes the file has extra lines of comments - can't anticipate, it changes...

	    getline(&line, &len, file_obs);
	    sscanf(line, "%s %s", text, text2);


	  }
	  sscanf(line, "%s %s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", year_obs, month_obs, day_obs, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &ap_obs);
	  /* printf("ddddd <%s> %s %s %f\n", year_obs, month_obs, day_obs, ap_obs); */

	if ( strcmp(year_obs, ":Product:") == 0 ){ // the end of the quarter of the year file is reached so skip header of new quarter
	  found_eoh = 0;
	  while ( found_eoh == 0 && !feof(file_obs)) {
	    getline(&line, &len, file_obs);
	    sscanf(line, "%s %s", text, text2);
	    //	    	      printf("<%s>\n", text2);
	    if (  strcmp( "Date", text2  ) == 0 )  {
	      found_eoh = 1;
	      getline(&line, &len, file_obs);
	      //		  printf("\n<%s>\n", line);
	      sscanf(line, "%s %s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", year_obs, month_obs, day_obs, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &ap_obs);
	      //		printf("%s %s %s %f\n", year_obs, month_obs, day_obs, ap_obs);
	    }
	  }
	}
	strcpy(time_obs, year_obs);
	strcat(time_obs, "-");
	strcat(time_obs, month_obs);
	strcat(time_obs, "-");
	strcat(time_obs, day_obs);
	strcat(time_obs, "T00:00:00.000");
	str2et_c(time_obs, &et_time_obs);
	if (fabs(et_time_obs - et_initial_epoch_midnight) < 0.01){// 0.01 for numerical reasons
	  start_saving_data = 1;
	}
	    
	// // // // start saving values only at the epoch start 
	if ( start_saving_data == 1){
	  if (feof(file_obs)){
	    break;
	  }

	  // // // // stop reading file at epoch end or one day before first prediction
	  if ( ifile == nb_file_obs_swpc_ap - 1 ){
	    if (fabs(et_time_obs - ( et_final_epoch_midnight + 2* 24 * 3600 ) ) < 0.01){// 0.01 for numerical reasons
	      break;
	    }
	    if (swpc_need_predictions == 1){
	      if ( et_time_obs > et_first_prediction - 0.01 ){
		break;
	      }
	    }
	  }
	  //	    printf("AAAAAAA = %d %d\n", index_nb_f107, nb_elements_in_file_f107);
	  x_ap_before_interpo[index_nb_ap] = et_time_obs;
	  ap_before_interpo[index_nb_ap] = ap_obs;
	  if ( index_nb_ap > 0 ){
	    if ( ap_before_interpo[index_nb_ap] > ( missing_data_value - 0.01 ) ){
	      ap_before_interpo[index_nb_ap] = ap_before_interpo[index_nb_ap-1];
	    }
	  }


	   
	  //	    	    printf("<%s> !!!!!!! %f || %d-%d\n", line, ap_before_interpo[index_nb_ap], index_nb_ap, nb_elements_in_file_ap);
	  ap_at_least_one_obs = 1;
	  index_nb_ap = index_nb_ap + 1;

	}

      }
      fclose(file_obs);
      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	printf("------ (load_options) (lin_interpolate_swpc) Done reading Ap file_obs.\n");
      }

    }

    if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
      printf("------ (load_options) (lin_interpolate_swpc) Done reading observation files.\n");
    }

    // below was commented because the situtation it covers can not happen since if the final epoch is more recent than the last date in the F10.7 file then we use predictions
    /* 	if (nb_elements_in_file_ap >= index_nb_ap){ // this happens when the final epoch is more recent than the last date in the F10.7 file. In this case, for all dates more recent than the final date in ap file, ap is set to constant equal to the last value in the ap file */

    /* 	  for (k = index_nb_ap; k < nb_elements_in_file_ap; k++){ */
    /* 	    ap_before_interpo[k] = ap_before_interpo[index_nb_ap-1]; */
    /* 	    x_ap_before_interpo[k] = x_ap_before_interpo[index_nb_ap-1]; */
	      
    /* 	  } */

    /* 	} */

    // end of below was commented because the situtation it covers can not happen since if the final epoch is more recent than the last date in the F10.7 file then we use predictions
  }


  //  Read prediction files
  if (swpc_need_predictions == 1){

    if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
      printf("------ (load_options) (lin_interpolate_swpc) Reading prediction files.\n");
    }

    // // open prediction file
    file_pred = fopen(filename_f107_ap_pred, "r");
    if (file_pred == NULL){
      print_error_any_iproc(iProc, "(load_options)(lin_interpolate_swpc) Could not fine filename_f107_ap_pred");
    }

    // // Ap block
    // // // skip ap header
    found_eoh = 0;
    while ( ( found_eoh == 0 ) && ( !feof(file_pred) ) ){
      getline(&line, &len, file_pred);
      sscanf(line, "%s", text);
      if (strcmp(text, "45-DAY") == 0){
	found_eoh = 1;
      }
    }


    // // // read predictions. !!! ASSUMPTION: 45 days predictions, 9 rows, 5 columns
    for (irow = 0; irow < nrow_pred; irow++){

      getline(&line, &len, file_pred);
      if (irow == 0){
	sscanf(line, "%s %lf %s %lf %s %lf %s %lf %s %lf", date_pred1, &float_trash, date_pred2, &float_trash, date_pred3, &float_trash, date_pred4, &float_trash, date_pred5, &float_trash);

	strcpy(date_pred_temp, "20");
	strncat(date_pred_temp, &date_pred1[0] + 5, 2);
	strcat(date_pred_temp, " ");
	strncat(date_pred_temp, &date_pred1[0] + 2, 3);
	strcat(date_pred_temp, " ");
	strncat(date_pred_temp, &date_pred1[0], 2);
	strcat(date_pred_temp, " 00:00:00.000");
	str2et_c(date_pred_temp, &et_date_pred_temp);
	if (swpc_need_observations == 1){
	  // //  Sometimes there is more than one day between the first pred and the last obs (if obs have not been updated by SWPC for instance). If it is the case, add "fake" observations from the last obs to the first pred, assuming ap constant. Ex: if last obs is 04/26 (ap = 12) and first predi s 04/29 thif an add fake obs for 04/27 and 04/28, assuming a ap constant for these 2 days: 12.

	  if ( (ap_at_least_one_obs == 1) && ( fabs(et_date_pred_temp - x_ap_before_interpo[index_nb_ap-1] ) > 25*3600 ) ){
	    ifakeobs = 0;
	    x_ap_before_interpo[index_nb_ap] = x_ap_before_interpo[index_nb_ap-1] + 24*3600;
	    ap_before_interpo[index_nb_ap] = ap_before_interpo[index_nb_ap - 1];
	    //	      et2utc_c(x_ap_before_interpo[index_nb_ap], "ISOC", 3, 255, time_x_ap_before_interpo);
	    //	      printf("%s %f (%d)\n", time_x_ap_before_interpo,  ap_before_interpo[index_nb_ap],  index_nb_ap);
	    while ( fabs(et_date_pred_temp - x_ap_before_interpo[index_nb_ap] ) > 25*3600 ){
	      index_nb_ap = index_nb_ap + 1;
	      x_ap_before_interpo[index_nb_ap] = x_ap_before_interpo[index_nb_ap-1] + 24*3600;
	      ap_before_interpo[index_nb_ap] = ap_before_interpo[index_nb_ap - 1];
	      //	      et2utc_c(x_ap_before_interpo[index_nb_ap], "ISOC", 3, 255, time_x_ap_before_interpo);
	      //	      printf("while: %s %f (%d)\n", time_x_ap_before_interpo,  ap_before_interpo[index_nb_ap],  index_nb_ap);

	    }
	    index_nb_ap = index_nb_ap + 1;
	  }
	}
 
      }

      nb_day_pred = index_nb_ap; // initalization but will be calculated after


      index_ap_first_pred = index_nb_ap;
      sscanf(line, "%s %lf %s %lf %s %lf %s %lf %s %lf", date_pred1, &ap_before_interpo_pred1, date_pred2, &ap_before_interpo_pred2, date_pred3, &ap_before_interpo_pred3, date_pred4, &ap_before_interpo_pred4, date_pred5, &ap_before_interpo_pred5);
     
      //	printf("%s %f %s %f %s %f %s %f %s %f\n", date_pred1, ap_before_interpo[index_nb_ap], date_pred2, ap_before_interpo[index_nb_ap+1], date_pred3, ap_before_interpo[index_nb_ap+2], date_pred4, ap_before_interpo[index_nb_ap+3], date_pred5, ap_before_interpo[index_nb_ap+4]);
      	
      // // // // convert date into seconds past J2000
      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred1[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred1[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred1[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_ap_before_interpo[index_nb_ap]);
      ap_before_interpo[index_nb_ap] = ap_before_interpo_pred1;
      if ((x_ap_before_interpo[index_nb_ap] >= et_final_epoch_midnight + 24 * 3600) || (index_nb_ap ==  nb_elements_in_file_ap + 4 - 1)){
	if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	  etprint(x_ap_before_interpo[index_nb_ap], "------ (load_options) (lin_interpolate_swpc) Last day of prediction");
	}
	nb_day_pred = index_nb_ap - nb_day_pred + 1 ; // nb_day_pred is equal to epoch end day - first value of prediction file + 1
	break;
      }

      //	  char date_pred_temp2[256];
      //	  et2utc_c(x_ap_before_interpo[index_nb_ap], "ISOC", 0, 255, date_pred_temp2);

      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred2[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred2[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred2[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_ap_before_interpo[index_nb_ap+1]);
      ap_before_interpo[index_nb_ap+1] = ap_before_interpo_pred2;
      if ( ( x_ap_before_interpo[index_nb_ap+1] >= et_final_epoch_midnight + 24 * 3600) || (index_nb_ap+1 ==  nb_elements_in_file_ap + 4 - 1)){
	if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	  etprint(x_ap_before_interpo[index_nb_ap+1], "------ (load_options) (lin_interpolate_swpc) Last day of prediction");
	}

	nb_day_pred = index_nb_ap - nb_day_pred + 1  + 1; // nb_day_pred is equal to epoch end day - first value of prediction file + 1

	break;
      }

      //	  et2utc_c(x_ap_before_interpo[index_nb_ap+1], "ISOC", 0, 255, date_pred_temp2);
      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred3[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred3[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred3[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_ap_before_interpo[index_nb_ap+2]);
      ap_before_interpo[index_nb_ap+2] = ap_before_interpo_pred3;
      if ((x_ap_before_interpo[index_nb_ap+2] >= et_final_epoch_midnight + 24 * 3600)  || (index_nb_ap+2 ==  nb_elements_in_file_ap + 4 - 1)){
	if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	  etprint(x_ap_before_interpo[index_nb_ap+2], "------ (load_options) (lin_interpolate_swpc) Last day of prediction");
	}

	nb_day_pred = index_nb_ap - nb_day_pred + 1  + 2; // nb_day_pred is equal to epoch end day - first value of prediction file + 1

	break;
      }

      //	  et2utc_c(x_ap_before_interpo[index_nb_ap+2], "ISOC", 0, 255, date_pred_temp2);

      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred4[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred4[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred4[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");

      str2et_c(date_pred_temp, &x_ap_before_interpo[index_nb_ap+3]);
	  
      ap_before_interpo[index_nb_ap+3] = ap_before_interpo_pred4;

	  
      if ((x_ap_before_interpo[index_nb_ap+3] >= et_final_epoch_midnight + 24 * 3600) || (index_nb_ap+3 ==  nb_elements_in_file_ap + 4 - 1)){

	if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	  etprint(x_ap_before_interpo[index_nb_ap+3], "------ (load_options) (lin_interpolate_swpc) Last day of prediction");
	}

	nb_day_pred = index_nb_ap - nb_day_pred + 1  + 3; // nb_day_pred is equal to epoch end day - first value of prediction file + 1
	//	    printf("%f %f %f %f\n", ap_before_interpo[index_nb_ap],ap_before_interpo[index_nb_ap+1],ap_before_interpo[index_nb_ap+2],ap_before_interpo[index_nb_ap+3]);
	break;
      }

      //	  et2utc_c(x_ap_before_interpo[index_nb_ap+3], "ISOC", 0, 255, date_pred_temp2);

      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred5[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred5[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred5[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_ap_before_interpo[index_nb_ap+4]);
      ap_before_interpo[index_nb_ap+4] = ap_before_interpo_pred5;
      if ((x_ap_before_interpo[index_nb_ap+4] >= et_final_epoch_midnight + 24 * 3600) || (index_nb_ap+4 ==  nb_elements_in_file_ap + 4 - 1)){
	if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	  etprint(x_ap_before_interpo[index_nb_ap+4], "------ (load_options) (lin_interpolate_swpc) Last day of prediction");
	}
	nb_day_pred = index_nb_ap - nb_day_pred + 1 +4; // nb_day_pred is equal to epoch end day - first value of prediction file + 1
	break;
      }
      //	  	  et2utc_c(x_ap_before_interpo[index_nb_ap+4], "ISOC", 0, 255, date_pred_temp2);

      //	printf("%f %f %f %f %f\n", ap_before_interpo[index_nb_ap],  ap_before_interpo[index_nb_ap+1],  ap_before_interpo[index_nb_ap+2],  ap_before_interpo[index_nb_ap+3],  ap_before_interpo[index_nb_ap+4] );
      mean_ap_45_days = mean_ap_45_days + ap_before_interpo[index_nb_ap] + ap_before_interpo[index_nb_ap+1] + ap_before_interpo[index_nb_ap+2]+ap_before_interpo[index_nb_ap+3]+ap_before_interpo[index_nb_ap+4]; // if the user chose swpc for the f107 and ap, and that the final epoch is beyong +45 days in the future (limit of swpc predictions) then values of ap and f107 beyong +45 days are set to constant, equal to the mean of the predicted values by SWPC (over the 45 days)

      index_nb_ap = index_nb_ap + 5;  // be careful if using index_nb_ap in the rest of the script




    }



    /* 	for (k = 0; k < ; k++){ */
	  
    /* 	  etprint( x_ap_before_interpo[i], "") */
    /* 	    printf("%f\n"ap_before_interpo[k]); */
    /* 	  } */
    /* 	exit(0); */


    mean_ap_45_days  = mean_ap_45_days / 45;


    swpc_et_last_observation = x_ap_before_interpo[index_ap_first_pred]-24*3600.;
    if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
      etprint(swpc_et_last_observation, "------ (load_options) (lin_interpolate_swpc) Last day of observation");
      pti(nb_day_pred, "------ (load_options) (lin_interpolate_swpc) Number of prediction days");
    }

    // // F107 block
    // // // skip f107 header
    found_eoh = 0;
    while ( ( found_eoh == 0 ) && ( !feof(file_pred) ) ){
      getline(&line, &len, file_pred);
      sscanf(line, "%s", text);
      if (strcmp(text, "45-DAY") == 0){
	found_eoh = 1;
      }
    }
    
    // // // read predictions. !!! ASSUMPTION: 45 days predictions, 9 rows, 5 columns
    for (irow = 0; irow < nrow_pred; irow++){
      getline(&line, &len, file_pred);


      //  See comment below
      if (irow == 0){
	sscanf(line, "%s %lf %s %lf %s %lf %s %lf %s %lf", date_pred1, &float_trash, date_pred2, &float_trash, date_pred3, &float_trash, date_pred4, &float_trash, date_pred5, &float_trash);
	strcpy(date_pred_temp, "20");
	strncat(date_pred_temp, &date_pred1[0] + 5, 2);
	strcat(date_pred_temp, " ");
	strncat(date_pred_temp, &date_pred1[0] + 2, 3);
	strcat(date_pred_temp, " ");
	strncat(date_pred_temp, &date_pred1[0], 2);
	strcat(date_pred_temp, " 00:00:00.000");
	str2et_c(date_pred_temp, &et_date_pred_temp);
	// //  Sometimes the first pred corresponds to the last obs so in that case we overwrite the last obs (as it is often flawed because it corresponds to the current day)
	if  (swpc_need_observations == 1){
	  // //  Sometimes there is more than one day between the first pred and the last obs (if obs have not been updated by SWPC for instance). If it is the case, add "fake" observations from the last obs to the first pred, assuming F10.7 constant. Ex: if last obs is 04/26 (F10.7 = 90) and first pred is 04/29 thif an add fake obs for 04/27 and 04/28, assuming a F10.7 constant for these 2 days: 90.
	  if ( (f107_at_least_one_obs == 1) && fabs(et_date_pred_temp - x_f107_before_interpo[index_nb_f107-1] ) > 25*3600 ){
	    ifakeobs = 0;
	    x_f107_before_interpo[index_nb_f107] = x_f107_before_interpo[index_nb_f107-1] + 24*3600;
	    f107_before_interpo[index_nb_f107] = f107_before_interpo[index_nb_f107 - 1];
	    /* char time_x_f107_before_interpo[256]; */
	    /*       printf("%f %d\n", f107_before_interpo[index_nb_f107], index_nb_f107); */
	    /* et2utc_c(x_f107_before_interpo[index_nb_f107], "ISOC", 3, 255, time_x_f107_before_interpo); */
	    /*      printf("%s %f (%d)\n", time_x_f107_before_interpo,  f107_before_interpo[index_nb_f107],  index_nb_f107); */
	    while ( fabs(et_date_pred_temp - x_f107_before_interpo[index_nb_f107] ) > 25*3600 ){
	      index_nb_f107 = index_nb_f107 + 1;
	      x_f107_before_interpo[index_nb_f107] = x_f107_before_interpo[index_nb_f107-1] + 24*3600;
	      f107_before_interpo[index_nb_f107] = f107_before_interpo[index_nb_f107 - 1];
	      /* char time_x_f107_before_interpo[256]; */
	      /* et2utc_c(x_f107_before_interpo[index_nb_f107], "ISOC", 3, 255, time_x_f107_before_interpo); */
	      /* 	      printf("while: %s %f (%d)\n", time_x_f107_before_interpo,  f107_before_interpo[index_nb_f107],  index_nb_f107); */

	    }
	    index_nb_f107 = index_nb_f107 + 1;

	  }
	}
      }

      /* et2utc_c(x_f107_before_interpo[index_nb_f107-1], "ISOC", 3, 255, time_x_f107_before_interpo); */
      /* printf("out %s %f (%d-%d)\n", time_x_f107_before_interpo, f107_before_interpo[index_nb_f107-1], index_nb_f107, nb_elements_in_file_f107); */
      /* exit(0); */
      sscanf(line, "%s %lf %s %lf %s %lf %s %lf %s %lf", date_pred1, &f107_before_interpo_pred1, date_pred2, &f107_before_interpo_pred2, date_pred3, &f107_before_interpo_pred3, date_pred4, &f107_before_interpo_pred4, date_pred5, &f107_before_interpo_pred5);

      //printf("%s %f %s %f %s %f %s %f %s %f\n", date_pred1, f107_before_interpo[index_nb_f107], date_pred2, f107_before_interpo[index_nb_f107+1], date_pred3, f107_before_interpo[index_nb_f107+2], date_pred4, f107_before_interpo[index_nb_f107+3], date_pred5, f107_before_interpo[index_nb_f107+4]);

      // // // // convert date into seconds past J2000
      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred1[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred1[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred1[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_f107_before_interpo[index_nb_f107]);
      f107_before_interpo[index_nb_f107] = f107_before_interpo_pred1;
      /* etprint(x_f107_before_interpo[index_nb_f107], "x_f107_before_interpo[index_nb_f107]"); */
      /* printf("%f %d\n",f107_before_interpo[index_nb_f107],index_nb_f107); */

      if ((x_f107_before_interpo[index_nb_f107] >= et_final_epoch_midnight + 24*3600) || (index_nb_f107 ==  nb_elements_in_file_f107 + 4 - 1)){

	break;
      }

      //	  char date_pred_temp2[256];
      //	  et2utc_c(x_f107_before_interpo[index_nb_f107], "ISOC", 0, 255, date_pred_temp2);

      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred2[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred2[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred2[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_f107_before_interpo[index_nb_f107+1]);
      f107_before_interpo[index_nb_f107+1] = f107_before_interpo_pred2;
      /* etprint(x_f107_before_interpo[index_nb_f107+1], "x_f107_before_interpo[index_nb_f107+1]"); */
      /* printf("%f %d\n",f107_before_interpo[index_nb_f107+1],index_nb_f107+1); */

      if ((x_f107_before_interpo[index_nb_f107+1] >= et_final_epoch_midnight + 24*3600) || (index_nb_f107+1 ==  nb_elements_in_file_f107 + 4 - 1)){
	//	    test_print("A");
	break;
      }

      //	  et2utc_c(x_f107_before_interpo[index_nb_f107+1], "ISOC", 0, 255, date_pred_temp2);

      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred3[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred3[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred3[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_f107_before_interpo[index_nb_f107+2]);
      f107_before_interpo[index_nb_f107+2] = f107_before_interpo_pred3;
      if ((x_f107_before_interpo[index_nb_f107+2] >= et_final_epoch_midnight + 24 * 3600) || (index_nb_f107+2 ==  nb_elements_in_file_f107 + 4 - 1)){
	break;
      }

      //	  et2utc_c(x_f107_before_interpo[index_nb_f107+2], "ISOC", 0, 255, date_pred_temp2);

      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred4[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred4[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred4[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_f107_before_interpo[index_nb_f107+3]);
      f107_before_interpo[index_nb_f107+3] = f107_before_interpo_pred4;
      if ((x_f107_before_interpo[index_nb_f107+3] >= et_final_epoch_midnight + 24 * 3600) || (index_nb_f107+3 ==  nb_elements_in_file_f107 + 4 - 1)){
	break;
      }

      //	  et2utc_c(x_f107_before_interpo[index_nb_f107+3], "ISOC", 0, 255, date_pred_temp2);

      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred5[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred5[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred5[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_f107_before_interpo[index_nb_f107+4]);
      f107_before_interpo[index_nb_f107+4] = f107_before_interpo_pred5;
      if ((x_f107_before_interpo[index_nb_f107+4] >= et_final_epoch_midnight + 24 * 3600) || (index_nb_f107+4 ==  nb_elements_in_file_f107 + 4 - 1)){
	break;
      }

      //	  et2utc_c(x_f107_before_interpo[index_nb_f107+4], "ISOC", 0, 255, date_pred_temp2)
      mean_f107_45_days = mean_f107_45_days + f107_before_interpo[index_nb_f107] + f107_before_interpo[index_nb_f107+1] + f107_before_interpo[index_nb_f107+2]+f107_before_interpo[index_nb_f107+3]+f107_before_interpo[index_nb_f107+4]; // if the user chose swpc for the f107 and ap, and that the final epoch is beyong +45 days in the future (limit of swpc predictions) then values of ap and f107 beyong +45 days are set to constant, equal to the mean of the predicted values by SWPC (over the 45 days)
      index_nb_f107 = index_nb_f107 + 5; // be careful if using index_nb_f107 in the rest of the script

    }

    fclose(file_pred);
    if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
      printf("------ (load_options) (lin_interpolate_swpc) Done reading prediction files\n");
    }

    //	printf("%f %f %f %f %d\n",f107_before_interpo[index_nb_f107],f107_before_interpo[index_nb_f107+1],f107_before_interpo[index_nb_f107+2],f107_before_interpo[index_nb_f107+3], index_nb_f107);

  }

  mean_f107_45_days = mean_f107_45_days / 45; // if the user chose swpc for the f107 and ap, and that the final epoch is beyong +45 days in the future (limit of swpc predictions) then values of ap and f107 beyong +45 days are set to constant, equal to the mean of the predicted values by SWPC (over the 45 days)

  /*   int ooo; */
  /* for (ooo = 0; ooo < nb_elements_in_file_ap; ooo++){ */
  /*   et2utc_c(x_ap_before_interpo[ooo], "ISOC", 0, 255, date_pred_temp); */
  /*   printf("%s %f\n", date_pred_temp, ap_before_interpo[ooo]); */
  /* } */
  /* exit(0); */
  if ( time_step_before_interpo > 23 * 3600 ){ // recall the assumption: the time step of the F10.7/Ap file has to be either one hour, or one day, or one month. If it is one day or one month, then we have NRLMSIS00e use daily Ap. Files from SWPC are assumed to be one day time step so Ap dauly is used, not ap hist
    *use_ap_hist = 0;
  }



  ///////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// LINEAR INTERPOLATION OF F10.7 //////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
    printf("----- (load_options) (lin_interpolate_swpc) Linear interpolating F10.7...\n");
  }

/*   i = 0; */
/*   x_after_interpo[i] = et_oldest_tle_epoch;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*   while (x_after_interpo[i] < et_initial){ */
/*     i = i + 1; */
/*     x_after_interpo[i] = et_oldest_tle_epoch + dt*i / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*     //   etprint(x_after_interpo[i], "tle"); */
	
/*   } // leave this loop when x_after_interpo gets newer than constellation epoch (et_initial) */
/*   int j = 0; */
/*   while (i<nb_time_steps_simu){ */
/*     x_after_interpo[i] = et_initial + dt*j / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*     //            etprint(x_after_interpo[i], "const"); */

/*     i = i +1; */
/*     j = j+1; */
/*   } */

  /* //  char times[256]; */
  /* //  pti(nb_time_steps_simu, "nb_time_steps_simu"); */
  /* for (i = 0; i < nb_time_steps_simu; i++){ */
  /* 	x_after_interpo[i] = et_initial + dt*i / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
  /* 	//    et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times); */
  /* 	//   pti(i, "i"); */
  /* } */

  /* /\\* Compute y_after_interpo *\\/ */
  //  a = malloc( ( nb_time_steps_simu + (int)(( 2 * 24 * 3600. ) / dt) ) * sizeof(double) );
  //  b = malloc( ( nb_time_steps_simu + (int)(( 2 * 24 * 3600. ) / dt) ) * sizeof(double) );
  /* /\\* Linear interpolate F10.7 *\\/ */

  int last_element_x_f107_before_interpo;
  if (swpc_final_epoch_more_recent_than_45_days_from_current_day == 1){

    last_element_x_f107_before_interpo = index_nb_f107-1;
    for(k = index_nb_f107-1; k < nb_elements_in_file_f107; k ++){
      x_f107_before_interpo[k] = x_f107_before_interpo[last_element_x_f107_before_interpo];
      f107_before_interpo[k] = mean_f107_45_days;// if the user chose swpc for the f107 and ap, and that the final epoch is beyong +45 days in the future (limit of swpc predictions) then values of ap and f107 beyong +45 days are set to constant, equal to the mean of the predicted values by SWPC (over the 45 days)
    }
    
  }
  else{
    last_element_x_f107_before_interpo = nb_elements_in_file_f107-1;
  }


  //  etprint(x_f107_before_interpo[index_nb_f107-1], "");
  for (i = 1; i < nb_time_steps_simu-1; i++){
    if ((swpc_final_epoch_more_recent_than_45_days_from_current_day == 1) && ( x_after_interpo[i] > x_f107_before_interpo[index_nb_f107-1]) ){// if the user chose swpc for the f107 and ap, and that the final epoch is beyong +45 days in the future (limit of swpc predictions) then values of ap and f107 beyong +45 days are set to constant, equal to the mean of the predicted values by SWPC (over the 45 days)
		
      f107_after_interpo[i] = mean_f107_45_days;
/* 			 	etprint(x_after_interpo[i] ,"a"); */
/* 			 	etprint(x_f107_before_interpo[x_min_index] ,"b"); */

    }
    else{
      //   printf("%d %d\n", i, nb_time_steps_simu-2);
      previous_index(&x_min_index, x_f107_before_interpo, x_after_interpo[i], nb_elements_in_file_f107);
      if (i == 1){
	save_first_index_in_f107_before_interpo = x_min_index;
      }
      /*  et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times);  */
      /* printf("\nx_after_interpo[%d] = %s \n", i, times);  */
      /* et2utc_c(x_f107_before_interpo[x_min_index], "C" ,3 ,255 , times);  */
      /* et2utc_c(x_f107_before_interpo[x_min_index+1], "C" ,3 ,255 , times2);  */
      /* printf("x_f107_before_interpo[x_min_index] = %s || %s (x_min_index = %d)\n",  times,  times2, x_min_index); exit(0); */
      x_min = x_f107_before_interpo[x_min_index];
      /* 	etprint(x_f107_before_interpo[last_element_x_f107_before_interpo] ,"file"); */

      /* 	etprint(x_after_interpo[i] ,"interpo"); */
      /* 	printf("%d\n",last_element_x_f107_before_interpo); */
      if (  x_after_interpo[i] - x_f107_before_interpo[last_element_x_f107_before_interpo]  > 0.01 ){
	f107_after_interpo[i] = f107_before_interpo[last_element_x_f107_before_interpo];
	//      	  printf("%f %d\n", f107_before_interpo[last_element_x_f107_before_interpo], last_element_x_f107_before_interpo);

      }
      else {
	//	print_test();

	y_min = f107_before_interpo[x_min_index];

	x_max = x_f107_before_interpo[x_min_index+1];
	y_max = f107_before_interpo[x_min_index+1];
	//	  printf("%f %f\n", f107_before_interpo[x_min_index+1],f107_before_interpo[x_min_index]);
	a[i- 1] = (y_max - y_min) / (x_max - x_min);
	b[i-1] = y_max - a[i-1]*x_max;
	//          printf("%f %d | %f\n", f107_before_interpo[x_min_index], x_min_index, f107_before_interpo[x_min_index+1]);
	f107_after_interpo[i] = a[i-1]*x_after_interpo[i] + b[i-1];
	
	//  	  printf("%f %d %f %f\n", f107_after_interpo[i], i,a[i-1],  b[i-1]);
      }
      /* 		etprint(x_after_interpo[i] ,"a");  */
      /* 		printf("%f\n", f107_before_interpo[x_min_index]); */


    }
/*     			 	etprint(x_after_interpo[i] ,"a"); */
/*     			 	etprint(x_f107_before_interpo[x_min_index] ,"b"); */
/*     				printf("%f %f %d\n", f107_after_interpo[i], f107_before_interpo[x_min_index], x_min_index); */
    //	printf("%f %f %f %f %d\n",f107_before_interpo[86],f107_before_interpo[86+1],f107_before_interpo[86+2],f107_before_interpo[86+3], 86);
  }
  //      exit(0);

 
  if( fabs( x_after_interpo[0] - x_f107_before_interpo[0] ) > 0.01 ) // if the first time in the driver file is different (by more than 0.01s here) from the first time of propagation then the first interpolated value of the driver is equal to the second one.// not correct but it;s ok (should not take index 0 for x_f107_before_interpo but index that corresponds to epoch initial)
    f107_after_interpo[0] = f107_after_interpo[1];
  else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    f107_after_interpo[0] = f107_before_interpo[0];
  if( fabs( x_after_interpo[nb_time_steps_simu-1] - x_f107_before_interpo[last_element_x_f107_before_interpo] ) > 0.01 ) // same comments as the two previous lines
    f107_after_interpo[nb_time_steps_simu-1] = f107_after_interpo[nb_time_steps_simu-2];
  else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    f107_after_interpo[nb_time_steps_simu-1] = f107_before_interpo[last_element_x_f107_before_interpo];
  /* /\\* Done with linear interpolating F10.7 *\\/ */
  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
    printf("----- (load_options) (lin_interpolate_swpc) Done linear interpolating F10.7.\n");
  }

  //  exit(0);
  /* DONE WITH F10.7 //////////////////////////////// */
  

  ///////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// LINEAR INTERPOLATION OF AP /////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  if (*use_ap_hist == 0){
    int last_element_x_ap_before_interpo;
    if (swpc_final_epoch_more_recent_than_45_days_from_current_day == 1){
      last_element_x_ap_before_interpo = index_nb_ap-1;
    }
    else{
      last_element_x_ap_before_interpo = nb_elements_in_file_ap-1;
    }


    if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
      printf("----- (load_options) (lin_interpolate_swpc) Linear interpolating Ap daily/monthly...\n");
    }

    // //   Linear interpolate Ap (if NRLMSIS00e uses daily Ap)
    for (i = 1; i < nb_time_steps_simu-1; i++){
      if ((swpc_final_epoch_more_recent_than_45_days_from_current_day == 1) && ( x_after_interpo[i] > x_ap_before_interpo[index_nb_ap-1])) {// if the user chose swpc for the f107 and ap, and that the final epoch is beyong +45 days in the future (limit of swpc predictions) then values of ap and f107 beyong +45 days are set to constant, equal to the mean of the predicted values by SWPC (over the 45 days)
		

	ap_after_interpo[i] = mean_ap_45_days;
      }
      else{
	previous_index(&x_min_index, x_ap_before_interpo, x_after_interpo[i], nb_elements_in_file_ap);
	/* /\\*     et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times);  *\\/ */
	/* /\\*     printf("\nx_after_interpo[%d] = %s \n", i, times); *\\/ */
	/* /\\*     et2utc_c(x_ap_before_interpo[x_min_index], "C" ,3 ,255 , times);  *\\/ */
	/* /\\*     et2utc_c(x_ap_before_interpo[x_min_index+1], "C" ,3 ,255 , times2);  *\\/ */
	/* /\\*     printf("x_ap_before_interpo[x_min_index] = %s || %s\n",  times,  times2); *\\/ */
	x_min = x_ap_before_interpo[x_min_index];
	if (  x_after_interpo[i] - x_ap_before_interpo[last_element_x_ap_before_interpo]  > 0.01 ){
	  ap_after_interpo[i] = ap_before_interpo[last_element_x_ap_before_interpo];
	}
	else {

	  y_min = ap_before_interpo[x_min_index];
	  x_max = x_ap_before_interpo[x_min_index+1];
	  y_max = ap_before_interpo[x_min_index+1];
	  a[i- 1] = (y_max - y_min) / (x_max - x_min);
	  b[i-1] = y_max - a[i-1]*x_max;

	  ap_after_interpo[i] = a[i-1]*x_after_interpo[i] + b[i-1];
	}

      }
	    
      //		 	etprint(x_ap_before_interpo[x_min_index] ,"b");
      /* 	etprint(x_after_interpo[i] ,"a"); */
      /* 				printf("%f %f (%d)\n", ap_after_interpo[i], ap_before_interpo[x_min_index], x_min_index); */
      /* 				printf("%f\n", ap_after_interpo[i]); */


      //      ap_after_interpo[i] = 2* ap_after_interpo[i];// !!!!!!!!!!!!!!!!!!!!!!!!!!! REMOVE THIE LINE
      //      printf("%f\n", ap_after_interpo[i]);
    }



 
    if( fabs( x_after_interpo[0] - x_ap_before_interpo[0] ) > 0.01 ) // if the first time in the driver file is different (by more than 0.01s here) from the first time of propagation then the first interpolated value of the driver is equal to the second one.
      ap_after_interpo[0] = ap_after_interpo[1];
    else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
      ap_after_interpo[0] = ap_before_interpo[0];
    if( fabs( x_after_interpo[nb_time_steps_simu-1] - x_ap_before_interpo[last_element_x_ap_before_interpo] ) > 0.01 ) // same comments as the two previous lines
      ap_after_interpo[nb_time_steps_simu-1] = ap_after_interpo[nb_time_steps_simu-2];
    else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
      ap_after_interpo[nb_time_steps_simu-1] = ap_before_interpo[last_element_x_ap_before_interpo];

    /* /\\* Done with linear interpolating Ap (if NRLMSIS00e uses daily Ap) *\\/ */
    if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
      printf("----- (load_options) (lin_interpolate_swpc) Done linear interpolating Ap daily/monthly.\n");
    }
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////// Done with Ap (if time step of Ap file is not one hour but one day or one month //////////////////


  else{

    printf("***! (load_options)(lin_interpolate_swpc) If using SWPC, the observations and predictions have to be every one day. The program will stop.\n!***"); MPI_Finalize();exit(0);
  }

  /* char tt[256]; */
  /* for (i = 0; i < nb_time_steps_simu; i++){ */
  /*   et2utc_c(x_after_interpo[i], "ISOC", 3, 255, tt); */
  /*   printf("%s: %f\n", tt, ap_after_interpo[i] ); */
    
  /* } */


  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////// Done with Ap ////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////// F10.7 AVERAGE ///////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  /////// NOTE: the statements below are taken form the function lin_interpolate and are adapted to this function. Therefore, a few things might look weird. But it works correctly.
  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
    printf("----- (load_options) (lin_interpolate_swpc) Starting to calculate and interpolate F10.7 81 days average.\n");
  }
  f107A_before_interpo = malloc( (nb_elements_in_file_f107 - save_first_index_in_f107_before_interpo + 1 ) * sizeof(double) );
  x_f107A_before_interpo = malloc( (nb_elements_in_file_f107 - save_first_index_in_f107_before_interpo + 1 ) * sizeof(double) );

  /////////////////////// CALCULATE f107A_before_interpo /////////////////////////////

  //  /\* Initialize sum_81_days as the sum of all F10.7 values from the initial epoch - 81 days to the initial epoch  *\/
  //    char current_time[256];
  int start_index_in_f107_before_interpo_to_calculate_f107_average = save_first_index_in_f107_before_interpo;
  int stop_index_in_f107_before_interpo_to_calculate_f107_average = save_first_index_in_f107_before_interpo + 1;
  int current_index_in_f107_average = 0; //save_first_index_in_f107_before_interpo;
  int current_index_in_f107_before_interpo = save_first_index_in_f107_before_interpo;
  //  // in f107_before_interpo, go back from epoch start to as far in the past as we can until reaching 81 days
  while ( ( ( save_first_index_in_f107_before_interpo - start_index_in_f107_before_interpo_to_calculate_f107_average ) * time_step_before_interpo < 81 * 24 * 3600. ) && ( start_index_in_f107_before_interpo_to_calculate_f107_average > 0 ) ) {
    sum_81_days = sum_81_days + f107_before_interpo[start_index_in_f107_before_interpo_to_calculate_f107_average];
    start_index_in_f107_before_interpo_to_calculate_f107_average = start_index_in_f107_before_interpo_to_calculate_f107_average - 1;
  }

  //    /\* Now that we have initialized sum_81_days, we just have to update it at each time step by removing from it the first element of f107_before_interpo and adding the last element of f107_before_interpo (by first and last I mean in the window time of 81 days) *\/

  f107A_before_interpo[current_index_in_f107_average] = sum_81_days / ( stop_index_in_f107_before_interpo_to_calculate_f107_average - start_index_in_f107_before_interpo_to_calculate_f107_average - 1 ); // here: current_index_in_f107_average = 0
  x_f107A_before_interpo[current_index_in_f107_average] = x_f107_before_interpo[current_index_in_f107_before_interpo];

  /*    char current_time[256];  */
  /* et2utc_c(x_f107A_before_interpo[current_index_in_f107_average], "C" ,3 ,255 , current_time);  */
  /* printf("%s: %f | %d\n", current_time, f107A_before_interpo[current_index_in_f107_average], current_index_in_f107_average);  */

  while (current_index_in_f107_average < (nb_elements_in_file_f107 - save_first_index_in_f107_before_interpo  -1 ) ){
    current_index_in_f107_before_interpo = current_index_in_f107_before_interpo + 1;
    if ( ( start_index_in_f107_before_interpo_to_calculate_f107_average >= 0 ) && ( ( current_index_in_f107_before_interpo - start_index_in_f107_before_interpo_to_calculate_f107_average - 1 ) * time_step_before_interpo >= 81 * 24 * 3600. ) ){ // we get in here if there is at least 81 days before the current index

      start_index_in_f107_before_interpo_to_calculate_f107_average = start_index_in_f107_before_interpo_to_calculate_f107_average + 1;
      sum_81_days = sum_81_days - f107_before_interpo[start_index_in_f107_before_interpo_to_calculate_f107_average];
    }

    if ( stop_index_in_f107_before_interpo_to_calculate_f107_average + 1 < nb_elements_in_file_f107 ){
      stop_index_in_f107_before_interpo_to_calculate_f107_average = stop_index_in_f107_before_interpo_to_calculate_f107_average + 1;
      sum_81_days = sum_81_days + f107_before_interpo[stop_index_in_f107_before_interpo_to_calculate_f107_average];
    }


    current_index_in_f107_average = current_index_in_f107_average + 1;
    x_f107A_before_interpo[current_index_in_f107_average] = x_f107_before_interpo[current_index_in_f107_before_interpo];
    f107A_before_interpo[current_index_in_f107_average] = sum_81_days / ( stop_index_in_f107_before_interpo_to_calculate_f107_average - start_index_in_f107_before_interpo_to_calculate_f107_average - 1 );

  }

 
  /*            char current_time[256]; */
  /*            int ppp; */
  /*         for (ppp = 0; ppp < (nb_elements_in_file_f107 - save_first_index_in_f107_before_interpo ); ppp++){ */
  /*         et2utc_c(x_f107A_before_interpo[ppp], "C" ,3 ,255 , current_time); */
  /*         printf("%s: %f\n", current_time, f107A_before_interpo[ppp]); */
  /*       } */


  /////////////////////// LINEAR INTERPOLATE f107A_before_interpo /////////////////////////////
  int nb_elements_in_file_f107A = (nb_elements_in_file_f107 - save_first_index_in_f107_before_interpo );
  //  Linear interpolate F10.7 average
  for (i = 1; i < nb_time_steps_simu-1; i++){
    previous_index(&x_min_index, x_f107A_before_interpo, x_after_interpo[i], nb_elements_in_file_f107A);
    x_min = x_f107A_before_interpo[x_min_index];
    if (  x_after_interpo[i] - x_f107A_before_interpo[nb_elements_in_file_f107A-1]  > 0.01 ){
      f107A_after_interpo[i] = f107A_before_interpo[nb_elements_in_file_f107A-1];
    }
    else {
      y_min = f107A_before_interpo[x_min_index];
      x_max = x_f107A_before_interpo[x_min_index+1];
      y_max = f107A_before_interpo[x_min_index+1];
      a[i- 1] = (y_max - y_min) / (x_max - x_min);
      b[i-1] = y_max - a[i-1]*x_max;
      f107A_after_interpo[i] = a[i-1]*x_after_interpo[i] + b[i-1];
    }

    //		 	etprint(x_f107A_before_interpo[x_min_index] ,"b");
    //			printf("%f %f (%d)\n", f107A_after_interpo[i], f107A_before_interpo[x_min_index], x_min_index);
    /* 		etprint(x_after_interpo[i] ,"a"); */

    /* 		printf("%f\n", f107A_after_interpo[i]); */

  }

  if( fabs( x_after_interpo[0] - x_f107A_before_interpo[0] ) > 0.01 ) // if the first time in the driver file is different (by more than 0.01s here) from the first time of propagation then the first interpolated value of the driver is equal to the second one.
    f107A_after_interpo[0] = f107A_after_interpo[1];
  else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    f107A_after_interpo[0] = f107A_before_interpo[0];
  if( fabs( x_after_interpo[nb_time_steps_simu-1] - x_f107A_before_interpo[nb_elements_in_file_f107A-1] ) > 0.01 ) // same comments as the two previous lines
    f107A_after_interpo[nb_time_steps_simu-1] = f107A_after_interpo[nb_elements_in_file_f107A-2];
  else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    f107A_after_interpo[nb_time_steps_simu-1] = f107A_before_interpo[nb_elements_in_file_f107A-1];
  //  Done with linear interpolating F10.7

  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
    printf("----- (load_options) (lin_interpolate_swpc) Done calculating and interpolating F10.7 81 days average.\n");
  }

  /* char tt[256]; */
  /* for (i = 0; i < nb_time_steps_simu; i++){ */
  /*   et2utc_c(x_after_interpo[i], "ISOC", 3, 255, tt); */
  /*   printf("%s: %f\n", tt, f107A_after_interpo[i] ); */
    
  /* } */


  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////// Done with F10.7 AVERAGE /////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// UNCERTAINTY FILE ON F10.7 AND AP ////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  if (swpc_need_predictions){


    if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
      printf("----- (load_options) (lin_interpolate) Reading and interpolating uncertainties in F10.7, Ap, and F10.7A.\n");
    }


    if ( nb_ensembles_density > 0 ){

      FILE *file_uncertainty = NULL;
      char filename_uncertainty[1100];
      // Read uncerainty file
      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	printf("------ (load_options) (lin_interpolate) Reading uncertainties in F10.7, Ap, and F10.7A.\n");
      }

      //newstructure
      /*       strcpy(filename_uncertainty, dir_input_density_msis); */
      strcpy(filename_uncertainty, "");
      //newstructure
      strcat(filename_uncertainty, "/sigma_f107_ap.txt");
      file_uncertainty = fopen(filename_uncertainty, "r");
      if (file_uncertainty == NULL){
	printf("***! (load_options)(lin_interpolate_swpc) The file for the uncertainty on F10.7 and Ap was not found. The program will stop. !***\n"); MPI_Finalize(); exit(0);
      }
      // // Skip header
      found_eoh = 0;
      while ( (found_eoh == 0) && (!feof(file_uncertainty)) ){
	getline(&line, &len, file_uncertainty);
	sscanf(line, "%s", text);
	if (strcmp(text, "#Time(day)") == 0){
	  found_eoh = 1;
	}
      }
      if (found_eoh == 0){
	printf("***! (load_options)(lin_interpolate_swpc) No data has been found in the file for the uncertainty on F10.7 and Ap was not found. The program will stop. !***\n"); MPI_Finalize(); exit(0);
      }

      int nb_day_pred_in_file_uncertainty = 0;
      // IMPORTANT: the uncertainty a midnight of the first day of pred is 0, the uncertainty at midnight of the second day of pred is the first value of the uncertainty file. And so on. Then linear interpolation.
      sigma_f107_before_interpo = malloc( 46 * sizeof(double) ); // ASSUMPTIONS: uncertainty file is 45 day predictions.
      sigma_ap_before_interpo = malloc( 46 * sizeof(double) ); // ASSUMPTIONS: uncertainty file is 45 day predictions
      x_sigma_f107_ap_before_interpo = malloc( 46 * sizeof(double) ); // ASSUMPTIONS: uncertainty file is 45 day predictions
	
      x_sigma_f107_ap_before_interpo[0] = 0; // the uncertainty a midnight of the first day of pred is 0
      int ipred=1;
      while ( (ipred < nb_day_pred + 1) && (!feof(file_uncertainty)) ){
	getline(&line, &len, file_uncertainty);
	nb_day_pred_in_file_uncertainty = nb_day_pred_in_file_uncertainty + 1;
	sscanf(line, "%lf %lf %lf", &x_sigma_f107_ap_before_interpo[ipred], &sigma_f107_before_interpo[ipred], &sigma_ap_before_interpo[ipred]);
	x_sigma_f107_ap_before_interpo[ipred] = x_sigma_f107_ap_before_interpo[ipred] * 24 * 3600;

	if (feof(file_uncertainty)){
	  printf("***! The number of days in the file with the uncertainty on F10.7 and Ap (%d) is smaller than the number of days in the future that the user wants to propagate for (%d). The program will stop. !***\n", nb_day_pred_in_file_uncertainty, nb_day_pred); MPI_Finalize(); exit(0);
	}
	//	      	printf("%f %f %f\n", x_sigma_f107_ap_before_interpo[ipred], sigma_f107_before_interpo[ipred], sigma_ap_before_interpo[ipred]);
	ipred = ipred + 1;
      }
      fclose(file_uncertainty);

      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	printf("------ (load_options) (lin_interpolate) Done reading uncertainties in F10.7, Ap, and F10.7A.\n");
      }
      // // linear interpolate the uncertainties
      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	printf("------ (load_options) (lin_interpolate) Interpolating uncertainties in F10.7, Ap, and F10.7A.\n");
      }

      str2et_c(final_epoch, &et_final_epoch);
      int nb_time_steps_in_prediction = ceil( ( ( et_final_epoch - et_first_prediction ) / dt * 2 ) )  + 1;
      for (i = 0; i < nb_time_steps_in_prediction; i++){

	x_sigma_f107_ap_after_interpo[i] =  dt*i / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method
	//       ptd(x_sigma_f107_ap_after_interpo[i], "time pred:");
      }
      //  exit(0);
      for (i = 1; i < nb_time_steps_in_prediction; i++){
	previous_index(&x_min_index, x_sigma_f107_ap_before_interpo, x_sigma_f107_ap_after_interpo[i], nb_day_pred+1);
	//     printf("%d %f\n", x_min_index, x_sigma_f107_ap_after_interpo[i]);
	// F10.7
	x_min = x_sigma_f107_ap_before_interpo[x_min_index];
	y_min = sigma_f107_before_interpo[x_min_index];
	x_max = x_sigma_f107_ap_before_interpo[x_min_index+1];
	y_max = sigma_f107_before_interpo[x_min_index+1];

	a[i- 1] = (y_max - y_min) / (x_max - x_min);
	b[i-1] = y_max - a[i-1]*x_max;
	sigma_f107_after_interpo[i] = a[i-1]*x_sigma_f107_ap_after_interpo[i] + b[i-1];
	//        printf("%f %f %f %f %d\n", sigma_f107_after_interpo[i], x_sigma_f107_ap_after_interpo[i], y_min, y_max, i);

	// Ap
	y_min = sigma_ap_before_interpo[x_min_index];
	y_max = sigma_ap_before_interpo[x_min_index+1];
	a[i- 1] = (y_max - y_min) / (x_max - x_min);
	b[i-1] = y_max - a[i-1]*x_max;

	sigma_ap_after_interpo[i] = a[i-1]*x_sigma_f107_ap_after_interpo[i] + b[i-1];
            

      }
      //            exit(0);
      sigma_f107_after_interpo[0] = sigma_f107_before_interpo[0]; // equal to 0
      sigma_ap_after_interpo[0] = sigma_ap_before_interpo[0]; // equal to 0
      /* sigma_f107_after_interpo[nb_time_steps_in_prediction-1] = sigma_f107_before_interpo[nb_day_pred]; */
      /* sigma_ap_after_interpo[nb_time_steps_in_prediction-1] = sigma_ap_before_interpo[nb_day_pred]; */


      free(sigma_f107_before_interpo);
      free(sigma_ap_before_interpo);


    }

    if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
      printf("------ (load_options) (lin_interpolate) Done interpolating uncertainties in F10.7, Ap, and F10.7A.\n");
    }

    if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
      printf("----- (load_options) (lin_interpolate) Done reading and interpolating uncertainties in F10.7, Ap, and F10.7A.\n");
    }


  }
  /*       /////////////////////////////////////////////////////////////////////////////////// */
  /*       //////////////////////////// DONE with UNCERTAINTY FILE ON F10.7 AND AP /////////// */
  /*       /////////////////////////////////////////////////////////////////////////////////// */

  free(a);
  free(b);
  free(x_f107_before_interpo);

  free(x_ap_before_interpo);
		  
  free(ap_before_interpo);
  free(f107_before_interpo);

  if ((iProc == 0) & ( iDebugLevel >= 3 ) ){
    printf("---- (load_options) (lin_interpolate_swpc) Done linear interpolating F10.7 and Ap.\n");
  }

  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
    printf("----- (load_options) (lin_interpolate_swpc) Just got out of lin_interpolate_swpc.\n");
  }

  return 0;


}




/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           lin_interpolate_swpc_mod
//  Purpose:        Linear interpolation of F10.7 and Ap, as well as derivation and linear interpolation of F10.7A 
//  Assumptions:    - Files are previously (automatically) downloaded from SWPC
//                  - The F10.7 and Ap files must have the same time step: 1 day
//                  - prediciton file has Ap and F10.7 in same file AND 45 days predictions, 9 rows, 5 columns
//                  - the modification file on F10.7 and Ap is assumed to be written before running the propagator. It has to be put in ./input/density/density_NRLMSIS00e/ and have the correct format 
//                  - the modification file on F10.7 and Ap is assumed to be one day time step, predictions for day+1, day+2, ..., day+45 
// Notes:           This function has been copy paster from lin_interpolate_swpc and modified from it. So there might be inconsitencies in the comments
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 11/17/2016    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////

int lin_interpolate_swpc_mod(double *f107_after_interpo,
			     double *f107A_after_interpo,
			     double *ap_after_interpo,
			     double *x_after_interpo,
			     double *mod_f107_after_interpo, 
			     double *mod_ap_after_interpo,
			     double *x_mod_f107_ap_after_interpo,
			     double *use_ap_hist,
			     char **list_filename_obs_swpc_f107,
			     char **list_filename_obs_swpc_ap,
			     char filename_f107_ap_pred[256],
			     int nb_time_steps_simu,
			     char initial_epoch[256],
			     char final_epoch[256],
			     double dt,
			     double missing_data_value,
			     int iDebugLevel,
			     int nb_file_obs_swpc_f107,
			     int nb_file_obs_swpc_ap,
			     int swpc_need_observations,
			     int swpc_need_predictions,
			     double *swpc_et_first_prediction, int iProc, char filename_f107_ap_mod[1000]
			     ){

  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
    printf( "----- (load_options) (lin_interpolate_swpc_mod) Just got in lin_interpolate_swpc_mod.\n");
  }


  /* Declarations */
    int nb_time_steps_in_prediction = 0;
  double swpc_et_last_observation;
  double et_final_epoch;
  int index_ap_first_pred ;
  int ap_at_least_one_obs = 0; 
  int f107_at_least_one_obs = 0;
  int ifakeobs;
  /* int ooo; */
  /*     char time_x_f107_before_interpo[256];	   char time_x_ap_before_interpo[256]; */
  double et_first_prediction = 0;
  double ap_before_interpo_pred1, ap_before_interpo_pred2, ap_before_interpo_pred3, ap_before_interpo_pred4, ap_before_interpo_pred5;
  double f107_before_interpo_pred1, f107_before_interpo_pred2, f107_before_interpo_pred3, f107_before_interpo_pred4, f107_before_interpo_pred5;
  int nb_day_pred = 0;
  double et_date_pred_temp;
  char date_pred1[256], date_pred2[256], date_pred3[256], date_pred4[256], date_pred5[256], date_pred_temp[256];
  int nrow_pred = 9; // you can't just chage this number, you also need to modify the code accordingly

  int irow;
  int ifile;
  FILE *file_obs=NULL, *file_pred=NULL;
  int  index_nb_f107, index_nb_ap;
  int start_saving_data;
  char text2[256];
  char year_obs[6], month_obs[6], day_obs[6], str_trash[100];
  double float_trash;
  double f107_obs, ap_obs;
  char time_obs[256];
  double et_time_obs;
  char initial_epoch_midnight_minus_81_days[256];


  double sum_81_days = 0 ;
  double time_step_before_interpo;
  int x_min_index, save_first_index_in_f107_before_interpo;
  double *f107_before_interpo = NULL;
  double *f107A_before_interpo = NULL;
  double *ap_before_interpo = NULL;
  double *x_f107_before_interpo = NULL;
  double *x_f107A_before_interpo = NULL;
  double *x_ap_before_interpo = NULL;
  double *mod_f107_before_interpo = NULL;
  double *mod_ap_before_interpo = NULL;
  double *x_mod_f107_ap_before_interpo = NULL;
  
  char *line = NULL;
  size_t len = 0;
  int found_eoh = 0;
  char text[256];

  int i;
  int nb_elements_in_file_f107=0, nb_elements_in_file_ap=0;
  double et_initial;  
  double *a = NULL;
  double *b = NULL;
  b = malloc(  ( nb_time_steps_simu + (int)(( 2 * 24 * 3600. ) / dt) )  * sizeof(double) );
  a = malloc( ( nb_time_steps_simu + (int)(( 2 * 24 * 3600. ) / dt) ) * sizeof(double) );


  double x_min, y_min, x_max, y_max;
  // J2000 time of intial and final epoch at midnight
  char initial_epoch_midnight[256], final_epoch_midnight[256];
  str2et_c(initial_epoch, &et_initial);
  strcpy(initial_epoch_midnight, "");  strcpy(final_epoch_midnight, "");
  strncat(initial_epoch_midnight, &initial_epoch[0], 10);

  strcat(initial_epoch_midnight, "T00:00:00.000");
  strncat(final_epoch_midnight, &final_epoch[0], 10);
  strcat(final_epoch_midnight, "T00:00:00.000");
  double et_final_epoch_midnight, et_initial_epoch_midnight;
  str2et_c(final_epoch_midnight, &et_final_epoch_midnight);  str2et_c(initial_epoch_midnight, &et_initial_epoch_midnight);

  double et_initial_epoch_midnight_minus_81_days = et_initial_epoch_midnight - 81. * 24 * 3600;
  et2utc_c(et_initial_epoch_midnight_minus_81_days, "ISOC", 0, 255, initial_epoch_midnight_minus_81_days);

  // number of days between epoch start day - 81 days and epoch end day
  int nb_days_initial_minus_81_days_to_final = ceil(( et_final_epoch_midnight - et_initial_epoch_midnight_minus_81_days )) / 3600. / 24 + 1;// ceil is necessary for numerical reasons
  
  int nb_days_initial_to_final = ceil(( et_final_epoch_midnight - et_initial_epoch_midnight )) / 3600. / 24 + 1; // ceil is necessary for numerical reasons



  time_step_before_interpo  = 24 * 3600; // the SWPC files  are assumed to be 1 day time step

  // The number of elements of x_f107_before_interpo and f107_before_interpo is equal to the number of days between epoch start day - 81 days and epoch end day
  nb_elements_in_file_f107 = nb_days_initial_minus_81_days_to_final + 1; // +1 because if only observations (no predictions), we want to have the day right after the peoch day so that we can linear interpolate values during end epoch day
  nb_elements_in_file_ap = nb_days_initial_to_final +1;

  x_f107_before_interpo = malloc( nb_elements_in_file_f107 * sizeof(double) );
  f107_before_interpo = malloc( nb_elements_in_file_f107 * sizeof(double) );
  x_ap_before_interpo = malloc( nb_elements_in_file_ap * sizeof(double) );
  ap_before_interpo = malloc( nb_elements_in_file_ap * sizeof(double) );


  // Read observation files
  if (swpc_need_observations == 1){
    if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
      printf("------ (load_options) (lin_interpolate_swpc_mod) Reading observation files\n");
    }

    // // First read the first day of the predction file: this is the day where to stop the observations (if pred start on 09-08 then stop obs at 09-07). We do that because   sometimes the first pred corresponds to the last obs (so basically we ignore the last obs as it is often flawed because it corresponds to the current day)
    // // // Read prediction files

    if (swpc_need_predictions == 1){
      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	printf("------ (load_options) (lin_interpolate_swpc_mod) First reading the first day of the predction file.\n");
      }

      // //  // // open prediction file
      //      printf("<%s>\n", filename_f107_ap_pred);
      file_pred = fopen(filename_f107_ap_pred, "r");
      if (file_pred == NULL){
	print_error_any_iproc(iProc, "(load_options)(lin_interpolate_swpc_mod) Could not find filename_f107_ap_pred");
      }
      // // // //Ap block
      // // //  // // skip ap header
      found_eoh = 0;
      while ( ( found_eoh == 0 ) && ( !feof(file_pred) ) ){
	getline(&line, &len, file_pred);
	sscanf(line, "%s", text);
	if (strcmp(text, "45-DAY") == 0){
	  found_eoh = 1;
	}
      }

      // // // // // // read first line of predictions !!! ASSUMPTION: 45 days predictions, 9 rows, 5 columns
      getline(&line, &len, file_pred);      
      //	printf("HH <%s>\n",line);
      sscanf(line, "%s %lf %s %lf %s %lf %s %lf %s %lf", date_pred1, &float_trash, date_pred2, &float_trash, date_pred3, &float_trash, date_pred4, &float_trash, date_pred5, &float_trash);

      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred1[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred1[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred1[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &et_date_pred_temp);	
      et_first_prediction = et_date_pred_temp; 

      fclose(file_pred);
      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	etprint(et_first_prediction, "------ (load_options) (lin_interpolate_swpc_mod) First day of prediction");
      }

      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	printf("------ (load_options) (lin_interpolate_swpc_mod) Done first reading the first day of the predction file.\n");
      }

    }
    else{
      et_first_prediction = -1;
    }
    *swpc_et_first_prediction = et_first_prediction;

    //    // Go over each F10.7 observation file
    index_nb_f107 = 0;
    start_saving_data = 0;

    for (ifile = 0; ifile < nb_file_obs_swpc_f107 ; ifile++){
      // // // open f10.7 observation file
      file_obs = fopen(list_filename_obs_swpc_f107[ifile], "r");
      if (file_obs == NULL){
	printf("***! (load_options)(lin_interpolate_swpc_mod) The F10.7 observation file was not found. The program will stop. !***\n"); MPI_Finalize(); exit(0);
      }
      //    printf("<%s>\n", list_filename_obs_swpc_f107[ifile]);
      // // // skip header
      
      found_eoh = 0;
      while ( found_eoh == 0 && !feof(file_obs)) {
	getline(&line, &len, file_obs);
	sscanf(line, "%s", text);
	strcpy(text2, "");
	strncat(text2, &text[0], 4);
	if (  strcmp( "#---", text2  ) == 0 )  {
	  found_eoh = 1;
	}
      }


      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	printf("------ (load_options) (lin_interpolate_swpc_mod) Reading f10.7 file_obs.\n");
      }

      // // // Read file
      while (!feof(file_obs)){
	getline(&line, &len, file_obs);

	if (feof(file_obs)){
	  /* printf("<%s>\n",list_filename_obs_swpc_f107[ifile]); */
	  /* 	  test_print("P"); */
	  break;
	}

	//	printf("\n<%s>\n", line);
	sscanf(line, "%s %s %s %lf %s", year_obs, month_obs, day_obs, &f107_obs, str_trash);
	if ( strcmp(year_obs, ":Product:") == 0 ){ // the end of the quarter of the year file is reached so skip header of new quarter

	  found_eoh = 0;
	  while ( found_eoh == 0 && !feof(file_obs)) {
	    getline(&line, &len, file_obs);
	    sscanf(line, "%s", text);
	    strcpy(text2, "");
	    strncat(text2, &text[0], 4);
	    if (  strcmp( "#---", text2  ) == 0 )  {
	      found_eoh = 1;
	      getline(&line, &len, file_obs);
	      //    printf("\n<%s>\n", line);
	      sscanf(line, "%s %s %s %lf %s", year_obs, month_obs, day_obs, &f107_obs, str_trash);
	    }
	  }
	}
	//	  printf("XXXXXXXXXXXXXXXXX <%s>", line);
	strcpy(time_obs, year_obs);
	strcat(time_obs, "-");
	strcat(time_obs, month_obs);
	strcat(time_obs, "-");
	strcat(time_obs, day_obs);
	strcat(time_obs, "T00:00:00.000");
	str2et_c(time_obs, &et_time_obs);



	if (fabs(et_time_obs - et_initial_epoch_midnight_minus_81_days) < 0.01){// 0.01 for numerical reasons
	  start_saving_data = 1;
	}

	// // // // start saving values only at the epoch start minus 81 days
	if ( start_saving_data == 1){
	  if (feof(file_obs)){
	    break;
	  }
	 

	  // // // // stop reading file at epoch end or one day before first prediction
	  if ( ifile == nb_file_obs_swpc_f107 - 1 ){
	    if (fabs(et_time_obs - ( et_final_epoch_midnight + 2* 24 * 3600 ) ) < 0.01){// 0.01 for numerical reasons
	      break;
	    }
	    if (swpc_need_predictions == 1){
	      if ( et_time_obs > et_first_prediction  + 0.01 ){
		break;
	      }
	    }
	  }
	  //	    printf("AAAAAAA = %d %d\n", index_nb_f107, nb_elements_in_file_f107);
	  x_f107_before_interpo[index_nb_f107] = et_time_obs;
	  f107_before_interpo[index_nb_f107] = f107_obs;
	  /* etprint(x_f107_before_interpo[index_nb_f107], "" ); */
	  /*      printf("%f %d\n", f107_before_interpo[index_nb_f107], index_nb_f107); */
	  if ( index_nb_f107 > 0 ){
	    if ( f107_before_interpo[index_nb_f107] > ( missing_data_value - 0.01 ) )
	      f107_before_interpo[index_nb_f107] = f107_before_interpo[index_nb_f107-1];
	  }


	  //	printf("<%s>\n", line);
	  f107_at_least_one_obs = 1;
	  index_nb_f107 = index_nb_f107 + 1;
	}
      }
      fclose(file_obs);
      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	printf("------ (load_options) (lin_interpolate_swpc_mod) Done reading f10.7 file_obs.\n");
      }

    } // end of going over f107 observaation files


      // Ap
    index_nb_ap = 0;
    start_saving_data = 0;
    //    // Go over each ap observation file
    for (ifile = 0; ifile < nb_file_obs_swpc_ap ; ifile++){
      // // // open ap observation file
      file_obs = fopen(list_filename_obs_swpc_ap[ifile], "r");
      if (file_obs == NULL){
	printf("***! (load_options)(lin_interpolate_swpc_mod) The Ap observation file was not found. The program will stop. !***\n"); MPI_Finalize(); exit(0);
      }
      //	  printf("<%s>\n", list_filename_obs_swpc_ap[ifile]);
      // // // skip header
      
      found_eoh = 0;
      while ( found_eoh == 0 && !feof(file_obs)) {
	getline(&line, &len, file_obs);
	sscanf(line, "%s %s", text, text2);
	if (  strcmp( "Date", text2  ) == 0 )  {
	  found_eoh = 1;
	}
      }


      // // // Read file
      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	printf("------ (load_options) (lin_interpolate_swpc_mod) Reading Ap file_obs.\n");
      }

      while (!feof(file_obs)){
	getline(&line, &len, file_obs);
	if (feof(file_obs)){
	  break;
	}

	//	printf("\n<%s>\n", line);
	sscanf(line, "%s %s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", year_obs, month_obs, day_obs, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &ap_obs);
	//		printf("%s %s %s %f\n", year_obs, month_obs, day_obs, ap_obs);
	if ( strcmp(year_obs, ":Product:") == 0 ){ // the end of the quarter of the year file is reached so skip header of new quarter

	  found_eoh = 0;
	  while ( found_eoh == 0 && !feof(file_obs)) {
	    getline(&line, &len, file_obs);
	    sscanf(line, "%s %s", text, text2);
	    //	      printf("<%s>\n", text2);
	    if (  strcmp( "Date", text2  ) == 0 )  {
	      found_eoh = 1;
	      getline(&line, &len, file_obs);
	      //		  printf("\n<%s>\n", line);
	      sscanf(line, "%s %s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", year_obs, month_obs, day_obs, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &float_trash, &ap_obs);
	      //		printf("%s %s %s %f\n", year_obs, month_obs, day_obs, ap_obs);
	    }
	  }
	}
	strcpy(time_obs, year_obs);
	strcat(time_obs, "-");
	strcat(time_obs, month_obs);
	strcat(time_obs, "-");
	strcat(time_obs, day_obs);
	strcat(time_obs, "T00:00:00.000");
	str2et_c(time_obs, &et_time_obs);

	if (fabs(et_time_obs - et_initial_epoch_midnight) < 0.01){// 0.01 for numerical reasons
	  start_saving_data = 1;
	}
	    
	// // // // start saving values only at the epoch start 
	if ( start_saving_data == 1){
	  if (feof(file_obs)){
	    break;
	  }

	  // // // // stop reading file at epoch end or one day before first prediction

	  if ( ifile == nb_file_obs_swpc_ap - 1 ){
	    if (fabs(et_time_obs - ( et_final_epoch_midnight + 2* 24 * 3600 ) ) < 0.01){// 0.01 for numerical reasons
	      break;
	    }
	    if (swpc_need_predictions == 1){
	      if ( et_time_obs > et_first_prediction - 0.01 ){
		break;
	      }
	    }
	  }
	  //	    printf("AAAAAAA = %d %d\n", index_nb_f107, nb_elements_in_file_f107);
	  x_ap_before_interpo[index_nb_ap] = et_time_obs;


	  ap_before_interpo[index_nb_ap] = ap_obs;

	  if ( index_nb_ap > 0 ){
	    if ( ap_before_interpo[index_nb_ap] > ( missing_data_value - 0.01 ) ){
	      ap_before_interpo[index_nb_ap] = ap_before_interpo[index_nb_ap-1];
	    }
	  }


	   
	  //	    	    printf("<%s> !!!!!!! %f || %d-%d\n", line, ap_before_interpo[index_nb_ap], index_nb_ap, nb_elements_in_file_ap);
	  ap_at_least_one_obs = 1;
	  index_nb_ap = index_nb_ap + 1;

	}

      }
      fclose(file_obs);
      if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	printf("------ (load_options) (lin_interpolate_swpc_mod) Done reading Ap file_obs.\n");
      }

    }

    if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
      printf("------ (load_options) (lin_interpolate_swpc_mod) Done reading observation files.\n");
    }

  }


  //  Read prediction files
  if (swpc_need_predictions == 1){
    if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
      printf("------ (load_options) (lin_interpolate_swpc_mod) Reading prediction files.\n");
    }

    // // open prediction file
    file_pred = fopen(filename_f107_ap_pred, "r");
    if (file_pred == NULL){
      print_error_any_iproc(iProc, "(load_options)(lin_interpolate_swpc_mod) Could not fine filename_f107_ap_pred");
    }

    // // Ap block
    // // // skip ap header
    found_eoh = 0;
    while ( ( found_eoh == 0 ) && ( !feof(file_pred) ) ){
      getline(&line, &len, file_pred);
      sscanf(line, "%s", text);
      if (strcmp(text, "45-DAY") == 0){
	found_eoh = 1;
      }
    }


    // // // read predictions. !!! ASSUMPTION: 45 days predictions, 9 rows, 5 columns
    for (irow = 0; irow < nrow_pred; irow++){
      getline(&line, &len, file_pred);
      if (irow == 0){
      	sscanf(line, "%s %lf %s %lf %s %lf %s %lf %s %lf", date_pred1, &float_trash, date_pred2, &float_trash, date_pred3, &float_trash, date_pred4, &float_trash, date_pred5, &float_trash);

	strcpy(date_pred_temp, "20");
	strncat(date_pred_temp, &date_pred1[0] + 5, 2);
	strcat(date_pred_temp, " ");
	strncat(date_pred_temp, &date_pred1[0] + 2, 3);
	strcat(date_pred_temp, " ");
	strncat(date_pred_temp, &date_pred1[0], 2);
	strcat(date_pred_temp, " 00:00:00.000");
	str2et_c(date_pred_temp, &et_date_pred_temp);
	if (swpc_need_observations == 1){
	  // //  Sometimes there is more than one day between the first pred and the last obs (if obs have not been updated by SWPC for instance). If it is the case, add "fake" observations from the last obs to the first pred, assuming ap constant. Ex: if last obs is 04/26 (ap = 12) and first predi s 04/29 thif an add fake obs for 04/27 and 04/28, assuming a ap constant for these 2 days: 12.

	  if ( (ap_at_least_one_obs == 1) && ( fabs(et_date_pred_temp - x_ap_before_interpo[index_nb_ap-1] ) > 25*3600 ) ){
	    ifakeobs = 0;
	    x_ap_before_interpo[index_nb_ap] = x_ap_before_interpo[index_nb_ap-1] + 24*3600;
	    ap_before_interpo[index_nb_ap] = ap_before_interpo[index_nb_ap - 1];
	    //	      et2utc_c(x_ap_before_interpo[index_nb_ap], "ISOC", 3, 255, time_x_ap_before_interpo);
	    //	      printf("%s %f (%d)\n", time_x_ap_before_interpo,  ap_before_interpo[index_nb_ap],  index_nb_ap);
	    while ( fabs(et_date_pred_temp - x_ap_before_interpo[index_nb_ap] ) > 25*3600 ){
  	      index_nb_ap = index_nb_ap + 1;
  	      x_ap_before_interpo[index_nb_ap] = x_ap_before_interpo[index_nb_ap-1] + 24*3600;
  	      ap_before_interpo[index_nb_ap] = ap_before_interpo[index_nb_ap - 1];
  	      //	      et2utc_c(x_ap_before_interpo[index_nb_ap], "ISOC", 3, 255, time_x_ap_before_interpo);
  	      //	      printf("while: %s %f (%d)\n", time_x_ap_before_interpo,  ap_before_interpo[index_nb_ap],  index_nb_ap);

	    }
	    index_nb_ap = index_nb_ap + 1;
	  }
	}
 
      }

      nb_day_pred = index_nb_ap; // initalization but will be calculated after
      index_ap_first_pred = index_nb_ap;
      sscanf(line, "%s %lf %s %lf %s %lf %s %lf %s %lf", date_pred1, &ap_before_interpo_pred1, date_pred2, &ap_before_interpo_pred2, date_pred3, &ap_before_interpo_pred3, date_pred4, &ap_before_interpo_pred4, date_pred5, &ap_before_interpo_pred5);
     
      //            	printf("%s %f %s %f %s %f %s %f %s %f\n", date_pred1, ap_before_interpo_pred1, date_pred2, ap_before_interpo_pred2, date_pred3, ap_before_interpo_pred3, date_pred4, ap_before_interpo_pred4, date_pred5, ap_before_interpo_pred5);

      // // // // convert date into seconds past J2000
      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred1[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred1[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred1[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_ap_before_interpo[index_nb_ap]);
      ap_before_interpo[index_nb_ap] = ap_before_interpo_pred1;

      if (x_ap_before_interpo[index_nb_ap] >= et_final_epoch_midnight + 24 * 3600 - dt/4.){ // - dt/4. for numerical reason
	if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	  etprint(x_ap_before_interpo[index_nb_ap], "------ (load_options) (lin_interpolate_swpc_mod) Last day of prediction");
	}
	nb_day_pred = index_nb_ap + 1 ; // used to be nb_day_pred = index_nb_ap - nb_day_pred + 1 ; // nb_day_pred is equal to epoch end day - first value of prediction file + 1
	break;
      }

      //	  char date_pred_temp2[256];
      //	  et2utc_c(x_ap_before_interpo[index_nb_ap], "ISOC", 0, 255, date_pred_temp2);

      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred2[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred2[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred2[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_ap_before_interpo[index_nb_ap+1]);
      ap_before_interpo[index_nb_ap+1] = ap_before_interpo_pred2;

      if (x_ap_before_interpo[index_nb_ap+1] >= et_final_epoch_midnight + 24 * 3600 - dt/4.){// - dt/4. for numerical reason
	if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	  etprint(x_ap_before_interpo[index_nb_ap+1], "------ (load_options) (lin_interpolate_swpc_mod) Last day of prediction");
	}

	
	nb_day_pred = index_nb_ap + 1  + 1;  // used to be index_nb_ap - nb_day_pred + 1  + 1; // nb_day_pred is equal to epoch end day - first value of prediction file + 1
	break;
      }

      //	  et2utc_c(x_ap_before_interpo[index_nb_ap+1], "ISOC", 0, 255, date_pred_temp2);

      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred3[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred3[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred3[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_ap_before_interpo[index_nb_ap+2]);
      ap_before_interpo[index_nb_ap+2] = ap_before_interpo_pred3;
      //      etprint(x_ap_before_interpo[index_nb_ap+2],"x_ap_before_interpo[index_nb_ap+2]");
      if (x_ap_before_interpo[index_nb_ap+2] >= et_final_epoch_midnight + 24 * 3600 - dt/4.){ // - dt/4. for numerical reason (sometimes x_ap_before_interpo[index_nb_ap+2] = et_final_epoch_midnight + 24 * 3600 but you don't get in this if because of numerical errors)
	if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	  etprint(x_ap_before_interpo[index_nb_ap+2], "------ (load_options) (lin_interpolate_swpc_mod) Last day of prediction");
	}

	nb_day_pred =index_nb_ap + 1  + 2; // used to be index_nb_ap - nb_day_pred + 1  + 2; // nb_day_pred is equal to epoch end day - first value of prediction file + 1

	break;
      }

      //	  et2utc_c(x_ap_before_interpo[index_nb_ap+2], "ISOC", 0, 255, date_pred_temp2);


      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred4[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred4[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred4[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      
      str2et_c(date_pred_temp, &x_ap_before_interpo[index_nb_ap+3]);
/*       etprint(et_final_epoch_midnight, "et_final_epoch_midnight"); */
/*       printf("sss %f %d\n", ap_before_interpo[0], index_nb_ap);exitf(); */
/*       etprint(x_ap_before_interpo[index_nb_ap+3] , "x_ap_before_interpo[index_nb_ap+3] "); */
/*       etprint(et_final_epoch_midnight + 24 * 3600 - dt/4., "et_final_epoch_midnight + 24 * 3600 - dt/4."); */
      ap_before_interpo[index_nb_ap+3] = ap_before_interpo_pred4;
      if (x_ap_before_interpo[index_nb_ap+3] >= et_final_epoch_midnight + 24 * 3600 - dt/4.){// - dt/4. for numerical reason
	if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	  etprint(x_ap_before_interpo[index_nb_ap+3], "------ (load_options) (lin_interpolate_swpc_mod) Last day of prediction");
	}
	nb_day_pred = index_nb_ap + 1  + 3; // used to be index_nb_ap - nb_day_pred + 1  + 3; // nb_day_pred is equal to epoch end day - first value of prediction file + 1

	break;
      }

      //	  et2utc_c(x_ap_before_interpo[index_nb_ap+3], "ISOC", 0, 255, date_pred_temp2);

      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred5[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred5[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred5[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_ap_before_interpo[index_nb_ap+4]);
      ap_before_interpo[index_nb_ap+4] = ap_before_interpo_pred5;
      if (x_ap_before_interpo[index_nb_ap+4] >= et_final_epoch_midnight + 24 * 3600 - dt/4.){// - dt/4. for numerical reason
	if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
	  etprint(x_ap_before_interpo[index_nb_ap+4], "------ (load_options) (lin_interpolate_swpc_mod) Last day of prediction");
	}
	nb_day_pred = index_nb_ap  + 1 +4; // used to be index_nb_ap - nb_day_pred + 1 +4; // nb_day_pred is equal to epoch end day - first value of prediction file + 1
	break;
      }
      //	  	  et2utc_c(x_ap_before_interpo[index_nb_ap+4], "ISOC", 0, 255, date_pred_temp2);

      //	printf("%f %f %f %f %f\n", ap_before_interpo[index_nb_ap],  ap_before_interpo[index_nb_ap+1],  ap_before_interpo[index_nb_ap+2],  ap_before_interpo[index_nb_ap+3],  ap_before_interpo[index_nb_ap+4] );
      index_nb_ap = index_nb_ap + 5;  // be careful if using index_nb_ap in the rest of the script

    }
    swpc_et_last_observation = x_ap_before_interpo[index_ap_first_pred]-24*3600.;
    if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
        
    etprint(swpc_et_last_observation, "------ (load_options) (lin_interpolate_swpc_mod) Last day of observation");
      pti(nb_day_pred, "------ (load_options) (lin_interpolate_swpc_mod) Number of prediction days");
    }
    /* pti(nb_day_pred,"nb_day_pred"); */
    /* exit(0); */

    // // F107 block
    // // // skip f107 header
    found_eoh = 0;
    while ( ( found_eoh == 0 ) && ( !feof(file_pred) ) ){
      getline(&line, &len, file_pred);
      sscanf(line, "%s", text);
      if (strcmp(text, "45-DAY") == 0){
	found_eoh = 1;
      }
    }
    
    // // // read predictions. !!! ASSUMPTION: 45 days predictions, 9 rows, 5 columns
    for (irow = 0; irow < nrow_pred; irow++){
      getline(&line, &len, file_pred);


      //  See comment below
      if (irow == 0){
  	sscanf(line, "%s %lf %s %lf %s %lf %s %lf %s %lf", date_pred1, &float_trash, date_pred2, &float_trash, date_pred3, &float_trash, date_pred4, &float_trash, date_pred5, &float_trash);
	strcpy(date_pred_temp, "20");
	strncat(date_pred_temp, &date_pred1[0] + 5, 2);
	strcat(date_pred_temp, " ");
	strncat(date_pred_temp, &date_pred1[0] + 2, 3);
	strcat(date_pred_temp, " ");
	strncat(date_pred_temp, &date_pred1[0], 2);
	strcat(date_pred_temp, " 00:00:00.000");
	str2et_c(date_pred_temp, &et_date_pred_temp);
	// //  Sometimes the first pred corresponds to the last obs so in that case we overwrite the last obs (as it is often flawed because it corresponds to the current day)
	if  (swpc_need_observations == 1){
	  // //  Sometimes there is more than one day between the first pred and the last obs (if obs have not been updated by SWPC for instance). If it is the case, add "fake" observations from the last obs to the first pred, assuming F10.7 constant. Ex: if last obs is 04/26 (F10.7 = 90) and first pred is 04/29 thif an add fake obs for 04/27 and 04/28, assuming a F10.7 constant for these 2 days: 90.
	  if ( (f107_at_least_one_obs == 1) && fabs(et_date_pred_temp - x_f107_before_interpo[index_nb_f107-1] ) > 25*3600 ){
	    ifakeobs = 0;
	    x_f107_before_interpo[index_nb_f107] = x_f107_before_interpo[index_nb_f107-1] + 24*3600;
	    f107_before_interpo[index_nb_f107] = f107_before_interpo[index_nb_f107 - 1];
	    /* char time_x_f107_before_interpo[256]; */
	    /*       printf("%f %d\n", f107_before_interpo[index_nb_f107], index_nb_f107); */
	    /* et2utc_c(x_f107_before_interpo[index_nb_f107], "ISOC", 3, 255, time_x_f107_before_interpo); */
	    /*      printf("%s %f (%d)\n", time_x_f107_before_interpo,  f107_before_interpo[index_nb_f107],  index_nb_f107); */
	    while ( fabs(et_date_pred_temp - x_f107_before_interpo[index_nb_f107] ) > 25*3600 ){
	      index_nb_f107 = index_nb_f107 + 1;
	      x_f107_before_interpo[index_nb_f107] = x_f107_before_interpo[index_nb_f107-1] + 24*3600;
	      f107_before_interpo[index_nb_f107] = f107_before_interpo[index_nb_f107 - 1];
	      /* char time_x_f107_before_interpo[256]; */
	      /* et2utc_c(x_f107_before_interpo[index_nb_f107], "ISOC", 3, 255, time_x_f107_before_interpo); */
	      /* 	      printf("while: %s %f (%d)\n", time_x_f107_before_interpo,  f107_before_interpo[index_nb_f107],  index_nb_f107); */

	    }
	    index_nb_f107 = index_nb_f107 + 1;

	  }
	}
      }
      /* print_test(); */
      /* et2utc_c(x_f107_before_interpo[index_nb_f107-1], "ISOC", 3, 255, time_x_f107_before_interpo); */
      /* printf("out %s %f (%d-%d)\n", time_x_f107_before_interpo, f107_before_interpo[index_nb_f107-1], index_nb_f107, nb_elements_in_file_f107); */
      /* exit(0); */
      sscanf(line, "%s %lf %s %lf %s %lf %s %lf %s %lf", date_pred1, &f107_before_interpo_pred1, date_pred2, &f107_before_interpo_pred2, date_pred3, &f107_before_interpo_pred3, date_pred4, &f107_before_interpo_pred4, date_pred5, &f107_before_interpo_pred5);

      //printf("%s %f %s %f %s %f %s %f %s %f\n", date_pred1, f107_before_interpo[index_nb_f107], date_pred2, f107_before_interpo[index_nb_f107+1], date_pred3, f107_before_interpo[index_nb_f107+2], date_pred4, f107_before_interpo[index_nb_f107+3], date_pred5, f107_before_interpo[index_nb_f107+4]);

      // // // // convert date into seconds past J2000
      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred1[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred1[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred1[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_f107_before_interpo[index_nb_f107]);
      f107_before_interpo[index_nb_f107] = f107_before_interpo_pred1;
      /* etprint(x_f107_before_interpo[index_nb_f107], "x_f107_before_interpo[index_nb_f107]"); */
      /* printf("%f %d\n",f107_before_interpo[index_nb_f107],index_nb_f107); */

      if (x_f107_before_interpo[index_nb_f107] >= et_final_epoch_midnight + 24*3600 - dt/4.){// - dt/4. for numerical reason
	//	    print_test();
	break;
      }

      //	  char date_pred_temp2[256];
      //	  et2utc_c(x_f107_before_interpo[index_nb_f107], "ISOC", 0, 255, date_pred_temp2);

      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred2[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred2[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred2[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_f107_before_interpo[index_nb_f107+1]);
      f107_before_interpo[index_nb_f107+1] = f107_before_interpo_pred2;
      /* etprint(x_f107_before_interpo[index_nb_f107+1], "x_f107_before_interpo[index_nb_f107+1]"); */
      /* printf("%f %d\n",f107_before_interpo[index_nb_f107+1],index_nb_f107+1); */

      if (x_f107_before_interpo[index_nb_f107+1] >= et_final_epoch_midnight + 24*3600 - dt/4.){// - dt/4. for numerical reason
	//	    test_print("A");
	break;
      }

      //	  et2utc_c(x_f107_before_interpo[index_nb_f107+1], "ISOC", 0, 255, date_pred_temp2);

      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred3[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred3[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred3[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_f107_before_interpo[index_nb_f107+2]);
      f107_before_interpo[index_nb_f107+2] = f107_before_interpo_pred3;
      if (x_f107_before_interpo[index_nb_f107+2] >= et_final_epoch_midnight + 24 * 3600 - dt/4.){// - dt/4. for numerical reason
	break;
      }

      //	  et2utc_c(x_f107_before_interpo[index_nb_f107+2], "ISOC", 0, 255, date_pred_temp2);

      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred4[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred4[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred4[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_f107_before_interpo[index_nb_f107+3]);
      f107_before_interpo[index_nb_f107+3] = f107_before_interpo_pred4;
      if (x_f107_before_interpo[index_nb_f107+3] >= et_final_epoch_midnight + 24 * 3600 - dt/4.){// - dt/4. for numerical reason
	break;
      }

      //	  et2utc_c(x_f107_before_interpo[index_nb_f107+3], "ISOC", 0, 255, date_pred_temp2);

      strcpy(date_pred_temp, "20");
      strncat(date_pred_temp, &date_pred5[0] + 5, 2);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred5[0] + 2, 3);
      strcat(date_pred_temp, " ");
      strncat(date_pred_temp, &date_pred5[0], 2);
      strcat(date_pred_temp, " 00:00:00.000");
      str2et_c(date_pred_temp, &x_f107_before_interpo[index_nb_f107+4]);
      f107_before_interpo[index_nb_f107+4] = f107_before_interpo_pred5;
      if (x_f107_before_interpo[index_nb_f107+4] >= et_final_epoch_midnight + 24 * 3600 - dt/4.){// - dt/4. for numerical reason
	break;
      }

      //	  et2utc_c(x_f107_before_interpo[index_nb_f107+4], "ISOC", 0, 255, date_pred_temp2);

      index_nb_f107 = index_nb_f107 + 5; // be careful if using index_nb_f107 in the rest of the script
    }

    fclose(file_pred);
    if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
      printf("------ (load_options) (lin_interpolate_swpc) Done reading prediction files\n");
    }

  }


/*     int ooo; */
/*   for (ooo = 0; ooo < nb_elements_in_file_ap; ooo++){ */
/*     et2utc_c(x_ap_before_interpo[ooo], "ISOC", 0, 255, date_pred_temp); */
/*     printf("%s %f\n", date_pred_temp, ap_before_interpo[ooo]); */
/*   } */
  /* exit(0); */
  if ( time_step_before_interpo > 23 * 3600 ){ // recall the assumption: the time step of the F10.7/Ap file has to be either one hour, or one day, or one month. If it is one day or one month, then we have NRLMSIS00e use daily Ap. Files from SWPC are assumed to be one day time step so Ap dauly is used, not ap hist
    *use_ap_hist = 0;
  }



  ///////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// LINEAR INTERPOLATION OF F10.7 //////////////////
  ///////////////////////////////////////////////////////////////////////////////////


  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
    printf("----- (load_options) (lin_interpolate_swpc) Linear interpolating F10.7...\n");
  }


/*   //  char times[256]; */
/*   //  pti(nb_time_steps_simu, "nb_time_steps_simu"); */
/*   for (i = 0; i < nb_time_steps_simu; i++){ */
/*     x_after_interpo[i] = et_initial + dt*i / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*     //    et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times); */
/*     //   pti(i, "i"); */
/*   } */

  /* /\\* Compute y_after_interpo *\\/ */
  //  a = malloc( ( nb_time_steps_simu + (int)(( 2 * 24 * 3600. ) / dt) ) * sizeof(double) );
  //  b = malloc( ( nb_time_steps_simu + (int)(( 2 * 24 * 3600. ) / dt) ) * sizeof(double) );
  /* /\\* Linear interpolate F10.7 *\\/ */
  for (i = 1; i < nb_time_steps_simu-1; i++){
    
    //   printf("%d %d\n", i, nb_time_steps_simu-2);
    previous_index(&x_min_index, x_f107_before_interpo, x_after_interpo[i], nb_elements_in_file_f107);
    if (i == 1){
      save_first_index_in_f107_before_interpo = x_min_index;
    }
    /*  et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times);  */
    /* printf("\nx_after_interpo[%d] = %s \n", i, times);  */
    /* et2utc_c(x_f107_before_interpo[x_min_index], "C" ,3 ,255 , times);  */
    /* et2utc_c(x_f107_before_interpo[x_min_index+1], "C" ,3 ,255 , times2);  */
    /* printf("x_f107_before_interpo[x_min_index] = %s || %s (x_min_index = %d)\n",  times,  times2, x_min_index); exit(0); */
    x_min = x_f107_before_interpo[x_min_index];
    if (  x_after_interpo[i] - x_f107_before_interpo[nb_elements_in_file_f107-1]  > 0.01 ){
      f107_after_interpo[i] = f107_before_interpo[nb_elements_in_file_f107-1];
      //      printf("%f %d\n", f107_before_interpo[nb_elements_in_file_f107-1], nb_elements_in_file_f107-1);
    }
    else {

      y_min = f107_before_interpo[x_min_index];
      x_max = x_f107_before_interpo[x_min_index+1];
      y_max = f107_before_interpo[x_min_index+1];
      a[i- 1] = (y_max - y_min) / (x_max - x_min);
      b[i-1] = y_max - a[i-1]*x_max;
      //          printf("%f %d | %f\n", f107_before_interpo[x_min_index], x_min_index, f107_before_interpo[x_min_index+1]);
      f107_after_interpo[i] = a[i-1]*x_after_interpo[i] + b[i-1];
      //       printf("%f %d\n", f107_after_interpo[i], i);
    }


  }


 
  if( fabs( x_after_interpo[0] - x_f107_before_interpo[0] ) > 0.01 ) // if the first time in the driver file is different (by more than 0.01s here) from the first time of propagation then the first interpolated value of the driver is equal to the second one.// not correct but it;s ok (should not take index 0 for x_f107_before_interpo but index that corresponds to epoch initial)
    f107_after_interpo[0] = f107_after_interpo[1];
  else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    f107_after_interpo[0] = f107_before_interpo[0];
  if( fabs( x_after_interpo[nb_time_steps_simu-1] - x_f107_before_interpo[nb_elements_in_file_f107-1] ) > 0.01 ) // same comments as the two previous lines
    f107_after_interpo[nb_time_steps_simu-1] = f107_after_interpo[nb_time_steps_simu-2];
  else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    f107_after_interpo[nb_time_steps_simu-1] = f107_before_interpo[nb_elements_in_file_f107-1];
  /* /\\* Done with linear interpolating F10.7 *\\/ */
  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
    printf("----- (load_options) (lin_interpolate_swpc) Done linear interpolating F10.7.\n");
  }

  //  exit(0);  
  /* DONE WITH F10.7 //////////////////////////////// */
  

  ///////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// LINEAR INTERPOLATION OF AP /////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  if (*use_ap_hist == 0){
    if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
      printf("----- (load_options) (lin_interpolate_swpc) Linear interpolating Ap daily/monthly...\n");
    }

    // //   Linear interpolate Ap (if NRLMSIS00e uses daily Ap)
    for (i = 1; i < nb_time_steps_simu-1; i++){
      previous_index(&x_min_index, x_ap_before_interpo, x_after_interpo[i], nb_elements_in_file_ap);
/*       /\\*     et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times);  *\\/ */
/*       /\\*     printf("\nx_after_interpo[%d] = %s \n", i, times); *\\/ */
/*       /\\*     et2utc_c(x_ap_before_interpo[x_min_index], "C" ,3 ,255 , times);  *\\/ */
/*       /\\*     et2utc_c(x_ap_before_interpo[x_min_index+1], "C" ,3 ,255 , times2);  *\\/ */
/*       /\\*     printf("x_ap_before_interpo[x_min_index] = %s || %s\n",  times,  times2); *\\/ */
      x_min = x_ap_before_interpo[x_min_index];
      if (  x_after_interpo[i] - x_ap_before_interpo[nb_elements_in_file_ap-1]  > 0.01 ){
	ap_after_interpo[i] = ap_before_interpo[nb_elements_in_file_ap-1];
      }
      else {

	y_min = ap_before_interpo[x_min_index];
	x_max = x_ap_before_interpo[x_min_index+1];
	y_max = ap_before_interpo[x_min_index+1];
	a[i- 1] = (y_max - y_min) / (x_max - x_min);
	b[i-1] = y_max - a[i-1]*x_max;

	ap_after_interpo[i] = a[i-1]*x_after_interpo[i] + b[i-1];
/* 	 etprint(x_after_interpo[i], "time");  */
/* 	 printf("%f %f %d\n", y_max, y_min, x_min_index); */
      }
/*       etprint(x_after_interpo[i], "time"); */
/*           printf("%f\n", ap_after_interpo[i]); */
    }

/*     printf("nb_time_steps_simu %d\n",nb_time_steps_simu); */
/*     exitf(); */
 
    if( fabs( x_after_interpo[0] - x_ap_before_interpo[0] ) > 0.01 ) // if the first time in the driver file is different (by more than 0.01s here) from the first time of propagation then the first interpolated value of the driver is equal to the second one.
      ap_after_interpo[0] = ap_after_interpo[1];
    else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
      ap_after_interpo[0] = ap_before_interpo[0];
    if( fabs( x_after_interpo[nb_time_steps_simu-1] - x_ap_before_interpo[nb_elements_in_file_ap-1] ) > 0.01 ) // same comments as the two previous lines
      ap_after_interpo[nb_time_steps_simu-1] = ap_after_interpo[nb_time_steps_simu-2];
    else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
      ap_after_interpo[nb_time_steps_simu-1] = ap_before_interpo[nb_elements_in_file_ap-1];

    /* /\\* Done with linear interpolating Ap (if NRLMSIS00e uses daily Ap) *\\/ */
    if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
      printf("----- (load_options) (lin_interpolate_swpc) Done linear interpolating Ap daily/monthly.\n");
    }
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////// Done with Ap (if time step of Ap file is not one hour but one day or one month //////////////////


  else{

    printf("***! (load_options)(lin_interpolate_swpc) If using SWPC, the observations and predictions have to be every one day. The program will stop.\n!***"); MPI_Finalize();exit(0);
  }

  /*   /\* char tt[256]; *\/ */
  /*   /\* for (i = 0; i < nb_time_steps_simu; i++){ *\/ */
  /*   /\*   et2utc_c(x_after_interpo[i], "ISOC", 3, 255, tt); *\/ */
  /*   /\*   printf("%s: %f\n", tt, ap_after_interpo[i] ); *\/ */
    
  /*   /\* } *\/ */


  /*   /////////////////////////////////////////////////////////////////////////////////// */
  /*   ///////////////////////////////////// Done with Ap //////////////////////////////// */
  /*   /////////////////////////////////////////////////////////////////////////////////// */

  /*   /////////////////////////////////////////////////////////////////////////////////// */
  /*   ///////////////////////////////////// F10.7 AVERAGE /////////////////////////////// */
  /*   /////////////////////////////////////////////////////////////////////////////////// */

  /////// NOTE: the statements below are taken form the function lin_interpolate and are adapted to this function. Therefore, a few things might look weird. But it works correctly.
  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
    printf("----- (load_options) (lin_interpolate_swpc) Starting to calculate and interpolate F10.7 81 days average.\n");
  }
  f107A_before_interpo = malloc( (nb_elements_in_file_f107 - save_first_index_in_f107_before_interpo + 1 ) * sizeof(double) );
  x_f107A_before_interpo = malloc( (nb_elements_in_file_f107 - save_first_index_in_f107_before_interpo + 1 ) * sizeof(double) );
  
  /////////////////////// CALCULATE f107A_before_interpo /////////////////////////////
  //  /\* Initialize sum_81_days as the sum of all F10.7 values from the initial epoch - 81 days to the initial epoch  *\/
  //    char current_time[256];
  int start_index_in_f107_before_interpo_to_calculate_f107_average = save_first_index_in_f107_before_interpo;
  int stop_index_in_f107_before_interpo_to_calculate_f107_average = save_first_index_in_f107_before_interpo + 1;
  int current_index_in_f107_average = 0; //save_first_index_in_f107_before_interpo;
  int current_index_in_f107_before_interpo = save_first_index_in_f107_before_interpo;
  //  // in f107_before_interpo, go back from epoch start to as far in the past as we can until reaching 81 days
  while ( ( ( save_first_index_in_f107_before_interpo - start_index_in_f107_before_interpo_to_calculate_f107_average ) * time_step_before_interpo < 81 * 24 * 3600. ) && ( start_index_in_f107_before_interpo_to_calculate_f107_average > 0 ) ) {
    sum_81_days = sum_81_days + f107_before_interpo[start_index_in_f107_before_interpo_to_calculate_f107_average];
    start_index_in_f107_before_interpo_to_calculate_f107_average = start_index_in_f107_before_interpo_to_calculate_f107_average - 1;
  }

  //    /\* Now that we have initialized sum_81_days, we just have to update it at each time step by removing from it the first element of f107_before_interpo and adding the last element of f107_before_interpo (by first and last I mean in the window time of 81 days) *\/

  f107A_before_interpo[current_index_in_f107_average] = sum_81_days / ( stop_index_in_f107_before_interpo_to_calculate_f107_average - start_index_in_f107_before_interpo_to_calculate_f107_average - 1 ); // here: current_index_in_f107_average = 0
  x_f107A_before_interpo[current_index_in_f107_average] = x_f107_before_interpo[current_index_in_f107_before_interpo];

  /*    char current_time[256];  */
  /* et2utc_c(x_f107A_before_interpo[current_index_in_f107_average], "C" ,3 ,255 , current_time);  */
  /* printf("%s: %f | %d\n", current_time, f107A_before_interpo[current_index_in_f107_average], current_index_in_f107_average);  */

  while (current_index_in_f107_average < (nb_elements_in_file_f107 - save_first_index_in_f107_before_interpo  -1 ) ){
    current_index_in_f107_before_interpo = current_index_in_f107_before_interpo + 1;
    if ( ( start_index_in_f107_before_interpo_to_calculate_f107_average >= 0 ) && ( ( current_index_in_f107_before_interpo - start_index_in_f107_before_interpo_to_calculate_f107_average - 1 ) * time_step_before_interpo >= 81 * 24 * 3600. ) ){ // we get in here if there is at least 81 days before the current index

      start_index_in_f107_before_interpo_to_calculate_f107_average = start_index_in_f107_before_interpo_to_calculate_f107_average + 1;
      sum_81_days = sum_81_days - f107_before_interpo[start_index_in_f107_before_interpo_to_calculate_f107_average];
    }

    if ( stop_index_in_f107_before_interpo_to_calculate_f107_average + 1 < nb_elements_in_file_f107 ){
      stop_index_in_f107_before_interpo_to_calculate_f107_average = stop_index_in_f107_before_interpo_to_calculate_f107_average + 1;
      sum_81_days = sum_81_days + f107_before_interpo[stop_index_in_f107_before_interpo_to_calculate_f107_average];
    }


    current_index_in_f107_average = current_index_in_f107_average + 1;
    x_f107A_before_interpo[current_index_in_f107_average] = x_f107_before_interpo[current_index_in_f107_before_interpo];
    f107A_before_interpo[current_index_in_f107_average] = sum_81_days / ( stop_index_in_f107_before_interpo_to_calculate_f107_average - start_index_in_f107_before_interpo_to_calculate_f107_average - 1 );

  }

 
  /*      char current_time[256];  */
  /*      int ppp; */
  /*   for (ppp = 0; ppp < (nb_elements_in_file_f107 - save_first_index_in_f107_before_interpo ); ppp++){ */
  /*   et2utc_c(x_f107A_before_interpo[ppp], "C" ,3 ,255 , current_time); */
  /*   printf("%s: %f\n", current_time, f107A_before_interpo[ppp]); */
  /* } */


/*   // !!!!!!!!!!!! REMOVE BLOCK BELOW collison paper 2*/ 
/*        int ppp; */
/*     for (ppp = 0; ppp < (nb_elements_in_file_f107 - save_first_index_in_f107_before_interpo ); ppp++){ */
/*       f107A_before_interpo[ppp] = 100.; */
/*   } */
/*   // !!!!!!!!!!!! end of REMOVE BLOCK BELOW */

  /////////////////////// LINEAR INTERPOLATE f107A_before_interpo /////////////////////////////
  int nb_elements_in_file_f107A = (nb_elements_in_file_f107 - save_first_index_in_f107_before_interpo );
  //  Linear interpolate F10.7 average
  for (i = 1; i < nb_time_steps_simu-1; i++){
    previous_index(&x_min_index, x_f107A_before_interpo, x_after_interpo[i], nb_elements_in_file_f107A);
    x_min = x_f107A_before_interpo[x_min_index];
    if (  x_after_interpo[i] - x_f107A_before_interpo[nb_elements_in_file_f107A-1]  > 0.01 ){
      f107A_after_interpo[i] = f107A_before_interpo[nb_elements_in_file_f107A-1];
    }
    else {
      y_min = f107A_before_interpo[x_min_index];
      x_max = x_f107A_before_interpo[x_min_index+1];
      y_max = f107A_before_interpo[x_min_index+1];
      a[i- 1] = (y_max - y_min) / (x_max - x_min);
      b[i-1] = y_max - a[i-1]*x_max;
      f107A_after_interpo[i] = a[i-1]*x_after_interpo[i] + b[i-1];
    }

  }

  if( fabs( x_after_interpo[0] - x_f107A_before_interpo[0] ) > 0.01 ) // if the first time in the driver file is different (by more than 0.01s here) from the first time of propagation then the first interpolated value of the driver is equal to the second one.
    f107A_after_interpo[0] = f107A_after_interpo[1];
  else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    f107A_after_interpo[0] = f107A_before_interpo[0];
  if( fabs( x_after_interpo[nb_time_steps_simu-1] - x_f107A_before_interpo[nb_elements_in_file_f107A-1] ) > 0.01 ) // same comments as the two previous lines
    f107A_after_interpo[nb_time_steps_simu-1] = f107A_after_interpo[nb_elements_in_file_f107A-2];
  else // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    f107A_after_interpo[nb_time_steps_simu-1] = f107A_before_interpo[nb_elements_in_file_f107A-1];
  //  Done with linear interpolating F10.7

  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
    printf("----- (load_options) (lin_interpolate_swpc) Done calculating and interpolating F10.7 81 days average.\n");
  }

  /* char tt[256]; */
  /* for (i = 0; i < nb_time_steps_simu; i++){ */
  /*   et2utc_c(x_after_interpo[i], "ISOC", 3, 255, tt); */
  /*   printf("%s: %f\n", tt, f107A_after_interpo[i] ); */
    
  /* } */


  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////// Done with F10.7 AVERAGE /////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// UNCERTAINTY FILE ON F10.7 AND AP ////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  if (swpc_need_predictions){


    if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
      printf("----- (load_options) (lin_interpolate) Reading and interpolating uncertainties in F10.7, Ap, and F10.7A.\n");
    }




    FILE *file_uncertainty = NULL;
    char filename_uncertainty[1100];
    // Read uncerainty file
    if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
      printf("------ (load_options) (lin_interpolate) Reading uncertainties in F10.7, Ap, and F10.7A.\n");
    }


    strcpy(filename_uncertainty, filename_f107_ap_mod);
    file_uncertainty = fopen(filename_uncertainty, "r");

    if (file_uncertainty == NULL){
      printf("***! (load_options)(lin_interpolate_swpc) The file for the uncertainty on F10.7 and Ap was not found. The program will stop. !***\n"); MPI_Finalize(); exit(0);
    }
    // // Skip header
    found_eoh = 0;
    while ( (found_eoh == 0) && (!feof(file_uncertainty)) ){
      getline(&line, &len, file_uncertainty);
      sscanf(line, "%s", text);
      if (strcmp(text, "#Time(day)") == 0){
	found_eoh = 1;
      }
    }
    if (found_eoh == 0){
      printf("***! (load_options)(lin_interpolate_swpc) No data has been found in the file for the uncertainty on F10.7 and Ap was not found. The program will stop. !***\n"); MPI_Finalize(); exit(0);
    }

    int nb_day_pred_in_file_uncertainty = 0;
    // IMPORTANT: the uncertainty a midnight of the first day of pred is 0, the uncertainty at midnight of the second day of pred is the first value of the uncertainty file. And so on. Then linear interpolation.
    mod_f107_before_interpo = malloc( 46 * sizeof(double) ); // ASSUMPTIONS: uncertainty file is 45 day predictions.
    mod_ap_before_interpo = malloc( 46 * sizeof(double) ); // ASSUMPTIONS: uncertainty file is 45 day predictions
    x_mod_f107_ap_before_interpo = malloc( 46 * sizeof(double) ); // ASSUMPTIONS: uncertainty file is 45 day predictions
	
    x_mod_f107_ap_before_interpo[0] = 0; // the uncertainty a midnight of the first day of pred is 0
    mod_f107_before_interpo[0] = 0;
    mod_ap_before_interpo[0] = 0;
    int ipred=1;

    while ( (ipred < nb_day_pred + 1) && (!feof(file_uncertainty)) ){
      getline(&line, &len, file_uncertainty);
      nb_day_pred_in_file_uncertainty = nb_day_pred_in_file_uncertainty + 1;
      sscanf(line, "%lf %lf %lf", &x_mod_f107_ap_before_interpo[ipred], &mod_f107_before_interpo[ipred], &mod_ap_before_interpo[ipred]);
      x_mod_f107_ap_before_interpo[ipred] = x_mod_f107_ap_before_interpo[ipred] * 24 * 3600;

      if (feof(file_uncertainty)){
	printf("***! The number of days in the file with the uncertainty on F10.7 and Ap (%d) is smaller than the number of days in the future that the user wants to propagate for (%d). The program will stop. !***\n", nb_day_pred_in_file_uncertainty, nb_day_pred); MPI_Finalize(); exit(0);
      }


      //      printf("oooooo %f %f %f | %d\n", x_mod_f107_ap_before_interpo[ipred], mod_f107_before_interpo[ipred], mod_ap_before_interpo[ipred], iProc);
      ipred = ipred + 1;
    }
    fclose(file_uncertainty);

    if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
      printf("------ (load_options) (lin_interpolate) Done reading uncertainties in F10.7, Ap, and F10.7A.\n");
    }
    // // linear interpolate the uncertainties
    if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
      printf("------ (load_options) (lin_interpolate) Interpolating uncertainties in F10.7, Ap, and F10.7A.\n");
    }

    str2et_c(final_epoch, &et_final_epoch);

      nb_time_steps_in_prediction= ceil( ( ( et_final_epoch - et_first_prediction ) / dt * 2 ) )  + 1;
    for (i = 0; i < nb_time_steps_in_prediction; i++){
      x_mod_f107_ap_after_interpo[i] =  dt*i / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method
      //            ptd(x_mod_f107_ap_after_interpo[i], "time pred:");
    }

   //  exit(0);


    for (i = 1; i < nb_time_steps_in_prediction; i++){
      previous_index(&x_min_index, x_mod_f107_ap_before_interpo, x_mod_f107_ap_after_interpo[i], nb_day_pred+1);
      //      printf("%d %f %f %f | %d\n", x_min_index, x_mod_f107_ap_after_interpo[i], mod_f107_before_interpo[x_min_index], mod_f107_before_interpo[x_min_index+1], iProc);
      // F10.7
      x_min = x_mod_f107_ap_before_interpo[x_min_index];
      y_min = mod_f107_before_interpo[x_min_index];
      x_max = x_mod_f107_ap_before_interpo[x_min_index+1];

      y_max = mod_f107_before_interpo[x_min_index+1];

      a[i- 1] = (y_max - y_min) / (x_max - x_min);
      b[i-1] = y_max - a[i-1]*x_max;
      mod_f107_after_interpo[i] = a[i-1]*x_mod_f107_ap_after_interpo[i] + b[i-1];
      //    printf("xxx %f %f %f %f %d\n", mod_f107_after_interpo[i], x_mod_f107_ap_after_interpo[i], y_min, y_max, i);
      //      printf("xxx %f | %d\n", mod_f107_after_interpo[i], iProc);

      // Ap
      y_min = mod_ap_before_interpo[x_min_index];
      y_max = mod_ap_before_interpo[x_min_index+1];
      a[i- 1] = (y_max - y_min) / (x_max - x_min);
      b[i-1] = y_max - a[i-1]*x_max;

      mod_ap_after_interpo[i] = a[i-1]*x_mod_f107_ap_after_interpo[i] + b[i-1];

      //      printf("xxx ap %f, x %f, y1 %f, y2 %f, x2 %f, x1 %f, a %f\n", mod_ap_after_interpo[i], x_mod_f107_ap_after_interpo[i], y_min, y_max, x_max,x_min, x_min_index);
    }
    //        exit(0);
    mod_f107_after_interpo[0] = mod_f107_before_interpo[0]; // equal to 0
    mod_ap_after_interpo[0] = mod_ap_before_interpo[0]; // equal to 0
    /* mod_f107_after_interpo[nb_time_steps_in_prediction-1] = mod_f107_before_interpo[nb_day_pred]; */
    /* mod_ap_after_interpo[nb_time_steps_in_prediction-1] = mod_ap_before_interpo[nb_day_pred]; */


    free(mod_f107_before_interpo);
    free(mod_ap_before_interpo);


      

    if ((iProc == 0) & ( iDebugLevel >= 5 ) ){
      printf("------ (load_options) (lin_interpolate) Done interpolating uncertainties in F10.7, Ap, and F10.7A.\n");
    }

    if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
      printf("----- (load_options) (lin_interpolate) Done reading and interpolating uncertainties in F10.7, Ap, and F10.7A.\n");
    }


  }
  /*       /////////////////////////////////////////////////////////////////////////////////// */
  /*       //////////////////////////// DONE with UNCERTAINTY FILE ON F10.7 AND AP /////////// */
  /*       /////////////////////////////////////////////////////////////////////////////////// */

  free(a);
  free(b);
  free(x_f107_before_interpo);

  free(x_ap_before_interpo);
		  
  free(ap_before_interpo);
  free(f107_before_interpo);

/*   int ipred = 0; */
/*   for (i = 0; i < nb_time_steps_simu; i++){ */


/*     etprint(x_after_interpo[i], "time"); */
/*     printf("%f %f %f\n", f107_after_interpo[i], f107A_after_interpo[i], ap_after_interpo[i]); */

/* /\*     if (x_after_interpo[i] > *swpc_et_first_prediction ){ *\/ */
/* /\*       printf("mod %f %f %d\n", mod_f107_after_interpo[ipred],       x_mod_f107_ap_after_interpo[ipred], ipred); *\/ */
/* /\*       ipred = ipred + 1; *\/ */
/* /\*   } *\/ */
/*   } */

/*   exitf(); */
  
  if ((iProc == 0) & ( iDebugLevel >= 3 ) ){
    printf("---- (load_options) (lin_interpolate_swpc) Done linear interpolating F10.7 and Ap.\n");
  }

  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
    printf("----- (load_options) (lin_interpolate_swpc) Just got out of lin_interpolate_swpc.\n");
  }

  //      print_test();
  return 0;


}





///////////////////////////////////////////////////////////////////////////////////////// */
//
//  Name:           lin_interpolate_ap_hist
//  Purpose:        Linear interpolation of ap 
//  Assumptions:    - the Ap file has to be 1 hour time step
//                  - the first element of the Ap file starts at hour 0 of the day
//                  - the last element of the Ap file ends at hour 23 of the day
//                  - the Ap file must be 3-hour average (for example 00:00:00 15, 01:00:00 15, 02:00:00 15,  03:00:00 9, 04:00:00 9, 05:00:00 9)  
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 08/14/2015    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////

int lin_interpolate_ap_hist(double **y_after_interpo,
			    double *x_after_interpo,
			    char filename[256],
			    char src_file[256],
			    int nb_time_steps_simu,
			    char initial_epoch[256], 
			    double et_oldest_tle_epoch,
			    double dt,
			    double missing_data_value, int iDebugLevel, int iProc){

  /* Declarations */

  int index_in_daily_ap;
  int count_day = 0;
  double x_min_index_save;
  double average_eight_3hr;
  int sss;
  double lat_sat, long_sat, alt_sat, loc_time,dens_norm400, dens_norm410, dens_uncert, ss;
  int nb_data_points;
  int x_min_index;
  double *y_before_interpo = NULL;
  double *x_before_interpo = NULL;
  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  int found_eoh = 0;
  char text[256];
  int line_num;
  char yy[256], doy[256], hh[256];
  double driver_temp;
  double x_before_interpo_temp;
  double *daily_ap=NULL;
  int i;
  int nb_elements_in_file=0; 
  double et_initial;
  double *a = NULL;
  double *b = NULL;
  double x_min, y_min, x_max, y_max;

  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
    printf("----- (load_options) (lin_interpolate_ap_hist) (lin_interpolate_ap_hist) Starting lin_interpolate_ap_hist.\n");
  }


  /* Algorithm */


  fp = fopen(filename, "r");
  if (fp == NULL){
    printf("***! (load_options) (lin_interpolate_ap_hist) (lin_interpolate_ap_hist) The file %s could not be found. The program will stop. ***!\n", filename); MPI_Finalize(); exit(0);
  }

  if ( ( strcmp(src_file, "omniweb") == 0 ) || ( strcmp(src_file, "dynamic_manual") == 0 ) ){    

    found_eoh = 0;

    while ( found_eoh == 0 && !feof(fp)) {
      getline(&line, &len, fp);
      sscanf(line, "%s", text);
      if (  strcmp( "YEAR", text  ) == 0 )  {
	found_eoh = 1;
      }
      if (  strcmp( "YYYY-MM-DD", text  ) == 0 )  {
	found_eoh = 1;
      }
    }

    if (strcmp( "YEAR", text  ) == 0){
      /* Calculates the x and y arrays to interpolate */
      if ( strcmp(src_file, "dynamic_manual") == 0 ){
	nb_elements_file(&nb_elements_in_file, filename, "YEAR","#ENDOFFILE");
      }
      if ( strcmp(src_file, "omniweb") == 0 ){
	nb_elements_file(&nb_elements_in_file, filename, "YEAR","</pre><hr><HR>");
      }
      x_before_interpo = malloc( nb_elements_in_file * sizeof(double) );
      y_before_interpo = malloc( nb_elements_in_file * sizeof(double) );
      daily_ap = malloc( nb_elements_in_file/24. * sizeof(double) );
      for (line_num = 0; line_num<nb_elements_in_file; line_num++){
	getline(&line, &len, fp);
	sscanf(line, "%s %s %s %lf", yy, doy, hh,&driver_temp);
	strcpy(text,yy);
	strcat(text,"-");
	strcat(text,doy);
	strcat(text,"/");
	strcat(text,hh);
	strcat(text,":");
	strcat(text,"00");
	str2et_c (text, &x_before_interpo_temp );
	x_before_interpo[line_num] = x_before_interpo_temp;
	y_before_interpo[line_num] = driver_temp;
	if ( line_num > 0 ){
	  if ( y_before_interpo[line_num] > ( missing_data_value - 0.01 ) ) 
	    y_before_interpo[line_num] = y_before_interpo[line_num-1];
	}
	// this assumes the file is one-hour time step AND that the first element of the file starts at hour 0 of the day AND the last element of the file ends at hour 23 of the day AND must be a 3-hour average 

	if ( fmod( count_day, 24.) == 0 ){

	  daily_ap[count_day / 24] = driver_temp;
	  count_day = count_day + 1;
	  //	  printf("XXXXXXXXXXX = %s | %f | %d\n", text, driver_temp, count_day / 24);
	}
	else{
	  //	  printf("YYYYYYYYYYYYYYYYYYYYYYYYYYY %f\n",daily_ap[count_day / 24 ]);
	  daily_ap[count_day / 24] = daily_ap[count_day / 24 ] + driver_temp;
	  count_day = count_day + 1;
	}

	if ( fmod( count_day, 24.) == 0 ){ // end of the day
	  daily_ap[(count_day-1 )/ 24] = nearbyint( daily_ap[(count_day-1) / 24] / 24. ); // we round it to the closest integer because STK does that too.
	}
      }

    }

    else if (  strcmp( "YYYY-MM-DD", text  ) == 0 ) {
      /* Calculates the x and y arrays to interpolate */
      nb_elements_file(&nb_elements_in_file, filename, "YYYY-MM-DD","#ENDOFFILE");
      x_before_interpo = malloc( nb_elements_in_file * sizeof(double) );
      y_before_interpo = malloc( nb_elements_in_file * sizeof(double) );

      for (line_num = 0; line_num<nb_elements_in_file; line_num++){
	getline(&line, &len, fp);
	sscanf(line, "%20[^\n] %lf", text, &driver_temp);
	str2et_c (text, &x_before_interpo_temp );
	x_before_interpo[line_num] = x_before_interpo_temp;
	y_before_interpo[line_num] = driver_temp;
	if ( line_num > 0 ){
	  if ( y_before_interpo[line_num] > ( missing_data_value - 0.01 ) ) 
	    y_before_interpo[line_num] = y_before_interpo[line_num-1];
	}

	// this assumes the file is one-hour time step AND that the first element of the file starts at hour 0 of the day AND the last element of the file ends at hour 23 of the day AND must be a 3-hour average 

	if ( fmod( count_day, 24.) == 0 ){
	  //  printf("%d, %f\n",nb_elements_in_file,nb_elements_in_file/24.);
	  daily_ap[count_day / 24] = driver_temp;
	  count_day = count_day + 1;
	  //	  printf("XXXXXXXXXXX = %s | %f | %d\n", text, driver_temp, count_day / 24);
	}
	else{
	  //	  printf("YYYYYYYYYYYYYYYYYYYYYYYYYYY %f\n",daily_ap[count_day / 24 ]);
	  daily_ap[count_day / 24] = daily_ap[count_day / 24 ] + driver_temp;
	  count_day = count_day + 1;
	}

	if ( fmod( count_day, 24.) == 0 ){ // end of the day
	  daily_ap[(count_day-1 )/ 24] = nearbyint( daily_ap[(count_day-1) / 24] / 24. ); // we round it to the closest integer because STK does that too.
	}

      }

    } 
  }

  else if ( strcmp(src_file, "raid3") == 0 ){
    /* Calculates the x and y arrays to interpolate */
    nb_elements_file(&nb_elements_in_file, filename, "Two-digit","#ENDOFFILE");
    x_before_interpo = malloc( nb_elements_in_file * sizeof(double) );
    y_before_interpo = malloc( nb_elements_in_file * sizeof(double) );

    getline(&line, &len, fp);
    for (line_num = 0; line_num<nb_elements_in_file; line_num++){
      getline(&line, &len, fp);
      sscanf(line, "%s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", yy, doy, &ss, &lat_sat, &long_sat, &alt_sat, &loc_time, &y_before_interpo[line_num], &dens_norm400, &dens_norm410, &dens_uncert, &nb_data_points );
      //      printf("%s |  %s | %f %f %f %f %f %e %e %e %e  %d",yy, doy, ss, lat_sat, long_sat, alt_sat, loc_time, y_before_interpo[line_num], dens_norm400, dens_norm410, dens_uncert, nb_data_points );
      strcpy(text,yy);
      strcat(text,"-");
      strcat(text,doy);
      strcat(text,"/");
      strcat(text,"00");
      strcat(text,":");
      strcat(text,"00");
      str2et_c (text, &x_before_interpo_temp );
      x_before_interpo_temp = x_before_interpo_temp + ss;
      x_before_interpo[line_num] = x_before_interpo_temp;
    }
  }

  free(line);
  fclose(fp);

  a = malloc( nb_time_steps_simu * sizeof(double) );
  b = malloc( nb_time_steps_simu * sizeof(double) );
  /* Algorithm */  
  /* Calculates the x array on which the interpolation is done */
  str2et_c(initial_epoch, &et_initial);  
  char times[300];
/*         i = 0; */
/*   x_after_interpo[i] = et_oldest_tle_epoch;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*   while (x_after_interpo[i] < et_initial){ */
/*     i = i + 1; */
/*     x_after_interpo[i] = et_oldest_tle_epoch + dt*i / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*     //   etprint(x_after_interpo[i], "tle"); */
	
/*   } // leave this loop when x_after_interpo gets newer than constellation epoch (et_initial) */
/*   int j = 0; */
/*   while (i<nb_time_steps_simu){ */
/*     x_after_interpo[i] = et_initial + dt*j / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*     //            etprint(x_after_interpo[i], "const"); */

/* 	i = i +1; */
/* 	j = j+1; */
/*   } */

/*   for (i = 0; i < nb_time_steps_simu; i++){ */
/*     x_after_interpo[i] = et_initial + dt*i / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*     // et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times); */
/*   } */
  for (i = 1; i < nb_time_steps_simu-1; i++){
    previous_index(&x_min_index, x_before_interpo, x_after_interpo[i], nb_elements_in_file);
    et2utc_c(x_after_interpo[i], "C" ,6 ,255 , times);
    x_min_index_save = x_min_index;

    if ( x_after_interpo[i] - x_before_interpo[nb_elements_in_file-1]  > 0.01 ){ // if the ap file ends before the epoch time

      y_after_interpo[0][i] = y_before_interpo[nb_elements_in_file-1];
      y_after_interpo[1][i] = y_before_interpo[nb_elements_in_file-1]; // that's not really rigorous
      y_after_interpo[2][i] = y_before_interpo[nb_elements_in_file-1]; // that's not really rigorous
      y_after_interpo[3][i] = y_before_interpo[nb_elements_in_file-1]; // that's not really rigorous
      y_after_interpo[4][i] = y_before_interpo[nb_elements_in_file-1]; // that's not really rigorous
      y_after_interpo[5][i] = y_before_interpo[nb_elements_in_file-1]; // that's not really rigorous
      y_after_interpo[6][i] = y_before_interpo[nb_elements_in_file-1]; // that's not really rigorous
    }

    else {
      /*   *   1 : 3 hr AP index for current time  */

      x_min = x_before_interpo[x_min_index];

      /* y_min = y_before_interpo[x_min_index]; */
      /* x_max = x_before_interpo[x_min_index+1]; */
      /* y_max = y_before_interpo[x_min_index+1]; */
      /* a[i- 1] = (y_max - y_min) / (x_max - x_min); */
      /* b[i-1] = y_max - a[i-1]*x_max; */
      /* y_after_interpo[1][i] = a[i-1]*x_after_interpo[i] + b[i-1]; */
      y_after_interpo[1][i] = y_before_interpo[x_min_index]; // I commented the 6 lines above because actually the ap shold be a step functino and so this interpolation was not correct

      /* *   2 : 3 hr AP index for 3 hrs before current time */
      x_min_index = x_min_index-3; // This assumes that the Ap file has to be 1 hour time step
      /* x_min = x_before_interpo[x_min_index]; */
      /* y_min = y_before_interpo[x_min_index]; */
      /* x_max = x_before_interpo[x_min_index+1]; */
      /* y_max = y_before_interpo[x_min_index+1]; */
      /* a[i- 1] = (y_max - y_min) / (x_max - x_min); */
      /* b[i-1] = y_max - a[i-1]*x_max; */
      /* y_after_interpo[ 2][i] = a[i-1]*(x_after_interpo[i]-3*3600.) + b[i-1]; */
      y_after_interpo[ 2][i] = y_before_interpo[x_min_index];// I commented the 6 lines above because actually the ap shold be a step functino and so this interpolation was not correct
      /* *   3 : 3 hr AP index for 6 hrs before current time */
      x_min_index = x_min_index-3; // This assumes that the Ap file has to be 1 hour time step
      /* x_min = x_before_interpo[x_min_index]; */
      /* y_min = y_before_interpo[x_min_index]; */
      /* x_max = x_before_interpo[x_min_index+1]; */
      /* y_max = y_before_interpo[x_min_index+1]; */
      /* a[i- 1] = (y_max - y_min) / (x_max - x_min); */
      /* b[i-1] = y_max - a[i-1]*x_max; */
      /* y_after_interpo[3][i] = a[i-1]*(x_after_interpo[i]-6*3600.) + b[i-1]; */
      y_after_interpo[3][i] =  y_before_interpo[x_min_index];// I commented the 6 lines above because actually the ap shold be a step functino and so this interpolation was not correct
      /* *   4 : 3 hr AP index for 9 hrs before current time */
      x_min_index = x_min_index-3; // This assumes that the Ap file has to be 1 hour time step
      x_min = x_before_interpo[x_min_index];
      y_min = y_before_interpo[x_min_index];
      x_max = x_before_interpo[x_min_index+1];
      y_max = y_before_interpo[x_min_index+1];
      a[i- 1] = (y_max - y_min) / (x_max - x_min);
      b[i-1] = y_max - a[i-1]*x_max;
      y_after_interpo[ 4][i] = a[i-1]*(x_after_interpo[i]-9*3600.) + b[i-1];
      y_after_interpo[ 4][i] = y_before_interpo[x_min_index];// I commented the 6 lines above because actually the ap shold be a step functino and so this interpolation was not correct
      /* *   5 : Average of eight 3 hr AP indicies from 12 to 33 hrs  */
      /* *           prior to current time  */
      average_eight_3hr = 0.0;
      for (sss = 0; sss < 8; sss++){
	x_min_index = x_min_index-3; // This assumes that the Ap file has to be 1 hour time step
	y_min = y_before_interpo[x_min_index];
	average_eight_3hr = average_eight_3hr + y_min;
      }
      average_eight_3hr = average_eight_3hr / 8.;
      y_after_interpo[5][i] = nearbyint( average_eight_3hr ); // no interpolation here (not sure it makes a lot more sense since we're already averaging over a bunch of indices)// we round it to the closest integer because STK does that too.
      /* *           prior to current time */
      /* *   6 : Average of eight 3 hr AP indicies from 36 to 57 hrs  */
      /* *           prior to current time  */
      average_eight_3hr = 0.0;
      for (sss = 0; sss < 8; sss++){
	x_min_index = x_min_index-3; // This assumes that the Ap file has to be 1 hour time step
	y_min = y_before_interpo[x_min_index];
	average_eight_3hr = average_eight_3hr + y_min;
      }
      average_eight_3hr = average_eight_3hr / 8.;
      y_after_interpo[ 6][i] = nearbyint( average_eight_3hr ) ; // no int
      /*   *   0 : daily AP */
      // BLOCK: THAT IS IF BY 'DAILY' IT MEANS THE AVERAGE OVER THE 24 HOURS OF THE SAME DAY (AND NOT OVER THE LAST 24 HOURS)

      index_in_daily_ap = (int) ( ( x_after_interpo[i] - x_before_interpo[0] ) / ( 3600. * 24) ); // this is the number of days between the current index and the first day of the file
      y_after_interpo[0][i] = daily_ap[index_in_daily_ap];  
      // END OF BLOCK: THAT IS IF BY 'DAILY' IT MEANS THE AVERAGE OVER THE 24 HOURS OF THE SAME DAY (AND NOT OVER THE LAST 24 HOURS)

      // BLOCK: THAT IS IF BY 'DAILY' IT MEANS THE AVERAGE OVER THE LAST 24 HOURS (AND NOT OVER THE 24 HOURS OF THE SAME DAY)
      /* // compute the average ap of the day, called daily ap */
      /* average_eight_3hr = 0.0;  */
      /* for (sss = 0; sss < 8; sss++){ */
      /* 	if (sss == 0){ */
      /* 	  x_min_index = x_min_index_save - 3; // This assumes that the Ap file has to be 1 hour time step */
      /* 	} */
      /* 	else{ */
      /* 	  x_min_index = x_min_index -3; */
      /* 	} */
      /* 	et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times); */

      /* 	//printf("%s | %f\n",times, y_after_interpo[0][i]); */

      /* 	y_min = y_before_interpo[x_min_index]; */
      /* 	average_eight_3hr = average_eight_3hr + y_min; */
      /* } */
      /* average_eight_3hr = average_eight_3hr / 8.; */
      /* y_after_interpo[0][i] = average_eight_3hr;  */
      // END OF BLOCK: THAT IS IF BY 'DAILY' IT MEANS THE AVERAGE OVER THE LAST 24 HOURS (AND NOT OVER THE 24 HOURS OF THE SAME DAY)

    }

/*     etprint(x_after_interpo[i], ""); */
/*     printf("%f %f %f %f %f %f %f %f\n", y_after_interpo[0][i], y_after_interpo[1][i],y_after_interpo[2][i],y_after_interpo[3][i],y_after_interpo[4][i],y_after_interpo[5][i],y_after_interpo[6][i]); */
  }

  // Below is pretty bad but it just says that for the first time step the ap_hist is constant (equal to the one at the second time step). That's not a huge approximation anyway 
  y_after_interpo[0][0] = y_after_interpo[0][1];
  y_after_interpo[1][0] = y_after_interpo[1][1]; 
  y_after_interpo[2][0] = y_after_interpo[2][1]; 
  y_after_interpo[3][0] = y_after_interpo[3][1]; 
  y_after_interpo[4][0] = y_after_interpo[4][1]; 
  y_after_interpo[5][0] = y_after_interpo[5][1]; 
  y_after_interpo[6][0] = y_after_interpo[6][1]; 
  // same comment but for the last time step
  y_after_interpo[0][nb_time_steps_simu-1] = y_after_interpo[0][nb_time_steps_simu-2];
  y_after_interpo[1][nb_time_steps_simu-1] = y_after_interpo[1][nb_time_steps_simu-2]; 
  y_after_interpo[2][nb_time_steps_simu-1] = y_after_interpo[2][nb_time_steps_simu-2]; 
  y_after_interpo[3][nb_time_steps_simu-1] = y_after_interpo[3][nb_time_steps_simu-2]; 
  y_after_interpo[4][nb_time_steps_simu-1] = y_after_interpo[4][nb_time_steps_simu-2]; 
  y_after_interpo[5][nb_time_steps_simu-1] = y_after_interpo[5][nb_time_steps_simu-2]; 
  y_after_interpo[6][nb_time_steps_simu-1] = y_after_interpo[6][nb_time_steps_simu-2]; 

  free(a); 
  free(b); 
  free(x_before_interpo);
  free(y_before_interpo);


  if ((iProc == 0) & ( iDebugLevel >= 4 ) ){
    printf("----- (load_options) (lin_interpolate_ap_hist) (lin_interpolate_ap_hist) Ended lin_interpolate_ap_hist.\n");
  }


  return 0;

}



/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           lin_interpolate_attitude
//  Purpose:        Linear interpolation of the attitude file
//  Assumptions:    if the file gives quaternions, there should not be more tahn 5 space character on the first line of the data (this is how SpOCK differentiates it from a pitch roll yaw order-picth order_roll order_yaw file)
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 10/25/2015    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////

int lin_interpolate_attitude(double **quaternion_after_interpo,
			     double *pitch_after_interpo,
			     double *roll_after_interpo,
			     double *yaw_after_interpo,
			     int *pitch_order_after_interpo,
			     int *roll_order_after_interpo,
			     int *yaw_order_after_interpo,
			     double *x_after_interpo,
			     char filename[256],
			     int nb_time_steps_simu,
			     char initial_epoch[256], 
			     double et_oldest_tle_epoch,
			     double dt,
			     int iProc,
			     int *file_is_quaternion,
			     int use_kalman){

  /* Declarations */
  char read_time[10]; char len_fake_text[10];    char read_all_line[100];
      double **quaternion_before_interpo=NULL;
  char fake_text[256];
  int x_min_index;
  double *pitch_before_interpo = NULL;
  double *roll_before_interpo = NULL;
  double *yaw_before_interpo = NULL;
  double *pitch_order_before_interpo = NULL;
  double *roll_order_before_interpo = NULL;
  double *yaw_order_before_interpo = NULL;
  double *x_before_interpo = NULL;
  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  int found_eoh = 0;
  char text[256];
  int line_num;
  double x_before_interpo_temp;
  int i;
  int nb_elements_in_file=0; 
  double et_initial;
  double *a = NULL;
  double *b = NULL;
  double x_min, x_max;
  double y_min_quaternion[4], y_max_quaternion[4];
  double y_min_pitch, y_max_pitch;
  double y_min_roll, y_max_roll;
  double y_min_yaw, y_max_yaw;

  /* Algorithm */

  str2et_c(initial_epoch, &et_initial);  
  /* Calculates the x and y arrays to interpolate */
  nb_elements_file(&nb_elements_in_file, filename, "#ENDOFHEADER","#ENDOFFILE");



  x_before_interpo = malloc( nb_elements_in_file * sizeof(double) );
  a = malloc( nb_time_steps_simu * sizeof(double) );
  b = malloc( nb_time_steps_simu * sizeof(double) );

  fp = fopen(filename, "r");
  
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#ENDOFHEADER", text  ) == 0 )  {
      found_eoh = 1;
    }
  }

  int found_start = 0; 
  ssize_t read;
      while ( (found_start == 0) && (read = getline(&line, &len, fp)) != -1 ) {
	if (strchr(&line[0],'2') != NULL)
	  found_start = 1;
      }

    sscanf(line, "%s", text);

     int count_bracket=0;
       int  icount_bracket;

      for (icount_bracket = 0; icount_bracket< strlen(line); icount_bracket++){
	if (line[icount_bracket]=='('){
	  count_bracket = count_bracket + 1;
	  }
      }


      //    for (count_bracket=0; text[count_bracket]; text[count_bracket]==' ' ? count_bracket++ : *text++);
    // if more than 5 spaces then the file is pitch roll yaw. otherwise it's quaternion

    if (count_bracket  >= 1){
      *file_is_quaternion = 0;
    }
    else{
      *file_is_quaternion = 1;
    }

    if (*file_is_quaternion == 0) {
  pitch_before_interpo = malloc( nb_elements_in_file * sizeof(double) );
  roll_before_interpo = malloc( nb_elements_in_file * sizeof(double) );
  yaw_before_interpo = malloc( nb_elements_in_file * sizeof(double) );
  pitch_order_before_interpo = malloc( nb_elements_in_file * sizeof(double) );
  roll_order_before_interpo = malloc( nb_elements_in_file * sizeof(double) );
  yaw_order_before_interpo = malloc( nb_elements_in_file * sizeof(double) );

  for (line_num = 0; line_num<nb_elements_in_file; line_num++){
    if (line_num > 0){
    getline(&line, &len, fp);
    }
    strcpy(text,"");
    sscanf(line, "%s", fake_text);
    strcpy(read_time, "%");
    sprintf(len_fake_text, "%lu", strlen(fake_text));
    strcat(read_time, len_fake_text);
    strcat(read_time, "[^\\n] ");
    strcpy(read_all_line, read_time);
    strcat(read_all_line, "(%lf;%lf;%lf) (%lf;%lf;%lf)");
    RemoveSpaces(line);
    sscanf(line, read_all_line, text, &pitch_before_interpo[line_num], &roll_before_interpo[line_num], &yaw_before_interpo[line_num], &pitch_order_before_interpo[line_num], &roll_order_before_interpo[line_num], &yaw_order_before_interpo[line_num]);  
    

    //    sscanf(line, "%17[^\n] (%lf;%lf;%lf) (%lf;%lf;%lf)", text, &pitch_before_interpo[line_num], &roll_before_interpo[line_num], &yaw_before_interpo[line_num], &pitch_order_before_interpo[line_num], &roll_order_before_interpo[line_num], &yaw_order_before_interpo[line_num]);  

    /* print_test(); */
/*     printf("<%s> || (%f %f %f) (%f %f %f)\n", text, pitch_before_interpo[line_num], roll_before_interpo[line_num], yaw_before_interpo[line_num], pitch_order_before_interpo[line_num], roll_order_before_interpo[line_num],  yaw_order_before_interpo[line_num]); */
    /*   exit(0); */

    str2et_c (text, &x_before_interpo_temp );

    x_before_interpo[line_num] = x_before_interpo_temp;
  }

  free(line);
  fclose(fp);


  /* Calculates the x array on which the interpolation is done */

/*   //char times[256]; */
  if (use_kalman != 1){
    i = 0;
  x_after_interpo[i] = et_oldest_tle_epoch;  // "/ 2.0" because of the Runge Kunta orfer 4 method
  while (x_after_interpo[i] < et_initial){
    i = i + 1;
    x_after_interpo[i] = et_oldest_tle_epoch + dt*i / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    //        etprint(x_after_interpo[i], "tle");

  } // leave this loop when x_after_interpo gets newer than constellation epoch (et_initial)
  int j = 0;
  //   printf("Nd %d\n",nb_time_steps_simu);exit(0);
  while (i<nb_time_steps_simu){
       x_after_interpo[i] = et_initial + dt*j / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method
       //          etprint(x_after_interpo[i], "const");
	i = i +1;
	j = j+1;
  }
  }


  /* for (i = 0; i < nb_time_steps_simu; i++){ */
  /*   x_after_interpo[i] = et_initial + dt*i / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */

  /*   //    et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times); */
  /* } */

  /* Compute y_after_interpo */

  //  FILE *fp_temp;
  //  fp_temp = fopen("att-test.txt", "w+");
  //   char times[256]; 
  /*   char times2[256]; */
  for (i = 1; i < nb_time_steps_simu-1; i++){
    if ( ( x_after_interpo[i] >= x_before_interpo[0] ) && (x_after_interpo[i] <= x_before_interpo[nb_elements_in_file-1] ) ) {
      previous_index(&x_min_index, x_before_interpo, x_after_interpo[i], nb_elements_in_file);
      /*     et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times);  */
      /*     printf("\nx_after_interpo[%d] = %s \n", i, times); */
      /*     et2utc_c(x_before_interpo[x_min_index], "C" ,3 ,255 , times);  */
      /*     et2utc_c(x_before_interpo[x_min_index+1], "C" ,3 ,255 , times2);  */
      /*     printf("x_before_interpo[x_min_index] = %s || %s\n",  times,  times2); */
      x_min = x_before_interpo[x_min_index];
      x_max = x_before_interpo[x_min_index+1];
      // pitch
      y_min_pitch = pitch_before_interpo[x_min_index];
      y_max_pitch = pitch_before_interpo[x_min_index+1];
      a[i- 1] = (y_max_pitch - y_min_pitch) / (x_max - x_min);
      b[i-1] = y_max_pitch - a[i-1]*x_max;
      pitch_after_interpo[i] = a[i-1]*x_after_interpo[i] + b[i-1];
      pitch_order_after_interpo[i] = pitch_order_before_interpo[x_min_index];
      // roll
      y_min_roll = roll_before_interpo[x_min_index];
      y_max_roll = roll_before_interpo[x_min_index+1];
      a[i- 1] = (y_max_roll - y_min_roll) / (x_max - x_min);
      b[i-1] = y_max_roll - a[i-1]*x_max;
      roll_after_interpo[i] = a[i-1]*x_after_interpo[i] + b[i-1];
      roll_order_after_interpo[i] = roll_order_before_interpo[x_min_index];
      // yaw
      y_min_yaw = yaw_before_interpo[x_min_index];
      y_max_yaw = yaw_before_interpo[x_min_index+1];
      a[i- 1] = (y_max_yaw - y_min_yaw) / (x_max - x_min);
      b[i-1] = y_max_yaw - a[i-1]*x_max;
      yaw_after_interpo[i] = a[i-1]*x_after_interpo[i] + b[i-1];
      yaw_order_after_interpo[i] = yaw_order_before_interpo[x_min_index];
    }
    else if ( x_after_interpo[i] < x_before_interpo[0] ){
      pitch_after_interpo[i] = pitch_before_interpo[0];
      roll_after_interpo[i] = roll_before_interpo[0];
      yaw_after_interpo[i] = yaw_before_interpo[0];
      pitch_order_after_interpo[i] = pitch_order_before_interpo[0];
      roll_order_after_interpo[i] = roll_order_before_interpo[0];
      yaw_order_after_interpo[i] = yaw_order_before_interpo[0];
    }

    else{
      pitch_after_interpo[i] = pitch_before_interpo[nb_elements_in_file-1];
      roll_after_interpo[i] = roll_before_interpo[nb_elements_in_file-1];
      yaw_after_interpo[i] = yaw_before_interpo[nb_elements_in_file-1];
      pitch_order_after_interpo[i] = pitch_order_before_interpo[nb_elements_in_file-1];
      roll_order_after_interpo[i] = roll_order_before_interpo[nb_elements_in_file-1];
      yaw_order_after_interpo[i] = yaw_order_before_interpo[nb_elements_in_file-1];
    }
/*     printf("\n"); */
/*     etprint(x_after_interpo[i], ""); */
/* printf("%f %f %f %d %d %d\n", pitch_after_interpo[i], roll_after_interpo[i], yaw_after_interpo[i], pitch_order_after_interpo[i], roll_order_after_interpo[i], yaw_order_after_interpo[i]); */
    //    et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times);  
    //    fprintf(fp_temp, "%s %f %f %f %d %d %d\n", times,  pitch_after_interpo[i], roll_after_interpo[i], yaw_after_interpo[i], pitch_order_after_interpo[i], roll_order_after_interpo[i], yaw_order_after_interpo[i]);
  }

  if( fabs( x_after_interpo[0] - x_before_interpo[0] ) > 0.01 ){ // if the first time in the driver file is different (by more than 0.01s here) from the first time of propagation then the first interpolated value of the driver is equal to the second one. 
    pitch_after_interpo[0] = pitch_after_interpo[1];
    roll_after_interpo[0] = roll_after_interpo[1];
    yaw_after_interpo[0] = yaw_after_interpo[1];
    pitch_order_after_interpo[0] = pitch_order_after_interpo[1];
    roll_order_after_interpo[0] = roll_order_after_interpo[1];
    yaw_order_after_interpo[0] = yaw_order_after_interpo[1];

  }
  else{ // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    pitch_after_interpo[0] = pitch_before_interpo[0];
    roll_after_interpo[0] = roll_before_interpo[0];
    yaw_after_interpo[0] = yaw_before_interpo[0];
    pitch_order_after_interpo[0] = pitch_order_before_interpo[0];
    roll_order_after_interpo[0] = roll_order_before_interpo[0];
    yaw_order_after_interpo[0] = yaw_order_before_interpo[0];
  }

  if( fabs( x_after_interpo[nb_time_steps_simu-1] - x_before_interpo[nb_elements_in_file-1] ) > 0.01 ){ // same comments as the two previous blocks
    pitch_after_interpo[nb_time_steps_simu-1] = pitch_after_interpo[nb_time_steps_simu-2];
    roll_after_interpo[nb_time_steps_simu-1] = roll_after_interpo[nb_time_steps_simu-2];
    yaw_after_interpo[nb_time_steps_simu-1] = yaw_after_interpo[nb_time_steps_simu-2];
    pitch_order_after_interpo[nb_time_steps_simu-1] = pitch_order_after_interpo[nb_time_steps_simu-2];
    roll_order_after_interpo[nb_time_steps_simu-1] = roll_order_after_interpo[nb_time_steps_simu-2];
    yaw_order_after_interpo[nb_time_steps_simu-1] = yaw_order_after_interpo[nb_time_steps_simu-2];

  }
  else{ // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    pitch_after_interpo[nb_time_steps_simu-1] = pitch_before_interpo[nb_elements_in_file-1];
    roll_after_interpo[nb_time_steps_simu-1] = roll_before_interpo[nb_elements_in_file-1];
    yaw_after_interpo[nb_time_steps_simu-1] = yaw_before_interpo[nb_elements_in_file-1];
    pitch_order_after_interpo[nb_time_steps_simu-1] = pitch_order_before_interpo[nb_elements_in_file-1];
    roll_order_after_interpo[nb_time_steps_simu-1] = roll_order_before_interpo[nb_elements_in_file-1];
    yaw_order_after_interpo[nb_time_steps_simu-1] = yaw_order_before_interpo[nb_elements_in_file-1];
  }




  free(pitch_before_interpo);
  free(roll_before_interpo);
  free(yaw_before_interpo);
  free(pitch_order_before_interpo);
  free(roll_order_before_interpo);
  free(yaw_order_before_interpo);
    }


    else{

      quaternion_before_interpo = malloc( nb_elements_in_file * sizeof(double*) );

      for (line_num = 0; line_num<nb_elements_in_file; line_num++){
	quaternion_before_interpo[line_num] = malloc( 4 * sizeof(double) );
	if (line_num > 0){
	  getline(&line, &len, fp);
	}

	sscanf(line, "%s %lf %lf %lf %lf", text,&quaternion_before_interpo[line_num][0], &quaternion_before_interpo[line_num][1], &quaternion_before_interpo[line_num][2], &quaternion_before_interpo[line_num][3]);
      

      str2et_c(text, &x_before_interpo[line_num]);  
      
      //etprint(x_before_interpo[line_num], "");
		 //      printf("%f %f %f %f\n",quaternion_before_interpo[line_num][0], quaternion_before_interpo[line_num][1], quaternion_before_interpo[line_num][2], quaternion_before_interpo[line_num][3]);
/*       printf("%d %d\n", line_num, nb_elements_in_file); */
    }


      //      exit(0);


  /* Calculates the x array on which the interpolation is done */

/*   //char times[256]; */
  if (use_kalman != 1){
    i = 0;
  x_after_interpo[i] = et_oldest_tle_epoch;  // "/ 2.0" because of the Runge Kunta orfer 4 method
  while (x_after_interpo[i] < et_initial){
    i = i + 1;
    x_after_interpo[i] = et_oldest_tle_epoch + dt*i / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method
    //        etprint(x_after_interpo[i], "tle");

  } // leave this loop when x_after_interpo gets newer than constellation epoch (et_initial)
  int j = 0;
  //   printf("Nd %d\n",nb_time_steps_simu);exit(0);
  while (i<nb_time_steps_simu){
       x_after_interpo[i] = et_initial + dt*j / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method
       //          etprint(x_after_interpo[i], "const");
	i = i +1;
	j = j+1;
  }
  }



/*   i = 0; */
/*   x_after_interpo[i] = et_oldest_tle_epoch;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*   while (x_after_interpo[i] < et_initial){ */
/*     i = i + 1; */
/*     x_after_interpo[i] = et_oldest_tle_epoch + dt*i / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*     //        etprint(x_after_interpo[i], "tle"); */
	
/*   } // leave this loop when x_after_interpo gets newer than constellation epoch (et_initial) */
/*   int j = 0; */
/*   while (i<nb_time_steps_simu){ */
/*     x_after_interpo[i] = et_initial + dt*j / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*     //        etprint(x_after_interpo[i], "const"); */
/* 	i = i +1; */
/* 	j = j+1; */
/*   } */


/*   for (i = 0; i < nb_time_steps_simu; i++){ */
/*     x_after_interpo[i] = et_initial + dt*i / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */

/*     //    et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times); */
/*   } */


  for (i = 1; i < nb_time_steps_simu-1; i++){
    if ( ( x_after_interpo[i] >= x_before_interpo[0] ) && (x_after_interpo[i] <= x_before_interpo[nb_elements_in_file-1] ) ) {
      previous_index(&x_min_index, x_before_interpo, x_after_interpo[i], nb_elements_in_file);
      /*     et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times);  */
      /*     printf("\nx_after_interpo[%d] = %s \n", i, times); */
      /*     et2utc_c(x_before_interpo[x_min_index], "C" ,3 ,255 , times);  */
      /*     et2utc_c(x_before_interpo[x_min_index+1], "C" ,3 ,255 , times2);  */
      /*     printf("x_before_interpo[x_min_index] = %s || %s\n",  times,  times2); */
      x_min = x_before_interpo[x_min_index];
      x_max = x_before_interpo[x_min_index+1];
      // quaternion
      y_min_quaternion[0] = quaternion_before_interpo[x_min_index][0];
      y_min_quaternion[1] = quaternion_before_interpo[x_min_index][1];
      y_min_quaternion[2] = quaternion_before_interpo[x_min_index][2];
      y_min_quaternion[3] = quaternion_before_interpo[x_min_index][3];

      y_max_quaternion[0] = quaternion_before_interpo[x_min_index+1][0];
      y_max_quaternion[1] = quaternion_before_interpo[x_min_index+1][1];
      y_max_quaternion[2] = quaternion_before_interpo[x_min_index+1][2];
      y_max_quaternion[3] = quaternion_before_interpo[x_min_index+1][3];


      a[i- 1] = (y_max_quaternion[0] - y_min_quaternion[0]) / (x_max - x_min);
      b[i-1] = y_max_quaternion[0] - a[i-1]*x_max;
      quaternion_after_interpo[i][0] = a[i-1]*x_after_interpo[i] + b[i-1];

      a[i- 1] = (y_max_quaternion[1] - y_min_quaternion[1]) / (x_max - x_min);
      b[i-1] = y_max_quaternion[1] - a[i-1]*x_max;
      quaternion_after_interpo[i][1] = a[i-1]*x_after_interpo[i] + b[i-1];

      a[i- 1] = (y_max_quaternion[2] - y_min_quaternion[2]) / (x_max - x_min);
      b[i-1] = y_max_quaternion[2] - a[i-1]*x_max;
      quaternion_after_interpo[i][2] = a[i-1]*x_after_interpo[i] + b[i-1];

      a[i- 1] = (y_max_quaternion[3] - y_min_quaternion[3]) / (x_max - x_min);
      b[i-1] = y_max_quaternion[3] - a[i-1]*x_max;
      quaternion_after_interpo[i][3] = a[i-1]*x_after_interpo[i] + b[i-1];

/*       etprint(x_after_interpo[i],""); */
/*       printf("%f %f %f %f\n\n", quaternion_after_interpo[i][0], quaternion_after_interpo[i][1], quaternion_after_interpo[i][2], quaternion_after_interpo[i][3]); */
    }
    else if ( x_after_interpo[i] < x_before_interpo[0] ){
      quaternion_after_interpo[i][0] = quaternion_before_interpo[0][0];
      quaternion_after_interpo[i][1] = quaternion_before_interpo[0][1];
      quaternion_after_interpo[i][2] = quaternion_before_interpo[0][2];
      quaternion_after_interpo[i][3] = quaternion_before_interpo[0][3];
    }

    else{
      quaternion_after_interpo[i][0] = quaternion_before_interpo[nb_elements_in_file-1][0];
      quaternion_after_interpo[i][1] = quaternion_before_interpo[nb_elements_in_file-1][1];
      quaternion_after_interpo[i][2] = quaternion_before_interpo[nb_elements_in_file-1][2];
      quaternion_after_interpo[i][3] = quaternion_before_interpo[nb_elements_in_file-1][3];
    }


       //       fprintf(fp_temp, "%s %f %f %f %d %d %d\n", times,  quaternion_after_interpo[i], roll_after_interpo[i], yaw_after_interpo[i], quaternion_order_after_interpo[i], roll_order_after_interpo[i], yaw_order_after_interpo[i]);
  }
  //  exit(0);
  if( fabs( x_after_interpo[0] - x_before_interpo[0] ) > 0.01 ){ // if the first time in the driver file is different (by more than 0.01s here) from the first time of propagation then the first interpolated value of the driver is equal to the second one. 
    quaternion_after_interpo[0][0] = quaternion_after_interpo[1][0];
    quaternion_after_interpo[0][1] = quaternion_after_interpo[1][1];
    quaternion_after_interpo[0][2] = quaternion_after_interpo[1][2];
    quaternion_after_interpo[0][3] = quaternion_after_interpo[1][3];
  }
  else{ // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    quaternion_after_interpo[0][0] = quaternion_before_interpo[0][0];
    quaternion_after_interpo[0][1] = quaternion_before_interpo[0][1];
    quaternion_after_interpo[0][2] = quaternion_before_interpo[0][2];
    quaternion_after_interpo[0][3] = quaternion_before_interpo[0][3];
  }

  if( fabs( x_after_interpo[nb_time_steps_simu-1] - x_before_interpo[nb_elements_in_file-1] ) > 0.01 ){ // same comments as the two previous blocks
    quaternion_after_interpo[nb_time_steps_simu-1][0] = quaternion_after_interpo[nb_time_steps_simu-2][0];
    quaternion_after_interpo[nb_time_steps_simu-1][1] = quaternion_after_interpo[nb_time_steps_simu-2][1];
    quaternion_after_interpo[nb_time_steps_simu-1][2] = quaternion_after_interpo[nb_time_steps_simu-2][2];
    quaternion_after_interpo[nb_time_steps_simu-1][3] = quaternion_after_interpo[nb_time_steps_simu-2][3];
  }
  else{ // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    quaternion_after_interpo[nb_time_steps_simu-1][0] = quaternion_before_interpo[nb_elements_in_file-1][0];
    quaternion_after_interpo[nb_time_steps_simu-1][1] = quaternion_before_interpo[nb_elements_in_file-1][1];
    quaternion_after_interpo[nb_time_steps_simu-1][2] = quaternion_before_interpo[nb_elements_in_file-1][2];
    quaternion_after_interpo[nb_time_steps_simu-1][3] = quaternion_before_interpo[nb_elements_in_file-1][3];
  }





      free(line);
      fclose(fp);
    }
    
    free(a);
    free(b);

    free(quaternion_before_interpo);
  free(x_before_interpo);

  //   fclose(fp_temp);



  return 0;

}




/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           nb_elements_file
//  Purpose:        Returns the number of lines after the header. header_end represents the last line of the file before nb_elements_in_file is incremented
//  Assumptions:    None
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 08/17/2015    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////


int nb_elements_file(int *nb_element_in_file,
		     char filename[256],
		     char header_end[256],
		     char file_end[256]){

  /* Declarations */
  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  int found_eoh = 0;
  char text[256];
  
  /* Algorithm */
  fp = fopen(filename, "r");
  // Count number of elements in file
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( text, file_end  ) == 0 )  {
      found_eoh = 1;
    }
        *nb_element_in_file = *nb_element_in_file+1;
  }
  rewind(fp);
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( header_end, text  ) == 0 )  {
      found_eoh = 1;
      *nb_element_in_file = *nb_element_in_file-1;
    }
    *nb_element_in_file = *nb_element_in_file-1;
  }

  return 0;
  

}


/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           nb_time_steps
//  Purpose:        Computes the number of time steps in the simulation.
//  Assumptions:    None
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 08/18/2015    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////


int nb_time_steps(int *nb_time_steps_simu,
		  double et_initial_epoch,
		  char final_epoch[256],
		  double dt){

  /* Declarations */

  double et_final;
  /* Algorithm */

  str2et_c(final_epoch, &et_final);
  *nb_time_steps_simu = (int)(ceil( ( et_final - et_initial_epoch ) / dt ) ) + 1; 

  if (*nb_time_steps_simu < 0){
    *nb_time_steps_simu = -(*nb_time_steps_simu);
  }
  
  return 0;

}



/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           previous_index
//  Purpose:        Among all the elements of the array x_before_interpo that are lower than value, computes the maximum index 
//  Assumptions:    x_before_interpo is sorted in the ascending order
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 08/18/2015    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////


int previous_index(int *x_min_index, 
		   double *x_before_interpo, 
		   double value,
		   double size_x_before_interpo){

  /* Declarations */
  int i_index = 0;

  /* Algorithm */

  while ( (x_before_interpo[i_index] <= value) && (i_index<size_x_before_interpo) ){
    i_index++;
  }
  *x_min_index = i_index-1;


       

  return 0;

}


/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           read_gps_tle
//  Purpose:        Read the GPS TLE file and return the name of each GPS satellite output file and the number of GPS in the TLE file
//  Assumptions:    None
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 08/18/2015    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////


int read_gps_tle(char gps_tle_filename[1000],
		 int *nb_gps,
		 char gps_name[N_SATS][1000]){

  /* Declarations */
  FILE *gps_tle_file = NULL;
  int i,j;//, count_nb_of__;
  char *line = NULL;
  char text[256];
  size_t len = 0;
  int first_space;
  ssize_t read;

  /* Algorithm */
  gps_tle_file = fopen(gps_tle_filename, "r");
  if (gps_tle_file == NULL){
    printf("!***!\nThe GPS TLE file:\n%s\ncould not be found. The program will stop.\n!***!\n", gps_tle_filename); MPI_Finalize(); exit(0);
  }
  // Reads the entire file one time to calcualte the number of GPS
  while ( (read = getline(&line, &len, gps_tle_file)) != -1 ){

    if (line[0] == '1'){
      *nb_gps = *nb_gps + 1 ;
    }
  }


  // Now 
  rewind(gps_tle_file);
  for (i = 0; i<*nb_gps; i++){
    first_space = 1;
    getline(&line, &len, gps_tle_file);
    //    printf("<%s> %d\n", line, *nb_gps);
    sscanf(line, "%255[^\n]s", text);  
    getline(&line, &len, gps_tle_file);
    getline(&line, &len, gps_tle_file);

    j = 0;
    while ((unsigned)j < strlen(text)){
      if ( (text[j] == ' ') && (first_space == 0) )  {
	//	strcat(gps_name[i],".txt");
	j = 1000;
      }
      if (j < 999)
	gps_name[i][j] = text[j];
      if ( (text[j] == ' ') && first_space == 1 )  {
	gps_name[i][j] = '_';
	first_space = 0;
      }
      j = j+1;

    }
  }

  fclose(gps_tle_file);

  return 0;
}




void RemoveSpaces(char* source)
{
  char* i = source;
  char* j = source;
  while(*j != 0)
    {
      *i = *j++;
      if(*i != ' ')
	i++;
    }
  *i = 0;
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           calculate_f107_average
//  Purpose:        Write a file of the F10.7 81-day average
//  Assumptions:    None
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 06/08/2016    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////

int calculate_f107_average(char filename_f107_to_calculate_f107_average[256], char filename_f107A[256],  char initial_epoch_wget[8], char final_epoch_wget[8], char correct_end_date_omniweb_temp[8], double et_initial_epoch, double et_final_epoch) //correct_end_date_omniweb_temp: last day of filename_f107_to_calculate_f107_average
{
  /* Declarations */

  int ierr;
  FILE *f107A;
  FILE *file_f107_to_calculate_f107_average;
  int found_eoh;  char *line = NULL;  size_t len = 0; char text[256];
  int nb_element_in_file_f107_to_calculate_f107_average = 0;
  double *f107_file;
  char year_file[8], doy_file[8];
  double  hr_file=0;
  int i;
  int nb_days_in_f107=0;
  int iloop = 0;
  int   iloop_ave;
  int nb_elts_for_average; 
  double f107_average;
  double sum_f107_for_average;
  char time_for_file_f107A[256];
  char h_str[10];
  double et_for_file_f107A;
  char time_for_file_f107A_temp[256];
  int h;
  char initial_epoch_wget_str[256];
  double et_initial_epoch_wget;
  char final_epoch_wget_str[256];
  double et_final_epoch_wget;
  char  year_final_epoch_wget_str[4]; char doy_final_epoch_wget_str[3];
  char correct_end_date_omniweb_temp_str[256];
  double et_correct_end_date_omniweb_temp;
  char  year_correct_end_date_omniweb_temp_str[4]; char doy_correct_end_date_omniweb_temp_str[3];

  /* Algorithm */
  // Open files
  file_f107_to_calculate_f107_average = fopen(filename_f107_to_calculate_f107_average, "r"); // read the F10.7 values for the computation of the 81-day average. This is a daily average file
  f107A = fopen(filename_f107A, "w+"); // file with the results: F10.7A as a function of time (from initial_epoch to final_epoch)
  fprintf(f107A, "YEAR DOY HR F10.7A\n");

  strcpy(initial_epoch_wget_str, "");
  strncat(initial_epoch_wget_str, initial_epoch_wget, 4); strcat(initial_epoch_wget_str, "-"); strncat(initial_epoch_wget_str, initial_epoch_wget+4, 2);  strcat(initial_epoch_wget_str, "-"); strncat(initial_epoch_wget_str, initial_epoch_wget+6, 2); strcat(initial_epoch_wget_str, "T00:00");
  str2et_c(initial_epoch_wget_str, &et_initial_epoch_wget);
  et2utc_c(et_initial_epoch_wget, "D" ,0 ,255 , initial_epoch_wget_str);  

  strcpy(final_epoch_wget_str, "");
  strncat(final_epoch_wget_str, final_epoch_wget, 4); strcat(final_epoch_wget_str, "-"); strncat(final_epoch_wget_str, final_epoch_wget+4, 2);  strcat(final_epoch_wget_str, "-"); strncat(final_epoch_wget_str, final_epoch_wget+6, 2); strcat(final_epoch_wget_str, "T00:00");
  str2et_c(final_epoch_wget_str, &et_final_epoch_wget);
  et2utc_c(et_final_epoch_wget, "D" ,0 ,255 , final_epoch_wget_str);  
  strcpy(year_final_epoch_wget_str, ""); strcpy(doy_final_epoch_wget_str, "");
  strncat(year_final_epoch_wget_str, final_epoch_wget_str, 4);
  strncat(doy_final_epoch_wget_str, final_epoch_wget_str+5,3);

  strcpy(correct_end_date_omniweb_temp_str, "");
  strncat(correct_end_date_omniweb_temp_str, correct_end_date_omniweb_temp, 4); strcat(correct_end_date_omniweb_temp_str, "-"); strncat(correct_end_date_omniweb_temp_str, correct_end_date_omniweb_temp+4, 2);  strcat(correct_end_date_omniweb_temp_str, "-"); strncat(correct_end_date_omniweb_temp_str, correct_end_date_omniweb_temp+6, 2); strcat(correct_end_date_omniweb_temp_str, "T00:00");
  str2et_c(correct_end_date_omniweb_temp_str, &et_correct_end_date_omniweb_temp);
  et2utc_c(et_correct_end_date_omniweb_temp, "D" ,0 ,255 , correct_end_date_omniweb_temp_str);  
  strcpy(year_correct_end_date_omniweb_temp_str, ""); strcpy(doy_correct_end_date_omniweb_temp_str, "");
  strncat(year_correct_end_date_omniweb_temp_str, correct_end_date_omniweb_temp_str, 4);
  strncat(doy_correct_end_date_omniweb_temp_str, correct_end_date_omniweb_temp_str+5,3);

  nb_days_in_f107 = (int)(ceil(( et_final_epoch - et_initial_epoch ) / 3600. / 24.));

  // skip header in file_f107_to_calculate_f107_average
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(file_f107_to_calculate_f107_average)) {
    getline(&line, &len, file_f107_to_calculate_f107_average);
    sscanf(line, "%s", text);
    if (  strcmp( "YEAR", text  ) == 0 )  {
      found_eoh = 1;
    }
  }

  nb_elements_file(&nb_element_in_file_f107_to_calculate_f107_average, filename_f107_to_calculate_f107_average, "YEAR", "</pre><hr><HR>");

  f107_file = malloc( nb_element_in_file_f107_to_calculate_f107_average * sizeof(double) );
  if (  f107_file == NULL ){
    printf("***! Could not allow memory space for f107_file.\nThe program will stop. !***\n");
    ierr =  MPI_Finalize();
    exit(0);
  }
  for (i = 0; i < nb_element_in_file_f107_to_calculate_f107_average; i++){
    getline(&line, &len, file_f107_to_calculate_f107_average);
    sscanf(line, "%s %s %lf %lf", year_file, doy_file, &hr_file, &f107_file[i]);
    /* if (stop_increment == 0){ */
    /* nb_days_in_f107 = nb_days_in_f107 + 1; */
    /*   } */
    /* if ( ( strcmp(year_file, year_final_epoch_wget_str) == 0 ) && ( strcmp(doy_file, doy_final_epoch_wget_str) == 0 ) ){ */
    /*   stop_increment = 1; */
    /* } */

  }

  while (iloop < nb_days_in_f107){ // when iloop reaches nb_days_in_f107 - 1, it corresponds to the final time step of the propagation, so the final step for which we compute F10.7A
    iloop_ave =0;
    sum_f107_for_average = 0;
    nb_elts_for_average = 0;

    while ( ( iloop_ave < 81 ) && ( iloop_ave + iloop < nb_element_in_file_f107_to_calculate_f107_average )  ){
      sum_f107_for_average = sum_f107_for_average + f107_file[iloop_ave + iloop];
      nb_elts_for_average = nb_elts_for_average + 1;
      iloop_ave = iloop_ave + 1;
    }
    f107_average = sum_f107_for_average / nb_elts_for_average;
    for (h = 0; h < 24; h++){
      et_for_file_f107A = et_initial_epoch_wget + iloop * 24 * 3600 + h * 3600.;// we start at et_initial_epoch_wget and increment by one hour because lin_interpolate has been buit for interpolation from one hour time step (and the file for F10.7 and Ap are one hour time step)
      et2utc_c(et_for_file_f107A, "ISOD" ,0 ,255 , time_for_file_f107A_temp);  
      strcpy(time_for_file_f107A, "");
      strncat(time_for_file_f107A,time_for_file_f107A_temp,4);
      strcat(time_for_file_f107A, " ");
      strncat(time_for_file_f107A,time_for_file_f107A_temp+5,3);
      strcat(time_for_file_f107A, " ");
      sprintf(h_str, "%d", h);
      strcat(time_for_file_f107A,h_str );
      fprintf(f107A, "%s %f\n", time_for_file_f107A, f107_average);
    }
    iloop = iloop + 1;
  }

  fclose(file_f107_to_calculate_f107_average);
  fprintf(f107A, "</pre><hr><HR>\n");
  fclose(f107A);

  return 0;
}



// You must free the result if result is non-NULL.
char *str_replace(char *orig, char *rep, char *with) {
  char *result; // the return string
  char *ins;    // the next insert point
  char *tmp;    // varies
  int len_rep;  // length of rep
  int len_with; // length of with
  int len_front; // distance between rep and end of last rep
  int count;    // number of replacements

  if (!orig)
    return NULL;
  if (!rep)
    rep = "";
  len_rep = strlen(rep);
  if (!with)
    with = "";
  len_with = strlen(with);

  ins = orig;
  for (count = 0; (tmp = strstr(ins, rep)); ++count) {
    ins = tmp + len_rep;
  }

  // first time through the loop, all the variable are set correctly
  // from here on,
  //    tmp points to the end of the result string
  //    ins points to the next occurrence of rep in orig
  //    orig points to the remainder of orig after "end of rep"
  tmp = result = malloc(strlen(orig) + (len_with - len_rep) * count + 1);

  if (!result)
    return NULL;

  while (count--) {
    ins = strstr(orig, rep);
    len_front = ins - orig;
    tmp = strncpy(tmp, orig, len_front) + len_front;
    tmp = strcpy(tmp, with) + len_with;
    orig += len_front + len_rep; // move to next "end of rep"
  }
  strcpy(tmp, orig);
  return result;
}


int generate_ensemble_f107_ap(OPTIONS_T *OPTIONS, int iDebugLevel, int iProc){

  // Declarations
  /* int aaa_count; */
  /* int eee;  */
  /* int aaa_sigma; */
  int aaa;
  //  char time[256];
  /* double *ensemble_f107_at_given_time = NULL; */
  /* double *ensemble_f107_at_given_time_sorted = NULL; */
  /* double *ensemble_ap_at_given_time = NULL; */
  /* double *ensemble_ap_at_given_time_sorted = NULL; */


  if ((iProc == 0) & ( iDebugLevel >= 3 ) ){
    printf("---- (load_options) (generate_ensemble_f107_ap) Generating ensembles for F10.7 and Ap.\n");
  }

  /* ensemble_f107_at_given_time = malloc( OPTIONS->nb_ensembles_density * sizeof( double ) );  */
  /* if (ensemble_f107_at_given_time == NULL){ */
  /*   printf("***! (load_options) (generate_ensemble_f107_ap) There is not enough memory for ensemble_f107_at_given_time. The program will stop. !***\n"); MPI_Finalize();exit(0); */
  /* } */
  /* ensemble_f107_at_given_time_sorted = malloc( OPTIONS->nb_ensembles_density *  sizeof( double ) );  */
  /* if (ensemble_f107_at_given_time_sorted == NULL){ */
  /*   printf("***! (load_options) (generate_ensemble_f107_ap) There is not enough memory for ensemble_f107_at_given_time_sorted. The program will stop. !***\n"); MPI_Finalize();exit(0); */
  /* } */

  OPTIONS->f107_ensemble = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double *) ) ; // OPTIONS->nb_time_steps * 2 is a mximum value, it's actaully not the correct nb of elements if there are observations in addition to predictions
  for (aaa = 0; aaa< OPTIONS->nb_time_steps * 2; aaa++){
    OPTIONS->f107_ensemble[aaa] = malloc( OPTIONS->nb_ensembles_density * sizeof(double) );
  }
  if ( OPTIONS->f107_ensemble  == NULL ) {
 printf("***! (load_options) (generate_ensemble_f107_ap) There is not enough memory OPTIONS->f107_ensemble. The program will stop. !***\n"); MPI_Finalize();exit(0);
  }

  /* ensemble_ap_at_given_time = malloc( OPTIONS->nb_ensembles_density * sizeof( double ) );  */
  /* if (ensemble_ap_at_given_time == NULL){ */
  /*   printf("***! (load_options) (generate_ensemble_ap_ap) There is not enough memory for ensemble_ap_at_given_time. The program will stop. !***\n"); MPI_Finalize();exit(0); */
  /* } */
  /* ensemble_ap_at_given_time_sorted = malloc( OPTIONS->nb_ensembles_density *  sizeof( double ) );  */
  /* if (ensemble_ap_at_given_time_sorted == NULL){ */
  /*   printf("***! (load_options) (generate_ensemble_ap_ap) There is not enough memory for ensemble_ap_at_given_time_sorted. The program will stop. !***\n"); MPI_Finalize();exit(0); */
  /* } */

  OPTIONS->Ap_ensemble = malloc( OPTIONS->nb_time_steps * 2 * sizeof(double *) ) ; // OPTIONS->nb_time_steps * 2 is a mximum value, it's actaully not the correct nb of elements if there are observations in addition to predictions
  for (aaa = 0; aaa< OPTIONS->nb_time_steps * 2; aaa++){
    OPTIONS->Ap_ensemble[aaa] = malloc( OPTIONS->nb_ensembles_density * sizeof(double) );
  }
  if ( OPTIONS->Ap_ensemble  == NULL ) {
 printf("***! (load_options) (generate_ensemble_ap_ap) There is not enough memory OPTIONS->ap_ensemble. The program will stop. !***\n"); MPI_Finalize();exit(0);
  }
  
  if (OPTIONS->et_interpo[0] > OPTIONS->swpc_et_first_prediction){
    previous_index( &OPTIONS->aaa_sigma, OPTIONS->et_sigma_f107_ap, OPTIONS->et_interpo[0] - OPTIONS->swpc_et_first_prediction, ( OPTIONS->nb_time_steps + (int)(( 2 * 24 * 3600. ) / OPTIONS->dt))); 

      }
  else{
    OPTIONS->aaa_sigma = 0;
  }

	       //    printf("OPTIONS->aaa_sigma: %d\n", OPTIONS->aaa_sigma);exit(0);
/*   aaa_count = 0; */

/* // go through time steps that corrresponds to the predictions only (not the observations)  */
/*   for (aaa = 0; aaa< OPTIONS->nb_time_steps * 2; aaa++){  // ("* 2.0" because of the Runge Kunta order 4 method) */
/*     if ( OPTIONS->et_interpo[aaa] >= OPTIONS->swpc_et_first_prediction ){ // only overwrite predictions of ensembles (not observations). OPTIONS->et_interpo[aaa] corresponds to the first say of predictions */
/*       //      printf("XXXXXXXXXXXXXX\nXXXXXXXXXXXXXXXXXXX\n"); */
/*       // // Generate nb_ensembles_density normal random values */

/*       for ( eee = 0; eee < OPTIONS->nb_ensembles_density; eee++){ */
/* 	ensemble_f107_at_given_time[eee]   = randn( OPTIONS->f107[aaa], OPTIONS->sigma_f107[aaa_sigma]); */
/* 	ensemble_ap_at_given_time[eee]   = randn(  OPTIONS->Ap[aaa], OPTIONS->sigma_ap[aaa_sigma]); */
/* 	/\* if (eee == 0){ *\/ */
/* 	/\* etprint(OPTIONS->et_interpo[aaa], "time"); *\/ */
/* 	/\* printf("Ap[%d]: %f | sigma_ap[%d]: %f \n",aaa, OPTIONS->Ap[aaa],aaa_sigma, OPTIONS->sigma_ap[aaa_sigma]); *\/ */
/* 	/\* } *\/ */
	
/*       } */

/*       // // Order values in ascending order */
/*       sort_asc_order(ensemble_f107_at_given_time_sorted,  ensemble_f107_at_given_time,OPTIONS->nb_ensembles_density); */
/*       sort_asc_order(ensemble_ap_at_given_time_sorted,  ensemble_ap_at_given_time,OPTIONS->nb_ensembles_density); */

/*       // // Save values */
/*       for ( eee = 0; eee < OPTIONS->nb_ensembles_density; eee++){ */
/* 	OPTIONS->f107_ensemble[aaa_count][eee] = ensemble_f107_at_given_time_sorted[eee]; */
/* 	OPTIONS->Ap_ensemble[aaa_count][eee] = ensemble_ap_at_given_time_sorted[eee]; */
/*       /\* 	if (OPTIONS->Ap_ensemble[aaa_count][eee] < 0){ *\/ */
/*       /\* etprint(OPTIONS->et_interpo[aaa], "\ntime"); *\/ */
/*       /\* printf("ap_ens[%d][%d]: %f - ap_ref[%d]: %f - sigma_ap[%d]: %f\n", aaa_count, eee,  OPTIONS->Ap_ensemble[aaa_count][eee], aaa, OPTIONS->Ap[aaa], aaa_sigma, OPTIONS->sigma_ap[aaa_sigma] ); *\/ */
/*       /\* //      printf("%f %d %d\n",  OPTIONS->f107_ensemble[aaa_count][eee] , eee, aaa_count ); *\/ */
/*       /\* /\\* if () *\\/ *\/ */
/*       /\* 	} *\/ */

/*       } */


/*       aaa_sigma = aaa_sigma + 1; */
/*       aaa_count = aaa_count + 1; */
/*     } */
/*   } */
  if ((iProc == 0) & ( iDebugLevel >= 3 ) ){
    printf("---- (load_options) (generate_ensemble_f107_ap) Done generating ensembles for F10.7 and Ap.\n");
  }


  return 0;
}



/*
 * C program to accept N numbers and arrange them in an ascending order
 */
#include <stdio.h>
 
void sort_asc_order( double *array_out, double *array_in, int n)
{
  int i, j;
    double a;

  for (i = 0; i < n; ++i){
    array_out[i] =  array_in[i];
  }

 
 for (i = 0; i < n; ++i)
    {
      for (j = i + 1; j < n; ++j)
        {
	  if (array_out[i] > array_out[j])
            {
	      a =  array_out[i];
	      array_out[i] = array_out[j];
	      array_out[j] = a;
            }
        }
    }
}



int ptd( double d, char str[256]){
  printf("%s: %.10f\n", str, d);
  return 0;
}

int pti( int i, char str[256]){
  printf("%s: %d\n", str, i);
  return 0;
}


// Assumtpions:
// - for VCM, assumes that  fr is 1 for prograde orbits, -1 for retrograde orbits. This is a convention, but there are others that assume fr to be 1 all the time. make sure that the VCM was generated using the same convention
int ini_collision( OPTIONS_T *OPTIONS, int iProc ){ // !!!!!! the collision input file must be in meters
  //GOOD REFERENCE: http://www.prepacom.net/HEC2/math/cours/Changement%20de%20bases.pdf
  int isc;
  int i;
  // Declarations
       /* double inva[36]; */
       /*     int s; */

  if (OPTIONS->coll_vcm == 0){ // if collison file is not a vcm

  FILE *file_in_collision= NULL;
  char filename_in_collision[256];
  strcpy(filename_in_collision, OPTIONS->filename_input_collision);
  file_in_collision = fopen(filename_in_collision, "r");
  if (file_in_collision == NULL){
    printf("***! The collision input file could not be found (name of the file: %s). The program will stop. !***\n", OPTIONS->filename_input_collision); MPI_Finalize(); exit(0);
  }
  int find_state_eci = 0;
  int find_covariance = 0;
  int find_nb_ensembles = 0;
  char *line = NULL;
  size_t len = 0;
  char text[256];

  // Find state ECI section
  while ( ( find_state_eci == 0 ) && ( !feof( file_in_collision ) ) ){
    getline( &line, &len, file_in_collision );
    sscanf( line, "%s", text);
    if ( strcmp( text, "#STATE_ECI" ) == 0 ){
      find_state_eci = 1;
    }
  }
  if ( feof( file_in_collision ) ){
    print_error(iProc, "No section #STATE_ECI was found in the collision input file");
  }
  
  // Read ECI state for each sc

  for ( isc = 0; isc < OPTIONS->nb_satellites_not_including_gps; isc++ ){
    getline( &line, &len, file_in_collision );
    if (line[0] == '\0'){
      print_error(iProc, "It looks like you did not write a state for each of the spacecraft in the collision input file");
    }      
    RemoveSpaces(line);
    sscanf(line, "(%lf;%lf;%lf)(%lf;%lf;%lf)", &OPTIONS->x_eci[isc], &OPTIONS->y_eci[isc], &OPTIONS->z_eci[isc], &OPTIONS->vx_eci[isc], &OPTIONS->vy_eci[isc], &OPTIONS->vz_eci[isc]);
    OPTIONS->x_eci[isc] = OPTIONS->x_eci[isc] / 1000.; // !!!!!!!!! the collission input file must be in meters
    OPTIONS->y_eci[isc] = OPTIONS->y_eci[isc] / 1000.;
    OPTIONS->z_eci[isc] = OPTIONS->z_eci[isc] / 1000.;
    OPTIONS->vx_eci[isc] = OPTIONS->vx_eci[isc] / 1000.;
    OPTIONS->vy_eci[isc] = OPTIONS->vy_eci[isc] / 1000.;
    OPTIONS->vz_eci[isc] = OPTIONS->vz_eci[isc] / 1000.;
  }

  // Find covariace section
  rewind(file_in_collision);
  while ( ( find_covariance == 0 ) && ( !feof( file_in_collision ) ) ){
    getline( &line, &len, file_in_collision );
    sscanf( line, "%s", text);
    if ( strcmp( text, "#COVARIANCE" ) == 0 ){
      find_covariance = 1;
    }
  }
  if ( feof( file_in_collision ) ){
    print_error(iProc, "No section #COVARIANCE was found in the collision input file");
  }
  
  // Read covariance matrix for each sc


  OPTIONS->covariance_matrix = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof( double **) );
  for ( isc = 0; isc < OPTIONS->nb_satellites_not_including_gps; isc++ ){
    OPTIONS->covariance_matrix[isc] = malloc( 6 * sizeof( double *) );
    for ( i = 0; i < 6; i++ ){
      OPTIONS->covariance_matrix[isc][i] = malloc( 6 * sizeof( double ) );
      if ( OPTIONS->covariance_matrix[isc][i] == NULL ){
	print_error(iProc, "Not enough memory for the covariance matrix");
      }
    }
    if ( OPTIONS->covariance_matrix[isc] == NULL ){
      print_error(iProc, "Not enough memory for the covariance matrix");
    }
  }
  if ( OPTIONS->covariance_matrix == NULL ){
    print_error(iProc, "Not enough memory for the covariance matrix");
  }

  for ( isc = 0; isc < OPTIONS->nb_satellites_not_including_gps; isc++ ){
    for ( i = 0; i < 6; i++ ){
    getline( &line, &len, file_in_collision );
    if (line[0] == '\0'){
      print_error(iProc, "It looks like you did not write a covariance matrix for each of the spacecraft in the collision input file or that you incorrectly wrote it");
    }      
    RemoveSpaces(line);
    if ( i == 0 ){
      sscanf(line, "((%lf;%lf;%lf;%lf;%lf;%lf);", &OPTIONS->covariance_matrix[isc][i][0], &OPTIONS->covariance_matrix[isc][i][1], &OPTIONS->covariance_matrix[isc][i][2], &OPTIONS->covariance_matrix[isc][i][3], &OPTIONS->covariance_matrix[isc][i][4], &OPTIONS->covariance_matrix[isc][i][5]); 
    }
    else if ( i == 5 ){
      sscanf(line, "(%lf;%lf;%lf;%lf;%lf;%lf))", &OPTIONS->covariance_matrix[isc][i][0], &OPTIONS->covariance_matrix[isc][i][1], &OPTIONS->covariance_matrix[isc][i][2], &OPTIONS->covariance_matrix[isc][i][3], &OPTIONS->covariance_matrix[isc][i][4], &OPTIONS->covariance_matrix[isc][i][5]); 
    }
    else{
      sscanf(line, "(%lf;%lf;%lf;%lf;%lf;%lf);", &OPTIONS->covariance_matrix[isc][i][0], &OPTIONS->covariance_matrix[isc][i][1], &OPTIONS->covariance_matrix[isc][i][2], &OPTIONS->covariance_matrix[isc][i][3], &OPTIONS->covariance_matrix[isc][i][4], &OPTIONS->covariance_matrix[isc][i][5]); 
    }
    }
  }


  // Find state #NB_ENSEMBLES_COLLISION section
  rewind(file_in_collision);
  while ( ( find_nb_ensembles == 0 ) && ( !feof( file_in_collision ) ) ){
    getline( &line, &len, file_in_collision );
    sscanf( line, "%s", text);
    if ( strcmp( text, "#NB_ENSEMBLES_COLLISION" ) == 0 ){
      find_nb_ensembles = 1;
    }
  }
  if ( feof( file_in_collision ) ){
    print_error(iProc, "No section #NB_ENSEMBLES_COLLISION was found in the collision input file");
  }
  getline( &line, &len, file_in_collision );
  OPTIONS->nb_ensembles = 0;
  sscanf(line,"%d", &OPTIONS->nb_ensembles);

  // Find state #MIN_DISTANCE_CLOSE_APPROACH section
  int find_min_dist_close_approach = 0;
  while ( ( find_min_dist_close_approach == 0 ) && ( !feof( file_in_collision ) ) ){
    getline( &line, &len, file_in_collision );
    sscanf( line, "%s", text);
    if ( strcmp( text, "#MIN_DISTANCE_CLOSE_APPROACH" ) == 0 ){
      find_min_dist_close_approach = 1;
    }
  }
  if ( feof( file_in_collision ) ){
    print_error(iProc, "No section #MIN_DISTANCE_CLOSE_APPROACH was found in the collision input file");
  }
  getline( &line, &len, file_in_collision );
  OPTIONS->min_dist_close_approach = 0;
  sscanf(line,"%lf", &OPTIONS->min_dist_close_approach);
  OPTIONS->min_dist_close_approach = OPTIONS->min_dist_close_approach / 1000.;  // !!!!!!!!! the collission input file must be in meters



  // Find state #MIN_DISTANCE_COLLISION section
  int find_min_dist_collision = 0;
  while ( ( find_min_dist_collision == 0 ) && ( !feof( file_in_collision ) ) ){
    getline( &line, &len, file_in_collision );
    sscanf( line, "%s", text);
    if ( strcmp( text, "#MIN_DISTANCE_COLLISION" ) == 0 ){
      find_min_dist_collision = 1;
    }
  }
  if ( feof( file_in_collision ) ){
    print_error(iProc, "No section #MIN_DISTANCE_COLLISION was found in the collision input file");
  }
  getline( &line, &len, file_in_collision );
  OPTIONS->min_dist_collision = 0;
  sscanf(line,"%lf", &OPTIONS->min_dist_collision);
  OPTIONS->min_dist_collision = OPTIONS->min_dist_collision / 1000.;  // !!!!!!!!! the collission input file must be in meters
  fclose(file_in_collision);
  } // end of if not a vcm
/*     else{ // if vcm, the information (r/v etc) has been read in the function read_vcm. Need to fill variables and convert equinoctial cov to eci cov */

  else{



/*   OPTIONS->covariance_matrix = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof( double **) ); */
/*   for ( isc = 0; isc < OPTIONS->nb_satellites_not_including_gps; isc++ ){ */
/*     OPTIONS->covariance_matrix[isc] = malloc( 6 * sizeof( double *) ); */
/*     for ( i = 0; i < 6; i++ ){ */
/*       OPTIONS->covariance_matrix[isc][i] = malloc( 6 * sizeof( double ) ); */
/*       if ( OPTIONS->covariance_matrix[isc][i] == NULL ){ */
/* 	print_error(iProc, "Not enough memory for the covariance matrix"); */
/*       } */
/*     } */
/*     if ( OPTIONS->covariance_matrix[isc] == NULL ){ */
/*       print_error(iProc, "Not enough memory for the covariance matrix"); */
/*     } */
/*   } */
/*   if ( OPTIONS->covariance_matrix == NULL ){ */
/*     print_error(iProc, "Not enough memory for the covariance matrix"); */
/*   } */

/*   for ( isc = 0; isc < OPTIONS->nb_satellites_not_including_gps; isc++ ){ */



/*  double rvec[3], vvec[3]; */
/*  rvec[0] = OPTIONS->x_eci[isc]; rvec[1] = OPTIONS->y_eci[isc]; rvec[2] = OPTIONS->z_eci[isc]; */
/*  vvec[0] = OPTIONS->vx_eci[isc]; vvec[1] = OPTIONS->vy_eci[isc]; vvec[2] = OPTIONS->vz_eci[isc]; */
/*  double mu = 398600.4418; // km^3/s^2 (ideally, should read from propagator.c -> load_params but the call to load_params is later) */
/*     double af, ag, lequin, nequin, chi, psi; */
/*  double fr; */
/*  // fr is 1 for prograde orbits, -1 for retrograde orbits. !!!! this is a convention, but there are others that assume fr to be 1 all the time. make sure that the VCM was generated using the same convention */
/*  ORBITAL_ELEMENTS_T oe_temp; */
/*  cart2kep(&oe_temp, rvec, vvec, OPTIONS->et_vcm[isc] , mu); */

/*  if (oe_temp.inclination <= M_PI/2.){ */
/*    fr = 1; */
/*  } */
/*  else{ */
/*    fr = -1; */
/*  } */
/*  // printf("%f %f %f %f %f %f %f\n", oe_temp.sma, oe_temp.eccentricity, oe_temp.inclination*180./M_PI, oe_temp.long_an*180./M_PI, oe_temp.w*180./M_PI, oe_temp.f*180./M_PI, fr); */
/*  double T_equin_to_cart[6][6]; */

/*     cart_to_equin( &af, &ag, &lequin, &nequin, &chi, &psi,  mu,  fr, rvec, vvec); // r and v in km km/s */
/*     double atest = pow(mu / (nequin*nequin),1./3); */
/*     //    printf("equin here %f %f %f %f %f %f\n", atest, af,ag, chi,psi, lequin*180./M_PI); */
/*  compute_T_deriv_equin_to_cart( T_equin_to_cart, af, ag, lequin, nequin, chi, psi, mu, fr); */




/*      double m_eq[6][6]; */

/*      // T_equin_to_cart is to transform equi to cart cov assuming equi cov is: n af ag chi psi l. But VCM gives cov in the order: af ag l n chi psi. Also, the covariance in the VCM is dienisionless -> all terms involvning the mean motion are divied by the mean motion (non-diag) or the squre of the mean motnion (diag). So, here need to unnormalized */
/*      m_eq[0][0] = OPTIONS->covariance_matrix_equinoctial[isc][3][3] * nequin * nequin; */
/*      m_eq[1][0] = OPTIONS->covariance_matrix_equinoctial[isc][3][0] * nequin; m_eq[1][1] = OPTIONS->covariance_matrix_equinoctial[isc][0][0]; */
/*      m_eq[2][0] = OPTIONS->covariance_matrix_equinoctial[isc][3][1] * nequin; m_eq[2][1] = OPTIONS->covariance_matrix_equinoctial[isc][1][0]; m_eq[2][2] = OPTIONS->covariance_matrix_equinoctial[isc][1][1]; */
/*      m_eq[3][0] = OPTIONS->covariance_matrix_equinoctial[isc][4][3] * nequin; m_eq[3][1] = OPTIONS->covariance_matrix_equinoctial[isc][4][0]; m_eq[3][2] = OPTIONS->covariance_matrix_equinoctial[isc][4][1]; m_eq[3][3] = OPTIONS->covariance_matrix_equinoctial[isc][4][4]; */
/*      m_eq[4][0] = OPTIONS->covariance_matrix_equinoctial[isc][5][3] * nequin; m_eq[4][1] = OPTIONS->covariance_matrix_equinoctial[isc][5][0]; m_eq[4][2] = OPTIONS->covariance_matrix_equinoctial[isc][5][1]; m_eq[4][3] = OPTIONS->covariance_matrix_equinoctial[isc][5][4]; m_eq[4][4] = OPTIONS->covariance_matrix_equinoctial[isc][5][5];  */
/*      m_eq[5][0] = OPTIONS->covariance_matrix_equinoctial[isc][3][2] * nequin; m_eq[5][1] = OPTIONS->covariance_matrix_equinoctial[isc][2][0]; m_eq[5][2] = OPTIONS->covariance_matrix_equinoctial[isc][2][1]; m_eq[5][3] = OPTIONS->covariance_matrix_equinoctial[isc][4][2]; m_eq[5][4] = OPTIONS->covariance_matrix_equinoctial[isc][5][2]; m_eq[5][5] = OPTIONS->covariance_matrix_equinoctial[isc][2][2];  */
/* 	     int j; */
/*     for (i = 0; i < 6; i++){ */
/*       for (j = 0; j < 6; j++){ */
/* 	if (j > i){ */
/* 	m_eq[i][j] = m_eq[j][i]; */
/* 	} */
/*       } */
/*     } */


/*     char time_test[260]; */
/*     strcpy(time_test, ""); */
/*     et2utc_c(OPTIONS->et_vcm[isc], "ISOC" ,3 ,255 , time_test); */
/*     char year[10], month[10], hour[10], min[10], sec[10], day[10]; */
/*     if (iProc == 0){ */
/*  printf("%f %f %f %f %f %f %f\n", oe_temp.sma, oe_temp.eccentricity, oe_temp.inclination*180./M_PI, oe_temp.long_an*180./M_PI, oe_temp.w*180./M_PI, oe_temp.f*180./M_PI, fr); */
/*     printf("<%s>\n\n\n", time_test); */
/*     } */
/*     strcpy(year, ""); strcpy(month, ""); strcpy(hour, ""); strcpy(min, ""); strcpy(sec, ""); strcpy(day, ""); */
/*     strncat(year, time_test, 4); */
/*     strncat(month, time_test+5, 2); */
/*     strncat(day, time_test+8, 2); */
/*     strncat(hour, time_test+11, 2); */
/*     strncat(min, time_test+14, 2); */
/*     strncat(sec, time_test+17, 6); */
/*     if (iProc == 0){ */
/*         printf("reci = [%.8f; %.8f; %.8f;];\n", OPTIONS->x_eci[isc], OPTIONS->y_eci[isc], OPTIONS->z_eci[isc]); */
/*     printf("veci = [%.8f; %.8f; %.8f;];\n", OPTIONS->vx_eci[isc], OPTIONS->vy_eci[isc], OPTIONS->vz_eci[isc]); */
/*     printf("aeci = [0.001;0.002;0.003];\n"); */
/*     printf("year = %s;\nmon = %s;\nday = %s;\nhr = %s;\nmin = %s;\nsec = %s;\n", year, month, day, hour, min, sec); */

/*     printf("dut1 =  0.10597;\ndat  = 32;\nxp   =  0.0;\nyp   =  0.0;\nlod  =  0.0;\nterms = 2;\ntimezone= 0;\norder = 106;\n\n"); */
    
/*               m_print6_temp(m_eq, "eqco"); */
/* 	      printf("\n\n"); */
/*     } */


/*      double m_cart[6][6], m_cart_temp[6][6]; */
/*      m_x_m6bis( m_cart_temp, T_equin_to_cart, m_eq ); */
/*      double T_equin_to_cart_trans[6][6]; */
/*      m_trans6(   T_equin_to_cart_trans, T_equin_to_cart); */
/*      m_x_m6bis( m_cart, m_cart_temp, T_equin_to_cart_trans); */
/*     if (iProc == 0){ */
/*       //      m_print6(T_equin_to_cart,"T_equin_to_cart"); */
/*         printf("\nCartesian\n"); */
/*     } */

/*     for (i = 0; i < 6; i++){ */
/*       for (j = 0; j < 6; j++){ */
/* 	m_cart[i][j] = m_cart[i][j] * 1.e6; // km2 to m2 */
/* 	OPTIONS->covariance_matrix[isc][i][j] = m_cart[i][j];  */
/* 	if (iProc == 0){ */
/* 	  	  	  printf("%e ",  OPTIONS->covariance_matrix[isc][i][j]); */
/* 	} */
/*       } */
/*       if (iProc == 0){ */
/* 	printf("\n"); */
/*       } */
/*     } */
/*     /\* // !!!!!! REMOVE BLOCK BELOW *\/ */
/*     /\* if (isc == 0){ *\/ */
/*     /\*   OPTIONS->covariance_matrix[isc][0][0] = 2.616647e+03; OPTIONS->covariance_matrix[isc][0][1] = 7.637742e+03; OPTIONS->covariance_matrix[isc][0][2] = -4.370318e+03; OPTIONS->covariance_matrix[isc][0][3] = 4.411248e+00; OPTIONS->covariance_matrix[isc][0][4] = -1.718072e+00; OPTIONS->covariance_matrix[isc][0][5] = -2.078885e+00 ; *\/ */
/*     /\*   OPTIONS->covariance_matrix[isc][1][0] = 7.637742e+03; OPTIONS->covariance_matrix[isc][1][1] = 2.308465e+04; OPTIONS->covariance_matrix[isc][1][2] = -1.303136e+04; OPTIONS->covariance_matrix[isc][1][3] = 1.313077e+01; OPTIONS->covariance_matrix[isc][1][4] = -5.425810e+00; OPTIONS->covariance_matrix[isc][1][5] = -6.115084e+00 ; *\/ */
/*     /\*   OPTIONS->covariance_matrix[isc][2][0] = -4.370318e+03; OPTIONS->covariance_matrix[isc][2][1] = -1.303136e+04; OPTIONS->covariance_matrix[isc][2][2] = 7.446650e+03; OPTIONS->covariance_matrix[isc][2][3] = -7.490541e+00; OPTIONS->covariance_matrix[isc][2][4] = 2.972651e+00; OPTIONS->covariance_matrix[isc][2][5] = 3.558777e+00 ; *\/ */
/*     /\*   OPTIONS->covariance_matrix[isc][3][0] = 4.411248e+00; OPTIONS->covariance_matrix[isc][3][1] = 1.313077e+01; OPTIONS->covariance_matrix[isc][3][2] = -7.490541e+00; OPTIONS->covariance_matrix[isc][3][3] = 7.561292e-03; OPTIONS->covariance_matrix[isc][3][4] = -2.980035e-03; OPTIONS->covariance_matrix[isc][3][5] = -3.582930e-03 ; *\/ */
/*     /\*   OPTIONS->covariance_matrix[isc][4][0] = -1.718072e+00; OPTIONS->covariance_matrix[isc][4][1] = -5.425810e+00; OPTIONS->covariance_matrix[isc][4][2] = 2.972651e+00; OPTIONS->covariance_matrix[isc][4][3] = -2.980035e-03; OPTIONS->covariance_matrix[isc][4][4] = 1.426881e-03; OPTIONS->covariance_matrix[isc][4][5] = 1.282848e-03 ; *\/ */
/*     /\*   OPTIONS->covariance_matrix[isc][5][0] = -2.078885e+00; OPTIONS->covariance_matrix[isc][5][1] = -6.115084e+00; OPTIONS->covariance_matrix[isc][5][2] = 3.558777e+00; OPTIONS->covariance_matrix[isc][5][3] = -3.582930e-03; OPTIONS->covariance_matrix[isc][5][4] = 1.282848e-03; OPTIONS->covariance_matrix[isc][5][5] = 1.821817e-03 ; *\/ */
/*     /\* } *\/ */
/*     /\* else{ *\/ */
/*     /\*   OPTIONS->covariance_matrix[isc][0][0] = 1.278413e+04; OPTIONS->covariance_matrix[isc][0][1] = -1.146899e+04; OPTIONS->covariance_matrix[isc][0][2] = -1.443598e+04; OPTIONS->covariance_matrix[isc][0][3] = -6.705461e+00; OPTIONS->covariance_matrix[isc][0][4] = 3.354092e+00; OPTIONS->covariance_matrix[isc][0][5] = -9.703897e+00 ; *\/ */
/*     /\*   OPTIONS->covariance_matrix[isc][1][0] = -1.146899e+04; OPTIONS->covariance_matrix[isc][1][1] = 1.039556e+04; OPTIONS->covariance_matrix[isc][1][2] = 1.296609e+04; OPTIONS->covariance_matrix[isc][1][3] = 6.054387e+00; OPTIONS->covariance_matrix[isc][1][4] = -3.024363e+00; OPTIONS->covariance_matrix[isc][1][5] = 8.695377e+00 ; *\/ */
/*     /\*   OPTIONS->covariance_matrix[isc][2][0] = -1.443598e+04; OPTIONS->covariance_matrix[isc][2][1] = 1.296609e+04; OPTIONS->covariance_matrix[isc][2][2] = 1.656856e+04; OPTIONS->covariance_matrix[isc][2][3] = 7.519223e+00; OPTIONS->covariance_matrix[isc][2][4] = -3.733794e+00; OPTIONS->covariance_matrix[isc][2][5] = 1.117779e+01 ; *\/ */
/*     /\*   OPTIONS->covariance_matrix[isc][3][0] = -6.705461e+00; OPTIONS->covariance_matrix[isc][3][1] = 6.054387e+00; OPTIONS->covariance_matrix[isc][3][2] = 7.519223e+00; OPTIONS->covariance_matrix[isc][3][3] = 3.564894e-03; OPTIONS->covariance_matrix[isc][3][4] = -1.766507e-03; OPTIONS->covariance_matrix[isc][3][5] = 5.033724e-03 ; *\/ */
/*     /\*   OPTIONS->covariance_matrix[isc][4][0] = 3.354092e+00; OPTIONS->covariance_matrix[isc][4][1] = -3.024363e+00; OPTIONS->covariance_matrix[isc][4][2] = -3.733794e+00; OPTIONS->covariance_matrix[isc][4][3] = -1.766507e-03; OPTIONS->covariance_matrix[isc][4][4] = 9.472415e-04; OPTIONS->covariance_matrix[isc][4][5] = -2.476999e-03 ; *\/ */
/*     /\*   OPTIONS->covariance_matrix[isc][5][0] = -9.703897e+00; OPTIONS->covariance_matrix[isc][5][1] = 8.695377e+00; OPTIONS->covariance_matrix[isc][5][2] = 1.117779e+01; OPTIONS->covariance_matrix[isc][5][3] = 5.033724e-03; OPTIONS->covariance_matrix[isc][5][4] = -2.476999e-03; OPTIONS->covariance_matrix[isc][5][5] = 7.584359e-03 ; *\/ */

/*     /\* } *\/ */

/*     /\* // !!!!!! end of REMOVE BLOCK BELOW *\/ */

/*       if (iProc == 0){ */
/* 		printf("\n"); */
/*       } */

/*       double T_inrtl_2_ntw_6by6[6][6], T_inrtl_2_ntw_6by6_trans[6][6]; */
/*       compute_T_inrtl_2_ntw_6by6( T_inrtl_2_ntw_6by6, rvec, vvec); */
/*       double m_ntw[6][6], m_ntw_temp[6][6]; */
/*            m_x_m6bis( m_ntw_temp, T_inrtl_2_ntw_6by6, m_cart); */
/*      m_trans6(   T_inrtl_2_ntw_6by6_trans, T_inrtl_2_ntw_6by6); */
/*      m_x_m6bis( m_ntw, m_ntw_temp, T_inrtl_2_ntw_6by6_trans); */
/*       if (iProc == 0){ */
/* 	//	     m_print6(T_inrtl_2_ntw_6by6, "T_inrtl_2_ntw_6by6"); */
/* 	     //	     m_print6(m_ntw_temp, "m_ntw_temp"); */
/* /\* 	     double T_inrtl_2_ntw_test[3][3]; *\/ */
/* /\* 	     compute_T_inrtl_2_ntw(T_inrtl_2_ntw_test, rvec, vvec); *\/ */
/* /\* 	     m_print(T_inrtl_2_ntw_test, "T_inrtl_2_ntw"); *\/ */
/* 	m_print6(m_ntw, "NTW"); */
/* 	     printf("radial %.4f, in-track %.4f, cross  %.4f (km)\nv_radial  %.4f, v_in-track  %.4f, v_cross  %.4f (km/s)\n", sqrt(m_ntw[0][0])/1000., sqrt(m_ntw[1][1])/1000., sqrt(m_ntw[2][2])/1000., sqrt(m_ntw[3][3])/1000., sqrt(m_ntw[4][4])/1000., sqrt(m_ntw[5][5])/1000.);  */
/*       } */
      
   //    exitf();

      //  exitf();
    //  } // end go over all sc


  } // end of if colision file has VCM format
  //  exitf();
  ////////////// DIAGNOALZATION OF COVARIANCES MATRICES (source: https://www.gnu.org/software/gsl/manual/html_node/Eigenvalue-and-Eigenvector-Examples.html)
  int j;
  double data[36];
  // Allocate memory for eigenvalues, eigenvectors, and rotation matrix
  OPTIONS->eigenvalue_covariance_matrix = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof(double *) );
  for (isc = 0; isc < OPTIONS->nb_satellites_not_including_gps; isc++){
    OPTIONS->eigenvalue_covariance_matrix[isc] = malloc( 6 * sizeof( double ) );
    if (  OPTIONS->eigenvalue_covariance_matrix[isc] == NULL ){
      print_error(iProc, "Not enough memory for eigenvalues of the covariance matrix for collision assessment");
    }
  }
  if ( OPTIONS->eigenvalue_covariance_matrix == NULL ){
    print_error(iProc, "Not enough memory for eigenvalues of the covariance matrix for collision assessment");
  }

  double ***eigenvector_covariance_matrix;
  eigenvector_covariance_matrix = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof(double **) );
  for (isc = 0; isc < OPTIONS->nb_satellites_not_including_gps; isc++){
    eigenvector_covariance_matrix[isc] = malloc( 6 * sizeof( double *) );
    for ( i = 0; i < 6; i++ ){
      eigenvector_covariance_matrix[isc][i] = malloc( 6 * sizeof( double) );
    if (  eigenvector_covariance_matrix[isc][i] == NULL ){
      print_error(iProc, "Not enough memory for eigenvectors of the covariance matrix for collision assessment");
    }

    }
    if (  eigenvector_covariance_matrix[isc] == NULL ){
      print_error(iProc, "Not enough memory for eigenvectors of the covariance matrix for collision assessment");
    }
  }
  if ( eigenvector_covariance_matrix == NULL ){
    print_error(iProc, "Not enough memory for eigenvectors of the covariance matrix for collision assessment");
  }


  OPTIONS->rotation_matrix_for_diagonalization = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof( double **) );
  for ( isc = 0; isc < OPTIONS->nb_satellites_not_including_gps; isc++ ){
    OPTIONS->rotation_matrix_for_diagonalization[isc] = malloc( 6 * sizeof( double *) );
    for ( i = 0; i < 6; i++ ){
      OPTIONS->rotation_matrix_for_diagonalization[isc][i] = malloc( 6 * sizeof( double ) );
      if ( OPTIONS->rotation_matrix_for_diagonalization[isc][i] == NULL ){
	print_error(iProc, "Not enough memory for the covariance matrix");
      }
    }
    if ( OPTIONS->rotation_matrix_for_diagonalization[isc] == NULL ){
      print_error(iProc, "Not enough memory for the covariance matrix");
    }
  }
  if ( OPTIONS->rotation_matrix_for_diagonalization  == NULL ){
    print_error(iProc, "Not enough memory for the covariance matrix");
  }



  OPTIONS->inverse_rotation_matrix_for_diagonalization = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof( double **) );
  for ( isc = 0; isc < OPTIONS->nb_satellites_not_including_gps; isc++ ){
    OPTIONS->inverse_rotation_matrix_for_diagonalization[isc] = malloc( 6 * sizeof( double *) );
    for ( i = 0; i < 6; i++ ){
      OPTIONS->inverse_rotation_matrix_for_diagonalization[isc][i] = malloc( 6 * sizeof( double ) );
      if ( OPTIONS->inverse_rotation_matrix_for_diagonalization[isc][i] == NULL ){
	print_error(iProc, "Not enough memory for the covariance matrix");
      }
    }
    if ( OPTIONS->inverse_rotation_matrix_for_diagonalization[isc] == NULL ){
      print_error(iProc, "Not enough memory for the covariance matrix");
    }
  }
  if ( OPTIONS->inverse_rotation_matrix_for_diagonalization   == NULL ){
    print_error(iProc, "Not enough memory for the covariance matrix");
  }



  /* OPTIONS->inverse_rotation_matrix_for_diagonalization_copy = malloc( OPTIONS->nb_satellites_not_including_gps * sizeof( double **) ); */
  /* for ( isc = 0; isc < OPTIONS->nb_satellites_not_including_gps; isc++ ){ */
  /*   OPTIONS->inverse_rotation_matrix_for_diagonalization_copy[isc] = malloc( 6 * sizeof( double *) ); */
  /*   for ( i = 0; i < 6; i++ ){ */
  /*     OPTIONS->inverse_rotation_matrix_for_diagonalization_copy[isc][i] = malloc( 6 * sizeof( double ) ); */
  /*     if ( OPTIONS->inverse_rotation_matrix_for_diagonalization_copy[isc][i] == NULL ){ */
  /* 	print_error(iProc, "Not enough memory for the covariance matrix"); */
  /*     } */
  /*   } */
  /*   if ( OPTIONS->inverse_rotation_matrix_for_diagonalization_copy[isc] == NULL ){ */
  /*     print_error(iProc, "Not enough memory for the covariance matrix"); */
  /*   } */
  /* } */
  /* if ( OPTIONS->inverse_rotation_matrix_for_diagonalization_copy   == NULL ){ */
  /*   print_error(iProc, "Not enough memory for the covariance matrix"); */
  /* } */


  // GOOD REFERENCE: http://www.prepacom.net/HEC2/math/cours/Changement%20de%20bases.pdf
  /* if (iProc == 0){ */
  /* m_print6(OPTIONS->covariance_matrix[0], "M"); */
  /* } */
  for ( isc = 0; isc < OPTIONS->nb_satellites_not_including_gps; isc++){
    // Convert covariance matrix in correct format for GSL library
    if ( strcmp(OPTIONS->type_orbit_initialisation, "collision" ) == 0 ){
      for (i = 0; i < 6; i++){
	for (j = 0; j < 6; j++){
	  data[i*6 + j] = OPTIONS->covariance_matrix[isc][i][j]; // !!!!!!!! should be data[i*6 + j] = OPTIONS->covariance_matrix[isc][i][j];      
	}
      }
    }
    else if ( strcmp(OPTIONS->type_orbit_initialisation, "collision_vcm" ) == 0 ){
      for (i = 0; i < 6; i++){
	for (j = 0; j < 6; j++){
	  data[i*6 + j] = OPTIONS->covariance_matrix_equinoctial_only_rv[isc][i][j]; // !!!!!!!! should be data[i*6 + j] = OPTIONS->covariance_matrix[isc][i][j];      
	}
      }
    }


    // Allocate variables for GSL libaries
  gsl_matrix_view m = gsl_matrix_view_array (data, 6, 6);
  gsl_vector_complex *eval = gsl_vector_complex_alloc (6);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (6, 6);
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (6);  

// Diagonalization
  gsl_eigen_nonsymmv (&m.matrix, eval, evec,  w); 
  // Write eigenvalues in OPTIONS->eigenvalue_covariance_matrix[isc] and eigenvectors in eigenvector_covariance_matrix[isc] for each sc
  //  print_test();
    for (i = 0; i < 6; i++){
      // Eigenvalues
        gsl_complex eval_i = gsl_vector_complex_get (eval, i);
	OPTIONS->eigenvalue_covariance_matrix[isc][i] = GSL_REAL(eval_i);
	/* if (iProc==0){ */
	/* /\* if (isc == 1){ *\/ */
	/*   //	  printf("lambda[%d] = %g\n",i, OPTIONS->eigenvalue_covariance_matrix[isc][i]); */
	/*   //	  printf("lambda[%d] = %g\n",i, eval[i]); */
	/* /\* } *\/ */
	/*   } */
	if  ( OPTIONS->eigenvalue_covariance_matrix[isc][i] < 0 ){
	  print_error(iProc, "The eigenvalue of the covariance matrix is negative");
	}
	/* printf("balbalbalbalbal %e %d %d\n", OPTIONS->eigenvalue_covariance_matrix[isc][i], isc, i); */
	/* exitf(); */
	// Eigenvectors
	gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);

	for ( j = 0; j < 6; j++ ){
	  gsl_complex z =  gsl_vector_complex_get(&evec_i.vector, j);
	  eigenvector_covariance_matrix[isc][i][j] = GSL_REAL(z);
	  OPTIONS->rotation_matrix_for_diagonalization[isc][j][i] = eigenvector_covariance_matrix[isc][i][j];
	  //	  OPTIONS->inverse_rotation_matrix_for_diagonalization_copy[isc][i][j] = eigenvector_covariance_matrix[isc][i][j];
	  /* if (iProc==0){ */
	  /*   if (isc == 1){ */
	      
	  /*     	  	    ptd(eigenvector_covariance_matrix[isc][i][j] , "e"); */
	  /*   } */
	  /* } */
	}
	//		printf("\n");
    }
    //    m_print6(OPTIONS->rotation_matrix_for_diagonalization[isc], "Rotation matrix" );
  

    /* for (i = 0; i < 6; i++){ */
    /*   for (j = 0; j < 6; j++){ */
    /* 	data[i*6 + j] = OPTIONS->rotation_matrix_for_diagonalization[isc][i][j]; */
    /*   } */
    /* } */

    // Allocate variables for GSL libaries
      /* gsl_matrix_view m_to_inv = gsl_matrix_view_array (data, 6, 6); */
      /* gsl_matrix_view inv = gsl_matrix_view_array(inva,6,6); */
      /* gsl_permutation * p = gsl_permutation_alloc (6); */

      /* gsl_linalg_LU_decomp (&m_to_inv.matrix, p, &s); */
      /* gsl_linalg_LU_invert (&m_to_inv.matrix, p, &inv.matrix); */

      /* for (i = 0; i < 6; ++i){ */
      /* 	for (j = 0; j < 6; ++j){ */
      /* 	  OPTIONS->inverse_rotation_matrix_for_diagonalization[isc][i][j] = gsl_matrix_get(&inv.matrix,i,j); */
      /* 	} */
      /* } */
      /* //      m_print6(OPTIONS->inverse_rotation_matrix_for_diagonalization[isc], "OPTIONS->inverse_rotation_matrix_for_diagonalization[isc]"); */
      /* gsl_permutation_free (p); */
    
    /* // exitall(); */
    /* /\*   if ( (iProc == 0) ){ *\/ */
    /* /\*        m_print6(OPTIONS->inverse_rotation_matrix_for_diagonalization[isc], "inverse rot mat"); *\/ */
    /* /\* } *\/ */

  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
  gsl_eigen_nonsymmv_free(w);



/*   print_test(); */
/*   m_print6_temp(OPTIONS->covariance_matrix_equinoctial_only_rv[isc], "eqco6"); */
/*   m_print6_temp(OPTIONS->rotation_matrix_for_diagonalization[isc], "rot6"); */
/*   printf("\nXXX\nEIGEN VALUES SC %d\n", isc); */
/*   for (i = 0; i < 6; i++){ */
/*         printf("%15.10e\n", OPTIONS->eigenvalue_covariance_matrix[isc][i]); */
/*   } */

/*   exitf(); */
  }


  //         exitall();
  /* // Check diagonalization */
  /* double **mat_temp = NULL, **mat_temp2 = NULL; */
  /* mat_temp = malloc( 6 * sizeof( double * ) ); */
  /* if (mat_temp == NULL){ */
  /*   print_error(iProc, "Not enough memory for check matrix mat_temp"); */
  /* } */
  /* for (i = 0; i < 6; i++){ */
  /*   mat_temp[i] = malloc( 6 * sizeof( double ) ); */
  /* if (mat_temp[i] == NULL){ */
  /*   print_error(iProc, "Not enough memory for check matrix mat_temp"); */
  /* } */

  /* } */
  /* mat_temp2 = malloc( 6 * sizeof( double * ) ); */
  /* if (mat_temp2 == NULL){ */
  /*   print_error(iProc, "Not enough memory for check matrix mat_temp2"); */
  /* } */

  /* for (i = 0; i < 6; i++){ */
  /*   mat_temp2[i] = malloc( 6 * sizeof( double ) ); */
  /* if (mat_temp2[i] == NULL){ */
  /*   print_error(iProc, "Not enough memory for check matrix mat_temp2"); */
  /* } */

  /* } */

  /* if ( (iProc == 0) ){ */
  /*   isc = 0; */

  /* /\* m_x_m6( mat_temp, OPTIONS->covariance_matrix[isc], OPTIONS->rotation_matrix_for_diagonalization[isc] ); *\/ */
  /* /\* m_x_m6( mat_temp2, OPTIONS->inverse_rotation_matrix_for_diagonalization[isc], mat_temp ); *\/ */
  /* m_x_m6(mat_temp2, OPTIONS->inverse_rotation_matrix_for_diagonalization[isc], OPTIONS->rotation_matrix_for_diagonalization[isc]); */
  /* m_print6(mat_temp2, "inverse sc 0"); */
  /* /\* m_x_m6(mat_temp2, OPTIONS->inverse_rotation_matrix_for_diagonalization_copy[isc], OPTIONS->rotation_matrix_for_diagonalization[isc]); *\/ */
  /* /\* m_print6(mat_temp2, "XXXXXXXXXXX"); *\/ */

  /*   isc = 1; */

  /* /\* m_x_m6( mat_temp, OPTIONS->covariance_matrix[isc], OPTIONS->rotation_matrix_for_diagonalization[isc] ); *\/ */
  /* /\* m_x_m6( mat_temp2, OPTIONS->inverse_rotation_matrix_for_diagonalization[isc], mat_temp ); *\/ */
  /* m_x_m6(mat_temp2, OPTIONS->inverse_rotation_matrix_for_diagonalization[isc], OPTIONS->rotation_matrix_for_diagonalization[isc]); */
  /* m_print6(mat_temp2, "inverse sc 1"); */

  /* /\* for (i = 0; i < 6; i++){ *\/ */
  /* /\*   printf("%15.10e\n", OPTIONS->eigenvalue_covariance_matrix[isc][i]); *\/ */
  /* /\* } *\/ */
  /* } */
  /* // END of Check diagonalization */


  //  exitall();




  return 0;
}



int print_error(int iProc,  char *error_message){

  if ( iProc == 0 ){
  printf("***! %s. The program will stop. !***\n", error_message); MPI_Finalize(); exit(0);
  }
  return 0;
}

int print_error_any_iproc( int iProc,  char *error_message){

  printf("***! %s. The program will stop. (iProc: %d) !***\n", error_message, iProc); MPI_Finalize(); exit(0);

  return 0;
}



int exitall(){
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  exit(0);
  return 0;
}


int print_oe(ORBITAL_ELEMENTS_T OE, PARAMS_T *PARAMS){

  double alt_apogee = OE.sma * ( 1 + OE.eccentricity ) - PARAMS->EARTH.radius;
  ptd(alt_apogee, "alt_apogee");
  ptd(OE.inclination * RAD2DEG, "inclination");
  ptd(OE.w * RAD2DEG, "arg_perigee");
  ptd(OE.long_an * RAD2DEG, "RAAN");
  ptd(OE.f * RAD2DEG, "true_ano");
  ptd(OE.eccentricity, "eccentricity");
  ptd(OE.sma, "sma");
  ptd( fmod((OE.w + OE.f) * RAD2DEG, 360) , "arg_perigee + true_anomaly");
  return 0;
}


int compute_time_interpo(OPTIONS_T *OPTIONS){

  double fake_d;
  double et_meas;
  int found_eoh;
  char *line = NULL;
  char text[300];
  int ierr;
  size_t len = 0;
  FILE *fp_meas = NULL;
  FILE *fp_kalman_init = NULL;
  char filename_meas[300];
  int nb_time_steps_temp;

  // Calculate the maximum number of time steps to allo memory for OPTIONS->et_interpo
  nb_time_steps(&nb_time_steps_temp, OPTIONS->et_oldest_tle_epoch, OPTIONS->final_epoch, OPTIONS->dt);// this numbe of steps is added to the number of measurements to get the max number of times steps
  // // get name of measurement file 
  fp_kalman_init = fopen(OPTIONS->filename_kalman_init, "r");  
  getline(&line, &len, fp_kalman_init);
  strtok(line, "\n");  strtok(line, "\r"); 
  strcpy(filename_meas, "");
  sscanf(line, "%s", filename_meas);
  fclose(fp_kalman_init);
  fp_meas = fopen(filename_meas, "r");
  // // skip header of measurment file
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp_meas)) {
    getline(&line, &len, fp_meas);
    sscanf(line, "%s", text);
    if (  strcmp( text, "#START" ) == 0  )  {
      found_eoh = 1;
    }
  }
  if (feof(fp_meas)){
    printf("***! No measurements were found. The program will stop. !***\n");
    ierr =  MPI_Finalize();
    exit(0);
  }
  int nb_time_steps_meas_total = 0;
  while (feof(fp_meas) == 0){
    getline(&line, &len, fp_meas);
    nb_time_steps_meas_total = nb_time_steps_meas_total + 1;
  }
  rewind(fp_meas);
  int nb_time_steps_max = nb_time_steps_meas_total + nb_time_steps_temp;
  OPTIONS->et_interpo = malloc( nb_time_steps_max * sizeof(double) * 2);
  double *et_interpo_temp;
  et_interpo_temp = malloc( nb_time_steps_max * sizeof(double) );
  
  // // skip header of measurment file (again since we used rewind)
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp_meas)) {
    getline(&line, &len, fp_meas);
    sscanf(line, "%s", text);
    if (  strcmp( text, "#START" ) == 0  )  {
      found_eoh = 1;
    }
  }
  if (feof(fp_meas)){
    printf("***! No measurements were found. The program will stop. !***\n");
    ierr =  MPI_Finalize();
    exit(0);
  }

  double et_initial_epoch, et_final_epoch;
  str2et_c(OPTIONS->initial_epoch, &et_initial_epoch);
  str2et_c(OPTIONS->final_epoch, &et_final_epoch);
  
  int i=0;
  double et_next_time_step;


  getline(&line, &len, fp_meas);
  if (feof(fp_meas) == 0){
    sscanf(line, "%s %lf %lf %lf %lf %lf %lf", text, &fake_d,  &fake_d,  &fake_d,  &fake_d,  &fake_d,  &fake_d);
    str2et_c(text, &et_meas);
  }
  // skip all observations before initial epoch
  while (et_meas < et_initial_epoch){
    getline(&line, &len, fp_meas);
    if (feof(fp_meas) != 0){
      printf("***! No measurements more recent than the initial epoch were found. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    else{
      sscanf(line, "%s %lf %lf %lf %lf %lf %lf", text, &fake_d,  &fake_d,  &fake_d,  &fake_d,  &fake_d,  &fake_d);
      str2et_c(text, &et_meas);
    }
  }
  if (et_meas > et_final_epoch){
      printf("***! No measurements older than the final epoch were found. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);

  }
  // the first time step is equal to the first measurement (since in the KF, we initialize the r/v as the first measuremed r/v)
  int already_propagated_until_time_step_right_after_meas;
  if ( et_meas > et_initial_epoch ){
    et_interpo_temp[0] = et_initial_epoch;
    while ( et_interpo_temp[i] < et_meas ){
      i = i+1;
      et_interpo_temp[i] = et_initial_epoch + OPTIONS->dt * i; // the /2 for Runge kutta is done after (in OPTIONS->et_interpo)
    }
  }
  et_interpo_temp[i] = et_meas;
  //  et_next_time_step = et_initial_epoch + (int)( ( et_meas - et_initial_epoch ) / OPTIONS->dt ) * OPTIONS->dt + OPTIONS->dt;
  double dt_temp = OPTIONS->dt;
  while ( feof(fp_meas) == 0 ){
    et_next_time_step = et_initial_epoch + (int)( ( et_meas - et_initial_epoch ) / OPTIONS->dt ) * OPTIONS->dt + OPTIONS->dt;
    //       etprint(et_next_time_step , "next");
    getline(&line, &len, fp_meas);
    sscanf(line, "%s %lf %lf %lf %lf %lf %lf", text, &fake_d,  &fake_d,  &fake_d,  &fake_d,  &fake_d,  &fake_d);
    //printf("<%s>\n", text);
    str2et_c(text, &et_meas);
    already_propagated_until_time_step_right_after_meas = 0;

    while ( et_interpo_temp[i] != et_meas ) {
      i = i+1;
      //      printf("i = %d | %d\n", i, nb_time_steps_max);
      if ( ( et_meas - et_interpo_temp[i-1]  ) < dt_temp ) {
	dt_temp = ( et_meas - et_interpo_temp[i-1] ) ; 
      }

      if ( ( ( et_next_time_step - et_interpo_temp[i-1]  ) < dt_temp ) && ( already_propagated_until_time_step_right_after_meas == 0 ) ) {
	dt_temp = ( et_next_time_step - et_interpo_temp[i-1] ) ;
	et_interpo_temp[i] = et_next_time_step;
	already_propagated_until_time_step_right_after_meas = 1;
      }
      else{
	et_interpo_temp[i] = et_interpo_temp[i-1] + dt_temp; 
      }

	dt_temp = OPTIONS->dt;// if got in the if condition above "if ( ( et_next_time_step - et_interpo_temp[i-1]  ) < dt_temp )" then we need to reset dt_temp to the initial value (OPTIONS->dt). 
       
      if ( et_interpo_temp[i] > et_final_epoch ){
	break;
      }
      //            etprint(et_interpo_temp[i], "interpo");
    }
      if ( et_interpo_temp[i] > et_final_epoch ){
	break;
      }

  }
  OPTIONS->nb_time_steps = i+1;

  int j;
  for (j = 0; j < OPTIONS->nb_time_steps; j ++){
    
    OPTIONS->et_interpo[2*j] = et_interpo_temp[j];
    if (j != (OPTIONS->nb_time_steps - 1)){
    OPTIONS->et_interpo[2*j+1] = (et_interpo_temp[j+1]+et_interpo_temp[j])/2.; 
    }
  }
  OPTIONS->et_interpo[2*OPTIONS->nb_time_steps-1] = et_interpo_temp[OPTIONS->nb_time_steps-1];
/*   for (j = 0; j < OPTIONS->nb_time_steps; j ++){ */
/*     etprint(et_interpo_temp[j], ""); */
/*   } */
/*   for (j = 0; j < OPTIONS->nb_time_steps*2; j ++){ */
/*     etprint(OPTIONS->et_interpo[j], ""); */
/*   } */
  
  /* fclose(fp_meas); */
  //       MPI_Finalize(); exit(0);
  return 0;
}



int exitf(){
  MPI_Finalize();
  exit(0);
  return 0;
  
}

// read_cdm: reads collision cdm file.
/* assumptions: */
// the order of object (names of VCM files) in SpOCK main input file is the same as the order of objects in CDM file
int read_cdm(char filename[1000], OPTIONS_T *OPTIONS){
    FILE *fp;
    char *line = NULL;
    size_t len = 0;

    fp = fopen(filename, "r");
    getline(&line, &len, fp);
    char *next, *start, *end;
    int found = 0;
    char text[500];
    char epoch[100];    char epoch1[100];     char epoch2[100];    
    strcpy(epoch, "");    strcpy(epoch1, "");     strcpy(epoch2, "");    
    double et;

    char dummy[100];

    int isc = 0;
    // BC
    double bc_cdm;
    while ( found == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s %s %lf %s", text, dummy, &bc_cdm, dummy);
    if (strcmp(text, "CD_AREA_OVER_MASS") == 0){
      OPTIONS->bc_cdm[isc] = bc_cdm;
      found = 1;
    }
    }

    // SRP
    rewind(fp);
    found = 0;
    double srp_cdm;
    while ( found == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s %s %lf %s", text, dummy, &srp_cdm, dummy);
    if (strcmp(text, "CR_AREA_OVER_MASS") == 0){
      OPTIONS->srp_cdm[isc] = srp_cdm;
      found = 1;
    }
    }

    // VARIANCE BC
    rewind(fp);
    found = 0;
    double bc_cdm_std;
    while ( found == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s %s %lf %s", text, dummy, &bc_cdm_std, dummy);
    if (strcmp(text, "CDRG_DRG") == 0){
      OPTIONS->bc_cdm_std[isc] = sqrt(bc_cdm_std);
      found = 1;
    }
    }

    // VARIANCE SRP
    rewind(fp);
    found = 0;
    double srp_cdm_std;
    while ( found == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s %s %lf %s", text, dummy, &srp_cdm_std, dummy);
    if (strcmp(text, "CSRP_SRP") == 0){
      OPTIONS->srp_cdm_std[isc] = sqrt(srp_cdm_std);
      found = 1;
    }
    }


     isc = 1;
     
    // BC
     found = 0;
    while ( found == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s %s %lf %s", text, dummy, &bc_cdm, dummy);
    if (strcmp(text, "CD_AREA_OVER_MASS") == 0){
      OPTIONS->bc_cdm[isc] = bc_cdm;
      found = 1;
    }
    }

    // SRP
    rewind(fp);
    found = 0;
    while ( found < 2 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s %s %lf %s", text, dummy, &srp_cdm, dummy);
    if (strcmp(text, "CR_AREA_OVER_MASS") == 0){
      OPTIONS->srp_cdm[isc] = srp_cdm;
      found = found + 1;
    }
    }

    // VARIANCE BC
    rewind(fp);
    found = 0;
    while ( found < 2 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s %s %lf %s", text, dummy, &bc_cdm_std, dummy);
    if (strcmp(text, "CDRG_DRG") == 0){
      OPTIONS->bc_cdm_std[isc] = sqrt(bc_cdm_std);
      found = found + 1;
    }
    }

    // VARIANCE SRP
    rewind(fp);
    found = 0;
    while ( found < 2 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s %s %lf %s", text, dummy, &srp_cdm_std, dummy);
    if (strcmp(text, "CSRP_SRP") == 0){
      OPTIONS->srp_cdm_std[isc] = sqrt(srp_cdm_std);
      found = found + 1;
    }
    }


/*     isc = 0; */
/*     printf("bc[%d] = %f m2/kg, variance = %e m**4/kg**2\n", isc, OPTIONS->bc_cdm[isc], OPTIONS->bc_cdm_std[isc] * OPTIONS->bc_cdm_std[isc] ); */
/*     printf("srp[%d] = %f m2/kg, variance = %e m**4/kg**2\n", isc, OPTIONS->srp_cdm[isc], OPTIONS->srp_cdm_std[isc] * OPTIONS->srp_cdm_std[isc] ); */
/*     isc = 1; */
/*     printf("bc[%d] = %f m2/kg, variance = %e m**4/kg**2\n", isc, OPTIONS->bc_cdm[isc], OPTIONS->bc_cdm_std[isc] * OPTIONS->bc_cdm_std[isc] ); */
/*     printf("srp[%d] = %f m2/kg, variance = %e m**4/kg**2\n", isc, OPTIONS->srp_cdm[isc], OPTIONS->srp_cdm_std[isc] * OPTIONS->srp_cdm_std[isc] ); */

/*     exitf(); */

}

// read_vcm: reads collision vcm files.
/* assumptions: */
/* - position and velocity of object are in km and akm/s. */
/* - position and velocity of objects are in j2000 frame of reference */
/* - cvoaraine atrix of object is in equinoctial elements */
/* - ballistic coeff and solar radiation pressure coeff are in m2/kg  */
int read_vcm(char filename[1000], OPTIONS_T *OPTIONS, int isc){

    FILE *fp;
    char *line = NULL;
    size_t len = 0;

    fp = fopen(filename, "r");
    getline(&line, &len, fp);
    char *next, *start, *end;
    int found = 0;
    char text[500];
    char epoch[100];    char epoch1[100];     char epoch2[100];    
    strcpy(epoch, "");    strcpy(epoch1, "");     strcpy(epoch2, "");    
    double et;
    // EPOCH. way below is not the smartestway to do it. i could simlpy do like i do for the position, velocity etc
    while ( found == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    strcpy(text, "");
    strncat(text, &line[0], 13);
    if (strcmp(text, "<> EPOCH TIME") == 0){
      next = &line[0];
      start = strstr(next, "(UTC):")+7;
      end = strstr(next, " EPOCH REV");
      strncat(epoch1, start,  (int)(end) - (int)(start));
      //      printf("<%s><||> <%s> || %d %d\n", line, epoch1, (int)(start - next), (int)(end-next) );

      // remove extra day and month
      next = &epoch1[0];
      start = next;
      end = strstr(next, "(");
      strncat(epoch2, start,  (int)(end) - (int)(start));
      RemoveSpaces(epoch2);
      strncat(epoch, epoch2, 4);
      strcat(epoch, "-");
      strncat(epoch, epoch2+4, (int)(&epoch2[-1] - (&epoch2[0]+4)));
      //            printf("epoch 2 <%s>\n", epoch2);
      //    printf("epoch  <%s>\n", epoch);

      strcat(epoch, "T");

      next = &epoch1[0];
      start = strrchr(next, ')') + 1;
      end = &next[-1];
      strncat(epoch, start,  (int)(end) - (int)(start));
      RemoveSpaces(epoch);

            /* printf("XXXX <%s>\n", epoch); */

	    /* //      strcpy(epoch, "2017-305T05:11:6.478"); // !!!!!! */
	                str2et_c(epoch, &OPTIONS->et_vcm[isc]);

	    /* char epoch_check[300]; */
	    /* strcpy(epoch_check, ""); */
	    /* et2utc_c(et, "ISOC" ,6 ,255 , epoch_check); */
	    /* printf("<%s> VS <%s>\n", epoch, epoch_check); */

    found = 1;
    }
    }
    // POSITION
    rewind(fp); found = 0;
    char dummy[100];
    double r[3];
    while ( found == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    strcpy(text, "");
    strncat(text, &line[0], 10);
    if (strcmp(text, "<> J2K POS") == 0){
      sscanf(line, "%s %s %s %s %lf %lf %lf", dummy, dummy,  dummy,  dummy,  &OPTIONS->x_eci[isc], &OPTIONS->y_eci[isc], &OPTIONS->z_eci[isc]);
      //      v_print(r, "r");

    found = 1;
    }
    }

    // VELOCOTY
    rewind(fp); found = 0;
    double v[3];
    while ( found == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    strcpy(text, "");
    strncat(text, &line[0], 10);
    if (strcmp(text, "<> J2K VEL") == 0){
      sscanf(line, "%s %s %s %s %lf %lf %lf", dummy, dummy,  dummy,  dummy,  &OPTIONS->vx_eci[isc], &OPTIONS->vy_eci[isc], &OPTIONS->vz_eci[isc]);
      //      v_print(v, "v");

    found = 1;
    }
    }

    // BALLISTIC COEFFICIENT
    rewind(fp); found = 0;

    while ( found == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    strcpy(text, "");
    strncat(text, &line[0], 17);
    if (strcmp(text, "<> BALLISTIC COEF") == 0){
      sscanf(line, "%s %s %s %s %lf %s", dummy, dummy,  dummy,  dummy,  &OPTIONS->bc_vcm[isc], dummy); 
      //      printf("bc %e\n",bc);
    found = 1;
    }
    }

    // SOLAR RADIATION PRESSURE COEFFICIENT
    rewind(fp); found = 0;

    while ( found == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    strcpy(text, "");
    strncat(text, &line[0], 24);
    if (strcmp(text, "<> SOLAR RAD PRESS COEFF") == 0){
      sscanf(line, "%s %s %s %s %s %s %lf %s", dummy, dummy,  dummy,  dummy, dummy, dummy , &OPTIONS->srp_vcm[isc], dummy);
      //      printf("srp %e\n",srp);
    found = 1;
    }
    }
     

    // COVARIANCE MATRIX. ONLY THE FIRST 21 ELEMENTS (ALL EQUINOCITAL ELEMENTS) + BC/BC and AGOM/AGOM
    rewind(fp); found = 0;
    int i,j;
    for (i = 0; i < 9; i++){
      for (j = 0; j < 9; j++){
	OPTIONS->covariance_matrix_equinoctial[isc][i][j] = 0;
      }
    }

    for (i = 0; i < 8; i++){
      for (j = 0; j < 8; j++){
	OPTIONS->covariance_matrix_equinoctial_no_bdot[isc][i][j] = 0;
      }
    }


    for (i = 0; i < 6; i++){
      for (j = 0; j < 6; j++){
	OPTIONS->covariance_matrix_equinoctial_no_bdot[isc][i][j] = 0;
      }
    }

    while ( found == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    strcpy(text, "");
    strncat(text, &line[0], 20);
    
    if (strcmp(text, "<> COVARIANCE MATRIX") == 0){
    getline(&line, &len, fp);
//    printf("<%s>\n", line);
    sscanf(line, "%s %lf %lf %lf %lf %lf", dummy, &OPTIONS->covariance_matrix_equinoctial[isc][0][0], &OPTIONS->covariance_matrix_equinoctial[isc][1][0], &OPTIONS->covariance_matrix_equinoctial[isc][1][1], &OPTIONS->covariance_matrix_equinoctial[isc][2][0], &OPTIONS->covariance_matrix_equinoctial[isc][2][1]);
    getline(&line, &len, fp);
    sscanf(line, "%s %lf %lf %lf %lf %lf", dummy, &OPTIONS->covariance_matrix_equinoctial[isc][2][2], &OPTIONS->covariance_matrix_equinoctial[isc][3][0], &OPTIONS->covariance_matrix_equinoctial[isc][3][1], &OPTIONS->covariance_matrix_equinoctial[isc][3][2], &OPTIONS->covariance_matrix_equinoctial[isc][3][3]);
    getline(&line, &len, fp);
    sscanf(line, "%s %lf %lf %lf %lf %lf", dummy, &OPTIONS->covariance_matrix_equinoctial[isc][4][0], &OPTIONS->covariance_matrix_equinoctial[isc][4][1], &OPTIONS->covariance_matrix_equinoctial[isc][4][2], &OPTIONS->covariance_matrix_equinoctial[isc][4][3], &OPTIONS->covariance_matrix_equinoctial[isc][4][4]);
    getline(&line, &len, fp);
    sscanf(line, "%s %lf %lf %lf %lf %lf", dummy, &OPTIONS->covariance_matrix_equinoctial[isc][5][0], &OPTIONS->covariance_matrix_equinoctial[isc][5][1], &OPTIONS->covariance_matrix_equinoctial[isc][5][2], &OPTIONS->covariance_matrix_equinoctial[isc][5][3], &OPTIONS->covariance_matrix_equinoctial[isc][5][4]);

    getline(&line, &len, fp); // 21 to 25
    sscanf(line, "%s %lf %lf %lf %lf %lf", dummy, &OPTIONS->covariance_matrix_equinoctial[isc][5][5], &OPTIONS->covariance_matrix_equinoctial[isc][6][0], &OPTIONS->covariance_matrix_equinoctial[isc][6][1], &OPTIONS->covariance_matrix_equinoctial[isc][6][2], &OPTIONS->covariance_matrix_equinoctial[isc][6][3]);

    getline(&line, &len, fp); // 26 to 30
    sscanf(line, "%s %lf %lf %lf %lf %lf", dummy, &OPTIONS->covariance_matrix_equinoctial[isc][6][4], &OPTIONS->covariance_matrix_equinoctial[isc][6][5], &OPTIONS->covariance_matrix_equinoctial[isc][6][6], &OPTIONS->covariance_matrix_equinoctial[isc][7][0], &OPTIONS->covariance_matrix_equinoctial[isc][7][1]);

    getline(&line, &len, fp); // 31 to 35
    sscanf(line, "%s %lf %lf %lf %lf %lf", dummy, &OPTIONS->covariance_matrix_equinoctial[isc][7][2], &OPTIONS->covariance_matrix_equinoctial[isc][7][3], &OPTIONS->covariance_matrix_equinoctial[isc][7][4], &OPTIONS->covariance_matrix_equinoctial[isc][7][5], &OPTIONS->covariance_matrix_equinoctial[isc][7][6]);

    getline(&line, &len, fp); // 36 to 40
    sscanf(line, "%s %lf %lf %lf %lf %lf", dummy, &OPTIONS->covariance_matrix_equinoctial[isc][7][7], &OPTIONS->covariance_matrix_equinoctial[isc][8][0], &OPTIONS->covariance_matrix_equinoctial[isc][8][1], &OPTIONS->covariance_matrix_equinoctial[isc][8][2], &OPTIONS->covariance_matrix_equinoctial[isc][8][3]);

    getline(&line, &len, fp); // 41 to 45
    sscanf(line, "%s %lf %lf %lf %lf %lf", dummy, &OPTIONS->covariance_matrix_equinoctial[isc][8][4], &OPTIONS->covariance_matrix_equinoctial[isc][8][5], &OPTIONS->covariance_matrix_equinoctial[isc][8][6], &OPTIONS->covariance_matrix_equinoctial[isc][8][7], &OPTIONS->covariance_matrix_equinoctial[isc][8][8]);

/*     OPTIONS->bc_vcm_std[isc] = sqrt(OPTIONS->covariance_matrix_equinoctial[isc][6][6] * OPTIONS->bc_vcm[isc] * OPTIONS->bc_vcm[isc] ); */
/*     OPTIONS->srp_vcm_std[isc] = sqrt(OPTIONS->covariance_matrix_equinoctial[isc][8][8] * OPTIONS->srp_vcm[isc] * OPTIONS->srp_vcm[isc] ); */
    /* printf("\nsc %d\n", isc); */
    /* printf("bc %e %e\n", bc_vcm_std2, OPTIONS->bc_vcm_std[isc]); */
    /* printf("srp %e %e\n", srp_vcm_std2, OPTIONS->srp_vcm_std[isc]); */


    // cov is symmetric
    for (i = 0; i < 9; i++){
      for (j = 0; j < 9; j++){
	if (j > i){
	OPTIONS->covariance_matrix_equinoctial[isc][i][j] = OPTIONS->covariance_matrix_equinoctial[isc][j][i];
	}
	if ((i < 7) && (j < 7)){
	OPTIONS->covariance_matrix_equinoctial_no_bdot[isc][i][j] = OPTIONS->covariance_matrix_equinoctial[isc][i][j];
	}


	if ((i < 6) && (j < 6)){
	OPTIONS->covariance_matrix_equinoctial_only_rv[isc][i][j] = OPTIONS->covariance_matrix_equinoctial[isc][i][j];
	}
      
      }
    }

      for (j = 0; j < 7; j++){
	OPTIONS->covariance_matrix_equinoctial_no_bdot[isc][7][j] = OPTIONS->covariance_matrix_equinoctial[isc][8][j];
      }
    	OPTIONS->covariance_matrix_equinoctial_no_bdot[isc][7][7] = OPTIONS->covariance_matrix_equinoctial[isc][8][8];
    for (i = 0; i < 7; i++){
	OPTIONS->covariance_matrix_equinoctial_no_bdot[isc][i][7] = OPTIONS->covariance_matrix_equinoctial[isc][i][8];
}


			       //            m_print9_temp(OPTIONS->covariance_matrix_equinoctial[isc], "eqco 9*9");
			       //            m_print8_temp(OPTIONS->covariance_matrix_equinoctial_no_bdot[isc], "eqco");
			       //exitf();
    /*   int ielt = 0; */
    /*   while (ielt < 21){ */
    /* 	if (mod(ielt , 5) == 0){ */
    /* getline(&line, &len, fp); */
    /* sscanf(line, "%lf %lf %lf %lf %lf", &OPTIONS->covariance_matrix_equinoctial[isc][][], */
    /* 	} */
    /*   printf("srp %e\n",srp); */
    /*   ielt = ielt + 1; */
    
    found = 1;
    }
    }
     
        fclose(fp);

  return 0;
}


int read_thrust(OPTIONS_T *OPTIONS){ // if section #THRUST exists in the main input file, this functions reads the file that contains information about the external thrust applied to the sc
  FILE *thrust_file = NULL;
  thrust_file = fopen(OPTIONS->thrust_filename, "r");
  if (thrust_file == NULL){
    printf("!***!\nThe thrust file:\n%s\ncould not be found. The program will stop.\n!***!\n", OPTIONS->thrust_filename); MPI_Finalize(); exit(0);
  }

     char *line = NULL;


     char text[256];
     size_t len = 0;
     getline(&line, &len, thrust_file);
     sscanf(line, "%s", text);
     strcpy(OPTIONS->thrust_start, text );
     getline(&line, &len, thrust_file);
     sscanf(line, "%s", text);
     strcpy(OPTIONS->thrust_stop, text );
     getline(&line, &len, thrust_file);
     sscanf(line, "%lf %lf %lf", &OPTIONS->thrust_accel_lvlh[0], &OPTIONS->thrust_accel_lvlh[1], &OPTIONS->thrust_accel_lvlh[2]);
     OPTIONS->thrust_accel_lvlh[0] = OPTIONS->thrust_accel_lvlh[0] / 1000.;// m/s2 to km/s2
     OPTIONS->thrust_accel_lvlh[1] = OPTIONS->thrust_accel_lvlh[1] / 1000.;// m/s2 to km/s2
     OPTIONS->thrust_accel_lvlh[2] = OPTIONS->thrust_accel_lvlh[2] / 1000.;// m/s2 to km/s2
     
     str2et_c(OPTIONS->thrust_start, &OPTIONS->et_thrust_start);
     str2et_c(OPTIONS->thrust_stop, &OPTIONS->et_thrust_stop); 

  

  fclose(thrust_file);
  


  return 0;

}
