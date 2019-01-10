#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include "time.h"
#include "fly_storm.h" 
#include "propagator.h" 
#include "moat_prototype.h"
#include "options.h"

int nProcs;
int iProc;

int debug_cbv = 0;// set this variable to 1 to print in the binary more variables than with debug_cbv = 0. It's different from the option iMin, it was added by cbv later (06-02-18). !!!!!!! MUST be the same as in find_specular_points.c

int main(int argc, char * argv[]) {

  // Declarations
  int ierr;

  OPTIONS_T       OPTIONS;
  PARAMS_T        PARAMS;
  char filename_input[300];
  char filename_storm_interpolated[300][N_STORM];
  int iDebugLevel = -1; // for now set it to -1
    char filename_input_no_path[300];
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &iProc);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
  if (iProc == 0) {
     printf("Interpolating the positions of the satellites...\n"); //flying over the storms.\n");
  }


  //  Load Options
  strcpy(filename_input, "./input/main_input/");
  strcat(filename_input, argv[1]);
    strcpy(filename_input_no_path, "");

    strcat(filename_input_no_path, argv[1]);
  load_options( &OPTIONS, filename_input, iProc, nProcs,  iDebugLevel, filename_input_no_path);
 

  // Load Params
  //newstructure
  load_params( &PARAMS,  iDebugLevel, OPTIONS.earth_fixed_frame, OPTIONS.use_ap_hist, iProc , OPTIONS.path_to_spice);
  //    load_params( &PARAMS, main_directory_location, iDebugLevel, OPTIONS.earth_fixed_frame, OPTIONS.use_ap_hist, iProc );
    //newstructure

  

  int compute_coverage = 0;
/*   if ( strcmp(argv[2], "1") == 0 ){ */
/*     compute_coverage = 1; */
/*   } */

  int iStart, iEnd, nCygEachPe, nCygLeft,i, iCygnss;
  char filenamella[300];
  char filenameecef[300];
  nCygEachPe = (OPTIONS.n_satellites - OPTIONS.nb_gps)/nProcs;
  nCygLeft = (OPTIONS.n_satellites - OPTIONS.nb_gps) - (nCygEachPe * nProcs);
  iStart = 0;
  for (i=0; i<iProc; i++) {
    iStart += nCygEachPe;
    if (i < nCygLeft && iProc > 0) iStart++;
  }
  iEnd = iStart+nCygEachPe;
  if (iProc < nCygLeft) iEnd++;
  for (iCygnss=iStart; iCygnss<iEnd; iCygnss++) {
    // Time linear interpolation of the positions of the satellites (ECEF and LLA). USED FOR VISUALIZATION
    strcpy(filenameecef, OPTIONS.dir_output_run_name_sat_name[iCygnss]);
    strcat(filenameecef, "/ECEF_");
    strcat(filenameecef, OPTIONS.filename_output[iCygnss]);
    strcpy(filenamella, OPTIONS.dir_output_run_name_sat_name[iCygnss]);
    strcat(filenamella, "/LLA_");
    strcat(filenamella, OPTIONS.filename_output[iCygnss]);
    //    printf("%d <%s> <%s>\n", iCygnss, filenamella, filenameecef);
    interpolate_position( filenamella, filenameecef,  OPTIONS.dt, 1.0, iCygnss, &OPTIONS, &PARAMS); // one second interpolation
      // End of time linear interpolation of the positions of the satellites (ECEF and LLA). USED FOR VISUALIZATION   
    // Read the bin files created by find_specular_points.c and compute the coverage of the specular points in the storm
    fly_storm(&OPTIONS, filename_storm_interpolated, iCygnss, compute_coverage, &PARAMS);
  }

  /* // Notify Exit */
  if (iProc == 0) {
     printf("Done interpolating the positions of the satellites.\n"); //flying over the storms.\n");
  }

  ierr = MPI_Finalize();    

  return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           fly_storm
//  Purpose:        Computes when a spacecraft fly over a tropcial storm
//  Assumptions:    Only one storm
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 10/01/2015    |   ---     | Initial implementation
//
/////////////////////////////////////////////////////////////////////////////////////////
int fly_storm( OPTIONS_T *OPTIONS,   char filename_storm_interpolated[300][N_STORM], int iCygnss, int compute_coverage, PARAMS_T *PARAMS)
{

  /* Declarations */
  //      char tstorm[300];
    char dist_to_name_storm[300], in_name_storm[300];
  double storm_ecef[3];
  int specular_in_storm ;
  //  double lat_storm, lon_storm, storm_radius_uncertainty, storm_34_radius;
  //  double ecef_x_storm, ecef_y_storm, ecef_z_storm;
  //  double et_storm;
  //  char time_storm[300];
  char *line = NULL;
  size_t len = 0;
  double distance_specular_point_to_center_of_storm;
  //  FILE *file_storm_interpolated;
  FILE *file_specular_position_in;
  char filename_specular_position_in[300];
  char filename_specular_position_out[300];
  FILE *file_specular_position_out;

  /* Algorithm */

  /* if (compute_coverage == 1){ */
  /*   file_storm_interpolated = fopen(filename_storm_interpolated[0], "r"); */
  /*   getline(&line,&len,file_storm_interpolated); */
  /*   sscanf(line, "%19[^\n] %lf %lf %lf %lf %lf %lf", time_storm, &lat_storm, &lon_storm, &ecef_x_storm, &ecef_y_storm, &ecef_z_storm, &storm_radius_uncertainty); */
  /*   str2et_c(time_storm, &et_storm); */
  /* } */
  // Read the storm files
  // // Define variables

  int istorm;
  FILE *file_storm;
  char ***storm_time,  **storm_basin, **storm_number;
  double **storm_forecast_period,**storm_lat, **storm_lon, **storm_cone_radius, **storm_34_radius; // storm_34_radius is the max of the four quadrants given in the NOAA storm file
 double **et_time_plus_forecast_period;
 int **wind_radii_code;
 int *nb_lines_in_storm_file;
  int max_nb_lines_storm_file = 45; // 45 because it should be 30 but we take a margin of 15 lines. 30 because 10 * 3. '10' because there are 10 possible forecast periods: 0 3 6 12 24 36 48 72 96 120. '3' because for each forecast period, there are maximum 3 lines (34, 50 and 64 kts wind radii)
  
  int iline;

  // // Allocate memory
 if (compute_coverage == 1){
   nb_lines_in_storm_file = malloc( OPTIONS->nb_storm * sizeof(int) );
  storm_basin  = malloc( OPTIONS->nb_storm * sizeof(char *) );
  storm_number  = malloc( OPTIONS->nb_storm * sizeof(char *) );
  storm_time  = malloc( OPTIONS->nb_storm * sizeof(char **) );
  storm_forecast_period  = malloc( OPTIONS->nb_storm * sizeof(double *) );
  et_time_plus_forecast_period  = malloc( OPTIONS->nb_storm * sizeof(double *) );
  storm_lat  = malloc( OPTIONS->nb_storm * sizeof(double *) );
  storm_lon  = malloc( OPTIONS->nb_storm * sizeof(double *) );
  storm_cone_radius  = malloc( OPTIONS->nb_storm * sizeof(double *) );
  storm_34_radius  = malloc( OPTIONS->nb_storm * sizeof(double *) );
  wind_radii_code  = malloc( OPTIONS->nb_storm * sizeof(int *) );
  for (istorm = 0; istorm < OPTIONS->nb_storm; istorm++){
    storm_basin[istorm] = malloc( 4 * sizeof(char));
    storm_number[istorm] = malloc( 4 * sizeof(char));
    storm_time[istorm] = malloc( max_nb_lines_storm_file * sizeof(char *));
    storm_forecast_period[istorm] = malloc( max_nb_lines_storm_file * sizeof(double));
    et_time_plus_forecast_period[istorm] = malloc( max_nb_lines_storm_file * sizeof(double));
    for ( iline = 0; iline < max_nb_lines_storm_file; iline++ ){
      storm_time[istorm][iline] = malloc( 16 * sizeof(char));
    }
    storm_lat[istorm] = malloc( max_nb_lines_storm_file * sizeof(double));
    storm_lon[istorm] = malloc( max_nb_lines_storm_file * sizeof(double));
    storm_cone_radius[istorm] = malloc( max_nb_lines_storm_file * sizeof(double));
    storm_34_radius[istorm] = malloc( max_nb_lines_storm_file * sizeof(double));
    wind_radii_code[istorm] = malloc( max_nb_lines_storm_file * sizeof(int));
  }
  if ( ( storm_basin == NULL ) ||   ( storm_number == NULL ) ||  ( storm_time == NULL ) || ( storm_forecast_period == NULL ) || ( et_time_plus_forecast_period == NULL ) ||   ( storm_lat == NULL ) ||   ( storm_lon == NULL ) ||   ( storm_cone_radius == NULL ) ||   ( storm_34_radius == NULL ) ||   ( wind_radii_code == NULL ) || ( nb_lines_in_storm_file == NULL) ){
    printf("***! (fly_storm) Could not allow memory space for at least one of the variables of the storm files \n. The program will stop. !***\n");
    MPI_Finalize();
    exit(0);
  } 
  char filename_storm_with_path[1200] ;
  //  int dumb;
  char dumb_str[10]; // dumb_str: for the comma

 double radius_first_quadrant, radius_second_quadrant, radius_third_quadrant, radius_fourth_quadrant;
	double et_storm_time; 	char storm_time_format_ok[16];
	//	char test_test[300], test_test2[300];
  // // read positions and radii of each storm
  for (istorm = 0; istorm < OPTIONS->nb_storm; istorm++){
    nb_lines_in_storm_file[istorm] = 0;
    strcpy(filename_storm_with_path, OPTIONS->dir_input_coverage_storm);
    strcat(filename_storm_with_path, "/");
    strcat(filename_storm_with_path, OPTIONS->filename_storm[istorm]);
  file_storm = fopen(filename_storm_with_path, "r");
  iline = 0;
  while( !feof(file_storm) ){
        getline(&line, &len, file_storm);
	if (!feof(file_storm)){
	  nb_lines_in_storm_file[istorm] = nb_lines_in_storm_file[istorm] + 1;
	sscanf(line, "%s %s %s %lf, %lf, %lf, %lf, %s %s %s %d, %s %lf, %lf, %lf, %lf", storm_basin[istorm], storm_number[istorm],  storm_time[istorm][iline], &storm_forecast_period[istorm][iline], &storm_lat[istorm][iline],  &storm_lon[istorm][iline], &storm_cone_radius[istorm][iline], dumb_str, dumb_str,dumb_str, &wind_radii_code[istorm][iline],dumb_str, &radius_first_quadrant, &radius_second_quadrant, &radius_third_quadrant, &radius_fourth_quadrant );
	strtok(storm_basin[istorm], ","); 
	strtok(storm_number[istorm], ","); 
	strtok(storm_time[istorm][iline], ",");

	// // // Convert time in seconds past J2000
	// // // // Convert storm_time to format to be read by str2et_c
	strcpy(storm_time_format_ok, "");
	strncat(storm_time_format_ok, &storm_time[istorm][iline][0],4);
	strcat(storm_time_format_ok, "-");
	strncat(storm_time_format_ok, &storm_time[istorm][iline][0]+4,2);
	strcat(storm_time_format_ok, "-");
	strncat(storm_time_format_ok, &storm_time[istorm][iline][0]+6,2);
	strcat(storm_time_format_ok, "T");
	strncat(storm_time_format_ok, &storm_time[istorm][iline][0]+8,2);
	strcat(storm_time_format_ok, ":00");
	str2et_c(storm_time_format_ok, &et_storm_time);

	// // // // //  just to check (can comment the 2 lines below)
	/* et2utc_c(et_storm_time, "ISOC", 0, 255, test_test);  */
	/*   //printf("<%s> -> <%s>\n",storm_time[istorm][iline], test_test);  */

	// // // Add this time to the forecast period (in seconds) to get et_time_plus_forecast_period
	et_time_plus_forecast_period[istorm][iline] = et_storm_time + storm_forecast_period[istorm][iline] * 3600.;
	// // // // //  just to check (can comment the 2 lines below)
	/* et2utc_c(et_time_plus_forecast_period[istorm][iline], "ISOC", 0, 255, test_test2); */
	/* printf("<%s> + %f hours: <%s>\n",test_test, storm_forecast_period[istorm][iline],test_test2); */

	//	printf("\nbasin: <%s> \nnumber: <%s> \ntime: <%s> \nperiod: <%f> \nlat: <%f> \nlon: <%f> \nradius cone: <%f> \nquadrants: <%f> <%f> <%f> <%f>\n", storm_basin[istorm], storm_number[istorm], storm_time[istorm][iline], storm_forecast_period[istorm][iline], storm_lat[istorm][iline], storm_lon[istorm][iline], storm_cone_radius[istorm][iline], radius_first_quadrant, radius_second_quadrant, radius_third_quadrant, radius_fourth_quadrant);

	// storm_34_radius[istorm][iline]  is the max of the 4 quadrant radii
	storm_34_radius[istorm][iline]  = radius_first_quadrant;
	if ( radius_second_quadrant > storm_34_radius[istorm][iline] ){
	  storm_34_radius[istorm][iline] = radius_second_quadrant;
	}
	if ( radius_third_quadrant > storm_34_radius[istorm][iline] ){
	  storm_34_radius[istorm][iline] = radius_third_quadrant;
	}	
	if ( radius_fourth_quadrant > storm_34_radius[istorm][iline] ){
	  storm_34_radius[istorm][iline] = radius_fourth_quadrant;
	}
	storm_34_radius[istorm][iline] =  storm_34_radius[istorm][iline] * 1.852; // nautic miles to km
	storm_cone_radius[istorm][iline] = storm_cone_radius[istorm][iline] * 1.852; // nautic miles to km 
	// printf("radius: <%f>\n", storm_34_radius[istorm][iline]);
	iline = iline + 1;
	}
  }

  fclose(file_storm);
 }
 }

  // read specular file
  strcpy(filename_specular_position_in, OPTIONS->dir_output_run_name_sat_name[iCygnss]);
  strcat(filename_specular_position_in, "/specular_");
  char *next; int find_file_name;      char sat_name[300];
  next = &OPTIONS->filename_output[iCygnss][0];
  find_file_name =  (int)(strchr(next, '.') - next);
  strncat(filename_specular_position_in, next, find_file_name);
  strcat(filename_specular_position_in, ".bin");
  file_specular_position_in = fopen(filename_specular_position_in, "r");

  // create output file 
  if ( compute_coverage == 1){
    strcpy(filename_specular_position_out, OPTIONS->dir_output_run_name_coverage_storm);
    strcat(filename_specular_position_out, "/coverage_specular_");
    strncat(filename_specular_position_out, next, find_file_name);
    strcpy(sat_name, "");
    strncat(sat_name, next, find_file_name);
    strcat(filename_specular_position_out, ".txt");
    file_specular_position_out = fopen(filename_specular_position_out, "w+");
    fprintf(file_specular_position_out, "This file shows the positions of the specular points for %s, as well as the distance from them to the storms.\nTIME LON LAT GAIN NAME_GPS ", sat_name);
    //    fprintf(file_specular_position_out, "This file shows the positions of the specular points for %s, as well as the distance from them to the storms.\nTIME ECEF_X ECEF_Y ECEF_Z LON LAT GAIN NAME_GPS ", sat_name);
  for (istorm = 0; istorm < OPTIONS->nb_storm; istorm++){

    strcpy(dist_to_name_storm, "DIST_"); 
    strcat(dist_to_name_storm, storm_basin[istorm]);
    strcat(dist_to_name_storm, storm_number[istorm]);    
    strcpy(in_name_storm, "IN_");
    strcat(in_name_storm, storm_basin[istorm]);
    strcat(in_name_storm, storm_number[istorm]);
    fprintf(file_specular_position_out, "%s %s ", dist_to_name_storm, in_name_storm);
  }
  fprintf(file_specular_position_out, "\n#START\n");
  }
  else{
    strcpy(filename_specular_position_out, OPTIONS->dir_output_run_name_sat_name[iCygnss]);
    strcat(filename_specular_position_out, "/specular_");
    strncat(filename_specular_position_out, next, find_file_name);
    strcpy(sat_name, "");
    strncat(sat_name, next, find_file_name);
    strcat(filename_specular_position_out, ".txt");
    file_specular_position_out = fopen(filename_specular_position_out, "w+");
           fprintf(file_specular_position_out, "This file shows the positions of the specular points for %s.\nNote: there is also the possibility of calculating the distance from each specular point to a storm.\nTIME LON_SPEC LAT_SPEC GAIN NAME_GPS NORM_POWER ECEF_SC_X ECEF_SC_Y ECEF_SC_Z ECEF_GPS_X ECEF_GPS_Y ECEF_GPS_Z ECEF_SPEC_X ECEF_SPEC_Y ECEF_SPEC_Z\n#START\n", sat_name);
	   //       fprintf(file_specular_position_out, "This file was generated by SpOCK. It shows the positions of the specular points for %s.\nFor questions, email cbv@umich.edu.\nTIME LON_SPEC LAT_SPEC\n#START\n", sat_name);

    //     fprintf(file_specular_position_out, "This file shows the positions of the specular points for %s.\nNote: there is also the possibility of calculating the distance from each specular point to a storm.\nTIME ECEF_SPEC_X ECEF_SPEC_Y ECEF_SPEC_Z LON_SPEC LAT_SPEC GAIN NAME_GPS ECEF_GPS_X ECEF_GPS_Y ECEF_GPS_Z ECEF_GPS_VX ECEF_GPS_VY ECEF_GPS_VZ ECEF_CYGNSS_X ECEF_CYGNSS_Y ECEF_CYGNSS_Z ECEF_CYGNSS_VX ECEF_CYGNSS_VY ECEF_CYGNSS_VZ\n#START\n", sat_name);
  }

  // Read the bin file created by find_specular_points.c
  int iGps, iPt;
    int iPtInner;
    float lon_spec, lat_spec,  ecef_x_spec, ecef_y_spec, ecef_z_spec, time_ymdhmsm[7];
    int8_t gain_spec,  normpower_spec;  // onboard algoirthm computes RCG as integer. normpower_spec (or gain_spec) represents RCG.
    int sss;
    char year_str[15], month_str[15], day_str[15], hour_str[15], minute_str[15], second_str[15], time_str[120];
    double et_spec, et_spec_save; char time_spec[300];
    float lat_sat, lon_sat, heading_sat;
  while(!feof(file_specular_position_in)){

    fread(&iGps,sizeof(iGps),1,file_specular_position_in);
    fread(&iPt,sizeof(iPt),1,file_specular_position_in);
    iPtInner = ( iPt % (int) (OPTIONS->dt_output) );

    if (iPtInner == 0){
      for (sss = 0; sss < 7; sss ++){
	fread(&time_ymdhmsm[sss],sizeof(time_ymdhmsm[sss]),1,file_specular_position_in);
      }
      sprintf(year_str, "%d", (int)(time_ymdhmsm[0]));
      sprintf(month_str, "%d", (int)(time_ymdhmsm[1]));
      sprintf(day_str, "%d", (int)(time_ymdhmsm[2]));
      sprintf(hour_str, "%d", (int)(time_ymdhmsm[3]));
      sprintf(minute_str, "%d", (int)(time_ymdhmsm[4]));
      sprintf(second_str, "%d", (int)(time_ymdhmsm[5]));
      strcpy(time_str, "");
      strcat(time_str, year_str); strcat(time_str, "-"); strcat(time_str, month_str); strcat(time_str, "-"); strcat(time_str, day_str); strcat(time_str, "T"); strcat(time_str, hour_str); strcat(time_str, ":"); strcat(time_str, minute_str); strcat(time_str, ":"); strcat(time_str, second_str);
      str2et_c(time_str, &et_spec);
      et_spec_save = et_spec;
    }
    else{
      et_spec = et_spec_save + iPtInner; // iPtInner always increments by one second in find_specular_points.c
    }

    et2utc_c(et_spec, "ISOC" ,0 ,255 , time_spec);
    fread(&lon_sat,sizeof(lon_sat),1,file_specular_position_in);
    fread(&lat_sat,sizeof(lat_sat),1,file_specular_position_in);
    fread(&heading_sat,sizeof(heading_sat),1,file_specular_position_in);

	      if ( debug_cbv == 1){
    fread(&ecef_x_spec,sizeof(ecef_x_spec),1,file_specular_position_in);
    fread(&ecef_y_spec,sizeof(ecef_y_spec),1,file_specular_position_in);
    fread(&ecef_z_spec,sizeof(ecef_z_spec),1,file_specular_position_in);
	      }
    // !!!!!!!!!!!!!!!! COMMENT THESE LINES BELOW!
    /* float ecef_x_gps, ecef_y_gps, ecef_z_gps, ecef_x_sat, ecef_y_sat, ecef_z_sat; */
    /* fread(&ecef_x_gps,sizeof(ecef_x_gps),1,file_specular_position_in); */
    /* fread(&ecef_y_gps,sizeof(ecef_y_gps),1,file_specular_position_in); */
    /* fread(&ecef_z_gps,sizeof(ecef_z_gps),1,file_specular_position_in); */
    /* fread(&ecef_x_sat,sizeof(ecef_x_sat),1,file_specular_position_in); */
    /* fread(&ecef_y_sat,sizeof(ecef_y_sat),1,file_specular_position_in); */
    /* fread(&ecef_z_sat,sizeof(ecef_z_sat),1,file_specular_position_in); */
    // !!!!!!!!!!!!!!!! END OF COMMENT THESE LINES BELOW!
    //    printf("AAAAAAAAAAAAAAAAAAAAAAA\n");
    fread(&lon_spec,sizeof(lon_spec),1,file_specular_position_in);
    fread(&lat_spec,sizeof(lat_spec),1,file_specular_position_in);
    fread(&gain_spec,sizeof(gain_spec),1,file_specular_position_in);

    // !!!!!!!!!!!!!!!! COMMENT THESE LINES BELOW!
    float ecef_x_gps, ecef_y_gps, ecef_z_gps, ecef_x_sat, ecef_y_sat, ecef_z_sat,ecef_vx_gps, ecef_vy_gps, ecef_vz_gps, ecef_vx_sat, ecef_vy_sat, ecef_vz_sat;
    //  if ( iPt % 60 == 0 ){
	      if ( debug_cbv == 1){
      fread(&ecef_x_gps,sizeof(ecef_x_gps),1,file_specular_position_in);
      fread(&ecef_y_gps,sizeof(ecef_y_gps),1,file_specular_position_in);
      fread(&ecef_z_gps,sizeof(ecef_z_gps),1,file_specular_position_in);
      fread(&ecef_vx_gps,sizeof(ecef_vx_gps),1,file_specular_position_in);
      fread(&ecef_vy_gps,sizeof(ecef_vy_gps),1,file_specular_position_in);
      fread(&ecef_vz_gps,sizeof(ecef_vz_gps),1,file_specular_position_in);
      fread(&ecef_x_sat,sizeof(ecef_x_sat),1,file_specular_position_in);
      fread(&ecef_y_sat,sizeof(ecef_y_sat),1,file_specular_position_in);
      fread(&ecef_z_sat,sizeof(ecef_z_sat),1,file_specular_position_in);
      fread(&ecef_vx_sat,sizeof(ecef_vx_sat),1,file_specular_position_in);
      fread(&ecef_vy_sat,sizeof(ecef_vy_sat),1,file_specular_position_in);
      fread(&ecef_vz_sat,sizeof(ecef_vz_sat),1,file_specular_position_in);
      fread(&normpower_spec,sizeof(normpower_spec),1,file_specular_position_in);
	      }
      int8_t elev_spec; // the on-board algorithm computes the azim and elev as integers
	int16_t azim_spec;// the on-board algorithm computes the azim and elev as integers
	float elev_gps_from_cyg;
	float azim_spec_not_int, elev_spec_not_int;
	      if ( debug_cbv == 1){
      fread(&elev_spec,sizeof(elev_spec),1,file_specular_position_in);
      fread(&azim_spec,sizeof(azim_spec),1,file_specular_position_in);
      fread(&elev_gps_from_cyg,sizeof(elev_gps_from_cyg),1,file_specular_position_in);
      fread(&elev_spec_not_int,sizeof(elev_spec_not_int),1,file_specular_position_in);
      fread(&azim_spec_not_int,sizeof(azim_spec_not_int),1,file_specular_position_in);

      if (azim_spec < 0){
	azim_spec = 360 + azim_spec;
      }


      if (azim_spec_not_int < 0){
	azim_spec_not_int = 360 + azim_spec_not_int;
      }

	      }
      // UNCOMMENT 6 LINES BELOW FOR ONE YEAR RUN FOR E2ES
     /*  if (!feof(file_specular_position_in)){ */
     /* if ( iPt % 60 == 0 ){ */
     /*  	fprintf(file_specular_position_out, " %f %f %f %f %f %f %f %f %f %f %f %f", ecef_x_gps, ecef_y_gps, ecef_z_gps,ecef_vx_gps, ecef_vy_gps, ecef_vz_gps, ecef_x_sat, ecef_y_sat, ecef_z_sat,ecef_vx_sat, ecef_vy_sat, ecef_vz_sat); */
     /* } */
     /*  } */
      ///}
    // !!!!!!!!!!!!!!!! END OF COMMENT THESE LINES BELOW!

      //      printf("BBBBBBBBBBBBBBBBBBBBBBBBBBBB\n");
    if (!feof(file_specular_position_in)){
      if ( debug_cbv == 1){

      fprintf(file_specular_position_out, "%s %f %f %d %s %d %f %f %f %f %f %f %f %f %f %d %d %f %f %f", time_spec, lon_spec, lat_spec, gain_spec, OPTIONS->gps_file_name[iGps], normpower_spec, ecef_x_sat, ecef_y_sat, ecef_z_sat, ecef_x_gps, ecef_y_gps, ecef_z_gps, ecef_x_spec, ecef_y_spec, ecef_z_spec, elev_spec, azim_spec, elev_gps_from_cyg, elev_spec_not_int, azim_spec_not_int);
      }
      else{
      fprintf(file_specular_position_out, "%s %f %f %d %s", time_spec, lon_spec, lat_spec, gain_spec, OPTIONS->gps_file_name[iGps]);
      }
      //   fprintf(file_specular_position_out, "%s %f %f", time_spec, lon_spec, lat_spec );
      //      fprintf(file_specular_position_out, "%s %f %f %f %f %f %f %s", time_spec, ecef_x_spec, ecef_y_spec, ecef_z_spec, lon_spec, lat_spec, gain_spec, OPTIONS->gps_file_name[iGps]);
    }

    //    printf("CCCCCCCCCCCCCCCCCCCCcc\n");

    // !!!!!!!!!!!!!!!! COMMENT THESE LINES BELOW! these are just when I want to look at the positions of the specular points with respect to the CYGNSS and GPS satellites
    /* if (!feof(file_specular_position_in)){ */
    /*   fprintf(file_specular_position_out, "%s %f %f %f %f %f %f %f %f %f", time_spec, ecef_x_spec, ecef_y_spec, ecef_z_spec, ecef_x_gps, ecef_y_gps, ecef_z_gps, ecef_x_sat, ecef_y_sat, ecef_z_sat); */
    /* } */
    // !!!!!!!!!!!!!!!! END OF COMMENT THESE LINES BELOW!

    if ( compute_coverage == 1 ){

  for (istorm = 0; istorm < OPTIONS->nb_storm; istorm++){
    /* et2utc_c(et_time_plus_forecast_period[istorm][0],"ISOC", 0, 255, tstorm); */
    /* printf("storm: <%s> || spec: <%s>\n", tstorm, time_spec); */
      iline= 0;
      if (  et_time_plus_forecast_period[istorm][0] - 0.0001 <= et_spec){
	while (et_time_plus_forecast_period[istorm][iline] - 0.0001 <= et_spec ){
	  iline = iline + 1;
	  if (iline >= nb_lines_in_storm_file[istorm]){
	    break;
	  }
	}
	if (iline < nb_lines_in_storm_file[istorm]){
	iline = iline - 1; // this is the line for which the forecast period is before the spec time and for which the next line is after the spec time

	// if the wind_radii_code = 50 or 64, we don't want to look at it so we go up to find the line with wind_radii_code = 34 at the same forecast period. If there is no line with 34 wind radius then we keep whatever line we have (50 or 64)
	if ( ( wind_radii_code[istorm][iline] != 34 )  && ( storm_forecast_period[istorm][iline-1] == storm_forecast_period[istorm][iline] ) ){
	  iline = iline-1;
	}
	// one more time in case now wind_radii_code[istorm][iline] = 50 (so oringally wind_radii_code[istorm][iline] was 64)
	if ( ( wind_radii_code[istorm][iline] != 34 )  && ( storm_forecast_period[istorm][iline-1] == storm_forecast_period[istorm][iline] ) ){
	  iline = iline-1;
	}

	//	printf("(%d) %s: %f %f %f\n", istorm, time_spec, storm_lat[istorm][iline], storm_lon[istorm][iline],storm_34_radius[istorm][iline] / 1.852);
	// convert the latitude/longitude of the storm at that time into ECEF coordinates
      geodetic_to_geocentric( PARAMS->EARTH.flattening, 0.0, storm_lat[istorm][iline]*DEG2RAD, storm_lon[istorm][iline]*DEG2RAD, PARAMS->EARTH.radius, storm_ecef);
      	distance_specular_point_to_center_of_storm = sqrt( ( ecef_x_spec - storm_ecef[0] ) * ( ecef_x_spec - storm_ecef[0]  )  +  ( ecef_y_spec - storm_ecef[1] ) * ( ecef_y_spec - storm_ecef[1]  ) +  ( ecef_z_spec - storm_ecef[2] ) * ( ecef_z_spec - storm_ecef[2] ) ) ;
	if ( distance_specular_point_to_center_of_storm <  storm_34_radius[istorm][iline] ){
	  specular_in_storm = 1;
	}
	else{
specular_in_storm = 0;
	}
    if (!feof(file_specular_position_in)){
	fprintf(file_specular_position_out, " %f %d", distance_specular_point_to_center_of_storm, specular_in_storm);
    }
	}
	else{
    if (!feof(file_specular_position_in)){
	fprintf(file_specular_position_out, " NO_DATA NO_DATA");
    }

	}
      }
      else{
    if (!feof(file_specular_position_in)){
	fprintf(file_specular_position_out, " NO_DATA NO_DATA");
    }
      }
  }
  fprintf(file_specular_position_out, "\n");   
 }
    else{
      fprintf(file_specular_position_out,"\n");
    }
  }
  fclose(file_specular_position_in);
  fclose(file_specular_position_out); 

 
  if (compute_coverage == 1){
    free(storm_basin); free(storm_forecast_period); free(storm_number); free(storm_lat); free(storm_lon); free(storm_cone_radius); free(storm_34_radius); free(storm_time); free(et_time_plus_forecast_period); free(wind_radii_code); free(nb_lines_in_storm_file);
  }

  return 0;

}


/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           initialize_fly_storm
//  Purpose:        Computes when a spacecraft fly over a tropcial storm
//  Assumptions:    There could be more than one storm, contrarily to fly_storm
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 10/01/2015    |   ---     | Initial implementation
//
/////////////////////////////////////////////////////////////////////////////////////////
int initialize_fly_storm( OPTIONS_T *OPTIONS, PARAMS_T *PARAMS,  char filename_storm_interpolated[300][N_STORM], char main_directory_location[300])

{
  /* Declarations */
  char text_filename_storm[300];
  double radius_uncertainty[2][8];
  double alt_test, lat_test,lon_test;
  double vect_ecef_test[3];
  double time_for_interpolation = 0;
  double storm_ecef_interpolated_x, storm_ecef_interpolated_y, storm_ecef_interpolated_z;
  double storm_ecef_x_previous_step = 0, storm_ecef_y_previous_step = 0, storm_ecef_z_previous_step = 0;
  double save_previous_et=0;
  int line_number[N_STORM];
  ssize_t read;
  char *line = NULL;
  int sss;
  size_t len = 0;
  char date_storm_string[300];
  int find_file_name;
  char *next;
  char str_useless[10];
  char text_output[300];
  char times[300]; 
  FILE *file_radius_uncertainty_interpolated;
  double number_day_prediction_storm_in_seconds;
  double time_for_interpolation_radius_uncertainty;
  double radius_uncertainty_interpolated;
  double et_radius;

  FILE *file_storm[N_STORM];
  FILE *file_storm_interpolated[N_STORM];
  char storm_name[300][N_STORM];
  double latitude[N_STORM]; 
  double longitude[N_STORM];
  double storm_ECEF[N_STORM][3];
  double et_storm[N_STORM];


  /* Algorithm */

  printf("Interpolating the storms files...\n");
  /* Radius of uncertainty */
  // initialize the radius of uncertainty
  // first column is the time 
  radius_uncertainty[0][0] = 0;
  radius_uncertainty[0][1] = 12;
  radius_uncertainty[0][2] = 24;
  radius_uncertainty[0][3] = 36;
  radius_uncertainty[0][4] = 48;
  radius_uncertainty[0][5] = 72;
  radius_uncertainty[0][6] = 96;
  radius_uncertainty[0][7] = 120;
  for (sss = 0; sss < 8; sss++){
    radius_uncertainty[0][sss] = radius_uncertainty[0][sss] * 3600.;
  }
  // second column is the uncertainty at the given time. "* 1.852" is to convert from nautical miles to km. "+ 300" is the radius of the hurricane (assumption here) // source for the uncertainty: http://www.nhc.noaa.gov/aboutcone.shtml
  radius_uncertainty[1][0] = 0 * 1.852 + 300.;
  radius_uncertainty[1][1] = 32 * 1.852 + 300.;
  radius_uncertainty[1][2] = 52 * 1.852 + 300.;
  radius_uncertainty[1][3] = 71 * 1.852 + 300.;
  radius_uncertainty[1][4] = 90 * 1.852 + 300.;
  radius_uncertainty[1][5] = 122 * 1.852 + 300.;
  radius_uncertainty[1][6] = 170 * 1.852 + 300.;
  radius_uncertainty[1][7] = 225 * 1.852 + 300.;
      
  // Now interpolate over the N days and write the result in a file 
  file_radius_uncertainty_interpolated = fopen("file_radius_uncertainty_interpolated.txt", "w+");
  number_day_prediction_storm_in_seconds = radius_uncertainty[0][7];
  for (sss = 0; sss < 7 ; sss++){
    time_for_interpolation_radius_uncertainty = radius_uncertainty[0][sss];
    while(time_for_interpolation_radius_uncertainty <= radius_uncertainty[0][sss+1]){
      radius_uncertainty_interpolated = radius_uncertainty[1][sss] + ( time_for_interpolation_radius_uncertainty - radius_uncertainty[0][sss] ) / ( radius_uncertainty[0][sss+1] - radius_uncertainty[0][sss] ) * ( radius_uncertainty[1][sss+1] - radius_uncertainty[1][sss] );
      fprintf(file_radius_uncertainty_interpolated, "%f %f \n", time_for_interpolation_radius_uncertainty, radius_uncertainty_interpolated);
      time_for_interpolation_radius_uncertainty = time_for_interpolation_radius_uncertainty + 1.0;
    }
  }

  /* STORMS */
  for (sss = 0; sss<OPTIONS->nb_storm; sss++){ 
    line_number[sss] = 0;
    rewind(file_radius_uncertainty_interpolated); 
    strcpy(text_filename_storm, OPTIONS->dir_input_coverage_storm);
    strcat(text_filename_storm, "/");
    strcat(text_filename_storm, OPTIONS->filename_storm[sss]);
    strtok(text_filename_storm, "\n");
    file_storm[sss] = fopen(text_filename_storm, "r");

    // interplation file 
    file_storm_interpolated[sss] = fopen(filename_storm_interpolated[sss], "w+");

    // name of the storm
    getline(&line, &len, file_storm[sss]);
    sscanf(line,"%s", storm_name[sss]);

    // month and year of the storm
    getline(&line, &len, file_storm[sss]);
    next = &line[0];
    strcpy(text_output," ");
    find_file_name =  (int)(strchr(next, ',') - next)-1;
    strncat(text_output, next, find_file_name+7);

    while ( (read = getline(&line, &len, file_storm[sss])) != -1 ) {

      // day, hour and minute of the storm
      sscanf(line, "%9[^\n] %lf %lf", str_useless, &latitude[sss], &longitude[sss]);
      strcpy(date_storm_string,"");
      strncat(date_storm_string,&str_useless[0],2);
      strncat(date_storm_string,&text_output[0],25);
      strncat(date_storm_string,&str_useless[4],3);
      strncat(date_storm_string,":",1);
      strncat(date_storm_string,&str_useless[7],2);

      // conversion to second past 2000
      str2et_c(date_storm_string, &et_storm[sss]);
    
      // conversion lat/lon/ to ECEF
      geodetic_to_geocentric( PARAMS->EARTH.flattening, 0.0, latitude[sss]*DEG2RAD, longitude[sss]*DEG2RAD, PARAMS->EARTH.radius, storm_ECEF[sss]);
      //  v_print(storm_ECEF[sss], "");

      if (line_number[sss] > 0){
  	time_for_interpolation = save_previous_et;
  	while ((int)(time_for_interpolation) <= (int)(et_storm[sss])){

  	  // read radius of uncertainty from the interpolated radius uncertainty file
  	  getline(&line, &len, file_radius_uncertainty_interpolated);
  	  sscanf(line, "%lf %lf", &et_radius, &radius_uncertainty_interpolated);

  	  storm_ecef_interpolated_x = storm_ecef_x_previous_step + ( ( time_for_interpolation - save_previous_et ) / ( et_storm[sss] - save_previous_et )  ) * ( storm_ECEF[sss][0] - storm_ecef_x_previous_step ) ;
  	  storm_ecef_interpolated_y = storm_ecef_y_previous_step + ( ( time_for_interpolation - save_previous_et ) / ( et_storm[sss] - save_previous_et ) ) * ( storm_ECEF[sss][1] - storm_ecef_y_previous_step ) ;
  	  storm_ecef_interpolated_z = storm_ecef_z_previous_step + ( ( time_for_interpolation - save_previous_et ) / ( et_storm[sss] - save_previous_et ) ) * ( storm_ECEF[sss][2] - storm_ecef_z_previous_step ) ;

  	  // write in the storm extrapolated file the interpolated positions
  	  et2utc_c(time_for_interpolation, "ISOC" ,0 ,255 , times);
  	  //	  fprintf(file_storm_interpolated[sss], "%s %f %f %f \n", times, storm_ecef_interpolated_x, storm_ecef_interpolated_y, storm_ecef_interpolated_z);

  	  time_for_interpolation = time_for_interpolation + 1.0; // the linear interpolation is every second because a bigger time step is not accurate enough (if one minute for instance, the satellites moves too fast anbd we could miss its fly over the storm)

  	  vect_ecef_test[0] = storm_ecef_interpolated_x;
  	  vect_ecef_test[1] = storm_ecef_interpolated_y;
  	  vect_ecef_test[2] = storm_ecef_interpolated_z;
	  geocentric_to_geodetic( vect_ecef_test,
				  &PARAMS->EARTH.radius,
				  &PARAMS->EARTH.flattening,
				  &alt_test,
				  &lat_test,
				  &lon_test);

	  fprintf(file_storm_interpolated[sss], "%s %f %f %f %f %f %f\n", times, lat_test*RAD2DEG, lon_test*RAD2DEG,  storm_ecef_interpolated_x, storm_ecef_interpolated_y, storm_ecef_interpolated_z, radius_uncertainty_interpolated);

  	}

      }

      save_previous_et = et_storm[sss];
      storm_ecef_x_previous_step = storm_ECEF[sss][0];
      storm_ecef_y_previous_step = storm_ECEF[sss][1];
      storm_ecef_z_previous_step = storm_ECEF[sss][2];

      line_number[sss] = line_number[sss] + 1;

    }

    free(line);

    fclose(file_radius_uncertainty_interpolated);
    fclose(file_storm[sss]);
    fclose(file_storm_interpolated[sss]);

  }

  return 0;

}

