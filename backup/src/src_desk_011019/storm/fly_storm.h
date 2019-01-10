
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SpiceUsr.h"
#include "propagator.h" 
#include "moat_prototype.h"
#include "options.h"



// Prototypes


int initialize_fly_storm( OPTIONS_T *OPTIONS, PARAMS_T *PARAMS, char filename_storm_interpolated[100][N_STORM], char main_directory_location[256]);

int fly_storm( OPTIONS_T *OPTIONS,   char filename_storm_interpolated[100][N_STORM], int iCygnss, int compute_coverage, PARAMS_T *PARAMS);



//int initialize_fly_storm( int *nb_storm,   char filename_storm_interpolated[100][N_STORM]);
/* int geodetic_to_geocentric(double flattening,            /\* IN:     flattening parameter     *\/ */
/* 			   double h,                     /\* IN:  M  height above ellipsoid  *\/ */
/* 			   double lat,                   /\* IN:  r  geodetic latitude       *\/ */
/* 			   double longitude,              /\* IN:  r  longitude               *\/ */
/* 			   double equatorial_radius,       /\* IN: equatorial radius *\/ */
/* 			   double  R_ecef_2_cg_ECEF[3]);   /\* OUT:     vector in ECEF           *\/ */


/* int geocentric_to_geodetic( */
/* 			   double  R_pt_wrt_ecef_ECEF[3], /\* vector in ECEF           *\/ */
/* 			   double *semimajor_axis,        /\* planetary radius         *\/ */
/* 			   double *flattening,            /\* flattening parameter     *\/ */
/* 			   double *h,                     /\* height above ellipsoid  *\/ */
/* 			   double *lat,                   /\* geodetic latitude       *\/ */
/* 			   double *longitude   );         /\* longitude               *\/ */


