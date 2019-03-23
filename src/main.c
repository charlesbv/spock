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
    char filename_input_raw[300];
    int ierr;
    CONSTELLATION_T *CONSTELLATION = malloc(sizeof(CONSTELLATION_T));

    ierr = MPI_Init(&argc, &argv);
     
    /* find out MY process ID, and how many processes were started. */
      
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &iProc);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
/*     char filename[1000]; */
/*     FILE *file = NULL; */
/*      strcpy(filename,"test_area_attitude_5deg_ascii.bin"); */
/*   file = fopen(filename, "r"); */
/*   double bla; */
/*     fread(&bla,sizeof(bla),1,file); */
/*     printf("%f\n", bla); */
/*     MPI_Finalize();exit(0); */
/*     char test_wget[1000]filename_specular_position_in; */
/*     strcpy(test_wget, "wget http://services.swpc.noaa.gov/text/45-day-ap-forecast.txt"); */
/*     system(test_wget); */
    
//    printf("enwofjoiewjfoiewjfew\neowifheiowjfewewj\nfeowfjhweoinfeiw\n");
    strcpy( filename_input, "./input/main_input/");
    strcat(filename_input, argv[1]);
    strcpy(filename_input_raw, "");
    strcat(filename_input_raw, argv[1]);
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

    load_options( &OPTIONS, filename_input, iProc, nProcs, iDebugLevel, filename_input_raw);

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
    int degree = (int)(OPTIONS.degree);
    load_params( &PARAMS,  iDebugLevel, OPTIONS.earth_fixed_frame, OPTIONS.use_ap_hist, iProc, OPTIONS.path_to_spice, degree, OPTIONS.gravity_map, OPTIONS.include_earth_pressure );
    //newstructure

/*     // !!!!!!!!!!!!!!REMOVE BLOCK BELOW */
    
/*      // EQNCPV ( ET, EPOCH, EQEL, RAPOL, DECPOL, STATE ) */
/*      char date_temp[100]; */
/*      double et_temp; */
/*      double r_eq[9], r_eq_raw[6]; */
/*      double mu = PARAMS.EARTH.GRAVITY.mu * 1e9; */
/*      //     strcpy(date_temp, "2000-12-15T16:58:50.208000"); r_eq_raw[0] = -0.0010197; r_eq_raw[1] = 0.0003038; r_eq_raw[2] = 228.5098015; r_eq_raw[3] = 0.0011110; r_eq_raw[4] = 1.1243593; r_eq_raw[5] = 0.2074336; // af ag meanlon n chi psi */
/*      //strcpy(date_temp, "2000-12-15T04:59:09.164000"); r_eq_raw[0] = -0.3624854; r_eq_raw[1] = 0.6441811; r_eq_raw[2] = 45.1201104; r_eq_raw[3] = 0.0001549; r_eq_raw[4] = -0.4194393; r_eq_raw[5] = -0.4308122; // af ag meanlon n chi psi */
/*      //     strcpy(date_temp, "2000-12-14T05:25:03.461000"); r_eq_raw[0] = -0.0002825 ; r_eq_raw[1] = -0.0002319; r_eq_raw[2] = 318.0245529; r_eq_raw[3] = 0.0011035; r_eq_raw[4] = -0.3137868; r_eq_raw[5] = -3.8490362; // af ag meanlon n chi psi */

/*     double m_eq[6][6]; //af ag meanlon rad n chi psi */
/*     //    strcpy(date_temp, "2000-12-15T16:58:50.208000"); r_eq_raw[0] = -0.0010197; r_eq_raw[1] = 0.0003038; r_eq_raw[2] = 228.5098015;r_eq_raw[3] = 0.0011110; r_eq_raw[4] = 1.1243593; r_eq_raw[5] = 0.2074336; */
/* /\*      m_eq[0][0] = 8.042040e-013; m_eq[0][1] = 7.419230e-013; m_eq[0][2] = -2.026230e-012; m_eq[0][3] = 3.787930e-016; m_eq[0][4] = -1.773020e-013; m_eq[0][5] = 2.483520e-013; *\/ */
/* /\*      m_eq[1][0] = 7.419230e-013; m_eq[1][1] = 2.190440e-012; m_eq[1][2] = -3.820340e-012; m_eq[1][3] = 1.190100e-015; m_eq[1][4] = -2.838440e-013; m_eq[1][5] = 3.679250e-013; *\/ */
/* /\*      m_eq[2][0] = -2.026230e-012; m_eq[2][1] = -3.820340e-012; m_eq[2][2] = 1.975700e-011; m_eq[2][3] = 9.354380e-015; m_eq[2][4] = 1.473860e-012; m_eq[2][5] = -3.015420e-012; *\/ */
/* /\*      m_eq[3][0] = 3.787930e-016; m_eq[3][1] = 1.190100e-015; m_eq[3][2] = 9.354380e-015; m_eq[3][3] = 1.562360e-017; m_eq[3][4] = -8.860930e-017; m_eq[3][5] = -4.569740e-016; *\/ */
/* /\*      m_eq[4][0] = -1.773020e-013; m_eq[4][1] = -2.838440e-013; m_eq[4][2] = 1.473860e-012; m_eq[4][3] = -8.860930e-017; m_eq[4][4] = 9.677970e-013; m_eq[4][5] = -7.230720e-013; *\/ */
/* /\*      m_eq[5][0] = 2.483520e-013; m_eq[5][1] = 3.679250e-013; m_eq[5][2] = -3.015420e-012; m_eq[5][3] = -4.569740e-016; m_eq[5][4] = -7.230720e-013; m_eq[5][5] = 1.841230e-012; *\/ */


/* // Vallado 30 */
/* /\* m_eq[0][0] = 1.307401e-013; m_eq[0][1] = 8.454837e-014; m_eq[0][2] = -1.025769e-013; m_eq[0][3] = 2.767180e-016; m_eq[0][4] = -2.243993e-014; m_eq[0][5] = -2.467241e-014; *\/ */
/* /\* m_eq[1][0] = 8.454837e-014; m_eq[1][1] = 8.622244e-014; m_eq[1][2] = -6.605814e-014; m_eq[1][3] = 2.301159e-016; m_eq[1][4] = -1.730851e-014; m_eq[1][5] = -1.900285e-014; *\/ */
/* /\* m_eq[2][0] = -1.025769e-013; m_eq[2][1] = -6.605814e-014; m_eq[2][2] = 1.021784e-013; m_eq[2][3] = -2.251061e-016; m_eq[2][4] = 3.247419e-014; m_eq[2][5] = 5.462808e-016; *\/ */
/* /\* m_eq[3][0] = 2.767180e-016; m_eq[3][1] = 2.301159e-016; m_eq[3][2] = -2.251061e-016; m_eq[3][3] = 7.173977e-019; m_eq[3][4] = -4.709651e-017; m_eq[3][5] = -5.180442e-017; *\/ */
/* /\* m_eq[4][0] = -2.243993e-014; m_eq[4][1] = -1.730851e-014; m_eq[4][2] = 3.247419e-014; m_eq[4][3] = -4.709651e-017; m_eq[4][4] = 1.685802e-014; m_eq[4][5] = -8.851445e-015; *\/ */
/* /\* m_eq[5][0] = -2.467241e-014; m_eq[5][1] = -1.900285e-014; m_eq[5][2] = 5.462808e-016; m_eq[5][3] = -5.180442e-017; m_eq[5][4] = -8.851445e-015; m_eq[5][5] = 2.129877e-014; *\/ */

/* /\*      double af = r_eq_raw[0]; double ag = r_eq_raw[1]; double lequin = r_eq_raw[2]*M_PI/180.;  *\/ */
/* /\* double nequin = r_eq_raw[3]; double chi = r_eq_raw[4]; double psi = r_eq_raw[5]; *\/ */

/* // Vallado 15       a = (mu / n^2)^(1/3)  */

/*  /\*    // !!!erase  belwo *\/ */
/*  /\*    double af, ag, lequin, nequin, chi, psi; *\/ */
/*  /\*    double fr = -1; *\/ */
/*  /\*    double rvec[3], vvec[3];   *\/ */
/*  /\*    rvec[0] = -605.7922166; rvec[1] = -5870.2295111; rvec[2] = 3493.0531990; vvec[0] = -1.568254290; vvec[1] = -3.702348910; vvec[2] = -6.479483950; *\/ */


/*  /\*    cart_to_equin( &af, &ag, &lequin, &nequin, &chi, &psi,  mu/1.e9,  fr, rvec, vvec); // mu in km3 *\/ */
/*  /\* printf("equin: %f %f %f %f %f %f\n", af, ag, chi, psi, lequin * 180./M_PI, nequin); *\/ */
/*  /\* exitf(); *\/ */
/*  /\*    // !!!end of errase  belwo *\/ */

/*     strcpy(date_temp, "2000-12-15T16:58:50.208000");      double af = 0.0010610; double ag = 0.0000800; double lequin = 69.4157838*M_PI/180.; */
/*     double nequin = sqrt(mu/pow(6860.7631490*1000, 3)); double chi = 0.8601197; double psi = 0.1586839; double fr = -1; */

/*     // n rad af ag chi psi meanlonM rad */
/*     m_eq[0][0] = 4.306971e-17; m_eq[0][1] = -6.624132e-15; m_eq[0][2] = -1.753261e-14; m_eq[0][3] = 4.779764e-17; m_eq[0][4] = -2.416461e-17; m_eq[0][5] = 7.156497e-17; */
/*     m_eq[1][0] = -6.624132e-15; m_eq[1][1] = 4.505630e-12; m_eq[1][2] = 1.756086e-12; m_eq[1][3] = -1.499603e-14; m_eq[1][4] = 7.462691e-15; m_eq[1][5] = -5.202264e-12; */
/*     m_eq[2][0] = -1.753261e-14; m_eq[2][1] = 1.756086e-12; m_eq[2][2] = 8.443544e-12; m_eq[2][3] = -2.467493e-14; m_eq[2][4] = 1.404855e-14; m_eq[2][5] = 1.909008e-12; */
/*     m_eq[3][0] = 4.779764e-17; m_eq[3][1] = -1.499603e-14; m_eq[3][2] = -2.467493e-14; m_eq[3][3] = 1.372028e-12; m_eq[3][4] = 1.069441e-13; m_eq[3][5] = -1.331001e-13; */
/*     m_eq[4][0] = -2.416461e-17; m_eq[4][1] = 7.462691e-15; m_eq[4][2] = 1.404855e-14; m_eq[4][3] = 1.069441e-13; m_eq[4][4] = 1.615648e-12; m_eq[4][5] = 1.550717e-12; */
/*     m_eq[5][0] = 7.156497e-17; m_eq[5][1] = -5.202264e-12; m_eq[5][2] = 1.909008e-12; m_eq[5][3] = -1.331001e-13; m_eq[5][4] = 1.550717e-12; m_eq[5][5] = 1.048619e-11; */

/* /\*     // a m af ag chi psi meanlonM rad *\/ */
/* /\*     m_eq[0][0] = 7.299847e+02; m_eq[0][1] = 2.727093e-05; m_eq[0][2] = 7.218012e-05; m_eq[0][3] = -1.967784e-07; m_eq[0][4] = 9.948340e-08; m_eq[0][5] = -2.946262e-07; *\/ */
/* /\*     m_eq[1][0] = 2.727093e-05; m_eq[1][1] = 4.505630e-12; m_eq[1][2] = 1.756086e-12; m_eq[1][3] = -1.499603e-14; m_eq[1][4] = 7.462691e-15; m_eq[1][5] = -5.202264e-12; *\/ */
/* /\*     m_eq[2][0] = 7.218012e-05; m_eq[2][1] = 1.756086e-12; m_eq[2][2] = 8.443544e-12; m_eq[2][3] = -2.467493e-14; m_eq[2][4] = 1.404855e-14; m_eq[2][5] = 1.909008e-12; *\/ */
/* /\*     m_eq[3][0] = -1.967784e-07; m_eq[3][1] = -1.499603e-14; m_eq[3][2] = -2.467493e-14; m_eq[3][3] = 1.372028e-12; m_eq[3][4] = 1.069441e-13; m_eq[3][5] = -1.331001e-13; *\/ */
/* /\*     m_eq[4][0] = 9.948340e-08; m_eq[4][1] = 7.462691e-15; m_eq[4][2] = 1.404855e-14; m_eq[4][3] = 1.069441e-13; m_eq[4][4] = 1.615648e-12; m_eq[4][5] = 1.550717e-12; *\/ */
/* /\*     m_eq[5][0] = -2.946262e-07; m_eq[5][1] = -5.202264e-12; m_eq[5][2] = 1.909008e-12; m_eq[5][3] = -1.331001e-13; m_eq[5][4] = 1.550717e-12; m_eq[5][5] = 1.048619e-11; *\/ */

/* /\*      double r_class[7]; *\/ */
/* /\*      r_class[0] = sma; r_class[1] = ecc; r_class[2] = inc*180/M_PI; r_class[3] = raan*180/M_PI; r_class[4] = arg_per*180/M_PI; r_class[5] = true_ano*180./M_PI, r_class[6] = mean_ano * 180./M_PI; *\/ */
/*      //m_x_v6bis(r_class, T_equin_2_class, r_eq_raw); */
/* /\*      v_print6(r_eq_raw, "equin"); *\/ */
/* /\*      v_print7(r_class, "class"); *\/ */

/* /\*      r_eq[0] = pow( mu / ( r_eq_raw[3]*r_eq_raw[3] ), 1/3.); // a = (mu / n^2)^(1/3) *\/ */
/* /\* r_eq[1] = r_eq_raw[1]; r_eq[2] = r_eq_raw[0]; r_eq[3] = r_eq_raw[2]*M_PI/180.; r_eq[4] = r_eq_raw[4]; r_eq[5] = r_eq_raw[5]; r_eq[6] = 0; r_eq[7] = r_eq_raw[3]; r_eq[8] = 0; // *\/ */

/* /\*      //     double r_out[6]; *\/ */
/* /\*      SpiceDouble        r_out[6]; *\/ */
/* /\*      eqncpv_c(et_temp, et_temp, r_eq, -M_PI/2, M_PI/2., r_out); *\/ */
/* /\*  printf("%f %f %f %f %f %f\n", r_out[0], r_out[1], r_out[2], r_out[3], r_out[4], r_out[5], r_out[6]); *\/ */


/*     char filename_vcm[1000]; */
/*     strcpy(filename_vcm, "data/41890_27434_20171105_043841/41890_20171101_051106.txt"); */
/*     int isc = 0; */
/*     read_vcm(filename_vcm, &OPTIONS, isc); */
/*     exitf(); */
/* double rvec[3], vvec[3]; */
/*  equin_to_cart( rvec,  vvec,   af,  ag,  lequin,  nequin,  chi,  psi, mu, fr); */
/*  v_print(rvec, "rvec"); */
/*  v_print(vvec, "vvec"); */
/*  double T_equin_to_cart[6][6]; */
/*  compute_T_deriv_equin_to_cart( T_equin_to_cart, af, ag, lequin, nequin, chi, psi, mu, fr); */

/* /\*  T_equin_to_cart[0][0] = -8.829808e-02; T_equin_to_cart[0][1] = -2.430153e+06; T_equin_to_cart[0][2] = 1.556428e+06; T_equin_to_cart[0][3] = -2.902619e+06; T_equin_to_cart[0][4] = -5.721404e+06; T_equin_to_cart[0][5] = -1.411582e+06; *\/ */
/* /\*  T_equin_to_cart[1][0] = -8.556234e-01; T_equin_to_cart[1][1] = -4.175372e+06; T_equin_to_cart[1][2] = 7.832891e+06; T_equin_to_cart[1][3] = -1.089294e+05; T_equin_to_cart[1][4] = 4.548598e+06; T_equin_to_cart[1][5] = -3.332477e+06; *\/ */
/* /\*  T_equin_to_cart[2][0] = 5.091348e-01; T_equin_to_cart[2][1] = -1.214958e+07; T_equin_to_cart[2][2] = 8.149317e+05; T_equin_to_cart[2][3] = -6.864554e+05; T_equin_to_cart[2][4] = 6.651869e+06; T_equin_to_cart[2][5] = -5.832170e+06; *\/ */
/* /\*  T_equin_to_cart[3][0] = 1.142915e-04; T_equin_to_cart[3][1] = 8.143420e+01; T_equin_to_cart[3][2] = -1.704747e+03; T_equin_to_cart[3][3] = 8.007979e+03; T_equin_to_cart[3][4] = -3.608485e+03; T_equin_to_cart[3][5] = 6.739326e+02; *\/ */
/* /\*  T_equin_to_cart[4][0] = 2.698205e-04; T_equin_to_cart[4][1] = 4.819935e+03; T_equin_to_cart[4][2] = -5.754917e+03; T_equin_to_cart[4][3] = -2.819928e+02; T_equin_to_cart[4][4] = -5.813753e+03; T_equin_to_cart[4][5] = 6.530521e+03; *\/ */
/* /\*  T_equin_to_cart[5][0] = 4.722131e-04; T_equin_to_cart[5][1] = -5.912878e+03; T_equin_to_cart[5][2] = -4.706736e+03; T_equin_to_cart[5][3] = -1.777072e+03; T_equin_to_cart[5][4] = 4.195328e+03; T_equin_to_cart[5][5] = -3.885957e+03; *\/ */

/*      double m_cart[6][6], m_cart_temp[6][6]; */
/*      m_x_m6bis( m_cart_temp, T_equin_to_cart, m_eq ); */
/*      double T_equin_to_cart_trans[6][6]; */
/*      m_trans6(   T_equin_to_cart_trans, T_equin_to_cart); */
/*      m_x_m6bis( m_cart, m_cart_temp, T_equin_to_cart_trans); */
/*      //m_print6(T_equin_to_cart,"T_equin_to_cart"); */
/* 	     m_print6(m_cart,"m_cart"); */

/*  exitf(); */
/* /\*      double T_equin_2_class_temp[6][6]; *\/ */
/* /\*      compute_T_deriv_equin_to_class( T_equin_2_class_temp,  af,  ag,  lequin,  nequin,  chi,  psi,   mu, fr); *\/ */

/* /\*      double m_class[6][6], m_class_temp[6][6]; *\/ */
/* /\*      double T_equin_2_class[6][6]; *\/ */
/* /\*      //     m_trans6(   T_equin_2_class, T_equin_2_class_temp); // !!!!!! *\/ */
/* /\*           m_copy6(   T_equin_2_class, T_equin_2_class_temp); // !!!!!! *\/ */
/* /\*      m_x_m6bis( m_class_temp, T_equin_2_class, m_eq ); *\/ */
/* /\*      double T_equin_2_class_trans[6][6]; *\/ */
/* /\*      m_trans6(   T_equin_2_class_trans, T_equin_2_class); *\/ */
/* /\*      m_x_m6bis( m_class, m_class_temp, T_equin_2_class_trans); *\/ */



/* /\*      // Vallado 03 *\/ */
/* /\* /\\*      m_class[0][0] = 1.215911e+001; m_class[0][1] = 8.212505e-007; m_class[0][2] = 1.988270e-007; m_class[0][3] = -1.526735e-007; m_class[0][4] = 1.159226e-003; m_class[0][5] = -1.158147e-003; *\\/ *\/ */
/* /\* /\\*      m_class[1][0] = 8.212505e-007; m_class[1][1] = 8.083254e-014; m_class[1][2] = 1.698441e-014; m_class[1][3] = -1.304184e-014; m_class[1][4] = 7.796741e-011; m_class[1][5] = -7.787492e-011; *\\/ *\/ */
/* /\* /\\*      m_class[2][0] = 1.988270e-007; m_class[2][1] = 1.698441e-014; m_class[2][2] = 1.040397e-014; m_class[2][3] = 5.668433e-015; m_class[2][4] = 2.215181e-011; m_class[2][5] = -2.212971e-011; *\\/ *\/ */
/* /\* /\\*      m_class[3][0] = -1.526735e-007; m_class[3][1] = -1.304184e-014; m_class[3][2] = 5.668433e-015; m_class[3][3] = 1.859767e-014; m_class[3][4] = -1.700668e-011; m_class[3][5] = 1.699277e-011; *\\/ *\/ */
/* /\* /\\*      m_class[4][0] = 1.159226e-003; m_class[4][1] = 7.796741e-011; m_class[4][2] = 2.215181e-011; m_class[4][3] = -1.700668e-011; m_class[4][4] = 1.202832e-007; m_class[4][5] = -1.201791e-007; *\\/ *\/ */
/* /\* /\\*      m_class[5][0] = -1.158147e-003; m_class[5][1] = -7.787492e-011; m_class[5][2] = -2.212971e-011; m_class[5][3] = 1.699277e-011; m_class[5][4] = -1.201791e-007; m_class[5][5] = 1.200752e-007; *\\/ *\/ */
/* /\*      // Vallado 15 *\/ */
/* /\* /\\*      m_class[0][0] = 7.299847e+02; m_class[0][1] = 3.262250e-05; m_class[0][2] = 1.988270e-07; m_class[0][3] = -1.526735e-07; m_class[0][4] = 6.571931e-02; m_class[0][5] = -6.571976e-02; *\\/ *\/ */
/* /\* /\\*      m_class[1][0] = 3.262250e-05; m_class[1][1] = 4.791316e-12; m_class[1][2] = 1.698441e-14; m_class[1][3] = -1.304184e-14; m_class[1][4] = 1.909371e-09; m_class[1][5] = -1.914428e-09; *\\/ *\/ */
/* /\* /\\*      m_class[2][0] = 1.988270e-07; m_class[2][1] = 1.698441e-14; m_class[2][2] = 1.821030e-12; m_class[2][3] = 1.857460e-13; m_class[2][4] = 2.217579e-11; m_class[2][5] = -2.216053e-11; *\\/ *\/ */
/* /\* /\\*      m_class[3][0] = -1.526735e-07; m_class[3][1] = -1.304184e-14; m_class[3][2] = 1.857460e-13; m_class[3][3] = 2.051627e-12; m_class[3][4] = -1.673598e-11; m_class[3][5] = 1.701644e-11; *\\/ *\/ */
/* /\* /\\*      m_class[4][0] = 6.571931e-02; m_class[4][1] = 1.909371e-09; m_class[4][2] = 2.217579e-11; m_class[4][3] = -1.673598e-11; m_class[4][4] = 7.206135e-06; m_class[4][5] = -7.203997e-06; *\\/ *\/ */
/* /\* /\\*      m_class[5][0] = -6.571976e-02; m_class[5][1] = -1.914428e-09; m_class[5][2] = -2.216053e-11; m_class[5][3] = 1.701644e-11; m_class[5][4] = -7.203997e-06; m_class[5][5] = 7.201868e-06; *\\/ *\/ */
/* /\*      double sma,  ecc,  inc,  raan,  arg_per,  mean_ano, true_ano; *\/ */
/* /\*      str2et_c(date_temp, &et_temp); *\/ */
/* /\*      equin_to_class( &sma,  &ecc,  &inc,  &raan,  &arg_per, &true_ano, &mean_ano,  af,  ag,  lequin,  nequin,  chi,  psi,   mu, et_temp, fr); *\/ */
/* /\*      printf("%f %f %f %f %f %f %f\n", sma,  ecc,  inc*180./M_PI,  raan*180./M_PI,  arg_per*180./M_PI,  true_ano*180./M_PI,  mean_ano * 180./M_PI); *\/ */
/* /\*      double T_class_to_cart_temp[6][6]; *\/ */
/* /\*      compute_T_deriv_class_to_cart(  T_class_to_cart_temp,  sma,  ecc,  inc,  raan,  arg_per,  true_ano,mu ); *\/ */
/* /\*      double T_class_to_cart[6][6]; *\/ */
/* /\*      m_trans6(T_class_to_cart, T_class_to_cart_temp); *\/ */

     
     
/* /\* /\\*      T_class_to_cart[0][0] = -8.829808e-002; T_class_to_cart[0][1] = 2.551516e+005; T_class_to_cart[0][2] = 3.435083e+006; T_class_to_cart[0][3] = 5.870230e+006; T_class_to_cart[0][4] = -1.409737e+006; T_class_to_cart[0][5] = -1.411582e+006; *\\/ *\/ */
/* /\* /\\*      T_class_to_cart[1][0] = -8.556234e-001; T_class_to_cart[1][1] = 2.472462e+006; T_class_to_cart[1][2] = -6.337403e+005; T_class_to_cart[1][3] = -6.057922e+005; T_class_to_cart[1][4] = -3.323832e+006; T_class_to_cart[1][5] = -3.332476e+006; *\\/ *\/ */
/* /\* /\\*      T_class_to_cart[2][0] = 5.091348e-001; T_class_to_cart[2][1] = -1.471227e+006; T_class_to_cart[2][2] = -4.692898e+005; T_class_to_cart[2][3] = 0.000000e+000; T_class_to_cart[2][4] = -5.830333e+006; T_class_to_cart[2][5] = -5.832170e+006; *\\/ *\/ */
/* /\* /\\*      T_class_to_cart[3][0] = 1.142915e-004; T_class_to_cart[3][1] = -1.269885e+003; T_class_to_cart[3][2] = -6.371951e+003; T_class_to_cart[3][3] = 3.702349e+003; T_class_to_cart[3][4] = 6.721174e+002; T_class_to_cart[3][5] = 6.739326e+002; *\\/ *\/ */
/* /\* /\\*      T_class_to_cart[4][0] = 2.698205e-004; T_class_to_cart[4][1] = -7.476387e+003; T_class_to_cart[4][2] = 1.175565e+003; T_class_to_cart[4][3] = -1.568254e+003; T_class_to_cart[4][4] = 6.524030e+003; T_class_to_cart[4][5] = 6.530521e+003; *\\/ *\/ */
/* /\* /\\*      T_class_to_cart[5][0] = 4.722131e-004; T_class_to_cart[5][1] = 8.010540e+002; T_class_to_cart[5][2] = 8.705152e+002; T_class_to_cart[5][3] = 0.000000e+000; T_class_to_cart[5][4] = -3.890477e+003; T_class_to_cart[5][5] = -3.885957e+003; *\\/ *\/ */

/* /\*      double m_cart[6][6], m_cart_temp[6][6]; *\/ */
/* /\*      m_x_m6bis( m_cart_temp, T_class_to_cart, m_class ); *\/ */
/* /\*      double T_class_to_cart_trans[6][6]; *\/ */
/* /\*      m_trans6(   T_class_to_cart_trans, T_class_to_cart); *\/ */
/* /\*      m_x_m6bis( m_cart, m_cart_temp, T_class_to_cart_trans); *\/ */

/* /\*      //m_print6(m_eq, "equinoctial"); *\/ */
/* /\*      //     m_print6(T_equin_2_class, "T_equin_2_class"); *\/ */
/* /\* /\\*      m_print6(T_equin_2_class_trans, "T_equin_2_class_trans"); *\\/ *\/ */
/* /\* //     m_print6(m_class, "classical"); *\/ */
/* /\*      //           m_print6(T_class_to_cart, "T_class_to_cart"); *\/ */
/* /\*      //   m_print6(T_class_to_cart_trans, "T_class_to_cart_trans"); *\/ */
/* /\*                    m_print6(m_cart, "cartesian"); *\/ */

/* /\*     print_test(); *\/ */
/* /\*     exitf(); *\/ */

/* /\*     // !!!!!!!!!!!!!!end of REMOVE BLOCK BELOW *\/ */



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

    //    exitf();
    //    printf("\nSSSSSSSSSSSS %d\n",OPTIONS.nb_time_steps);


    //  Create a 3d map of the gravitational potential derivatives dUdr, dUdlat, and dUdlong. These are then used in compute_gravity to compute the acceleration due to the Earth gravity

    /* if (OPTIONS.gravity_map == 1){ */
    /*   gravity_map(CONSTELLATION, PARAMS.EARTH.GRAVITY, degree, iProc); */
    /*   printf("Done building the 3D gravity map.\n"); */
    /* } */

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

/*  LocalWords:  eq
 */
