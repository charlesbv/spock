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

#include "kalman_9state_new_with_tau.h"

int kalman_filt( MEAS_T *meas, KALMAN_T *kf, PARAMS_T *params, OPTIONS_T *OPTIONS, GROUND_STATION_T *GROUND_STATION, CONSTELLATION_T  *CONSTELLATION, int iProc, int iDebugLevel, int nProcs ) {

      double et_time_step_right_after_first_meas;
  // read measurements

      read_meas(meas, meas->fp_meas, kf->init_complete, params, OPTIONS);
  //  etprint(meas->et, "newmeas");
       if ((meas->et > OPTIONS->et_final_epoch) && (kf->sc.et >= OPTIONS->et_final_epoch) ){ // the second condition (kf->sc.et >= OPTIONS->et_final_epoch) is because if not verified then the sc still needs to be propagated to the next time step of propagation (et_next_time_step)
    return 0;
  }  

  if (feof(meas->fp_meas) == 0){
    // Initialize the first time KF is called
    if ( kf->init_complete == 0 ){
      // first propagate from the propagation intial epoch (first line of section #TIME in main input file) until the first measurement. It's only a state propagation (like if no KF) since the KF starts only at the first measurement. This step is useless for the rest of the process (by "rest of the process" we mean starting at the first measurement) but is just done here to be consistent with the initial epoch chosen in the main input file
      kf->fp_kalman = fopen(kf->sc.filenamekalman, "w+");
      propagate_spaecraft_from_initial_epoch_until_first_meas( meas, kf, params, OPTIONS, GROUND_STATION,CONSTELLATION, iProc, iDebugLevel, nProcs );

      kalman_init( meas, kf, CONSTELLATION,  params, OPTIONS, iProc, iDebugLevel );
      kf->fp_kalman_meas = fopen(kf->sc.filenamekalman_meas, "w+");
      kf->is_obs = 1;
      kalman_write_out( kf, kf->fp_kalman, meas) ;


      return 0;
    }
   
    

    // Propagate State & Cov to Measurement time
    kalman_prop( meas, kf, params, OPTIONS, GROUND_STATION, CONSTELLATION, iProc, iDebugLevel, nProcs );

    if ( kf->is_obs != 0 ){ // is_obs can be equal to 0 for the last time step: if the final epoch doesn't occur at the same time as a measurement then we still want to propagate (in kalman_prop) until the final epoch (that's why there is the condition (kf->sc.et_next_time_step == OPTIONS->et_final_epoch) in the while loop of kalm_prop). So in that case, and only that case, we get out of kalm_prop not at the time of a measurement (but at the final epoch) so we should not use the KF equations for this time (since there is no measurment for that time). So in this case: is_obs = 0
    // Process GPS
    kalman_process_GPS( meas, kf);

    // Write Output

    kalman_write_out( kf, kf->fp_kalman, meas) ;
    }

  }

  return 0;
}

int kalman_joseph_update(  KALMAN_T *kf, double H[10], int meas_i) {

  // Declarations
  double HP[10] = {0.0};
  double PHt[10]= {0.0};
  double HPHt;
  double HPHtR;
  double K[10] = {0.0};
  double I_KH[10][10]= {{0.0}};
  double Pnew[10][10]= {{0.0}};
  double Pnew2[10][10]= {{0.0}};
  double KRK[10][10] = {{0.0}};
  int ii;
  int jj;

  for (ii = 0; ii < 10; ii++ ){
    HP[ii]  = 0.0;
    PHt[ii] = 0.0;
    for (jj = 0; jj < 10; jj++) {
      HP[ii]  += H[jj] * kf->P[jj][ii];
      PHt[ii] += kf->P[ii][jj] * H[jj];
    }
  }

  HPHt = 0.0;
  for (ii = 0; ii < 10; ii++ ) {
    HPHt += HP[ii] * H[ii];
  }

  HPHtR = HPHt + kf->R[meas_i];

  kf->sy[meas_i] = sqrt(kf->y[meas_i]*kf->y[meas_i] / HPHtR);
  //     HPHtR = HPHtR * (1.0 + kf->uw ); // !!!!!! cbv: not sure Joel's calculation is ok over here . More details in http://sites.utexas.edu/renato/files/2017/04/Underweight_v4.pdf. cbv chooses to comment this line for now (Joel had set uw to 0 so it doesn't really cahnge anything if commented or not anyway)

  for (ii = 0; ii < 10; ii++ ){
    K[ii]       = PHt[ii] / HPHtR;
    kf->dX[ii]  = kf->update[ii] * K[ii] * kf->y[meas_i];
    //    printf("K[%d] = %e\n",ii, K[ii]);
  }

  // K*H is outer product
  for (ii = 0; ii < 10; ii ++ ){
    for (jj = 0; jj < 10; jj++ ) {
      I_KH[ii][jj] = -K[ii]*H[jj];
      if (ii==jj) {
	I_KH[ii][jj] += 1.0;
      }
    }
  }

  for (ii = 0; ii < 10; ii ++ ){
    for (jj = 0; jj < 10; jj++ ) {
      KRK[ii][jj] = K[ii]*K[jj] * kf->R[meas_i] ; //!!!! cbv: is that outer product too?
    }
  }


  // Lazy filter update
  m_x_m10_for_kalman( Pnew, I_KH, kf->P);
  m_x_mt10_for_kalman( Pnew2, Pnew, I_KH);
  m_add10_for_kalman( kf->P , Pnew2, KRK);

  /* // !!!!!! cbv: why did Joel do that? This is basically using the typical updarta of coavria,ce, wehreas right above it's using Joseph. */
  /* for (ii = 0 ; ii < 9; ii++ ){ */
  /* 	 for (jj = 0; jj < 9; jj++ ) { */
  /* 	     kf->P[ii][jj] = Pnew[ii][jj]; */
  /* 	 } */
  /* } */
  /* // !!!!!! end of cbv: why did Joel do that? */


  return 0;
}


int kalman_process_GPS( MEAS_T *meas, KALMAN_T *kf ) {

  // Declarations
  int ii;
  double H[10] = {0.0};

  // Cycle through the measurments
  for (ii = 0; ii < 6; ii++ ){

    // Compute the residuals
    kalman_compute_resids( meas, kf, 0, ii );

    // Compute H
    kalman_computeH_gps( H, kf, ii );

    // Update State & Cov
    kalman_joseph_update( kf, H , ii);
    
    // Update State
    kalman_update_state( kf );

    // Post fit residuals
    kalman_compute_resids( meas, kf, 1, ii );

  }

  return 0;

}

int kalman_update_state(KALMAN_T *kf ){

  int ii;

  for (ii = 0; ii < 3; ii++ ){
    kf->sc.r_i2cg_INRTL[ii]   += kf->dX[ii];
    kf->sc.v_i2cg_INRTL[ii]   += kf->dX[ii+3];
    kf->sc.a_i2cg_kalman[ii]  += kf->dX[ii+6];
    if (ii == 0){
      kf->tau += kf->dX[ii+9];
    }
  }

  return 0;

}


int kalman_computeH_gps( double H[10], KALMAN_T *kf, int ii ){

  // Zero everything out
  int jj;
  for (jj = 0; jj < 10; jj++ ){
    H[jj] = 0.0;
  }
  H[ii] = 1.0;

  return 0;
}

int kalman_compute_resids( MEAS_T *meas, KALMAN_T *kf, int flag, int index ) {

  // GPS position resids
  if (flag == 0 ){ // residual before estimation = measurement - propagated state = measurement  - reference trajectory. For an EXTENDED Kalman filter, the propagated state is the propagation of the state estimated at the previous measurement (since in EKF, the reference trajectory (about which the linearization is done) is reset at each measurement to be equal to the estimated state)
    if (index < 3 ) {
      kf->y[index] = meas->r_i2cg_INRTL[index] - kf->sc.r_i2cg_INRTL[index];
    } else {
      kf->y[index] = meas->v_i2cg_INRTL[index-3] - kf->sc.v_i2cg_INRTL[index-3];
    }
  } else {  // post fit residual = measurement - estimated state. For an EXTENDED Kalman filter, the reference trajectory will be reset to be equal to this new estimated state.
    if (index < 3 ) {
      kf->y_pf[index]  = meas->r_i2cg_INRTL[index] - kf->sc.r_i2cg_INRTL[index];
    } else {
      kf->y_pf[index]  = meas->v_i2cg_INRTL[index-3] - kf->sc.v_i2cg_INRTL[index-3];
    }
  }

  return 0;

}


int kalman_init( MEAS_T *meas, KALMAN_T *kf, CONSTELLATION_T  *CONSTELLATION,       PARAMS_T *PARAMS, OPTIONS_T *OPTIONS, int iProc, int iDebugLevel) {
  char *line = NULL;
  char text[256];
  size_t len = 0;


  if ((strcmp(CONSTELLATION->spacecraft[0][0].INTEGRATOR.attitude.attitude_profile, "ensemble_angular_velocity") == 0) || (strcmp(CONSTELLATION->spacecraft[0][0].INTEGRATOR.attitude.attitude_profile, "ensemble_initial_attitude") == 0)){ // !!!!!!! all sc
    print_error(iProc, "The Kalman filter can't be run with any kind of ensembles on the attitude");
  }

  if (CONSTELLATION->spacecraft[0][0].INTEGRATOR.file_is_quaternion != 0){
    printf("\n***! There used to be an error message saying that the Kalman filter can't be run with quaternions. cbv removed this error message on 02-15-18 but maybe cbv shoudn't have done so. !***\n\n\n");
    //    print_error(iProc, "The Kalman filter can't be run with quaternions");
  }
	  
 if (OPTIONS->new_cd == 1){
   print_error(iProc, "the kalman filter can't be run with the accommodation coefficient algorithm (for the computation of cd). cd needs to be input in the geometry file, not the accommodation coefficient");
 }
  
  //

  int ii;
  int jj;
  double T[3][3]      = {{0.0}};
  double T9[9][9]     = {{0.0}};
  double Ptemp9[9][9] = {{0.0}};
  double Plvlh[9][9]  = {{0.0}};
  double P9[9][9]  = {{0.0}};

  // Zero everything out
  //     memset ( kf, 0, sizeof(&kf) ); // !!!!! cbv commented

  // Copy over the GPS
  kf->sc.et = meas->et;
     
  v_copy( kf->sc.r_i2cg_INRTL, meas->r_i2cg_INRTL);
  v_copy( kf->sc.v_i2cg_INRTL, meas->v_i2cg_INRTL);
  // !!!!!!!! TO REMOVE
  // true (6888.1370000000; 0.0000000000; 0.0000000000) (-0.0000000000; 6.2313553648; 4.3632419997)
  // noise (6888.1303268012; 0.4617380056; 0.2147756887) (-0.0007906759; 6.2294657357; 4.3617893191)
  /* 	       kf->sc.r_i2cg_INRTL[0] = 6888.1303268012; kf->sc.r_i2cg_INRTL[1] = 0.4617380056; kf->sc.r_i2cg_INRTL[2] = 0.2147756887; */
  /* kf->sc.v_i2cg_INRTL[0] = -0.0007906759; kf->sc.v_i2cg_INRTL[1] = 6.2294657357; kf->sc.v_i2cg_INRTL[2] = 4.3617893191; */

  /* kf->sc.r_i2cg_INRTL[0] = 6888.1370000000; kf->sc.r_i2cg_INRTL[1] = 0.0; kf->sc.r_i2cg_INRTL[2] = 0.0; */
  /* kf->sc.v_i2cg_INRTL[0] = -0.0; kf->sc.v_i2cg_INRTL[1] = 6.2313553648; kf->sc.v_i2cg_INRTL[2] = 4.3632419997; */
  // !!!!!!!! end of TO REMOVE
     
  // Load Params // !!!!!!!! cbv commented
  /* kf->sc.INTEGRATOR.Ap        = 6; */
  /* kf->sc.INTEGRATOR.f107      = 116;    // Daily average of F10.9 flux */
  /* kf->sc.INTEGRATOR.f107A     = 116;    // 81 day average of F10.9 flux */
  /* kf->sc.INTEGRATOR.dt        = 1;      // Integrator timestep !!!!!!! was 1 */
  /* kf->sc.INTEGRATOR.surface[0].Cd        = 2.2;      // Coefficient of drag // !!!! do it for all surfaces  */
  /* kf->sc.INTEGRATOR.mass      = 25.0;     // Mass of spacecraft */
  /* kf->sc.INTEGRATOR.surface[0].area      = 0*0.06; // !!!! do it for all surfaces */
  /* kf->sc.INTEGRATOR.degree    = 4;        // Gravity degree */
  /* kf->sc.INTEGRATOR.order     = 4;        // Gravity order */

  //    strcpy(kf->sc.filenamekalman, CONSTELLATION->spacecraft[0][0].filenamekalman); // !!!!!! all spacecraft

  //load_params( params );
  // Initialize P

  for (ii = 0; ii < 6; ii++){
    getline(&line, &len, kf->fp_kalman_init);
    RemoveSpaces(line);  strtok(line, "\n");  strtok(line, "\r"); 

    if (ii == 0){
      sscanf(line, "((%lf;%lf;%lf;%lf;%lf;%lf);", &Plvlh[ii][0], &Plvlh[ii][1], &Plvlh[ii][2], &Plvlh[ii][3], &Plvlh[ii][4], &Plvlh[ii][5]);
    }
    else if (ii == 5){
      sscanf(line, "(%lf;%lf;%lf;%lf;%lf;%lf))", &Plvlh[ii][0], &Plvlh[ii][1], &Plvlh[ii][2], &Plvlh[ii][3], &Plvlh[ii][4], &Plvlh[ii][5]);
    }
    else{
      sscanf(line, "(%lf;%lf;%lf;%lf;%lf;%lf);", &Plvlh[ii][0], &Plvlh[ii][1], &Plvlh[ii][2], &Plvlh[ii][3], &Plvlh[ii][4], &Plvlh[ii][5]);
    }
  }

  /*      for (ii = 0; ii < 6; ii++){ */
  /*        for(jj = 0; jj < 6; jj++ ) { */
	 
  /* 	 printf("%e ", Plvlh[ii][jj] ); */
  /*        } */
  /*        printf("\n"); */
  /*      } */



  Plvlh[0][5] = -0.95*sqrt(Plvlh[0][0]*Plvlh[5][5]);
  Plvlh[5][0] = Plvlh[0][5];

  compute_T_inrtl_2_lvlh(T, kf->sc.r_i2cg_INRTL, kf->sc.v_i2cg_INRTL);

  for (ii = 0; ii < 3; ii++){
    for(jj = 0; jj < 3; jj++ ) {

      T9[ii][jj]     = T[ii][jj];
      T9[ii+3][ii+3] = T[ii][jj]; // !!!!!!!!! cbv: mistake? shouldn't it be: T9[ii+3][jj+3] = T[ii][jj]; ?

    }
  }

  T9[6][6] = 1.0;

  m_x_mt9_for_kalman( Ptemp9, Plvlh, T9);
  m_x_m9_for_kalman( P9, T9, Ptemp9);
  for (ii = 0; ii < 9; ii++){
    for(jj = 0; jj < 9; jj++ ) {
      kf->P[ii][jj] = P9[ii][jj];
    }
  }
  //     m_x_m9_for_kalman( kf->P, T9, Ptemp9);

  kf->P[6][6] = 1.0e-24;
  kf->P[7][7] = 1.0e-24;
  kf->P[8][8] = 1.0e-24;

  kf->P[9][9] = 10;

  /*      // !!!!!!!!!!!!!!! to remove */
  /*      kf->P[6][6] = 1.0e-18; */
  /*      kf->P[7][7] = 1.0e-18; */
  /*      kf->P[8][8] = 1.0e-18; */
  /*      // !!!!!!!!!!!!!!! end to remove */
     
  /*      // !!!!!!!!!!!!!! to remove */
  /*      for (ii = 0; ii < 9; ii++){ */
  /*        for(jj = 0; jj < 9; jj++ ) { */
	 
  /*      	 kf->P[ii][jj] = 0.0; */
  /*        } */
  /*      } */
  /*      // !!!!!!!!!!!!!! end of to remove */


  /*      for (ii = 0; ii < 9; ii++){ */
  /*        for(jj = 0; jj < 9; jj++ ) { */
	 
  /* 	 printf("%e ", kf->P[ii][jj] ); */
  /*        } */
  /*        printf("\n"); */
  /*      } */



  // R
  /* kf->R[0] = 0.001 * 0.001; // 0.002 * 0.002 / (10000.0); */
  /* kf->R[1] = 0.001 * 0.001; // 0.002 * 0.002 / (10000.0); */
  /* kf->R[2] = 0.001 * 0.001; // 0.002 * 0.002 / (10000.0); */
  /* kf->R[3] = 0.00001 * 0.00001; // 0.002 * 0.002 / (10000.0); */
  /* kf->R[4] = 0.00001 * 0.00001; // 0.002 * 0.002 / (10000.0); */
  /* kf->R[5] = 0.00001 * 0.00001; // 0.002 * 0.002 / (10000.0); */


  getline(&line, &len, kf->fp_kalman_init);
  getline(&line, &len, kf->fp_kalman_init);
  RemoveSpaces(line);  strtok(line, "\n");  strtok(line, "\r"); 

  sscanf(line, "(%lf;%lf;%lf;%lf;%lf;%lf)", &kf->R[0], &kf->R[1], &kf->R[2], &kf->R[3], &kf->R[4], &kf->R[5]);


  kf->uw = 0;

/*   // calculate the acceleration */
/*   double v_dummy[3]; */
/*   double test_density; */
/*   double a_inrtl[3]; */
/*   double et_initial_epoch; */
/*   str2et_c(OPTIONS->initial_epoch, &et_initial_epoch); */
/*   compute_dxdt( v_dummy,a_inrtl, &kf->sc.et,kf->sc.r_i2cg_INRTL,kf->sc.v_i2cg_INRTL, PARAMS, &kf->sc.INTEGRATOR, et_initial_epoch, OPTIONS->et_oldest_tle_epoch, &test_density, kf->sc.INTEGRATOR.index_in_attitude_interpolated, kf->sc.INTEGRATOR.index_in_driver_interpolated, CONSTELLATION, OPTIONS, iProc, iDebugLevel, &kf->sc); */

  //       m_x_v(kf->sc.a_i2cg_kalman, T, kf->sc.a_i2cg_INRTL_drag);
  kf->sc.a_i2cg_kalman[0] = 0;
  kf->sc.a_i2cg_kalman[1] = 0;
  kf->sc.a_i2cg_kalman[2] = 0;



  // Tau and sigma !!!!!!! joel ussed to call his sigma my sigma**2
  getline(&line, &len, kf->fp_kalman_init);
  getline(&line, &len, kf->fp_kalman_init);
  sscanf(line, "%lf", &kf->tau);
  getline(&line, &len, kf->fp_kalman_init);
  sscanf(line, "%lf", &kf->sigma);
  getline(&line, &len, kf->fp_kalman_init);
  sscanf(line, "%lf", &kf->sigma_tau);

  /* kf->tau = 10.0*60.0; */

  /* // Sigma */
  /* kf->sigma = 1e-19; */
     
  // Update array
  kf->update[0] = 1;
  kf->update[1] = 1;
  kf->update[2] = 1;
  kf->update[3] = 1;
  kf->update[4] = 1;
  kf->update[5] = 1;
  kf->update[6] = 1;
  kf->update[7] = 1;
  kf->update[8] = 1;
  kf->update[9] = 1;
     

  // Finished completion
  kf->init_complete = 1;

  // Q is set to 1 to get in the function computeQ
  //kf->read_q = 1; 

  return 0;
}


int kalman_cov_prop( KALMAN_T *kf, PARAMS_T *PARAMS, OPTIONS_T *OPTIONS,  CONSTELLATION_T *CONSTELLATION, int iProc, int iDebugLevel ) {

  // Declarations
  double PhiP[N_STATES][N_STATES] = {{0.0}};
  double Pnew[N_STATES][N_STATES] = {{0.0}};

  // Compute STM
  kalman_computeSTM( kf->Phi, &kf->sc, PARAMS, kf->sc.INTEGRATOR.dt, kf->tau, OPTIONS,  CONSTELLATION,iProc, iDebugLevel );

  // Compute Q
  kalman_computeQ( kf->Q, kf->Phi, kf->sc.INTEGRATOR.dt, kf->tau, kf->sigma, kf->fp_kalman_init, kf->sigma_tau);
     
  // P = Phi * P * Phi^t + Q
  m_x_m10_for_kalman( PhiP, kf->Phi, kf->P);
  m_x_mt10_for_kalman( Pnew, PhiP, kf->Phi);
  m_add10_for_kalman( kf->P, Pnew, kf->Q ); 

  return 0;
}

int kalman_computeA( double A[6][6], SPACECRAFT_T *sc, PARAMS_T *PARAMS , double dt, double tau, OPTIONS_T *OPTIONS,  CONSTELLATION_T *CONSTELLATION, int iProc, int iDebugLevel) {

  // Declarations
  double r2;
  double r3;
  double r4;
  double r5;
  double rmag;
  int ii;
  int jj;
  double T[3][3]    = {{0.0}};
  double Tinv[3][3] = {{0.0}};
  double a1[3]  = {0.0};
  double a2[3]  = {0.0};
  double a[3]   = {0.0};
  double r[3]   = {0.0};
  double v[3]   = {0.0};
  double scale_h = 200.0; // in km, from table A-1 of cruisckshank98
  double r0 = 7298.1450; // in km, from table A-1 of cruisckshank98
  double rho_0 = 0.0004; // kg/km^3 , from table A-1 of cruisckshank98
  double rho_exp;
  double omega_e = 0.00007292158553; // rad/s
  double va[3]; // speed of sc relative to rotating atmo
  double vamag;
  double x = sc->r_i2cg_INRTL[0];
  double y = sc->r_i2cg_INRTL[1];
  double z = sc->r_i2cg_INRTL[2];
  double xdot = sc->v_i2cg_INRTL[0];
  double ydot = sc->v_i2cg_INRTL[1];
  double zdot = sc->v_i2cg_INRTL[2];
  double re2 = pow(PARAMS->EARTH.radius, 2);
  double j2 = PARAMS->EARTH.GRAVITY.j2;
  double mu = PARAMS->EARTH.GRAVITY.mu;
  double m = sc->INTEGRATOR.mass;
  double x2 = pow(x, 2);
  double y2 = pow(y, 2);
  double z2 = pow(z, 2);

  // Intermediate
  v_mag( &rmag, sc->r_i2cg_INRTL );
  r2 = pow( rmag, 2 );
  r3 = pow( rmag, 3 );
  r4 = pow( rmag, 4 );
  r5 = pow( rmag, 5 );
  compute_T_inrtl_2_lvlh(T, sc->r_i2cg_INRTL, sc->v_i2cg_INRTL);
  m_trans(Tinv, T);

  for (ii = 0; ii < 6 ; ii++ ){
    for (jj = 0; jj < 6; jj++ ){
      A[ii][jj] = 0.0;
    }
  }

  // eq 230 of cruisckshank98
  A[0][3] = 1;
  A[1][4] = 1;
  A[2][5] = 1;

  /*   A[jj+3][ii] = (3.0 * PARAMS->EARTH.GRAVITY.mu * sc->r_i2cg_INRTL[ii] * sc->r_i2cg_INRTL[jj] ) /r5; */
  /*   if (ii == jj ) { */
  /*     A[jj+3][ii] -= PARAMS->EARTH.GRAVITY.mu / r3; */

  va[0] = sc->v_i2cg_INRTL[0] + omega_e * sc->r_i2cg_INRTL[1];
  va[1] = sc->v_i2cg_INRTL[1] - omega_e * sc->r_i2cg_INRTL[0];
  va[2] = sc->v_i2cg_INRTL[2];
  v_mag( &vamag, va );
  rho_exp = rho_0 * exp( - ( rmag - r0 ) / scale_h);

/*   // calculate the sum of cd*A*cos(angle(va, normal to surface)) over all surfaces */
/*   sc->INTEGRATOR.sum_cd_a_cos = 0; */
/*   int sss; */
/*   double v_angle[3]; */
/*   int order_rotation[3]; */
/*   v_angle[0] = sc->INTEGRATOR.attitude.pitch[sc->INTEGRATOR.index_in_attitude_interpolated]; */
/*   v_angle[1] = sc->INTEGRATOR.attitude.roll[sc->INTEGRATOR.index_in_attitude_interpolated]; */
/*   v_angle[2] = sc->INTEGRATOR.attitude.yaw[sc->INTEGRATOR.index_in_attitude_interpolated]; */
/*   order_rotation[0] = sc->INTEGRATOR.attitude.order_pitch[sc->INTEGRATOR.index_in_attitude_interpolated]; */
/*   order_rotation[1] = sc->INTEGRATOR.attitude.order_roll[sc->INTEGRATOR.index_in_attitude_interpolated]; */
/*   order_rotation[2] = sc->INTEGRATOR.attitude.order_yaw[sc->INTEGRATOR.index_in_attitude_interpolated]; */
/*   double T_sc_to_lvlh[3][3]; */
/*   compute_T_sc_to_lvlh( T_sc_to_lvlh, v_angle, order_rotation, sc->INTEGRATOR.attitude.attitude_profile, &sc->et,  sc->r_i2cg_INRTL, sc->v_i2cg_INRTL); // body to lvlh */
/*   double normal_in_lvlh[3]; */
/*   double normal_in_inertial[3]; */
/*   double T_lvlh_2_inrtl[3][3]; */
/*   m_trans(T_lvlh_2_inrtl, T); */
/*   double va_dot_normal; */
/*   double cos_va_normal; */
/*   double normal_in_inertial_norm; */
/*   for (sss = 0; sss < sc->INTEGRATOR.nb_surfaces; sss++){ */
/* 	/\* SC to LVLH *\/ */
/* 	  m_x_v(normal_in_lvlh, T_sc_to_lvlh, sc->INTEGRATOR.surface[sss].normal ); */
/* 	  /\* LVLH to inertial *\/ */
/* 	  m_x_v(normal_in_inertial, T_lvlh_2_inrtl, normal_in_lvlh); */
/* 	  // cos(angle(va, normal to surface)) */
/* 	  v_dot(&va_dot_normal, va, normal_in_inertial); */
/* 	  v_mag(&normal_in_inertial_norm, normal_in_inertial);// just in case it's not normalized  */
/* 	  cos_va_normal = va_dot_normal / ( vamag * normal_in_inertial_norm ); */
/* 	  if (cos_va_normal > 5e-3){ // 1e-3 for numerical reasons (should be 0 theoritically) */
/* 	    sc->INTEGRATOR.sum_cd_a_cos = sc->INTEGRATOR.sum_cd_a_cos + sc->INTEGRATOR.surface[sss].Cd * sc->INTEGRATOR.surface[sss].area * cos_va_normal; */
/* 	  } */
/*   } */

  // dax / dx
  A[3][0] = -mu / r3 *
    ( 1 - 3 * x2 / r2 + 15. / 2 * j2 * re2 * z2 / r4 * ( 7 * x2 / r2 - 1 ) + 3./2 * j2 * re2 / r2 * ( -5 * x2 / r2 + 1 ) )
    + 1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp * ( xdot + omega_e * y) * ( omega_e * ( ydot - omega_e * x ) / vamag + vamag * x / ( scale_h * rmag ) );

  // dax / dy
  A[3][1] = 3 * mu * x * y / r5 *
    ( 1 - 35./2 * j2 * re2 * z2 / r4 + 5./2 * j2 * re2/r2 )
    -1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp * ( vamag * omega_e + ( xdot + omega_e * y ) * ( ( xdot + omega_e * y ) * omega_e / vamag - vamag * y / ( scale_h * rmag) ) );

  // dax / dz
  A[3][2] = 3 * mu * x * z / r5 *
    ( 1 - 35./2 * j2 * re2 * z2 / r4 + 15./2 * j2 * re2/r2 )
    + 1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp * vamag * z * ( xdot + omega_e * y ) / ( scale_h * rmag);

  // dax / dvx
  A[3][3] =  -1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp * ( vamag + ( xdot + omega_e * y ) * ( xdot + omega_e * y ) / vamag );

  // dax / dvy
  A[3][4] = -1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp / vamag * ( xdot + omega_e * y ) * ( ydot - omega_e * x );

  // dax / dvz
  A[3][5] = -1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp * zdot / vamag * ( xdot + omega_e * y ) ;

  // day / dx
  A[4][0] = 3 * mu * x * y / r5 *
    ( 1 - 35./2 * j2 * re2 * z2 / r4 + 5./2 * j2 * re2/r2 )
    -1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp * ( -vamag * omega_e - ( ydot - omega_e * x ) * ( ( ydot - omega_e * x ) * omega_e / vamag + vamag * x / ( scale_h * rmag) ) );
  ;

  // day / dy
  A[4][1] =  -mu / r3 *
    ( 1 - 3 * y2 / r2 + 15. / 2 * j2 * re2 * z2 / r4 * ( 7 * y2 / r2 - 1 ) + 3./2 * j2 * re2 / r2 * ( -5 * y2 / r2 + 1 ) )
    - 1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp * ( ydot - omega_e * x ) * ( omega_e * ( xdot + omega_e * y ) / vamag - vamag * y / ( scale_h * rmag ) );

  // day / dz
  A[4][2] =  3 * mu * y * z / r5 *
    ( 1 - 35./2 * j2 * re2 * z2 / r4 + 15./2 * j2 * re2/r2 )
    + 1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp * vamag * z * ( ydot - omega_e * x ) / ( scale_h * rmag);

  // day / dvx
  A[4][3] =  -1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp / vamag * ( xdot + omega_e * y ) * ( ydot - omega_e * x );

  // day / dvy
  A[4][4] = -1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp * ( vamag + ( ydot - omega_e * x ) * ( ydot - omega_e * x ) / vamag );

  // day / dvz
  A[4][5] = -1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp * zdot / vamag * ( ydot - omega_e * x ) ;;

  // daz / dx
  A[5][0] = 3 * mu * x * z / r5 *
    ( 1 + 15./2 * j2 * re2/r2 - 35./2 * j2 * re2 * z2 / r4  )
    + 1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp * zdot * ( omega_e * ( ydot - omega_e * x ) / vamag + vamag * x / ( scale_h * rmag) ) ;

  // daz / dy
  A[5][1] = 3 * mu * y * z / r5 *
    ( 1 + 15./2 * j2 * re2/r2 - 35./2 * j2 * re2 * z2 / r4  )
    - 1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp * zdot * ( omega_e * ( xdot + omega_e * y ) / vamag - vamag * y / ( scale_h * rmag ) ) ;

  // daz / dz
  A[5][2] =  -mu / r3 *
    ( 1 - 3 * z2 / r2 + 3. * j2 * re2 / r2 * ( 3./2 - 15. * z2 / r2 + 35./2 * z2 / r4 ) )
    + 1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp * vamag * z * zdot / ( scale_h * rmag ) ;

  // daz / dvx
  A[5][3] =  -1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp * zdot * (  xdot + omega_e * y  ) / vamag;

  // daz / dvy
  A[5][4] = -1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp * zdot * ( ydot - omega_e * x ) / vamag;

  // daz / dvz
  A[5][5] =  -1./2 * sc->INTEGRATOR.sum_cd_a_cos / m * rho_exp * ( zdot * zdot / vamag + vamag ) ;
     
     
  /* printf(" --------- A --------- \n"); */
  /* for (jj = 0; jj < 6; jj++ ){ */
  /*   printf("   %e %e %e %e %e %e %e %e %e\n", A[jj][0], A[jj][1],A[jj][2],A[jj][3],A[jj][4],A[jj][5], A[jj][6],A[jj][7],A[jj][8] ); */
  /* } */

  return 0;
}

int kalman_computeSTM( double Phi[10][10], SPACECRAFT_T *sc, PARAMS_T *PARAMS, double dt, double tau,OPTIONS_T *OPTIONS,  CONSTELLATION_T *CONSTELLATION, int iProc, int iDebugLevel ) {


  int ii, jj;
  for (ii = 0; ii < 10 ; ii++ ){
    for (jj = 0; jj < 10; jj++ ){
      Phi[ii][jj] = 0.0;
    }
  }

  // use equation 19 of cruisckshank98: determine first phif from A of appendix a of cruisckshank98 (A is calculated in kalman_computeA using appendix A of cruisckshank98). The calculate phiw,  M, phitau, and N (using notation of eq 19 of cruisckshank98)
  double phif[6][6];
  for (ii = 0; ii < 6 ; ii++ ){
    for (jj = 0; jj < 6; jj++ ){
      phif[ii][jj] = 0.0;
    }
  }

  // Compute A and from A get phif
  double A[6][6] = {{0.0}};
  double A2[6][6] = {{0.0}};
  kalman_computeA( A , sc, PARAMS, dt, tau, OPTIONS,  CONSTELLATION, iProc, iDebugLevel );
  m_x_m6_for_kalman( A2, A, A );

  // phif = I + A*dt
  for (ii = 0; ii < 6; ii++) {
    for (jj = 0; jj<6; jj++) {
      phif[ii][jj] = A[ii][jj] * dt + A2[ii][jj] * dt *dt * 0.5;
      if (ii == jj){
	phif[ii][jj] += 1.0;
      }
    }
  }

  // Compute phiw
  double phiwv;
  double phiwp;
  double beta;
  beta = 1 / tau;
  phiwp = 1 / ( beta*beta ) * ( exp( -beta * dt ) - 1 ) + 1 / beta * dt;
  phiwv = 1 / beta * ( 1 - exp( -beta * dt ) );
  double phiw[6][3];
  for (ii = 0; ii < 6; ii++) {
    for (jj = 0; jj<3; jj++) {
      phiw[ii][jj] = 0;
    }
  }
  phiw[0][0] = phiwp;
  phiw[1][1] = phiwp;
  phiw[2][2] = phiwp;

  phiw[3][0] = phiwv;
  phiw[4][1] = phiwv;
  phiw[5][2] = phiwv;

  // Compute M
  double M[3][3];
  for (ii = 0; ii < 3; ii++) {
    for (jj = 0; jj<3; jj++) {
      M[ii][jj] = 0;
    }
  }
  M[0][0] = exp( - beta * dt );
  M[1][1] = exp( - beta * dt );
  M[2][2] = exp( - beta * dt );

  // Compute phitau
  double wx0 = sc->a_i2cg_kalman[0];
  double wy0 = sc->a_i2cg_kalman[1];
  double wz0 = sc->a_i2cg_kalman[2];
  double phitau[6];
  phitau[0] = 2 * wx0 / beta * ( exp( -beta * dt ) - 1 ) + wx0 * dt * ( exp( -beta * dt ) + 1 );
  phitau[1] = 2 * wy0 / beta * ( exp( -beta * dt ) - 1 ) + wy0 * dt * ( exp( -beta * dt ) + 1 );
  phitau[2] = 2 * wz0 / beta * ( exp( -beta * dt ) - 1 ) + wz0 * dt * ( exp( -beta * dt ) + 1 );
  phitau[3] = wx0 * ( 1 - exp( -beta * dt ) ) - wx0 * beta * dt * exp( -beta * dt );
  phitau[4] = wy0 * ( 1 - exp( -beta * dt ) ) - wy0 * beta * dt * exp( -beta * dt );
  phitau[5] = wz0 * ( 1 - exp( -beta * dt ) ) - wz0 * beta * dt * exp( -beta * dt );
     
  // Compute N
  double N[3];
  N[0]  = wx0 * beta * beta * dt * exp( -beta * dt );
  N[1]  = wy0 * beta * beta * dt * exp( -beta * dt );
  N[2]  = wz0 * beta * beta * dt * exp( -beta * dt );
     
  // Combine phif, phiw, M, phitau, and N to get Phi
  for (ii = 0; ii < 6; ii++) {
    for (jj = 0; jj < 6; jj++) {
      Phi[ii][jj] = phif[ii][jj];
    }
  }
  for (ii = 0; ii < 6; ii++) {
    for (jj = 6; jj < 9; jj++) {
      Phi[ii][jj] = phiw[ii][jj-6];
    }
  }
  for (ii = 6; ii < 9; ii++) {
    for (jj = 6; jj < 9; jj++) {
      Phi[ii][jj] = M[ii-6][jj-6];
    }
  }

  for (ii = 0; ii < 6; ii++) {
    Phi[ii][9] = phitau[ii-6];
  }

  for (ii = 6; ii < 9; ii++) {
    Phi[ii][9] = N[ii-6];
  }

  Phi[9][9] = 1;
     
  return 0;
}

int kalman_computeQ( double Q[10][10], double Phi[N_STATES][N_STATES], double dt, double tau, double sigma, FILE *fp_kalman_init, double sigma_tau) {
  char *line = NULL;
  char text[256];
  size_t len = 0;

  int ii;
  int jj;

  for (ii = 0; ii < 10; ii++ ){
    for (jj = 0; jj < 10; jj++ ){

      Q[ii][jj] = 0.0;
    }
  }
     
  // eq 6.51 ad 6.52 of lee05, equiavalent to eq 69 and 70  of cruisckshank98
  double lambda_scal = tau * sigma * sigma / 2. * ( 1 - exp( -2 * dt / tau ) );
  for (ii = 0 ; ii < 3 ; ii++ ){
    Q[ii][ii]     = pow( dt, 4) / 4.0 * lambda_scal;
    Q[ii][ii+3]   = pow( dt, 3) / 2.0 * lambda_scal;
    Q[ii+3][ii]   = pow( dt, 3) / 2.0 * lambda_scal;
    Q[ii+3][ii+3] = pow( dt, 2) * lambda_scal;
    Q[ii+3][ii+6] = pow( dt, 2) / 2. * lambda_scal;
    Q[ii+6][ii] = pow( dt, 2) / 2. * lambda_scal;
    Q[ii+3][ii+6] = dt * lambda_scal;
    Q[ii+6][ii+3] = dt * lambda_scal;
    Q[ii+6][ii+6] =  lambda_scal;
  }
  Q[9][9] = sigma_tau *  sigma_tau * dt;



  /*      double Qo[N_STATES][N_STATES]    = {{0.0}}; */
  /*      double Qtemp[N_STATES][N_STATES] = {{0.0}}; */
     
  /*      for (ii = 0; ii < N_STATES; ii++ ){ */
  /* 	 for (jj = 0; jj < N_STATES; jj++ ){ */

  /* 	     Qo[ii][jj] = 0.0; */
  /* 	 } */
  /*      } */

  /*      for (ii = 0 ; ii < 3 ; ii++ ){ */

  /* 	 /\* Qo[ii][ii]       = sigma * pow( dt, 4) / 4.0; *\/ */
  /* 	 /\* Qo[ii][ii+3]     = sigma * pow( dt, 3) / 2.0; *\/ */
  /* 	 /\* Qo[ii+3][ii]     = sigma * pow( dt, 3) / 2.0; *\/ */
  /* 	 /\* Qo[ii+3][ii+3]   = sigma * pow( dt, 2); *\/ */

  /*      } */


  /*      //          printf("q = %e | %e | %e\n\n", sigma * ( 1.0 - exp( -2.0 * dt/ tau)), exp( -2.0 * dt/ tau), exp( - dt/ tau)); */
  /*      /\* Qo[6][6] = sigma*sigma * tau/2. * ( 1.0 - exp( -2.0 * dt/ tau)); *\/ */
  /*      /\* Qo[7][7] = sigma*sigma * tau/2. * ( 1.0 - exp( -2.0 * dt/ tau)); *\/ */
  /*      /\* Qo[8][8] = sigma*sigma * tau/2. * ( 1.0 - exp( -2.0 * dt/ tau)); *\/ */

  /*      Qo[6][6] = sigma * ( 1.0 - exp( -2.0 * dt/ tau)); */
  /*      Qo[7][7] = sigma * ( 1.0 - exp( -2.0 * dt/ tau)); */
  /*      Qo[8][8] = sigma * ( 1.0 - exp( -2.0 * dt/ tau)); */

  /*      /\* Qo[6][6] = sigma * sqrt( 1.0 - exp( -2.0 * dt/ tau)); *\/ */
  /*      /\* Qo[7][7] = sigma * sqrt( 1.0 - exp( -2.0 * dt/ tau)); *\/ */
  /*      /\* Qo[8][8] = sigma * sqrt( 1.0 - exp( -2.0 * dt/ tau)); *\/ */

  /*      m_x_mt6_for_kalman( Qtemp, Qo, Phi); */
  /*      m_x_m6_for_kalman( Q, Phi, Qtemp ); */



  return 0;
}


int kalman_prop(  MEAS_T *meas, KALMAN_T *kf, PARAMS_T *PARAMS, OPTIONS_T *OPTIONS, GROUND_STATION_T *GROUND_STATION,CONSTELLATION_T  *CONSTELLATION, int iProc, int iDebugLevel, int nProcs ) {
  
  // Declarations
  double dt;
int       already_propagated_until_time_step_right_after_meas = 0;
 int need_to_reinitialize_dt = 1;
  // Save off dt in case we have to update it
  dt = kf->sc.INTEGRATOR.dt;

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
  int ii;
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
  int ccc;
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


  double density;
  double starttime;
  str2et_c(OPTIONS->initial_epoch, &starttime);

  while ( ( kf->sc.et != meas->et ) || ( (kf->sc.et_next_time_step == OPTIONS->et_final_epoch) && (meas->et>OPTIONS->et_final_epoch)) ){// get in here to propagate until the next measurment or propagate one more time from the previous measurement until the final epoch


    if ( ( meas->et - kf->sc.et  ) < kf->sc.INTEGRATOR.dt ){
      kf->sc.INTEGRATOR.dt = meas->et - kf->sc.et;
    }
    //        etprint(kf->sc.et, "kf->sc.et");

    if ( ( ( kf->sc.et_next_time_step - kf->sc.et  ) < kf->sc.INTEGRATOR.dt ) && ( already_propagated_until_time_step_right_after_meas == 0 ) ) {
      kf->sc.INTEGRATOR.dt = kf->sc.et_next_time_step - kf->sc.et;
      already_propagated_until_time_step_right_after_meas = 1;
    }

    
    // Propagate covariance matrix
    kalman_cov_prop( kf, PARAMS, OPTIONS,  CONSTELLATION ,  iProc,  iDebugLevel );

    // Propagate state
    // Position and velocity 
    propagate_spacecraft( &kf->sc, PARAMS,starttime, OPTIONS->et_oldest_tle_epoch, &density, GROUND_STATION, OPTIONS, CONSTELLATION, iProc, iDebugLevel,  start_ensemble, array_sc );
    
    // Gauss Markov process acceleration
    kf->sc.a_i2cg_kalman[0] = exp( -kf->sc.INTEGRATOR.dt / kf->tau ) * kf->sc.a_i2cg_kalman[0];
    kf->sc.a_i2cg_kalman[1] = exp( -kf->sc.INTEGRATOR.dt / kf->tau ) * kf->sc.a_i2cg_kalman[1];
    kf->sc.a_i2cg_kalman[2] = exp( -kf->sc.INTEGRATOR.dt / kf->tau ) * kf->sc.a_i2cg_kalman[2];

    if (  kf->sc.et != meas->et ){ // without this if, at each measrurement, the output file would have the prediciton and the estimate at this same time (so two lines for the same time, one predicted, the other one estimated). With this if, the predicted state at the measuremnt is not printed in the kalman output file`
/*       etprint(kf->sc.et, "kf->sc.et"); */
/*       etprint(meas->et, "meas->et"); */
/*       printf("%d\n", kf->sc.et != meas->et ); */
	
      kf->is_obs = 0;
      kalman_write_out( kf, kf->fp_kalman, meas) ;
    }
    // Gauss Markov process correlation time
    //    kf->tau = kf->tau; // the law of propagation of tau is that tau is constant

    kf->sc.INTEGRATOR.dt = dt; // if got in the if condition above "if ( ( ( kf->sc.et_next_time_step - kf->sc.et  ) < kf->sc.INTEGRATOR.dt ) " then we need to reset kf->sc.INTEGRATOR.dt to its intial value (dt). 


    if (kf->sc.et >= OPTIONS->et_final_epoch){
      break;
    }

  }

    kf->is_obs = 1;

  if (kf->sc.et != meas->et){ // this could happen for the last time step: if the final epoch doesn't occur at the same time as a measurement then we still want to propagate until the final epoch (that's why there is the condition (kf->sc.et_next_time_step == OPTIONS->et_final_epoch) in the previous while loop). So in that case, and only that case, we get out of kalm_prop not at the time of a measurement (but at the final epoch) so we should not use the KF equations for this time (since there is no measurment for that time)
    kf->is_obs = 0;
  }




  kf->sc.INTEGRATOR.dt = dt;

  return 0;

}


int kalman_write_out( KALMAN_T *kf, FILE *fp, MEAS_T *meas) {

  // Declarations
  double r[3] = {0.0};
  double v[3] = {0.0};
  double et;
  double ad;
  double rho;
  double drag_sigma = 0.0;
  double dd_dx[N_STATES] = {0.0};
  double term;
  double Area;
  double accel = 0.0;
  double v2;
  double prod[N_STATES] = {0.0};
  int ii;
  int jj;
  double T[3][3], T_inv[3][3];
  double v_estimate_lvlh[3];
  double amag_estimate, vmag_estimate;
  char text2[256];
  double rho_estimate;
  
  et = kf->sc.et;
  v_copy( r, kf->sc.r_i2cg_INRTL );
  v_copy( v, kf->sc.v_i2cg_INRTL );
  ad  = kf->sc.a_i2cg_kalman[0];

  // cbv did that

    double a_drag_inrtl_estimate[3];
  v_add(a_drag_inrtl_estimate, kf->sc.a_i2cg_INRTL_drag,kf->sc.a_i2cg_kalman);
  double a_drag_inrtl_estimate_mag;
  v_mag(&a_drag_inrtl_estimate_mag, a_drag_inrtl_estimate);
  double vmag;
  v_mag(&vmag, kf->sc.v_i2cg_INRTL);


  double va[3];
  double vamag;
  double omega_e = 0.00007292158553; // rad/s
  va[0] = kf->sc.v_i2cg_INRTL[0] + omega_e * kf->sc.r_i2cg_INRTL[1];
  va[1] = kf->sc.v_i2cg_INRTL[1] - omega_e * kf->sc.r_i2cg_INRTL[0];
  va[2] = kf->sc.v_i2cg_INRTL[2];
  v_mag( &vamag, va );

  
  rho_estimate = a_drag_inrtl_estimate_mag /  ( 1./2 * kf->sc.INTEGRATOR.sum_cd_a_cos / kf->sc.INTEGRATOR.mass * vamag*vamag ) / ( 1000*1000*1000 ) ;// in kg/m^3 
  
  //    printf("%e %e %e\n", a_drag_inrtl_estimate_mag, kf->sc.INTEGRATOR.sum_cd_a_cos, vamag);


  
  // end of cbv did that

  // cbv doesn't use that (Joel used it though (cbv made modifications to it but then stopped using it))
  /* rho = kf->sc.rho;  */
  /* v2 = kf->sc.v_i2cg_INRTL[0]*kf->sc.v_i2cg_INRTL[0] + kf->sc.v_i2cg_INRTL[1]*kf->sc.v_i2cg_INRTL[1] + kf->sc.v_i2cg_INRTL[2]*kf->sc.v_i2cg_INRTL[2]; */
  /* Area = kf->sc.INTEGRATOR.surface[0].area / ( 1000.0 * 1000 ); // km^2 to m^2 !!!!!! all surfaces  */
  /* compute_T_inrtl_2_lvlh(T, kf->sc.r_i2cg_INRTL, kf->sc.v_i2cg_INRTL); */
  /* m_x_v(v_estimate_lvlh, T, kf->sc.v_i2cg_INRTL); */
       
  /* v_mag(&amag_estimate, kf->sc.a_i2cg_kalman); // a_i2cg_kalman is supposed ot represent the unmodeled accleration */
  /* v_mag(&vmag_estimate, v_estimate_lvlh); */
  /* //     double a_drag_estimate[3];// a_drag_estimate  = a_total - a_gravity (only 2 forces: drag and gravity) */
  /* rho_estimate = fabs(kf->sc.a_i2cg_kalman[0]/(1/2. * v_estimate_lvlh[0] * vmag_estimate * kf->sc.INTEGRATOR.surface[0].Cd * kf->sc.INTEGRATOR.surface[0].area / kf->sc.INTEGRATOR.mass))/(1e9); // in kg/m^3 */
  /* //     rho_estimate = amag_estimate/(1/2. * v2 * kf->sc.INTEGRATOR.surface[0].Cd * kf->sc.INTEGRATOR.surface[0].area / kf->sc.INTEGRATOR.mass); */
  /* //printf("%e %e %e %e %e %e\n", rho_estimate, amag_estimate, v2, kf->sc.INTEGRATOR.surface[0].Cd, kf->sc.INTEGRATOR.surface[0].area, kf->sc.INTEGRATOR.mass); */
 
  /* // cbv: accel = ax (without the minus sign) */
  /* accel = 0.5 * rho * v2 * kf->sc.INTEGRATOR.surface[0].Cd * Area / kf->sc.INTEGRATOR.mass;// !!!!!! all surfaces */
  /* term = 2.0 * kf->sc.INTEGRATOR.mass / (kf->sc.INTEGRATOR.surface[0].Cd * Area );// !!!!!! all surfaces */
  /* // cbv: density = term * accel / v^2 */
  /* // cbv: dd_dx seems to be d(rho)/dx (= 0?), d(rho)/dy (= 0?), d(rho)/dz (= 0?), d(rho)/dvx = (see dd_dx[3] below), d(rho)/dvy = (see dd_dx[4] below), d(rho)/dvz = (see dd_dx[5] below), d(rho)/dax = (see dd_dx[6] below, note that ax = accel), d(rho)/day = 0, d(rho)/daz = 0 */
    
  /* dd_dx[3] = -term * accel / (v2 * v2 ) * 2 * kf->sc.v_i2cg_INRTL[0]; */
  /* dd_dx[4] = -term * accel/ (v2 * v2 ) * 2 * kf->sc.v_i2cg_INRTL[1]; */
  /* dd_dx[5] = -term * accel/ (v2 * v2 ) * 2 * kf->sc.v_i2cg_INRTL[2]; */
    
  /* dd_dx[6] = term / v2; */
    
  /* for (ii = 0; ii < N_STATES; ii++ ){ */
  /*   prod[ii] = 0.0; */
  /*   for (jj = 0; jj< N_STATES; jj++ ) { */
  /*     prod[ii] += dd_dx[jj]*kf->P[jj][ii]; */
  /*   } */
  /* } */

  /* drag_sigma = 0.0; */
  /* for (ii = 0; ii < N_STATES; ii++ ){ */
    
  /*   drag_sigma += prod[ii] * dd_dx[ii]; */
    
  /* } */
  // end of cbv doesn't use that (Joel used it though (cbv made modifications to it but then stopped using it))

  et2utc_c(et, "ISOC", 6, 255, text2);

  fprintf(fp, "%s %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f",
	  text2, r[0], r[1], r[2], v[0], v[1], v[2]);
  fprintf(fp, " %e %e %e",
	  kf->sc.a_i2cg_INRTL_drag[0], kf->sc.a_i2cg_INRTL_drag[1], kf->sc.a_i2cg_INRTL_drag[2]);
  fprintf( fp, " %e", kf->tau ); 
  fprintf(fp, " %e %e %e",
	  kf->sc.a_i2cg_kalman[0], kf->sc.a_i2cg_kalman[1], kf->sc.a_i2cg_kalman[2]);



  fprintf(fp," %e",
	  rho_estimate);

  fprintf(fp, " %e", 
	  kf->sc.INTEGRATOR.sum_cd_a_cos );

      fprintf(fp, " %6.3f %6.3f %6.3f",
	    kf->sc.INTEGRATOR.attitude.pitch_current,
	    kf->sc.INTEGRATOR.attitude.roll_current,
	      kf->sc.INTEGRATOR.attitude.yaw_current);

      //      printf("%s %e\n",text2,kf->sc.INTEGRATOR.A_ref_tot);
  if (kf->is_obs == 1){ // if this time corresponds to a time where there is a measurement
    fprintf(fp, " obs");
  }
  else{// otherwise the state is not estimated but only predicted (estimated is when we use the equations of the Kalman filter to use the measrument to get the estimate)
    fprintf(fp, " pred");
  }






  // Measurement file converted to the same unit as the state
  if (kf->is_obs == 1){ // if this time corresponds to a time where there is a measurement
  char text_meas[300];
    et2utc_c(meas->et, "ISOC", 6, 255, text_meas);
  fprintf(kf->fp_kalman_meas, "%s %e %e %e %e %e %e\n",
	  text_meas, meas->r_i2cg_INRTL[0], meas->r_i2cg_INRTL[1], meas->r_i2cg_INRTL[2], meas->v_i2cg_INRTL[0], meas->v_i2cg_INRTL[1], meas->v_i2cg_INRTL[2] );
  fprintf(fp, " %e %e %e %e %e %e",
	  kf->y_pf[0], kf->y_pf[1], kf->y_pf[2], kf->y_pf[3], kf->y_pf[4], kf->y_pf[5]);

  }

  fprintf(fp, " \n");

  // OLD
/*   fprintf(fp, " %e %e %e %e %e %e", */
/* 	  kf->P[0][0], kf->P[1][1], kf->P[2][2], kf->P[3][3], kf->P[4][4], kf->P[5][5]); */
/*   fprintf(fp, " %e %e %e %e %e %e", */
/* 	  kf->dX[0], kf->dX[1], kf->dX[2], kf->dX[3], kf->dX[4], kf->dX[5]); */
/*   fprintf(fp, " %e %e %e", */
/* 	  kf->y[0], kf->y[1], kf->y[2]); */
/*   fprintf(fp, " %e %e %e", */
/* 	  kf->sy[0], kf->sy[1], kf->sy[2]); */
/*   fprintf(fp, " %e %e %e", */
/* 	  kf->y_pf[0], kf->y_pf[1], kf->y_pf[2]); */
/*   fprintf(fp, " %e %e %e", */
/* 	  kf->sc.a_i2cg_kalman[0], kf->sc.a_i2cg_kalman[1], kf->sc.a_i2cg_kalman[2]); */
/*   fprintf(fp," %e %e %e %e %e", */
/* 	  ad, rho_estimate, kf->dX[6], kf->P[6][6], sqrt(drag_sigma));  */
/*   fprintf( fp, " %e", kf->tau );  */
/*   fprintf(fp, " %e %e %e", */
/* 	  kf->sc.a_i2cg_INRTL_drag[0], kf->sc.a_i2cg_INRTL_drag[1], kf->sc.a_i2cg_INRTL_drag[2]); */
  // end of OLD



/* char  time_temp_str[256]; */
/*   double time_temp; */
/*   strcpy(time_temp_str, "2017-05-04T00:00:40.999261"); */
/*   str2et_c(time_temp_str, &time_temp); */
/*   if (kf->sc.et >= time_temp){ */
/*     MPI_Finalize(); exit(0); */
/*   } */

  return 0;

}


int read_meas(MEAS_T *meas, FILE *fp_meas, int init_complete, PARAMS_T *PARAMS, OPTIONS_T *OPTIONS){
  SpiceDouble       xform[6][6];
  double estate[6], jstate[6];
    double fake_d;
  int found_eoh;
  char *line = NULL;
  char text[300];
  int ierr;
  size_t len = 0;
  double meas_r_ecef[3];
  double meas_v_ecef[3];
  if ( init_complete == 0 ){ // first time measurement file is open so skip header
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
  // skip all observations before initial epoch

  while (meas->et < OPTIONS->et_initial_epoch){
    getline(&line, &len, fp_meas);

      sscanf(line, "%s %lf %lf %lf %lf %lf %lf", text, &fake_d,  &fake_d,  &fake_d,  &fake_d,  &fake_d,  &fake_d);
      str2et_c(text, &meas->et);

  }

  } // end of if ( init_complete == 0 )



  else{
    getline(&line, &len, fp_meas);
  }
  if (feof(fp_meas) == 0){
/*     // here the measurment file is directly in ECI r/v */
/*     sscanf(line, "%s %lf %lf %lf %lf %lf %lf", text, &meas->r_i2cg_INRTL[0],  &meas->r_i2cg_INRTL[1],  &meas->r_i2cg_INRTL[2], &meas->v_i2cg_INRTL[0],  &meas->v_i2cg_INRTL[1],  &meas->v_i2cg_INRTL[2]); */
/*     str2et_c(text, &meas->et); */
/*     // end of here the measurment file is directly in ECI r/v  */

      // !!!!!!!!! the file I use from the netcdf data is in ECEF so need to convert to ECI
   
      sscanf(line, "%s %lf %lf %lf %lf %lf %lf", text, &meas_r_ecef[0],  &meas_r_ecef[1],  &meas_r_ecef[2], &meas_v_ecef[0],  &meas_v_ecef[1],  &meas_v_ecef[2]);
    str2et_c(text, &meas->et);
    // // Convert ECEF to ECI


    	    estate[0] = meas_r_ecef[0];estate[1] = meas_r_ecef[1];estate[2] = meas_r_ecef[2];
    	    estate[3] = meas_v_ecef[0];estate[4] = meas_v_ecef[1];estate[5] = meas_v_ecef[2];
    	    sxform_c (  PARAMS->EARTH.earth_fixed_frame,  "J2000", meas->et,    xform  );
    	    mxvg_c   (  xform,       estate,   6,  6, jstate );
    	    meas->r_i2cg_INRTL[0] = jstate[0]; meas->r_i2cg_INRTL[1] = jstate[1]; meas->r_i2cg_INRTL[2] = jstate[2];
    	    meas->v_i2cg_INRTL[0] = jstate[3]; meas->v_i2cg_INRTL[1] = jstate[4]; meas->v_i2cg_INRTL[2] = jstate[5];




/*       char text2[256]; */
/*      	et2utc_c(meas->et, "ISOC", 6, 255, text2); */
/*       printf("<%s> <%s>\n", text, text2); */

  
  }


  
  return 0;
}

int m_x_m9_for_kalman( double M[9][9], double M1[9][9], double M2[9][9] ){

  // Declarations
  int ii, jj, kk;

  for (ii = 0 ; ii < 9; ii++ ){

    for (jj = 0; jj < 9; jj++ ) {
      M[ii][jj] = 0;
      for (kk = 0; kk < 9; kk++) {
	M[ii][jj] += (M1[ii][kk] * M2[kk][jj]);
      }
    }
  }

  return 0;
}

int m_x_mt9_for_kalman( double M[9][9], double M1[9][9], double M2[9][9] ){

  // Declarations
  int ii, jj, kk;

  for (ii = 0 ; ii < 9; ii++ ){

    for (jj = 0; jj < 9; jj++ ) {
      M[ii][jj] = 0;
      for (kk = 0; kk < 9; kk++) {
	M[ii][jj] += (M1[ii][kk] * M2[jj][kk]);
      }
    }
  }

  return 0;
}

int m_add9_for_kalman( double M[9][9], double M1[9][9], double M2[9][9] ) {

  // Declarations
  int ii, jj;

  for (ii = 0 ; ii < 9; ii++ ){
    for (jj = 0; jj < 9; jj++ ) {

      M[ii][jj] = M1[ii][jj] + M2[ii][jj];

    }
  }


  return 0;
}


int m_x_m10_for_kalman( double M[10][10], double M1[10][10], double M2[10][10] ){

  // Declarations
  int ii, jj, kk;

  for (ii = 0 ; ii < 10; ii++ ){

    for (jj = 0; jj < 10; jj++ ) {
      M[ii][jj] = 0;
      for (kk = 0; kk < 10; kk++) {
	M[ii][jj] += (M1[ii][kk] * M2[kk][jj]);
      }
    }
  }

  return 0;
}

int m_x_mt10_for_kalman( double M[10][10], double M1[10][10], double M2[10][10] ){

  // Declarations
  int ii, jj, kk;

  for (ii = 0 ; ii < 10; ii++ ){

    for (jj = 0; jj < 10; jj++ ) {
      M[ii][jj] = 0;
      for (kk = 0; kk < 10; kk++) {
	M[ii][jj] += (M1[ii][kk] * M2[jj][kk]);
      }
    }
  }

  return 0;
}

int m_add10_for_kalman( double M[10][10], double M1[10][10], double M2[10][10] ) {

  // Declarations
  int ii, jj;

  for (ii = 0 ; ii < 10; ii++ ){
    for (jj = 0; jj < 10; jj++ ) {

      M[ii][jj] = M1[ii][jj] + M2[ii][jj];

    }
  }


  return 0;
}


int m_x_m6_for_kalman( double M[6][6], double M1[6][6], double M2[6][6] ){

  // Declarations
  int ii, jj, kk;

  for (ii = 0 ; ii < 6; ii++ ){

    for (jj = 0; jj < 6; jj++ ) {
      M[ii][jj] = 0;
      for (kk = 0; kk < 6; kk++) {
	M[ii][jj] += (M1[ii][kk] * M2[kk][jj]);
      }
    }
  }

  return 0;
}

int m_x_mt6_for_kalman( double M[6][6], double M1[6][6], double M2[6][6] ){

  // Declarations
  int ii, jj, kk;

  for (ii = 0 ; ii < 6; ii++ ){

    for (jj = 0; jj < 6; jj++ ) {
      M[ii][jj] = 0;
      for (kk = 0; kk < 6; kk++) {
	M[ii][jj] += (M1[ii][kk] * M2[jj][kk]);
      }
    }
  }

  return 0;
}

int m_add6_for_kalman( double M[6][6], double M1[6][6], double M2[6][6] ) {

  // Declarations
  int ii, jj;
  for (ii = 0 ; ii < 6; ii++ ){
    for (jj = 0; jj < 6; jj++ ) {
      M[ii][jj] = M1[ii][jj] + M2[ii][jj];
    }
  }
  return 0;
}


int propagate_spaecraft_from_initial_epoch_until_first_meas( MEAS_T *MEAS,  KALMAN_T *kf, PARAMS_T *PARAMS, OPTIONS_T *OPTIONS, GROUND_STATION_T *GROUND_STATION, CONSTELLATION_T  *CONSTELLATION, int iProc, int iDebugLevel, int nProcs ) {


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
  int ii;
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
  int ccc;
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
  double density;
      kalman_write_out( kf, kf->fp_kalman, MEAS) ;
      while (kf->sc.et + OPTIONS->dt < MEAS->et){
       
    propagate_spacecraft( &kf->sc, PARAMS, OPTIONS->et_initial_epoch, OPTIONS->et_oldest_tle_epoch, &density, GROUND_STATION, OPTIONS, CONSTELLATION, iProc, iDebugLevel,  start_ensemble, array_sc );
      kf->is_obs = 0;
      kalman_write_out( kf, kf->fp_kalman, MEAS) ;

      }
      // propagate one more time until first meas
      double dt_save  = kf->sc.INTEGRATOR.dt;
      if ( (( MEAS->et - kf->sc.et  ) < kf->sc.INTEGRATOR.dt ) && (fabs(MEAS->et - kf->sc.et) > 1e-6)) { // (fabs(MEAS->et - kf->sc.et) > 1e-6) is to eliminate the case where MEAS->et = kf->sc.et (taking into account numerical errors)
      kf->sc.INTEGRATOR.dt = MEAS->et - kf->sc.et;
          propagate_spacecraft( &kf->sc, PARAMS, OPTIONS->et_initial_epoch, OPTIONS->et_oldest_tle_epoch, &density, GROUND_STATION, OPTIONS, CONSTELLATION, iProc, iDebugLevel,  start_ensemble, array_sc );
      kf->is_obs = 0;
      //      kalman_write_out( kf, kf->fp_kalman, MEAS) ;

    }
	  kf->sc.INTEGRATOR.dt = dt_save;

/* 	  	  etprint(kf->sc.et, "first prop done"); */
/* 	  print_test(); */
  return 0;

}
