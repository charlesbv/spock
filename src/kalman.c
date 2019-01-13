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

#include "kalman.h"


int kalman_filt( MEAS_T *meas, KALMAN_T *kf, PARAMS_T *params, OPTIONS_T *OPTIONS, GROUND_STATION_T *GROUND_STATION, CONSTELLATION_T  *CONSTELLATION, int iProc, int iDebugLevel, int nProcs ) {

  // read measurements

  read_meas(meas, meas->fp_meas, kf->init_complete, params);
  if (feof(meas->fp_meas) == 0){
    // Initialize if necessary
    if ( kf->init_complete == 0 ){
   

      kalman_init( meas, kf, CONSTELLATION );
      printf("Iniitalizing Nav \n");
      
      char fp_name[256];

      kf->fp_kalman = fopen(kf->sc.filenamekalman, "w+");

      kalman_write_out( kf, kf->fp_kalman) ;

	       return 0;
        
    }

    
    
    

    // Propagate State & Cov to Measurement time
    kalman_prop( meas, kf, params, OPTIONS, GROUND_STATION, CONSTELLATION, iProc, iDebugLevel, nProcs );
    
    // Process GP
    kalman_process_GPS( meas, kf);

     // Write Output
     kalman_write_out( kf, kf->fp_kalman) ;

  }

     return 0;
 }

 int kalman_joseph_update(  KALMAN_T *kf, double H[9], int meas_i) {

     // Declarations
     double HP[9] = {0.0};
     double PHt[9]= {0.0};
     double HPHt;
     double HPHtR;
     double K[9] = {0.0};
     double I_KH[9][9]= {{0.0}};
     double Pnew[9][9]= {{0.0}};
     double Pnew2[9][9]= {{0.0}};
     double KRK[9][9] = {{0.0}};
     int ii;
     int jj;

     for (ii = 0; ii < 9; ii++ ){
	 HP[ii]  = 0.0;
	 PHt[ii] = 0.0;
	 for (jj = 0; jj < 9; jj++) {
	     HP[ii]  += H[jj] * kf->P[jj][ii];
	     PHt[ii] += kf->P[ii][jj] * H[jj];
	 }
     }

     HPHt = 0.0;
     for (ii = 0; ii < 9; ii++ ) {

	 HPHt += HP[ii] * H[ii];

     }

     HPHtR = HPHt + kf->R[meas_i];

     kf->sy[meas_i] = sqrt(kf->y[meas_i]*kf->y[meas_i] / HPHtR);

     //     HPHtR = HPHtR * (1.0 + kf->uw ); // !!!!!! cbv: not sure Joel's calculation is ok over here . More details in http://sites.utexas.edu/renato/files/2017/04/Underweight_v4.pdf. cbv chooses to comment this line for now (Joel had set uw to 0 so it doesn't really cahnge anything if commented or not anyway)

     for (ii = 0; ii < 9; ii++ ){

	 K[ii]       = PHt[ii] / HPHtR;
	 kf->dX[ii]  = kf->update[ii] * K[ii] * kf->y[meas_i];

     }
     
     /* printf("-------------K---------------\n"); */
     /* printf("%e %e %e %e %e %e %e %e %e\n", K[0], K[1],K[2],K[3],K[4],K[5],K[6], K[7], K[8]); */
     // K*H is outer product
     for (ii = 0; ii < 9; ii ++ ){
	 for (jj = 0; jj < 9; jj++ ) {

	     I_KH[ii][jj] = -K[ii]*H[jj];

	     if (ii==jj) {

		 I_KH[ii][jj] += 1.0;

	     }

	 }
     }

     for (ii = 0; ii < 9; ii ++ ){
	 for (jj = 0; jj < 9; jj++ ) {

	   KRK[ii][jj] = K[ii]*K[jj] * kf->R[meas_i] ; //!!!! cbv: is that outer product too?

	 }
     }


     // Lazy filter update
     m_x_m9( Pnew, I_KH, kf->P);
     m_x_mt9( Pnew2, Pnew, I_KH);
     m_add9( kf->P , Pnew2, KRK);


     // !!!!!! cbv: why did Joel do that? This is basically using the typical updarta of coavria,ce, wehreas right above it's using Joseph.
     for (ii = 0 ; ii < 9; ii++ ){
	 for (jj = 0; jj < 9; jj++ ) {
	     kf->P[ii][jj] = Pnew[ii][jj];
	 }
     }
     // !!!!!! end of cbv: why did Joel do that?

     /*
     printf("-------------dx---------------\n");
     printf("%e %e %e %e %e %e %e \n", kf->dX[0], kf->dX[1],kf->dX[2],kf->dX[3],kf->dX[4],kf->dX[5],kf->dX[6]);
     */

     return 0;
 }


 int kalman_process_GPS( MEAS_T *meas, KALMAN_T *kf ) {

     // Declarations
     int ii;
     double H[9] = {0.0};



     // Cycle through the measurments
     for (ii = 0; ii < 6; ii++ ){

	 // Compute the residuals
	 kalman_compute_resids( meas, kf, 0, ii );

	 kalman_computeH_gps( H, kf, ii );

	 // Update State & Cov
	 kalman_joseph_update( kf, H , ii);

	 // Update State
	 kalman_update_state( kf );

	 // Post fit residuals
	 kalman_compute_resids( meas, kf, 1, ii );

     }



     //printf(" Acc %e \n ",  kf->sc.a_i2cg_LVLH[0] );


     return 0;

 }

 int kalman_update_state(KALMAN_T *kf ){

     //
     int ii;
     for (ii = 0; ii < 3; ii++ ){
	 kf->sc.r_i2cg_INRTL[ii]   += kf->dX[ii];
	 kf->sc.v_i2cg_INRTL[ii]   += kf->dX[ii+3];
	 kf->sc.a_i2cg_LVLH[ii]    += kf->dX[ii+6];
     }



     return 0;
 }


 int kalman_computeH_gps( double H[9], KALMAN_T *kf, int ii ){

     // Zero everything out
     int jj;
     for (jj = 0; jj < 9; jj++ ){
	 H[jj] = 0.0;
     }
     H[ii] = 1.0;

     /*
     printf("-------------K---------------\n");
     printf("H[ %i ] = %e %e %e %e %e %e %e \n", ii, H[0], H[1],H[2],H[3],H[4],H[5],H[6]);
     */

     return 0;
 }

 int kalman_compute_resids( MEAS_T *meas, KALMAN_T *kf, int flag, int index ) {

     // GPS position resids
     if (flag == 0 ){
	 if (index < 3 ) {
	     kf->y[index] = meas->r_i2cg_INRTL[index] - kf->sc.r_i2cg_INRTL[index];
/* 	     if (index == 1){ */
/* 	       printf("%e %e\n", kf->y[index]); */
	       
/* 	       //	       MPI_Finalize();exit(0); */
/* 	       } */
	 } else {
	     kf->y[index] = meas->v_i2cg_INRTL[index-3] - kf->sc.v_i2cg_INRTL[index-3];
	 }
     } else {
	 if (index < 3 ) {
	     kf->y_pf[index]  = meas->r_i2cg_INRTL[index] - kf->sc.r_i2cg_INRTL[index];
	 } else {
	     kf->y_pf[index]  = meas->v_i2cg_INRTL[index-3] - kf->sc.v_i2cg_INRTL[index-3];
	 }
     }




     /*
     printf(" ------------- Resids -----------\n");
     printf(" %g %g %g \n", kf->y[0], kf->y[1], kf->y[2]);
     */
     return 0;

 }


int kalman_init( MEAS_T *meas, KALMAN_T *kf, CONSTELLATION_T  *CONSTELLATION) {

     //
     int ii;
     int jj;
     double T[3][3]      = {{0.0}};
     double T9[9][9]     = {{0.0}};
     double Ptemp9[9][9] = {{0.0}};
     double Plvlh[9][9]  = {{0.0}};

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
     Plvlh[0][0] = 0.001*0.001;
     Plvlh[1][1] = 0.001*0.001;
     Plvlh[2][2] = 0.001*0.001;
     Plvlh[3][3] = 0.00001*0.00001;
     Plvlh[4][4] = 0.00001*0.00001;
     Plvlh[5][5] = 0.00001*0.00001;

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

     m_x_mt9( Ptemp9, Plvlh, T9);
     m_x_m9( kf->P, T9, Ptemp9);

     kf->P[6][6] = 1.0e-24;
     kf->P[7][7] = 1.0e-24;
     kf->P[8][8] = 1.0e-24;

     // !!!!!!!!!!!!!! to remove
     for (ii = 0; ii < 9; ii++){
       for(jj = 0; jj < 9; jj++ ) {
	 
	 kf->P[ii][jj] = 0.0;
       }
     }
     // !!!!!!!!!!!!!! end of to remove

     // Tau
     kf->tau =  95.0*60.0;

     // Sigma
     kf->sigma = 1e-18;

     // R
     kf->R[0] = 0.001 * 0.001; // 0.002 * 0.002 / (10000.0);
     kf->R[1] = 0.001 * 0.001; // 0.002 * 0.002 / (10000.0);
     kf->R[2] = 0.001 * 0.001; // 0.002 * 0.002 / (10000.0);
     kf->R[3] = 0.00001 * 0.00001; // 0.002 * 0.002 / (10000.0);
     kf->R[4] = 0.00001 * 0.00001; // 0.002 * 0.002 / (10000.0);
     kf->R[5] = 0.00001 * 0.00001; // 0.002 * 0.002 / (10000.0);
     /* // !!!!!!!!!! to remove */
     /* kf->R[0] = 0; // 0.002 * 0.002 / (10000.0); */
     /* kf->R[1] = 0; // 0.002 * 0.002 / (10000.0); */
     /* kf->R[2] = 0; // 0.002 * 0.002 / (10000.0); */
     /* kf->R[3] = 0; // 0.002 * 0.002 / (10000.0); */
     /* kf->R[4] = 0; // 0.002 * 0.002 / (10000.0); */
     /* kf->R[5] = 0; // 0.002 * 0.002 / (10000.0); */
     /* // !!!!!!!!!! end of to remove */


     kf->uw = 0;

     // Update array
     kf->update[0] = 1;
     kf->update[1] = 1;
     kf->update[2] = 1;
     kf->update[3] = 1;
     kf->update[4] = 1;
     kf->update[5] = 1;
     kf->update[6] = 1;
     kf->update[7] = 0;
     kf->update[8] = 0;
     
     // Finished completion
     kf->init_complete = 1;
     
     return 0;
 }

 int m_x_m9( double M[9][9], double M1[9][9], double M2[9][9] ){

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

 int m_x_mt9( double M[9][9], double M1[9][9], double M2[9][9] ){

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

 int m_add9( double M[9][9], double M1[9][9], double M2[9][9] ) {

     // Declarations
     int ii, jj;

     for (ii = 0 ; ii < 9; ii++ ){
	 for (jj = 0; jj < 9; jj++ ) {

	     M[ii][jj] = M1[ii][jj] + M2[ii][jj];

	 }
     }


     return 0;
 }

 int kalman_cov_prop( KALMAN_T *kf, PARAMS_T *PARAMS ) {

     // Declarations
     double PhiP[N_STATES][N_STATES] = {{0.0}};
     double Pnew[N_STATES][N_STATES] = {{0.0}};

     // Compute STM
     kalman_computeSTM( kf->Phi, &kf->sc, PARAMS, kf->sc.INTEGRATOR.dt, kf->tau );

     // Compute Q
     kalman_computeQ( kf->Q, kf->Phi, kf->sc.INTEGRATOR.dt, kf->tau, kf->sigma);

     // P = Phi * P * Phi^t + Q
     m_x_m9( PhiP, kf->Phi, kf->P);
     m_x_mt9( Pnew, PhiP, kf->Phi);
     m_add9( kf->P, Pnew, kf->Q ); // OK

     /*
     int jj;
     printf(" --------- P --------- \n");
     for (jj = 0; jj < 9; jj++ ){
	 printf("   %e %e %e %e %e %e %e \n", kf->P[jj][0] , kf->P[jj][1], kf->P[jj][2], kf->P[jj][3], kf->P[jj][4], kf->P[jj][5], kf->P[jj][6] );
     }
     */
     return 0;
 }

 int kalman_computeA( double A[9][9], SPACECRAFT_T *sc, PARAMS_T *PARAMS , double dt, double tau) {

     // Declarations
     double r3;
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

     // Intermediate
     v_mag( &rmag, sc->r_i2cg_INRTL );
     r3 = pow( rmag, 3 );
     r5 = pow( rmag, 5 );
     compute_T_inrtl_2_lvlh(T, sc->r_i2cg_INRTL, sc->v_i2cg_INRTL);
     m_trans(Tinv, T);

     for (ii = 0; ii < 9 ; ii++ ){
	 for (jj = 0; jj < 9; jj++ ){
	     A[ii][jj] = 0.0;
	 }
     }


     // Compute A
     for (ii = 0; ii < 3; ii++ ){

	 // dv/dv
	 A[ii][ii+3] = 1.0;

	 // dv/da
	 //A[ii][6] = Tinv[ii][0] * dt ;

	 // da/da
	 A[ii+3][6] = Tinv[ii][0];
	 A[ii+3][7] = Tinv[ii][1];
	 A[ii+3][8] = Tinv[ii][2];

	 // da/dr
	 //for (jj = 0; jj < 3 ; jj++ ){

	     //kf->A[jj+3][ii] = (3.0 * PARAMS->EARTH.GRAVITY.mu * kf->sc.r_i2cg_INRTL[ii] * kf->sc.r_i2cg_INRTL[jj] ) /r5;
	//     if (ii == jj ) {

		 //kf->A[jj+3][ii] -= PARAMS->EARTH.GRAVITY.mu / r3;

	//     }

	// }

	 // Cheating!!!  TODO: Fix this
	 v_copy( r, sc->r_i2cg_INRTL);
	 r[ii] = r[ii] + 0.0005;
	 compute_gravity( a2, r, sc->et, &PARAMS->EARTH.GRAVITY, sc->INTEGRATOR.degree, PARAMS->EARTH.earth_fixed_frame, PARAMS->EARTH.flattening, PARAMS->EARTH.radius);

	 v_copy( r, sc->r_i2cg_INRTL);
	 r[ii] = r[ii] - 0.0005;
	 //	 printf("%f\n", sc->INTEGRATOR.degree);
	 compute_gravity( a1, r, sc->et, &PARAMS->EARTH.GRAVITY, sc->INTEGRATOR.degree, PARAMS->EARTH.earth_fixed_frame, PARAMS->EARTH.flattening, PARAMS->EARTH.radius);

	 v_sub( a, a2, a1);

	 for (jj = 0; jj < 3; jj++ ){

	     a[jj] = a[jj] / 0.001;

	 }

	 A[ii+3][0] = a[0];
	 A[ii+3][1] = a[1];
	 A[ii+3][2] = a[2];


     }

     A[6][6] = -1 / tau;
     A[7][7] = -1 / tau;
     A[8][8] = -1 / tau;

     /*
     printf(" --------- A --------- \n");
     for (jj = 0; jj < 9; jj++ ){
	 printf("   %e %e %e %e %e %e %e \n", kf->A[jj][0], kf->A[jj][1],kf->A[jj][2],kf->A[jj][3],kf->A[jj][4],kf->A[jj][5],kf->A[jj][6] );
     }
     */
     return 0;
 }

 int kalman_computeSTM( double Phi[9][9], SPACECRAFT_T *sc, PARAMS_T *PARAMS, double dt, double tau ) {

     // Declarations
     int ii, jj;
     double A[9][9] = {{0.0}};
     double A2[9][9] = {{0.0}};

     for (ii = 0; ii < 9 ; ii++ ){
	 for (jj = 0; jj < 9; jj++ ){
	     Phi[ii][jj] = 0.0;
	 }
     }


     // Compute A
     kalman_computeA( A , sc, PARAMS, dt, tau );
     m_x_m9( A2, A, A );

     // Phi = I + A*dt
     for (ii = 0; ii < 9; ii++) {
	 for (jj = 0; jj<9; jj++) {

	     Phi[ii][jj] = A[ii][jj] * dt + A2[ii][jj] * dt *dt * 0.5;

	     if (ii == jj){

		 Phi[ii][jj] += 1.0;

	     }

	 }
     }




     /*
     printf(" --------- STM --------- \n");
     for (jj = 0; jj < 9; jj++ ){
	 printf("   %e %e %e %e %e %e %e \n", kf->Phi[jj][0], kf->Phi[jj][1],kf->Phi[jj][2],kf->Phi[jj][3],kf->Phi[jj][4],kf->Phi[jj][5],kf->Phi[jj][6] );
     }
     */
     return 0;
 }

 int kalman_computeQ( double Q[9][9], double Phi[N_STATES][N_STATES], double dt, double tau, double sigma) {

     int ii;
     int jj;

     double Qo[N_STATES][N_STATES]    = {{0.0}};
     double Qtemp[N_STATES][N_STATES] = {{0.0}};

     for (ii = 0; ii < N_STATES; ii++ ){
	 for (jj = 0; jj < N_STATES; jj++ ){

	     Qo[ii][jj] = 0.0;
	 }
     }

     for (ii = 0 ; ii < 3 ; ii++ ){

	 //Qo[ii][ii]       = sigma * pow( dt, 4) / 4.0;
	 //Qo[ii][ii+3]     = sigma * pow( dt, 3) / 2.0;
	 //Qo[ii+3][ii]     = sigma * pow( dt, 3) / 2.0;
	 //Qo[ii+3][ii+3]   = sigma * pow( dt, 2);

     }

     Qo[6][6] = sigma * ( 1.0 - exp( -2.0 * dt/ tau));
     Qo[7][7] = sigma * ( 1.0 - exp( -2.0 * dt/ tau));
     Qo[8][8] = sigma * ( 1.0 - exp( -2.0 * dt/ tau));

     m_x_mt9( Qtemp, Qo, Phi);
     m_x_m9( Q, Phi, Qtemp );

     /* // !!!!!!!!!!!!!! to remove */
     /* for (ii = 0; ii < 9; ii++){ */
     /*   for(jj = 0; jj < 9; jj++ ) { */
	 
     /* 	 Q[ii][jj] = 0.0; */
     /*   } */
     /* } */
     /* // !!!!!!!!!!!!!! end of to remove */


     return 0;
 }


int kalman_prop(  MEAS_T *meas, KALMAN_T *kf, PARAMS_T *PARAMS, OPTIONS_T *OPTIONS, GROUND_STATION_T *GROUND_STATION,CONSTELLATION_T  *CONSTELLATION, int iProc, int iDebugLevel, int nProcs ) {


     // Declarations
     double dt;

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



	 double starttime;
	 str2et_c(OPTIONS->initial_epoch, &starttime);
	 double density;
     while ( kf->sc.et != meas->et ) {

	 if (( meas->et - kf->sc.et  ) < kf->sc.INTEGRATOR.dt) {
	     kf->sc.INTEGRATOR.dt = meas->et - kf->sc.et;

	 }

	 kalman_cov_prop( kf, PARAMS );
	 /* etprint(kf->sc.et, "before"); */
	 /* v_print(kf->sc.r_i2cg_INRTL, "before"); */
	 /* v_print(kf->sc.v_i2cg_INRTL, "before"); */
	 
	 propagate_spacecraft( &kf->sc, PARAMS,starttime, &density, GROUND_STATION, OPTIONS, CONSTELLATION, iProc, iDebugLevel,  start_ensemble, array_sc );
	 /* etprint(kf->sc.et, "after"); */
	 /* v_print(kf->sc.r_i2cg_INRTL, "after"); */
	 /* v_print(kf->sc.v_i2cg_INRTL, "after"); */
	 //	 MPI_Finalize();exit(0);
	 //	 print_test();
	 //kf->sc.a_i2cg_LVLH[0] = exp( -kf->sc.INTEGRATOR.dt / kf->tau ) * kf->sc.a_i2cg_LVLH[0];
	 //kf->sc.a_i2cg_LVLH[1] = exp( -kf->sc.INTEGRATOR.dt / kf->tau ) * kf->sc.a_i2cg_LVLH[1];
	 //kf->sc.a_i2cg_LVLH[2] = exp( -kf->sc.INTEGRATOR.dt / kf->tau ) * kf->sc.a_i2cg_LVLH[2];

     }

     kf->sc.INTEGRATOR.dt = dt;


     return 0;

 }


 int kalman_write_out( KALMAN_T *kf, FILE *fp) {

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

     et = kf->sc.et;
     v_copy( r, kf->sc.r_i2cg_INRTL );
     v_copy( v, kf->sc.v_i2cg_INRTL );
     ad  = kf->sc.a_i2cg_LVLH[0];
     rho = kf->sc.rho;

     v2 = kf->sc.v_i2cg_INRTL[0]*kf->sc.v_i2cg_INRTL[0] + kf->sc.v_i2cg_INRTL[1]*kf->sc.v_i2cg_INRTL[1] + kf->sc.v_i2cg_INRTL[2]*kf->sc.v_i2cg_INRTL[2];

     Area = kf->sc.INTEGRATOR.surface[0].area / ( 1000.0 * 1000 ); // !!!!!! all surfaces !!!!! units ok?

     // cbv: accel = ax (without the minus sign)
     accel = 0.5 * rho * v2 * kf->sc.INTEGRATOR.surface[0].Cd * Area / kf->sc.INTEGRATOR.mass;// !!!!!! all surfaces 


     term = 2.0 * kf->sc.INTEGRATOR.mass / (kf->sc.INTEGRATOR.surface[0].Cd * Area );// !!!!!! all surfaces
     // cbv: density = term * accel / v^2
     // cbv: dd_dx seems to be d(rho)/dx (= 0?), d(rho)/dy (= 0?), d(rho)/dz (= 0?), d(rho)/dvx = (see dd_dx[3] below), d(rho)/dvy = (see dd_dx[4] below), d(rho)/dvz = (see dd_dx[5] below), d(rho)/dax = (see dd_dx[6] below, note that ax = accel), d(rho)/day = 0, d(rho)/daz = 0
    dd_dx[3] = -term * accel / (v2 * v2 ) * 2 * kf->sc.v_i2cg_INRTL[0];
    dd_dx[4] = -term * accel/ (v2 * v2 ) * 2 * kf->sc.v_i2cg_INRTL[1];
    dd_dx[5] = -term * accel/ (v2 * v2 ) * 2 * kf->sc.v_i2cg_INRTL[2];
    
    dd_dx[6] = term / v2;
    
    for (ii = 0; ii < N_STATES; ii++ ){
        prod[ii] = 0.0;
        for (jj = 0; jj< N_STATES; jj++ ) {
            prod[ii] += dd_dx[jj]*kf->P[jj][ii];
        }
    }

    drag_sigma = 0.0;
    for (ii = 0; ii < N_STATES; ii++ ){
    
        drag_sigma += prod[ii] * dd_dx[ii];
    
    }
    //printf("Sigma Drag %e \n", sqrt(drag_sigma) ) ;
  char text2[256];
 	et2utc_c(et, "ISOC", 6, 255, text2);


    fprintf(fp, "%s %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %e %e",
        text2, r[0], r[1], r[2], v[0], v[1], v[2], ad, rho);

/*     fprintf(fp, "%s %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f", */
/*         text2, r[0], r[1], r[2], v[0], v[1], v[2]); */

    fprintf(fp, " %e %e %e %e %e %e %e",
        kf->P[0][0], kf->P[1][1], kf->P[2][2], kf->P[3][3], kf->P[4][4], kf->P[5][5], kf->P[6][6]);
    
    fprintf(fp, " %e %e %e %e %e %e %e",
        kf->dX[0], kf->dX[1], kf->dX[2], kf->dX[3], kf->dX[4], kf->dX[5], kf->dX[6]);
    
    fprintf(fp, " %e %e %e",
        kf->y[0], kf->y[1], kf->y[2]);
    
    fprintf(fp, " %e %e %e",
        kf->sy[0], kf->sy[1], kf->sy[2]);
    
    fprintf(fp, " %e %e %e",
        kf->y_pf[0], kf->y_pf[1], kf->y_pf[2]);
    
    fprintf( fp, " %e", sqrt(drag_sigma) );

    
    fprintf(fp, " \n");


    return 0;

}






int read_meas(MEAS_T *meas, FILE *fp_meas, int init_complete, PARAMS_T *PARAMS){
  SpiceDouble       xform[6][6];
  double estate[6], jstate[6];

  int found_eoh;
  char *line = NULL;
  char text[256];
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
  } // end of if ( init_complete == 0 )

  getline(&line, &len, fp_meas);
  if (feof(fp_meas) == 0){
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




/*   char text2[256]; */
/*  	et2utc_c(meas->et, "ISOC", 6, 255, text2); */
/*   printf("<%s> <%s>\n", text, text2); */

  
  }

  
  return 0;
}
