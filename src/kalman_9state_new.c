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

#include "kalman_9state_new.h"


int kalman_filt( MEAS_T *meas, KALMAN_T *kf, PARAMS_T *params, OPTIONS_T *OPTIONS, GROUND_STATION_T *GROUND_STATION, CONSTELLATION_T  *CONSTELLATION, int iProc, int iDebugLevel, int nProcs ) {

  // read measurements
 
  read_meas(meas, meas->fp_meas, kf->init_complete, params);
  if (feof(meas->fp_meas) == 0){
    // Initialize if necessary
    if ( kf->init_complete == 0 ){
   

      kalman_init( meas, kf, CONSTELLATION,  params, OPTIONS, iProc, iDebugLevel );
      //      printf("Iniitalizing Nav \n");
      
      //      char fp_name[256];

      kf->fp_kalman = fopen(kf->sc.filenamekalman, "w+");

            kalman_write_out( kf, kf->fp_kalman) ;

	       return 0;
        
    }

    
    
     

    // Propagate State & Cov to Measurement time
    kalman_prop( meas, kf, params, OPTIONS, GROUND_STATION, CONSTELLATION, iProc, iDebugLevel, nProcs );

    //    kf->read_q = 0; // !!!!!! if Q is not to be recomputed at each time step (constant Q)

    // Process GP
    kalman_process_GPS( meas, kf);
     /* 	 int jj; */
     /* 	 printf(" --------- P after prop, after Kalman --------- \n"); */
     /* for (jj = 0; jj < 9; jj++ ){ */
     /* 	 printf("   %e %e %e %e %e %e \n", kf->P[jj][0] , kf->P[jj][1], kf->P[jj][2], kf->P[jj][3], kf->P[jj][4], kf->P[jj][5] ); */
     /* } */

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
	 //	 printf("K[%d] = %e\n", ii, K[ii]);

     }
     
     /* printf("------------- K ---------------\n"); */
     /* printf("%e %e %e %e %e %e\n", K[0], K[1],K[2],K[3],K[4],K[5]); */
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
	   //	   printf("K[%d] = %e * K[%d] = %e * R[%d] = %e -> %e\n",ii, K[ii],jj, K[jj],meas_i, kf->R[meas_i], KRK[ii][jj]);
	 }
     }


     // Lazy filter update
     m_x_m9_for_kalman( Pnew, I_KH, kf->P);
     m_x_mt9_for_kalman( Pnew2, Pnew, I_KH);
     m_add9_for_kalman( kf->P , Pnew2, KRK);

     
     /* 	 printf(" --------- P joseph - P not joseph --------- \n"); */
     /* for (jj = 0; jj < 9; jj++ ){ */
     /* 	 printf("   %e %e %e %e %e %e \n", kf->P[jj][0] - Pnew[jj][0], kf->P[jj][1] - Pnew[jj][1], kf->P[jj][2] - Pnew[jj][2], kf->P[jj][3] - Pnew[jj][3], kf->P[jj][4] - Pnew[jj][4], kf->P[jj][5] - Pnew[jj][5] ); */
     /* } */


     /* // !!!!!! cbv: why did Joel do that? This is basically using the typical updarta of coavria,ce, wehreas right above it's using Joseph. */
     /* for (ii = 0 ; ii < 9; ii++ ){ */
     /* 	 for (jj = 0; jj < 9; jj++ ) { */
     /* 	     kf->P[ii][jj] = Pnew[ii][jj]; */
     /* 	 } */
     /* } */
     /* // !!!!!! end of cbv: why did Joel do that? */

     /* 	 printf(" --------- I_KH --------- \n"); */
     /* for (jj = 0; jj < 9; jj++ ){ */
     /* 	 printf("   %e %e %e %e %e %e \n", I_KH[jj][0], I_KH[jj][1], I_KH[jj][2] , I_KH[jj][3], I_KH[jj][4] , I_KH[jj][5]  ); */
     /* } */


     
     /* 	 printf(" --------- KRK --------- \n"); */
     /* for (jj = 0; jj < 9; jj++ ){ */
     /* 	 printf("   %e %e %e %e %e %e \n", KRK[jj][0], KRK[jj][1], KRK[jj][2] , KRK[jj][3], KRK[jj][4] , KRK[jj][5]  ); */
     /* } */


     
     /*
     printf("-------------dx---------------\n");
     printf("%e %e %e %e %e %e %e \n", kf->dX[0], kf->dX[1],kf->dX[2],kf->dX[3],kf->dX[4],kf->dX[5],kf->dX[9]);
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
	 kalman_update_state( kf,ii );

	 // Post fit residuals
	 kalman_compute_resids( meas, kf, 1, ii );

     }



     //printf(" Acc %e \n ",  kf->sc.a_i2cg_LVLH[0] );


     return 0;

 }

int kalman_update_state(KALMAN_T *kf, int imeas ){

     //
     int ii;
/*      for (ii = 3; ii < 6; ii++){ */
/*        //       printf("dX[%d]: %e %e | ", ii, kf->dX[ii], kf->sc.a_i2cg_LVLH[ii-6]); */
/*        printf("dX[%d]: %e %e %d| ", ii, kf->dX[ii], kf->sc.v_i2cg_INRTL[ii-3], imeas); */
/*      } */
/*      printf("\n"); */

     for (ii = 0; ii < 3; ii++ ){
	 kf->sc.r_i2cg_INRTL[ii]   += kf->dX[ii];
	 kf->sc.v_i2cg_INRTL[ii]   += kf->dX[ii+3];
	 kf->sc.a_i2cg_LVLH[ii]    += kf->dX[ii+6];
/* 	 printf("r[%d]: %e | ",ii, kf->dX[ii]/kf->sc.r_i2cg_INRTL[ii]); */
/* 	   printf("v[%d]: %e | ",ii, kf->dX[ii+3]/kf->sc.v_i2cg_INRTL[ii]); */
//	 printf("a[%d]: %e (%e) | ",ii, kf->dX[ii+6]/kf->sc.a_i2cg_INRTL[ii], kf->dX[ii+6]);
     }
     //       printf("\n");
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


int kalman_init( MEAS_T *meas, KALMAN_T *kf, CONSTELLATION_T  *CONSTELLATION,       PARAMS_T *PARAMS, OPTIONS_T *OPTIONS, int iProc, int iDebugLevel) {
  char *line = NULL;
  char text[256];
  size_t len = 0;

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
     m_x_m9_for_kalman( kf->P, T9, Ptemp9);

     kf->P[6][6] = 1.0e-24;
     kf->P[7][7] = 1.0e-24;
     kf->P[8][8] = 1.0e-24;

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

       // calculate the acceleration
       double v_dummy[3];
       double test_density;
       double a_inrtl[3];
       double et_initial_epoch;
       str2et_c(OPTIONS->initial_epoch, &et_initial_epoch);
       compute_dxdt( v_dummy,a_inrtl, &kf->sc.et,kf->sc.r_i2cg_INRTL,kf->sc.v_i2cg_INRTL, PARAMS, &kf->sc.INTEGRATOR, et_initial_epoch, OPTIONS->et_oldest_tle_epoch, &test_density, kf->sc.INTEGRATOR.index_in_attitude_interpolated, kf->sc.INTEGRATOR.index_in_driver_interpolated, CONSTELLATION, OPTIONS, iProc, iDebugLevel, &kf->sc);

       m_x_v(kf->sc.a_i2cg_LVLH, T, kf->sc.a_i2cg_INRTL_drag);
       //     kf->sc.a_i2cg_LVLH[0] = 0; kf->sc.a_i2cg_LVLH[1] = 0; kf->sc.a_i2cg_LVLH[2] = 0;  // a_i2cg_LVLH is supposed ot represent the unmodeled accleration, not the total acceleration. total acceleration = a_gravity + a_unmodeled and kf->sc.a_i2cg_LVLH is an estimation of a_unmodeled



     // Tau and sigma !!!!!!! joel ussed to call his sigma my sigma**2
                     getline(&line, &len, kf->fp_kalman_init);
              getline(&line, &len, kf->fp_kalman_init);
       sscanf(line, "%lf", &kf->tau);
              getline(&line, &len, kf->fp_kalman_init);
       sscanf(line, "%lf", &kf->sigma);

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
     

     // Finished completion
     kf->init_complete = 1;

     // Q is set to 1 to get in the function computeQ
     //kf->read_q = 1; 

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

int kalman_cov_prop( KALMAN_T *kf, PARAMS_T *PARAMS, OPTIONS_T *OPTIONS,  CONSTELLATION_T *CONSTELLATION, int iProc, int iDebugLevel ) {

     // Declarations
     double PhiP[N_STATES][N_STATES] = {{0.0}};
     double Pnew[N_STATES][N_STATES] = {{0.0}};


     // Compute STM

     kalman_computeSTM( kf->Phi, &kf->sc, PARAMS, kf->sc.INTEGRATOR.dt, kf->tau, OPTIONS,  CONSTELLATION,iProc, iDebugLevel );

     // Compute Q

     //     if (kf->read_q == 1){
      
     kalman_computeQ( kf->Q, kf->Phi, kf->sc.INTEGRATOR.dt, kf->tau, kf->sigma, kf->fp_kalman_init);
     //     }
     /* int jj; */
     /* printf(" --------- Q  --------- \n"); */
     /* for (jj = 0; jj < 9; jj++ ){ */
     /* 	 printf("   %e %e %e %e %e %e \n", kf->Q[jj][0] , kf->Q[jj][1], kf->Q[jj][2], kf->Q[jj][3], kf->Q[jj][4], kf->Q[jj][5] ); */
     /* } */


     
     // P = Phi * P * Phi^t + Q
     m_x_m9_for_kalman( PhiP, kf->Phi, kf->P);
     m_x_mt9_for_kalman( Pnew, PhiP, kf->Phi);
     m_add9_for_kalman( kf->P, Pnew, kf->Q ); // OK

     
     /*      int jj; */
     /* printf(" --------- P after prop, before Kalman --------- \n"); */
     /* for (jj = 0; jj < 9; jj++ ){ */
     /* 	 printf("   %e %e %e %e %e %e \n", kf->P[jj][0] , kf->P[jj][1], kf->P[jj][2], kf->P[jj][3], kf->P[jj][4], kf->P[jj][5] ); */
     /* } */
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
     double scale_h = 200.0; // in km, from table A-1 of cruisckshank98
     double r0 = 7298.1450; // in km, from table A-1 of cruisckshank98
     double rho_0 = 0.0004; // kg/km^3 , from table A-1 of cruisckshank98
     double rho_exp;
     double omega_e = 0.00007292158553; // rad/s
     double va[3]; // speed of sc relative to rotating atmo
     va[0] = sc->v_i2cg_INRTL[0] + omega_e * sc->r_i2cg_INRTL[1];
     va[1] = sc->v_i2cg_INRTL[1] - omega_e * sc->r_i2cg_INRTL[0];
     va[2] = sc->v_i2cg_INRTL[2];
     double vamag;
     v_mag( &vamag, va );
     rho_exp = rho_0 * exp( - ( rmag - r0 ) / scale_h);
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

     // dax / dx
     A[3][0] = -mu / r3 *
       ( 1 - 3 * x2 / r2 + 15. / 2 * j2 * re2 * z2 / r4 * ( 7 * x2 / r2 - 1 ) + 3./2 * j2 * re2 / r2 * ( -5 * x2 / r2 + 1 ) )
       + 1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp * ( xdot + omega_e * y) * ( omega_e * ( ydot - omega_e * x ) / vamag + vamag * x / ( scale_h * rmag ) );

     // dax / dy
     A[3][1] = 3 * mu * x * y / r5 *
       ( 1 - 35./2 * j2 * re2 * z2 / r4 + 5./2 * j2 * re2/r2 )
       -1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp * ( vamag * omega_e + ( xdot + omega_e * y ) * ( ( xdot + omega_e * y ) * omega_e / vamag - vamag * y / ( scale_h * rmag) ) );

     // dax / dz
     A[3][2] = 3 * mu * x * z / r5 *
       ( 1 - 35./2 * j2 * re2 * z2 / r4 + 15./2 * j2 * re2/r2 )
       + 1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp * vamag * z * ( xdot + omega_e * y ) / ( scale_h * rmag);

     // dax / dvx
     A[3][3] =  -1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp * ( vamag + ( xdot + omega_e * y ) * ( xdot + omega_e * y ) / vamag );

     // dax / dvy
     A[3][4] = -1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp / vamag * ( xdot + omega_e * y ) * ( ydot - omega_e * x );

     // dax / dvz
     A[3][5] = -1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp * zdot / vamag * ( xdot + omega_e * y ) ;

     // day / dx
     A[4][0] = 3 * mu * x * y / r5 *
       ( 1 - 35./2 * j2 * re2 * z2 / r4 + 5./2 * j2 * re2/r2 )
       -1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp * ( -vamag * omega_e - ( ydot - omega_e * x ) * ( ( ydot - omega_e * x ) * omega_e / vamag + vamag * x / ( scale_h * rmag) ) );
     ;

     // day / dy
     A[4][1] =  -mu / r3 *
       ( 1 - 3 * y2 / r2 + 15. / 2 * j2 * re2 * z2 / r4 * ( 7 * y2 / r2 - 1 ) + 3./2 * j2 * re2 / r2 * ( -5 * y2 / r2 + 1 ) )
       - 1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp * ( ydot - omega_e * x ) * ( omega_e * ( xdot + omega_e * y ) / vamag - vamag * y / ( scale_h * rmag ) );

          // day / dz
     A[4][2] =  3 * mu * y * z / r5 *
       ( 1 - 35./2 * j2 * re2 * z2 / r4 + 15./2 * j2 * re2/r2 )
       + 1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp * vamag * z * ( ydot - omega_e * x ) / ( scale_h * rmag);

          // day / dvx
     A[4][3] =  -1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp / vamag * ( xdot + omega_e * y ) * ( ydot - omega_e * x );

               // day / dvy
     A[4][4] = -1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp * ( vamag + ( ydot - omega_e * x ) * ( ydot - omega_e * x ) / vamag );

          // day / dvz
     A[4][5] = -1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp * zdot / vamag * ( ydot - omega_e * x ) ;;

               // daz / dx
     A[5][0] = 3 * mu * x * z / r5 *
       ( 1 + 15./2 * j2 * re2/r2 - 35./2 * j2 * re2 * z2 / r4  )
       + 1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp * zdot * ( omega_e * ( ydot - omega_e * x ) / vamag + vamag * x / ( scale_h * rmag) ) ;

     // daz / dy
     A[5][1] = 3 * mu * y * z / r5 *
       ( 1 + 15./2 * j2 * re2/r2 - 35./2 * j2 * re2 * z2 / r4  )
       - 1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp * zdot * ( omega_e * ( xdot + omega_e * y ) / vamag - vamag * y / ( scale_h * rmag ) ) ;

     // daz / dz
     A[5][2] =  -mu / r3 *
       ( 1 - 3 * z2 / r2 + 3. * j2 * re2 / r2 * ( 3./2 - 15. * z2 / r2 + 35./2 * z2 / r4 ) )
       + 1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp * vamag * z * zdot / ( scale_h * rmag ) ;

     // daz / dvx
     A[5][3] =  -1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp * zdot * (  xdot + omega_e * y  ) / vamag;

     // daz / dvy
     A[5][4] = -1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp * zdot * ( ydot - omega_e * x ) / vamag;

     // daz / dvz
     A[5][5] =  -1./2 * sc->INTEGRATOR.surface[0].Cd * sc->INTEGRATOR.surface[0].area / m * rho_exp * ( zdot * zdot / vamag + vamag ) ;
     
     
     /* printf(" --------- A --------- \n"); */
     /* for (jj = 0; jj < 6; jj++ ){ */
     /*   printf("   %e %e %e %e %e %e %e %e %e\n", A[jj][0], A[jj][1],A[jj][2],A[jj][3],A[jj][4],A[jj][5], A[jj][6],A[jj][7],A[jj][8] ); */
     /* } */
     return 0;
 }

int kalman_computeSTM( double Phi[9][9], SPACECRAFT_T *sc, PARAMS_T *PARAMS, double dt, double tau,OPTIONS_T *OPTIONS,  CONSTELLATION_T *CONSTELLATION, int iProc, int iDebugLevel ) {


     int ii, jj;
     for (ii = 0; ii < 9 ; ii++ ){
	 for (jj = 0; jj < 9; jj++ ){
	     Phi[ii][jj] = 0.0;
	 }
     }

     // use equation 19 of cruisckshank98: determine first phif from A of appendix a of cruisckshank98 (A is calculated in kalman_computeA using appendix A of cruisckshank98). The calculate phiw and M (using notation of eq 19 of cruisckshank98)
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

     // Combine phif, phiw, and M to get Phi
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

     
     /* printf(" --------- A --------- \n"); */
     /* for (jj = 0; jj < 6; jj++ ){ */
     /* 	 printf("   %e %e %e %e %e %e \n", A[jj][0], A[jj][1],A[jj][2],A[jj][3],A[jj][4],A[jj][5] ); */
     /* } */
     
     /* printf(" --------- A2 --------- \n"); */
     /* for (jj = 0; jj < 6; jj++ ){ */
     /* 	 printf("   %e %e %e %e %e %e \n", A2[jj][0], A2[jj][1],A2[jj][2],A2[jj][3],A2[jj][4],A2[jj][5] ); */
     /* } */
     
     /* printf(" --------- STM --------- \n"); */
     /* for (jj = 0; jj < 9; jj++ ){ */
     /*   printf("   %e %e %e %e %e %e %e %e %e\n", Phi[jj][0], Phi[jj][1],Phi[jj][2],Phi[jj][3],Phi[jj][4],Phi[jj][5],Phi[jj][6],Phi[jj][7],Phi[jj][8] ); */
     /* } */

     return 0;
 }

int kalman_computeQ( double Q[9][9], double Phi[N_STATES][N_STATES], double dt, double tau, double sigma, FILE *fp_kalman_init) {
  char *line = NULL;
  char text[256];
  size_t len = 0;

     int ii;
     int jj;

     for (ii = 0; ii < 9; ii++ ){
	 for (jj = 0; jj < 9; jj++ ){

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

/*      /\*        getline(&line, &len, fp_kalman_init); *\/ */
/*      /\* for (ii = 0; ii < 9; ii++){ *\/ */
/*      /\*   getline(&line, &len, fp_kalman_init); *\/ */
/*      /\*   RemoveSpaces(line);  strtok(line, "\n");  strtok(line, "\r");  *\/ */

/*      /\*   if (ii == 0){ *\/ */
/*      /\* 	 sscanf(line, "((%lf;%lf;%lf;%lf;%lf;%lf);", &Q[ii][0], &Q[ii][1], &Q[ii][2], &Q[ii][3], &Q[ii][4], &Q[ii][5]); *\/ */
/*      /\*   } *\/ */
/*      /\*   else if (ii == 5){ *\/ */
/*      /\* 	 sscanf(line, "(%lf;%lf;%lf;%lf;%lf;%lf))", &Q[ii][0], &Q[ii][1], &Q[ii][2], &Q[ii][3], &Q[ii][4], &Q[ii][5]); *\/ */
/*      /\*   } *\/ */
/*      /\*   else{ *\/ */
/*      /\* 	 sscanf(line, "(%lf;%lf;%lf;%lf;%lf;%lf);", &Q[ii][0], &Q[ii][1], &Q[ii][2], &Q[ii][3], &Q[ii][4], &Q[ii][5]); *\/ */
/*      /\*   } *\/ */
/*      /\* } *\/ */
/*      /\* // !!!!!!!!!!!!!! to remove *\/ */
/*      /\* for (ii = 0; ii < 9; ii++){ *\/ */
/*      /\*   for(jj = 0; jj < 9; jj++ ) { *\/ */
	 
/*      /\* 	 printf("%e ", Q[ii][jj]); *\/ */
/*      /\*   } *\/ */
/*      /\*   printf("\n"); *\/ */
/*      /\* } *\/ */
/*      /\* // !!!!!!!!!!!!!! end of to remove *\/ */



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



/*      // !!!!!!!!! Comment block below */
/* /\*      double S[3][3]; *\/ */
/* /\*      for (ii = 0; ii < 3; ii++ ){ *\/ */
/* /\*        for (jj = 0; jj < 3; jj++ ){ *\/ */
/* /\*      	 S[ii][jj] = 0; *\/ */
/* /\*        } *\/ */
/* /\*      } *\/ */
/* /\*      S[0][0] = sigma * ( 1.0 - exp( -2.0 * dt/ tau)); *\/ */
/* /\*      S[1][1] = sigma * ( 1.0 - exp( -2.0 * dt/ tau)); *\/ */
/* /\*      S[2][2] = sigma * ( 1.0 - exp( -2.0 * dt/ tau)); *\/ */

/* /\*      for (ii = 0; ii < 3; ii++){ *\/ */
/* /\*        for (jj = 0; jj < 3; jj++ ){ *\/ */
/* /\*      	 Q[ii][jj] = 1./4 * dt*dt*dt*dt * S[ii][jj]; *\/ */
/* /\*      	 Q[ii][jj+3] = 1./2 * dt*dt*dt * S[ii][jj]; *\/ */
/* /\*      	 Q[ii][jj+6] =  1./2 * dt*dt * S[ii][jj]; *\/ */
/* /\*      	 Q[ii+3][jj] = 1./2 * dt*dt*dt * S[ii][jj]; *\/ */
/* /\*      	 Q[ii+3][jj+3] = dt*dt * S[ii][jj]; *\/ */
/* /\*      	 Q[ii+3][jj+6] = dt * S[ii][jj]; *\/ */
/* /\*      	 Q[ii+6][jj] = 1./2 * dt*dt * S[ii][jj]; *\/ */
/* /\*      	 Q[ii+6][jj+3] = dt * S[ii][jj]; *\/ */
/* /\*      	 Q[ii+6][jj+6] = S[ii][jj]; *\/ */
/* /\*        } *\/ */
/* /\*      } *\/ */
/*      // !!!!!!!!! end of comment block below */



     
/*      /\* // !!!!!!!!!!!!!! to remove *\/ */
/*      /\* for (ii = 0; ii < 9; ii++){ *\/ */
/*      /\*   for(jj = 0; jj < 9; jj++ ) { *\/ */
	 
/*      /\* 	 printf("%e ", Q[ii][jj]); *\/ */
/*      /\*   } *\/ */
/*      /\*   printf("\n"); *\/ */
/*      /\* } *\/ */
/*      /\* printf("\n"); *\/ */
/*      /\* // !!!!!!!!!!!!!! end of to remove *\/ */



/* /\*      Q[0][0] = (10./1000)*(10./1000); *\/ */
/* /\*      Q[1][1] = (10./1000)*(10./1000); *\/ */
/* /\*      Q[2][2] = (10./1000)*(10./1000); *\/ */
/* /\*      Q[3][3] = (1./100000)*(1./100000); *\/ */
/* /\*      Q[4][4] = (1./100000)*(1./100000); *\/ */
/* /\*      Q[5][5] = (1./100000)*(1./100000); *\/ */

/* /\*      /\\* // !!!!!!!!!!!!!! to remove *\\/ *\/ */
/* /\*      /\\* for (ii = 0; ii < 9; ii++){ *\\/ *\/ */
/* /\*      /\\*   for(jj = 0; jj < 9; jj++ ) { *\\/ *\/ */
	 
/* /\*      /\\* 	 Q[ii][jj] = 0.0; *\\/ *\/ */
/* /\*      /\\*   } *\\/ *\/ */
/* /\*      /\\* } *\\/ *\/ */
/* /\*      /\\* // !!!!!!!!!!!!!! end of to remove *\\/ *\/ */


     return 0;
 }


int kalman_prop(  MEAS_T *meas, KALMAN_T *kf, PARAMS_T *PARAMS, OPTIONS_T *OPTIONS, GROUND_STATION_T *GROUND_STATION,CONSTELLATION_T  *CONSTELLATION, int iProc, int iDebugLevel, int nProcs ) {

  double a_i2cg_LVLH_save[3];
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

	 kalman_cov_prop( kf, PARAMS, OPTIONS,  CONSTELLATION ,  iProc,  iDebugLevel );
	 /* etprint(kf->sc.et, "before"); */
	 /* v_print(kf->sc.r_i2cg_INRTL, "before"); */
	 /* v_print(kf->sc.v_i2cg_INRTL, "before"); */
	 v_copy(a_i2cg_LVLH_save, kf->sc.a_i2cg_LVLH); // we don't wannt kf->sc.a_i2cg_LVLH to be calculated in the propagator because it's an unmodeled acceleration
	 propagate_spacecraft( &kf->sc, PARAMS,starttime, OPTIONS->et_oldest_tle_epoch, &density, GROUND_STATION, OPTIONS, CONSTELLATION, iProc, iDebugLevel,  start_ensemble, array_sc );

	 /* etprint(kf->sc.et, "after"); */
	 /* v_print(kf->sc.r_i2cg_INRTL, "after"); */
	 /* v_print(kf->sc.v_i2cg_INRTL, "after"); */
	 //	 MPI_Finalize();exit(0);
	 //	 print_test();


	 kf->sc.a_i2cg_LVLH[0] = exp( -kf->sc.INTEGRATOR.dt / kf->tau ) * a_i2cg_LVLH_save[0];
	 kf->sc.a_i2cg_LVLH[1] = exp( -kf->sc.INTEGRATOR.dt / kf->tau ) * a_i2cg_LVLH_save[1];
	 kf->sc.a_i2cg_LVLH[2] = exp( -kf->sc.INTEGRATOR.dt / kf->tau ) * a_i2cg_LVLH_save[2];
	 
/* 	 kf->sc.a_i2cg_LVLH[0] = exp( -kf->sc.INTEGRATOR.dt / kf->tau ) * kf->sc.a_i2cg_LVLH[0]; */
/* 	 kf->sc.a_i2cg_LVLH[1] = exp( -kf->sc.INTEGRATOR.dt / kf->tau ) * kf->sc.a_i2cg_LVLH[1]; */
/* 	 kf->sc.a_i2cg_LVLH[2] = exp( -kf->sc.INTEGRATOR.dt / kf->tau ) * kf->sc.a_i2cg_LVLH[2]; */
	 //	 	 v_print(kf->sc.a_i2cg_LVLH, "kf->sc.a_i2cg_LVLH");
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
     double rho_estimate;


     ad  = kf->sc.a_i2cg_LVLH[0];
     rho = kf->sc.rho;

     v2 = kf->sc.v_i2cg_INRTL[0]*kf->sc.v_i2cg_INRTL[0] + kf->sc.v_i2cg_INRTL[1]*kf->sc.v_i2cg_INRTL[1] + kf->sc.v_i2cg_INRTL[2]*kf->sc.v_i2cg_INRTL[2];


     Area = kf->sc.INTEGRATOR.surface[0].area / ( 1000.0 * 1000 ); // !!!!!! all surfaces !!!!!!! units ok?

     double T[3][3], T_inv[3][3];
     double v_estimate_lvlh[3];
     compute_T_inrtl_2_lvlh(T, kf->sc.r_i2cg_INRTL, kf->sc.v_i2cg_INRTL);
     m_x_v(v_estimate_lvlh, T, kf->sc.v_i2cg_INRTL);
       

     double amag_estimate, vmag_estimate;
     v_mag(&amag_estimate, kf->sc.a_i2cg_LVLH); // a_i2cg_LVLH is supposed ot represent the unmodeled accleration
     v_mag(&vmag_estimate, v_estimate_lvlh);
     //     double a_drag_estimate[3];// a_drag_estimate  = a_total - a_gravity (only 2 forces: drag and gravity)
     rho_estimate = fabs(kf->sc.a_i2cg_LVLH[0]/(1/2. * v_estimate_lvlh[0] * vmag_estimate * kf->sc.INTEGRATOR.surface[0].Cd * kf->sc.INTEGRATOR.surface[0].area / kf->sc.INTEGRATOR.mass))/(1e9); // in kg/m^3
     //     rho_estimate = amag_estimate/(1/2. * v2 * kf->sc.INTEGRATOR.surface[0].Cd * kf->sc.INTEGRATOR.surface[0].area / kf->sc.INTEGRATOR.mass);
     //printf("%e %e %e %e %e %e\n", rho_estimate, amag_estimate, v2, kf->sc.INTEGRATOR.surface[0].Cd, kf->sc.INTEGRATOR.surface[0].area, kf->sc.INTEGRATOR.mass);
 
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


    /* fprintf(fp, "%s %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %e %e", */
    /*     text2, r[0], r[1], r[2], v[0], v[1], v[2], ad, rho); */

	fprintf(fp, "%s %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f",
        text2, r[0], r[1], r[2], v[0], v[1], v[2]);

/*     fprintf(fp, "%s %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f", */
/*         text2, r[0], r[1], r[2], v[0], v[1], v[2]); */

        fprintf(fp, " %e %e %e %e %e %e",
        kf->P[0][0], kf->P[1][1], kf->P[2][2], kf->P[3][3], kf->P[4][4], kf->P[5][5]);

    /* fprintf(fp, " %e %e %e %e %e %e %e", */
    /*     kf->P[0][0], kf->P[1][1], kf->P[2][2], kf->P[3][3], kf->P[4][4], kf->P[5][5], kf->P[6][6]); */

    fprintf(fp, " %e %e %e %e %e %e",
        kf->dX[0], kf->dX[1], kf->dX[2], kf->dX[3], kf->dX[4], kf->dX[5]);

	
    /* fprintf(fp, " %e %e %e %e %e %e %e", */
    /*     kf->dX[0], kf->dX[1], kf->dX[2], kf->dX[3], kf->dX[4], kf->dX[5], kf->dX[6]); */
    
    fprintf(fp, " %e %e %e",
        kf->y[0], kf->y[1], kf->y[2]);
    
    fprintf(fp, " %e %e %e",
        kf->sy[0], kf->sy[1], kf->sy[2]);
    
    fprintf(fp, " %e %e %e",
        kf->y_pf[0], kf->y_pf[1], kf->y_pf[2]);

    fprintf(fp, " %e %e %e",
	    kf->sc.a_i2cg_LVLH[0], kf->sc.a_i2cg_LVLH[1], kf->sc.a_i2cg_LVLH[2]);
    //        v_print(kf->sc.a_i2cg_LVLH,"");

    fprintf(fp," %e %e %e %e %e", ad, rho_estimate, kf->dX[6], kf->P[6][6], sqrt(drag_sigma)); 

    //   fprintf( fp, " %e", sqrt(drag_sigma) );

    
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
  /* double meas_r_ecef[3]; */
  /* double meas_v_ecef[3]; */
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
      sscanf(line, "%s %lf %lf %lf %lf %lf %lf", text, &meas->r_i2cg_INRTL[0],  &meas->r_i2cg_INRTL[1],  &meas->r_i2cg_INRTL[2], &meas->v_i2cg_INRTL[0],  &meas->v_i2cg_INRTL[1],  &meas->v_i2cg_INRTL[2]);
  str2et_c(text, &meas->et);
 
  /*   // !!!!!!!!! the file I use from the netcdf data is in ECEF so need to convert to ECI */
   
  /*   sscanf(line, "%s %lf %lf %lf %lf %lf %lf", text, &meas_r_ecef[0],  &meas_r_ecef[1],  &meas_r_ecef[2], &meas_v_ecef[0],  &meas_v_ecef[1],  &meas_v_ecef[2]); */
  /* str2et_c(text, &meas->et);  */
  /* // // Convert ECEF to ECI */


  /* 	    estate[0] = meas_r_ecef[0];estate[1] = meas_r_ecef[1];estate[2] = meas_r_ecef[2]; */
  /* 	    estate[3] = meas_v_ecef[0];estate[4] = meas_v_ecef[1];estate[5] = meas_v_ecef[2]; */
  /* 	    sxform_c (  PARAMS->EARTH.earth_fixed_frame,  "J2000", meas->et,    xform  );  */
  /* 	    mxvg_c   (  xform,       estate,   6,  6, jstate );  */
  /* 	    meas->r_i2cg_INRTL[0] = jstate[0]; meas->r_i2cg_INRTL[1] = jstate[1]; meas->r_i2cg_INRTL[2] = jstate[2]; */
  /* 	    meas->v_i2cg_INRTL[0] = jstate[3]; meas->v_i2cg_INRTL[1] = jstate[4]; meas->v_i2cg_INRTL[2] = jstate[5]; */




/*   char text2[256]; */
/*  	et2utc_c(meas->et, "ISOC", 6, 255, text2); */
/*   printf("<%s> <%s>\n", text, text2); */

  
  }


  
  return 0;
}
