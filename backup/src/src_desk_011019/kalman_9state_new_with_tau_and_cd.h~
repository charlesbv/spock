#ifndef __kalman_h__
#define __kalman_h__

#include "propagator.h"
#include "options.h"
#define N_STATES 10

typedef struct {
  double r_i2cg_INRTL[3];
  double v_i2cg_INRTL[3];
  double et;
  FILE *fp_meas;
} MEAS_T;

typedef struct {  

  double X[10];
  double P[10][10];
  double dX[10];
  double update[10];

  double Phi[10][10];
  double Q[10][10];

    double tau;
 double sigma_tau;
    double y[6];
    double y_pf[6];
    double sy[6];
    
    double sigma;
    
    int    init_complete;

    double R[6];
    double uw;
    
    SPACECRAFT_T sc;
     FILE *fp_kalman;
  FILE *fp_kalman_init;    
  int read_q;
} KALMAN_T;

#endif

int kalman_cov_prop( KALMAN_T *kf, PARAMS_T *PARAMS, OPTIONS_T *OPTIONS,  CONSTELLATION_T *CONSTELLATION, int iProc, int iDebugLevel );
int kalman_filt( MEAS_T *meas, KALMAN_T *kf, PARAMS_T *params, OPTIONS_T *OPTIONS, GROUND_STATION_T *GROUND_STATION, CONSTELLATION_T  *CONSTELLATION, int iProc, int iDebugLevel, int nProcs );
int kalman_init( MEAS_T *meas, KALMAN_T *kf, CONSTELLATION_T  *CONSTELLATION    ,   PARAMS_T *PARAMS, OPTIONS_T *OPTIONS, int iProc, int iDebugLevel);
int kalman_prop( MEAS_T *meas, KALMAN_T *kf, PARAMS_T *params, OPTIONS_T *OPTIONS,GROUND_STATION_T *GROUND_STATION, CONSTELLATION_T  *CONSTELLATION, int iProc, int iDebugLevel, int nProcs );
int kalman_compute_resids( MEAS_T *meas, KALMAN_T *kf, int flag, int index );
int kalman_computeA( double A[6][6], SPACECRAFT_T *sc, PARAMS_T *PARAMS , double dt, double tau, OPTIONS_T *OPTIONS,  CONSTELLATION_T *CONSTELLATION, int iProc, int iDebugLevel);
int kalman_computeSTM( double Phi[10][10], SPACECRAFT_T *sc, PARAMS_T *PARAMS, double dt, double tau,OPTIONS_T *OPTIONS,  CONSTELLATION_T *CONSTELLATION, int iProc, int iDebugLevel );
int kalman_computeQ( double Q[10][10], double Phi[N_STATES][N_STATES], double dt, double tau, double sigma, FILE *fp_kalman_init, double sigma_tau );
int kalman_process_GPS( MEAS_T *meas, KALMAN_T *kf ) ;
int kalman_computeH_gps( double H[10], KALMAN_T *kf, int ii );
int kalman_joseph_update(  KALMAN_T *kf, double H[10], int meas_i);
int kalman_update_state(KALMAN_T *kf);
int kalman_write_out( KALMAN_T *kf, FILE *fp);

int m_x_m10_for_kalman( double M[10][10], double M1[10][10], double M2[10][10] );
int m_x_mt10_for_kalman( double M[10][10], double M1[10][10], double M2[10][10] );
int m_add10_for_kalman
( double M[10][10], double M1[10][10], double M2[10][10] );

int m_x_m9_for_kalman( double M[9][9], double M1[9][9], double M2[9][9] );
int m_x_mt9_for_kalman( double M[9][9], double M1[9][9], double M2[9][9] );
int m_add9_for_kalman
( double M[9][9], double M1[9][9], double M2[9][9] );


int m_x_m6_for_kalman( double M[6][6], double M1[6][6], double M2[6][6] );
int m_x_mt6_for_kalman( double M[6][6], double M1[6][6], double M2[6][6] );
int m_add6_for_kalman
( double M[6][6], double M1[6][6], double M2[6][6] );

int read_meas(MEAS_T *meas, FILE *fp_meas, int init_complete, PARAMS_T *PARAMS);
