#ifndef __kalman_h__
#define __kalman_h__

#include "propagator.h"
#include "options.h"
#define N_STATES 6

typedef struct {
  double r_i2cg_INRTL[3];
  double v_i2cg_INRTL[3];
  double et;
  FILE *fp_meas;
} MEAS_T;

typedef struct {  

  double X[6];
  double P[6][6];
  double dX[6];
  double update[6];

  double Phi[6][6];
  double Q[6][6];

    double tau;
    
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

int kalman_cov_prop( KALMAN_T *kf, PARAMS_T *PARAMS );
int kalman_filt( MEAS_T *meas, KALMAN_T *kf, PARAMS_T *params, OPTIONS_T *OPTIONS, GROUND_STATION_T *GROUND_STATION, CONSTELLATION_T  *CONSTELLATION, int iProc, int iDebugLevel, int nProcs );
int kalman_init( MEAS_T *meas, KALMAN_T *kf, CONSTELLATION_T  *CONSTELLATION);
int kalman_prop( MEAS_T *meas, KALMAN_T *kf, PARAMS_T *params, OPTIONS_T *OPTIONS,GROUND_STATION_T *GROUND_STATION, CONSTELLATION_T  *CONSTELLATION, int iProc, int iDebugLevel, int nProcs );
int kalman_compute_resids( MEAS_T *meas, KALMAN_T *kf, int flag, int index );
int kalman_computeA( double A[6][6], SPACECRAFT_T *sc, PARAMS_T *PARAMS , double dt, double tau);
int kalman_computeSTM( double Phi[6][6], SPACECRAFT_T *sc, PARAMS_T *PARAMS, double dt, double tau );
int kalman_computeQ( double Q[6][6], double Phi[N_STATES][N_STATES], double dt, double tau, double sigma, FILE *fp_kalman_init );
int kalman_process_GPS( MEAS_T *meas, KALMAN_T *kf ) ;
int kalman_computeH_gps( double H[6], KALMAN_T *kf, int ii );
int kalman_joseph_update(  KALMAN_T *kf, double H[6], int meas_i);
int kalman_update_state(KALMAN_T *kf );
int kalman_write_out( KALMAN_T *kf, FILE *fp);

int m_x_m6_for_kalman( double M[6][6], double M1[6][6], double M2[6][6] );
int m_x_mt6_for_kalman( double M[6][6], double M1[6][6], double M2[6][6] );
int m_add6_for_kalman
( double M[6][6], double M1[6][6], double M2[6][6] );

int read_meas(MEAS_T *meas, FILE *fp_meas, int init_complete, PARAMS_T *PARAMS);
