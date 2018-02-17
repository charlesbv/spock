#ifndef __PROP_MATH_H__
#define __PROP_MATH_H__
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "propagator.h"
#define SMALL_NUM 1.0e-12
#define DEG2RAD 0.01745329251994
#define RAD2DEG 57.29577951308233
#define NELEMS(x)  (sizeof(x) / sizeof((x)[0]))

#define  CONST_MAT          ( ConstSpiceDouble   (*) [3] )


// Prototypes
int m_print(double m_to_print[3][3],
	    char name[256]);

int v_norm_print(double v_to_print_norm[3],
		 char name[256]);

int v_print( double v_to_print[3],
	     char name[256]);

int v_mag( double *v_mag,
           double v_in[3]);

int v_norm( double u_out[3],
            double v_in[3]);

int v_dot(  double *dot,
            double v1[3],
            double v2[3]);

int v_cross(    double v_cross[3],
                double v1[3],
                double v2[3]);

int v_scale(    double v_out[3],
                double v_in[3],
                double scale);

int v_sub(  double v_out[3],
            double v1[3],
            double v2[3]);

int v_copy( double v_out[3],
            double v_in[3]);

int m_copy( double m_out[3][3],
            double m_in[3][3]);

int m_x_v(  double v_out[3],
            double m_in[3][3],
            double v_in[3]);

int m_x_m( double m_out[3][3],
	   double m_in1[3][3],
	   double m_in2[3][3] );


int m_trans(    double m_out[3][3],
                double m_in[3][3]);

int v_add(  double v_out[3],
            double v1[3],
            double v2[3]);

int compute_T_inrtl_2_ntw( double T_inrtl_2_ntw[3][3],
			   double r[3],
			   double v[3]);

int compute_T_inrtl_2_lvlh( double T_inrtl_2_lvlh[3][3],
                            double r[3],
                            double v[3]);

int compute_T_sc_to_lvlh(double T_sc_to_lvlh[3][3], 
			 double v_angle[3], 
			 int order_rotation[3], 
			 char attitude_profile[256], 
			 double *et,  
			 double r_i2cg_INRTL[3], 
			 double v_i2cg_INRTL[3],
			 int file_is_quaternion,
			 double quaternion[4]);


int compute_T_enu_to_ecef( double T_enu_to_ecef[3][3],
			   double geodetic_latitude, 
			   double longitude, 
			   double flattening);
			   //	   double equatorial_radius );

int etprint( double et_to_print, char str_print[256] );


int m_x_m6( double **m_out,
	   double **m_in1,
	    double **m_in2 );
int m_print6(double **m_to_print,
	     char name[256]);
int v_norm6( double u_out[6],
	     double v_in[6]);
int v_mag6( double *v_mag,
	    double v_in[6]);
int m_x_v6(  double v_out[6],
            double **m_in,
	     double v_in[6]);

/* int compute_T_inrtl_2_sc( double T_inrtl_2_sc[3][3], */
/*                             double r[3], */
/* 			  double v[3], */
/* 			  double v_angle[3], */
/* 			  int order_rotation[3], */
/* 			  char attitude_profile[256], */
/* 			  double *et); */
int q_copy(double q_out[4], double q_in[4]);


#endif
