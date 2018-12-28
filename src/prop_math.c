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
// These are math functions


#include "prop_math.h"
#include "propagator.h"

int m_print(double m_to_print[3][3],
	    char name[256])
{
  int c, d;
  
  printf("%s:\n", name);
  for (c = 0; c < 3; c++) {
    for (d = 0; d < 3; d++)
      printf("%f\t", m_to_print[c][d]);
    printf("\n");
  }
  

  return 0;
}

int m_print6(double **m_to_print,
	    char name[256])
{
  int c, d;
  
  printf("%s:\n", name);
  for (c = 0; c < 6; c++) {
    for (d = 0; d < 6; d++)
      printf("%e\t", m_to_print[c][d]);
    printf("\n");
  }
  

  return 0;
}

int v_norm_print(double v_to_print_norm[3],
		 char name[256])
{
  double norm = 0.0;
  v_mag(&norm,v_to_print_norm);
  printf("%s %.10e\n",name, norm);

  return 0;
}

int v_print( double v_to_print[3],
	     char name[256])
{
  printf("%s: (%.10e, %.10e, %.10e)\n",name, v_to_print[0], v_to_print[1], v_to_print[2]);
  return 0;
}

// Compute the vector mag
int v_mag( double *v_mag,
           double v_in[3])

{

    v_mag[0] = sqrt( v_in[0]*v_in[0] + v_in[1]*v_in[1] + v_in[2]*v_in[2]);
    return 0;
}


// Compute the vector mag
int v_mag6( double *v_mag,
           double v_in[6])

{

  v_mag[0] = sqrt( v_in[0]*v_in[0] + v_in[1]*v_in[1] + v_in[2]*v_in[2] + v_in[3]*v_in[3] + v_in[4]*v_in[4] + v_in[5]*v_in[5] );
    return 0;
}

// Compute the normal vector
int v_norm( double u_out[3],
            double v_in[3])

{
    double v_in_mag;
    int ii;
    
    v_mag( &v_in_mag, v_in );
    
    for (ii = 0; ii < 3; ii++) {
    
        u_out[ii] = v_in[ii] / v_in_mag;
    
    }
    
    return 0;
}


// Compute the normal vector
int v_norm6( double u_out[6],
            double v_in[6])

{
    double v_in_mag;
    int ii;
    
    v_mag6( &v_in_mag, v_in );
    
    for (ii = 0; ii < 6; ii++) {
    
        u_out[ii] = v_in[ii] / v_in_mag;
    
    }
    
    return 0;
}

// Compute the dot product
int v_dot(  double *dot,
            double v1[3],
            double v2[3])

{
    
    dot[0] = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    return  0;
}

// Compute the cross product
int v_cross(    double v_cross[3],
                double v1[3],
                double v2[3])

{
    
    v_cross[0] = v1[1]*v2[2] - v1[2]*v2[1];
    v_cross[1] = v1[2]*v2[0] - v1[0]*v2[2];
    v_cross[2] = v1[0]*v2[1] - v1[1]*v2[0];
    
    return  0;
}

// Scale a vector
int v_scale(    double v_out[3],
                double v_in[3],
                double scale)

{
    
    v_out[0] = v_in[0] * scale;
    v_out[1] = v_in[1] * scale;
    v_out[2] = v_in[2] * scale;
    
    return  0;
}

// Subtract vectors
int v_sub(  double v_out[3],
            double v1[3],
            double v2[3])

{
    
    v_out[0] = v1[0] - v2[0];
    v_out[1] = v1[1] - v2[1];
    v_out[2] = v1[2] - v2[2];
    
    return  0;
}

// Add vectors
int v_add(  double v_out[3],
            double v1[3],
            double v2[3])

{
    
    v_out[0] = v1[0] + v2[0];
    v_out[1] = v1[1] + v2[1];
    v_out[2] = v1[2] + v2[2];
    
    return  0;
}


// Copy vectors
int v_copy( double v_out[3],
            double v_in[3])

{
    
    v_out[0] = v_in[0];
    v_out[1] = v_in[1];
    v_out[2] = v_in[2];
    
    return  0;
}


// Copy matrices
int m_copy( double m_out[3][3],
            double m_in[3][3])

{
    
  int i,j;
  for (i = 0; i < 3; i++){
    for (j = 0; j < 3; j++){
      m_out[i][j] = m_in[i][j];
    }
  }
    
    return  0;
}



// Matrix x Vector
int m_x_v(  double v_out[3],
            double m_in[3][3],
            double v_in[3])
{
    v_out[0] = m_in[0][0]*v_in[0] + m_in[0][1]*v_in[1] + m_in[0][2]*v_in[2];
    v_out[1] = m_in[1][0]*v_in[0] + m_in[1][1]*v_in[1] + m_in[1][2]*v_in[2];
    v_out[2] = m_in[2][0]*v_in[0] + m_in[2][1]*v_in[1] + m_in[2][2]*v_in[2];
    
    return 0;
}


// Matrix x Vector
int m_x_v6(  double v_out[6],
            double **m_in,
            double v_in[6])
{
  v_out[0] = m_in[0][0]*v_in[0] + m_in[0][1]*v_in[1] + m_in[0][2]*v_in[2] + m_in[0][3]*v_in[3] + m_in[0][4]*v_in[4] + m_in[0][5]*v_in[5];
    v_out[1] = m_in[1][0]*v_in[0] + m_in[1][1]*v_in[1] + m_in[1][2]*v_in[2] + m_in[1][3]*v_in[3] + m_in[1][4]*v_in[4] +	m_in[1][5]*v_in[5];
    v_out[2] = m_in[2][0]*v_in[0] + m_in[2][1]*v_in[1] + m_in[2][2]*v_in[2] + m_in[2][3]*v_in[3] + m_in[2][4]*v_in[4] +	m_in[2][5]*v_in[5];
    v_out[3] = m_in[3][0]*v_in[0] + m_in[3][1]*v_in[1] + m_in[3][2]*v_in[2] + m_in[3][3]*v_in[3] + m_in[3][4]*v_in[4] +	m_in[3][5]*v_in[5];
    v_out[4] = m_in[4][0]*v_in[0] + m_in[4][1]*v_in[1] + m_in[4][2]*v_in[2] + m_in[4][3]*v_in[3] + m_in[4][4]*v_in[4] +	m_in[4][5]*v_in[5];
    v_out[5] = m_in[5][0]*v_in[0] + m_in[5][1]*v_in[1] + m_in[5][2]*v_in[2] + m_in[5][3]*v_in[3] + m_in[5][4]*v_in[4] +	m_in[5][5]*v_in[5];    
    return 0;
}



// Matrix x Matrix (CBV 07/24/2015)

int m_x_m( double m_out[3][3],
	   double m_in1[3][3],
	   double m_in2[3][3] )
{

  int c, d, k;
  double sum = 0;

  for (c = 0; c < 3; c++) {
    for (d = 0; d < 3; d++) {
      for (k = 0; k < 3; k++) {
	sum = sum + m_in1[c][k]*m_in2[k][d];
      }
      m_out[c][d] = sum;
      sum = 0;
    }
  }

  return 0;
}

// Matrix x Matrix (CBV 10/04/2016)

int m_x_m6( double **m_out,
	   double **m_in1,
	   double **m_in2 )
{

  int c, d, k;
  double sum = 0;
  for (c = 0; c < 6; c++) {
    for (d = 0; d < 6; d++) {
      for (k = 0; k < 6; k++) {
	sum = sum + m_in1[c][k]*m_in2[k][d];
      }
      m_out[c][d] = sum;
      sum = 0;
    }
  }
  return 0;
}

// Matrix Transpose
int m_trans(    double m_out[3][3],
                double m_in[3][3])
{
        m_out[0][1] = m_in[1][0];
        m_out[0][2] = m_in[2][0];
    
        m_out[1][0] = m_in[0][1];
        m_out[1][2] = m_in[2][1];
    
        m_out[2][0] = m_in[0][2];
        m_out[2][1] = m_in[1][2];
    
        m_out[0][0] = m_in[0][0];
        m_out[1][1] = m_in[1][1];
        m_out[2][2] = m_in[2][2];

        return 0;
}
// Inertial to NTW (CBV) (source: Vallado3 page 172)

int compute_T_inrtl_2_ntw( double T_inrtl_2_ntw[3][3],
			   double r[3],
			   double v[3])
{

  // Declarations
  double t[3], w[3], n[3]; // t = in-track; w = cross-track; n = vector in the orbital plane normal to the velocity (directed outwards from the Earth)
  int col;

  // Algorithm
  /* In-track (in the direction of the velocity even if the eccentricity is not 0 (contrarily to the along-track direction in the LVLH frame))*/
  v_norm(t,v); 

  /* Cross-track */
  v_cross(w,r,v);
  v_norm(w,w);

  /* In-track cross cross-track */
  v_cross(n,t,w);

  for (col = 0; col < 3; col++){
    T_inrtl_2_ntw[0][col] = t[col]; // in-track
    T_inrtl_2_ntw[1][col] = w[col]; // cross-track
    T_inrtl_2_ntw[2][col] = n[col]; // in-track cross cross-track
  }

    return 0;
}


/* int compute_T_inrtl_2_sc( double T_inrtl_2_sc[3][3], */
/*                             double r[3], */
/* 			  double v[3], */
/* 			  double v_angle[3], */
/* 			  int order_rotation[3], */
/* 			  char attitude_profile[256], */
/* 			  double *et) */
			  
/* { */


/*   double T_inrtl_2_lvlh[3][3]; */
/*     double T_sc_to_lvlh[3][3]; */
/*   double T_lvlh_2_sc[3][3]; */


/*     compute_T_inrtl_2_lvlh(T_inrtl_2_lvlh, r, v); */
/*     compute_T_sc_to_lvlh(T_sc_to_lvlh,  v_angle,  order_rotation,  attitude_profile, et,  r, v); */
/*     m_trans(T_lvlh_2_sc, T_sc_to_lvlh); */
      
/*   m_x_m(T_inrtl_2_sc, T_inrtl_2_lvlh, T_lvlh_2_sc); */


/*   return 0; */
/* } */

int q_copy(double q_out[4], double q_in[4]){
  
  q_out[0] = q_in[0];
  q_out[1] = q_in[1];
  q_out[2] = q_in[2];
  q_out[3] = q_in[3];

  return 0;
}

//  Inertial to LVLH
int compute_T_inrtl_2_lvlh( double T_inrtl_2_lvlh[3][3],
                            double r[3],
                            double v[3])
{
    // Declarations
    double ir[3];
    double iv[3];
    double in[3];
    
    // Algorithm
 
    /* Radial (1/2) */
    v_norm(ir, r);
    v_scale(ir, ir, -1.0);    

    /* Cross-track */
    v_cross(in, ir, v);
    v_norm(in, in);
    v_scale(in, in, -1.0);
    
    /* Along-track (different from the velocity direction if the eccentricity is not 0)*/
    v_cross(iv, in, ir);
    v_norm(iv, iv);
    v_scale(iv, iv, -1.0);
    
    /* Radial (2/2) */
    v_scale(ir, ir, -1.0);   // this lign has been added by CBV on 07-29-2015 so that the basis (iv, in, ir) is oriented in the right hand direction

    /* Along-track */
    T_inrtl_2_lvlh[0][0] = iv[0]; 
    T_inrtl_2_lvlh[0][1] = iv[1];
    T_inrtl_2_lvlh[0][2] = iv[2];
    
    /* Cross-track */
    T_inrtl_2_lvlh[1][0] = in[0];
    T_inrtl_2_lvlh[1][1] = in[1];
    T_inrtl_2_lvlh[1][2] = in[2];
    
    /* Radial */
    T_inrtl_2_lvlh[2][0] = ir[0]; 
    T_inrtl_2_lvlh[2][1] = ir[1]; 
    T_inrtl_2_lvlh[2][2] = ir[2]; 
    

    return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           Compute_T_Inrtl_2_Lvlh
//  Purpose:        Convert a vector in the SC reference system into the LVLH reference system 
//  Assumptions:    None
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 08/01/2015    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                                                                                                                                                     
int compute_T_sc_to_lvlh(double T_sc_to_lvlh[3][3], double v_angle[3], int order_rotation[3], char attitude_profile[256], double *et,  double r_i2cg_INRTL[3], double v_i2cg_INRTL[3]	,int file_is_quaternion,  double quaternion[4]){

  /* Declarations */
  int row;
  double x[6]; 
  double lt;
  double r_earth2sun_J2000[3]; 
  double r_cg2sun_J2000[3];
  double r_cg2sun_J2000_normalized[3]; 
  double e_z_body_in_inrtl[3];
  double T_inrtl_2_lvlh[3][3];
  double e_z_body_in_lvlh[3];
  /* double random_vect_not_colinear_to_z[3];  */
  /* double e_y_body_in_lvlh[3];  */
  /* double e_x_body_in_lvlh[3];  */
  double e_z_body_in_lvlh_normalized[3]; 
  /* double e_y_body_in_lvlh_normalized[3]; */
  /* double e_x_body_in_lvlh_normalized[3]; */
  int i;
  double theta, phi, psi;
  double pitch_mat[3][3], roll_mat[3][3], yaw_mat[3][3];
  double first_mat[3][3], second_mat[3][3], third_mat[3][3];
  double T_sc_to_lvlh_temp[3][3];
  //  double T_sc_intermediary_to_lvlh[3][3];
  double v_angle_rad[3];

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  ///////////////////// QUATERNION /////////////////////////////////////////////

  double *p=malloc(sizeof(double)), *r=malloc(sizeof(double)), *y=malloc(sizeof(double));
  if (file_is_quaternion == 1){ // if attitude is set using quaternions     
    q2m_c(quaternion,T_sc_to_lvlh);
/*             m_print(T_sc_to_lvlh, "quat T_sc_to_lvlh"); */
/* 	    m2eul_c ( T_sc_to_lvlh, 2, 1, 3, p, r, y ); */
/* 	    	    printf("%f %f %f\n",p[0]*180./M_PI, r[0]*180./M_PI, y[0]*180./M_PI ); */
	      
  }

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  ///////////////////// NOT QUATERNION /////////////////////////////////////////////
  else{ //no quaternion


  //////////////////////////////////////////////////////////////////////////////
  ///////////////////// ANY ATTITUDE EXCEPT SUN POINTED ////////////////////////
  if (  strcmp( attitude_profile, "sun_pointed"  ) != 0 )  { // if the attitude profile is anything except sun pointed
  for (i = 0; i < 3; i++){
    v_angle_rad[i] = v_angle[i] *  DEG2RAD;
  }
  
  theta = v_angle_rad[0];// pitch
  phi = v_angle_rad[1]; // roll
  psi = v_angle_rad[2]; // yaw

  // Pitch rotation
  pitch_mat[0][0] = cos(theta);   pitch_mat[0][1] = 0;   pitch_mat[0][2] = sin(theta);
  pitch_mat[1][0] = 0;   pitch_mat[1][1] = 1;   pitch_mat[1][2] = 0;
  pitch_mat[2][0] = -sin(theta);   pitch_mat[2][1] = 0;   pitch_mat[2][2] = cos(theta);

  // Roll rotation
  roll_mat[0][0] = 1;   roll_mat[0][1] = 0;   roll_mat[0][2] = 0;
  roll_mat[1][0] = 0;   roll_mat[1][1] = cos(phi);   roll_mat[1][2] = -sin(phi);
  roll_mat[2][0] = 0;   roll_mat[2][1] = sin(phi);   roll_mat[2][2] = cos(phi);

  // Yaw rotation
  yaw_mat[0][0] = cos(psi);   yaw_mat[0][1] = -sin(psi);   yaw_mat[0][2] = 0;
  yaw_mat[1][0] = sin(psi);   yaw_mat[1][1] = cos(psi);   yaw_mat[1][2] = 0;
  yaw_mat[2][0] = 0;   yaw_mat[2][1] = 0;   yaw_mat[2][2] = 1;

  // Order of rotation
  // // First matrix
/*     etprint(et[0], ""); */
/*     printf("%d %d %d\n", order_rotation[0], order_rotation[1], order_rotation[2]); */

//  printf("%d %d %d\n", order_rotation[0], order_rotation[1], order_rotation[2]);

  if (order_rotation[0] == 1)
    m_copy(first_mat, pitch_mat);
  else if (order_rotation[1] == 1)
    m_copy(first_mat, roll_mat);
  else if (order_rotation[2] == 1)
    m_copy(first_mat, yaw_mat);
  else{
    printf("***! (compute_T_sc_to_lvlh) You did not correctly choose the first rotation. The program will stop. !***\n"); MPI_Finalize(); 
    exit(0);
  }
  // // Second matrix
  if (order_rotation[0] == 2)
    m_copy(second_mat, pitch_mat);
  else if (order_rotation[1] == 2)
    m_copy(second_mat, roll_mat);
  else if (order_rotation[2] == 2)
    m_copy(second_mat, yaw_mat);
  else{
    printf("***! (compute_T_sc_to_lvlh) You did not correctly choose the second rotation. The program will stop. !***\n"); MPI_Finalize(); 
    exit(0);
  }
  // // Third matrix
  if (order_rotation[0] == 3)
    m_copy(third_mat, pitch_mat);
  else if (order_rotation[1] == 3)
    m_copy(third_mat, roll_mat);
  else if (order_rotation[2] == 3)
    m_copy(third_mat, yaw_mat);
  else{
    printf("***! (compute_T_sc_to_lvlh) You did not correctly choose the third rotation. The program will stop. !***\n"); MPI_Finalize(); 
    exit(0);
  }
  // // Final matrix: T_sc_to_lvlh
  m_x_m(T_sc_to_lvlh_temp, second_mat, first_mat);
  m_x_m(T_sc_to_lvlh, third_mat, T_sc_to_lvlh_temp);
  //            m_print(T_sc_to_lvlh, "euler T_sc_to_lvlh");
  }

  //////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// SUN POINTED ///////////////////////////////////
  if (  strcmp( attitude_profile, "sun_pointed"  ) == 0 )  {
    /* Sun-pointed */
    /* Express e_z_body in Inertial frame (J2000) */
    spkez_c(10, et[0], "J2000", "NONE", 399, x, &lt); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration.
    r_earth2sun_J2000[0] = x[0];
    r_earth2sun_J2000[1] = x[1];
    r_earth2sun_J2000[2] = x[2];


    v_sub(r_cg2sun_J2000, r_earth2sun_J2000, r_i2cg_INRTL);
    v_norm(r_cg2sun_J2000_normalized, r_cg2sun_J2000);
    v_copy(e_z_body_in_inrtl, r_cg2sun_J2000_normalized);

    /* Convert e_z_body in LVLH frame */
    compute_T_inrtl_2_lvlh(T_inrtl_2_lvlh, r_i2cg_INRTL, v_i2cg_INRTL);
    m_x_v(e_z_body_in_lvlh, T_inrtl_2_lvlh, e_z_body_in_inrtl);
    v_norm(e_z_body_in_lvlh_normalized,e_z_body_in_lvlh); // useless cause e_z_body_in_lvlh should already be normalized


    /************************************** BLOCK A: UNCOMMENT BELOW TO SET X_BODY PERPENDICULAR TO NADIR **********************************/
    // // IMPORTANT: IF YOU UNCOMMENT THIS BLOCK THEN COMMENT BLOCK B
    // Sun pointed gives the direction of z_body. The sc can still yaw around this vector sc to Sun. By default, we assume that the yaw angle is so that the x_body is perpendicular to nadir. In other words, x_body is the cross product of nadir with sat_to_sun vector (both expressed in lvlh coordinates)
    double lvlh_z_in_lvlh[3], new_e_x_body_in_lvlh[3], new_e_y_body_in_lvlh_normalized[3], new_e_x_body_in_lvlh_normalized[3];
    lvlh_z_in_lvlh[0] = 0; lvlh_z_in_lvlh[1] = 0; lvlh_z_in_lvlh[2] = 1; // lvlh_z is in the nadir (oriented away from the Earth)
    v_cross(new_e_x_body_in_lvlh, lvlh_z_in_lvlh, e_z_body_in_lvlh_normalized);
    v_norm(new_e_x_body_in_lvlh_normalized, new_e_x_body_in_lvlh);
    
    v_cross(new_e_y_body_in_lvlh_normalized, e_z_body_in_lvlh_normalized, new_e_x_body_in_lvlh_normalized);


    // // !!!!!!! BELOW IS TO ROTATE X_BODY BY 10 DEGREES IN THE PLANE PERPENDICULAR TO SAT_TO_SUN VECTOR. THEN X_BODY IS 10 DEGREES FROM BEING PERPENDICULAR TO NADIR (BUT STILL PERPENDICULAR TO Z_BODY (WHICH IS EQUAL TO SAT_TO_SUN VECTOR)). PARDON THE LONG NAMES
    // // !!!!!! THEREFORE, COMMENT THE BLOCK BELOW (UNLESS YOU KNOW WHAT YOU ARE DOING)
    double new_e_x_body_in_lvlh_normalized_rotated_10_degrees[3], new_e_x_body_in_lvlh_normalized_rotated_10_degrees_component_along_new_e_x_body_in_lvlh_normalized[3], new_e_x_body_in_lvlh_normalized_rotated_10_degrees_component_along_new_e_y_body_in_lvlh_normalized[3],new_e_x_body_in_lvlh_normalized_rotated_10_degrees_component_along_e_z_body_in_lvlh_normalized[3], new_e_x_body_in_lvlh_normalized_rotated_10_degrees_temp[3], new_e_x_body_in_lvlh_normalized_rotated_10_degrees_normalized[3];
    v_scale(new_e_x_body_in_lvlh_normalized_rotated_10_degrees_component_along_new_e_x_body_in_lvlh_normalized, new_e_x_body_in_lvlh_normalized, cos( 10. * DEG2RAD ));
    v_scale( new_e_x_body_in_lvlh_normalized_rotated_10_degrees_component_along_new_e_y_body_in_lvlh_normalized, new_e_y_body_in_lvlh_normalized, sin( 10. * DEG2RAD ));
    new_e_x_body_in_lvlh_normalized_rotated_10_degrees_component_along_e_z_body_in_lvlh_normalized[0] = 0;  new_e_x_body_in_lvlh_normalized_rotated_10_degrees_component_along_e_z_body_in_lvlh_normalized[1] = 0;  new_e_x_body_in_lvlh_normalized_rotated_10_degrees_component_along_e_z_body_in_lvlh_normalized[2] = 0;

    v_add( new_e_x_body_in_lvlh_normalized_rotated_10_degrees_temp, new_e_x_body_in_lvlh_normalized_rotated_10_degrees_component_along_new_e_x_body_in_lvlh_normalized, new_e_x_body_in_lvlh_normalized_rotated_10_degrees_component_along_new_e_y_body_in_lvlh_normalized );
    v_add( new_e_x_body_in_lvlh_normalized_rotated_10_degrees, new_e_x_body_in_lvlh_normalized_rotated_10_degrees_temp, new_e_x_body_in_lvlh_normalized_rotated_10_degrees_component_along_e_z_body_in_lvlh_normalized );
    v_norm( new_e_x_body_in_lvlh_normalized_rotated_10_degrees_normalized, new_e_x_body_in_lvlh_normalized_rotated_10_degrees);
    new_e_x_body_in_lvlh_normalized[0] = 0;  new_e_x_body_in_lvlh_normalized[1] = 0;  new_e_x_body_in_lvlh_normalized[2] = 0;
    v_copy( new_e_x_body_in_lvlh_normalized, new_e_x_body_in_lvlh_normalized_rotated_10_degrees_normalized );
    new_e_y_body_in_lvlh_normalized[0] = 0;  new_e_y_body_in_lvlh_normalized[1] = 0;  new_e_y_body_in_lvlh_normalized[2] = 0;
    v_cross( new_e_y_body_in_lvlh_normalized, e_z_body_in_lvlh_normalized, new_e_x_body_in_lvlh_normalized );

    // // !!!!!! END OF 'THEREFORE, COMMENT THE BLOCK BELOW (UNLESS YOU KNOW WHAT YOU ARE DOING)'

    /************************************** END OF 'BLOCK A: UNCOMMENT BELOW TO SET X_BODY PERPENDICULAR TO NADIR' **********************************/

    /************************************** BLOCK B: UNCOMMENT BELOW TO SET X_BODY IN THE DIRECTION OF LVLH_X **********************************/
    /* /\* // // IMPORTANT: IF YOU UNCOMMENT THIS BLOCK THEN COMMENT BLOCK A *\/ */
    /* /\* // Sun pointed gives the direction of z_body. The sc can still roll around this vector sc to Sun. By default, we assume that the roll angle is so that the x_body points towards lvlh_x. In other words, the projection of lvlv_x on the plane perpendicular to the Sun give x_body *\/ */
    /* double v_int[3], lvlh_x_in_lvlh[3], new_e_x_body_in_lvlh[3], new_e_y_body_in_lvlh_normalized[3], v_int_norm[3], new_e_x_body_in_lvlh_normalized[3]; */
    
    /* lvlh_x_in_lvlh[0] = 1; lvlh_x_in_lvlh[1] = 0; lvlh_x_in_lvlh[2] = 0; */
    /* v_cross(v_int, e_z_body_in_lvlh_normalized, lvlh_x_in_lvlh); */
    /* v_norm(v_int_norm, v_int); */
    /* v_cross(new_e_x_body_in_lvlh, v_int_norm, e_z_body_in_lvlh_normalized); */
    /* v_norm(new_e_x_body_in_lvlh_normalized, new_e_x_body_in_lvlh); */
    /* v_cross(new_e_y_body_in_lvlh_normalized, e_z_body_in_lvlh_normalized, new_e_x_body_in_lvlh_normalized); */

    /* /\*   /\\* /\\\* // Another way: *\\\/ *\\/ *\/ */
    /* /\*   /\\* double u_dot_n, u_dot_n_scale[3]; *\\/ *\/ */
    /* /\*   /\\* v_dot(&u_dot_n, lvlh_x_in_lvlh, e_z_body_in_lvlh_normalized); *\\/ *\/ */
    /* /\*   /\\* v_scale(u_dot_n_scale, e_z_body_in_lvlh_normalized,u_dot_n); *\\/ *\/ */
    /* /\*   /\\* v_sub(new_e_x_body_in_lvlh_normalized, lvlh_x_in_lvlh, u_dot_n_scale); *\\/ *\/ */
    /* /\*   /\\* v_norm(new_e_x_body_in_lvlh_normalized, new_e_x_body_in_lvlh_normalized); *\\/ *\/ */
    /* /\*   /\\* v_print(new_e_x_body_in_lvlh_normalized, "new_e_x_body_in_lvlh_normalized"); *\\/ *\/ */
    /************************************** END OF 'BLOCK B: UNCOMMENT BELOW TO SET X_BODY IN THE DIRECTION OF LVLH_X' **********************************/

    /* // We now have the new body frame, expressed in the LVLH frame: the z axis points towards the Sun, and the x axis is the projection of LVLH_X on the plane perdicular to the satellite-to-Sun vector. */




    for (row = 0; row<3; row++){
      T_sc_to_lvlh[row][0] = new_e_x_body_in_lvlh_normalized[row];
      T_sc_to_lvlh[row][1] = new_e_y_body_in_lvlh_normalized[row];
      T_sc_to_lvlh[row][2] = e_z_body_in_lvlh_normalized[row];
    }
    

  }
  } // not quaternion


  return 0;

  // // OLD STUFF FOR SUN POINTED:
    /* /\* Now the z axis of the body is fixed, set the rotation around this z axis to complete the (x, y, z)body *\/ */
    /* /\* For now, we just take any "random" vector orthogonal to e_z_body to be e_y_body. This means that the x and y axes have a "random" orientation in the plane orthonal to the satellite-Sun direction. *\/ */
    /* /\* Creation of a vector not colinear to e_z_body *\/ */
    /* if (e_z_body_in_lvlh[0] != 0){ */
    /*   random_vect_not_colinear_to_z[1] = 1.0; */
    /*   random_vect_not_colinear_to_z[0] = 0.0; */
    /*   random_vect_not_colinear_to_z[2] = 0.0; */
    /* } */
    /* else if (e_z_body_in_lvlh[1] != 0){ */
    /*   random_vect_not_colinear_to_z[0] = 1.0; */
    /*   random_vect_not_colinear_to_z[1] = 0.0; */
    /*   random_vect_not_colinear_to_z[2] = 0.0; */
    /* } */
    /* else if (e_z_body_in_lvlh[2] != 0){ */
    /*   random_vect_not_colinear_to_z[1] = 1.0; */
    /*   random_vect_not_colinear_to_z[2] = 0.0; */
    /*   random_vect_not_colinear_to_z[0] = 0.0; */
    /* } */
    /* else{ */
    /*   printf("It seems that the vector satellite-to-Sun is 0. The program will stop.\n"); */
    /*   exit(0); */
    /* } */

    /* /\* e_y_body is the cross product between e_z_body and this "random" non-colinear vector to e_z_body. Consequently, e_y_body is a "random" vector orthogonal to e_z_body *\/ */
    /* v_cross(e_y_body_in_lvlh, e_z_body_in_lvlh, random_vect_not_colinear_to_z); */
    /* v_cross(e_x_body_in_lvlh, e_y_body_in_lvlh, e_z_body_in_lvlh); */

    /* /\* Normalization of the SC basis *\/ */
    /* v_norm(e_x_body_in_lvlh_normalized,e_x_body_in_lvlh); */
    /* v_norm(e_y_body_in_lvlh_normalized,e_y_body_in_lvlh); */
    /* v_norm(e_z_body_in_lvlh_normalized,e_z_body_in_lvlh); // useless cause e_z_body_in_lvlh should already be normalized */

    /* /\* T_sc_intermediary_to_lvlh is the rotation from body to lvlh *\/ */
    /* for (row = 0; row<3; row++){ */
    /*   T_sc_intermediary_to_lvlh[row][0] = e_x_body_in_lvlh_normalized[row]; */
    /*   T_sc_intermediary_to_lvlh[row][1] = e_y_body_in_lvlh_normalized[row]; */
    /*   T_sc_intermediary_to_lvlh[row][2] = e_z_body_in_lvlh_normalized[row]; */
    /* } */

}




/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           compute_T_enu_to_ecef
//  Purpose:        Compute the rotation matrix from ENU to ECEF (taking into account the oblateness of the Earth)
//  Assumptions:    None
//  References      Vallado third edition section 3.2
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 06/02/2016    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                                                                                                                            
int compute_T_enu_to_ecef( double T_enu_to_ecef[3][3], // out: rotation matrix ENU to ECEF
			   double geodetic_latitude, // in: geodetic latitude of the place at the surface of the Earth
			   double longitude, // in: longitude of the place at the surface of the Earth
			   double flattening) /* IN:     flattening parameter     */
			   //			   double equatorial_radius ) /* IN: equatorial radius  (km) */

{
  /* Declarations */
  double T_enu_to_east_north_center_earth_to_local[3][3], T_east_north_center_earth_to_local_to_ecef[3][3];
  double geocentric_latitude;    
  /* Algorithm */

  // Compute the geocentric latitude (equation 3.11 of Vallado)
   geocentric_latitude = atan( ( 1 - flattening ) * ( 1 - flattening ) * tan( geodetic_latitude ) );

  // Compute the rotation matrix from ENU to the same frame as ENU except that we rotate the ENU of an angle (geodetic_latitude - geocentric_latitude) around the East axis 
  double rot_angle;
   rot_angle = geodetic_latitude - geocentric_latitude; 
   T_enu_to_east_north_center_earth_to_local[0][0] = 1; T_enu_to_east_north_center_earth_to_local[0][1] = 0; T_enu_to_east_north_center_earth_to_local[0][2] = 0;   
   T_enu_to_east_north_center_earth_to_local[1][0] = 0; T_enu_to_east_north_center_earth_to_local[1][1] = cos(rot_angle); T_enu_to_east_north_center_earth_to_local[1][2] = sin(rot_angle); 
   T_enu_to_east_north_center_earth_to_local[2][0] = 0; T_enu_to_east_north_center_earth_to_local[2][1] = -sin(rot_angle); T_enu_to_east_north_center_earth_to_local[2][2] = cos(rot_angle);


   // Compute the rotation matrix from the local east_north_center_earth_to_local of the place at the surface of the earth to ECEF. This is different from ENU because of the oblateness of the Earth. This frame uses the geocentric latitude (because it goes through the center of the Earth, like the ECEF). ENU uses the geodetic latitude. Therefore, to convert from ENU to ECEF, first use compute_T_enu_to_east_north_center_earth_to_local then compute_t_east_north_center_earth_to_local_to_ecef
   T_east_north_center_earth_to_local_to_ecef[0][0] = -sin(longitude); T_east_north_center_earth_to_local_to_ecef[0][1] = -sin(geocentric_latitude)*cos(longitude); T_east_north_center_earth_to_local_to_ecef[0][2] = cos(geocentric_latitude)*cos(longitude);   
   T_east_north_center_earth_to_local_to_ecef[1][0] = cos(longitude); T_east_north_center_earth_to_local_to_ecef[1][1] = -sin(geocentric_latitude)*sin(longitude); T_east_north_center_earth_to_local_to_ecef[1][2] = cos(geocentric_latitude)*sin(longitude); 
   T_east_north_center_earth_to_local_to_ecef[2][0] = 0; T_east_north_center_earth_to_local_to_ecef[2][1] = cos(geocentric_latitude); T_east_north_center_earth_to_local_to_ecef[2][2] = sin(geocentric_latitude);

   // Rotation matrix from ENU to ECEF: first rotate from ENU to east_north_center_earth_to_local. Then rotation from east_north_center_earth_to_local to ECEF
  m_x_m(T_enu_to_ecef, T_east_north_center_earth_to_local_to_ecef, T_enu_to_east_north_center_earth_to_local);


  return 0;
}


int etprint( double et_to_print, char str_print[256] ){
  char str[256];
  et2utc_c(et_to_print, "ISOC", 6, 255, str);
  printf("%s: %s\n", str_print,str);
  
  return 0;
}




/* ///////////////////////////////////////////////////////////////////////////////////////// */
/* // */
/* //  Name:           Compute_T_Inrtl_2_Lvlh */
/* //  Purpose:        Convert a vector in the SC reference system into the LVLH reference system given the components of the along-track axis (lvlh_alongtrack_in_body_cartesian) and the cross-track axis (lvlh_crosstrack_in_body_cartesian) of the LVLH frame in the SC reference system under the cartesian representation */
/* //  Assumptions:    None. */
/* //  References       */
/* // */
/* //  Change Log: */
/* //      |   Developer   |       Date    |   SCR     |   Notes */
/* //      | --------------|---------------|-----------|------------------------------- */
/* //      | C. Bussy-Virat| 07/24/2015    |   ---     | Initial Implementation */
/* // */
/* ///////////////////////////////////////////////////////////////////////////////////////// */

/* int compute_T_sc_to_lvlh( double T_sc_to_lvlh[3][3], */
/* 			  double lvlh_alongtrack_in_body_cartesian[3], // Along-track vector in the SC reference system */
/* 			  double lvlh_crosstrack_in_body_cartesian[3], // Cross-track vector in the SC reference system */
/* 			  double *et, */
/* 			  double r_i2cg_INRTL[3], */
/* 			  double v_i2cg_INRTL[3], */
/* 			  INTEGRATOR_T    *INTEGRATOR)  */
/* { */

/*   /\* Variable declarations *\/ */

/*   int column, row; */
/*   double lvlh_radial_in_body_cartesian[3], lvlh_radial_in_body_cartesian_normalized[3]; // Radial vector in the SC reference system */
/*   double lvlh_alongtrack_in_body_cartesian_normalized[3],lvlh_crosstrack_in_body_cartesian_normalized[3]; */
/*   double lvlh_alongtrack_in_body_cartesian_DOT_lvlh_crosstrack_in_body_cartesian = 0; */
/*   double x[6]; */
/*   double lt; */
/*   double r_earth2sun_J2000[3]; */
/*   double r_cg2sun_J2000[3]; */
/*   double r_cg2sun_J2000_normalized[3]; */
/*   double e_z_body_in_inrtl[3]; */
/*   double e_z_body_in_lvlh[3]; */
/*   double e_y_body_in_lvlh[3]; */
/*   double e_x_body_in_lvlh[3]; */
/*   double e_z_body_in_lvlh_normalized[3]; */
/*   double e_y_body_in_lvlh_normalized[3]; */
/*   double e_x_body_in_lvlh_normalized[3]; */
/*   double T_inrtl_2_lvlh[3][3]; */
/*   double random_vect_not_colinear_to_z[3]; */

 
/*   if (  strcmp( INTEGRATOR->attitude.attitude_profile, "sun_pointed"  ) == 0 )  { */
/*     /\* Sun-pointed *\/ */

/*     /\* Express e_z_body in Inertial frame (J2000) *\/ */
/*     spkez_c(10, et[0], "J2000", "NONE", 399, x, &lt); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration.  */
/*     r_earth2sun_J2000[0] = x[0]; */
/*     r_earth2sun_J2000[1] = x[1]; */
/*     r_earth2sun_J2000[2] = x[2]; */
  
/*     v_sub(r_cg2sun_J2000, r_earth2sun_J2000, r_i2cg_INRTL); */
/*     v_norm(r_cg2sun_J2000_normalized, r_cg2sun_J2000); */
/*     v_copy(e_z_body_in_inrtl, r_cg2sun_J2000_normalized); */
/*     v_scale(e_z_body_in_inrtl,e_z_body_in_inrtl,-1); */

/*     /\* Convert e_z_body in LVLH frame *\/ */
/*     compute_T_inrtl_2_lvlh(T_inrtl_2_lvlh, r_i2cg_INRTL, v_i2cg_INRTL); */
/*     m_x_v(e_z_body_in_lvlh, T_inrtl_2_lvlh, e_z_body_in_inrtl); */

/*     /\* Now the z axis of the body is fixed, set the rotation around this z axis to complete the (x, y, z)body *\/ */
/*     /\* For now, we just take any "random" vector orthogonal to e_z_body to be e_y_body. This means that the x and y axes have a "random" orientation in the plane orthonal to the satellite-Sun direction. Then we add the possibility of rotating around the z axis with a rate called rate_rot_z_body (in radian / s) *\/ */

/*     /\* Creation of a vector not colinear to e_z_body *\/ */
/*     if (e_z_body_in_lvlh[0] != 0){ */
/*       random_vect_not_colinear_to_z[1] = 1.0;  */
/*       random_vect_not_colinear_to_z[0] = 0.0;  */
/*       random_vect_not_colinear_to_z[2] = 0.0;  */
/*     } */
/*     else if (e_z_body_in_lvlh[1] != 0){ */
/*       random_vect_not_colinear_to_z[0] = 1.0;  */
/*       random_vect_not_colinear_to_z[1] = 0.0;  */
/*       random_vect_not_colinear_to_z[2] = 0.0;  */
/*     } */
/*     else if (e_z_body_in_lvlh[2] != 0){ */
/*       random_vect_not_colinear_to_z[1] = 1.0;  */
/*       random_vect_not_colinear_to_z[2] = 0.0;  */
/*       random_vect_not_colinear_to_z[0] = 0.0;  */
/*     } */
/*     else{ */
/*       printf("It seems that the vector satellite-to-Sun is 0. The program will stop.\n"); */
/*       exit(0); */
/*     } */

/*     /\* e_y_body is the cross product between e_z_body and this "random" non-colinear vector to e_z_body. Consequently, e_y_body is a "random" vector orthogonal to e_z_body *\/ */
/*     v_cross(e_y_body_in_lvlh, e_z_body_in_lvlh, random_vect_not_colinear_to_z); */
/*     v_cross(e_x_body_in_lvlh, e_y_body_in_lvlh, e_z_body_in_lvlh); */

/*     /\* Normalization of the SC basis *\/ */
/*     v_norm(e_x_body_in_lvlh_normalized,e_x_body_in_lvlh); */
/*     v_norm(e_y_body_in_lvlh_normalized,e_y_body_in_lvlh); */
/*     v_norm(e_z_body_in_lvlh_normalized,e_z_body_in_lvlh); // useless cause e_z_body_in_lvlh should already be normalize  */

/*     /\* Now we can allow the satellite to rotate around the z axis with a constant angular rate rate_rot_z_body *\/ */
/*     //////////////////// TO DO //////////////// */

/*     /\* T_sc_to_lvlh is the rotation from body to lvlh *\/ */
/*     for (row = 0; row<3; row++){ */
/*       T_sc_to_lvlh[row][0] = e_x_body_in_lvlh_normalized[row]; */
/*       T_sc_to_lvlh[row][1] = e_y_body_in_lvlh_normalized[row]; */
/*       T_sc_to_lvlh[row][2] = e_z_body_in_lvlh_normalized[row]; */
/*     } */

/* /\*     printf("\n T_sc_to_lvlh for sun_pointed:\n"); *\/ */
/* /\*     int c, d; *\/ */
/* /\*     for (c = 0; c < 3; c++) { *\/ */
/* /\*       for (d = 0; d < 3; d++) *\/ */
/* /\*         printf("%f\t", T_sc_to_lvlh[c][d]); *\/ */
/* /\*        printf("\n"); *\/ */
/* /\*     } *\/ */

/*   } */

/*   /\* any attitude profile but sun_pointed *\/ */
/*   else{ */
/*     /\* v_print(lvlh_alongtrack_in_body_cartesian,"along"); *\/ */
/* /\*     v_print(lvlh_crosstrack_in_body_cartesian,"cross"); *\/ */
    
/*     v_dot(&lvlh_alongtrack_in_body_cartesian_DOT_lvlh_crosstrack_in_body_cartesian,lvlh_alongtrack_in_body_cartesian,lvlh_crosstrack_in_body_cartesian); */

/*     if (floor(fabs(lvlh_alongtrack_in_body_cartesian_DOT_lvlh_crosstrack_in_body_cartesian)*1000000.0) != 0){ */
/*       printf("\nError: the along-track and cross-track directions (LVLH) are not perpendicular to each other. The program will stop.\n"); */
/*       exit(0); */
/*     } */

/*     /\* Radial vector (= along-track cross cross-track) *\/ */
/*     v_cross(lvlh_radial_in_body_cartesian, lvlh_alongtrack_in_body_cartesian, lvlh_crosstrack_in_body_cartesian);  */
  
/*     /\* Normalization of the LVLH basis (in case the user enters not normalized vectors in the input file) *\/ */
/*     v_norm(lvlh_radial_in_body_cartesian_normalized,lvlh_radial_in_body_cartesian); */
/*     v_norm(lvlh_alongtrack_in_body_cartesian_normalized,lvlh_alongtrack_in_body_cartesian); */
/*     v_norm(lvlh_crosstrack_in_body_cartesian_normalized,lvlh_crosstrack_in_body_cartesian); */

/*     for (column = 0; column<3; column++){ */
/*       T_sc_to_lvlh[0][column] = lvlh_alongtrack_in_body_cartesian_normalized[column]; */
/*       T_sc_to_lvlh[1][column] = lvlh_crosstrack_in_body_cartesian_normalized[column]; */
/*       T_sc_to_lvlh[2][column] = lvlh_radial_in_body_cartesian_normalized[column]; */
/*     } */
/*   } */
/*   /\*   printf("\n T_sc_to_lvlh:\n"); *\/ */
/*   /\*   int c, d; *\/ */
/*   /\*   for (c = 0; c < 3; c++) { *\/ */
/*   /\*     for (d = 0; d < 3; d++) *\/ */
/*   /\*       printf("%f\t", T_sc_to_lvlh[c][d]); *\/ */
/*   /\*      printf("\n"); *\/ */
/*   /\*   } *\/ */
  

  
/*   /\*      printf("\n                  SC               ->               LVLH CBV\n"); *\/ */
/*   /\*      int c; *\/ */
/*   /\*      for (c = 0; c < 3; c++) { *\/ */
/*   /\*        printf("%30.20f    %30.20f\n",v_sc[c],v_lvlh[c]); *\/ */
/*   /\*      } *\/ */
/*   /\*      printf("\n"); *\/ */



/*   return 0; */

/* } /\* ---------- end of compute_T_sc_to_lvlh ----------*\/ */
