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

int m_print6(double m_to_print[6][6],
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
int m_print6_temp(double m_to_print[6][6],
	    char name[256])
{
  int c, d;
  
  printf("%s = [ ...\n", name);
  for (c = 0; c < 6; c++) {
    if (c != 5){
    for (d = 0; d < 6; d++){
      if (d != 5){
      printf("%e,\t", m_to_print[c][d]);
      }
      else{
      printf("%e; ...\n", m_to_print[c][d]);
      }
    }
    }
    else{
    for (d = 0; d < 6; d++){
      if (d != 5){
      printf("%e,\t", m_to_print[c][d]);
      }
      else{
      printf("%e];\n", m_to_print[c][d]);
      }
    }
    }
  }
  

  return 0;
}






int m_print9_temp(double m_to_print[9][9],
	    char name[256])
{
  int c, d;
  
  printf("%s = [ ...\n", name);
  for (c = 0; c < 9; c++) {
    if (c != 8){
    for (d = 0; d < 9; d++){
      if (d != 8){
      printf("%e,\t", m_to_print[c][d]);
      }
      else{
      printf("%e; ...\n", m_to_print[c][d]);
      }
    }
    }
    else{
    for (d = 0; d < 9; d++){
      if (d != 8){
      printf("%e,\t", m_to_print[c][d]);
      }
      else{
      printf("%e];\n", m_to_print[c][d]);
      }
    }
    }
  }
  

  return 0;
}


int m_print8_temp(double m_to_print[8][8],
	    char name[256])
{
  int c, d;
  
  printf("%s = [ ...\n", name);
  for (c = 0; c < 8; c++) {
    if (c != 7){
    for (d = 0; d < 8; d++){
      if (d != 7){
      printf("%e,\t", m_to_print[c][d]);
      }
      else{
      printf("%e; ...\n", m_to_print[c][d]);
      }
    }
    }
    else{
    for (d = 0; d < 8; d++){
      if (d != 7){
      printf("%e,\t", m_to_print[c][d]);
      }
      else{
      printf("%e];\n", m_to_print[c][d]);
      }
    }
    }
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
  printf("%s: (%.8e, %.8e, %.8e)\n",name, v_to_print[0], v_to_print[1], v_to_print[2]);
  return 0;
}

int v_print6( double v_to_print[6],
	     char name[256])
{
  printf("%s: (%.6e, %.6e, %.6e) (%.6e, %.6e, %.6e)\n",name, v_to_print[0], v_to_print[1], v_to_print[2],v_to_print[3], v_to_print[4], v_to_print[5]);
  return 0;
}

int v_print7( double v_to_print[7],
	     char name[256])
{
  printf("%s: (%.7f, %.7f, %.7f) (%.7f, %.7f, %.7f) %.7f\n",name, v_to_print[0], v_to_print[1], v_to_print[2],v_to_print[3], v_to_print[4], v_to_print[5], v_to_print[6]);
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

// Copy matrices
int m_copy6( double m_out[6][6],
            double m_in[6][6])

{
    
  int i,j;
  for (i = 0; i < 6; i++){
    for (j = 0; j < 6; j++){
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

int m_x_v6bis(  double v_out[6],
            double m_in[6][6],
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




// Matrix x Vector
int m_x_v9(  double v_out[9],
            double **m_in,
            double v_in[9])
{


  int i,j;

  for (i = 0; i < 9; i++){
    v_out[i]= 0;
  }

  for (i = 0; i < 9; i++){
    v_out[i] = m_in[i][0]*v_in[0];
    for (j = 1 ; j < 9; j++){
      v_out[i] = v_out[i] + m_in[i][j]*v_in[j];
    }
  }

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

int m_x_m6bis( double m_out[6][6],
	   double m_in1[6][6],
	   double m_in2[6][6] )
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


// Matrix Transpose
int m_trans6(    double m_out[6][6],
                double m_in[6][6])
{
  int i,j;
  for (i = 0; i< 6; i++){
  for (j = 0; j< 6; j++){
    m_out[i][j] = m_in[j][i];

  }
  }

        return 0;
}

int equin_to_class(double *sma, double *ecc, double *inc, double *raan, double *arg_per, double *true_ano, double *mean_ano, double af, double ag, double l, double n, double chi, double psi, double mu, double et, double fr){ // chi looks like a x

  *sma = pow(mu/(n*n), 1./3);
  *ecc = sqrt(af*af + ag*ag);
  double p,q;

  *inc = M_PI * (1 - fr)/2. + 2*fr*atan(sqrt(chi*chi + psi*psi)); //Vallado 03: 2*atan(sqrt(chi*chi + psi*psi));
  *raan = atan2(chi,psi);
  *arg_per = atan2(ag, af) - fr * atan2(chi, psi); //Vallado 03:atan2(ag,af) - atan2(chi,psi);

  *mean_ano = l - atan2(ag, af); //Vallado 03:fmod(l - *arg_per - *raan, 2*M_PI);

  //to get the true ano, convert to j2000 intertial the  use cart2kep
  // equin to j2000
  double        r_cart[6];
  double r_eq[9];
  r_eq[0] = *sma*1e-3; // a = (mu / n^2)^(1/3) in km for eqncpv_c
  r_eq[1] = ag; r_eq[2] = af; r_eq[3] = l; r_eq[4] = chi; r_eq[5] = psi; r_eq[6] = 0; r_eq[7] = n; r_eq[8] = 0; //
  eqncpv_c(et, et, r_eq, -M_PI/2, M_PI/2., r_cart);
  double r_cart_r[3], r_cart_v[3];
  r_cart_r[0] = r_cart[0]; r_cart_r[1] = r_cart[1]; r_cart_r[2] = r_cart[2];
  r_cart_v[0] = r_cart[3]; r_cart_v[1] = r_cart[4]; r_cart_v[2] = r_cart[5];
  v_print(r_cart_r, "r_cart_r");
  v_print(r_cart_v, "r_cart_v");
  // j2000 to classical
  ORBITAL_ELEMENTS_T  OE;
  cart2kep(&OE, r_cart_r, r_cart_v, et, mu/1e9);// !!!!!! mu here needs to be in km3
  *true_ano = OE.f;
  //  printf("OE %f %e %f %f %f %f\n", OE.sma, OE.eccentricity, OE.inclination*180/M_PI, OE.long_an*180/M_PI, OE.w*180/M_PI, OE.f*180/M_PI);

  return 0;
}

// derivative matrix equinoctial to classical orbital elements (page 32 of http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.112.5780&rep=rep1&type=pdf)
// Vallado 2003
/* int compute_T_deriv_equin_to_class(double T_equin_2_class[6][6], double af, double ag, double l, double n, double chi, double psi, double mu){ // chi looks like a x */
/*   double den1 = af*af + ag*ag; */
/*   double den2 = (1 + chi*chi + psi*psi) * sqrt(chi*chi + psi*psi); */
/*     double den3 = chi*chi + psi*psi; */
/*     T_equin_2_class[0][0] = 0; T_equin_2_class[0][1] = af/sqrt(den1); T_equin_2_class[0][2] = 0; T_equin_2_class[0][3] = 0; T_equin_2_class[0][4] = -ag/den1; T_equin_2_class[0][5] = ag/den1; */
/*     T_equin_2_class[1][0] = 0; T_equin_2_class[1][1] = ag/sqrt(den1); T_equin_2_class[1][2] = 0; T_equin_2_class[1][3] = 0; T_equin_2_class[1][4] = af/den1; T_equin_2_class[1][5] = -af/den1; */
/*   T_equin_2_class[2][0] = 0; T_equin_2_class[2][1] = 0; T_equin_2_class[2][2] = 0; T_equin_2_class[2][3] = 0; T_equin_2_class[2][4] = 0; T_equin_2_class[2][5] = 1; */
/*   T_equin_2_class[3][0] = -2./(3*n)*pow(mu/(n*n),1./3); T_equin_2_class[3][1] = 0; T_equin_2_class[3][2] = 0; T_equin_2_class[3][3] = 0; T_equin_2_class[3][4] = 0; T_equin_2_class[3][5] = 0; */
/*   T_equin_2_class[4][0] = 0; T_equin_2_class[4][1] = 0; T_equin_2_class[4][2] = 2*chi/den2; T_equin_2_class[4][3] = psi/den3; T_equin_2_class[4][4] = -psi/den3; T_equin_2_class[4][5] = 0; */
/*   T_equin_2_class[5][0] = 0; T_equin_2_class[5][1] = 0; T_equin_2_class[5][2] = 2*psi/den2; T_equin_2_class[5][3] = -chi/den3; T_equin_2_class[5][4] = chi/den3; T_equin_2_class[5][5] = 0; */
  

/*   return 0; */
/* } */
// Vallado and Alfano 2015 "Updated Analytical Partials for Covariance Transformations and Optimization" (a few typos) and Danielson et al. 1995 "Semianalytic Satellite Theory" (that's how i found the typos in Alfano. Alfano numerical results are correct though but there are some typos in the expresssions of the equations
int compute_T_deriv_equin_to_class(double T_equin_2_class[6][6], double af, double ag, double l, double n, double chi, double psi, double mu, double fr){ // chi looks like a x
  double den1 = af*af + ag*ag;
  double den2 = (1 + chi*chi + psi*psi) * sqrt(chi*chi + psi*psi);
    double den3 = chi*chi + psi*psi;
    T_equin_2_class[0][0] = -2./(3*n)*pow(mu/(n*n),1./3); T_equin_2_class[0][1] = 0; T_equin_2_class[0][2] = 0; T_equin_2_class[0][3] = 0; T_equin_2_class[0][4] = 0; T_equin_2_class[0][5] = 0; // da/d
    T_equin_2_class[1][0] = 0; T_equin_2_class[1][1] = af/sqrt(den1); T_equin_2_class[1][2] = ag/sqrt(den1); T_equin_2_class[1][3] = 0; T_equin_2_class[1][4] = 0; T_equin_2_class[1][5] = 0; // de/d
    T_equin_2_class[2][0] = 0; T_equin_2_class[2][1] = 0; T_equin_2_class[2][2] = 0; T_equin_2_class[2][3] = 2*fr*chi/den2; T_equin_2_class[2][4] = 2*fr*psi/den2; T_equin_2_class[2][5] = 0; // di/d
    T_equin_2_class[3][0] = 0; T_equin_2_class[3][1] = 0; T_equin_2_class[3][2] = 0; T_equin_2_class[3][3] = psi/den3; T_equin_2_class[3][4] = -chi/den3; T_equin_2_class[3][5] = 0; // do/d
    T_equin_2_class[4][0] = 0; T_equin_2_class[4][1] = -ag/den1; T_equin_2_class[4][2] = af/den1; T_equin_2_class[4][3] = -fr*psi/den3; T_equin_2_class[4][4] = fr*chi/den3; T_equin_2_class[4][5] = 0; // dw/d
    T_equin_2_class[5][0] = 0; T_equin_2_class[5][1] = ag/den1; T_equin_2_class[5][2] = -af/den1; T_equin_2_class[5][3] = 0; T_equin_2_class[5][4] = 0; T_equin_2_class[5][5] = 1; // dM/d
  

  return 0;
}


// r v cartensian ECI to  equinoctial elements (Vallado and Alfano 2015 "Updated Analytical Partials for Covariance Transformations and Optimization" (a few typos) and Danielson et al. 1995 "Semianalytic Satellite Theory" (that's how i found the typos in Alfano. Alfano numerical results are correct though but there are some typos in the expresssions of the equations)
int cart_to_equin( double *af, double *ag, double *l, double *n, double *chi, double *psi, double mu,  double fr, double rvec[3], double vvec[3]){ // rvec and vvec in km and km/s. mu in km3

  // Begin algorithm
  // a, n, e vector, and w 
  double rmag;
  double vmag;
  double rdotv;
  v_mag( &rmag, rvec);
  v_mag( &vmag, vvec);
  v_dot(&rdotv, rvec, vvec);
  double a = 1.0 / ( (2.0/rmag) - ( (vmag*vmag)/mu ) );
  *n = sqrt(mu / (a*a*a)); //  a = pow(mu/(n*n), 1./3);

  double tempv1[3];
  double tempv2[3];

  v_scale(tempv1, vvec, rdotv);
    
  double coeff = vmag*vmag - (mu/rmag);
      double e_vector[3];
  v_scale(tempv2, rvec, coeff);
  v_sub( e_vector , tempv2, tempv1);

  v_scale(e_vector, e_vector, 1/mu);  

  double w[3], wtemp[3];
  v_cross(wtemp, rvec, vvec);
  v_norm(w, wtemp);

  // chi and psi (_e is [0], _q is [1], _w is [2])
  *chi = w[0] / (1 + fr * w[2]);
  *psi = -w[1] / (1 + fr*w[2]);

  // f and g
  double f[3], g[3];
  coeff = 1/(1 + (*chi)*(*chi) + (*psi)*(*psi));
  f[0] = (1-(*chi)*(*chi) + (*psi)*(*psi))*coeff; f[1] = 2*(*chi)*(*psi)*coeff; f[2] = -2*fr*(*chi)*coeff;
  g[0] = 2*fr*(*chi)*(*psi)*coeff; g[1] = (1+(*chi)*(*chi)-(*psi)*(*psi))*fr*coeff; g[2] = 2*(*psi)*coeff;

  // af and ag
  v_dot(&(*af), e_vector, f);
  v_dot(&(*ag), e_vector, g);

  // l (lamba_m in vallado 15)  
  double x, y;
  v_dot(&x, rvec, f);// rvec in km but then in sinf and cosf a is in km so need to to keep km here too
  v_dot(&y, rvec, g);
  double b = 1./(1 + sqrt(1 - (*ag)*(*ag) - (*af)*(*af)));
  double sinf = (*ag) + ((1 - (*ag)*(*ag)*b)*y - (*ag)*(*af)*b*x) / (a*sqrt(1 - (*ag)*(*ag) - (*af)*(*af)));
  double cosf = (*af) + ((1 - (*af)*(*af)*b)*x - (*ag)*(*af)*b*y) / (a*sqrt(1 - (*ag)*(*ag) - (*af)*(*af)));
  double F = atan2(sinf, cosf);
  *l = F + (*ag)*cosf -(*af)*sinf;

  return 0;
}

// equinoctial elements to r v cartensian ECI (Vallado and Alfano 2015 "Updated Analytical Partials for Covariance Transformations and Optimization" (a few typos) and Danielson et al. 1995 "Semianalytic Satellite Theory" (that's how i found the typos in Alfano. Alfano numerical results are correct though but there are some typos in the expresssions of the equations)


int equin_to_cart(double rvec[3], double vvec[3],  double af, double ag, double l, double n, double chi, double psi, double mu,  double fr){

  // solve for F
  double F;
  double fi, fiplus;
  double lguess;
  double err_max = 1e-15; 
  double err;
  int i, iter;
  fi = l;
  lguess = fi + ag * cos(fi) - af*sin(fi);
  err = fabs((l - lguess)/l);
  iter = 0;
  //  printf("XXXXXXXXXXX error[%d]: %.2f | %f %f %f %f %f %f %f %f\n", iter, err*100,  af,ag,  l,  n,  chi, psi, mu,  fr);
  while (err > err_max){
    fiplus = fi - (fi + ag*cos(fi) - af*sin(fi) - l) / (1 - ag*sin(fi) - af*cos(fi));
    fi = fiplus;
    lguess = fi + ag * cos(fi) - af*sin(fi);
    err = fabs((l - lguess)/l);
    iter = iter+1;
    //  printf("error[%d]: %e > %e | %d\n", iter, err, err_max, (err > err_max) );
  }
  F = fi;
  //    printf("F = %e | %f %f\n", F, l*180./M_PI, lguess*180./M_PI);

  double f[3], g[3], w[3];
  double fgwfac = 1 / (1 + chi*chi + psi*psi);
  f[0] = 1 - chi*chi + psi*psi; f[1] = 2*chi*psi; f[2] = -2*fr*chi;
  g[0] = 2*fr*chi*psi; g[1] = (1 + chi*chi - psi*psi)*fr; g[2] = 2*psi;
  w[0] = 2*chi; w[1] = -2*psi; w[2] = (1-chi*chi-psi*psi)*fr;

  for (i = 0;i < 3; i++){
    f[i] = f[i] * fgwfac;
    g[i] = g[i] * fgwfac;
    w[i] = w[i] * fgwfac;
  }
  double a = pow(mu/(n*n), 1./3);

  double B = sqrt(1 - ag*ag - af*af);
  double b = 1 / (1 + B);

  double sinl = ( (1 - af*af*b) * sin(F) + ag*af*b*cos(F) - ag  ) / (1 - ag*sin(F) - af*cos(F));
  double cosl = ( (1 - ag*ag*b) * cos(F) + ag*af*b*sin(F) - af  ) / (1 - ag*sin(F) - af*cos(F));



/*   r = a*(1 - af*sin(F) - ag*cos(F)); */
/*   double  */
    
  double x = a*( (1-ag*ag*b)*cos(F) + ag*af*b*sin(F) - af );
  double y = a*( (1-af*af*b)*sin(F) + ag*af*b*cos(F) - ag );

  double xdot = -n*a*(ag + sinl) / B;
  double ydot = n*a*(af + cosl) / B;


  for (i = 0; i < 3; i++){
    rvec[i] = x*f[i] + y*g[i];
    vvec[i] = xdot*f[i] + ydot*g[i];
  }

  return 0;
}
// derivative matrix equinoctial elements to cartesian position velocity (Vallado and Alfano 2015 "Updated Analytical Partials for Covariance Transformations and Optimization" (a few typos) and Danielson et al. 1995 "Semianalytic Satellite Theory" (that's how i found the typos in Alfano. Alfano numerical results are correct though but there are some typos in the expresssions of the equations)
int compute_T_deriv_equin_to_cart( double T_equin_to_cart[6][6], double af, double ag, double l, double n, double chi, double psi, double mu, double fr){

  double rvec[3], vvec[3];

  // solve for F
  double F;
  double fi, fiplus;
  double lguess;
  double err_max = 1e-50; 
  double err;
  int i, iter;
  fi = l;
  lguess = fi + ag * cos(fi) - af*sin(fi);
  err = fabs((l - lguess)/l);
  iter = 0;
  ///  printf("error[%d]: %.2f\n", iter, err*100);
  while (err > err_max){
    fiplus = fi - (fi + ag*cos(fi) - af*sin(fi) - l) / (1 - ag*sin(fi) - af*cos(fi));
    fi = fiplus;
    lguess = fi + ag * cos(fi) - af*sin(fi);
    err = fabs((l - lguess)/l);
    iter = iter+1;
    // printf("error[%d]: %e > %e | %d\n", iter, err, err_max, (err > err_max) );
  }
  F = fi;
  //  printf("F = %e | %f %f\n", F, l*180./M_PI, lguess*180./M_PI);

  double f[3], g[3], w[3];
  double fgwfac = 1 / (1 + chi*chi + psi*psi);
  f[0] = 1 - chi*chi + psi*psi; f[1] = 2*chi*psi; f[2] = -2*fr*chi;
  g[0] = 2*fr*chi*psi; g[1] = (1 + chi*chi - psi*psi)*fr; g[2] = 2*psi;
  w[0] = 2*chi; w[1] = -2*psi; w[2] = (1-chi*chi-psi*psi)*fr;

  for (i = 0;i < 3; i++){
    f[i] = f[i] * fgwfac;
    g[i] = g[i] * fgwfac;
    w[i] = w[i] * fgwfac;
  }
  double a = pow(mu/(n*n), 1./3);

  double B = sqrt(1 - ag*ag - af*af);
  double b = 1 / (1 + B);

  double sinl = ( (1 - af*af*b) * sin(F) + ag*af*b*cos(F) - ag  ) / (1 - ag*sin(F) - af*cos(F));
  double cosl = ( (1 - ag*ag*b) * cos(F) + ag*af*b*sin(F) - af  ) / (1 - ag*sin(F) - af*cos(F));



/*   r = a*(1 - af*sin(F) - ag*cos(F)); */
/*   double  */
    
  double x = a*( (1-ag*ag*b)*cos(F) + ag*af*b*sin(F) - af );
  double y = a*( (1-af*af*b)*sin(F) + ag*af*b*cos(F) - ag );

  double xdot = -n*a*(ag + sinl) / B;
  double ydot = n*a*(af + cosl) / B;


  for (i = 0; i < 3; i++){
    rvec[i] = x*f[i] + y*g[i];
    vvec[i] = xdot*f[i] + ydot*g[i];
  }
 
  double rx = rvec[0]; double ry = rvec[1]; double rz = rvec[2];
  double vx = vvec[0]; double vy = vvec[1]; double vz = vvec[2];
  double fe = f[0]; double fq = f[1]; double fw = f[2];
  double ge = g[0]; double gq = g[1]; double gw = g[2];
  double we = w[0]; double wq = w[1]; double ww = w[2];

  double r, v;
  v_mag(&r, rvec);// !!!!! not sure check with danielson95
  v_mag(&v, vvec);// !!!!! not sure check with danielson95


  double A = n * a*a;

  double C = 1 + chi*chi + psi*psi;
  double dxdaf = ag * xdot / (n * (1 + B)) + a * y * xdot / (A*B) - a;
  double dydaf = ag * ydot / (n * (1 + B)) - a * x * xdot / (A*B);
  double dxdag = - af * xdot / (n * (1 + B)) + a * y * ydot / (A*B);
  double dydag = - af * ydot / (n * (1 + B)) - a * x * ydot / (A*B) - a;
  double r3 = r*r*r;
  double a3 = a*a*a;
  double dxdotdaf = a*xdot*ydot/(A*B) - A/r3 * (a*ag*x/(1+B) + x*y/B);
  double dydotdaf = -a*xdot*xdot/(A*B) - A/r3 * (a*ag*y/(1+B) - x*x/B);
  double dxdotdag = a*ydot*ydot/(A*B) + A/r3 * (a*af*x/(1+B) - y*y/B);
  double dydotdag = -a*xdot*ydot/(A*B) + A/r3 * (a*af*y/(1+B) + x*y/B);

/*   // if the convention is to use derivative with respect to the sma (and not the mean motion) */
/*   T_equin_to_cart[0][0] = rx/a; T_equin_to_cart[0][1] = dxdaf*fe + dydaf*ge; T_equin_to_cart[0][2] = dxdag*fe + dydag*ge; T_equin_to_cart[0][3] = 2*(fr*psi*(y*fe - x*ge) - x*we)/C; T_equin_to_cart[0][4] = 2*fr*(chi*(x*ge - y*fe) + y*we)/C; T_equin_to_cart[0][5] = vx/n; // drx/d */
/*   T_equin_to_cart[1][0] = ry/a; T_equin_to_cart[1][1] = dxdaf*fq + dydaf*gq; T_equin_to_cart[1][2] = dxdag*fq + dydag*gq; T_equin_to_cart[1][3] = 2*(fr*psi*(y*fq - x*gq) - x*wq)/C; T_equin_to_cart[1][4] = 2*fr*(chi*(x*gq - y*fq) + y*wq)/C; T_equin_to_cart[1][5] = vy/n; // drx/d */
/*   T_equin_to_cart[2][0] = rz/a; T_equin_to_cart[2][1] = dxdaf*fw + dydaf*gw; T_equin_to_cart[2][2] = dxdag*fw + dydag*gw; T_equin_to_cart[2][3] = 2*(fr*psi*(y*fw - x*gw) - x*ww)/C; T_equin_to_cart[2][4] = 2*fr*(chi*(x*gw - y*fw) + y*ww)/C; T_equin_to_cart[2][5] = vz/n; // drx/d */
/*   T_equin_to_cart[3][0] = -vx/(2*a); T_equin_to_cart[3][1] = dxdotdaf*fe + dydotdaf*ge; T_equin_to_cart[3][2] = dxdotdag*fe + dydotdag*ge; T_equin_to_cart[3][3] = 2*(fr*psi*(ydot*fe - xdot*ge) - xdot*we)/C; T_equin_to_cart[3][4] = 2*fr*(chi*(xdot*ge - ydot*fe) + ydot*we)/C; T_equin_to_cart[3][5] = -n*a3*rx/r3; // drx/d */
/*   T_equin_to_cart[4][0] = -vy/(2*a); T_equin_to_cart[4][1] = dxdotdaf*fq + dydotdaf*gq; T_equin_to_cart[4][2] = dxdotdag*fq + dydotdag*gq; T_equin_to_cart[4][3] = 2*(fr*psi*(ydot*fq - xdot*gq) - xdot*wq)/C; T_equin_to_cart[4][4] = 2*fr*(chi*(xdot*gq - ydot*fq) + ydot*wq)/C; T_equin_to_cart[4][5] = -n*a3*ry/r3; // drx/d */
/*   T_equin_to_cart[5][0] = -vz/(2*a); T_equin_to_cart[5][1] = dxdotdaf*fw + dydotdaf*gw; T_equin_to_cart[5][2] = dxdotdag*fw + dydotdag*gw; T_equin_to_cart[5][3] = 2*(fr*psi*(ydot*fw - xdot*gw) - xdot*ww)/C; T_equin_to_cart[5][4] = 2*fr*(chi*(xdot*gw - ydot*fw) + ydot*ww)/C; T_equin_to_cart[5][5] = -n*a3*rz/r3; // drx/d */





  // if the convention is to use derivative with respect to the mean motion (and not the sma)
  T_equin_to_cart[0][0] = -2*rx/(3*n); T_equin_to_cart[0][1] = dxdaf*fe + dydaf*ge; T_equin_to_cart[0][2] = dxdag*fe + dydag*ge; T_equin_to_cart[0][3] = 2*(fr*psi*(y*fe - x*ge) - x*we)/C; T_equin_to_cart[0][4] = 2*fr*(chi*(x*ge - y*fe) + y*we)/C; T_equin_to_cart[0][5] = vx/n; // drx/d
  T_equin_to_cart[1][0] = -2*ry/(3*n); T_equin_to_cart[1][1] = dxdaf*fq + dydaf*gq; T_equin_to_cart[1][2] = dxdag*fq + dydag*gq; T_equin_to_cart[1][3] = 2*(fr*psi*(y*fq - x*gq) - x*wq)/C; T_equin_to_cart[1][4] = 2*fr*(chi*(x*gq - y*fq) + y*wq)/C; T_equin_to_cart[1][5] = vy/n; // drx/d
  T_equin_to_cart[2][0] = -2*rz/(3*n); T_equin_to_cart[2][1] = dxdaf*fw + dydaf*gw; T_equin_to_cart[2][2] = dxdag*fw + dydag*gw; T_equin_to_cart[2][3] = 2*(fr*psi*(y*fw - x*gw) - x*ww)/C; T_equin_to_cart[2][4] = 2*fr*(chi*(x*gw - y*fw) + y*ww)/C; T_equin_to_cart[2][5] = vz/n; // drx/d
  T_equin_to_cart[3][0] = vx/(3*n); T_equin_to_cart[3][1] = dxdotdaf*fe + dydotdaf*ge; T_equin_to_cart[3][2] = dxdotdag*fe + dydotdag*ge; T_equin_to_cart[3][3] = 2*(fr*psi*(ydot*fe - xdot*ge) - xdot*we)/C; T_equin_to_cart[3][4] = 2*fr*(chi*(xdot*ge - ydot*fe) + ydot*we)/C; T_equin_to_cart[3][5] = -n*a3*rx/r3; // drx/d
  T_equin_to_cart[4][0] = vy/(3*n); T_equin_to_cart[4][1] = dxdotdaf*fq + dydotdaf*gq; T_equin_to_cart[4][2] = dxdotdag*fq + dydotdag*gq; T_equin_to_cart[4][3] = 2*(fr*psi*(ydot*fq - xdot*gq) - xdot*wq)/C; T_equin_to_cart[4][4] = 2*fr*(chi*(xdot*gq - ydot*fq) + ydot*wq)/C; T_equin_to_cart[4][5] = -n*a3*ry/r3; // drx/d
  T_equin_to_cart[5][0] = vz/(3*n); T_equin_to_cart[5][1] = dxdotdaf*fw + dydotdaf*gw; T_equin_to_cart[5][2] = dxdotdag*fw + dydotdag*gw; T_equin_to_cart[5][3] = 2*(fr*psi*(ydot*fw - xdot*gw) - xdot*ww)/C; T_equin_to_cart[5][4] = 2*fr*(chi*(xdot*gw - ydot*fw) + ydot*ww)/C; T_equin_to_cart[5][5] = -n*a3*rz/r3; // drx/d


  return 0;
}

// derivative matrix classical orbital elements to cartesian position velocity (page 28-31 of http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.112.5780&rep=rep1&type=pdf)
int compute_T_deriv_class_to_cart( double T_class_to_cart[6][6], double a, double e, double i, double o, double w, double nu, double mu ){ //o i raan, w is arg_per, nu is true ano

  //  a = 6860.7631*1000; e = 0.0010640; i = 97.65184*M_PI/180; o = 79.54701*M_PI/180; w = 83.86041*M_PI/180; nu = 65.21303*M_PI/180; // !!!!! remoee

  double ri[6]; // it's not ri, just a notation
  double c1 = (1-e*e) / (1+e*cos(nu));

  double c2 = -(2*a*e + a*cos(nu) + a*e*e*cos(nu)) / ((1+e*cos(nu))*(1+e*cos(nu)));
  double c3 = c1 / (1+e*cos(nu));
  ri[0] = c1 * ( cos(nu)*(cos(o)*cos(w) - sin(o)*sin(w)*cos(i)) - sin(nu)*(cos(o)*sin(w) + sin(o)*cos(w)*cos(i)) );
  ri[1] = c2 * ( cos(nu)*(cos(o)*cos(w) - sin(o)*sin(w)*cos(i)) - sin(nu)*(cos(o)*sin(w) + sin(o)*cos(w)*cos(i)) );
  ri[2] = a*c1*sin(o)*sin(i)*( cos(nu)*sin(w) + sin(nu)*cos(w) );
  ri[3] = a*c1*(-cos(nu)*sin(o)*cos(w) - cos(nu)*cos(o)*sin(w)*cos(i) + sin(nu)*sin(o)*sin(w) - sin(nu)*cos(o)*cos(w)*cos(i));
  ri[4] = a*c1*(-cos(nu)*cos(o)*sin(w) - cos(nu)*sin(o)*cos(w)*cos(i) - sin(nu)*cos(o)*cos(w) + sin(nu)*sin(o)*sin(w)*cos(i));
  ri[5] = a*c3*( -sin(nu) * (cos(o)*cos(w) - sin(o)*sin(w)*cos(i)) + (e + cos(nu))*(-cos(o)*sin(w) - sin(o)*cos(w)*cos(i)) );

  double rj[6]; // it's not rj, just a notation
  double c4 = sqrt(mu / (a*(1-e*e)));
  rj[0] = c1 * ( cos(nu)*sin(o)*cos(w) + cos(nu)*cos(o)*sin(w)*cos(i) - sin(nu)*sin(o)*sin(w) + sin(nu)*cos(o)*cos(w)*cos(i) );
  rj[1] = c2 * ( cos(nu)*(sin(o)*cos(w) + cos(o)*sin(w)*cos(i)) + sin(nu)*(-sin(o)*sin(w) + cos(o)*cos(w)*cos(i)) );
  rj[2] = -a*c1*cos(o)*sin(i)*( cos(nu)*sin(w) + sin(nu)*cos(w) );
  rj[3] = a*c1*(cos(nu)*cos(o)*cos(w) - cos(nu)*sin(o)*sin(w)*cos(i) - sin(nu)*cos(o)*sin(w) - sin(nu)*sin(o)*cos(w)*cos(i));
  // line below from Vallado 2015 
  //rj[3] = a*c1*(cos(nu)*(sin(o)*cos(w) + cos(o)*sin(w)*cos(i)) + sin(nu)*(-sin(o)*sin(w) + cos(o)*cos(w)*cos(i)));
  rj[4] = a*c1*(-cos(nu)*sin(o)*sin(w) + cos(nu)*cos(o)*cos(w)*cos(i) - sin(nu)*sin(o)*cos(w) - sin(nu)*cos(o)*sin(w)*cos(i));
  rj[5] = a*c3*( -sin(nu) * (sin(o)*cos(w) + cos(o)*sin(w)*cos(i)) + (e + cos(nu))*(-sin(o)*sin(w) + cos(o)*cos(w)*cos(i)) );

  double rk[6]; // it's not rk, just a notation
  rk[0] = c1*(cos(nu)*sin(w)*sin(i) + sin(nu)*cos(w)*sin(i));
  rk[1] = c2*(cos(nu)*sin(w)*sin(i) + sin(nu)*cos(w)*sin(i));
  rk[2] = a*c1*cos(i)*(cos(nu)*sin(w) + sin(nu)*cos(w));
  rk[3] = 0;
  rk[4] = a*c1*sin(i)*(cos(nu)*cos(w) - sin(nu)*sin(w));
  rk[5] = a*c3*(-sin(nu)*sin(w)*sin(i) + (e+cos(nu))*cos(w)*sin(i));


  double vi[6]; // it's not vi, just a notation
  c3 = sqrt(mu/(a*(1-e*e)));
  vi[0] = 1/(2*a)*c3 * ( sin(nu)*(cos(o)*cos(w) - sin(o)*sin(w)*cos(i)) - (e+cos(nu)) * (-cos(o)*sin(w) - sin(o)*cos(w)*cos(i) ));//!!!! looks like there s a typo in the equation online cos((w*cos(i)))
  vi[1] = 1 / (1-e*e) * c3 * ( -e*sin(nu)*(cos(o)*cos(w) - sin(o)*sin(w)*cos(i)) + (1+e*cos(nu)) * (-cos(o)*sin(w) - sin(o)*cos(w)*cos(i) ));//!!!! looks like there s a typo in the equation online cos((w*cos(i)))
  vi[2] = c3*sin(o)*( sin(nu)*sin(w)*sin(i) + (e+cos(nu))*cos(w)*sin(i) );
  vi[3] = c3*(sin(nu)*sin(o)*cos(w) + sin(nu)*cos(o)*sin(w)*cos(i) + (e+cos(nu)) * (sin(o)*sin(w) - cos(o)*cos(w)*cos(i) ));
  vi[4] = c3*(sin(nu)*cos(o)*sin(w) + sin(nu)*sin(o)*cos(w)*cos(i) + (e+cos(nu)) * (-cos(o)*cos(w) + sin(o)*sin(w)*cos(i) ));
  vi[5] = c3*(-cos(nu)*(cos(o)*cos(w)-sin(o)*sin(w)*cos(i)) - sin(nu)*(-cos(o)*sin(w) - sin(o)*cos(w)*cos(i)));

  double vj[6]; // it's not vj, just a notation
  vj[0] = 1/(2*a)*c3 * ( sin(nu)*(sin(o)*cos(w) + cos(o)*sin(w)*cos(i)) - (e+cos(nu)) * (-sin(o)*sin(w) + cos(o)*cos(w)*cos(i) ));//!!!! looks like there s a typo in the equation online cos((w*cos(i)))
  vj[1] = 1 / (1-e*e) * c3 * ( -e*sin(nu)*(sin(o)*cos(w) + cos(o)*sin(w)*cos(i)) + (1+e*cos(nu)) * (-sin(o)*sin(w) + cos(o)*cos(w)*cos(i) ));//!!!! looks like there s a typo in the equation online cos((w*cos(i)))
  vj[2] = c3*( sin(nu)*cos(o)*sin(w)*sin(i) - (e+cos(nu))*cos(o)*cos(w)*sin(i) );
  vj[3] = c3*(-sin(nu)*cos(o)*cos(w) + sin(nu)*sin(o)*sin(w)*cos(i) + (e+cos(nu)) * (-cos(o)*sin(w) - sin(o)*cos(w)*cos(i) ));
  vj[4] = c3*(sin(nu)*sin(o)*sin(w) - sin(nu)*cos(o)*cos(w)*cos(i) + (e+cos(nu)) * (-sin(o)*cos(w) - cos(o)*sin(w)*cos(i) ));
  vj[5] = c3*(-cos(nu)*(sin(o)*cos(w)+cos(o)*sin(w)*cos(i)) - sin(nu)*(-sin(o)*sin(w) + cos(o)*cos(w)*cos(i)));


  double vk[6]; // it's not vk, just a notation
  vk[0] = 1/(2*a)*c3 * ( sin(nu)*sin(w)*sin(i) - (e+cos(nu)) * cos(w)*sin(i) );
  vk[1] = 1 / (1-e*e) * c3 * ( -e*sin(nu)*sin(w)*sin(i) + (1+e*cos(nu)) * cos(w)*sin(i));//!!!! looks like there s a typo in the equation online cos((w*cos(i)))
  vk[2] = c3*( -sin(nu)*sin(w)*cos(i) + (e+cos(nu))*cos(w)*cos(i) );
  vk[3] = 0;
  vk[4] = c3*(-sin(nu)*cos(w)*sin(i) - (e+cos(nu)) * (sin(w)*sin(i) ));
  vk[5] = c3*(-cos(nu)*sin(w)*sin(i) - sin(nu)*cos(w)*sin(i));


/*   v_print6(ri,"ri");  v_print6(rj,"rj");  v_print6(rk,"rk"); */
/*   v_print6(vi,"vi");  v_print6(vj,"vj");  v_print6(vk,"vk"); */
    int j;
  for (j = 0; j< 6; j++){
    T_class_to_cart[j][0] = ri[j];
    T_class_to_cart[j][1] = rj[j];
    T_class_to_cart[j][2] = rk[j];
    T_class_to_cart[j][3] = vi[j];
    T_class_to_cart[j][4] = vj[j];
    T_class_to_cart[j][5] = vk[j];
  }


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

// compute transformation matrix from inertial to ntw but 6 by 6 dimensinos. this is used to convrt inertial covariance matcies *position, velocity) from inertial to ntw coordinates.
// important: contraritly to compute_T_inrtl_2_ntw, the outputs is ntw (twn in compute_T_inrtl_2_ntw). this is to compare to Vallado matlab's script covct2o2.m (in Code/cygnss/collision/vallado/matlab or email from Vallado o July 9, 2018)
// source: Vallado matlab's script covct2o2.m (in Code/cygnss/collision/vallado/matlab or email from Vallado o July 9, 2018)
int compute_T_inrtl_2_ntw_6by6( double T_inrtl_2_ntw_6by6[6][6],
			   double r[3],
			   double v[3])
{



  // Declarations


  double t[3], w[3], n[3]; // t = in-track; w = cross-track; n = vector in the orbital plane normal to the velocity (directed outwards from the Earth)
  int col;

  // if lvlh
    /* Radial (1/2) */
    v_norm(n, r);
    v_scale(n, n, -1.0);

    /* Cross-track */
    v_cross(w, n, v);
    v_norm(w, w);
    v_scale(w, w, -1.0);
    
    /* Along-track (different from the velocity direction if the eccentricity is not 0)*/
    v_cross(t, w, n);
    v_norm(t, t);
    v_scale(t, t, -1.0);
    
    /* Radial (2/2) */
    v_scale(n, n, -1.0);   // this lign has been added by CBV on 07-29-2015 so that the basis (iv, w, n) is oriented w the right hand direction

/*   //if ntw */

/*   // Algorithm */
/*   /\* In-track (in the direction of the velocity even if the eccentricity is not 0 (contrarily to the along-track direction in the LVLH frame))*\/ */
/*   v_norm(t,v); */

/*   /\* Cross-track *\/ */
/*   v_cross(w,r,v); */
/*   v_norm(w,w); */

/*   /\* In-track cross cross-track *\/ */
/*   v_cross(n,t,w); */
 
  int i,j;
    for (i = 0; i < 6; i++){
      for (j = 0; j < 6; j++){
	T_inrtl_2_ntw_6by6[i][j] = 0;
      }
    }

        T_inrtl_2_ntw_6by6[0][0] = n[0];
        T_inrtl_2_ntw_6by6[0][1] = n[1];
        T_inrtl_2_ntw_6by6[0][2] = n[2];
        T_inrtl_2_ntw_6by6[1][0] = t[0];
        T_inrtl_2_ntw_6by6[1][1] = t[1];
        T_inrtl_2_ntw_6by6[1][2] = t[2];
        T_inrtl_2_ntw_6by6[2][0] = w[0];
        T_inrtl_2_ntw_6by6[2][1] = w[1];
        T_inrtl_2_ntw_6by6[2][2] = w[2];

        T_inrtl_2_ntw_6by6[3][3] = n[0];
        T_inrtl_2_ntw_6by6[3][4] = n[1];
        T_inrtl_2_ntw_6by6[3][5] = n[2];
        T_inrtl_2_ntw_6by6[4][3] = t[0];
        T_inrtl_2_ntw_6by6[4][4] = t[1];
        T_inrtl_2_ntw_6by6[4][5] = t[2];
        T_inrtl_2_ntw_6by6[5][3] = w[0];
        T_inrtl_2_ntw_6by6[5][4] = w[1];
        T_inrtl_2_ntw_6by6[5][5] = w[2];



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




   /* Inertial to Earth pressure frame. The Earth pressure frame is defined as follow: */
  /*  Notations: */
  /* - O center of the Earth */
  /* - C position of the satellite */
  /* - H projection of the satellite lcoation on the surface of the Earth (ie sub-satellite point) */
   /* - S position of the Sun */
   // - the z vector is the direction OH
  // - the y vector is the vector in the plane OHS and in the direction of the Sun
  // - the x vector completes the orthogonal basis
int compute_T_inrtl_2_earth_pres_frame( double T_inrtl_2_earth_pres_frame[3][3],
                            double r[3],
					double et)
{
    // Declarations
  double xep[3], yep[3], zep[3]; // "ep": Earth Pressure frame -> coordiantes of the vectors of the EP frame in the inertial refererence frame
  double xep_temp[3];
  double rmag;
    double x[6];
  double lt;
  double r_earth2sun_J2000[3];
  v_mag(&rmag, r);  	  
  spkez_c(10, et, "J2000", "NONE", 399, x, &lt); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration.
  r_earth2sun_J2000[0] = x[0];
  r_earth2sun_J2000[1] = x[1];
  r_earth2sun_J2000[2] = x[2];
  double r_earth2sun_J2000_norm[3];
  v_norm(r_earth2sun_J2000_norm, r_earth2sun_J2000);
  v_norm(zep, r);
  v_cross(xep_temp, zep, r_earth2sun_J2000_norm);
  double xep_temp_norm[3];
  v_norm(xep_temp_norm, xep_temp);
  double yep_temp[3];
  v_cross(yep_temp, zep, xep_temp_norm);
  double  cos_yep_temp_earth2sun;
  v_dot(&cos_yep_temp_earth2sun, yep_temp, r_earth2sun_J2000_norm);
  if (cos_yep_temp_earth2sun < 0){// we want yep to be in the direction of the Sun.
    v_scale(yep, yep_temp, -1);
  }
  else{
    v_scale(yep, yep_temp, 1);
  }
  v_cross(xep, yep, zep);

  // vecotr completing the orthogonal basis
    T_inrtl_2_earth_pres_frame[0][0] = xep[0]; 
    T_inrtl_2_earth_pres_frame[0][1] = xep[1];
    T_inrtl_2_earth_pres_frame[0][2] = xep[2];
    
    // vector in the plane OCS and in the direction of the Sun
    T_inrtl_2_earth_pres_frame[1][0] = yep[0];
    T_inrtl_2_earth_pres_frame[1][1] = yep[1];
    T_inrtl_2_earth_pres_frame[1][2] = yep[2];
    
    // vector Earth to satellite
    T_inrtl_2_earth_pres_frame[2][0] = zep[0]; 
    T_inrtl_2_earth_pres_frame[2][1] = zep[1]; 
    T_inrtl_2_earth_pres_frame[2][2] = zep[2]; 
    

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
                                                                                                                                                                                                                     
int compute_T_sc_to_lvlh(double T_sc_to_lvlh[3][3], double v_angle[3], int order_rotation[3], char attitude_profile[256], double *et,  double r_i2cg_INRTL[3], double v_i2cg_INRTL[3]	,int file_is_quaternion,  double quaternion[4], PARAMS_T *PARAMS){

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
  // the quaternion code was made to read the quaternion files of swarm at https://swarm-diss.eo.esa.int/#swarm%2FLevel1b%2FEntire_mission_data%2FSTRxATT
  double *p=malloc(sizeof(double)), *r=malloc(sizeof(double)), *y=malloc(sizeof(double));
  if (file_is_quaternion == 1){ // if attitude is set using quaternions
    // the quaternion I get from swarm are the rotation from sc to ecef so we need to convert from ecef to lvlh now
    /// // first calculate r_ecef2cg_ECEF and v_ecef2cg_ECEF from tehir inertial coordinates
    SpiceDouble       xform[6][6];
    double estate[6], jstate[6];
    double r_ecef2cg_ECEF[3], v_ecef2cg_ECEF[3];
    estate[0] = r_i2cg_INRTL[0];estate[1] = r_i2cg_INRTL[1];estate[2] = r_i2cg_INRTL[2];
    estate[3] = v_i2cg_INRTL[0];estate[4] = v_i2cg_INRTL[1];estate[5] = v_i2cg_INRTL[2];
    sxform_c (  "J2000", PARAMS->EARTH.earth_fixed_frame,  *et,    xform  );
    mxvg_c   (  xform,       estate,   6,  6, jstate );
    r_ecef2cg_ECEF[0] = jstate[0]; r_ecef2cg_ECEF[1] = jstate[1]; r_ecef2cg_ECEF[2] = jstate[2];
    v_ecef2cg_ECEF[0] = jstate[3]; v_ecef2cg_ECEF[1] = jstate[4]; v_ecef2cg_ECEF[2] = jstate[5];

    // // then compute the transformation matrix ecef to lvlh (use the equations of the transformation eci to lvlh but it's the same equations, it's just the inputs r and v that are now ecef instead of eci)
    double T_ecef_2_lvlh[3][3];
    compute_T_inrtl_2_lvlh( T_ecef_2_lvlh, r_ecef2cg_ECEF, v_ecef2cg_ECEF);
    /* double T_ecef_2_ntw[3][3]; */
    /* compute_T_inrtl_2_ntw( T_ecef_2_ntw, r_ecef2cg_ECEF, v_ecef2cg_ECEF); */

    // // then compute the transformation matrix from sc to ecef from the quaternions sc to ecef
    double T_ecef_to_sc_temp[3][3];
    q2m_c(quaternion,T_ecef_to_sc_temp);//!!!!!!q2m_c(quaternion,T_sc_to_ecef);

    double T_ecef_to_sc[3][3];
    for (i = 0; i < 3; i++){ // this i sbecause in swarm the x body is the same as the x body of spock but the y body swarm is minus the y body spock and the z body swarm is minus the z body spock
      T_ecef_to_sc[0][i] = T_ecef_to_sc_temp[0][i];
      T_ecef_to_sc[1][i] = -T_ecef_to_sc_temp[1][i];
      T_ecef_to_sc[2][i] = -T_ecef_to_sc_temp[2][i];
      }
    
    double T_sc_to_ecef[3][3];
    m_trans(T_sc_to_ecef, T_ecef_to_sc);
    // // finally, compute the transformation matrix sc to lvlh
   
    m_x_m(T_sc_to_lvlh, T_ecef_2_lvlh, T_sc_to_ecef);// don't ask me why it's not  m_x_m(T_sc_to_lvlh, T_sc_to_ecef, T_ecef_2_lvlh)...

    /* etprint(*et, "time"); */
    /*     m_print(T_sc_to_lvlh, "T_sc_to_lvlh"); */
    /* // TEST */
    /* double x_sat_in_ecef[3], x_sat_in_body[3], r_ecef2cg_ECEF_normalized[3], r_ecef2cg_ECEF_normalized_in_lvlh[3], z_sat_in_body[3], z_sat_in_lvlh[3], z_lvlh_in_lvlh[3], x_lvlh_in_lvlh[3], x_sat_in_lvlh[3], y_sat_in_lvlh[3], y_lvlh_in_lvlh[3], y_sat_in_body[3], z_sat_in_ecef[3], v_ecef2cg_ECEF_normalized[3], v_ecef2cg_ECEF_normalized_in_lvlh[3]; */
    /* x_sat_in_body[0] = 1; x_sat_in_body[1] = 0; x_sat_in_body[2] = 0; */
    /* z_sat_in_body[0] = 0; z_sat_in_body[1] = 0; z_sat_in_body[2] = 1; */
    /* y_sat_in_body[0] = 0; y_sat_in_body[1] = 1; y_sat_in_body[2] = 0; */
    /* z_lvlh_in_lvlh[0] = 0; z_lvlh_in_lvlh[1] = 0; z_lvlh_in_lvlh[2] = 1; */
    /* x_lvlh_in_lvlh[0] = 1; x_lvlh_in_lvlh[1] = 0; x_lvlh_in_lvlh[2] = 0; */
    /* y_lvlh_in_lvlh[0] = 0; y_lvlh_in_lvlh[1] = 1; y_lvlh_in_lvlh[2] = 0; */
    /* m_x_v(x_sat_in_ecef, T_sc_to_ecef, x_sat_in_body); */
    /* m_x_v(z_sat_in_ecef, T_sc_to_ecef, z_sat_in_body); */
    /* v_print(x_sat_in_ecef, "x_sat_in_ecef"); */
    /* v_norm(r_ecef2cg_ECEF_normalized, r_ecef2cg_ECEF); */
    /* v_norm(v_ecef2cg_ECEF_normalized, v_ecef2cg_ECEF); */
    /* v_print(r_ecef2cg_ECEF_normalized, "r_ecef2cg_ECEF_normalized"); */
    /* v_print(z_sat_in_ecef, "z_sat_in_ecef should be similar to r_ecef2cg_ECEF_normalized"); */
    /* m_x_v(r_ecef2cg_ECEF_normalized_in_lvlh, T_ecef_2_lvlh, r_ecef2cg_ECEF_normalized); */
    /* m_x_v(v_ecef2cg_ECEF_normalized_in_lvlh, T_ecef_2_lvlh, v_ecef2cg_ECEF_normalized); */
    /* v_print(r_ecef2cg_ECEF_normalized_in_lvlh, "r_ecef2cg_ECEF_normalized_in_lvlh should be 0 0 1"); */
    /* v_print(v_ecef2cg_ECEF_normalized_in_lvlh, "r_ecef2cg_ECEF_normalized_in_lvlh should close (not equal) to 1 0 0"); */
    /* m_x_v(z_sat_in_lvlh, T_sc_to_lvlh, z_sat_in_body); */
    /* v_print(z_sat_in_lvlh, "z_sat_in_lvlh");// if sc was nadir (ie pitch roll yaw wrt lvlh are 0 then this should be 0 0 1) */
    /* double z_sat_in_lvlh_dot_z_lvlh_in_lvlh, x_sat_in_lvlh_dot_x_lvlh_in_lvlh; */
    /* v_dot(&z_sat_in_lvlh_dot_z_lvlh_in_lvlh, z_sat_in_lvlh, z_lvlh_in_lvlh); */
    /* printf("dot z %f\n", z_sat_in_lvlh_dot_z_lvlh_in_lvlh); */
    /* m_x_v(x_sat_in_lvlh, T_sc_to_lvlh, x_sat_in_body); */
    /* m_x_v(y_sat_in_lvlh, T_sc_to_lvlh, y_sat_in_body); */
    /* v_print(x_sat_in_lvlh, "x_sat_in_lvlh"); */
    /* v_print(y_sat_in_lvlh, "y_sat_in_lvlh"); */
    /* v_dot(&x_sat_in_lvlh_dot_x_lvlh_in_lvlh, x_sat_in_lvlh, x_lvlh_in_lvlh); */
    /* printf("dot x %f\n", x_sat_in_lvlh_dot_x_lvlh_in_lvlh); */
    // end of TEST

    // OTHER OLDER TEST BELOW
    /* double ex_ntw_in_ecef_coord[3], ex_ntw_in_ntw_coord[3]; */
    /* v_norm(ex_ntw_in_ecef_coord, v_ecef2cg_ECEF); */
    /* v_print(ex_ntw_in_ecef_coord, "ex_ntw_in_ecef_coord"); */
    /* m_x_v(ex_ntw_in_ntw_coord, T_ecef_2_ntw, ex_ntw_in_ecef_coord); */
    /* v_print(ex_ntw_in_ntw_coord, "ex_ntw_in_ntw_coord"); */
    /* double T_sc_to_ntw[3][3]; */
    /* m_x_m(T_sc_to_ntw, T_sc_to_ecef, T_ecef_2_ntw); */
    
    /* double x_sat_in_body[3], x_ecef1[3], x_sat_in_ntw[3]; */
    /* x_sat_in_body[0] = 1; x_sat_in_body[1] = 0; x_sat_in_body[2] = 0; */
    /*  m_x_v(x_ecef1, T_sc_to_ecef, x_sat_in_body); */
    /* m_x_v(x_sat_in_ntw, T_sc_to_ntw, x_sat_in_body); */
    /* double x_sat_in_ntw2[3]; */
    /* m_x_v(x_sat_in_ntw2, T_ecef_2_ntw, x_ecef1); */

    /* m_print(T_body_to_ntw, "T_body_to_ntw"); */
    /* m_print(T_sc_to_ntw, "T_sc_to_ntw"); */
    /* double x_sat_in_ntw3[3]; */
    /*     double T_body_to_ntw[3][3]; */
    /* m_x_m(T_body_to_ntw, T_ecef_2_ntw, T_sc_to_ecef); */

    /* m_x_v(x_sat_in_ntw3, T_body_to_ntw, x_sat_in_body); */
    /* v_print(x_sat_in_ntw, "x_sat_in_ntw"); */
    /* //        v_print(x_ecef1, "x_ecef1"); */
    /* v_print(x_sat_in_ntw2, "x_sat_in_ntw2"); */
    /* v_print(x_sat_in_ntw3, "x_sat_in_ntw3"); */
    /* double x_sat_dot_x_ntw; */
    /* v_dot(&x_sat_dot_x_ntw, x_sat_in_ntw3, ex_ntw_in_ntw_coord); */
    /* printf("dot %f\n", x_sat_dot_x_ntw); */
    /* //    v_print(v_ecef2cg_ECEF, "v_ecef2cg_ECEF"); */
    /* //    printf("%f %f %f %f\n", quaternion[0], quaternion[1], quaternion[2], quaternion[3]); */
    /* //m_print(T_sc_to_ecef, "T_sc_to_ecef"); */
    /*   //        m_print(T_sc_to_lvlh, "T_sc_to_lvlh");// ++++++ */
    /* /\* /\\* m_print(T_sc_to_lvlh_bis, "quat T_sc_to_lvlh"); *\\/ *\/ */
    // END OF  OTHER OLDER TEST BELOW
    //      print_test();exitf();

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
