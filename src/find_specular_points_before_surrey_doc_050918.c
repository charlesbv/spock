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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "stdlib.h"
#include <pwd.h>
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h>
//#include "SpiceUsr.h"
int iDebug;

int nProcs;
int iProc;

int nSpecularPointsMax = 6;

double re    = 6378000.0;
double rpole = 6378000.0;
double req   = 6378000.0;
double dtor  =       0.017453292;
double pi    =       3.1415927;
double twopi =       6.2831855;

void lam_darwin_malloc_linker_hack(){};

typedef struct {

  double          *pitch;
  double          *roll;
  double          *yaw;
  int             *order_pitch, *order_roll, *order_yaw;
  double          *et_interpo;     // for linear interpolation
  int nb_time_steps;
  double dt_prop;
  char initial_epoch[300];
  char final_epoch[300];
  double et_initial;
  double et_final;
  char attitude_profile[500];
  double **quaternion;
  int file_is_quaternion;
} ATTITUDE_T;

struct specular_point_info {

  double delta;
  double tau;
  double beta;
  double sigma;
  double theta;
  double a;
  double e;
  double d;
  double diff;
  double latitude;
  double longitude;
  double radius;
  int    nIters;

  double elev; // cbv 04-27-18
  double azim; // cbv 04-27-18

  double DistToTangentCygnss;
  double DistToTangentGps;
  double PercentAlongTangent;

  double CygnssToSpecularPoint;
  double SpecularPointToGPS;

  // This is the X, Y, Z coordinate of the Specular Point
  // in the space craft coordinate system:
  double SPscX, SPscY, SPscZ;

  double CygnssLOSVel;
  double GpsLOSVel;

  double gain;
  double power;
  double NormPower;
  double Incidence;

  float CygX;
  float CygY;
  float CygZ;

  float GpsX;
  float GpsY;
  float GpsZ;

  float SpX;
  float SpY;
  float SpZ;

};

enum {nPtsGrid=1000};

struct specular_point_info CourseGrid[nPtsGrid];

enum {nTimesMax=604800}; // used to be: 1209600    (14 days)
enum {nGps=31};
enum {nCygnss=16};

struct SatelliteDetails {

  int IsGood;
  int nPts;


   double lat[nTimesMax];
  double lon[nTimesMax];
  double radius[nTimesMax];
  double x[nTimesMax];
  double y[nTimesMax];
  double z[nTimesMax];
  double vx[nTimesMax];
  double vy[nTimesMax];
  double vz[nTimesMax];
  double time[nTimesMax];
  double yaw[nTimesMax];
  int nYawTimes;
  char name[200] ; // CBV 05/08/16
  int time_ymdhmsm[nTimesMax][7];  // CBV 05/09/16

};

enum {nThetaAntenna = 18}; // cbv new antenna onboard. before: 36
enum {nPhiAntenna   = 24}; // cbv new antenna onboard. before: 73
enum {nWinds        = 70};
enum {nAngles       = 90};

struct AntennaInfoStructure {

  float ptx;
  float gtx;
  float a0;
  float lambda;
  float noisefloor;
  float factor;
  float phi[nPhiAntenna];
  float theta[nThetaAntenna];
  float gain[2][nPhiAntenna][nThetaAntenna]; // cbv new antenna onboard. 2 for 2 antennas (port and starboard)
  float dTheta, dPhi;

  float Sigma0[nWinds][nAngles];
  float dIncidence;

  char AntennaFileNameStarboard[150]; // cbv new antenna onboard. before: char AntennaFileName
  char AntennaFileNamePort[150]; // cbv new antenna onboard. before: nothing
  char Sigma0FileName[150];
  char AntennaInfoFileName[150];

  int nAnt;
  float rolls[10];
  float yaws[10];
  float pitches[10];

};

struct AntennaInfoStructure AntennaInfo;


/*----------------------------------------------------------------------------

  returns the middle of a string.

----------------------------------------------------------------------------*/

char* strmid(char *instr, int start, int length) {

  char *out, *startptr;
  int i;

  out = malloc(sizeof(char)*(length+2));
  startptr = out;

  for (i=0; i<start; i++) instr++;

  for (i=0; i<length; i++) {
    *out = *instr;
    out++; instr++;
  }

  *out = '\0';

  return startptr;

}

//-----------------------------------------------------------------------
//
//-----------------------------------------------------------------------

char* c_int_str(int value, int length) {

  int ndig, dec, sign;
  char *pstr, *buff, *start, *temp;
  int c, strl;

  ndig = 1;
  dec = 2;
  sign = 0;

  pstr = fcvt((double) value, ndig, &dec, &sign);
  pstr = strmid(pstr,0,strlen(pstr)-1);

  if (dec > length) length = dec;

  if ( (buff=malloc(sizeof(char)*(length+2))) == NULL ) {
    printf("Not enough memory for the conversion!\n");
    start = pstr; }
  else {

    start = buff;
    c = 0;

    if (strlen(pstr) == 0) {
      *buff = '0';
      buff++; c++;
    } else {
      if (value == 0) c = strlen(pstr);
    }

    if (sign) {
      *buff = '-';
      buff++; c++;
    }

    while (*pstr || c < length-dec) {
      if (c < length-dec) {
	*buff = '0';
	buff++; c++;
      } else {
	*buff = *pstr;
	buff++; pstr++; c++;
      }
    }

    *buff = '\0';
  }

  strl = strlen(start);
  if (strl < length) {
    temp=malloc(sizeof(char)*(length+2));
    strcpy(temp,"0");
    strcat(temp,start);
    strcpy(start,temp);
    free(temp);
  }

  return start;

}

//-----------------------------------------------------------------------
//
//-----------------------------------------------------------------------

int read_yaw_file(char YawFileName[150],
		  struct SatelliteDetails CYGNSS[nCygnss]) {

  // This is super limited at this time:
  // 1. Assumes only a single satellite.
  // 2. Assumes that dt in satellite files and attitude file are the same
  // 3. Start times are identical.
  // 4. Attitude repeats over and over again.

  FILE *fpYaw;
  int TmpTime, t0, dt;
  float TmpYaw;
  int iError, i, TmpX, TmpY;
  char line[600];

  iError = 0;

  fpYaw = fopen(YawFileName,"r");
  if (!fpYaw) {
    printf("Can not read yaw file : %s\n",YawFileName);
    iError = 1;
    //newstructure
// return;
//newstructure
return 0;
  }

  i = 0;
  while (fscanf(fpYaw, "%d %f %f %f", &TmpTime, &TmpYaw, &TmpX, &TmpY) > 0) {
    //printf("TmpTime: %d %f\n",TmpTime,TmpYaw);
    CYGNSS[0].yaw[i] = TmpYaw;
    if (i == 0) t0 = TmpTime;
    if (i == 1) dt = TmpTime-t0;
    i++;
  }
  CYGNSS[0].nYawTimes = i;

  close(fpYaw);

  return iError;

}


// cbv new antenna onboard. function read_antenna_info replaced the one below it
//-----------------------------------------------------------------------
//
//-----------------------------------------------------------------------

int read_antenna_info() {

  FILE *fpAntenna, *fpSigma0;
  int i, j;
  int iError;
  float tmp;
  char line[600];

  iError = 0;

  AntennaInfo.ptx=25.0;
  AntennaInfo.gtx=20.0;
  AntennaInfo.a0=6.25e8;
  AntennaInfo.lambda=0.1904;
  AntennaInfo.noisefloor=5.8e-19;

  char *line_file = NULL; size_t len = 0; char text[1000];
  int iant;
  AntennaInfo.dTheta = 90.0/nThetaAntenna;
  AntennaInfo.dPhi = 360.0/nPhiAntenna;

  // starboard antenna
  iant = 0;
  fpAntenna = fopen(AntennaInfo.AntennaFileNameStarboard,"r");
  if (!fpAntenna) {
    printf("Can not read antenna file : %s\n",AntennaInfo.AntennaFileNameStarboard);
    iError = 1;
    //newstructure
    // return;
    //newstructure
    return 0;
  }


  for (j=0; j<nThetaAntenna; j++) { // row theta
    getline(&line_file, &len, fpAntenna);
    //    RemoveSpaces(line_file);
    strtok(line_file, "\n");strtok(line_file, "\r");

    //   printf("\n\nhh <%s>\n",line_file);
    sscanf(line_file, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", &AntennaInfo.gain[iant][0][j], &AntennaInfo.gain[iant][1][j], &AntennaInfo.gain[iant][2][j], &AntennaInfo.gain[iant][3][j], &AntennaInfo.gain[iant][4][j], &AntennaInfo.gain[iant][5][j], &AntennaInfo.gain[iant][6][j], &AntennaInfo.gain[iant][7][j], &AntennaInfo.gain[iant][8][j], &AntennaInfo.gain[iant][9][j], &AntennaInfo.gain[iant][10][j], &AntennaInfo.gain[iant][11][j], &AntennaInfo.gain[iant][12][j], &AntennaInfo.gain[iant][13][j], &AntennaInfo.gain[iant][14][j], &AntennaInfo.gain[iant][15][j], &AntennaInfo.gain[iant][16][j], &AntennaInfo.gain[iant][17][j], &AntennaInfo.gain[iant][18][j], &AntennaInfo.gain[iant][19][j], &AntennaInfo.gain[iant][20][j], &AntennaInfo.gain[iant][21][j], &AntennaInfo.gain[iant][22][j], &AntennaInfo.gain[iant][23][j]);


    AntennaInfo.theta[j] = 0 + j * AntennaInfo.dTheta;
    //printf("%f | ", AntennaInfo.theta[j]);
    for (i=0; i<nPhiAntenna; i++) { // column phi */
      AntennaInfo.phi[i] = 0 + i * AntennaInfo.dPhi; 
      //      printf("%f ", AntennaInfo.phi[i]);
	}
    //    printf("\n\n");

  }

  close(fpAntenna);


  // port antenna
  iant = 1;
  fpAntenna = fopen(AntennaInfo.AntennaFileNamePort,"r");
  if (!fpAntenna) {
    printf("Can not read antenna file : %s\n",AntennaInfo.AntennaFileNamePort);
    iError = 1;
    //newstructure
    // return;
    //newstructure
    return 0;
  }


  for (j=0; j<nThetaAntenna; j++) { // row theta
    getline(&line_file, &len, fpAntenna);
    //    RemoveSpaces(line_file);
    strtok(line_file, "\n");strtok(line_file, "\r");

    //   printf("\n\nhh <%s>\n",line_file);
    sscanf(line_file, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", &AntennaInfo.gain[iant][0][j], &AntennaInfo.gain[iant][1][j], &AntennaInfo.gain[iant][2][j], &AntennaInfo.gain[iant][3][j], &AntennaInfo.gain[iant][4][j], &AntennaInfo.gain[iant][5][j], &AntennaInfo.gain[iant][6][j], &AntennaInfo.gain[iant][7][j], &AntennaInfo.gain[iant][8][j], &AntennaInfo.gain[iant][9][j], &AntennaInfo.gain[iant][10][j], &AntennaInfo.gain[iant][11][j], &AntennaInfo.gain[iant][12][j], &AntennaInfo.gain[iant][13][j], &AntennaInfo.gain[iant][14][j], &AntennaInfo.gain[iant][15][j], &AntennaInfo.gain[iant][16][j], &AntennaInfo.gain[iant][17][j], &AntennaInfo.gain[iant][18][j], &AntennaInfo.gain[iant][19][j], &AntennaInfo.gain[iant][20][j], &AntennaInfo.gain[iant][21][j], &AntennaInfo.gain[iant][22][j], &AntennaInfo.gain[iant][23][j]);
    //    printf("%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", AntennaInfo.gain[iant][0][j], AntennaInfo.gain[iant][1][j], AntennaInfo.gain[iant][2][j], AntennaInfo.gain[iant][3][j], AntennaInfo.gain[iant][4][j], AntennaInfo.gain[iant][5][j], AntennaInfo.gain[iant][6][j], AntennaInfo.gain[iant][7][j], AntennaInfo.gain[iant][8][j], AntennaInfo.gain[iant][9][j], AntennaInfo.gain[iant][10][j], AntennaInfo.gain[iant][11][j], AntennaInfo.gain[iant][12][j], AntennaInfo.gain[iant][13][j], AntennaInfo.gain[iant][14][j], AntennaInfo.gain[iant][15][j], AntennaInfo.gain[iant][16][j], AntennaInfo.gain[iant][17][j], AntennaInfo.gain[iant][18][j], AntennaInfo.gain[iant][19][j], AntennaInfo.gain[iant][20][j], AntennaInfo.gain[iant][21][j], AntennaInfo.gain[iant][22][j], AntennaInfo.gain[iant][23][j]);


  }

  close(fpAntenna);



  
  AntennaInfo.factor = AntennaInfo.ptx * AntennaInfo.gtx *
    AntennaInfo.a0 * (AntennaInfo.lambda * AntennaInfo.lambda) /
    (4*4*pi*pi);

  fpSigma0 = fopen(AntennaInfo.Sigma0FileName,"r");
  if (!fpSigma0) {
    printf("Can not read Sigma0 file : %s\n",AntennaInfo.Sigma0FileName);
    iError = 1;
    //newstructure
    // return;
    //newstructure
    return 0;
  }

  fread(&i,sizeof(int), 1, fpSigma0);
  fread(&j,sizeof(int), 1, fpSigma0);
  //printf("%d %d\n",i,j);

  for (i=0; i<nWinds; i++) {
    for (j=0; j<nAngles; j++) {
      fread(&tmp,sizeof(float),1,fpSigma0);
      AntennaInfo.Sigma0[i][j]=tmp;
      //printf("%d %d %f\n",i,j,tmp);
    }
  }
  close(fpSigma0);

  AntennaInfo.dIncidence = 90.0/nAngles;

  fpAntenna = fopen(AntennaInfo.AntennaInfoFileName,"r");


  if (!fpAntenna) {
    printf("Can not read antenna file : %s\n",AntennaInfo.AntennaInfoFileName);
    iError = 1;
    //newstructure
    // return;
    //newstructure
    return 0;
  }

  int IsDone = 0;

  while (!IsDone) {

    if (fgets(line,500,fpAntenna) == NULL) IsDone = 1;

    if (strstr(line,"#N") != NULL) {
      fscanf(fpAntenna, "%d", &AntennaInfo.nAnt);
      //printf("%d\n",AntennaInfo.nAnt);
    }
    if (strstr(line,"#YAW") != NULL) {
      for (i=0;i<AntennaInfo.nAnt;i++) {
	fscanf(fpAntenna, "%f", &AntennaInfo.yaws[i]);
	//printf("%d %f\n",i,AntennaInfo.yaws[i]);
      }
    }
    if (strstr(line,"#ROLL") != NULL) {
      for (i=0;i<AntennaInfo.nAnt;i++) {
	fscanf(fpAntenna, "%f", &AntennaInfo.rolls[i]);
	//printf("%d %f\n",i,AntennaInfo.rolls[i]);
      }
    }
    if (strstr(line,"#PITCH") != NULL) {
      for (i=0;i<AntennaInfo.nAnt;i++) {
	fscanf(fpAntenna, "%f", &AntennaInfo.pitches[i]);
	//printf("%d %f\n",i,AntennaInfo.pitches[i]);
      }
    }

  }

  close(fpAntenna);



  return iError;

}

// cbv new antenna onboard. function below read_antenna_info below

/* //----------------------------------------------------------------------- */
/* // */
/* //----------------------------------------------------------------------- */

/* int read_antenna_info() { */

/*   FILE *fpAntenna, *fpSigma0; */
/*   int i, j; */
/*   int iError; */
/*   float tmp; */
/*   char line[600]; */

/*   iError = 0; */

/*   AntennaInfo.ptx=25.0; */
/*   AntennaInfo.gtx=20.0; */
/*   AntennaInfo.a0=6.25e8; */
/*   AntennaInfo.lambda=0.1904; */
/*   AntennaInfo.noisefloor=5.8e-19; */

/*   fpAntenna = fopen(AntennaInfo.AntennaFileName,"r"); */

/*   if (!fpAntenna) { */
/*     printf("Can not read antenna file : %s\n",AntennaInfo.AntennaFileName); */
/*     iError = 1; */
/*     //newstructure */
/* // return; */
/* //newstructure */
/* return 0; */
/*   } */

/*   fread(&i,sizeof(int), 1, fpAntenna); */
/*   fread(&j,sizeof(int), 1, fpAntenna); */
/*   //printf("%d %d\n",i,j); */

/*   for (i=0; i<nPhiAntenna; i++) { */
/*     for (j=0; j<nThetaAntenna; j++) { */
/*       fread(&tmp,sizeof(float),1,fpAntenna); */
/*       AntennaInfo.phi[i] = tmp; */
/*       fread(&tmp,sizeof(float),1,fpAntenna); */
/*       AntennaInfo.theta[j] = tmp; */
/*       fread(&tmp,sizeof(float),1,fpAntenna); */
/*       AntennaInfo.gain[i][j] = tmp; */
/*       //fread(AntennaInfo.phi[i][j],sizeof(float),1,fpAntenna); */
/*       //fread(AntennaInfo.theta[i][j],sizeof(float),1,fpAntenna); */
/*       //fread(AntennaInfo.gain[i][j],sizeof(float),1,fpAntenna); */

/*       //if (AntennaInfo.gain[i][j] > 20.0) printf("%d %d %f %f %f\n",i,j,AntennaInfo.phi[i], AntennaInfo.theta[j], AntennaInfo.gain[i][j]); */

/*     } */
/*   } */
/*   close(fpAntenna); */
/*   AntennaInfo.dTheta = 180.0/nThetaAntenna; */
/*   AntennaInfo.dPhi = AntennaInfo.dTheta; */

  
/*   AntennaInfo.factor = AntennaInfo.ptx * AntennaInfo.gtx * */
/*     AntennaInfo.a0 * (AntennaInfo.lambda * AntennaInfo.lambda) / */
/*     (4*4*pi*pi); */

/*   fpSigma0 = fopen(AntennaInfo.Sigma0FileName,"r"); */
/*   if (!fpSigma0) { */
/*     printf("Can not read Sigma0 file : %s\n",AntennaInfo.Sigma0FileName); */
/*     iError = 1; */
/*     //newstructure */
/* // return; */
/* //newstructure */
/* return 0; */
/*   } */

/*   fread(&i,sizeof(int), 1, fpSigma0); */
/*   fread(&j,sizeof(int), 1, fpSigma0); */
/*   //printf("%d %d\n",i,j); */

/*   for (i=0; i<nWinds; i++) { */
/*     for (j=0; j<nAngles; j++) { */
/*       fread(&tmp,sizeof(float),1,fpSigma0); */
/*       AntennaInfo.Sigma0[i][j]=tmp; */
/*       //printf("%d %d %f\n",i,j,tmp); */
/*     } */
/*   } */
/*   close(fpSigma0); */

/*   AntennaInfo.dIncidence = 90.0/nAngles; */

/*   fpAntenna = fopen(AntennaInfo.AntennaInfoFileName,"r"); */


/*   if (!fpAntenna) { */
/*     printf("Can not read antenna file : %s\n",AntennaInfo.AntennaInfoFileName); */
/*     iError = 1; */
/*     //newstructure */
/* // return; */
/* //newstructure */
/* return 0; */
/*   } */

/*   int IsDone = 0; */

/*   while (!IsDone) { */

/*     if (fgets(line,500,fpAntenna) == NULL) IsDone = 1; */

/*     if (strstr(line,"#N") != NULL) { */
/*       fscanf(fpAntenna, "%d", &AntennaInfo.nAnt); */
/*       //printf("%d\n",AntennaInfo.nAnt); */
/*     } */
/*     if (strstr(line,"#YAW") != NULL) { */
/*       for (i=0;i<AntennaInfo.nAnt;i++) { */
/* 	fscanf(fpAntenna, "%f", &AntennaInfo.yaws[i]); */
/* 	//printf("%d %f\n",i,AntennaInfo.yaws[i]); */
/*       } */
/*     } */
/*     if (strstr(line,"#ROLL") != NULL) { */
/*       for (i=0;i<AntennaInfo.nAnt;i++) { */
/* 	fscanf(fpAntenna, "%f", &AntennaInfo.rolls[i]); */
/* 	//printf("%d %f\n",i,AntennaInfo.rolls[i]); */
/*       } */
/*     } */
/*     if (strstr(line,"#PITCH") != NULL) { */
/*       for (i=0;i<AntennaInfo.nAnt;i++) { */
/* 	fscanf(fpAntenna, "%f", &AntennaInfo.pitches[i]); */
/* 	//printf("%d %f\n",i,AntennaInfo.pitches[i]); */
/*       } */
/*     } */

/*   } */

/*   close(fpAntenna); */

/*   return iError; */

/* } */


// sigma0
//  openw,1,file+'.bin'
//  writeu,1,nWinds,nAngles
//  for i=0,nWinds-1 do for j=0,nAngles-1 do writeu,1,Sigma0(i,j)
//  close,1

// antenna
//  openw,1,'antenna.bin'
//  writeu,1,nPhiN,nThetaN
//  for i=0,nPhiN-1 do for j=0,nThetaN-1 do $
//     writeu,1,phiN(i,j),thetaN(i,j),gainN(i,j)
//  close,1

//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

double radius_at_point(float lat) {

  double r;

  r = fabs(lat)/90.0 * (rpole-req) + req;
  
  return r;

}

//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

double c_i_to_r(int itime[7]) {

  int dayofmon[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
  double timereal;
  int nyear, nleap, nmonth, nday, nhour, nmin, nsec, nmilli, i;

  if (itime[0] % 4 == 0) dayofmon[1]++;

  timereal = 0.0;

  if (itime[0] > 1900) {
    nyear = itime[0] - 1965;
  } else {
    if (itime[0] > 65) {
      nyear = itime[0] - 65;
    } else {
      nyear = itime[0] + 100 - 65;
    }
  }

  nleap = nyear/4;

  nmonth = itime[1] - 1;

  nday = 0;

  for (i=0; i < nmonth; i++) nday = nday + dayofmon[i];

  nday = nday + itime[2] - 1;
  nhour  = itime[3];
  nmin   = itime[4];
  nsec   = itime[5];
  nmilli = itime[6];

  timereal =
    ((double) nmilli) / ((double) (1000.0)) +
    ((double) nsec) * ((double) (1.0)) +
    ((double) nmin) * ((double) (60.0)) +
    ((double) nhour) * ((double) (60.0*60.0)) +
    ((double) nday) * ((double) (24.0*60.0*60.0)) +
    ((double) nleap) * ((double) (24.0*60.0*60.0)) +
    ((double) nyear) * ((double) (365.0*24.0*60.0*60.0));
  
  return timereal;

}

//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

int compare(char a[], char b[]) {

   int c = 0;
   int la, lb;

   la = strlen(a);
   lb = strlen(b);

   if (la > 10) la = 10;
   if (lb < la) la = lb;

   // Only compare the first 10 (or less) characters
 
   while( a[c] == b[c] && c < la) {
      if( a[c] == '\0' || b[c] == '\0' )
         break;
      c++;
   }
   if (c == la) return 1;
   else return 0;

}

//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

int getint(char line[1000], int iStart, int nChars) {

  char linecopy[1000];
  int i, result;

  result = 0;

  for (i=0; i<nChars; i++) linecopy[i] = line[i+iStart];
  linecopy[i] = '\0';

  sscanf(linecopy, "%d", &result);

  return result;

}

//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

int getmonth(char line[1000], int iStart, int nChars) {

  char linecopy[1000];
  int i, result;

  result = 0;

  for (i=0; i<nChars; i++) linecopy[i] = line[i+iStart];
  linecopy[i] = '\0';

  if (compare(linecopy, "Jan")) result = 1;
  if (compare(linecopy, "Feb")) result = 2;
  if (compare(linecopy, "Mar")) result = 3;
  if (compare(linecopy, "Apr")) result = 4;
  if (compare(linecopy, "May")) result = 5;
  if (compare(linecopy, "Jun")) result = 6;
  if (compare(linecopy, "Jul")) result = 7;
  if (compare(linecopy, "Aug")) result = 8;
  if (compare(linecopy, "Sep")) result = 9;
  if (compare(linecopy, "Oct")) result = 10;
  if (compare(linecopy, "Nov")) result = 11;
  if (compare(linecopy, "Dec")) result = 12;

  return result;

}

//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

float getdouble(char line[1000], int nSkip) {

  char linecopy[1000];
  int i, iStart, nChars, nSpaces;
  double result;

  result = 0;

  iStart = 25;
  nSpaces = 0;
  nChars = 0;
  while (nSpaces < nSkip+1) {
    if (nSpaces == nSkip) {
      nChars++;
    } else {
      iStart++;
    }
    if (line[iStart+nChars] == 32 || line[iStart+nChars] == '\0') nSpaces++;

    //printf("%d %d %d %d ->%c<-\n",nSkip,nSpaces,iStart,nChars,line[iStart+nChars]);

  }

  //printf("line ::%s\n",line);

  //printf("istart nchars %d %d\n",iStart, nChars);

  for (i=0; i<nChars; i++) linecopy[i] = line[i+iStart];
  linecopy[i] = '\0';

  //printf("linecopy ::%s\n",linecopy);


  sscanf(linecopy, "%lf", &result);

  //printf("getdouble : %lf\n", result);

  return result;

}

//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

int strcompress(char line[600]) {

  char linecopy[600];
  int i, j, IsSpace, result;

  result = 0;

  j = 0;
  for (i=0; i<strlen(line); i++) {
    // turn commas into spaces:
    if (line[i] == 44) line[i] = 32;
    linecopy[j] = line[i];
    if (line[i] == 32) {
      IsSpace++;
    } else {
      IsSpace = 0;
    }
    if (IsSpace <= 1) j++;
  }
  linecopy[j] = '\0';

  strcpy(line, linecopy);

  return result;

}

//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

int read_stk(char filename[256], struct SatelliteDetails *sat) {

  FILE *fpSTK = NULL;
  char line[600], line2[600], testline[40], testline2[40];
  int IsFound, IsDone, IsDoneInner, IsPropagatorFile;
  int iSat, i, res, err;
  int iTime[7], nSats, nPts;

  double x, y, z, r, xy;

  err = 0;

  fpSTK = fopen(filename,"r");
        if (fpSTK == NULL){
	  printf("***! Could not find the file %s. The program will stop. !***\n", filename); MPI_Finalize();exit(0);
      }

  strcpy(testline,"Time (UTCG)");
  strcpy(testline2,"#START");
  IsPropagatorFile = 0;

  IsDone = 0;

  iSat = 0;

  nSats = 0;

  while (!IsDone) {

    if (fgets(line,500,fpSTK) == NULL) IsDone = 1;

    //printf("line-->%s<--\n",line);

    IsFound = 0;
    if (strstr(line,testline2) != NULL) {
      IsPropagatorFile = 1;
    }

    if (strstr(line,testline)  != NULL ||
	strstr(line,testline2) != NULL) {

      // files with commas in them don't have an extra line of space.

      if ( (strstr(line,",") == NULL) && (IsPropagatorFile == 0)) fgets(line,500,fpSTK);
      IsDoneInner = 0;

      nPts = 0;

      while (!IsDoneInner) {

	/* printf("here first!\n"); */

	/* printf("nPts : %d \n",nPts); */

	if (fgets(line,500,fpSTK) == NULL) IsDoneInner = 1;

	if (strlen(line) > 20) {

	  err = strcompress(line);
	  //	  printf("line-->%s<--\n",line);
	  //printf("nPts : %d \n",nPts);

	  // if the day is less than 10, then there is often only a
	  // single diget, which is bad.  So, fix this:
	  if (line[1] == 32) {
	    strcpy(line2,line);
	    strcpy(line," ");
	    strcat(line, line2);
	  }

	  if (IsPropagatorFile) {
	    iTime[0] = getint(line,0,4);
	    iTime[1] = getint(line,5,2);
	    iTime[2] = getint(line,8,2);
	    iTime[3] = getint(line,11,2);
	    iTime[4] = getint(line,14,2);
	    iTime[5] = getint(line,17,2);
	    iTime[6] = 0; //getint(line,20,3);
	    // the propagator files have a time-stamp that is 4 characters too short,
	    // so just add a space to the front....
	    strcpy(line2,line);
	    strcpy(line,"      ");
	    strcat(line, line2);
	  } else {
	    iTime[0] = getint(line,7,4);
	    iTime[1] = getmonth(line,3,3);
	    iTime[2] = getint(line,0,2);
	    iTime[3] = getint(line,12,2);
	    iTime[4] = getint(line,15,2);
	    iTime[5] = getint(line,18,2);
	    iTime[6] = getint(line,21,3);
	  }
	  // CBV 05/09/16
	  int sss;
	  for (sss = 0; sss < 7; sss++){
	    sat[nSats].time_ymdhmsm[nPts][sss] = iTime[sss];
	  }
	  // END CBV 05/09/16

	  sat[nSats].time[nPts] = c_i_to_r(iTime);

	  //printf("line-->%s<--\n",line);
	  x = getdouble(line,0)*1000.0;
	  y = getdouble(line,1)*1000.0;
	  z = getdouble(line,2)*1000.0;

	  //printf("xyz : %lf %lf %lf\n", x, y, z);

	  sat[nSats].x[nPts] = x;
	  sat[nSats].y[nPts] = y;
	  sat[nSats].z[nPts] = z;

	  r  = sqrt(x*x + y*y + z*z);
	  xy = sqrt(x*x + y*y);

	  sat[nSats].lat[nPts] = asin(z/r)/dtor;
	  sat[nSats].lon[nPts] = acos(x/xy)/dtor;
	  if (y < 0) sat[nSats].lon[nPts] = 360.0-sat[nSats].lon[nPts];
	  sat[nSats].radius[nPts] = r;

	  //printf("lat lon r : %lf %lf %lf\n", sat[nSats].lat[nPts], sat[nSats].lon[nPts], sat[nSats].radius[nPts]);

	  sat[nSats].vx[nPts] = getdouble(line,3);
	  //printf("here\n");
	  sat[nSats].vy[nPts] = getdouble(line,4);
	  //printf("here %d %d %s\n",nSats, nPts, line);
	  sat[nSats].vz[nPts] = getdouble(line,5);
	  //printf("here\n");

	  nPts++;

	  //printf("here last!\n");
	  
	} else {
	  IsDoneInner = 1;
	}
      }

      //      if (iProc == 0) printf("Done reading in sat %2d with %d points\n",nSats, nPts);
      sat[nSats].IsGood = 1;
      sat[nSats].nPts = nPts;
      nSats++;
    }

  }

  close(fpSTK);

  return err;

}

//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

double compute_diff(double theta, double r, double r0, double r1, double d, double delta) {

  double tau, a, e, beta, beta2, rat, sigma, diff;

  tau = theta - delta;

  a = sqrt(r*r + r0*r0 - 2 * r * r0 * cos(delta));
  e = sqrt(r*r + r1*r1 - 2 * r * r1 * cos(tau));

  beta  = pi - asin(r0/a * sin(delta));
  if (beta > twopi) beta -= twopi;
  beta2 = pi - asin(r1/e * sin(tau));
  if (beta2 > twopi) beta2 -= twopi;
  rat   = (a*a + e*e - d*d)/(2*a*e);
  if (rat > 1.0) rat = 1.0;
  sigma = acos(rat);
  diff  = twopi - sigma - 2*beta;

  //printf("compute: %lf %lf %lf %lf %lf %lf\n", theta, delta, tau, sigma, beta, diff);

  return diff;

}


//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

int specular_point(float CygLat, float CygLon, float CygR,
		   float GpsLat, float GpsLon, float GpsR,
		   int iMethod,
		   struct specular_point_info *info) {

  double r, cct, gct, max_angle;
  double CygLatr, CygLonr, GpsLatr, GpsLonr;
  double CygX, CygY, CygZ, GpsX, GpsY, GpsZ;
  double cosTheta, theta, delta, deltaOld, deltaNew;
  double CygR_2, GpsR_2, r_2, rguess, rguessold, rg_2;
  double d;
  double diff, oldDiff, maxDiff;
  int i, n, nMax, nOuter;
  double tau, a, e, beta, beta2, rat, sigma;

  double cl, gl, d1, d2, percent, clp, glp;
  double CygX_tp, CygY_tp, CygZ_tp;
  double GpsX_tp, GpsY_tp, GpsZ_tp;
  double x_sp, y_sp, z_sp;
  double lat_sp, lon_sp, r_sp, xy;

  double dlat, dr, dr_threshold, x;

  double diffLeft, diffCenter, diffRight;
  double deltaLeft, deltaCenter, deltaRight;

  float cX, cY, cZ, cXY, cRot1, cXn, cYn, cZn;
  float cXnn, cYnn, cZnn, cXZ, cRot2;

  float gX, gY, gZ;
  float gXn, gYn, gZn;
  float gXnn, gYnn, gZnn, gZY, t, p;

  float sX, sY, sZ;
  float sXn, sYn, sZn;
  float sXnn, sYnn, sZnn;

  int dummy;
  double timing;

  maxDiff = 0.000001;
  dr_threshold = 0.25;

  r = req;

  // cygnss - center of earth - tangent point
  cct = acos(r/CygR);
  // gps - center of earth - tangent point
  gct = acos(r/GpsR);

  max_angle = cct + gct;

  if (iDebug > 0) {
    printf("------------------------------------------------------\n");
    printf("%f %f %f\n", CygLat, CygLon, CygR);
    printf("%f %f %f\n", GpsLat, GpsLon, GpsR);
  }

  dlat = GpsLat-CygLat;

  CygLatr = CygLat*dtor;
  CygLonr = CygLon*dtor;

  GpsLatr = GpsLat*dtor;
  GpsLonr = GpsLon*dtor;

  CygX = cos(CygLatr)*cos(CygLonr);
  CygY = cos(CygLatr)*sin(CygLonr);
  CygZ = sin(CygLatr);

  GpsX = cos(GpsLatr)*cos(GpsLonr);
  GpsY = cos(GpsLatr)*sin(GpsLonr);
  GpsZ = sin(GpsLatr);

  cosTheta = CygX*GpsX + CygY*GpsY + CygZ*GpsZ;
  theta = acos(cosTheta);

  if (theta > max_angle*0.95) return -1;

  //timing = MPI_Wtime();
  //for (dummy=0;dummy<1000000;dummy++) {

  CygR_2 = CygR*CygR;
  GpsR_2 = GpsR*GpsR;

  d = sqrt(CygR_2 + GpsR_2 - 2*CygR*GpsR*cosTheta);

  // totally a guess at the latitude of the specular point
  // just to get a guess at the radius of the Earth
  lat_sp = (GpsLat-CygLat)*(0.08+0.2*theta/pi) + CygLat;
  rguess = radius_at_point(lat_sp);

  delta = 0.0;

  dr = 100.0;

  nMax = 100;

  if (iMethod == 0) {

    deltaLeft   = 0.0;
    deltaRight  = theta * 0.3;
    deltaCenter = (deltaLeft+deltaRight)/2.0;

  } else {

    delta = -1.0;

    i = 2;
    while (CourseGrid[i].theta < theta && CourseGrid[i].delta > 0.0 && i < nPtsGrid-1) i++;

    x = (CourseGrid[i].theta - theta) / (CourseGrid[i].theta - CourseGrid[i-1].theta);
    deltaLeft   = (1.0 - x) * CourseGrid[i-1].delta + x * CourseGrid[i-2].delta;
    deltaCenter = (1.0 - x) * CourseGrid[i  ].delta + x * CourseGrid[i-1].delta;
    deltaRight  = (1.0 - x) * CourseGrid[i+1].delta + x * CourseGrid[i  ].delta;

  }

  nOuter = 0;

  while (dr > dr_threshold) {

    rguessold = rguess;

    rg_2  = rguess*rguess;

    diff = 1.0e32;

    n = 0;

    if (nOuter > 0) {
      deltaCenter = delta;
      deltaLeft   = delta - (deltaRight-deltaLeft)*2.0;
      deltaRight  = delta + (delta - deltaLeft);
    }

    nOuter++;

    //printf("delta : %lf %lf %lf %lf\n", delta, deltaLeft, deltaCenter, deltaRight);

    while (fabs(diff) > maxDiff && n < nMax) {

      if (iDebug > 0) printf("diff : %lf\n",diff);

      diffLeft   = compute_diff(theta, rguess, CygR, GpsR, d, deltaLeft);
      diffCenter = compute_diff(theta, rguess, CygR, GpsR, d, deltaCenter);
      diffRight  = compute_diff(theta, rguess, CygR, GpsR, d, deltaRight);

      diff = fabs(diffLeft);
      if (fabs(diffCenter) < diff) diff = fabs(diffCenter);
      if (fabs(diffRight)  < diff) diff = fabs(diffRight);
    
      if (iMethod == 1) {
	maxDiff = deltaCenter;
	dr_threshold = 1.0e32;
      }

      if (fabs(diffLeft) < maxDiff) {
	delta = deltaLeft;
	diff  = diffLeft;
      }

      if (fabs(diffRight) < maxDiff) {
	delta = deltaRight;
	diff  = diffRight;
      }

      if (fabs(diffCenter) < maxDiff) {
	delta = deltaCenter;
	diff  = diffCenter;
      }

      //printf("delta : %lf %lf %lf %lf %lf\n", delta, diff, deltaLeft, deltaCenter, deltaRight);

      // Use bifrication method to determine point.

      if (fabs(diff) > maxDiff) {
	if (diffLeft*diffCenter < 0.0) {
	  deltaRight = deltaCenter;
	} else {
	  if (diffRight*diffCenter < 0.0) {
	    deltaLeft = deltaCenter;
	  } else {
	    if (deltaRight > 0.0 && deltaLeft > 0.0) {
	      deltaRight = deltaLeft;
	      deltaLeft = 0.0;
	    } else {
	      deltaRight = deltaRight*1.1;
	    }
	  }
	}
	deltaCenter = (deltaLeft + deltaRight)/2.0;
      }

      n++;
      //printf("diff %d %lf %lf %lf\n",n,diffLeft, diffCenter, diffRight);

    }

    //printf("n : %d %d %lf\n",nOuter, n, delta);

    tau = theta - delta;

    a = sqrt(rg_2 + CygR_2 - 2 * rguess * CygR * cos(delta));
    e = sqrt(rg_2 + GpsR_2 - 2 * rguess * GpsR * cos(tau));

    beta  = pi - asin(CygR/a * sin(delta));
    if (beta > twopi) beta -= twopi;
    beta2 = pi - asin(GpsR/e * sin(tau));
    if (beta2 > twopi) beta2 -= twopi;
    rat   = (a*a + e*e - d*d)/(2*a*e);
    if (fabs(rat) < 1.0) {
      sigma = acos(rat);
    } else {
      if (rat <= -1.0) {
	sigma = -pi;
      } else {
	sigma = 0.0;
      }
    }

    if (iDebug > 0) printf("tau : %lf %lf %lf %lf\n", tau, beta, beta2, rat);

    diff  = twopi - sigma - 2*beta;

    // we need to make sure that we have the correct radius of the Earth,
    // so we calculate the longitude quickly, and get the radius:
    lat_sp = dlat * delta/theta + CygLat;
    rguess = radius_at_point(lat_sp);
    dr = fabs(rguess-rguessold);

    //printf("lat_sp : %lf\n",lat_sp);

    //if (dr >= dr_threshold) delta = 0.99*delta;

  }

//   }
//timing = MPI_Wtime() - timing;
//printf("timing : %lf\n",timing);

  // In order to get the longitude of the specular point, it is a bit more
  // math.....

  // distance along CoE-Cygnss line to tangent point through SP
  cl = rguess / cos(delta);
  // distance along CoE-GPS line to tangent point through SP
  gl = rguess / cos(tau);

  //printf("rguess: %lf %lf %lf\n", rguess, cl/r0, gl/GpsR);

  // distance from Cygnss tangent point to SP
  d1 = cl * sin(delta);
  // distance from Gps tangent point to SP
  d2 = gl * sin(tau);

  //printf("d1 : %lf %lf\n", d1, sqrt(cl*cl - d1*d1));

  percent = d1/(d1+d2);

  //printf("cl: %f %f %f\n", cl/1000.0, gl/1000.0, percent);

  clp = cl/CygR;
  glp = gl/GpsR;

  CygX_tp = CygR * CygX * clp;
  CygY_tp = CygR * CygY * clp;
  CygZ_tp = CygR * CygZ * clp;

  GpsX_tp = GpsR * GpsX * glp;
  GpsY_tp = GpsR * GpsY * glp;
  GpsZ_tp = GpsR * GpsZ * glp;

//  cX = r0 * x0;
//  cY = r0 * y0;
//  cZ = r0 * z0;
//
//  gX = x1;
//  gY = y1;
//  gZ = z1;
//
//  cXY = sqrt(cX*cX + cY*cY);
//  cRot1 = acos(cX/cXY);
//  if (cY<0) cRot1 = twopi - cRot1;
//
//  cRot1 = -cRot1;
//
//  cXn = cX*cos(cRot1) - cY*sin(cRot1);
//  cYn = cX*sin(cRot1) + cY*cos(cRot1);
//  cZn = cZ;
//
//  gXn = gX*cos(cRot1) - gY*sin(cRot1);
//  gYn = gX*sin(cRot1) + gY*cos(cRot1);
//  gZn = gZ;
//
//  cXZ = sqrt(cXn*cXn + cZn*cZn);
//  cRot2 = acos(cXn/cXZ);
//  if (cZn<0) cRot2 = twopi - cRot2;
//
//  cRot2 = -cRot2;
//
//  cXnn = cXn*cos(cRot2) - cZn*sin(cRot2);
//  cYnn = cYn;
//  cZnn = cXn*sin(cRot2) + cZn*cos(cRot2);
//
//  gXnn = gXn*cos(cRot2) - gZn*sin(cRot2);
//  gYnn = gYn;
//  gZnn = gXn*sin(cRot2) + gZn*cos(cRot2);
//
//  //printf("gxnn : %lf %lf %lf\n", gXnn, gYnn, gZnn);
//
//  r = sqrt(gXnn*gXnn + gYnn*gYnn + gZnn*gZnn);
//  gZY = sqrt(gZnn*gZnn + gYnn*gYnn);
//  t = acos(gXnn / r);
//  p = acos(gZnn/gZY);
//  if (gYnn < 0.0) p = twopi-p;
//
//  //printf("t : %lf %lf %lf\n", t/dtor, theta/dtor, p/dtor);
//
//  // Specular point is on the surface of the Earth....
//
//  sX = rguess;
//  sY = 0.0;
//  sZ = 0.0;
//
//  // Rotated by an angle of delta from CYGNSS....
//
//  sXn = sX * cos(delta) - sZ * sin(delta);
//  sYn = sY;
//  sZn = sX * sin(delta) + sZ * cos(delta);
//
//  //printf("delta : %lf %lf\n", delta, theta);
//
//  //sXn = sX * cos(theta) - sZ * sin(theta);
//  //sYn = sY;
//  //sZn = sX * sin(theta) + sZ * cos(theta);
//
//  // And about an angle p from CYGNSS....
//
//  sXnn = sXn;
//  sYnn = sZn * sin(p) + sYn * cos(p);
//  sZnn = sZn * cos(p) - sYn * sin(p);
//
//  //printf("sxnn : %lf %lf %lf\n", sXnn, sYnn, sZnn);
//
//  // now we have to put the coordinates back to the real frame....
//
//  cRot2 = -cRot2;
//
//  sXn = sXnn*cos(cRot2) - sZnn*sin(cRot2);
//  sYn = sYnn;
//  sZn = sXnn*sin(cRot2) + sZnn*cos(cRot2);
//
//  cRot1 = -cRot1;
//
//  sX = sXn*cos(cRot1) - sYn*sin(cRot1);
//  sY = sXn*sin(cRot1) + sYn*cos(cRot1);
//  sZ = sZn;
//
//  r_sp = sqrt(sX*sX + sY*sY + sZ*sZ);
//  lat_sp = asin(sZ/r_sp);

  //printf("cx: %lf %lf %lf\n", x0*r0/1000.0, y0*r0/1000.0, z0*r0/1000.0);
  //printf("gx: %lf %lf %lf\n", x1*GpsR/1000.0, y1*GpsR/1000.0, z1*GpsR/1000.0);
  //printf("sp: %lf %lf %lf\n", sX/1000.0, sY/1000.0, sZ/1000.0);
  //printf("lat_sp : %lf %lf\n", lat_sp/dtor, r_sp);

  //printf("x0: %lf %lf %lf\n", x0_tp/1000.0, y0_tp/1000.0, z0_tp/1000.0);
  //printf("x1: %lf %lf %lf\n", x1_tp/1000.0, y1_tp/1000.0, z1_tp/1000.0);

  info->DistToTangentCygnss = cl;
  info->DistToTangentGps    = gl;
  info->PercentAlongTangent = percent;

  x_sp = CygX_tp + (GpsX_tp - CygX_tp)*percent;
  y_sp = CygY_tp + (GpsY_tp - CygY_tp)*percent;
  z_sp = CygZ_tp + (GpsZ_tp - CygZ_tp)*percent;

  r_sp = sqrt(x_sp*x_sp + y_sp*y_sp + z_sp*z_sp);
  lat_sp = asin(z_sp/r_sp);
  xy = sqrt(x_sp*x_sp + y_sp*y_sp);
  if (xy > 0) {
    lon_sp = acos(x_sp/xy);
    if (y_sp < 0) lon_sp = twopi - lon_sp;
  } else {
    lon_sp = 0.0;
  }

  info->delta     = delta; // Cygnss-Center of Earth (CoE)-Specular Point (SP)
  info->tau       = tau; // GPS - CoE - SP
  info->beta      = beta; // Cygnss - SP - CoE
  info->sigma     = sigma; // Cygnss - SP - GPS
  info->Incidence = (sigma/2)/dtor;
  if (iDebug > 0) printf("sigma : %lf %lf %lf\n", sigma, (sigma/2)/dtor, theta/dtor);
  info->theta     = theta; // Cygnss - CoE - GPS
  info->a         = a; // Cygnss - SP
  info->e         = e; // GPS - SP
  info->d         = d; // Cygnss - GPS
  info->diff      = diff; // Error in calculation
  info->latitude  = lat_sp/dtor;
  info->longitude = lon_sp/dtor;
  info->radius    = radius_at_point(lat_sp);

  info->CygX = CygR*CygX;
  info->CygY = CygR*CygY;
  info->CygZ = CygR*CygZ;

  info->GpsX = GpsR*GpsX;
  info->GpsY = GpsR*GpsY;
  info->GpsZ = GpsR*GpsZ;

  info->SpX = x_sp;
  info->SpY = y_sp;
  info->SpZ = z_sp;

  info->nIters = n;

  return 0;

}

//-----------------------------------------------------------------------
//
//-----------------------------------------------------------------------

int find_gain(float CygX,  float CygY,  float CygZ,
	      float CygVx, float CygVy, float CygVz,
	      float GpsX,  float GpsY,  float GpsZ,
	      float GpsVx, float GpsVy, float GpsVz,
	      float SpX,   float SpY,   float SpZ,
	      float AntennaTilt, float RotationAngle, float yawoverall,  char gps_name[200], int printVar, int iGps, int iPt, int iPtInner, int nPtsInner, ATTITUDE_T ATTITUDE,

	      struct specular_point_info *info) {
/*   // !!!!!!!!!! remove below */
/*   SpX = -857.732*1000.; */
/*   SpY = 5616.403*1000.; */
/*   SpZ = -2888.636*1000.; */
/*   // !!!!!!!!!! end of remove below */


  double C2Sx, C2Sy, C2Sz;
  double S2Gx, S2Gy, S2Gz;
  double SPscX, SPscY, SPscZ;
  double stosXr, stosYr, stosZr;
  double stosXrm, stosYrm, stosZrm;
  double CygLos, GpsLos;

  double C2Smag, S2Gmag;

  //double x_sp, y_sp, z_sp;

  double xi, xj, xk;
  double yi, yj, yk;
  double zi, zj, zk;
  double mag;

  double t, t2, ct, st, theta_sat, theta_satm, xy, phi_sat, phi_satm;

  int ii, jj;

  double fx, fy, gainR, gainL;

  if (iDebug > 0) printf("entering find_gain\n");

  mag = sqrt(CygVx*CygVx + CygVy*CygVy + CygVz*CygVz);
  xi = CygVx/mag;
  xj = CygVy/mag;
  xk = CygVz/mag;

  mag = sqrt(CygX*CygX + CygY*CygY + CygZ*CygZ);
  zi = -CygX/mag;
  zj = -CygY/mag;
  zk = -CygZ/mag;

  yi = -(xj*zk - xk*zj);
  yj =  (xi*zk - xk*zi);
  yk = -(xi*zj - xj*zi);
  // Y should be magnitude 1.0, but just in case....
  mag = sqrt(yi*yi + yj*yj + yk*yk);
  if (mag == 0) {
    printf("position and velocity of satellite are in the same line!\n");
    printf("This make the Y-axis undefined! Must stop!\n");
    //newstructure
// return;
//newstructure
return 0;
  }
  yi = yi/mag;
  yj = yj/mag;
  yk = yk/mag;

  C2Sx = -(CygX - SpX);
  C2Sy = -(CygY - SpY);
  C2Sz = -(CygZ - SpZ);
  C2Smag = sqrt(C2Sx*C2Sx + C2Sy*C2Sy + C2Sz*C2Sz);


  info->CygnssToSpecularPoint = C2Smag;
  // cbv 04-27-18
  double C2S[3];
  C2S[0] = C2Sx ;  C2S[1] = C2Sy ;  C2S[2] = C2Sz ;

  double cygnss_r[3];
  double cygnss_v[3];
  cygnss_r[0] = CygX;  cygnss_r[1] = CygY;  cygnss_r[2] = CygZ;
  cygnss_v[0] = CygVx;  cygnss_v[1] = CygVy;  cygnss_v[2] = CygVz;

  double T_ecef_2_ntw[3][3];
  compute_T_ecef_2_ntw_not_like_prop_math( T_ecef_2_ntw,
			 cygnss_r,
			 cygnss_v);

 double C2S_orbit[3];
 m_x_v(C2S_orbit, T_ecef_2_ntw, C2S);

 int index_attitude_interpo;
 index_attitude_interpo = iPt * nPtsInner + iPtInner;
 double T_sc_to_ntw[3][3];
 double v_angle[3]; // 0 pitch, 1 roll, 2 yaw. in degrees
 int order_rotation[3]; // 0 pitch, 1 roll, 2 yaw. equal to 1 means first rotation, 2 means second rtoation, 3 means third rotation. ex: order_rotation[0] = 2, order_rotation[1] = 1, order_rotation[2] = 3 means roll then pitch then yaw
 order_rotation[0] = ATTITUDE.order_pitch[index_attitude_interpo]; order_rotation[1] = ATTITUDE.order_roll[index_attitude_interpo]; order_rotation[2] = ATTITUDE.order_yaw[index_attitude_interpo];  // in  UM Doc. 148-0336, CYGNSS Science Data Processing Coordinate Systems D -> by default pitch then roll then yaw
 v_angle[0] = ATTITUDE.pitch[index_attitude_interpo]; v_angle[1] = ATTITUDE.roll[index_attitude_interpo]; v_angle[2] = ATTITUDE.yaw[index_attitude_interpo];
/*  etprint(ATTITUDE.et_interpo[index_attitude_interpo], "time"); */
/*  printf("<%s> %d %d\n", gps_name, iPt, iPtInner); */
/*  printf("i = %d\n", index_attitude_interpo); */
/*  printf("%f %f %f\n", v_angle[0], v_angle[1], v_angle[2]); */
 // printf("%d %d %d\n", order_rotation[0], order_rotation[1], order_rotation[2]);
 
 compute_T_sc_to_ntw(T_sc_to_ntw, v_angle, order_rotation);
 double T_ntw_to_sc[3][3];
 m_trans(T_ntw_to_sc, T_sc_to_ntw);
 double C2S_body[3];
 m_x_v(C2S_body, T_ntw_to_sc, C2S_orbit);
 


 double C2S_body_xy[3];
 C2S_body_xy[0] = C2S_body[0];
 C2S_body_xy[1] = C2S_body[1];
 C2S_body_xy[2] = 0;

 double C2S_body_xy_norm[3];
 v_norm(C2S_body_xy_norm, C2S_body_xy);


  if (C2S_body_xy_norm[1] < 0){
    info->azim = 360. - acos(C2S_body_xy_norm[0])/dtor;
  }
  else{
    info->azim = acos(C2S_body_xy_norm[0])/dtor;
  }
    double C2S_body_norm[3];
    v_norm(C2S_body_norm, C2S_body);
    info->elev = acos(C2S_body_norm[2])/dtor;

/*     if (printVar == 1){ */
/*           if ( strcmp(gps_name, "PRN_16") == 0){ */
/* 	    v_print(C2S_body, "C2S_body"); */
/* 	    printf("info->azim %f, elev %f\n", info->azim, info->elev); */
/* 	  } */
/*     } */


  // end of cbv 04-27-18

  C2Sx = C2Sx/C2Smag;
  C2Sy = C2Sy/C2Smag;
  C2Sz = C2Sz/C2Smag;




  //  if (iDebug > 0) {
  
/*     printf("c: %f %f %f\n",CygX/1000.0, CygY/1000.0, CygZ/1000.0); */
/*     printf("g: %f %f %f\n",GpsX/1000.0, GpsY/1000.0, GpsZ/1000.0); */
/*     printf("s: %f %f %f\n",SpX/1000.0, SpY/1000.0, SpZ/1000.0); */
/*     printf("c2s: %f %f %f\n",C2Sx, C2Sy, C2Sz); */
    //  }
    //    MPI_Finalize();exit(0);
  S2Gx = GpsX - SpX;
  S2Gy = GpsY - SpY;
  S2Gz = GpsZ - SpZ;
  S2Gmag = sqrt(S2Gx*S2Gx + S2Gy*S2Gy + S2Gz*S2Gz);
  S2Gx = S2Gx/S2Gmag;
  S2Gy = S2Gy/S2Gmag;
  S2Gz = S2Gz/S2Gmag;
  info->SpecularPointToGPS = S2Gmag;

  if (iDebug > 0) printf("find_gain: before los\n");

  CygLos = C2Sx*CygVx + C2Sy*CygVy + C2Sz*CygVz;
  GpsLos = S2Gx*GpsVx + S2Gy*GpsVy + S2Gz*GpsVz;

  info->CygnssLOSVel = CygLos;
  info->GpsLOSVel    = GpsLos;


  SPscX = -(C2Sx*xi + C2Sy*xj + C2Sz*xk);
  SPscY = -(C2Sx*yi + C2Sy*yj + C2Sz*yk);
  SPscZ = -(C2Sx*zi + C2Sy*zj + C2Sz*zk);

  mag = sqrt(SPscX*SPscX + SPscY*SPscY + SPscZ*SPscZ);
  SPscX = SPscX/mag;
  SPscY = SPscY/mag;
  SPscZ = SPscZ/mag;

  if (iDebug > 0) printf("find_gain: before rotations\n");

  // This is the position in the tilted coordinate system of the antenna

  float yaw, cy, sy;
  float pitch, cp, sp;
  float tmpX, tmpY, tmpZ;
  int iAnt;

  float GainHighest=-1.0e32;


  for (iAnt=0; iAnt < AntennaInfo.nAnt; iAnt++) {


    yaw = (AntennaInfo.yaws[iAnt] + yawoverall)*dtor;
    cy = cos(yaw);
    sy = sin(yaw);
    //    printf("yaw %f %f\n", AntennaInfo.yaws[iAnt], yawoverall);

    tmpX = SPscX * cy - SPscY * sy;
    tmpY = SPscY * cy + SPscX * sy;
    tmpZ = SPscZ;

    t = 0;//AntennaInfo.rolls[iAnt]*dtor;
    //    printf("RotationAngle[%d] %e %f\n", iAnt, AntennaInfo.rolls[iAnt], RotationAngle);

    t2 = RotationAngle*dtor;
    ct = cos(t+t2);
    st = sin(t+t2);

    stosXr =  tmpX;
    stosYr =  tmpY * ct - tmpZ * st;
    stosZr =  tmpZ * ct + tmpY * st;

    pitch = AntennaInfo.pitches[iAnt]*dtor;
    cp = cos(pitch);
    sp = sin(pitch);

    tmpX = stosXr * cp - stosZr * sp;
    tmpY = stosYr;
    tmpZ = stosZr * cp + stosXr * sp;

    stosXr = tmpX;
    stosYr = tmpY;
    stosZr = tmpZ;



    if (iDebug > 0) printf("find_gain: before theta_sat\n");

    theta_sat  = info->elev;//  cbv  modfiied, used to be acos(stosZr )/dtor;
    //theta_satm = acos(stosZrm)/dtor;

    xy = sqrt(stosXr*stosXr + stosYr*stosYr) + 1.0e-7;

    phi_sat = info->azim; // cbv modfified, used to be     phi_sat = acos(stosXr/xy)/dtor + 90.0;

/*     if (printVar == 1){ */
/*           if ( strcmp(gps_name, "PRN_16") == 0){ */
/* 	    if (iAnt == 0){ */
/* /\* 	      printf("%f %f %f %f %f\n", acos(stosXr)/dtor, acos(stosXr)/dtor + 90.0, acos(stosXr/xy)/dtor, acos(stosXr/xy)/dtor + 90.0, azim); *\/ */
/* /\*     printf("c: %f %f %f\n",CygX/1000.0, CygY/1000.0, CygZ/1000.0); *\/ */
/* /\*     printf("g: %f %f %f\n",GpsX/1000.0, GpsY/1000.0, GpsZ/1000.0); *\/ */
/* /\*     printf("s: %f %f %f\n",SpX/1000.0, SpY/1000.0, SpZ/1000.0); *\/ */
/* /\*     printf("c2s: %f %f %f\n",C2Sx, C2Sy, C2Sz); *\/ */
/*     //    printf("vDotC2S %f, crossDotC2S %f | %f\n", vDotC2S, crossDotC2S, acos(vDotC2S)/dtor); */
    
/* 	    } */

/*     } */
/*     } */

// cbv commented:    if (phi_sat > 360.0) phi_sat -= 360.0;

    xy = sqrt(stosXrm*stosXrm + stosYrm*stosYrm) + 1.0e-7;
    //    printf("<%s>\n", gps_name);

    phi_satm = acos(stosXrm/xy)/dtor; //+ 90.0;
    if (stosYrm < 0.0) phi_satm = 360.0-phi_satm;
    phi_satm += 90.0;
    if (phi_satm > 360.0) phi_satm -= 360.0;

    
    if (iDebug > 0) printf("find_gain: before gain lookup\n");


    ii = (int)(phi_sat / AntennaInfo.dPhi);

    jj = (int)(theta_sat / AntennaInfo.dTheta);


    fx = (  phi_sat - AntennaInfo.phi[ii])  /AntennaInfo.dPhi;
    fy = (theta_sat - AntennaInfo.theta[jj])/AntennaInfo.dTheta;
    



/*     info->elev = theta_sat; */
/*     info->azim = phi_sat; */


    // cbv new antenna onboard. before: block below this block
                  gainR = AntennaInfo.gain[iAnt][ii  ][jj  ] ;
    
/* 	      gainR = (1-fx)*(1-fy) * AntennaInfo.gain[iAnt][ii  ][jj  ] + */
/*             (  fx)*(1-fy) * AntennaInfo.gain[iAnt][ii+1][jj  ] + */
/*             (1-fx)*(  fy) * AntennaInfo.gain[iAnt][ii  ][jj+1] + */
/*             (  fx)*(  fy) * AntennaInfo.gain[iAnt][ii+1][jj+1]; */

/* 	      gainR = round(gainR); */

    if (printVar == 1){
    if ( strcmp(gps_name, "PRN_29") == 0){
      printf("%f %f | %f %f | %f || %f %f %f %f | %d %d\n", theta_sat, AntennaInfo.theta[jj], phi_sat, AntennaInfo.phi[ii] ,gainR, AntennaInfo.gain[iAnt][ii  ][jj  ], AntennaInfo.gain[iAnt][ii+1][jj  ] ,  AntennaInfo.gain[iAnt][ii  ][jj+1],  AntennaInfo.gain[iAnt][ii+1][jj+1], ii, jj);


    }
    }

/* 	printf("\n iAnt %d\n", iAnt); */
/* 	printf("phi %f %f %d %f\n", phi_sat, AntennaInfo.phi[ii], ii, AntennaInfo.dPhi); */
/* 	printf("theta %f %f %d %f\n", theta_sat, AntennaInfo.theta[jj], jj, AntennaInfo.dTheta); */
/* 	printf("gain %f %f\n", AntennaInfo.gain[iAnt][ii  ][jj  ], gainR); */
    //    printf("%d %f %f %d %d\n", iAnt, gainR, AntennaInfo.gain[iAnt][ii  ][jj  ], ii, jj);
    // cbv new antenna onboard. before: block below
/*     gainR = (1-fx)*(1-fy) * AntennaInfo.gain[ii  ][jj  ] + */
/*             (  fx)*(1-fy) * AntennaInfo.gain[ii+1][jj  ] + */
/*             (1-fx)*(  fy) * AntennaInfo.gain[ii  ][jj+1] + */
/*             (  fx)*(  fy) * AntennaInfo.gain[ii+1][jj+1]; */
    // end of cbv new antenna onboard. before: block below


    if (gainR > GainHighest) GainHighest = gainR;

    
  }
  //  printf("\n\n") ;

  //	MPI_Finalize(); exit(0);
//  ii = (int) (phi_satm / AntennaInfo.dPhi);
//  jj = (int) (theta_satm / AntennaInfo.dTheta);
//
//  fx = (  phi_satm - AntennaInfo.phi[ii])  /AntennaInfo.dPhi;
//  fy = (theta_satm - AntennaInfo.theta[jj])/AntennaInfo.dTheta;
//
//  gainL = (1-fx)*(1-fy) * AntennaInfo.gain[ii  ][jj  ] +
//          (  fx)*(1-fy) * AntennaInfo.gain[ii+1][jj  ] +
//          (1-fx)*(  fy) * AntennaInfo.gain[ii  ][jj+1] +
//          (  fx)*(  fy) * AntennaInfo.gain[ii+1][jj+1];

  if (iDebug > 0) printf("find_gain: before filling in info->gain\n");

  info->gain = GainHighest;
  //if (gainL > gainR) info->gain = gainL;

  if (iDebug > 0) printf("%f %f %f %f\n",info->gain, AntennaInfo.factor, C2Smag, S2Gmag);

  // cbv new antenna onboard. before: before was the block below this one
  //  info->NormPower =  info->gain ;
    
  // cbv new antenna onboard. before: block below
  info->NormPower =
    info->gain;// *
  /*   AntennaInfo.factor / */
/*     (C2Smag*C2Smag) / (S2Gmag*S2Gmag); */
  // end of cbv new antenna onboard. before: block below
  if (iDebug > 0) printf("gain : %f %e\n",info->gain,info->NormPower);

  if (iDebug > 0) printf("%f %f\n",info->Incidence, AntennaInfo.dIncidence);

  // determine the closest Sigma0 to get given the incidence angle:
  jj = (int) (info->Incidence/AntennaInfo.dIncidence + 0.5);

  if (iDebug > 0) printf("find_gain jj = %d\n",jj);

  // Estimate the power given a 10 m/s wind:
  info->power = info->NormPower * AntennaInfo.Sigma0[9][jj];

  if (iDebug > 0) printf("exiting find_gain\n");

}


//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

int fill_specular_point_grid(double rCygnss, double rGps) {

  int i, err;
  double cct, gct, MaxAngle, angle;
  struct specular_point_info specularpoint;

  double cLat, cLon, cR, gLat, gLon, gR;

  // cygnss - center of earth - tangent point
  cct = acos(re/rCygnss)/dtor;
  // gps - center of earth - tangent point
  gct = acos(re/rGps)/dtor;

  MaxAngle = cct + gct;

  for (i = 0; i < nPtsGrid; i++) {

    cLat = 0.0;
    cLon = 0.0;
    cR = rCygnss;

    gLat = 0.0;
    gLon = MaxAngle/nPtsGrid * (i+1);
    gR = rGps;

    err = specular_point(cLat, cLon, cR,
			 gLat, gLon, gR,
			 0,
			 &specularpoint);

    CourseGrid[i].DistToTangentCygnss = specularpoint.DistToTangentCygnss;
    CourseGrid[i].DistToTangentGps    = specularpoint.DistToTangentGps;
    CourseGrid[i].PercentAlongTangent = specularpoint.PercentAlongTangent;
    CourseGrid[i].delta               = specularpoint.delta;
    CourseGrid[i].theta               = specularpoint.theta;

  }

  return 0;

}


//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

int find_specular_points( struct SatelliteDetails GPS[nGps],
			  struct SatelliteDetails CYGNSS[nCygnss],
			  char filehead[],
			  int lonshift, float rotate, float yawoverall, int iUseMin,
			  int LimitToOneDay,
			  int nGps_prop, int nCygnss_prop, ATTITUDE_T ATTITUDE) { // CBV 05/08/16

  // On my computer, I can't seem to allocate enough memory to actually
  // store everything.  So, I will just write it out to a file.
  double gps_x_inner, gps_y_inner, gps_z_inner, gps_vx_inner, gps_vy_inner, gps_vz_inner;
  double cygnss_x_inner, cygnss_y_inner, cygnss_z_inner, cygnss_vx_inner, cygnss_vy_inner, cygnss_vz_inner;

  FILE *fpOutput;
  char filename[400], sCygnss[3];
  int iCygnss, iGps, iPt, nPts, err, nPtsInner, iPtInner, iPtTotal;
  int iStart, iEnd, nCygEachPe, nCygLeft, i, iYawPoint, iYawPointBefore;
  int nVars;

  double cLat, cLon, cR;
  double gLat, gLon, gR;
  double cVx, cVy, cVz;
  double gVx, gVy, gVz;

  double cDt, gDt;

  double x, lon1, lon2;

  float tmp;

  float LonShift;

  struct specular_point_info specularpoint;

  struct specular_point_info SpecularPointsAllGps[nGps];
  double NormPowerAllGps[nGps];

  float gLonAllGps[nGps];
  float gLatAllGps[nGps];
  float gRAllGps[nGps];

  float MaxNormPower;
  int iMaxNormPower;
  int iOrder[nGps];

  int iNumber, iIncidence, iWind;

  float power, ratio;
  char iRatio[nWinds];

  LonShift = (float) lonshift;

  //printf("rotate: %f\n",rotate);

  nCygEachPe = nCygnss_prop/nProcs;

  nCygLeft = nCygnss_prop - (nCygEachPe * nProcs);

  iStart = 0;
  for (i=0; i<iProc; i++) {
    iStart += nCygEachPe;
    if (i < nCygLeft && iProc > 0) iStart++;
  }

  iEnd = iStart+nCygEachPe;
  if (iProc < nCygLeft) iEnd++;

  
  for (iCygnss=iStart; iCygnss<iEnd; iCygnss++) {


    // CBV 05/08/16
    /* strcpy(filename,filehead); */
    /* strcat(filename,".cyg"); */
    /* strcat(filename, c_int_str(iCygnss,2)); */
    /* if (lonshift > 0) { */
    /*   strcat(filename,".lon"); */
    /*   strcat(filename, c_int_str(lonshift,3)); */
    /* } */
    strcpy(filename,filehead);
    strcat(filename, CYGNSS[iCygnss].name);
    strcat(filename, "/specular_");
    strcat(filename, CYGNSS[iCygnss].name);

    // END CBV 05/08/16

    strcat(filename,".bin");
    fpOutput = fopen(filename,"w");

    if (iUseMin > 0) {
      nVars = 5;
    } else {
      nVars = 15;
    }

    //    if (iCygnss == 0) fwrite(&nVars,sizeof(nVars),1,fpOutput); // CBV 05/09/16

    //        if (iProc == 0) printf("iCygnss : %d %d\n",iCygnss, CYGNSS[iCygnss].IsGood);


    if (CYGNSS[iCygnss].IsGood == 1) {

      cDt = CYGNSS[iCygnss].time[1] - CYGNSS[iCygnss].time[0];
      nPtsInner = (int) cDt;

      //      nPts = (CYGNSS[iCygnss].nPts/1400) * 1440;
      //      if (nPts > 1440 && LimitToOneDay) nPts = 1440;
      nPts = CYGNSS[iCygnss].nPts-1;

      //printf("nPts : %d\n",nPts);

      for (iPt=0; iPt < nPts; iPt++) { //      for (iPt=1; iPt < nPts; iPt++) { // CBV 05/09/16

	if (iCygnss == 3){
		  if (iPt == 2 ){
	  int sss;
	  for (sss = 0; sss < 7; sss++){
	    printf("%d ", CYGNSS[iCygnss].time_ymdhmsm[iPt][sss]);
	  }
	  printf("\n");
	}
	}


	//		if (iProc == 0) printf("iCygnss: %d (of %d) iPt: %d (of %d)\n",iCygnss+1, iEnd, iPt, nPts);

	for (iPtInner=0; iPtInner < nPtsInner; iPtInner++) {

	  iPtTotal = iPt*nPtsInner + iPtInner; //iPtTotal = (iPt-1)*nPtsInner + iPtInner; // CBV 05/09/16

	  x = ((double) iPtInner) / (double) nPtsInner;
	  cLat = x * CYGNSS[iCygnss].lat[iPt+1] + (1-x) * CYGNSS[iCygnss].lat[iPt];
	  cR   = x * CYGNSS[iCygnss].radius[iPt+1] + (1-x) * CYGNSS[iCygnss].radius[iPt];
	  
	  lon1 = CYGNSS[iCygnss].lon[iPt+1];

	  lon2 = CYGNSS[iCygnss].lon[iPt];
	  if (lon1 > 300.0 && lon2 < 60.0) lon1 = lon1-360.0;
	  if (lon2 > 300.0 && lon1 < 60.0) lon2 = lon2-360.0;
	  cLon = x * lon1 + (1-x) * lon2 + lonshift;
	  if (cLon < 0  ) cLon += 360.0;
	  if (cLon > 360) cLon -= 360.0;

	  for (iGps=0; iGps<nGps_prop; iGps++) {
	    NormPowerAllGps[iGps] = -1.0e32;
	    iOrder[iGps] = nGps_prop;
	  }
	  MaxNormPower = -1.0e32;
	  iMaxNormPower = -1;
	  iNumber = 0;

	  for (iGps=0; iGps<nGps_prop; iGps++) {
	    //for (iGps=0; iGps<1; iGps++) {

	    if (GPS[iGps].IsGood == 1) {

	      // CBV 05/09/16
	      /* gLat = x * GPS[iGps].lat[iPt] + (1-x) * GPS[iGps].lat[iPt-1]; */
	      /* gR   = x * GPS[iGps].radius[iPt] + (1-x) * GPS[iGps].radius[iPt-1]; */

	      /* lon1 = GPS[iGps].lon[iPt]; */
	      /* lon2 = GPS[iGps].lon[iPt-1]; */
	      gLat = x * GPS[iGps].lat[iPt+1] + (1-x) * GPS[iGps].lat[iPt];
	      gR   = x * GPS[iGps].radius[iPt+1] + (1-x) * GPS[iGps].radius[iPt];
	      
	      lon1 = GPS[iGps].lon[iPt+1];
	      lon2 = GPS[iGps].lon[iPt];


	      // END CBV 05/09/16
	      if (lon1 > 300.0 && lon2 < 60.0) lon1 = lon1-360.0;
	      if (lon2 > 300.0 && lon1 < 60.0) lon2 = lon2-360.0;
	      gLon = x * lon1 + (1-x) * lon2;
	      if (gLon < 0) gLon += 360.0;

	      err = specular_point(cLat, cLon, cR,
				   gLat, gLon, gR,
				   0,
				   &specularpoint);


	      if (err == 0) {

		cVx = x * CYGNSS[iCygnss].vx[iPt+1] + (1-x) * CYGNSS[iCygnss].vx[iPt];
		cVy = x * CYGNSS[iCygnss].vy[iPt+1] + (1-x) * CYGNSS[iCygnss].vy[iPt];
		cVz = x * CYGNSS[iCygnss].vz[iPt+1] + (1-x) * CYGNSS[iCygnss].vz[iPt];

		// CBV 05/09/16
		/* gVx = x * GPS[iGps].vx[iPt] + (1-x) * GPS[iGps].vx[iPt-1]; */
		/* gVy = x * GPS[iGps].vy[iPt] + (1-x) * GPS[iGps].vy[iPt-1]; */
		/* gVz = x * GPS[iGps].vz[iPt] + (1-x) * GPS[iGps].vz[iPt-1]; */

		gVx = x * GPS[iGps].vx[iPt+1] + (1-x) * GPS[iGps].vx[iPt];
		gVy = x * GPS[iGps].vy[iPt+1] + (1-x) * GPS[iGps].vy[iPt];
		gVz = x * GPS[iGps].vz[iPt+1] + (1-x) * GPS[iGps].vz[iPt];
		// END CBV 05/09/16

		// Find the yaw of the satellite:
		if (CYGNSS[0].nYawTimes > 0) {
		  iYawPoint = iPt % CYGNSS[0].nYawTimes;
		  iYawPointBefore = iYawPoint-1;
		  if (iYawPointBefore == -1) iYawPointBefore = CYGNSS[0].nYawTimes-1;

		  yawoverall = x * CYGNSS[0].yaw[iYawPoint] + (1-x) * CYGNSS[0].yaw[iYawPointBefore];
		  //printf("Changed Yaw : %f %d %d\n", yawoverall, iPt, iYawPoint);
		}

		if (iDebug > 0)
		  printf("find gain! iCygnss, iGps, iPts : %d %d %d\n", iCygnss, iGps, iPt);


int printVar = 0;
		  if (iPtInner == 0){
		  if (iPt == 2 ){
		    if (iCygnss == 3){
printVar = 1;
		    }
}
}


		  err = find_gain(specularpoint.CygX,specularpoint.CygY,specularpoint.CygZ,
				cVx, cVy, cVz,
				specularpoint.GpsX,specularpoint.GpsY,specularpoint.GpsZ,
				gVx, gVy, gVz,
				specularpoint.SpX,specularpoint.SpY,specularpoint.SpZ,
				  45.0, rotate, yawoverall, GPS[iGps].name, printVar, iGps, iPt, iPtInner,nPtsInner, ATTITUDE,
				&specularpoint);

		
		// Save all of the specular point and GPS information:

		// cbv 04-27-18
		SpecularPointsAllGps[iGps].elev = specularpoint.elev;
		SpecularPointsAllGps[iGps].azim = specularpoint.azim;
		// end of cbv 04-27-18

		// CBV 05/09/16
		SpecularPointsAllGps[iGps].SpX    = specularpoint.SpX;
		SpecularPointsAllGps[iGps].SpY    = specularpoint.SpY;
		SpecularPointsAllGps[iGps].SpZ    = specularpoint.SpZ;

    /* // !!!!!!!!!!!!!!!! COMMENT THESE LINES BELOW! These lines are just when I want to look at the positions of the specular points with respect to the CYGNSS and GPS satellites */
    /* 	      SpecularPointsAllGps[iGps].GpsX = specularpoint.GpsX; */
    /* 	      SpecularPointsAllGps[iGps].GpsY = specularpoint.GpsY; */
    /* 	      SpecularPointsAllGps[iGps].GpsZ = specularpoint.GpsZ; */
    /* 	      SpecularPointsAllGps[iGps].CygX = specularpoint.CygX; */
    /* 	      SpecularPointsAllGps[iGps].CygY = specularpoint.CygY; */
    /* 	      SpecularPointsAllGps[iGps].CygZ = specularpoint.CygZ; */
    /* // !!!!!!!!!!!!!!!! END OF COMMENT THESE LINES BELOW! */
// END CBV 05/09/16
		SpecularPointsAllGps[iGps].longitude    = specularpoint.longitude;
		SpecularPointsAllGps[iGps].latitude     = specularpoint.latitude;


		if (iDebug > 0)
		  printf("lon, lat : c: %f %f %f g: %f %f %f s: %f %f np: %e %f los: %f %f\n",
			 cLon, cLat, cR, gLon, gLat, gR,
			 specularpoint.longitude,specularpoint.latitude,
			 specularpoint.power, specularpoint.gain,
			 specularpoint.CygnssLOSVel,specularpoint.GpsLOSVel );

		SpecularPointsAllGps[iGps].radius       = specularpoint.radius;
		SpecularPointsAllGps[iGps].CygnssLOSVel = specularpoint.CygnssLOSVel;
		SpecularPointsAllGps[iGps].GpsLOSVel    = specularpoint.GpsLOSVel;
		SpecularPointsAllGps[iGps].gain         = specularpoint.gain;
		SpecularPointsAllGps[iGps].power        = specularpoint.power;
		SpecularPointsAllGps[iGps].NormPower    = specularpoint.NormPower;
		SpecularPointsAllGps[iGps].Incidence    = specularpoint.Incidence;

		gLonAllGps[iGps] = (float) gLon;
		gLatAllGps[iGps] = (float) gLon;
		gRAllGps[iGps]   = (float) gR;
		
		// Here we are trying to find the highest normalized gain (i.e., don't know the
		// Sigma0, which is wind speed dependent):
		NormPowerAllGps[iGps] = specularpoint.NormPower;
		if (NormPowerAllGps[iGps] > MaxNormPower) {
		  MaxNormPower = NormPowerAllGps[iGps];
		  iMaxNormPower = iGps;
		}

		//		printf("in err %d\n", iGps);
	      } // end of if err == 0 (err is result of specular_point functino)
	      //		printf("in isgood %d\n", iGps);
	    } // end of if (GPS[iGps].IsGood == 1)


	    //		printf("normpower[%d]: %e\n",iGps, NormPowerAllGps[iGps]);


	    
	  } // loop over gps



	  // Store the highest normalized gain GSP satellite:
	  if (iMaxNormPower > -1) {
	    iOrder[iMaxNormPower] = iNumber;
	    NormPowerAllGps[iMaxNormPower] = -1.1e32;
	    iNumber++;
	  }

	  // Determine the order of the other GPS satellites:
	  while (MaxNormPower > -1.0e32) {
	    MaxNormPower = -1.0e32;
	    iMaxNormPower = -1;
	    for (iGps=0; iGps<nGps_prop; iGps++) {
	      if (NormPowerAllGps[iGps] > MaxNormPower) {
		MaxNormPower = NormPowerAllGps[iGps];
		iMaxNormPower = iGps;
	      }
	    }
	    if (iMaxNormPower > -1) {
	      iOrder[iMaxNormPower] = iNumber;
	      NormPowerAllGps[iMaxNormPower] = -1.1e32;
	      iNumber++;
	    }
	    //	    printf("iNumber %d\n",iNumber);
	      }

/* 	  int sss; */
/* 	  if (iPtInner == 0){ */
/* 	  printf("\n"); */
/* 	  for ( sss = 0; sss < 7; sss++){ */
/* 	  printf("%d ",  GPS[iGps].time_ymdhmsm[iPt][sss]); */
/* 	  } */

/* 	  printf("\n"); */
/* 	  } */
	  for (iGps=0; iGps<nGps_prop; iGps++) {

	    if (iOrder[iGps] < nSpecularPointsMax) {
	      if (iCygnss == 3){
		if (iPtInner == 0){
		  if (iPt == 2 ){//> nPts - 4 ){
		    //		       if ( strcmp(GPS[iGps].name, "PRN_13") == 0){
		    printf("GPS[%d]: <%s> | normPower %.2f | x %f | elev %f, azim %f | %d\n", iGps, GPS[iGps].name, SpecularPointsAllGps[iGps].NormPower, GPS[iGps].x[iPt]/1000., SpecularPointsAllGps[iGps].elev, SpecularPointsAllGps[iGps].azim, iOrder[iGps]);
		    //}

		  }
		}
	      }
	      //	      fwrite(&iCygnss,  sizeof(iCygnss),  1,fpOutput); // CBV 05/09/16
	      fwrite(&iGps,     sizeof(iGps),     1,fpOutput);

	      fwrite(&iPtTotal, sizeof(iPtTotal), 1,fpOutput);

	      if (iUseMin == 0) {

		tmp = (float) cLon; fwrite(&tmp,sizeof(tmp),1,fpOutput);
		tmp = (float) cLat; fwrite(&tmp,sizeof(tmp),1,fpOutput);
		tmp = (float) cR;   fwrite(&tmp,sizeof(tmp),1,fpOutput);

		//tmp = (float) cVx; fwrite(&tmp,sizeof(tmp),1,fpOutput);
		//tmp = (float) cVy; fwrite(&tmp,sizeof(tmp),1,fpOutput);
		//tmp = (float) cVz; fwrite(&tmp,sizeof(tmp),1,fpOutput);

		tmp = gLonAllGps[iGps]; fwrite(&tmp,sizeof(tmp),1,fpOutput);
		tmp = gLatAllGps[iGps]; fwrite(&tmp,sizeof(tmp),1,fpOutput);
		tmp = gRAllGps[iGps];   fwrite(&tmp,sizeof(tmp),1,fpOutput);

	      }

	      //tmp = (float) gVx; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      //tmp = (float) gVy; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      //tmp = (float) gVz; fwrite(&tmp,sizeof(tmp),1,fpOutput);

		// CBV 05/09/16
	      int sss;
	      if (iPtInner == 0){
	      for ( sss = 0; sss < 7; sss++){
	      	tmp = (float) GPS[iGps].time_ymdhmsm[iPt][sss]; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      }
	      }


	      tmp =  SpecularPointsAllGps[iGps].SpX/1000.; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  SpecularPointsAllGps[iGps].SpY/1000.; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  SpecularPointsAllGps[iGps].SpZ/1000.; fwrite(&tmp,sizeof(tmp),1,fpOutput);

    /* // !!!!!!!!!!!!!!!! COMMENT THESE LINES BELOW! These lines are just when I want to look at the positions of the specular points with respect to the CYGNSS and GPS satellites */
    /* 	      tmp =  SpecularPointsAllGps[iGps].GpsX/1000.; fwrite(&tmp,sizeof(tmp),1,fpOutput); */
    /* 	      tmp =  SpecularPointsAllGps[iGps].GpsY/1000.; fwrite(&tmp,sizeof(tmp),1,fpOutput); */
    /* 	      tmp =  SpecularPointsAllGps[iGps].GpsZ/1000.; fwrite(&tmp,sizeof(tmp),1,fpOutput); */
    /* 	      tmp =  SpecularPointsAllGps[iGps].CygX/1000.; fwrite(&tmp,sizeof(tmp),1,fpOutput); */
    /* 	      tmp =  SpecularPointsAllGps[iGps].CygY/1000.; fwrite(&tmp,sizeof(tmp),1,fpOutput); */
    /* 	      tmp =  SpecularPointsAllGps[iGps].CygZ/1000.; fwrite(&tmp,sizeof(tmp),1,fpOutput); */
    /* // !!!!!!!!!!!!!!!! END OF COMMENT THESE LINES BELOW! */

		// END CBV 05/09/16
	      tmp = (float) SpecularPointsAllGps[iGps].longitude; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp = (float) SpecularPointsAllGps[iGps].latitude;  fwrite(&tmp,sizeof(tmp),1,fpOutput);

	      if (iUseMin == 0) {
		tmp = (float) SpecularPointsAllGps[iGps].radius;  fwrite(&tmp,sizeof(tmp),1,fpOutput);
		tmp = (float) SpecularPointsAllGps[iGps].CygnssLOSVel;
		fwrite(&tmp,sizeof(tmp),1,fpOutput);
		tmp = (float) SpecularPointsAllGps[iGps].GpsLOSVel;
		fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      }
	      tmp = (float) SpecularPointsAllGps[iGps].gain;
	      fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      // CBV 05/09/16
	      //	      if ((iPtInner == 0) && (iPt %2 == 0)){// THAT IS TO OUTPUT THE POSITION AND VELOCITY OF THE GPS/CYGNSS EVERY 60 SECONDS (WHEN THE PROPAGATOR OUTPUTS EVERY 30 SECONDS (THIS EXPLAINS THE iPt %2 == 0))

	      gps_x_inner = x * GPS[iGps].x[iPt+1] + (1-x) * GPS[iGps].x[iPt];
	      gps_y_inner = x * GPS[iGps].y[iPt+1] + (1-x) * GPS[iGps].y[iPt];
	      gps_z_inner = x * GPS[iGps].z[iPt+1] + (1-x) * GPS[iGps].z[iPt];
	      tmp =  (float) gps_x_inner/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  (float) gps_y_inner/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  (float) gps_z_inner/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput);

	      gps_vx_inner = x * GPS[iGps].vx[iPt+1] + (1-x) * GPS[iGps].vx[iPt];
	      gps_vy_inner = x * GPS[iGps].vy[iPt+1] + (1-x) * GPS[iGps].vy[iPt];
	      gps_vz_inner = x * GPS[iGps].vz[iPt+1] + (1-x) * GPS[iGps].vz[iPt];
	      tmp =  (float) gps_vx_inner; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  (float) gps_vy_inner; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  (float) gps_vz_inner; fwrite(&tmp,sizeof(tmp),1,fpOutput);

	      cygnss_x_inner = x * CYGNSS[iCygnss].x[iPt+1] + (1-x) * CYGNSS[iCygnss].x[iPt];
	      cygnss_y_inner = x * CYGNSS[iCygnss].y[iPt+1] + (1-x) * CYGNSS[iCygnss].y[iPt];
	      cygnss_z_inner = x * CYGNSS[iCygnss].z[iPt+1] + (1-x) * CYGNSS[iCygnss].z[iPt];
	      tmp =  (float) cygnss_x_inner/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  (float) cygnss_y_inner/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  (float) cygnss_z_inner/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput);

	      cygnss_vx_inner = x * CYGNSS[iCygnss].vx[iPt+1] + (1-x) * CYGNSS[iCygnss].vx[iPt];
	      cygnss_vy_inner = x * CYGNSS[iCygnss].vy[iPt+1] + (1-x) * CYGNSS[iCygnss].vy[iPt];
	      cygnss_vz_inner = x * CYGNSS[iCygnss].vz[iPt+1] + (1-x) * CYGNSS[iCygnss].vz[iPt];
	      tmp =  (float) cygnss_vx_inner; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  (float) cygnss_vy_inner; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  (float) cygnss_vz_inner; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      // used to be below but the problem is that it doesnt show the r/v interpolated
/* 	      tmp =  (float) GPS[iGps].x[iPt]/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput); */
/* 	      tmp =  (float) GPS[iGps].y[iPt]/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput); */
/* 	      tmp =  (float) GPS[iGps].z[iPt]/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput); */

/* 	      tmp =  (float) GPS[iGps].vx[iPt]; fwrite(&tmp,sizeof(tmp),1,fpOutput); */
/* 	      tmp =  (float) GPS[iGps].vy[iPt]; fwrite(&tmp,sizeof(tmp),1,fpOutput); */
/* 	      tmp =  (float) GPS[iGps].vz[iPt]; fwrite(&tmp,sizeof(tmp),1,fpOutput); */

/* 	      tmp =  (float) CYGNSS[iCygnss].x[iPt]/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput); */
/* 	      tmp =  (float) CYGNSS[iCygnss].y[iPt]/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput); */
/* 	      tmp =  (float) CYGNSS[iCygnss].z[iPt]/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput); */

/* 	      tmp =  (float) CYGNSS[iCygnss].vx[iPt]; fwrite(&tmp,sizeof(tmp),1,fpOutput); */
/* 	      tmp =  (float) CYGNSS[iCygnss].vy[iPt]; fwrite(&tmp,sizeof(tmp),1,fpOutput); */
/* 	      tmp =  (float) CYGNSS[iCygnss].vz[iPt]; fwrite(&tmp,sizeof(tmp),1,fpOutput); */
	      // end of used to be beliow
	      //	      }
	      /* /\* if ((iPtInner == 0) && (iCygnss == 0) && (iPt %2 == 0)){ *\/ */
	      /* /\* 	printf("%f %f %f\n", GPS[iGps].x[iPt]/1000, GPS[iGps].y[iPt]/1000, GPS[iGps].z[iPt]/1000); *\/ */
	      /* } */


	      /* tmp = (float) SpecularPointsAllGps[iGps].power; */
	      /* fwrite(&tmp,sizeof(tmp),1,fpOutput); */
	      tmp = (float) SpecularPointsAllGps[iGps].NormPower;
	      fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      // END CBV 05/09/16
	      if (iUseMin == 0) {
		tmp = (float) SpecularPointsAllGps[iGps].Incidence;
		fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      }

	      iIncidence = (int) (SpecularPointsAllGps[iGps].Incidence+0.5);
	      if (iIncidence < 0 ) iIncidence = 0;
	      if (iIncidence > 89) iIncidence = 89;

	      for (iWind = 0; iWind < nWinds; iWind++) {
		power = SpecularPointsAllGps[iGps].NormPower * AntennaInfo.Sigma0[iWind][iIncidence];
		ratio = power/AntennaInfo.noisefloor;
		if (ratio > 127) ratio = 127;
		if (ratio < 0  ) ratio = 0;
		iRatio[iWind] = (char) ratio;
	      }
	      // CBV 05/09/16
	      // fwrite(iRatio,sizeof(char),nWinds,fpOutput);
	      // END CBV 05/09/16
	    }

	  }
	  //	  MPI_Finalize();exit(0);


	}


      }

    }

  }

  close(fpOutput);

  return err;

}

//-----------------------------------------------------------------------
//
//-----------------------------------------------------------------------

double calc_mean_radius(struct SatelliteDetails *sat, int nPts) {

  int    i, j, k, nPtsTotal;
  double ave;

  ave = 0.0;
  nPtsTotal = 0;

  for (i=0; i<nPts; i++) {
    if (sat[i].IsGood == 1) {
      for (j=0; j < sat[i].nPts; j++) {
	ave = ave + sat[i].radius[j];
	nPtsTotal++;
      }
    }
  }

  ave = ave/nPtsTotal;

}





// CBV 05/08
void RemoveSpaces(char* source)
{
  char* i = source;
  char* j = source;
  while(*j != 0)
    {
      *i = *j++;
      if(*i != ' ')
	i++;
    }
  *i = 0;
}
// END CBV 05/08




int v_print( double v_to_print[3],
	     char name[256])
{
  printf("%s: (%f, %f, %f)\n",name, v_to_print[0], v_to_print[1], v_to_print[2]);
  return 0;
}



int v_dot(  double *dot,
            double v1[3],
            double v2[3])

{
    
    dot[0] = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
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

// different from the one in prop_math.c (in prop_math it's inrtl not ecef, and the y is port not startboard and the z in zenith not nadir)

int compute_T_ecef_2_ntw_not_like_prop_math( double T_ecef_2_ntw[3][3],
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
  v_cross(w,v,r);
  v_norm(w,w);

  /* In-track cross cross-track */
  v_cross(n,t,w);

  for (col = 0; col < 3; col++){
    T_ecef_2_ntw[0][col] = t[col]; // in-track
    T_ecef_2_ntw[1][col] = w[col]; // cross-track
    T_ecef_2_ntw[2][col] = n[col]; // in-track cross cross-track
  }

    return 0;
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


// Compute the vector mag
int v_mag( double *v_mag,
           double v_in[3])

{

    v_mag[0] = sqrt( v_in[0]*v_in[0] + v_in[1]*v_in[1] + v_in[2]*v_in[2]);
    return 0;
}



// an adaptation of the function compute_T_sc_to_lvlh  in prop_math.c. by adaptation: no sun pionted, no quaternion.also for the spec calculations, pitch toll andy= yaw are with respoect to the NTW refce system, not LVLH. The diffce between NTW and LVLH is taht x is velocity in NTW. In lvlh z is position. So NTW is different from LVLH when the eccentricity is different from 0
int compute_T_sc_to_ntw(double T_sc_to_ntw[3][3], double v_angle[3], int order_rotation[3]){

  /* Declarations */
  int row;
  double x[6]; 
  double lt;
  double r_earth2sun_J2000[3]; 
  double r_cg2sun_J2000[3];
  double r_cg2sun_J2000_normalized[3]; 
  double e_z_body_in_inrtl[3];
  double T_inrtl_2_ntw[3][3];
  double e_z_body_in_ntw[3];
  /* double random_vect_not_colinear_to_z[3];  */
  /* double e_y_body_in_ntw[3];  */
  /* double e_x_body_in_ntw[3];  */
  double e_z_body_in_ntw_normalized[3]; 
  /* double e_y_body_in_ntw_normalized[3]; */
  /* double e_x_body_in_ntw_normalized[3]; */
  int i;
  double theta, phi, psi;
  double pitch_mat[3][3], roll_mat[3][3], yaw_mat[3][3];
  double first_mat[3][3], second_mat[3][3], third_mat[3][3];
  double T_sc_to_ntw_temp[3][3];
  //  double T_sc_intermediary_to_ntw[3][3];
  double v_angle_rad[3];


/*   double *p=malloc(sizeof(double)), *r=malloc(sizeof(double)), *y=malloc(sizeof(double)); */
/*   if (file_is_quaternion == 1){ // if attitude is set using quaternions      */
/*     q2m_c(quaternion,T_sc_to_ntw); */
/* /\*             m_print(T_sc_to_ntw, "quat T_sc_to_ntw"); *\/ */
/* /\* 	    m2eul_c ( T_sc_to_ntw, 2, 1, 3, p, r, y ); *\/ */
/* /\* 	    	    printf("%f %f %f\n",p[0]*180./M_PI, r[0]*180./M_PI, y[0]*180./M_PI ); *\/ */
	      
/*   } */

/*   ////////////////////////////////////////////////////////////////////////////// */
/*   ////////////////////////////////////////////////////////////////////////////// */
/*   ///////////////////// NOT QUATERNION ///////////////////////////////////////////// */
/*   else{ //no quaternion */


  //////////////////////////////////////////////////////////////////////////////
  ///////////////////// ANY ATTITUDE EXCEPT SUN POINTED ////////////////////////
  //  if (  strcmp( attitude_profile, "sun_pointed"  ) != 0 )  { // if the attitude profile is anything except sun pointed
  for (i = 0; i < 3; i++){
    v_angle_rad[i] = v_angle[i] *  dtor;
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
    printf("***! (compute_T_sc_to_ntw) You did not correctly choose the first rotation. The program will stop. !***\n"); MPI_Finalize(); 
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
    printf("***! (compute_T_sc_to_ntw) You did not correctly choose the second rotation. The program will stop. !***\n"); MPI_Finalize(); 
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
    printf("***! (compute_T_sc_to_ntw) You did not correctly choose the third rotation. The program will stop. !***\n"); MPI_Finalize(); 
    exit(0);
  }
  // // Final matrix: T_sc_to_ntw
  m_x_m(T_sc_to_ntw_temp, second_mat, first_mat);
  m_x_m(T_sc_to_ntw, third_mat, T_sc_to_ntw_temp);
  //            m_print(T_sc_to_ntw, "euler T_sc_to_ntw");
  //  }

  //////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// SUN POINTED ///////////////////////////////////
/*   if (  strcmp( attitude_profile, "sun_pointed"  ) == 0 )  { */
/*     /\* Sun-pointed *\/ */
/*     /\* Express e_z_body in Inertial frame (J2000) *\/ */
/*     spkez_c(10, et[0], "J2000", "NONE", 399, x, &lt); //   Return the state (position and velocity) of a target body relative to an observing body, optionally corrected for light time (planetary aberration) and stellar aberration. */
/*     r_earth2sun_J2000[0] = x[0]; */
/*     r_earth2sun_J2000[1] = x[1]; */
/*     r_earth2sun_J2000[2] = x[2]; */


/*     v_sub(r_cg2sun_J2000, r_earth2sun_J2000, r_i2cg_INRTL); */
/*     v_norm(r_cg2sun_J2000_normalized, r_cg2sun_J2000); */
/*     v_copy(e_z_body_in_inrtl, r_cg2sun_J2000_normalized); */

/*     /\* Convert e_z_body in NTW frame *\/ */
/*     compute_T_inrtl_2_ntw(T_inrtl_2_ntw, r_i2cg_INRTL, v_i2cg_INRTL); */
/*     m_x_v(e_z_body_in_ntw, T_inrtl_2_ntw, e_z_body_in_inrtl); */
/*     v_norm(e_z_body_in_ntw_normalized,e_z_body_in_ntw); // useless cause e_z_body_in_ntw should already be normalized */


/*     /\************************************** BLOCK A: UNCOMMENT BELOW TO SET X_BODY PERPENDICULAR TO NADIR **********************************\/ */
/*     // // IMPORTANT: IF YOU UNCOMMENT THIS BLOCK THEN COMMENT BLOCK B */
/*     // Sun pointed gives the direction of z_body. The sc can still yaw around this vector sc to Sun. By default, we assume that the yaw angle is so that the x_body is perpendicular to nadir. In other words, x_body is the cross product of nadir with sat_to_sun vector (both expressed in ntw coordinates) */
/*     double ntw_z_in_ntw[3], new_e_x_body_in_ntw[3], new_e_y_body_in_ntw_normalized[3], new_e_x_body_in_ntw_normalized[3]; */
/*     ntw_z_in_ntw[0] = 0; ntw_z_in_ntw[1] = 0; ntw_z_in_ntw[2] = 1; // ntw_z is in the nadir (oriented away from the Earth) */
/*     v_cross(new_e_x_body_in_ntw, ntw_z_in_ntw, e_z_body_in_ntw_normalized); */
/*     v_norm(new_e_x_body_in_ntw_normalized, new_e_x_body_in_ntw); */
    
/*     v_cross(new_e_y_body_in_ntw_normalized, e_z_body_in_ntw_normalized, new_e_x_body_in_ntw_normalized); */


/*     // // !!!!!!! BELOW IS TO ROTATE X_BODY BY 10 DEGREES IN THE PLANE PERPENDICULAR TO SAT_TO_SUN VECTOR. THEN X_BODY IS 10 DEGREES FROM BEING PERPENDICULAR TO NADIR (BUT STILL PERPENDICULAR TO Z_BODY (WHICH IS EQUAL TO SAT_TO_SUN VECTOR)). PARDON THE LONG NAMES */
/*     // // !!!!!! THEREFORE, COMMENT THE BLOCK BELOW (UNLESS YOU KNOW WHAT YOU ARE DOING) */
/*     double new_e_x_body_in_ntw_normalized_rotated_10_degrees[3], new_e_x_body_in_ntw_normalized_rotated_10_degrees_component_along_new_e_x_body_in_ntw_normalized[3], new_e_x_body_in_ntw_normalized_rotated_10_degrees_component_along_new_e_y_body_in_ntw_normalized[3],new_e_x_body_in_ntw_normalized_rotated_10_degrees_component_along_e_z_body_in_ntw_normalized[3], new_e_x_body_in_ntw_normalized_rotated_10_degrees_temp[3], new_e_x_body_in_ntw_normalized_rotated_10_degrees_normalized[3]; */
/*     v_scale(new_e_x_body_in_ntw_normalized_rotated_10_degrees_component_along_new_e_x_body_in_ntw_normalized, new_e_x_body_in_ntw_normalized, cos( 10. * dtor )); */
/*     v_scale( new_e_x_body_in_ntw_normalized_rotated_10_degrees_component_along_new_e_y_body_in_ntw_normalized, new_e_y_body_in_ntw_normalized, sin( 10. * dtor )); */
/*     new_e_x_body_in_ntw_normalized_rotated_10_degrees_component_along_e_z_body_in_ntw_normalized[0] = 0;  new_e_x_body_in_ntw_normalized_rotated_10_degrees_component_along_e_z_body_in_ntw_normalized[1] = 0;  new_e_x_body_in_ntw_normalized_rotated_10_degrees_component_along_e_z_body_in_ntw_normalized[2] = 0; */

/*     v_add( new_e_x_body_in_ntw_normalized_rotated_10_degrees_temp, new_e_x_body_in_ntw_normalized_rotated_10_degrees_component_along_new_e_x_body_in_ntw_normalized, new_e_x_body_in_ntw_normalized_rotated_10_degrees_component_along_new_e_y_body_in_ntw_normalized ); */
/*     v_add( new_e_x_body_in_ntw_normalized_rotated_10_degrees, new_e_x_body_in_ntw_normalized_rotated_10_degrees_temp, new_e_x_body_in_ntw_normalized_rotated_10_degrees_component_along_e_z_body_in_ntw_normalized ); */
/*     v_norm( new_e_x_body_in_ntw_normalized_rotated_10_degrees_normalized, new_e_x_body_in_ntw_normalized_rotated_10_degrees); */
/*     new_e_x_body_in_ntw_normalized[0] = 0;  new_e_x_body_in_ntw_normalized[1] = 0;  new_e_x_body_in_ntw_normalized[2] = 0; */
/*     v_copy( new_e_x_body_in_ntw_normalized, new_e_x_body_in_ntw_normalized_rotated_10_degrees_normalized ); */
/*     new_e_y_body_in_ntw_normalized[0] = 0;  new_e_y_body_in_ntw_normalized[1] = 0;  new_e_y_body_in_ntw_normalized[2] = 0; */
/*     v_cross( new_e_y_body_in_ntw_normalized, e_z_body_in_ntw_normalized, new_e_x_body_in_ntw_normalized ); */

/*     // // !!!!!! END OF 'THEREFORE, COMMENT THE BLOCK BELOW (UNLESS YOU KNOW WHAT YOU ARE DOING)' */

/*     /\************************************** END OF 'BLOCK A: UNCOMMENT BELOW TO SET X_BODY PERPENDICULAR TO NADIR' **********************************\/ */

/*     /\************************************** BLOCK B: UNCOMMENT BELOW TO SET X_BODY IN THE DIRECTION OF NTW_X **********************************\/ */
/*     /\* /\\* // // IMPORTANT: IF YOU UNCOMMENT THIS BLOCK THEN COMMENT BLOCK A *\\/ *\/ */
/*     /\* /\\* // Sun pointed gives the direction of z_body. The sc can still roll around this vector sc to Sun. By default, we assume that the roll angle is so that the x_body points towards ntw_x. In other words, the projection of lvlv_x on the plane perpendicular to the Sun give x_body *\\/ *\/ */
/*     /\* double v_int[3], ntw_x_in_ntw[3], new_e_x_body_in_ntw[3], new_e_y_body_in_ntw_normalized[3], v_int_norm[3], new_e_x_body_in_ntw_normalized[3]; *\/ */
    
/*     /\* ntw_x_in_ntw[0] = 1; ntw_x_in_ntw[1] = 0; ntw_x_in_ntw[2] = 0; *\/ */
/*     /\* v_cross(v_int, e_z_body_in_ntw_normalized, ntw_x_in_ntw); *\/ */
/*     /\* v_norm(v_int_norm, v_int); *\/ */
/*     /\* v_cross(new_e_x_body_in_ntw, v_int_norm, e_z_body_in_ntw_normalized); *\/ */
/*     /\* v_norm(new_e_x_body_in_ntw_normalized, new_e_x_body_in_ntw); *\/ */
/*     /\* v_cross(new_e_y_body_in_ntw_normalized, e_z_body_in_ntw_normalized, new_e_x_body_in_ntw_normalized); *\/ */

/*     /\* /\\*   /\\\* /\\\\* // Another way: *\\\\/ *\\\/ *\\/ *\/ */
/*     /\* /\\*   /\\\* double u_dot_n, u_dot_n_scale[3]; *\\\/ *\\/ *\/ */
/*     /\* /\\*   /\\\* v_dot(&u_dot_n, ntw_x_in_ntw, e_z_body_in_ntw_normalized); *\\\/ *\\/ *\/ */
/*     /\* /\\*   /\\\* v_scale(u_dot_n_scale, e_z_body_in_ntw_normalized,u_dot_n); *\\\/ *\\/ *\/ */
/*     /\* /\\*   /\\\* v_sub(new_e_x_body_in_ntw_normalized, ntw_x_in_ntw, u_dot_n_scale); *\\\/ *\\/ *\/ */
/*     /\* /\\*   /\\\* v_norm(new_e_x_body_in_ntw_normalized, new_e_x_body_in_ntw_normalized); *\\\/ *\\/ *\/ */
/*     /\* /\\*   /\\\* v_print(new_e_x_body_in_ntw_normalized, "new_e_x_body_in_ntw_normalized"); *\\\/ *\\/ *\/ */
/*     /\************************************** END OF 'BLOCK B: UNCOMMENT BELOW TO SET X_BODY IN THE DIRECTION OF NTW_X' **********************************\/ */

/*     /\* // We now have the new body frame, expressed in the NTW frame: the z axis points towards the Sun, and the x axis is the projection of NTW_X on the plane perdicular to the satellite-to-Sun vector. *\/ */




/*     for (row = 0; row<3; row++){ */
/*       T_sc_to_ntw[row][0] = new_e_x_body_in_ntw_normalized[row]; */
/*       T_sc_to_ntw[row][1] = new_e_y_body_in_ntw_normalized[row]; */
/*       T_sc_to_ntw[row][2] = e_z_body_in_ntw_normalized[row]; */
/*     } */
    

/*   } */
  //  } // not quaternion


  return 0;

  // // OLD STUFF FOR SUN POINTED:
    /* /\* Now the z axis of the body is fixed, set the rotation around this z axis to complete the (x, y, z)body *\/ */
    /* /\* For now, we just take any "random" vector orthogonal to e_z_body to be e_y_body. This means that the x and y axes have a "random" orientation in the plane orthonal to the satellite-Sun direction. *\/ */
    /* /\* Creation of a vector not colinear to e_z_body *\/ */
    /* if (e_z_body_in_ntw[0] != 0){ */
    /*   random_vect_not_colinear_to_z[1] = 1.0; */
    /*   random_vect_not_colinear_to_z[0] = 0.0; */
    /*   random_vect_not_colinear_to_z[2] = 0.0; */
    /* } */
    /* else if (e_z_body_in_ntw[1] != 0){ */
    /*   random_vect_not_colinear_to_z[0] = 1.0; */
    /*   random_vect_not_colinear_to_z[1] = 0.0; */
    /*   random_vect_not_colinear_to_z[2] = 0.0; */
    /* } */
    /* else if (e_z_body_in_ntw[2] != 0){ */
    /*   random_vect_not_colinear_to_z[1] = 1.0; */
    /*   random_vect_not_colinear_to_z[2] = 0.0; */
    /*   random_vect_not_colinear_to_z[0] = 0.0; */
    /* } */
    /* else{ */
    /*   printf("It seems that the vector satellite-to-Sun is 0. The program will stop.\n"); */
    /*   exit(0); */
    /* } */

    /* /\* e_y_body is the cross product between e_z_body and this "random" non-colinear vector to e_z_body. Consequently, e_y_body is a "random" vector orthogonal to e_z_body *\/ */
    /* v_cross(e_y_body_in_ntw, e_z_body_in_ntw, random_vect_not_colinear_to_z); */
    /* v_cross(e_x_body_in_ntw, e_y_body_in_ntw, e_z_body_in_ntw); */

    /* /\* Normalization of the SC basis *\/ */
    /* v_norm(e_x_body_in_ntw_normalized,e_x_body_in_ntw); */
    /* v_norm(e_y_body_in_ntw_normalized,e_y_body_in_ntw); */
    /* v_norm(e_z_body_in_ntw_normalized,e_z_body_in_ntw); // useless cause e_z_body_in_ntw should already be normalized */

    /* /\* T_sc_intermediary_to_ntw is the rotation from body to ntw *\/ */
    /* for (row = 0; row<3; row++){ */
    /*   T_sc_intermediary_to_ntw[row][0] = e_x_body_in_ntw_normalized[row]; */
    /*   T_sc_intermediary_to_ntw[row][1] = e_y_body_in_ntw_normalized[row]; */
    /*   T_sc_intermediary_to_ntw[row][2] = e_z_body_in_ntw            _normalized[row]; */
    /* } */

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


int etprint( double et_to_print, char str_print[256] ){
  char str[256];
  et2utc_c(et_to_print, "ISOC", 6, 255, str);
  printf("%s: %s\n", str_print,str);
  
  return 0;
}


int nb_time_steps_f(int *nb_time_steps_simu,
		  double et_initial_epoch,
		  char final_epoch[256],
		  double dt){

  /* Declarations */

  double et_final;
  /* Algorithm */

  str2et_c(final_epoch, &et_final);
  //  printf("%f %f %f %f %f %f %d\n",et_final, et_initial_epoch, ( et_final - et_initial_epoch ), dt, ( et_final - et_initial_epoch ) / dt, ceil( ( et_final - et_initial_epoch ) / dt ), (int)(ceil( ( et_final - et_initial_epoch ) / dt ) ));
  *nb_time_steps_simu = (int)(ceil( ( et_final - et_initial_epoch ) / dt ) ) + 1; 

  return 0;

}






/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           load_attitude
//  Purpose:        Load the attitude file (nadir, sun_pointed or manual (= name of external file))
//  Assumptions:    Same attitude configuration for each satellite of the constellation
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 08/12/2015    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////


int load_attitude( ATTITUDE_T *ATTITUDE,
		    FILE *input_file, int nProcs, double ang_velo[3], int iProc){

  /* Declarations */
  int j;
  int ierr;
  int iia;
  char *line = NULL;
  size_t len = 0;
  int found_eoh = 0;
  char text_location[256];
  char temp_copy[256];
  char text[256];


  //char times[256];
  /* Algorithm */
  /* nadir */


  if (strcmp(ATTITUDE->attitude_profile, "nadir") == 0){
    //    printf("here nb_time_steps %d\n", ATTITUDE->nb_time_steps);
    ATTITUDE->et_interpo = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );

    ATTITUDE->pitch = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );
    ATTITUDE->roll = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );
    ATTITUDE->yaw = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );
    ATTITUDE->order_pitch = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );
    ATTITUDE->order_roll = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );
    ATTITUDE->order_yaw = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );

    if (  ATTITUDE->pitch == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->pitch \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  ATTITUDE->roll == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->roll \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  ATTITUDE->yaw == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->yaw \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }

    if (  ATTITUDE->order_pitch == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->order_pitch \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  ATTITUDE->order_roll == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->order_roll \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  ATTITUDE->order_yaw == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->order_yaw \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }


    if ( ATTITUDE->et_interpo == NULL ){
      printf("***! Could not allow memory space for ATTITUDE->et_interpo  \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);

    }


  iia = 0;

  ATTITUDE->et_interpo[iia] = ATTITUDE->et_initial;  

      ATTITUDE->pitch[iia] = 0;
      ATTITUDE->roll[iia] = 0;
      ATTITUDE->yaw[iia] = 0;
      ATTITUDE->order_pitch[iia] = 1;
      ATTITUDE->order_roll[iia] = 2;
      ATTITUDE->order_yaw[iia] = 3;

  while (ATTITUDE->et_interpo[iia] < ATTITUDE->et_initial){
    iia = iia + 1;

    ATTITUDE->et_interpo[iia] = ATTITUDE->et_initial + 1.0*iia;

    //        etprint(ATTITUDE->et_interpo[iia], "tle");
      ATTITUDE->pitch[iia] = 0;
      ATTITUDE->roll[iia] = 0;
      ATTITUDE->yaw[iia] = 0;
      ATTITUDE->order_pitch[iia] = 1;
      ATTITUDE->order_roll[iia] = 2;
      ATTITUDE->order_yaw[iia] = 3;
	
  } // leave this loop when x_after_interpo gets newer than constellation epoch (ATTITUDE->et_initial)
  j = 0;
  while (iia<ATTITUDE->nb_time_steps ){

    ATTITUDE->et_interpo[iia] = ATTITUDE->et_initial + 1.0*j;

    //        etprint(ATTITUDE->et_interpo[iia], "const");
      ATTITUDE->pitch[iia] = 0;
      ATTITUDE->roll[iia] = 0;
      ATTITUDE->yaw[iia] = 0;
      ATTITUDE->order_pitch[iia] = 1;
      ATTITUDE->order_roll[iia] = 2;
      ATTITUDE->order_yaw[iia] = 3;

/* 	etprint(ATTITUDE->et_interpo[iia], "time"); */
/* 		printf("ATTITUDE->yaw[%d] %f\n", iia, ATTITUDE->yaw[iia]); */

	iia = iia +1;
	j = j+1;
  }


/*     for (iia = 0; iia < ATTITUDE->nb_time_steps ; iia++){ */
/*       ATTITUDE->et_interpo[iia] = ATTITUDE->et_initial + ATTITUDE->dt_prop*iia / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*       ATTITUDE->pitch[iia] = 0; */
/*       ATTITUDE->roll[iia] = 0; */
/*       ATTITUDE->yaw[iia] = 0; */
/*       ATTITUDE->order_pitch[iia] = 1; */
/*       ATTITUDE->order_roll[iia] = 2; */
/*       ATTITUDE->order_yaw[iia] = 3; */
/*     } */

  }

  /* sun_pointed */
  else if (strcmp(ATTITUDE->attitude_profile, "sun_pointed") == 0){
    printf("***! Sun pointed is not a valide attitude to compute specular point locations.\n. The program will stop. !***\n"); MPI_Finalize(); exit(0);
  
  }

  else if (strcmp(ATTITUDE->attitude_profile, "ensemble_angular_velocity") == 0){
    printf("***! Ensembles on the attitude cannot be run to compute specular point locations.\n. The program will stop. !***\n"); MPI_Finalize(); exit(0);

  }

  else if (strcmp(ATTITUDE->attitude_profile, "ensemble_initial_attitude") == 0){
    printf("***! Ensembles on the attitude cannot be run to compute specular point locations.\n. The program will stop. !***\n"); MPI_Finalize(); exit(0);

  }

  else if (strcmp(ATTITUDE->attitude_profile, "angular_velocity") == 0){

      ATTITUDE->et_interpo = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );

    ATTITUDE->pitch = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );
    ATTITUDE->roll = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );
    ATTITUDE->yaw = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );
    ATTITUDE->order_pitch = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );
    ATTITUDE->order_roll = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );
    ATTITUDE->order_yaw = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );

    if (  ATTITUDE->pitch == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->pitch \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  ATTITUDE->roll == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->roll \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  ATTITUDE->yaw == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->yaw \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }

    if (  ATTITUDE->order_pitch == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->order_pitch \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  ATTITUDE->order_roll == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->order_roll \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  ATTITUDE->order_yaw == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->order_yaw \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }


    if ( ATTITUDE->et_interpo == NULL ){
      printf("***! Could not allow memory space for ATTITUDE->et_interpo  \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }


  iia = 0;

  ATTITUDE->et_interpo[iia] = ATTITUDE->et_initial; 

      ATTITUDE->pitch[iia] = iia * ang_velo[3] * 1.0  + ang_velo[0];
      ATTITUDE->roll[iia] = iia * ang_velo[4] * 1.0  + ang_velo[1];
      ATTITUDE->yaw[iia] = iia * ang_velo[5] * 1.0  + ang_velo[2] ;
/*       ATTITUDE->pitch[iia] = iia * ang_velo[0] * ATTITUDE->dt_prop / 2.; */
/*       ATTITUDE->roll[iia] = iia * ang_velo[1] * ATTITUDE->dt_prop / 2.; */
/*       ATTITUDE->yaw[iia] = iia * ang_velo[2] * ATTITUDE->dt_prop / 2.; */
      ATTITUDE->order_pitch[iia] = 1;
      ATTITUDE->order_roll[iia] = 2;
      ATTITUDE->order_yaw[iia] = 3;

  while (ATTITUDE->et_interpo[iia] < ATTITUDE->et_initial){
    iia = iia + 1;

    ATTITUDE->et_interpo[iia] = ATTITUDE->et_initial + 1.0*iia ;  

      ATTITUDE->pitch[iia] = iia * ang_velo[3] * 1.0  + ang_velo[0];
      ATTITUDE->roll[iia] = iia * ang_velo[4] * 1.0  + ang_velo[1];
      ATTITUDE->yaw[iia] = iia * ang_velo[5] * 1.0 + ang_velo[2] ;
/*       ATTITUDE->pitch[iia] = iia * ang_velo[0] * ATTITUDE->dt_prop / 2.; */
/*       ATTITUDE->roll[iia] = iia * ang_velo[1] * ATTITUDE->dt_prop / 2.; */
/*       ATTITUDE->yaw[iia] = iia * ang_velo[2] * ATTITUDE->dt_prop / 2.; */
      ATTITUDE->order_pitch[iia] = 1;
      ATTITUDE->order_roll[iia] = 2;
      ATTITUDE->order_yaw[iia] = 3;

    //        etprint(ATTITUDE->et_interpo[iia], "tle");
	
  } // leave this loop when x_after_interpo gets newer than constellation epoch (ATTITUDE->et_initial)
  j = 0;
  while (iia<ATTITUDE->nb_time_steps ){

    ATTITUDE->et_interpo[iia] = ATTITUDE->et_initial + 1.0*j ;

      ATTITUDE->pitch[iia] = iia * ang_velo[3] * 1.0 + ang_velo[0];
      ATTITUDE->roll[iia] = iia * ang_velo[4] * 1.0  + ang_velo[1];
      ATTITUDE->yaw[iia] = iia * ang_velo[5] * 1.0 + ang_velo[2] ;
/*       ATTITUDE->pitch[iia] = iia * ang_velo[0] * ATTITUDE->dt_prop / 2.; */
/*       ATTITUDE->roll[iia] = iia * ang_velo[1] * ATTITUDE->dt_prop / 2.; */
/*       ATTITUDE->yaw[iia] = iia * ang_velo[2] * ATTITUDE->dt_prop / 2.; */
      ATTITUDE->order_pitch[iia] = 1;
      ATTITUDE->order_roll[iia] = 2;
      ATTITUDE->order_yaw[iia] = 3;

/*       etprint(ATTITUDE->et_interpo[iia], "const"); */
/*       printf("iia %f %f %f\n", ATTITUDE->pitch[iia] , ATTITUDE->roll[iia] , ATTITUDE->yaw[iia] ); */
    //        etprint(ATTITUDE->et_interpo[iia], "const");
	iia = iia +1;
	j = j+1;
  }


/*     for (iia = 0; iia < ATTITUDE->nb_time_steps ; iia++){ */
/*       ATTITUDE->et_interpo[iia] = ATTITUDE->et_initial + ATTITUDE->dt_prop*iia / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*       ATTITUDE->pitch[iia] = iia * ang_velo[3] * ATTITUDE->dt_prop / 2. + ang_velo[0]; */
/*       ATTITUDE->roll[iia] = iia * ang_velo[4] * ATTITUDE->dt_prop / 2. + ang_velo[1]; */
/*       ATTITUDE->yaw[iia] = iia * ang_velo[5] * ATTITUDE->dt_prop / 2. + ang_velo[2] ; */
/* /\*       ATTITUDE->pitch[iia] = iia * ang_velo[0] * ATTITUDE->dt_prop / 2.; *\/ */
/* /\*       ATTITUDE->roll[iia] = iia * ang_velo[1] * ATTITUDE->dt_prop / 2.; *\/ */
/* /\*       ATTITUDE->yaw[iia] = iia * ang_velo[2] * ATTITUDE->dt_prop / 2.; *\/ */
/*       ATTITUDE->order_pitch[iia] = 1; */
/*       ATTITUDE->order_roll[iia] = 2; */
/*       ATTITUDE->order_yaw[iia] = 3; */
/*     } */

  }

  // if the attitude representation is give in a file
  else{


    strcpy(temp_copy, ATTITUDE->attitude_profile);
    //newstructure
/*     strcpy(text_location, ATTITUDE->dir_input_attitude); */
/*     strcat(text_location, "/"); */
    strcpy(text_location, "");
    //newstructure
    strcpy(ATTITUDE->attitude_profile, text_location);
    strcat(ATTITUDE->attitude_profile,temp_copy);

    ATTITUDE->et_interpo = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );

    ATTITUDE->pitch = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );
    ATTITUDE->roll = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );
    ATTITUDE->yaw = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );
    ATTITUDE->order_pitch = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );
    ATTITUDE->order_roll = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );
    ATTITUDE->order_yaw = malloc( ATTITUDE->nb_time_steps  * sizeof(double) );
    ATTITUDE->quaternion = malloc( ATTITUDE->nb_time_steps  * sizeof(double*) );
    for (iia = 0; iia < ATTITUDE->nb_time_steps ; iia++){
    ATTITUDE->quaternion[iia] = malloc( 4 * sizeof(double) );
    }
    if (  ATTITUDE->pitch == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->pitch \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  ATTITUDE->roll == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->roll \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  ATTITUDE->yaw == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->yaw \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }

    if (  ATTITUDE->order_pitch == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->order_pitch \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  ATTITUDE->order_roll == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->order_roll \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }
    if (  ATTITUDE->order_yaw == NULL ){
      printf("***! Could not allow memory space for  ATTITUDE->order_yaw \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }


    if ( ATTITUDE->et_interpo == NULL ){
      printf("***! Could not allow memory space for ATTITUDE->et_interpo  \n. The program will stop. !***\n");
      ierr =  MPI_Finalize();
      exit(0);
    }

    lin_interpolate_attitude(ATTITUDE->quaternion,ATTITUDE->pitch, ATTITUDE->roll, ATTITUDE->yaw, ATTITUDE->order_pitch, ATTITUDE->order_roll, ATTITUDE->order_yaw, ATTITUDE->et_interpo, ATTITUDE ->attitude_profile, ATTITUDE->nb_time_steps , ATTITUDE->initial_epoch, ATTITUDE->et_initial,1.0        , iProc, &ATTITUDE->file_is_quaternion); // 1.0 because inteprolation is every second for the copmutation of the specular points



  }
  //  MPI_Finalize();exit(0);


  return 0;

}




/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           nb_elements_file
//  Purpose:        Returns the number of lines after the header. header_end represents the last line of the file before nb_elements_in_file is incremented
//  Assumptions:    None
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 08/17/2015    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////


int nb_elements_file(int *nb_element_in_file,
		     char filename[256],
		     char header_end[256],
		     char file_end[256]){

  /* Declarations */
  FILE *fp=NULL;
  char *line = NULL;
  size_t len = 0;
  int found_eoh = 0;
  char text[256];
  
  /* Algorithm */
  fp = fopen(filename, "r");
  if (fp == NULL){
    printf("***! Could not find the file %s. The program will stop. !***\n", filename); MPI_Finalize();exit(0);
  }

  // Count number of elements in file
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( text, file_end  ) == 0 )  {
      found_eoh = 1;
    }
    *nb_element_in_file = *nb_element_in_file+1;
  }
  rewind(fp);
  
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( header_end, text  ) == 0 )  {
      found_eoh = 1;
      *nb_element_in_file = *nb_element_in_file-1;
    }
    *nb_element_in_file = *nb_element_in_file-1;
  }

  return 0;
  

}


/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           previous_index
//  Purpose:        Among all the elements of the array x_before_interpo that are lower than value, computes the maximum index 
//  Assumptions:    x_before_interpo is sorted in the ascending order
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 08/18/2015    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////


int previous_index(int *x_min_index, 
		   double *x_before_interpo, 
		   double value,
		   double size_x_before_interpo){

  /* Declarations */
  int i_index = 0;

  /* Algorithm */

  while ( (x_before_interpo[i_index] <= value) && (i_index<size_x_before_interpo) ){
    i_index++;
  }
  *x_min_index = i_index-1;


       

  return 0;

}





/////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           lin_interpolate_attitude
//  Purpose:        Linear interpolation of the attitude file
//  Assumptions:    if the file gives quaternions, there should not be more tahn 5 space character on the first line of the data (this is how SpOCK differentiates it from a pitch roll yaw order-picth order_roll order_yaw file)
//  References      None
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | C. Bussy-Virat| 10/25/2015    |   ---     | Initial Implementation 
//
/////////////////////////////////////////////////////////////////////////////////////////

int lin_interpolate_attitude(double **quaternion_after_interpo,
			     double *pitch_after_interpo,
			     double *roll_after_interpo,
			     double *yaw_after_interpo,
			     int *pitch_order_after_interpo,
			     int *roll_order_after_interpo,
			     int *yaw_order_after_interpo,
			     double *x_after_interpo,
			     char filename[256],
			     int nb_time_steps_simu,
			     char initial_epoch[256], 
			     double et_oldest_tle_epoch,
			     double dt,
			     int iProc,
			     int *file_is_quaternion){

  /* Declarations */
  char read_time[10]; char len_fake_text[10];    char read_all_line[100];
      double **quaternion_before_interpo=NULL;
  char fake_text[256];
  int x_min_index;
  double *pitch_before_interpo = NULL;
  double *roll_before_interpo = NULL;
  double *yaw_before_interpo = NULL;
  double *pitch_order_before_interpo = NULL;
  double *roll_order_before_interpo = NULL;
  double *yaw_order_before_interpo = NULL;
  double *x_before_interpo = NULL;
  FILE *fp = NULL;
  char *line = NULL;
  size_t len = 0;
  int found_eoh = 0;
  char text[256];
  int line_num;
  double x_before_interpo_temp;
  int i;
  int nb_elements_in_file=0; 
  double et_initial;
  double *a = NULL;
  double *b = NULL;
  double x_min, x_max;
  double y_min_quaternion[4], y_max_quaternion[4];
  double y_min_pitch, y_max_pitch;
  double y_min_roll, y_max_roll;
  double y_min_yaw, y_max_yaw;

  /* Algorithm */

  str2et_c(initial_epoch, &et_initial);  
  /* Calculates the x and y arrays to interpolate */
  nb_elements_file(&nb_elements_in_file, filename, "#ENDOFHEADER","#ENDOFFILE");



  x_before_interpo = malloc( nb_elements_in_file * sizeof(double) );
  a = malloc( nb_time_steps_simu * sizeof(double) );
  b = malloc( nb_time_steps_simu * sizeof(double) );

  fp = fopen(filename, "r");

  found_eoh = 0;
  while ( found_eoh == 0 && !feof(fp)) {
    getline(&line, &len, fp);
    sscanf(line, "%s", text);
    if (  strcmp( "#ENDOFHEADER", text  ) == 0 )  {
      found_eoh = 1;
    }
  }

  int found_start = 0; 
  ssize_t read;
      while ( (found_start == 0) && (read = getline(&line, &len, fp)) != -1 ) {
	if (strchr(&line[0],'2') != NULL)
	  found_start = 1;
      }

    sscanf(line, "%s", text);

     int count_bracket=0;
       int  icount_bracket;

      for (icount_bracket = 0; icount_bracket< strlen(line); icount_bracket++){
	if (line[icount_bracket]=='('){
	  count_bracket = count_bracket + 1;
	  }
      }


      //    for (count_bracket=0; text[count_bracket]; text[count_bracket]==' ' ? count_bracket++ : *text++);
    // if more than 5 spaces then the file is pitch roll yaw. otherwise it's quaternion

    if (count_bracket  >= 1){
      *file_is_quaternion = 0;
    }
    else{
      *file_is_quaternion = 1;
    }

    if (*file_is_quaternion == 0) {
  pitch_before_interpo = malloc( nb_elements_in_file * sizeof(double) );
  roll_before_interpo = malloc( nb_elements_in_file * sizeof(double) );
  yaw_before_interpo = malloc( nb_elements_in_file * sizeof(double) );
  pitch_order_before_interpo = malloc( nb_elements_in_file * sizeof(double) );
  roll_order_before_interpo = malloc( nb_elements_in_file * sizeof(double) );
  yaw_order_before_interpo = malloc( nb_elements_in_file * sizeof(double) );

  for (line_num = 0; line_num<nb_elements_in_file; line_num++){
    if (line_num > 0){
    getline(&line, &len, fp);
    }
    strcpy(text,"");
    sscanf(line, "%s", fake_text);
    strcpy(read_time, "%");
    sprintf(len_fake_text, "%lu", strlen(fake_text));
    strcat(read_time, len_fake_text);
    strcat(read_time, "[^\\n] ");
    strcpy(read_all_line, read_time);
    strcat(read_all_line, "(%lf;%lf;%lf) (%lf;%lf;%lf)");
    RemoveSpaces(line);
    sscanf(line, read_all_line, text, &pitch_before_interpo[line_num], &roll_before_interpo[line_num], &yaw_before_interpo[line_num], &pitch_order_before_interpo[line_num], &roll_order_before_interpo[line_num], &yaw_order_before_interpo[line_num]);  
    

    //    sscanf(line, "%17[^\n] (%lf;%lf;%lf) (%lf;%lf;%lf)", text, &pitch_before_interpo[line_num], &roll_before_interpo[line_num], &yaw_before_interpo[line_num], &pitch_order_before_interpo[line_num], &roll_order_before_interpo[line_num], &yaw_order_before_interpo[line_num]);  

    /* print_test(); */
/*     printf("<%s> || (%f %f %f) (%f %f %f)\n", text, pitch_before_interpo[line_num], roll_before_interpo[line_num], yaw_before_interpo[line_num], pitch_order_before_interpo[line_num], roll_order_before_interpo[line_num],  yaw_order_before_interpo[line_num]); */
    /*   exit(0); */

    str2et_c (text, &x_before_interpo_temp );

    x_before_interpo[line_num] = x_before_interpo_temp;
  }

  free(line);
  fclose(fp);


  /* Calculates the x array on which the interpolation is done */

/*   //char times[256]; */
  
    i = 0;
  x_after_interpo[i] = et_oldest_tle_epoch;  // "/ 2.0" because of the Runge Kunta orfer 4 method
  while (x_after_interpo[i] < et_initial){
    i = i + 1;
    x_after_interpo[i] = et_oldest_tle_epoch + dt*i;
    //        etprint(x_after_interpo[i], "tle");

  } // leave this loop when x_after_interpo gets newer than constellation epoch (et_initial)
  int j = 0;
  //   printf("Nd %d\n",nb_time_steps_simu);exit(0);
  while (i<nb_time_steps_simu){
    x_after_interpo[i] = et_initial + dt*j;
       //          etprint(x_after_interpo[i], "const");
	i = i +1;
	j = j+1;
  }
  


  /* for (i = 0; i < nb_time_steps_simu; i++){ */
  /*   x_after_interpo[i] = et_initial + dt*i / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */

  /*   //    et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times); */
  /* } */

  /* Compute y_after_interpo */

  //  FILE *fp_temp;
  //  fp_temp = fopen("att-test.txt", "w+");
  //   char times[256]; 
  /*   char times2[256]; */
  for (i = 1; i < nb_time_steps_simu-1; i++){
    if ( ( x_after_interpo[i] >= x_before_interpo[0] ) && (x_after_interpo[i] <= x_before_interpo[nb_elements_in_file-1] ) ) {
      previous_index(&x_min_index, x_before_interpo, x_after_interpo[i], nb_elements_in_file);
      /*     et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times);  */
      /*     printf("\nx_after_interpo[%d] = %s \n", i, times); */
      /*     et2utc_c(x_before_interpo[x_min_index], "C" ,3 ,255 , times);  */
      /*     et2utc_c(x_before_interpo[x_min_index+1], "C" ,3 ,255 , times2);  */
      /*     printf("x_before_interpo[x_min_index] = %s || %s\n",  times,  times2); */
      x_min = x_before_interpo[x_min_index];
      x_max = x_before_interpo[x_min_index+1];
      // pitch
      y_min_pitch = pitch_before_interpo[x_min_index];
      y_max_pitch = pitch_before_interpo[x_min_index+1];
      a[i- 1] = (y_max_pitch - y_min_pitch) / (x_max - x_min);
      b[i-1] = y_max_pitch - a[i-1]*x_max;
      pitch_after_interpo[i] = a[i-1]*x_after_interpo[i] + b[i-1];
      pitch_order_after_interpo[i] = pitch_order_before_interpo[x_min_index];
      // roll
      y_min_roll = roll_before_interpo[x_min_index];
      y_max_roll = roll_before_interpo[x_min_index+1];
      a[i- 1] = (y_max_roll - y_min_roll) / (x_max - x_min);
      b[i-1] = y_max_roll - a[i-1]*x_max;
      roll_after_interpo[i] = a[i-1]*x_after_interpo[i] + b[i-1];
      roll_order_after_interpo[i] = roll_order_before_interpo[x_min_index];
      // yaw
      y_min_yaw = yaw_before_interpo[x_min_index];
      y_max_yaw = yaw_before_interpo[x_min_index+1];
      a[i- 1] = (y_max_yaw - y_min_yaw) / (x_max - x_min);
      b[i-1] = y_max_yaw - a[i-1]*x_max;
      yaw_after_interpo[i] = a[i-1]*x_after_interpo[i] + b[i-1];
      yaw_order_after_interpo[i] = yaw_order_before_interpo[x_min_index];
    }
    else if ( x_after_interpo[i] < x_before_interpo[0] ){
      pitch_after_interpo[i] = pitch_before_interpo[0];
      roll_after_interpo[i] = roll_before_interpo[0];
      yaw_after_interpo[i] = yaw_before_interpo[0];
      pitch_order_after_interpo[i] = pitch_order_before_interpo[0];
      roll_order_after_interpo[i] = roll_order_before_interpo[0];
      yaw_order_after_interpo[i] = yaw_order_before_interpo[0];
    }

    else{
      pitch_after_interpo[i] = pitch_before_interpo[nb_elements_in_file-1];
      roll_after_interpo[i] = roll_before_interpo[nb_elements_in_file-1];
      yaw_after_interpo[i] = yaw_before_interpo[nb_elements_in_file-1];
      pitch_order_after_interpo[i] = pitch_order_before_interpo[nb_elements_in_file-1];
      roll_order_after_interpo[i] = roll_order_before_interpo[nb_elements_in_file-1];
      yaw_order_after_interpo[i] = yaw_order_before_interpo[nb_elements_in_file-1];
    }
/*     printf("\n"); */
/*     etprint(x_after_interpo[i], ""); */
/* printf("%f %f %f %d %d %d\n", pitch_after_interpo[i], roll_after_interpo[i], yaw_after_interpo[i], pitch_order_after_interpo[i], roll_order_after_interpo[i], yaw_order_after_interpo[i]); */
/*        et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times);   */
/*        fprintf(fp_temp, "%s %f %f %f %d %d %d\n", times,  pitch_after_interpo[i], roll_after_interpo[i], yaw_after_interpo[i], pitch_order_after_interpo[i], roll_order_after_interpo[i], yaw_order_after_interpo[i]); */
  }

  if( fabs( x_after_interpo[0] - x_before_interpo[0] ) > 0.01 ){ // if the first time in the driver file is different (by more than 0.01s here) from the first time of propagation then the first interpolated value of the driver is equal to the second one. 
    pitch_after_interpo[0] = pitch_after_interpo[1];
    roll_after_interpo[0] = roll_after_interpo[1];
    yaw_after_interpo[0] = yaw_after_interpo[1];
    pitch_order_after_interpo[0] = pitch_order_after_interpo[1];
    roll_order_after_interpo[0] = roll_order_after_interpo[1];
    yaw_order_after_interpo[0] = yaw_order_after_interpo[1];

  }
  else{ // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    pitch_after_interpo[0] = pitch_before_interpo[0];
    roll_after_interpo[0] = roll_before_interpo[0];
    yaw_after_interpo[0] = yaw_before_interpo[0];
    pitch_order_after_interpo[0] = pitch_order_before_interpo[0];
    roll_order_after_interpo[0] = roll_order_before_interpo[0];
    yaw_order_after_interpo[0] = yaw_order_before_interpo[0];
  }

  if( fabs( x_after_interpo[nb_time_steps_simu-1] - x_before_interpo[nb_elements_in_file-1] ) > 0.01 ){ // same comments as the two previous blocks
    pitch_after_interpo[nb_time_steps_simu-1] = pitch_after_interpo[nb_time_steps_simu-2];
    roll_after_interpo[nb_time_steps_simu-1] = roll_after_interpo[nb_time_steps_simu-2];
    yaw_after_interpo[nb_time_steps_simu-1] = yaw_after_interpo[nb_time_steps_simu-2];
    pitch_order_after_interpo[nb_time_steps_simu-1] = pitch_order_after_interpo[nb_time_steps_simu-2];
    roll_order_after_interpo[nb_time_steps_simu-1] = roll_order_after_interpo[nb_time_steps_simu-2];
    yaw_order_after_interpo[nb_time_steps_simu-1] = yaw_order_after_interpo[nb_time_steps_simu-2];

  }
  else{ // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    pitch_after_interpo[nb_time_steps_simu-1] = pitch_before_interpo[nb_elements_in_file-1];
    roll_after_interpo[nb_time_steps_simu-1] = roll_before_interpo[nb_elements_in_file-1];
    yaw_after_interpo[nb_time_steps_simu-1] = yaw_before_interpo[nb_elements_in_file-1];
    pitch_order_after_interpo[nb_time_steps_simu-1] = pitch_order_before_interpo[nb_elements_in_file-1];
    roll_order_after_interpo[nb_time_steps_simu-1] = roll_order_before_interpo[nb_elements_in_file-1];
    yaw_order_after_interpo[nb_time_steps_simu-1] = yaw_order_before_interpo[nb_elements_in_file-1];
  }




  free(pitch_before_interpo);
  free(roll_before_interpo);
  free(yaw_before_interpo);
  free(pitch_order_before_interpo);
  free(roll_order_before_interpo);
  free(yaw_order_before_interpo);
    }


    else{

      quaternion_before_interpo = malloc( nb_elements_in_file * sizeof(double*) );

      for (line_num = 0; line_num<nb_elements_in_file; line_num++){
	quaternion_before_interpo[line_num] = malloc( 4 * sizeof(double) );
	if (line_num > 0){
	  getline(&line, &len, fp);
	}

	sscanf(line, "%s %lf %lf %lf %lf", text,&quaternion_before_interpo[line_num][0], &quaternion_before_interpo[line_num][1], &quaternion_before_interpo[line_num][2], &quaternion_before_interpo[line_num][3]);
      

      str2et_c(text, &x_before_interpo[line_num]);  
      
      //           etprint(x_before_interpo[line_num], "");
/*       printf("%f %f %f %f\n",quaternion_before_interpo[line_num][0], quaternion_before_interpo[line_num][1], quaternion_before_interpo[line_num][2], quaternion_before_interpo[line_num][3]); */
/*       printf("%d %d\n", line_num, nb_elements_in_file); */
    }


      //      exit(0);


  /* Calculates the x array on which the interpolation is done */

/*   //char times[256]; */
  
    i = 0;
  x_after_interpo[i] = et_oldest_tle_epoch;  // "/ 2.0" because of the Runge Kunta orfer 4 method
  while (x_after_interpo[i] < et_initial){
    i = i + 1;
    x_after_interpo[i] = et_oldest_tle_epoch + dt*i;
    //        etprint(x_after_interpo[i], "tle");

  } // leave this loop when x_after_interpo gets newer than constellation epoch (et_initial)
  int j = 0;
  //   printf("Nd %d\n",nb_time_steps_simu);exit(0);
  while (i<nb_time_steps_simu){
    x_after_interpo[i] = et_initial + dt*j;
       //          etprint(x_after_interpo[i], "const");
	i = i +1;
	j = j+1;
  }
  



/*   i = 0; */
/*   x_after_interpo[i] = et_oldest_tle_epoch;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*   while (x_after_interpo[i] < et_initial){ */
/*     i = i + 1; */
/*     x_after_interpo[i] = et_oldest_tle_epoch + dt*i / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*     //        etprint(x_after_interpo[i], "tle"); */
	
/*   } // leave this loop when x_after_interpo gets newer than constellation epoch (et_initial) */
/*   int j = 0; */
/*   while (i<nb_time_steps_simu){ */
/*     x_after_interpo[i] = et_initial + dt*j / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */
/*     //        etprint(x_after_interpo[i], "const"); */
/* 	i = i +1; */
/* 	j = j+1; */
/*   } */


/*   for (i = 0; i < nb_time_steps_simu; i++){ */
/*     x_after_interpo[i] = et_initial + dt*i / 2.0;  // "/ 2.0" because of the Runge Kunta orfer 4 method */

/*     //    et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times); */
/*   } */


  for (i = 1; i < nb_time_steps_simu-1; i++){
    if ( ( x_after_interpo[i] >= x_before_interpo[0] ) && (x_after_interpo[i] <= x_before_interpo[nb_elements_in_file-1] ) ) {
      previous_index(&x_min_index, x_before_interpo, x_after_interpo[i], nb_elements_in_file);
      /*     et2utc_c(x_after_interpo[i], "C" ,3 ,255 , times);  */
      /*     printf("\nx_after_interpo[%d] = %s \n", i, times); */
      /*     et2utc_c(x_before_interpo[x_min_index], "C" ,3 ,255 , times);  */
      /*     et2utc_c(x_before_interpo[x_min_index+1], "C" ,3 ,255 , times2);  */
      /*     printf("x_before_interpo[x_min_index] = %s || %s\n",  times,  times2); */
      x_min = x_before_interpo[x_min_index];
      x_max = x_before_interpo[x_min_index+1];
      // quaternion
      y_min_quaternion[0] = quaternion_before_interpo[x_min_index][0];
      y_min_quaternion[1] = quaternion_before_interpo[x_min_index][1];
      y_min_quaternion[2] = quaternion_before_interpo[x_min_index][2];
      y_min_quaternion[3] = quaternion_before_interpo[x_min_index][3];

      y_max_quaternion[0] = quaternion_before_interpo[x_min_index+1][0];
      y_max_quaternion[1] = quaternion_before_interpo[x_min_index+1][1];
      y_max_quaternion[2] = quaternion_before_interpo[x_min_index+1][2];
      y_max_quaternion[3] = quaternion_before_interpo[x_min_index+1][3];


      a[i- 1] = (y_max_quaternion[0] - y_min_quaternion[0]) / (x_max - x_min);
      b[i-1] = y_max_quaternion[0] - a[i-1]*x_max;
      quaternion_after_interpo[i][0] = a[i-1]*x_after_interpo[i] + b[i-1];

      a[i- 1] = (y_max_quaternion[1] - y_min_quaternion[1]) / (x_max - x_min);
      b[i-1] = y_max_quaternion[1] - a[i-1]*x_max;
      quaternion_after_interpo[i][1] = a[i-1]*x_after_interpo[i] + b[i-1];

      a[i- 1] = (y_max_quaternion[2] - y_min_quaternion[2]) / (x_max - x_min);
      b[i-1] = y_max_quaternion[2] - a[i-1]*x_max;
      quaternion_after_interpo[i][2] = a[i-1]*x_after_interpo[i] + b[i-1];

      a[i- 1] = (y_max_quaternion[3] - y_min_quaternion[3]) / (x_max - x_min);
      b[i-1] = y_max_quaternion[3] - a[i-1]*x_max;
      quaternion_after_interpo[i][3] = a[i-1]*x_after_interpo[i] + b[i-1];

      etprint(x_after_interpo[i],"");
      printf("%f %f %f %f\n\n", quaternion_after_interpo[i][0], quaternion_after_interpo[i][1], quaternion_after_interpo[i][2], quaternion_after_interpo[i][3]);
    }
    else if ( x_after_interpo[i] < x_before_interpo[0] ){
      quaternion_after_interpo[i][0] = quaternion_before_interpo[0][0];
      quaternion_after_interpo[i][1] = quaternion_before_interpo[0][1];
      quaternion_after_interpo[i][2] = quaternion_before_interpo[0][2];
      quaternion_after_interpo[i][3] = quaternion_before_interpo[0][3];
    }

    else{
      quaternion_after_interpo[i][0] = quaternion_before_interpo[nb_elements_in_file-1][0];
      quaternion_after_interpo[i][1] = quaternion_before_interpo[nb_elements_in_file-1][1];
      quaternion_after_interpo[i][2] = quaternion_before_interpo[nb_elements_in_file-1][2];
      quaternion_after_interpo[i][3] = quaternion_before_interpo[nb_elements_in_file-1][3];
    }


       //       fprintf(fp_temp, "%s %f %f %f %d %d %d\n", times,  quaternion_after_interpo[i], roll_after_interpo[i], yaw_after_interpo[i], quaternion_order_after_interpo[i], roll_order_after_interpo[i], yaw_order_after_interpo[i]);
  }
  //  exit(0);
  if( fabs( x_after_interpo[0] - x_before_interpo[0] ) > 0.01 ){ // if the first time in the driver file is different (by more than 0.01s here) from the first time of propagation then the first interpolated value of the driver is equal to the second one. 
    quaternion_after_interpo[0][0] = quaternion_after_interpo[1][0];
    quaternion_after_interpo[0][1] = quaternion_after_interpo[1][1];
    quaternion_after_interpo[0][2] = quaternion_after_interpo[1][2];
    quaternion_after_interpo[0][3] = quaternion_after_interpo[1][3];
  }
  else{ // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    quaternion_after_interpo[0][0] = quaternion_before_interpo[0][0];
    quaternion_after_interpo[0][1] = quaternion_before_interpo[0][1];
    quaternion_after_interpo[0][2] = quaternion_before_interpo[0][2];
    quaternion_after_interpo[0][3] = quaternion_before_interpo[0][3];
  }

  if( fabs( x_after_interpo[nb_time_steps_simu-1] - x_before_interpo[nb_elements_in_file-1] ) > 0.01 ){ // same comments as the two previous blocks
    quaternion_after_interpo[nb_time_steps_simu-1][0] = quaternion_after_interpo[nb_time_steps_simu-2][0];
    quaternion_after_interpo[nb_time_steps_simu-1][1] = quaternion_after_interpo[nb_time_steps_simu-2][1];
    quaternion_after_interpo[nb_time_steps_simu-1][2] = quaternion_after_interpo[nb_time_steps_simu-2][2];
    quaternion_after_interpo[nb_time_steps_simu-1][3] = quaternion_after_interpo[nb_time_steps_simu-2][3];
  }
  else{ // if both first times are the same then the first interpolated value of the driver is directly equal to the one on the original driver file
    quaternion_after_interpo[nb_time_steps_simu-1][0] = quaternion_before_interpo[nb_elements_in_file-1][0];
    quaternion_after_interpo[nb_time_steps_simu-1][1] = quaternion_before_interpo[nb_elements_in_file-1][1];
    quaternion_after_interpo[nb_time_steps_simu-1][2] = quaternion_before_interpo[nb_elements_in_file-1][2];
    quaternion_after_interpo[nb_time_steps_simu-1][3] = quaternion_before_interpo[nb_elements_in_file-1][3];
  }





      free(line);
      fclose(fp);
    }
    
    free(a);
    free(b);

    free(quaternion_before_interpo);
  free(x_before_interpo);

  //   fclose(fp_temp);



  return 0;

}



//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

int main(int argc, char *argv[]) {
    ATTITUDE_T        ATTITUDE;
   char sat_nb_str[10];
  int err, i, lonshift, iUseMin, LimitToOneDay = 0;
  float rotate=0, yawoverall=0;
  char filename[150];
  char CygnssFilename[256], GpsFilename[256], YawFileName[150];
  float dt;
  double rCygnss, rGps;
  strcpy(YawFileName, "");
  iDebug = 0;

  lonshift = 0;
  iUseMin = 0;

  struct SatelliteDetails *GPS;
  GPS = malloc(nGps * sizeof(struct SatelliteDetails));
  if (GPS == NULL){
    printf("***! (find_specular_points.c)(main) Not enough memory for GPS. The program wills stop. !***\n");
    MPI_Finalize();
    exit(0);
  }

  struct SatelliteDetails *CYGNSS;
  CYGNSS = malloc(nCygnss * sizeof(struct SatelliteDetails));
  if (CYGNSS == NULL){
    printf("***! (find_specular_points.c)(main) Not enough memory for CYGNSS. The program wills stop. !***\n");
    MPI_Finalize();
    exit(0);
  }


 CYGNSS[0].nYawTimes = 0;

 MPI_Init(&argc, &argv);
 MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
 MPI_Comm_rank(MPI_COMM_WORLD, &iProc);

 if (iProc == 0) printf("Computing the specular point locations...\n");

 // CBV 05/08/16
 char filename_input_propagator[150];
 char filename_without_extension[150];

 //newstructure
 // strcpy(filename_input_propagator, "./input/main_input/";)
 strcpy(filename_input_propagator, "");
 //newstructure
 strcat(filename_input_propagator,argv[1]);
    char filename_input_no_path[256];
    strcpy(filename_input_no_path, "");
    strcat(filename_input_no_path, argv[1]);

 FILE *file_input_propagator;
 file_input_propagator = fopen(filename_input_propagator, "r");
 int nCygnss_prop, nGps_prop = 0;
 int found_eoh ; char *line = NULL; size_t len = 0; char text[256];
 found_eoh = 0;

 while ( (found_eoh == 0 ) && (!feof(file_input_propagator)) ) {
   getline(&line, &len, file_input_propagator);
   sscanf(line, "%s", text);
   if (  strcmp( "#SPACECRAFT", text  ) == 0 )  {
     found_eoh = 1;
   }
 }

 
 if (feof(file_input_propagator)){
   printf("***! No section #SPACECRAFT was found in %s. The program will stop.***!\n", filename_input_propagator); MPI_Finalize(); exit(0);
 }
 // Number of CYGNSS
 getline(&line,&len,file_input_propagator);
 sscanf(line,"%d",&nCygnss_prop);
 //  Number of GPS: read the TLE file
 char gps_tle_filename[256]; FILE *gps_tle_file;ssize_t read;
 getline(&line,&len,file_input_propagator);
 sscanf(line,"%s",text);
 if ( strcmp(text, "0") != 0 ){
   //newstructure
   //   strcpy(gps_tle_filename, "./input/tle/constellation_gps_tle/");
   strcpy(gps_tle_filename, "");
   //newstructure
   strcat(gps_tle_filename,text);
   gps_tle_file = fopen(gps_tle_filename, "r");
   while ( (read = getline(&line, &len, gps_tle_file)) != -1 ){
     if (line[0] == '1'){
       nGps_prop = nGps_prop + 1 ;
     }
   }

  // Now 
  rewind(gps_tle_file);
  int first_space;
  int j;
  for (i = 0; i<nGps_prop; i++){
    first_space = 1;
    getline(&line, &len, gps_tle_file);
    sscanf(line, "%255[^\n]s", text);  
    getline(&line, &len, gps_tle_file);
    getline(&line, &len, gps_tle_file);
    strcpy(GPS[i].name, "");
    j = 0;
    while ((unsigned)j < strlen(text)){
      if ( (text[j] == ' ') && (first_space == 0) )  {
	//	strcat(GPS[i],".txt");
	j = 1000;
      }
      if (j < 999)
	GPS[i].name[j] = text[j];
      if ( (text[j] == ' ') && first_space == 1 )  {
	GPS[i].name[j] = '_';
	first_space = 0;
      }
      j = j+1;

    }
    //        printf("<%s>\n", GPS[i].name);
  }

   fclose(gps_tle_file);
 }
 else{
   nGps_prop = 0;
 }
 // MPI_Finalize();exit(0);
 //newstructure
 //  SPICE path
 int sss; char text2[256];int find_file_name;char *next; int find_file_name2;
 char path_to_spice[256];

 found_eoh = 0;
 rewind(file_input_propagator);
 while ( (found_eoh == 0 ) && (!feof(file_input_propagator)) ) {
   getline(&line, &len, file_input_propagator);
   sscanf(line, "%s", text);
   if (  strcmp( "#SPICE", text  ) == 0 )  {
     found_eoh = 1;
   }
 }
 if (feof(file_input_propagator)){
   printf("***! No section #SPICE was found in %s. The program will stop.***!\n", filename_input_propagator); MPI_Finalize(); exit(0);
 }

  getline(&line, &len, file_input_propagator);

  strcpy(path_to_spice, "");
  sscanf(line, "%s", path_to_spice);// ASSUMPTION: the path of spice specified in the main inoput file in the section #SPICE must be the same as the path of SPICE installation in the Makefile, with /data at the end. So if in the Makefile the SPICE directory is hello/hi/ then the path of spice in the main input file in section #SPICE must be hello/hi/data/. This is in this path that the input files of find_specular_points (antenna.bin, antenna.info, sigma0_table.bin)have to be (so in this example in hello/hi/data/)

  strcat(path_to_spice, "/"); // don't put the '/' at the end of the path in section #SPICE

 //newstructure




//newstructure

  char            leap_sec_file[256]; // SPICE leap seconds file
  char leap_sec_file_path[256];


  // // Now load the leap second file. ASSUMPTION: it has to be in directory_to_cspice/data/ (where directory_to_cspice is indicated in the Makefile at SPICE_DIR)
  strcpy(leap_sec_file_path,path_to_spice);
  strcat(leap_sec_file_path, "naif0012.tls");
  strcpy(leap_sec_file, leap_sec_file_path );
  furnsh_c( leap_sec_file);
  // end of sectio SPICE until we get back to it after reading section #TIME
//newstructure



  found_eoh = 0;
  rewind(file_input_propagator);
  while ( found_eoh == 0 && !feof(file_input_propagator)) {
    getline(&line, &len, file_input_propagator);
    sscanf(line, "%s", text);
    if (  strcmp( "#TIME", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  if (feof(file_input_propagator)){
    printf("***! No section #TIME found in %s. The program will stop. !***\n", filename);
   MPI_Finalize();
    exit(0);

  }

  // Date of initial epoch
  getline(&line, &len, file_input_propagator);
  sscanf(line, "%s", text);
  

    strcpy(ATTITUDE.initial_epoch, text );
    // Date of final epoch
    getline(&line, &len, file_input_propagator);
    sscanf(line, "%s", text);
    strcpy(ATTITUDE.final_epoch, text );


  // Time step of simulation
  getline(&line, &len, file_input_propagator);
  RemoveSpaces(line);
  sscanf( line, "%lf", &ATTITUDE.dt_prop);

  // Re-evaluate final epoch and ATTITUDE.et_final so that final epoch is a multiple of ATTITUDE.dt_prop + inital epoch (if originally it is not, then take the closest multiple just before the orginal value of final epoch (ex: if initial epoch is at 12:00:00, ATTITUDE.dt_prop = 20s, and final epoch as written in main input file is 12:00:50 then re-evaluate final epoch to be at 12:00:40))
  
  str2et_c(ATTITUDE.final_epoch, &ATTITUDE.et_final);
  str2et_c(ATTITUDE.initial_epoch, &ATTITUDE.et_initial);
  //  printf("aaaaa <%s> %f\n", ATTITUDE.initial_epoch, ATTITUDE.et_initial);

/*   //  if ( fabs(ATTITUDE.dt_prop -  fmod( ATTITUDE.et_final - ATTITUDE.et_initial, ATTITUDE.dt_prop ) ) > 0.01 ){ // 0.01 for numerical reasons */
/*   if ( fabs(fmod( ATTITUDE.et_final - ATTITUDE.et_initial, ATTITUDE.dt_prop ) ) > 1e-6 ){ // 0.01 for numerical reasons */
/*     /\* printf("X = %f\n", fabs(ATTITUDE.dt_prop - fmod( ATTITUDE.et_final - ATTITUDE.et_initial, ATTITUDE.dt_prop ) )); *\/ */
/*     /\* printf("old: <%s>\n",ATTITUDE.final_epoch); *\/ */
/*     ATTITUDE.et_final = (int) ( (ATTITUDE.et_final - ATTITUDE.et_initial) / ATTITUDE.dt_prop ) * ATTITUDE.dt_prop + ATTITUDE.et_initial; */
/*         et2utc_c(ATTITUDE.et_final, "ISOC", 6, 300, ATTITUDE.final_epoch); */

/*     if (iProc == 0){ */
/*       //      printf("***! The final epoch has been reset to %s (so it is equal to initial epoch + a multiple of ATTITUDE.dt_prop). !***\n", ATTITUDE.final_epoch); */
/*       } */
/*     //    printf("new: <%s>\n",ATTITUDE.final_epoch); */
/* } */

  //    etprint(ATTITUDE.et_final   , "new again");

  
  /* ATTITUDE */


  rewind(file_input_propagator);
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(file_input_propagator)) {
    getline(&line, &len, file_input_propagator);
    sscanf(line, "%s", text);
    if (  strcmp( "#ATTITUDE", text  ) == 0 )  {
      found_eoh = 1;
    }
  }


  double ang_velo[6];\
  //  printf("ATTITUDE.et_initial %f\n", ATTITUDE.et_initial);
/*     etprint(ATTITUDE.et_initial, "ATTITUDE.et_initial"); */
/*     printf("final %s\n", ATTITUDE.final_epoch); */
/*     printf("ATTITUDE.dt_prop %f\n", ATTITUDE.dt_prop); */
    nb_time_steps_f(&ATTITUDE.nb_time_steps, ATTITUDE.et_initial, ATTITUDE.final_epoch , 1.0);// time step of 1 second for the interpolation so here number of times teps is the number of time steps in the interpolation
/*         printf("ATTITUDE.nb_time_steps %d\n", ATTITUDE.nb_time_steps); */
/*               MPI_Finalize();exit(0); */


  if (feof(file_input_propagator)){ // if the user does not include a section #ATTITUDE, then the attitude is set to nadir pointing by default
    strcpy(ATTITUDE.attitude_profile, "nadir");
     load_attitude( &ATTITUDE, file_input_propagator, nProcs, ang_velo, iProc );
  }
  else{
    getline(&line,&len,file_input_propagator);
    sscanf(line,"%s", text);
    strcpy(ATTITUDE.attitude_profile, text);


    // if the user put a vector to represent the angular velocity of the satellite
    if (ATTITUDE.attitude_profile[0] == '('){
    RemoveSpaces(line);  strtok(line, "\n");  strtok(line, "\r"); 
      sscanf(line, "(%lf;%lf;%lf)(%lf;%lf;%lf)", &ang_velo[0], &ang_velo[1], &ang_velo[2], &ang_velo[3], &ang_velo[4], &ang_velo[5]); // first 3 is the inital atitude, last 3 is the (constant) angular velocity. in degree and degree per sec
      //      sscanf(line, "(%lf;%lf;%lf)", &ang_velo[0], &ang_velo[1], &ang_velo[2]); // in degree per sec
      strcpy(ATTITUDE.attitude_profile, "angular_velocity");
    }

 load_attitude( &ATTITUDE, file_input_propagator, nProcs, ang_velo, iProc );


  }



 // NAME OF OUTPUT FILES
 //newstructure
 // int sss; char text2[256];int find_file_name;char *next; int find_file_name2;
 //newstructure
    /* OUTPUT */
    char dir_output_run_name_temp[1000], path_output_run_name_temp[100];
      char text_temp[256];
        struct passwd *pw = getpwuid(getuid());
        const char *homedir = pw->pw_dir;
	    double et_text_output_dir_now;
  char text_output_dir_now[256];
  int ccc;

  rewind(file_input_propagator);
  found_eoh = 0;
  while ( found_eoh == 0 && !feof(file_input_propagator)) {
    getline(&line, &len, file_input_propagator);
    sscanf(line, "%s", text);
    if (  strcmp( "#OUTPUT", text  ) == 0 )  {
      found_eoh = 1;
    }
  }
  if (feof(file_input_propagator)){
    printf("***! No section #OUTPUT found in %s. The program will stop. !***\n", filename);
    MPI_Finalize();
    exit(0);
  }

  // Names of output files
  getline(&line,&len,file_input_propagator);
        RemoveSpaces(line);  strtok(line, "\n");  strtok(line, "\r");
  sscanf(line, "%s", text);
  if (text[(int)(strlen(text))-1] == '/'){ // if the user wrote a '/' at the very end then remove it
    strcpy(text_temp, text);
    strcpy(text,"");
      strncat(text, &text_temp[0], (int)(strlen(text_temp))-1);
  }

  next = &text[0];

  int find_tild=0;
    find_tild = (int)(strchr(next, '~') - next);
  if (text[0] == '~'){
    strcpy(text_temp, homedir);
    strcat(text_temp, "/");
    strncat(text_temp,next + find_tild+2, next + (int)(strlen(text)));
    strcpy(text,text_temp);
  next = &text[0];
    }

    int find_dir=0;
    find_dir = (int)(strrchr(next, '/') - next);
    if ((find_dir < 0) || (find_dir > 10000)){// no / in the name so the user wants the otuput folder in the directory where he runs the command
    strcpy(text_temp, text);
    strcpy(text,"./");
      strncat(text, &text_temp[0], (int)(strlen(text_temp)));
    find_dir = (int)(strrchr(next, '/') - next);
    //      find_dir = -1;

    }

    strcpy(path_output_run_name_temp, "");
    strncat(path_output_run_name_temp, next, find_dir +1);
    //strncat(text_surface, next, find_file_name);
    strcpy(dir_output_run_name_temp, "");
    strncat(dir_output_run_name_temp, next + find_dir +1, next + (int)(strlen(text)));

    /* printf("<%s>\n",path_output_run_name_temp); */

 
    if (strcmp(dir_output_run_name_temp, "out") == 0 ){
      strcpy(filename_without_extension, "");
      strncat(filename_without_extension, &filename_input_no_path[0], (int)(strrchr(&filename_input_no_path[0], '.') - &filename_input_no_path[0]));
      strcpy(dir_output_run_name_temp, filename_without_extension);

    }
/*     printf("<%s>\n",dir_output_run_name_temp); */
/*     MPI_Finalize();exit(0); */


        /*     printf("<%s>\n",path_output_run_name_temp); */
        /* printf("<%s>\n",dir_output_run_name_temp); */
        /* MPI_Finalize();exit(0); */
         
  for (sss = 0; sss<nCygnss_prop; sss++){
    if (sss == 0){
      // CYGNSS CONSTELLATION FILE
      //newstructure
      //      strcpy(CygnssFilename, "");
      //     strcpy(CygnssFilename, "./output/");
     strcpy(CygnssFilename, path_output_run_name_temp);
     //newstructure
          strcat(CygnssFilename, dir_output_run_name_temp);
     strcat(CygnssFilename, "/CONSTELLATION_CYGNSS_for_run_");
     strcat(CygnssFilename, dir_output_run_name_temp);
     strcat(CygnssFilename, "1.txt");
     // GPS CONSTELLATION FILE
     strcpy(GpsFilename, path_output_run_name_temp);

          strcat(GpsFilename, dir_output_run_name_temp);
     strcat(GpsFilename, "/CONSTELLATION_GPS_for_run_");
     strcat(GpsFilename, dir_output_run_name_temp);
     strcat(GpsFilename, "1.txt");
     // Names for the CYGNSS and GPS satellites
     strcpy(CYGNSS[sss].name, dir_output_run_name_temp);
      sprintf(sat_nb_str, "%d", sss+1);
      strcat(CYGNSS[sss].name, sat_nb_str);
    }
    else{
      strcpy(CYGNSS[sss].name, dir_output_run_name_temp);
      sprintf(sat_nb_str, "%d", sss+1);
      strcat(CYGNSS[sss].name, sat_nb_str);
    }
    //    printf("CYGNSS[%d].name = <%s>\n",sss, CYGNSS[sss].name);
  }


 fclose(file_input_propagator);


  /* for (sss = 0; sss<OPTIONS->n_satellites; sss++){ */
  /*   if ( (OPTIONS->nb_gps > 0 ) &&  sss >= OPTIONS->nb_satellites_not_including_gps ) { */
  /*     strcpy( OPTIONS->filename_output[sss], OPTIONS->gps_file_name[sss-OPTIONS->nb_satellites_not_including_gps]); */
  /*     strcpy(OPTIONS->name_sat[sss],OPTIONS->filename_output[sss]); */
  /*     strcat( OPTIONS->filename_output[sss],".txt" ); */
  /*   } */
  /*   else{ */
  /*     strcpy(OPTIONS->filename_output[sss], dir_output_run_name_temp); */
  /*     sprintf(sat_nb_str, "%d", sss+1); */
  /*     strcat(OPTIONS->filename_output[sss], sat_nb_str); */
  /*     strcpy(OPTIONS->name_sat[sss],OPTIONS->filename_output[sss]); */
  /*     strcat(OPTIONS->filename_output[sss], ".txt"); */
  /*   } */
  /* } */

/*  rewind(file_input_propagator); */
/*  while ( (found_eoh == 0 ) && (!feof(file_input_propagator)) ) { */
/*    getline(&line, &len, file_input_propagator); */
/*    sscanf(line, "%s", text); */
/*    if (  strcmp( "#OUTPUT", text  ) == 0 )  { */
/*      found_eoh = 1; */
/*    } */
/*  } */
/*  if (feof(file_input_propagator)){ */
/*    printf("***! No section #OUTPUT was found in %s. The program will stop.***!\n", filename_input_propagator); MPI_Finalize(); exit(0); */
/*  } */
/*  getline(&line,&len,file_input_propagator); */
/*  sscanf(line, "%s", text); */
/*      if (strcmp(text, "out") == 0){ */
/*       strcpy(filename_without_extension, ""); */
/*       strncat(filename_without_extension, &filename_input_no_path[0], (int)(strchr(&filename_input_no_path[0], '.') - &filename_input_no_path[0])); */
/*        strcpy(text, filename_without_extension); */
/*      } */
     
/*   for (sss = 0; sss<nCygnss_prop; sss++){ */
/*     if (sss == 0){ */
/*       // CYGNSS CONSTELLATION FILE */
/*       //newstructure */
/*       //      strcpy(CygnssFilename, ""); */
/*       //     strcpy(CygnssFilename, "./output/"); */
/*      strcpy(CygnssFilename, "./spock/"); */
/*      //newstructure */
/*      strcat(CygnssFilename, text); */

/*      strcat(CygnssFilename, "/CONSTELLATION_CYGNSS_for_run_"); */
/*      strcat(CygnssFilename, text); */
/*      strcat(CygnssFilename, "1.txt"); */
/*      // GPS CONSTELLATION FILE */
/* /\*       strcpy(GpsFilename, ""); *\/ */
/* /\*      strcpy(GpsFilename, "./output/"); *\/ */
/*      strcpy(GpsFilename, "./spock/"); */
/*      //newstructure */
/*      strcat(GpsFilename, text); */
/*      strcat(GpsFilename, "/CONSTELLATION_GPS_for_run_"); */
/*      strcat(GpsFilename, text); */
/*      strcat(GpsFilename, "1.txt"); */
     
/*      // Names for the CYGNSS and GPS satellites */
/*      strcpy(CYGNSS[sss].name, ""); */
/*      strcpy(CYGNSS[sss].name, text); */
/*       sprintf(sat_nb_str, "%d", sss+1); */
/*       strcat(CYGNSS[sss].name, sat_nb_str); */
/*     } */
/*     else{ */
/*       strcpy(CYGNSS[sss].name, ""); */
/*       strcpy(CYGNSS[sss].name, text); */
/*       sprintf(sat_nb_str, "%d", sss+1); */
/*       strcat(CYGNSS[sss].name, sat_nb_str); */
/*     } */
/*     //    printf("CYGNSS[%d].name = <%s>\n",sss, CYGNSS[sss].name); */
/*   } */

/*  fclose(file_input_propagator); */

 /* for (sss = 0; sss < nCygnss_prop; sss++){ */
 /*   if (sss == 0){ */
 /*     next = &line[0]; */
 /*     find_file_name = (int)(strchr(next, '.') - next); */
 /*     strcpy(text, ""); */
 /*     strncat(text, next, find_file_name); */
 /*     strcpy(CygnssFilename, ""); */
 /*     strcpy(CygnssFilename, "./output/run_"); */
 /*     strcat(CygnssFilename, text); */
 /*     strcat(CygnssFilename, "/CONSTELLATION_CYGNSS_for_run_"); */
 /*     strcat(CygnssFilename, text); */
 /*     strcat(CygnssFilename, ".txt"); */
 /*     strcpy(GpsFilename, ""); */
 /*     strcpy(GpsFilename, "./output/run_"); */
 /*     strcat(GpsFilename, text); */
 /*     strcat(GpsFilename, "/CONSTELLATION_GPS_for_run_"); */
 /*     strcat(GpsFilename, text); */
 /*     strcat(GpsFilename, ".txt"); */
 /*     strcpy(CYGNSS[sss].name, ""); */
 /*     strcpy(CYGNSS[sss].name, text); */
 /*   } */
  /*  else{ */
 /*     next =  next + find_file_name + 1; */
 /*     find_file_name2 = (int)(strchr(next, ',') - next); */
 /*     find_file_name = (int)(strchr(next, '.') - next); */
 /*     strcpy(text2, ""); */
 /*     strncat(text2, next+find_file_name2+1, find_file_name - find_file_name2-1); */
 /*     RemoveSpaces(text2); */
 /*     strtok(text2, "\n"); */
 /*     strcpy(CYGNSS[sss].name, ""); */
 /*     strcpy(CYGNSS[sss].name, text2); */
 /*   } */
 /* } */
  // fclose(file_input_propagator);

 // END CBV 05/08/16


  for (i=1; i< argc; i++) {
    //printf("arg%d=%s\n", i, argv[i]);
    if (strstr(argv[i],"-lon=") != NULL) {
      lonshift = atoi(strmid(argv[i],5,3));
      //printf("lonshift : %d\n",lonshift);
    }
    if (strstr(argv[i],"-rot=") != NULL) {
      rotate = (float) atoi(strmid(argv[i],5,3));
      //printf("rotate : %f\n",rotate);
    }
    if (strstr(argv[i],"-yaw=") != NULL) {

      yawoverall = (float) atoi(strmid(argv[i],5,3));
      //printf("yawoverall : %f\n",yawoverall);
    }
    if (strstr(argv[i],"-nspm=") != NULL) {
      nSpecularPointsMax =  atoi(strmid(argv[i],6,3));
      printf("nSpecularPointsMax : %d\n",nSpecularPointsMax);
    }
    if (strstr(argv[i],"-yawfile=") != NULL) {
      strcpy(YawFileName,strmid(argv[i],9,71));
      printf("YawFile : %s\n",YawFileName);
      err = read_yaw_file(YawFileName, CYGNSS);
      if (err != 0) {
	printf("Error! Must Stop!\n");
	//newstructure
// return;
//newstructure
return 0;
      }
    }
    if (strstr(argv[i],"-day") != NULL) {
      LimitToOneDay = 1;
      printf("LimitToOneDay : %d\n",LimitToOneDay);
    }
    if (strstr(argv[i],"-min") != NULL) iUseMin = 1;
  }

  if (argc < 4) {

 // CBV 05/08/16
    if (iProc == 0) {
      printf("Usage: \n");
      printf("mpirun -np N a.out input_file_propagator\n");
      printf("                   -lon=XXX (longitude shift to add to each Cygnss lon)\n");
      printf("                   -rot=XXX (rotation of the satellite)\n");
      printf("                   -min (output minimal amounts of data for coverage only)\n");
    }
    /* if (iProc == 0) { */
    /*   printf("Usage: \n"); */
    /*   printf("mpirun -np N a.out CygnssOrbitSTK.txt GpsOrbitSTK.txt OutputFile\n"); */
    /*   printf("                   -lon=XXX (longitude shift to add to each Cygnss lon)\n"); */
    /*   printf("                   -rot=XXX (rotation of the satellite)\n"); */
    /*   printf("                   -min (output minimal amounts of data for coverage only)\n"); */
    /*   printf("                   -day (limit simulation to only one day)\n"); */
    /* } */
 // END CBV 05/08/16
  } else {

    for (i=0; i<nGps_prop; i++) GPS[i].IsGood = 0;
    for (i=0; i<nCygnss_prop; i++) CYGNSS[i].IsGood = 0;

    //newstructure
/*     strcpy(AntennaInfo.AntennaFileName,"./input/specular_points/antenna.bin"); */
/*     strcpy(AntennaInfo.AntennaInfoFileName,"./input/specular_points/antenna.info"); */
/*     strcpy(AntennaInfo.Sigma0FileName,"./input/specular_points/sigma0_table.bin"); */
    strcpy(AntennaInfo.AntennaFileNameStarboard, path_to_spice); // cbv new antenna onboard. before: strcpy(AntennaInfo.AntennaFileName, path_to_spice);
    strcpy(AntennaInfo.AntennaFileNamePort, path_to_spice); // cbv new antenna onboard. before: nothing
    strcpy(AntennaInfo.AntennaInfoFileName, path_to_spice);
    strcpy(AntennaInfo.Sigma0FileName, path_to_spice);

    strcat(AntennaInfo.AntennaFileNameStarboard,"ant_1_starboard_ddmi_v1.agm"); // cbv new antenna onboard. before: strcat(AntennaInfo.AntennaFileName,"antenna.bin");
    strcat(AntennaInfo.AntennaFileNamePort,"ant_1_port_ddmi_v1.agm"); // cbv new antenna onboard. before: nothing
    strcat(AntennaInfo.AntennaInfoFileName,"antenna.info");
    strcat(AntennaInfo.Sigma0FileName,"sigma0_table.bin");

    //newstructure

    err = read_antenna_info();

    //        printf("Reading CYGNSS...\n");
    //    strcpy(CygnssFilename,argv[1]); // CBV 05/08/16


    err = read_stk(CygnssFilename, CYGNSS);

    if (err != 0 ) {
      MPI_Finalize();
      //newstructure
// return;
//newstructure
return 0;
    }

    //    printf("Reading GPS...\n");
    //    strcpy(GpsFilename,argv[2]); // CBV 05/08/16


    err = read_stk(GpsFilename, GPS);

    if (err != 0 ) {
      MPI_Finalize();
      //newstructure
// return;
//newstructure
return 0;
    }


    rCygnss = calc_mean_radius(CYGNSS,nCygnss_prop);
    rGps    = calc_mean_radius(GPS,   nGps_prop);

    err = fill_specular_point_grid(rCygnss, rGps);
    
    if (err != 0 ) {
      MPI_Finalize();
      //newstructure
// return;
//newstructure
return 0;
    }







    // CBV 05/08/16
   //    strcpy(filename,argv[3]);
    //newstructure
/*     strcpy(filename, ""); */
/*     strcpy(filename, "./output/"); */
    /* strcpy(filename, "./spock/"); */
    /* //newstructure */
    /* strcat(filename, text); */
    /* strcat(filename, "/"); */

    strcpy(filename,path_output_run_name_temp);

    strcat(filename, dir_output_run_name_temp);

    strcat(filename, "/");
    // END CBV 05/08/16


    err =  find_specular_points( GPS, CYGNSS, filename, lonshift, rotate, yawoverall, iUseMin, LimitToOneDay, nGps_prop, nCygnss_prop, ATTITUDE); // CBV 05/08/16

 if (iProc == 0) printf("Done computing the specular point locations.\n");

    if (err != 0 ) {
      MPI_Finalize();
      //newstructure
// return;
//newstructure
return 0;
    }

  }


  free(GPS);
  free(CYGNSS);
  MPI_Finalize();

}



