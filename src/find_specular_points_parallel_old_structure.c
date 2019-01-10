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

int iDebug;

int nProcs;
int iProc;

int nSpecularPointsMax = 4;

double re    = 6378000.0;
double rpole = 6378000.0;
double req   = 6378000.0;
double dtor  =       0.017453292;
double pi    =       3.1415927;
double twopi =       6.2831855;

void lam_darwin_malloc_linker_hack(){};

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

enum {nTimesMax=51840000};
enum {nGps=31};
enum {nCygnss=8};

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
  char name[200] ;; // CBV 05/08/16
  int time_ymdhmsm[nTimesMax][7];  // CBV 05/09/16

};

enum {nThetaAntenna = 36};
enum {nPhiAntenna   = 73};
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
  float gain[nPhiAntenna][nThetaAntenna];
  float dTheta, dPhi;

  float Sigma0[nWinds][nAngles];
  float dIncidence;

  char AntennaFileName[150];
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
    return;
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

  fpAntenna = fopen(AntennaInfo.AntennaFileName,"r");

  if (!fpAntenna) {
    printf("Can not read antenna file : %s\n",AntennaInfo.AntennaFileName);
    iError = 1;
    return;
  }

  fread(&i,sizeof(int), 1, fpAntenna);
  fread(&j,sizeof(int), 1, fpAntenna);
  //printf("%d %d\n",i,j);

  for (i=0; i<nPhiAntenna; i++) {
    for (j=0; j<nThetaAntenna; j++) {
      fread(&tmp,sizeof(float),1,fpAntenna);
      AntennaInfo.phi[i] = tmp;
      fread(&tmp,sizeof(float),1,fpAntenna);
      AntennaInfo.theta[j] = tmp;
      fread(&tmp,sizeof(float),1,fpAntenna);
      AntennaInfo.gain[i][j] = tmp;
      //fread(AntennaInfo.phi[i][j],sizeof(float),1,fpAntenna);
      //fread(AntennaInfo.theta[i][j],sizeof(float),1,fpAntenna);
      //fread(AntennaInfo.gain[i][j],sizeof(float),1,fpAntenna);

      //if (AntennaInfo.gain[i][j] > 20.0) printf("%d %d %f %f %f\n",i,j,AntennaInfo.phi[i], AntennaInfo.theta[j], AntennaInfo.gain[i][j]);

    }
  }
  close(fpAntenna);
  AntennaInfo.dTheta = 180.0/nThetaAntenna;
  AntennaInfo.dPhi = AntennaInfo.dTheta;

  AntennaInfo.factor = AntennaInfo.ptx * AntennaInfo.gtx *
    AntennaInfo.a0 * (AntennaInfo.lambda * AntennaInfo.lambda) /
    (4*4*pi*pi);

  fpSigma0 = fopen(AntennaInfo.Sigma0FileName,"r");
  if (!fpSigma0) {
    printf("Can not read Sigma0 file : %s\n",AntennaInfo.Sigma0FileName);
    iError = 1;
    return;
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
    return;
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

  FILE *fpSTK;
  char line[600], line2[600], testline[40], testline2[40];
  int IsFound, IsDone, IsDoneInner, IsPropagatorFile;
  int iSat, i, res, err;
  int iTime[7], nSats, nPts;

  double x, y, z, r, xy;

  err = 0;

  fpSTK = fopen(filename,"r");
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

	//printf("here first!\n");

	//printf("nPts : %d \n",nPts);

	if (fgets(line,500,fpSTK) == NULL) IsDoneInner = 1;

	if (strlen(line) > 20) {

	  err = strcompress(line);
	  //printf("line-->%s<--\n",line);
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
	      float AntennaTilt, float RotationAngle, float yawoverall,
	      struct specular_point_info *info) {

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
    return;
  }
  yi = yi/mag;
  yj = yj/mag;
  yk = yk/mag;

  C2Sx = CygX - SpX;
  C2Sy = CygY - SpY;
  C2Sz = CygZ - SpZ;
  C2Smag = sqrt(C2Sx*C2Sx + C2Sy*C2Sy + C2Sz*C2Sz);
  C2Sx = C2Sx/C2Smag;
  C2Sy = C2Sy/C2Smag;
  C2Sz = C2Sz/C2Smag;

  info->CygnssToSpecularPoint = C2Smag;

  if (iDebug > 0) {
    printf("c: %f %f %f\n",CygX/1000.0, CygY/1000.0, CygZ/1000.0);
    printf("g: %f %f %f\n",GpsX/1000.0, GpsY/1000.0, GpsZ/1000.0);
    printf("s: %f %f %f\n",SpX/1000.0, SpY/1000.0, SpZ/1000.0);
    printf("c2s: %f %f %f\n",C2Sx, C2Sy, C2Sz);
  }

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


    tmpX = SPscX * cy - SPscY * sy;
    tmpY = SPscY * cy + SPscX * sy;
    tmpZ = SPscZ;

    t = AntennaInfo.rolls[iAnt]*dtor;
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

    theta_sat  = acos(stosZr )/dtor;
    //theta_satm = acos(stosZrm)/dtor;

    xy = sqrt(stosXr*stosXr + stosYr*stosYr) + 1.0e-7;
    phi_sat = acos(stosXr/xy)/dtor + 90.0;
    if (phi_sat > 360.0) phi_sat -= 360.0;

    //xy = sqrt(stosXrm*stosXrm + stosYrm*stosYrm) + 1.0e-7;
    //phi_satm = acos(stosXrm/xy)/dtor; //+ 90.0;
    //if (stosYrm < 0.0) phi_satm = 360.0-phi_satm;
    //phi_satm += 90.0;
    //if (phi_satm > 360.0) phi_satm -= 360.0;
  
    if (iDebug > 0) printf("find_gain: before gain lookup\n");

    ii = (int) (phi_sat / AntennaInfo.dPhi);
    jj = (int) (theta_sat / AntennaInfo.dTheta);

    fx = (  phi_sat - AntennaInfo.phi[ii])  /AntennaInfo.dPhi;
    fy = (theta_sat - AntennaInfo.theta[jj])/AntennaInfo.dTheta;

    gainR = (1-fx)*(1-fy) * AntennaInfo.gain[ii  ][jj  ] +
            (  fx)*(1-fy) * AntennaInfo.gain[ii+1][jj  ] +
            (1-fx)*(  fy) * AntennaInfo.gain[ii  ][jj+1] +
            (  fx)*(  fy) * AntennaInfo.gain[ii+1][jj+1];

    if (gainR > GainHighest) GainHighest = gainR;

  }


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

  info->NormPower =
    info->gain *
    AntennaInfo.factor /
    (C2Smag*C2Smag) / (S2Gmag*S2Gmag);

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
			  int nGps_prop, int nCygnss_prop) { // CBV 05/08/16

  // On my computer, I can't seem to allocate enough memory to actually
  // store everything.  So, I will just write it out to a file.

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

        if (iProc == 0) printf("iCygnss : %d %d\n",iCygnss, CYGNSS[iCygnss].IsGood);

    if (CYGNSS[iCygnss].IsGood == 1) {

      cDt = CYGNSS[iCygnss].time[1] - CYGNSS[iCygnss].time[0];
      nPtsInner = (int) cDt;

      //      nPts = (CYGNSS[iCygnss].nPts/1400) * 1440;
      //      if (nPts > 1440 && LimitToOneDay) nPts = 1440;
      nPts = CYGNSS[iCygnss].nPts-1;

      //printf("nPts : %d\n",nPts);

      for (iPt=0; iPt < nPts; iPt++) { //      for (iPt=1; iPt < nPts; iPt++) { // CBV 05/09/16

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

		err = find_gain(specularpoint.CygX,specularpoint.CygY,specularpoint.CygZ,
				cVx, cVy, cVz,
				specularpoint.GpsX,specularpoint.GpsY,specularpoint.GpsZ,
				gVx, gVy, gVz,
				specularpoint.SpX,specularpoint.SpY,specularpoint.SpZ,
				45.0, rotate, yawoverall,
				&specularpoint);

		// Save all of the specular point and GPS information:

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

	      }

	    }

	  }


	  // Store the highest normalized gain GSP satellite:
	  if (iMaxNormPower > -1) {
	    iOrder[iMaxNormPower] = iNumber;
	    NormPowerAllGps[iMaxNormPower] = -1.0e32;
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
	      NormPowerAllGps[iMaxNormPower] = -1.0e32;
	      iNumber++;
	    }
	  }

	  for (iGps=0; iGps<nGps_prop; iGps++) {

	    if (iOrder[iGps] < nSpecularPointsMax) {

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
	      tmp =  (float) GPS[iGps].x[iPt]/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  (float) GPS[iGps].y[iPt]/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  (float) GPS[iGps].z[iPt]/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput);

	      tmp =  (float) GPS[iGps].vx[iPt]; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  (float) GPS[iGps].vy[iPt]; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  (float) GPS[iGps].vz[iPt]; fwrite(&tmp,sizeof(tmp),1,fpOutput);

	      tmp =  (float) CYGNSS[iCygnss].x[iPt]/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  (float) CYGNSS[iCygnss].y[iPt]/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  (float) CYGNSS[iCygnss].z[iPt]/1000; fwrite(&tmp,sizeof(tmp),1,fpOutput);

	      tmp =  (float) CYGNSS[iCygnss].vx[iPt]; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  (float) CYGNSS[iCygnss].vy[iPt]; fwrite(&tmp,sizeof(tmp),1,fpOutput);
	      tmp =  (float) CYGNSS[iCygnss].vz[iPt]; fwrite(&tmp,sizeof(tmp),1,fpOutput);

	      //	      }
	      /* /\* if ((iPtInner == 0) && (iCygnss == 0) && (iPt %2 == 0)){ *\/ */
	      /* /\* 	printf("%f %f %f\n", GPS[iGps].x[iPt]/1000, GPS[iGps].y[iPt]/1000, GPS[iGps].z[iPt]/1000); *\/ */
	      /* } */


	      /* tmp = (float) SpecularPointsAllGps[iGps].power; */
	      /* fwrite(&tmp,sizeof(tmp),1,fpOutput); */
	      /* tmp = (float) SpecularPointsAllGps[iGps].NormPower; */
	      /* fwrite(&tmp,sizeof(tmp),1,fpOutput); */
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

//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

int main(int argc, char *argv[]) {

   char sat_nb_str[10];
  int err, i, lonshift, iUseMin, LimitToOneDay = 0;
  float rotate, yawoverall;
  char filename[150], sProc[3];
  char CygnssFilename[256], GpsFilename[256], YawFileName[150];
  float dt;
  double rCygnss, rGps;

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



 // CBV 05/08/16
 char filename_input_propagator[150];
 char filename_without_extension[150];

 strcpy(filename_input_propagator, "./input/main_input/");
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
   strcpy(gps_tle_filename, "./input/tle/constellation_gps_tle/");
   strcat(gps_tle_filename,text);
   gps_tle_file = fopen(gps_tle_filename, "r");
   while ( (read = getline(&line, &len, gps_tle_file)) != -1 ){
     if (line[0] == '1'){
       nGps_prop = nGps_prop + 1 ;
     }
   }
   fclose(gps_tle_file);
 }
 else{
   nGps_prop = 0;
 }

 // NAME OF OUTPUT FILES
 int sss; char text2[256];int find_file_name;char *next; int find_file_name2;
 found_eoh = 0;
 rewind(file_input_propagator);
 while ( (found_eoh == 0 ) && (!feof(file_input_propagator)) ) {
   getline(&line, &len, file_input_propagator);
   sscanf(line, "%s", text);
   if (  strcmp( "#OUTPUT", text  ) == 0 )  {
     found_eoh = 1;
   }
 }
 if (feof(file_input_propagator)){
   printf("***! No section #OUTPUT was found in %s. The program will stop.***!\n", filename_input_propagator); MPI_Finalize(); exit(0);
 }
 getline(&line,&len,file_input_propagator);
 sscanf(line, "%s", text);
     if (strcmp(text, "out") == 0){
      strcpy(filename_without_extension, "");
      strncat(filename_without_extension, &filename_input_no_path[0], (int)(strchr(&filename_input_no_path[0], '.') - &filename_input_no_path[0]));
       strcpy(text, filename_without_extension);
     }
     
  for (sss = 0; sss<nCygnss_prop; sss++){
    if (sss == 0){
      // CYGNSS CONSTELLATION FILE
      strcpy(CygnssFilename, "");
     strcpy(CygnssFilename, "./output/");
     strcat(CygnssFilename, text);

     strcat(CygnssFilename, "/CONSTELLATION_CYGNSS_for_run_");
     strcat(CygnssFilename, text);
     strcat(CygnssFilename, "1.txt");
     // GPS CONSTELLATION FILE
      strcpy(GpsFilename, "");
     strcpy(GpsFilename, "./output/");
     strcat(GpsFilename, text);
     strcat(GpsFilename, "/CONSTELLATION_GPS_for_run_");
     strcat(GpsFilename, text);
     strcat(GpsFilename, "1.txt");
     
     // Names for the CYGNSS and GPS satellites
     strcpy(CYGNSS[sss].name, "");
     strcpy(CYGNSS[sss].name, text);
      sprintf(sat_nb_str, "%d", sss+1);
      strcat(CYGNSS[sss].name, sat_nb_str);
    }
    else{
      strcpy(CYGNSS[sss].name, "");
      strcpy(CYGNSS[sss].name, text);
      sprintf(sat_nb_str, "%d", sss+1);
      strcat(CYGNSS[sss].name, sat_nb_str);
    }
    //    printf("CYGNSS[%d].name = <%s>\n",sss, CYGNSS[sss].name);
  }

 fclose(file_input_propagator);

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
	return;
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

    strcpy(AntennaInfo.AntennaFileName,"./input/specular_points/antenna.bin");
    strcpy(AntennaInfo.AntennaInfoFileName,"./input/specular_points/antenna.info");
    strcpy(AntennaInfo.Sigma0FileName,"./input/specular_points/sigma0_table.bin");
    err = read_antenna_info();

    //        printf("Reading CYGNSS...\n");
    //    strcpy(CygnssFilename,argv[1]); // CBV 05/08/16


    err = read_stk(CygnssFilename, CYGNSS);

    if (err != 0 ) {
      MPI_Finalize();
      return;
    }

    //    printf("Reading GPS...\n");
    //    strcpy(GpsFilename,argv[2]); // CBV 05/08/16


    err = read_stk(GpsFilename, GPS);

    if (err != 0 ) {
      MPI_Finalize();
      return;
    }

    rCygnss = calc_mean_radius(CYGNSS,nCygnss_prop);
    rGps    = calc_mean_radius(GPS,   nGps_prop);

    err = fill_specular_point_grid(rCygnss, rGps);
    
    if (err != 0 ) {
      MPI_Finalize();
      return;
    }

    strcpy(sProc,c_int_str(iProc,2));

    // CBV 05/08/16
   //    strcpy(filename,argv[3]);
    strcpy(filename, "");
    strcpy(filename, "./output/");
    strcat(filename, text);
    strcat(filename, "/");
    // END CBV 05/08/16

    err =  find_specular_points( GPS, CYGNSS, filename, lonshift, rotate, yawoverall, iUseMin, LimitToOneDay, nGps_prop, nCygnss_prop); // CBV 05/08/16
    if (err != 0 ) {
      MPI_Finalize();
      return;
    }

  }

  free(GPS);
  free(CYGNSS);
  MPI_Finalize();

}
