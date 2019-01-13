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

// STILL NEED TO CHECK THAT FUNCTION IS CORRECT

 /////////////////////////////////////////////////////////////////////////////////////////
//
//  Name:           load_gravity
//  Purpose:        Loads EGM96 gravity coefficients
//  Assumptions:    This isn't agnositic to file format, don't change current file
//  References      BMW
//
//  Change Log:
//      |   Developer   |       Date    |   SCR     |   Notes
//      | --------------|---------------|-----------|-------------------------------
//      | J. Getchius   | 05/17/2015    |   ---     | Initial Implementation
//      | C. Bussy-Virat| 10/14/2015    |   ---     | Change the inputs
//
/////////////////////////////////////////////////////////////////////////////////////////
//newstructure
//int load_gravity( GRAVITY_T *GRAVITY, char main_directory_location[256])
int load_gravity( GRAVITY_T *GRAVITY,  char path_to_spice[256])
//newstructure
{
  FILE *fp = NULL;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  int l = 0;
  int m = 0;
  float C;
  float S;
  float Csig;
  float Ssig;
  char text_location[256];
      double fac, facnum, facden, kfac;

  //  strcpy(text_location, main_directory_location);
  //newstructure
  //  egm96_to360_not_norm.txt is put with the SPICE files path_to_spice. ASSUMPTION: the path of spice specified in the main inoput file in the sectioN #SPICE must be the same as the path of SPICE installation in the Makefile, with /data at the end. So if in the Makefile the SPICE directory is hello/hi/ then the path of spice in the main input file in section #SPICE must be hello/hi/data/
  strcpy(text_location,path_to_spice);
  strcat(text_location, "egm96_to360.txt");//_not_norm.txt");
    //  strcpy(text_location, "input/egm96_to360_not_norm.txt");
  //newstructure
  // !!!!!!!!!! to read normalized coefficients then put ./code/egm96_to360.txt
  fp = fopen(text_location, "r");
  if (fp) {
    while ( (read = getline(&line, &len, fp)) != -1 ) {

      sscanf(line, "%i %i %e %e %e %e" , &l, &m, &C, &S, &Csig, &Ssig);
      GRAVITY->Clm[l][m] = (double)C; // l is the degree (CBV)
      GRAVITY->Slm[l][m] = (double)S; // m is the order (CBV)


      if (m == 0){
	kfac = 1;
      }
      else if (m > 0){
	kfac = 2;
      }
      else
	{
	  print_test();
	  exitf();
	}

      facnum = kfac * (2*l+1)*(double)factorial(l-m);
      facden = (double)factorial(l+m);
      //printf("%d %d %e %e\n", l,m,facnum, facden);
      fac = sqrt(facnum/facden);
      GRAVITY->Clm[l][m]  = GRAVITY->Clm[l][m] * fac;
      GRAVITY->Slm[l][m] = GRAVITY->Slm[l][m] * fac;
    }
    
    free(line);
    
    fclose(fp);

  } else {

    printf("The file input/egm96_to360_not_norm.txt has not been opened. The program will stop.\n");
    MPI_Finalize();
    exit(0);

  }

  return 0;
    
}
