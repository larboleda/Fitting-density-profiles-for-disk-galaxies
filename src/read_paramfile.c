#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "globalvars.c"

int readparamfile(char *paramfile)

{

  FILE *pfParam = NULL;

  pfParam = fopen(paramfile,"r");

  //returnRead = fscanf(pfParam,"%s",infile);
  returnRead = fscanf(pfParam,"%d",&component);
  returnRead = fscanf(pfParam,"%d",&symmetry);
  returnRead = fscanf(pfParam,"%lf",&factor);
  returnRead = fscanf(pfParam,"%lf",&Rmax);
  returnRead = fscanf(pfParam,"%lf",&zmax);
  returnRead = fscanf(pfParam,"%d",&minParticles);
  returnRead = fscanf(pfParam,"%d",&particlesPerBin);
  returnRead = fscanf(pfParam,"%d",&binType);
  returnRead = fscanf(pfParam,"%d",&fitOn);
  returnRead = fscanf(pfParam,"%d",&flag);
  returnRead = fscanf(pfParam,"%d %d %d %d %d %d",&Nhost[0],&Nhost[1],&Nhost[2],&Nhost[3],&Nhost[4],&Nhost[5]);
  returnRead = fscanf(pfParam,"%d %d %d %d %d %d",&Nsat[0],&Nsat[1],&Nsat[2],&Nsat[3],&Nsat[4],&Nsat[5]);
  returnRead = fscanf(pfParam,"%lf %lf %lf",&factor_comp[0], &factor_comp[1], &factor_comp[2]);
  returnRead = fscanf(pfParam,"%lf %lf %lf",&r_scalength[0],&r_scalength[1],&r_scalength[2]);
  returnRead = fscanf(pfParam,"%lf",&Deltat);
  


  fclose(pfParam);

  return 0;
}
