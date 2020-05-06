/////////////////////
// HEADERS //
////////////////////
#include<stdio.h>
#include<stdlib.h>  
#include<string.h>  
#include<malloc.h>
#include<math.h>
#include<time.h>

//CABECERAS DE GSL                                                                  
#include<gsl/gsl_errno.h>
#include<gsl/gsl_sort.h>
#include <gsl/gsl_permutation.h>
#include<gsl/gsl_sort_uint.h>
#include<gsl/gsl_fit.h>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>

//===============
// ADITIONAL CODE 
//===============

//=========
// UNITS
//=========

#define G 43007.1

//=====
// MAIN
//=====

int main(int argc, char *argv[])
{

  if(argc != 3)

    {
      printf("incorrect number of arguments, use as:\n");
      printf("%s \t data_file \t nlines\n",argv[0]);
    }

  int i,nlines;

  //command line arguments

  char *infile;
  
  FILE *pf=NULL, *wfile=NULL;

  double *tot_mass, *header_time, *Sfr;
  double deltat;
  
  infile = argv[1];

  nlines = atoi(argv[2]);
  
  pf = fopen(infile,"r");
  wfile =fopen("tasa_formacion_estelar.dat","w");
  
  
  tot_mass = (double *) malloc((size_t) nlines*sizeof(double));
  header_time = (double *) malloc((size_t) nlines*sizeof(double));
  Sfr = (double *) malloc((size_t) nlines*sizeof(double));
  
  //=========================
  // read integrated_sfr file
  //=========================

  for(i=0;i<nlines;i++)
    
    {
      fscanf(pf,"%lf %lf\n",&tot_mass[i],&header_time[i]);
    }


  deltat = 0.002;

  Sfr[0] = tot_mass[0]/deltat;
  fprintf(wfile,"%16.8e %16.8e %16.8e\n",header_time[0], Sfr[0], tot_mass[0]);
  printf("%16.8e %16.8e %16.8e\n",header_time[0], Sfr[0], tot_mass[0]);
  
  for(i=1;i<nlines;i++)
    
    {
      Sfr[i] = (tot_mass[i] - tot_mass[i-1])/deltat;
      fprintf(wfile,"%16.8e %16.8e %16.8e\n",header_time[i],Sfr[i],tot_mass[i]);
      printf("%16.8e %16.8e %16.8e\n",header_time[i], Sfr[i], tot_mass[i]);
    }
  
  return 0;
}
