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
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>


////////////////////
// MACROS //
///////////////////

#define X 0
#define Y 1
#define Z 2


//=========
// UNITS
//=========

#define G 43007.1

/*constants CGS*/
#define mp  1.6726e-24
#define KB  1.3806e-16
#define gama  (5.0/3.0)
#define GAMMA_MINUS1  (gama-1)
#define UnitLength_in_cm 3.085678e21 //kpc/h
#define UnitVelocity_in_cm_per_s 1.0e5 //1km /s
#define UnitTime_in_s  (UnitLength_in_cm / UnitVelocity_in_cm_per_s)
#define UnitMass_in_g 1.989e43 
#define UnitEnergy_in_cgs  (UnitMass_in_g * pow(UnitLength_in_cm, 2) / pow(UnitTime_in_s, 2))
//#define unitFactor  1.0e10
#define unitFactor  (UnitEnergy_in_cgs/UnitMass_in_g)
#define xH  0.76

//#define mu 1.1967658843732112192e-24 //(4mp/(1+3xH+(4xHXe)))
#define mu  ((4.0*mp)/(1.0 + 3.0 * xH))

////////////////////////////////////
// DATA STRUCTURES //
///////////////////////////////////

//******************
// For Gadget data
typedef struct 
{
  int type;
  int ID_generation;

#ifdef LONGIDS
  unsigned long long id;
#else
  unsigned int id;
  unsigned int ID_child_number;
#endif
  float pos[3];
  float vel[3];
  float mass;
#ifdef STELLARAGE
  float stellar_age;
#endif
#ifdef METALS
  float metallicity;
#endif
#ifdef OUTPUTPOTENTIAL
  float pot;
#endif
#ifdef OUTPUTACCELERATION
  float acce[3];
#endif
#ifdef OUTPUTTIMESTEP
  float timestep;
#endif
}particulas;

typedef struct
{
  float U;
  float rho;
#ifdef COOLING 
  float Ne;
  float Nh;
#endif
  float h;
#ifdef SFR
  float sfr;
#endif
#ifdef OUTPUTCHANGEOFENTROPY
  float ecr;
#endif
#ifdef OUTPUTCOOLRATE
  float coolr;
#endif
#ifdef MAGNETIC_FIELD
  //float B[3];
  double B[3];
  float divB;
  float Phi;
  float GradPhi[3];
#endif
#ifdef OUTPUT_DIV_CURL
  float div;
  float curl;
#endif
#ifdef OUTPUT_VORTICITY
  float vorticity[3];
#endif
#ifdef WINDS
  float delaytime;
#endif
} gas_properties;

typedef struct
{
  int Npart[6];
  double mass[6];
  double time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned int npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  int flag_stellarage;
  int flag_metals;
  unsigned int npartTotalHighWord[6];
  int  flag_entropy_instead_u;
  char     fill[256 - 6*4 - 6*8 - 2*8 - 2*4 - 6*4 - 2*4 - 4*8 - 2*4 - 6*4 - 1*4];  /* fills to 256 Bytes */
} io_header;

//*********************

typedef struct
{
  int id;
  double cm[3];
  double vcm[3];
  double lcm[3];
  double M;
} CM;

typedef struct
{
  double xini;
  double xfin;
  double x;
  size_t indexIni;
  size_t indexFin;
}Bin;

typedef struct
{
  unsigned int id;
  float pos[3];
  float vel[3];
  float mass;
  int type;
  float metallicity;
  float pot;
  
} galaxies;

typedef struct
{
  float U;
  float rho;
  float B[3];
  float divB;
  
#ifdef COOLING 
  float Ne;
  float Nh;
#endif
#ifdef SFR
  float sfr;
#endif
#ifdef OUTPUTCOOLRATE
  float coolr;
#endif
#ifdef WINDS
  float delaytime;
#endif
} galaxy_gas_properties;

//=================
// GLOBAL VARIABLES
//=================

particulas *particles;
galaxies *galaxy;
gas_properties *gaspro;
galaxy_gas_properties *galaxy_gaspro;

io_header Header;
CM cmdummy;
Bin *bins;

int N_part_total,N_part[6],N_min,N_max;
int n0,n1,n2,n3,n4,n5;
int returnRead;
char arr_name[4];

double Mass_total,Mass_tot[6];

float *U;

char binning[500];
int indexmin, indexmax, nDisk;

//for individual comp of meerger
int N_galaxy_tot;

//===================
//epicyclic variables
//===================

//double theta;
double *potential, *temperatures,*vsound;

double *cylindrical_radii, *Sfr;
double *v_disp_R, *v_radial,*v_disp_Z;
//double *v_disp_theta;
double *surf_dens, *dens_stars;
double *pot_soft;
double *r_soft;
int npot; //numer of points in new soft potential func.
int nfinalBins;

//================
// DATA STRUCTURES 
//================

typedef struct
{
  double R;
  double Sigma;
  double variance;
  int accepted;
}Profile;

//===================
//paramfile variables
//===================

int component, symmetry; //component in gadget
double factor; // factor to select bins number 
double Rmax, zmax; 
int minParticles, particlesPerBin, binType, fitOn, flag;  //  minimun particles number per bin, Type of binning,  particles number per bin
int Nhost[6], Nsat[6], flag;
double factor_comp[3];
double Deltat;
double r_scalength[3]; //entrada cero es rh_host, segunda entrada rh_sat


//=========
// ROUTINES
//=========

#include<epicyclic.c>
#include<medley.c>
#include<Bsplines.c>
#include<gas_temperature.c>
#include<surface_sfr_analysis.c>
#include<write_data.c>
#include<mass_profile.c>
#include<computing_index.c>
#include<integrated_sfr_time.c>
#include<translatios_rotations_split_merger.c>

#include<read_paramfile.c>

#include<split_merger.c>

#ifdef INPUT_OUTPUT_GADGET
#include<input_output_data.c>
#endif

#ifdef TRANSLATIONS_ROTATION
#include<translatios_rotations.c>
#endif

#ifdef BINS
#include<makes_bins.c>
#endif

/*
#ifdef MEDLEY
#include<../useful_routines/medley.c>
#endif
*/

#ifdef ANALYSIS
#include<analysis.c>
#endif


/*
/////////////////////////////////////////
//MACRO for gadget 2 format
/////////////////////////////////////////

#ifdef GADGET2_FORMAT
#define BLOCK_NAME_G2 \
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    \
  returnRead = fread(&arr_name,sizeof(char),4,fdata);	    \
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    \
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    \
  printf("------we've read \t %s\n",arr_name);)
#endif
*/
