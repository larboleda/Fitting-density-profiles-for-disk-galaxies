//===============
// ADITIONAL CODE 
//===============

#include<useful_routines.c>

//=========
// ROUTINES
//=========

int fitExponentialProfile(int nBins, double x[nBins], double y[nBins], double variance[nBins], FILE *pfLog, FILE *pfFit);
int binningType(int binType, int nBins, double minR, double maxR, int nDisk, double R[nDisk], int deltaN);

//=====
// MAIN
//=====

int main(int argc, char *argv[])
{

  if(argc != 3)

    {
      printf("incorrect number of arguments, use as:\n");
      printf("%s \t snap \t parameter file \n",argv[0]);
    }

  //int i; 
  //int *ind;
 
  CM cmDisk, cmdummy;

  //command line arguments

  char *infile=NULL, *paramfile=NULL;
  
  infile = argv[1];
  paramfile = argv[2];
  
  //===============
  // read paramfile
  //===============
  
  readparamfile(paramfile);
  
  //====================
  // read/write snapshot
  //====================

  read_gadget(infile);

  
#ifdef PRINT_HEADER
  print_header(infile);
#endif
  
#ifdef PRINT_ASCII_ONE_TYPE
  print_data_ascii(infile,component);
#endif
  
#ifdef PRINT_ASCII_ALL_TYPES
  print_data_ascii_all_types(infile);
#endif

//=============================================
//to split the merger into individual galaxies:
//host or sat according to flag in param file
//=============================================
  
#ifdef MERGER
  
  split_merger();

  //===================
  // WRITTE OPTIONS
  //===================
  
#ifdef PRINT_DATA
  
  print_data_toascii();
  
#endif

  
#endif

  //=============================================
  //to orientate the galaxy we will use the disk 
  //=============================================

  indices(2);
  
  //=====================================================
  //Computing center of mass and angular momentum of disk
  //=====================================================
  
  printf("\n 1. Computing center of mass and angular momentum of disk\n");
  if(flag==0)
    {
      centerMass(NULL, nDisk, indexmin, indexmax, &cmDisk);
      compute_angmom(NULL, nDisk, indexmin, indexmax, &cmDisk);
    }
  else
    {
      centerMass_mergercomp(NULL, nDisk, indexmin, indexmax, &cmDisk);
      compute_angmom_mergercomp(NULL, nDisk, indexmin, indexmax, &cmDisk);
      //printf("\n (2.****BEFORE TRANSLATE: cmDisk) cmx = %16.8f \t cmy = %16.8f \t cmz = %16.8f\n",cmDisk.cm[0], cmDisk.cm[1], cmDisk.cm[2]);
    }
  
  //===============================================
  //Traslating all galaxy to center of mass of disk
  //===============================================

  printf("\n 2. Traslating galaxy to center of mass of disk\n");
  
  if(flag==0)
    {
      totalTranslationMinus(&cmDisk,N_part_total);
      centerMass(NULL, nDisk, indexmin, indexmax, &cmdummy);
      compute_angmom(NULL, nDisk, indexmin, indexmax, &cmdummy);
    }
  else
    {
      traslate_mergercomp(&cmDisk, N_galaxy_tot );
      centerMass_mergercomp(NULL, nDisk, indexmin, indexmax, &cmdummy);
      //printf("\n (3.****AFTER TRANSLATE: DUMMY) cmx = %16.8f \t cmy = %16.8f \t cmz = %16.8f\n",cmdummy.cm[0], cmdummy.cm[1], cmdummy.cm[2]);
      compute_angmom_mergercomp(NULL, nDisk, indexmin, indexmax, &cmdummy);
    }

  //============================================================
  //Orientating the galaxy with the angular momentum of the disk
  //============================================================

  /*ind = (int *)malloc((size_t)nDisk*sizeof(int));
  if(ind == NULL)
    {
      printf("Allocation of ind failed\n");
      exit(0);
    }
  
  for( i=indexmin; i<indexmax; i++)
    {
      //    ind[i-indexmin] =  particles[i].id;
      ind[i-indexmin] =  i;
      }*/
  
  float posL[3];

  posL[0] = cmdummy.lcm[0];
  posL[1] = cmdummy.lcm[1];
  posL[2] = cmdummy.lcm[2];
  
  if(flag==0)
    {
      printf("\n 3. Orientating/rotating galaxy with disk's angular momentum\n");
      printf("using %f %f %f\n", cmdummy.lcm[0], cmdummy.lcm[1], cmdummy.lcm[2]);
      rotate(posL, N_part_total);

      printf("\n 4. Computing galaxy's new center of mass and angular momentum\n");
      centerMass(NULL, nDisk, indexmin, indexmax, &cmdummy);
      compute_angmom(NULL, nDisk, indexmin, indexmax, &cmdummy);

      
      printf("\n===========================\n");
      printf("Isolated galaxy simulation:\n");
      printf("===========================\n");
      printf("NEW RCM : %e %e %e\n",cmdummy.cm[0],cmdummy.cm[1],cmdummy.cm[2]);
      printf("NEW VCM : %e %e %e\n",cmdummy.vcm[0],cmdummy.vcm[1],cmdummy.vcm[2]);
      printf("NEW LCM : %e %e %e\n",cmdummy.lcm[0],cmdummy.lcm[1],cmdummy.lcm[2]); 
      printf("Done.\n");
      printf("\n");
      
      //free(ind);
    }
  else
    {
      printf("\n 3. Orientating galaxy with disk's angular momentum\n");
      printf("using %f %f %f\n", cmdummy.lcm[0], cmdummy.lcm[1], cmdummy.lcm[2]);
      rotate_mergercomp(posL, N_galaxy_tot);

      printf("\n 4. Computing galaxy's new center of mass and angular momentum\n");
      //centerMass_mergercomp(ind, nDisk, 0, 0, &cmdummy);
      centerMass_mergercomp(NULL, nDisk, indexmin, indexmax, &cmdummy);
      //printf("\n (4.****NEW CM OF DISK: DUMMY) cmx = %16.8f \t cmy = %16.8f \t cmz = %16.8f\n",cmdummy.cm[0], cmdummy.cm[1], cmdummy.cm[2]);
      compute_angmom_mergercomp(NULL, nDisk, indexmin, indexmax, &cmdummy);
      
      printf("\n======================================\n");
      printf("Individual galaxy component of merger:\n");
      printf("======================================\n");
      printf("NEW RCM : %e %e %e\n",cmdummy.cm[0],cmdummy.cm[1],cmdummy.cm[2]);
      printf("NEW VCM : %e %e %e\n",cmdummy.vcm[0],cmdummy.vcm[1],cmdummy.vcm[2]);
      printf("NEW LCM : %e %e %e\n",cmdummy.lcm[0],cmdummy.lcm[1],cmdummy.lcm[2]); 
      printf("Done.\n");
      printf("\n");
      
      //free(ind);
    }
  
  //==================================
  //estimating surface density profile
  //==================================
  
#ifdef SURF_DENS_PROFILE
  
  //============================================
  //updating indexmax & indexmin to analyse just
  //the component from paramfile
  //============================================
  
  nDisk = indexmin = indexmax = 0;
  
  indices(component);
  
  if(symmetry == 0)  //cylindrical
    {
      printf("\n --> Computing surface density profile for component %d\n",component);
      surface_density_profile(infile,indexmin,indexmax,nDisk,binning);
    }
  
  if(symmetry == 1)//spherical
    {
      printf("\n --> Computing volumetric density profile for component %d\n",component);
      volumetric_density(indexmin, indexmax, nDisk, binning);
    }
#endif
  

 
  //===================================================
  //estimating epyciclic frecuency & Toomre's parameter
  //===================================================
  
#ifdef TOOMRE

  //============================================
  //updating indexmax & indexmin to analyse just
  //the component from paramfile
  //============================================
  
  nDisk = indexmin = indexmax = 0;
  
  indices(component);

  //===========================
  //smothing potential function
  //===========================
  
#ifdef SOFT_POTENTIAL
  
  printf("\n --> Softening potential function with Bspline for component %d\n",component);
  
  Bspline(nfinalBins,cylindrical_radii,potential);
  
#endif

  if(component == 0) //gas
    {
      
      printf("\n --> Computing gas temperature\n");
      
      //==========================
      //estimating gas temperature
      //==========================
      
      temperature(indexmin,indexmax,nDisk,binning);
      
      printf("\n --> Computing epyciclic frecuency and Toomre's parameter for component %d\n",component);
      
#ifdef SOFT_POTENTIAL
      epicyclic_gsl_generic_gas(npot,r_soft,pot_soft,surf_dens,v_disp_R, v_disp_Z, temperatures,vsound, gsl_interp_cspline);
#else
      epicyclic_gsl_generic_gas(nfinalBins,cylindrical_radii,potential,surf_dens,v_disp_R, v_disp_Z, temperatures,vsound, gsl_interp_cspline);
#endif
      
      free(temperatures);
    }
  
  if(component == 2) //stars
    {
      printf("\n --> Computing epyciclic frecuency and Toomre's parameter for component %d\n",component);
      
#ifdef SOFT_POTENTIAL
      epicyclic_gsl_generic(npot,r_soft,pot_soft,surf_dens,v_disp_R, v_disp_Z,gsl_interp_cspline);
#else
      epicyclic_gsl_generic(nfinalBins,cylindrical_radii,potential,surf_dens,v_disp_R, v_disp_Z,gsl_interp_cspline);
#endif
    }
  
#endif
 
  //=======================
  //Star Formation Analysis
  //=======================

#ifdef SF_ANALYSIS
  sfr(infile, binning);
#endif
  
  //=======================
  //Mass profile 
  //=======================
  
#ifdef MASS_PROFILE

  //============================================
  //updating indexmax & indexmin to analyse just
  //the component from paramfile
  //============================================
  
  nDisk = indexmin = indexmax = 0;
  
  indices(component);
  
  printf("\n --> Computing mass profile at different radii\n");
  double rcomp;
  for(i=0; i<3; i++)
    {
      if(flag==1)
	{
	  rcomp  = factor_comp[i] * r_scalength[0];
	}
      if(flag==2)
	{
	  rcomp  = factor_comp[i] * r_scalength[1];
	}
      printf("rcomp = %lf\n",rcomp);
      massprofile(rcomp, indexmin, indexmax,i);
    }
 
#endif
  
  //=================================================
  //SFR integrated over time (in each snap) SATELLITE
  //=================================================
  
#ifdef INTEGRATED_SFR
  
  indexmin = indexmax = nDisk = 0;
  
  indices(4);

  printf("indexmin = %d \t indexmax = %d\n",indexmin,indexmax);
  
  
  printf("\n --> Computing integrated SFR at the given evolution time\n");
  
  integrated_sfr(indexmin, indexmax);
  
#endif

#ifdef GAS_SFR

  indexmin = indexmax = nDisk = 0;
  
  indices(0);

  printf("indexmin = %d \t indexmax = %d\n",indexmin,indexmax);
  
  
  printf("\n --> Computing gas SFR at the given evolution time\n");
  
  integrated_sfr_gas(indexmin, indexmax);
  
#endif
  
  free(particles);
  
#ifdef SURF_DENS_PROFILE
  free(potential);
  free(cylindrical_radii);
  free(v_disp_R);
  free(v_radial);
  free(v_disp_Z);
  free(surf_dens);
#endif
  
#ifdef SF_ANALYSIS
  free(dens_stars);
  free(Sfr);
#endif
  
#ifdef SOFT_POTENTIAL
  free(pot_soft);
  free(r_soft);
#endif
  
  return 0;
}

int fitExponentialProfile(int nBins, double x[nBins], double y[nBins], double variance[nBins],FILE *pfLog, FILE *pfFit)
{
  int i;
  double a, b, cov00, cov01, cov11, chisq;
  double w[nBins],xf,yf,variancef;
  double h, Sigma0, diskMass;
  
  for( i=0; i<nBins; i++ )
    {
      y[i] = log10(y[i]);
      w[i] = 1.0/variance[i];
      //printf("%d %lf %lf %lf\n",i,x[i],y[i],variance[i]);
    }

  fprintf(pfLog,"Fit made with gsl_fit_wlinear\n");   

  gsl_fit_wlinear(x, 1, w, 1, y, 1, nBins, 
		  &b, &a, &cov00, &cov01, &cov11, 
		  &chisq);

  fprintf(pfLog,"best fit: Log10(Sigma) = %lf R + %lf\n", a, b);
  fprintf(pfLog,"Covariance matrix:\n");
  fprintf(pfLog,"[ %lf, %lf\n  %lf, %lf]\n",cov00, cov01, cov01, cov11);
  fprintf(pfLog,"Chi square = %lf\n", chisq);

  for(xf = 0.0; xf < x[nBins-1]; xf = xf + 0.1)
    {
      
      gsl_fit_linear_est (xf, 
                          b, a, 
                          cov00, cov01, cov11, 
                          &yf, &variancef);

      fprintf(pfFit,"%lf %lf %lf %lf\n", xf, yf, -variancef, variancef);
    }
  
  fclose(pfFit);

  h = (-1.0/a)*log10(M_E);
  //r_scalength = h;
  Sigma0 = pow(10.0,b);
  diskMass = 2.0*M_PI*h*h*pow(10.0,b);

  fprintf(pfLog,"Disk scale length h = %lf    [kpc]\n",h);
  fprintf(pfLog,"Disk Sigma0 = %lf    [10 to 10 M_sol / kpc to 3]\n",Sigma0);
  fprintf(pfLog,"Disk mass = %lf     [10 to 10 M_sol]\n\n",diskMass);

  printf("\nDisk scale length h = %lf\n",h);
  printf("Disk Sigma0 = %lf\n",Sigma0);
  printf("Disk mass = %lf\n\n",diskMass);

  return 0;
}

int binningType(int binType, int nBins, double minR, double maxR, int nDisk, double R[nDisk], int deltaN)
{

  if( binType == 1 )
    {
      radialBins( minR, maxR, nBins);
      sprintf(binning,"radialBins");
    }
  if( binType == 2 )
    {
      logRadialBins( minR, maxR, nBins);
      sprintf(binning,"logradialBins");
    }
  if( binType == 3 )
    {
      nFixedBins(nDisk, deltaN, R);
      sprintf(binning,"nfixedBins");
    }
  
  return 0;
}
