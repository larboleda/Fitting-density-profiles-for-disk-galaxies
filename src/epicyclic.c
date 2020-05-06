int epicyclic_med_point(double *x, double *fx)

{
  int i;
  double *k2;
  double *logR;
  double *dfx; //potential first derivative respect to R
  double *d2fx; //potential second derivative respect to R

  //R is the radial cylindrical coordinate
  
  logR= (double *)malloc((size_t)nfinalBins*sizeof(double));
  if(logR == NULL)
    {
      printf("Allocation of potInBin failed\n");
      exit(0);
    }

  dfx= (double *)malloc((size_t)nfinalBins*sizeof(double));
  if(dfx == NULL)
    {
      printf("Allocation of potInBin failed\n");
      exit(0);
    }
  
  d2fx= (double *)malloc((size_t)nfinalBins*sizeof(double));
  if(d2fx == NULL)
    {
      printf("Allocation of potInBin failed\n");
      exit(0);
    }

  k2 = (double *)malloc((size_t)nfinalBins*sizeof(double));
  if(k2 == NULL)
    {
      printf("Allocation of potInBin failed\n");
      exit(0);
    }

  //=====
  // logR
  //=====
  
  for(i=0;i<nfinalBins;i++)
    {
      logR[i] = log10(x[i]);
    }
  
  
  //==================================================
  //first derivative of potential with respect to logR
  //==================================================
  
  firstDerivative(nfinalBins, logR, potential, dfx);
  
  //===================================================
  //second derivative of potential with respect to logR
  //===================================================
  
  secondDerivative(nfinalBins, logR, potential, d2fx);

  //===================
  //epicyclic frecuency
  //===================

  for(i=0;i<nfinalBins;i++)
    {
      k2[i] = (3.0/(x[i]*x[i]))*dfx[i] + (1.0/x[i])*d2fx[i];
      printf("k2 = %lf\n",k2[i]);
    }
  return 0;
}

//int epicyclic_gsl_generic(int Npoints, double *x, double *f, double *sdens, double *vdispR, double *vdispTheta, double *vdispZ,const gsl_interp_type *a) 
int epicyclic_gsl_generic(int Npoints, double *x, double *f, double *sdens, double *vdispR, double *vdispZ, const gsl_interp_type *a) 
{

  double xi, fi, dfx, d2fx, vdispR_i, vdispZ_i, sdens_i;
  //double vdispTheta_i ;
  int i;
  double k, k2, R, Q;

  //===========
  //output file
  //===========
  
  char  out[500];
  FILE *pf;

  sprintf(out,"sp_toomre_%s_type%d.dat",binning,component);

  pf = fopen(out,"w");
 
  //gsl_sort2(x, 1, f, 1, Npoints);
  //gsl_sort(x, 1, Npoints);
  
  //==========================================================
  // gsl interpolation, derivatives of a interpolated function
  //==========================================================
  
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();

  gsl_interp *interpolator  = gsl_interp_alloc(a, (Npoints)); //pote&derivatives
  gsl_interp *interp1  = gsl_interp_alloc(a, (Npoints)); //surface density
  gsl_interp *interp2  = gsl_interp_alloc(a, (Npoints)); //R-velocity dispersion
  gsl_interp *interp4  = gsl_interp_alloc(a, (Npoints)); //Z-velocity dispersion

  /*
  gsl_interp *interp3  = gsl_interp_alloc(a, (Npoints)); //theta-velocity dispersion
   */

  gsl_interp_init(interpolator, x, f, (Npoints));
  gsl_interp_init(interp1, x, sdens, (Npoints));
  gsl_interp_init(interp2, x, vdispR, (Npoints));
  // gsl_interp_init(interp3, x, vdispTheta, (Npoints));
  gsl_interp_init(interp4, x, vdispZ, (Npoints));

  double logr;
  
  for(xi = x[0]; xi < x[Npoints-1]; xi += 0.1)
    {

      logr = log10(xi); //x = R
      
      //=================
      //interp. potential
      //=================
      
      fi=gsl_interp_eval(interpolator, x, f, xi, acc);
      
      //=====================
      //first pot. derivative
      //=====================
      
      dfx = gsl_interp_eval_deriv(interpolator, x, f, xi, acc);

      //======================
      //second pot. derivative
      //======================
      
      d2fx = gsl_interp_eval_deriv2(interpolator, x, f, xi, acc);

      //===================
      //epicyclic frecuency
      //===================
      
      k2 = (3.0/xi)*(dfx) + (d2fx);
      k = sqrt(k2);
      
      //====================================
      //interpolation of velocity dispersion
      //and surface density
      //====================================

      sdens_i = gsl_interp_eval(interp1, x, sdens, xi, acc);
      vdispR_i =gsl_interp_eval(interp2, x, vdispR, xi, acc);
      // vdispTheta_i =gsl_interp_eval(interp3, x, vdispTheta, xi, acc);
      vdispZ_i =gsl_interp_eval(interp4, x, vdispZ, xi, acc);
      
      //==================
      //Toomre's parameter
      //==================

      Q = (k*vdispR_i)/(3.36*G*sdens_i);
      
      // printf(pf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",xi,sdens_i,vdispR_i,fi,dfx,d2fx, Q, k, vdispTheta_i, vdispZ_i);
      fprintf(pf,"%16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf\n",xi,sdens_i,vdispR_i,fi,dfx,d2fx, Q, k, vdispZ_i);
    }

  
  printf("\t using %s interpolation\n",gsl_interp_name(interpolator));
  
  gsl_interp_free (interpolator);
  gsl_interp_accel_free (acc);
  
  fclose(pf);
  
  return 0;
}


//for gas, Q depends on r in general: speeds of sound depends on temperature.

int epicyclic_gsl_generic_gas(int Npoints, double *x, double *f, double *sdens, double *vdispR, double *vdispZ, double *temp, double *cs, const gsl_interp_type *a) 
{

  double xi, fi, dfx, d2fx, vdispR_i, vdispZ_i,sdens_i, temp_i, cs_i;
  int i;
  double k, k2, R, Q;

  //===========
  //output file
  //===========
  
  char  out[500];
  FILE *pf;

  sprintf(out,"sp_toomre_%s_type%d.dat",binning,component);

  pf = fopen(out,"w");
 
  //gsl_sort2(x, 1, f, 1, Npoints);
  //gsl_sort(x, 1, Npoints);
  
  //==========================================================
  // gsl interpolation, derivatives of a interpolated function
  //==========================================================
  
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();

  gsl_interp *interpolator  = gsl_interp_alloc(a, (Npoints)); //pote&derivatives
  gsl_interp *interp1  = gsl_interp_alloc(a, (Npoints)); //surface density
  gsl_interp *interp2  = gsl_interp_alloc(a, (Npoints)); //R-velocity dispersion
  gsl_interp *interp3  = gsl_interp_alloc(a, (Npoints)); //temperature
  gsl_interp *interp4  = gsl_interp_alloc(a, (Npoints)); //sound speed
  gsl_interp *interp5  = gsl_interp_alloc(a, (Npoints)); //Z-velocity dispersion
 
 
  gsl_interp_init(interpolator, x, f, (Npoints));
  gsl_interp_init(interp1, x, sdens, (Npoints));
  gsl_interp_init(interp2, x, vdispR, (Npoints));
  gsl_interp_init(interp3, x, temp, (Npoints));
  gsl_interp_init(interp4, x, cs, (Npoints));
  gsl_interp_init(interp5, x, vdispZ, (Npoints));
 
  double logr;
  
  for(xi = x[0]; xi < x[Npoints-1]; xi += 0.1)
    {

      logr = log10(xi); //x = R
      
      //=================
      //interp. potential
      //=================
      
      fi=gsl_interp_eval(interpolator, x, f, xi, acc);
      
      //=====================
      //first pot. derivative
      //=====================
      
      dfx = gsl_interp_eval_deriv(interpolator, x, f, xi, acc);

      //======================
      //second pot. derivative
      //======================
      
      d2fx = gsl_interp_eval_deriv2(interpolator, x, f, xi, acc);

      //===================
      //epicyclic frecuency
      //===================
      
      k2 = (3.0/xi)*(dfx) + (d2fx);
      k = sqrt(k2);
      
      //====================================
      //interpolation of velocity dispersion
      //and surface density
      //====================================

      sdens_i = gsl_interp_eval(interp1, x, sdens, xi, acc);
      vdispR_i = gsl_interp_eval(interp2, x, vdispR, xi, acc);
      temp_i = gsl_interp_eval(interp3, x, temp, xi, acc);
      cs_i = gsl_interp_eval(interp4, x, cs, xi, acc);
      vdispZ_i = gsl_interp_eval(interp5, x, vdispZ, xi, acc);

      //==================
      //Toomre's parameter
      //==================

      Q = (k*cs_i)/(M_PI*G*sdens_i);

      //printf("%16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf\n",xi,sdens_i,vdispR_i,fi,dfx,d2fx, Q, k2,temp_i,cs_i);
      fprintf(pf,"%16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf\n",xi,sdens_i,vdispR_i,fi,dfx,d2fx, Q, k,temp_i,cs_i,vdispZ_i);
    }

  
  printf("\t using %s interpolation\n",gsl_interp_name(interpolator));
  
  gsl_interp_free (interpolator);
  gsl_interp_accel_free (acc);
  
  fclose(pf);
  
  return 0;
}
