int Bspline(int N, double *R, double *POT)
{

  int NCOEFFS = 25; //number of fit coefficients 
  int NBREAK = NCOEFFS -2; //nbreak = ncoeffs + 2 - k = ncoeffs - 2 since k = 4 
  
  const size_t n = N; // number of data points to fit 
  const size_t ncoeffs = NCOEFFS;
  const size_t nbreak = NBREAK;
  size_t i, j;
  gsl_bspline_workspace *bw;
  gsl_vector *B;
  double dy;
  gsl_rng *r;
  gsl_vector *c, *w;
  gsl_vector *x, *y;
  gsl_matrix *XX, *cov;
  gsl_multifit_linear_workspace *mw;
  double chisq, Rsq, dof, tss, xmin, xmax;
  
    
  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);
  
  //==========================================
  //allocate a cubic bspline workspace (k = 4)
  //==========================================
  
  bw = gsl_bspline_alloc(4, nbreak);
  B = gsl_vector_alloc(ncoeffs);
  
  x = gsl_vector_alloc(n);
  y = gsl_vector_alloc(n);
  XX = gsl_matrix_alloc(n, ncoeffs);
  c = gsl_vector_alloc(ncoeffs);
  w = gsl_vector_alloc(n);
  cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
  mw = gsl_multifit_linear_alloc(n, ncoeffs);

  //===============================================================
  //data to be fitted: arrays R and POT, used to define gsl vectors 
  //===============================================================
  
  for(i=0; i<n; ++i)
    {
      double sigma; 
      
      sigma = 1;
      
      if(i==0)
	{
          xmin=R[i];
        }
      if(i==(n-1))
	xmax=R[i];
      
      gsl_vector_set(x, i, R[i]);
      gsl_vector_set(y, i, POT[i]);
      gsl_vector_set(w, i, 1.0 / (sigma * sigma));
      
    }

  

  //==================================
  //use uniform breakpoints on [0, 15]
  //==================================
  
  gsl_bspline_knots_uniform(xmin, xmax, bw);

  
  //============================
  //construct the fit matrix XX
  //============================
  
  for(i=0; i<n; ++i)
    {
      double xi = gsl_vector_get(x, i);
      
      /* compute B_j(xi) for all j */
      gsl_bspline_eval(xi, B, bw);
      
      /* fill in row i of X */
      for(j=0; j<ncoeffs; ++j)
        {
          double Bj = gsl_vector_get(B, j);
          gsl_matrix_set(XX, i, j, Bj);
        }
    }
 
  //===========
  // do the fit
  //===========
  
  gsl_multifit_wlinear(XX, w, y, c, cov, &chisq, mw);
  
  dof = n - ncoeffs;
  tss = gsl_stats_wtss(w->data, 1, y->data, 1, y->size);
  Rsq = 1.0 - chisq / tss;
  
  fprintf(stderr, "\t chisq/dof = %e, Rsq = %f\n", chisq / dof, Rsq);
  printf("\n\n");
  
  double dx=0.1;
 
  npot = (int)((xmax - xmin) / dx);
  
  pot_soft = (double *) malloc((size_t)npot*sizeof(double)); //global var.
  r_soft   = (double *) malloc((size_t)npot*sizeof(double)); //global var.

 
  //==========================
  // output the smoothed curve
  //==========================
  {
    double xi, yi, yprev, ynext, yerr;
    double min1, min2;
    int FLAG;
    
    
    yprev=-1.0e10;
    FLAG=0;
    int count = 0;
    for(xi=xmin; count<npot; xi+=dx)
      {
	
        gsl_bspline_eval(xi, B, bw);
	gsl_multifit_linear_est(B, c, cov, &yi, &yerr);

	if(count >= npot)
	  printf("pendejo!\n");
	
	pot_soft[count] = yi;
	r_soft[count] = xi;
	
	//printf("%lf %lf\n", r_soft[count], pot_soft[count]);

	count += 1;
	
      }
  }
  
  gsl_rng_free(r);
  gsl_bspline_free(bw);
  gsl_vector_free(B);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_matrix_free(XX);
  gsl_vector_free(c);
  gsl_vector_free(w);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(mw);
  
  return 0;
  
}
