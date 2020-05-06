int sfr(char * infile, char *binning)
{
  int indexmin, indexmax, i;
  
  //Index according to number of particles in each
  //component (sat or hots) of the galaxy being used
  
  //==============================
  //isolated galaxy from snap
  //==============================
  if(flag==0)
    {
      if(N_part[4]>0)
	{
	  int nStars;
	  indexmin = indexmax = 0;
	  
	  nStars = N_part[4]; 
	  
	  //actualizing index for counting just the new stars type
	  for( i=0; i<4; i++)
	    {
	      indexmin = indexmin + N_part[i];
	    }
	  indexmax = indexmin + N_part[4];
	  
	  printf("\n 7. Computing star formation rate\n");
	  surface_sfr(infile, indexmin, indexmax, nStars,binning);
       } 
    }
  
  //==============================
  //host galaxy from merger
  //==============================
  if(flag==1)
    {
      if(Nhost[4]>0) //sfr
       {
	 int nStars;
	 indexmin = indexmax = 0;
	 
	 nStars = Nhost[4];
	 
	 for( i=0; i<4; i++)
	   {
	     //indexmin = indexmax = indexmin + N_part[i];
	     indexmin = indexmin + Nhost[i];
	   }
	 indexmax = indexmin + Nhost[4];
	 
	 printf("\n 8. Computing star formation rate\n");
	 surface_sfr(infile, indexmin, indexmax, nStars,binning);
	 
       }
     else
       {
	 printf("There'r no new formed stars in host galaxy\n");
       } 
    }
  
  //==============================
  // Sat galaxy from merger
  //==============================
  
  if(flag==2)
    {
     if(Nsat[4]>0) //sfr
       {
	 int nStars;
	 indexmin = indexmax = 0;
	 
	 nStars = Nsat[4];
	 
	 for( i=0; i<4; i++)
	   {
	     //indexmin = indexmax = indexmin + N_part[i];
	     indexmin = indexmin + Nsat[i];
	   }
	 
	 indexmax = indexmin + Nsat[4];
	 
	 //printf("Nsat[4] = %d \t indexmin = %d \t indexmax = %d\n",Nsat[4],indexmin,indexmax);
	 //exit(0);
	 
	 printf("\n 8. Computing star formation rate\n");
	 surface_sfr(infile, indexmin, indexmax, nStars,binning);
       }
     else
       {
	 printf("There'r no new formed stars in satellite galaxy\n");
       }
    }

  return 0;
}

int surface_sfr(char *infile, int indexmin, int indexmax, int nDisk, char *binning)

{
   size_t *index; // Index for particles
  int i, j, counter;  // disk particles number 
  int nBins, nParticlesBin; // Bins numbers, number of particles per bin
  int usedBins = 0; // number of bins with more than minParticles
  double *Radii, area, meanMass; // Radius for particles, area per bin, mean mass per bin
  double *masses_stars_InBin, mass_stars_InBin; // mass of particles within a bin, total mass in bin
  double *R, *Sigma, *variance; // vectors to fit the profile
  FILE *pfLog, *pfFit;
  FILE *pfProfile;
  Profile *profile; //used when fittig

  char outfiles[500];
  
  pfLog = fopen("sfr.log","w");
  
  //====================
  //ALLOCATION OF MEMORY
  //====================
  
  Radii = (double *)malloc((size_t)nDisk*sizeof(double));
  if(Radii == NULL)
    {
      printf("Allocation of Radii failed\n");
      exit(0);
    }
  
  masses_stars_InBin = (double *)malloc((size_t)nDisk*sizeof(double));
  if(masses_stars_InBin == NULL)
    {
      printf("Allocation of masses_stars_InBin failed\n");
      exit(0);
    }
  
  index = (size_t *)malloc((size_t)nDisk*sizeof(size_t));
  if(index == NULL)
    {
      printf("Allocation of index failed\n");
      exit(0);
    }
  

  //=======================================
  //estimating surface density of new stars
  //=======================================


  printf("\t indexmin =%d \t indexmax =%d\n",indexmin,indexmax);
  
  for( i=indexmin; i<indexmax; i++ )
    {
      if(flag==0)
	{
	  Radii[i-indexmin] = sqrt(particles[i].pos[X]*particles[i].pos[X] +
				   particles[i].pos[Y]*particles[i].pos[Y] );
	}
      else
	{
	   Radii[i-indexmin] = sqrt(galaxy[i].pos[X]*galaxy[i].pos[X] +
				   galaxy[i].pos[Y]*galaxy[i].pos[Y] );
	}
    }
  
  gsl_sort_index(index,Radii,1,nDisk);

  /*
  for(i=0; i<nDisk;i++)
    {
      printf("Radii[%d] = %16.8lf\n",i,Radii[i]);
    }
  
    exit(0);*/

  nBins = ceil(sqrt(1.0*nDisk)/factor);
  

  printf("\t nBins = %d \t nDisk = %d \t factor = %lf\n",nBins,nDisk,factor);
  
  if( binType == 3 ) 
    {
      nBins = ceil( nDisk/particlesPerBin );
      if( (nDisk/particlesPerBin) > 1.0*nBins )
	nBins = nBins +1;
    }

  profile = (Profile *)malloc((size_t)nBins*sizeof(Profile));
  if(profile == NULL)
    {
      printf("Allocation of profile failed\n");
      exit(0);
    }
  
  dens_stars = (double *)malloc((size_t)nBins*sizeof(double));
  if(surf_dens == NULL)
    {
      printf("Allocation of dens_stars failed\n");
      exit(0);
    }
  
  Sfr = (double *)malloc((size_t)nBins*sizeof(double));
  if(Sfr == NULL)
    {
      printf("Allocation of Sfr failed\n");
      exit(0);
    }

  if( binType == 1)
    fprintf(pfLog,"Binning type = fixed delta_R\n");
  if( binType == 2)
    fprintf(pfLog,"Binning type = Logaritmic delta_R\n");
  if( binType == 3)
    fprintf(pfLog,"Binning type = fixed delta_N\n");


  if( Radii[index[nDisk-1]] < Rmax )
    {
      Rmax = Radii[index[nDisk-1]];
      printf("----------Rmax = %lf\n",Rmax);
	
    }
  if(Radii[index[0]] > Rmax)
    {
      printf("\t Radii[index[0]] = %lf \t Rmax = %lf\n",Radii[index[0]],Rmax);
      exit(0);
    }
  
  fprintf(pfLog,"Used maximum radius  = %lf\n",Rmax);
  
  binningType(binType, nBins, Radii[index[0]], Rmax, nDisk, Radii, particlesPerBin);

  
  sprintf(outfiles,"disk_sfr_%s.output",binning);

  pfProfile = fopen(outfiles,"w");

  for( i=0; i<nBins; i++)
    {
      mass_stars_InBin = 0.0;
      nParticlesBin = 0;

      for( j=indexmin; j<indexmax; j++ )
	{
	  
	  if( (Radii[j-indexmin]>=bins[i].xini) && (Radii[j-indexmin]<bins[i].xfin) )
	    {
	      //==============================
	      //mass assignation of new stars
	      //==============================

	      if(flag==0)
		{
		  mass_stars_InBin = mass_stars_InBin + particles[j].mass;
		  masses_stars_InBin[nParticlesBin] =  particles[j].mass;
		  nParticlesBin++;
		}
	      else
		{
		  mass_stars_InBin = mass_stars_InBin + galaxy[j].mass;
		  masses_stars_InBin[nParticlesBin] =  galaxy[j].mass;
		  nParticlesBin++;
		}
	      
	    }
	}
      
      area=M_PI * (bins[i].xfin*bins[i].xfin - bins[i].xini*bins[i].xini);
      
      meanMass = mass_stars_InBin/nParticlesBin; //gass mass converted into new stars
      
      //=======================
      //surface density profile
      //=======================
      
      profile[i].R = bins[i].x; 
      profile[i].Sigma = mass_stars_InBin/area;
      profile[i].variance = meanMass;

      dens_stars[i] = mass_stars_InBin/area;

      Sfr[i] = dens_stars[i]/Header.time;
      
      if(nParticlesBin>minParticles)
	{
	  profile[i].accepted = 1;
	  fprintf(pfProfile,"%lf %lf %lf %lf %d\n",profile[i].R,profile[i].Sigma,profile[i].variance,mass_stars_InBin,nParticlesBin);

	  usedBins ++;
	}
    }

  //nfinalBins = usedBins;
 
  if( fitOn == 1  )
    {   
      sprintf(outfiles,"disk_surface_profile_%s_fit_type%d.output",binning,component);

      pfFit = fopen(outfiles,"w");
      
      R = (double *)malloc((size_t)usedBins*sizeof(double));
      if(R == NULL){
	printf("Allocation of R failed\n");
	exit(0);
      }
      
      Sigma = (double *)malloc((size_t)usedBins*sizeof(double));
      if(Sigma == NULL){
	printf("Allocation of Sigma failed\n");
	exit(0);
      }
      
      variance = (double *)malloc((size_t)usedBins*sizeof(double));
      if(variance == NULL){
	printf("Allocation of variance failed\n");
	exit(0);
      }
      
      counter = 0;
      for( i=0; i<nBins; i++)
	{
	  if(profile[i].accepted == 1)
	    {
	      R[counter] = profile[i].R;
	      Sigma[counter] = profile[i].Sigma;
	      variance[counter] = profile[i].variance;
	      counter ++;
	    }
	}
      
      fitExponentialProfile(usedBins, R, Sigma, variance, pfLog, pfFit);
      
      fclose(pfProfile);
      free(R);
      free(Sigma);
      free(variance);
    }
  
  
  
  fclose(pfLog);
  free(Radii);
  free(masses_stars_InBin);

  free(index);

  
  return 0;
}

