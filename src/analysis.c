int surface_density_profile(char *infile, int indexmin, int indexmax, int nDisk, char *binning)

{
  size_t *index; // Index for particles
  int i, j, counter;  // disk particles number 
  int nBins, nParticlesBin; // Bins numbers, number of particles per bin
  int usedBins = 0; // number of bins with more than minParticles
  double *Radii, area, meanMass, meanPotentialInBin; // Radius for particles, area per bin, mean mass per bin
  double *massesInBin,*potentialsInBin, *vR_InBin,*vZ_InBin,*vtheta_InBin, massInBin, pot; // mass of particles within a bin, total mass in bin
  double sumvR,sumvtheta,sumvZ;
  double *R, *Sigma, *variance; // vectors to fit the profile
  FILE *pfLog, *pfFit;
  FILE *pfProfile;
  Profile *profile; //used when fittig

  
  char  outfiles[500];
                                                                                                                                                    
  sprintf(outfiles,"disk_surface_profile_type%d.log",component);

  pfLog = fopen(outfiles,"w");
  
  fprintf(pfLog,"Surface density profile for : %s\n",infile);
  fprintf(pfLog,"Type Gadget = %d\n",component);
  if( symmetry == 0 )
    fprintf(pfLog,"Symmetry = cylindrical\n");
  if( symmetry == 1 )
    fprintf(pfLog,"Symmetry = spherical\n");
  fprintf(pfLog,"Factor bins number = %lf\n",factor);
  fprintf(pfLog,"Considered maximum radius  = %lf\n",Rmax);
  fprintf(pfLog,"Considered maximum height  = %lf\n",zmax);
  fprintf(pfLog,"Min particles per bin = %d\n",minParticles);
  fprintf(pfLog,"Particles number per bin = %d\n",particlesPerBin);
  fprintf(pfLog,"Binning type = %d\n",binType);
  fprintf(pfLog,"Flag to do fit = %d\n",fitOn);
  fprintf(pfLog,"Particles number = %d\n",nDisk);
  
  //====================
  //ALLOCATION OF MEMORY
  //====================
  
  Radii = (double *)malloc((size_t)nDisk*sizeof(double));
  if(Radii == NULL)
    {
      printf("Allocation of Radii failed\n");
      exit(0);
    }
  
  massesInBin = (double *)malloc((size_t)nDisk*sizeof(double));
  if(massesInBin == NULL)
    {
      printf("Allocation of massesInBin failed\n");
      exit(0);
    }

  potentialsInBin = (double *)malloc((size_t)nDisk*sizeof(double));
  if(potentialsInBin == NULL)
    {
      printf("Allocation of potInBin failed\n");
      exit(0);
    }

  vR_InBin = (double *)malloc((size_t)nDisk*sizeof(double));
  if(vR_InBin == NULL)
    {
      printf("Allocation of massesInBin failed\n");
      exit(0);
    }
  /*
  vtheta_InBin = (double *)malloc((size_t)nDisk*sizeof(double));
  if(vtheta_InBin == NULL)
    {
      printf("Allocation of vtheta in bin failed\n");
      exit(0);
    }*/

  vZ_InBin = (double *)malloc((size_t)nDisk*sizeof(double));
  if(vZ_InBin == NULL)
    {
      printf("Allocation of vz in bin failed\n");
      exit(0);
      }
  
  index = (size_t *)malloc((size_t)nDisk*sizeof(size_t));
  if(index == NULL)
    {
      printf("Allocation of index failed\n");
      exit(0);
    }
  
  
  fprintf(pfLog,"Index range = [ %d , %d )\n",indexmin,indexmax);

  //==========================
  //estimating surface density
  //==========================

  FILE  *pfile=NULL;
  char buf[200];
  sprintf(buf,"pot_type%d.dat",component);
  pfile = fopen(buf,"w");

  printf("\t indexmin =%d \t indexmax =%d\n",indexmin,indexmax);
  
  for( i=indexmin; i<indexmax; i++ )
    {
      if(flag==0)
	{
	  Radii[i-indexmin] = sqrt(particles[i].pos[X]*particles[i].pos[X] +
				   		       particles[i].pos[Y]*particles[i].pos[Y] );
	  fprintf(pfile,"%d %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf\n",particles[i].id, particles[i].pos[X], particles[i].pos[Y], particles[i].pos[Z],Radii[i-indexmin],particles[i].pot);
	}
      else
	{
	  Radii[i-indexmin] = sqrt(galaxy[i].pos[X]*galaxy[i].pos[X] +	
				   galaxy[i].pos[Y]*galaxy[i].pos[Y] );
	  fprintf(pfile,"%d %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf\n",galaxy[i].id, galaxy[i].pos[X], galaxy[i].pos[Y], galaxy[i].pos[Z],Radii[i-indexmin],galaxy[i].pot);
	}
    }
  
  gsl_sort_index(index,Radii,1,nDisk);

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

   potential = (double *)malloc((size_t)nBins*sizeof(double));
  if(potential == NULL)
    {
      printf("Allocation of potential failed\n");
      exit(0);
    }

  v_radial = (double *)malloc((size_t)nBins*sizeof(double));
  if(v_radial == NULL)
    {
      printf("Allocation of radial velocity failed\n");
      exit(0);
    }

   v_disp_R = (double *)malloc((size_t)nBins*sizeof(double));
  if(v_disp_R == NULL)
    {
      printf("Allocation of potential failed\n");
      exit(0);
    }
  /*
  v_disp_theta = (double *)malloc((size_t)nBins*sizeof(double));
  if(v_disp_theta == NULL)
    {
      printf("Allocation of ailed\n");
      exit(0);
    }*/

  v_disp_Z = (double *)malloc((size_t)nBins*sizeof(double));
  if(v_disp_Z == NULL)
    {
      printf("Allocation of vertical velocity dispersion failed\n");
      exit(0);
      }

   surf_dens = (double *)malloc((size_t)nBins*sizeof(double));
  if(surf_dens == NULL)
    {
      printf("Allocation of potential failed\n");
      exit(0);
    }
    cylindrical_radii = (double *)malloc((size_t)nBins*sizeof(double));
  if(cylindrical_radii == NULL)
    {
      printf("Allocation of cylindrical_radii failed\n");
      exit(0);
    }

  if( binType == 1)
    fprintf(pfLog,"Binning type = fixed delta_R\n");
  if( binType == 2)
    fprintf(pfLog,"Binning type = Logaritmic delta_R\n");
  if( binType == 3)
    fprintf(pfLog,"Binning type = fixed delta_N\n");

  fprintf(pfLog,"Tested Bins = %d\n",nBins);

  if( Radii[index[nDisk-1]] < Rmax )
    {
      Rmax = Radii[index[nDisk-1]];
      printf("----------Rmax = %lf\n",Rmax);
	
    }
  
  fprintf(pfLog,"Used maximum radius  = %lf\n",Rmax);

  printf("\t Radii[index[0]] = %lf \t Rmax = %lf\n",Radii[index[0]],Rmax);
  
  binningType(binType, nBins, Radii[index[0]], Rmax, nDisk, Radii, particlesPerBin);

  
  sprintf(outfiles,"disk_surface_profile_%s_type%d.output",binning,component);

  pfProfile = fopen(outfiles,"w");

  //file for raw data
  FILE *file;
  char  name[200];
  sprintf(name,"raw_type%d.dat",component);
  file = fopen(name,"w");
      
      
  for( i=0; i<nBins; i++)
    {

      //printf("----------xini = %lf\t xfin = %lf \n",bins[i].xini,bins[i].xfin);
      massInBin = 0.0;
      pot = 0.0;
      nParticlesBin = 0;
      sumvR = 0.0;

      for( j=indexmin; j<indexmax; j++ )
	{
	  
	  //printf("----------R =%lf \t xini = %lf\t xfin = %lf \n",Radii[j-indexmin],bins[i].xini,bins[i].xfin);
	  
	  if( (Radii[j-indexmin]>=bins[i].xini) && (Radii[j-indexmin]<bins[i].xfin) )
	    {
	      //============================
	      //mass assignation & potential
	      //============================
	      if(flag==0)
		{
		  massInBin = massInBin + particles[j].mass;
		  massesInBin[nParticlesBin] =  particles[j].mass;
		  pot = pot + particles[j].pot;
		  potentialsInBin[nParticlesBin] = particles[j].pot;
		  //=========================
	      //mean R velocity component
	      //=========================

	      vR_InBin[nParticlesBin] =   (1.0/Radii[j-indexmin]) * (particles[j].pos[0]*particles[j].vel[0] + particles[j].pos[1]*particles[j].vel[1] );

	      sumvR = sumvR + vR_InBin[nParticlesBin];
	      
	      /*
	      //=============================
	      //mean theta velocity component
	      //=============================
	      
	      theta = arctan2(particles[j].pos[1], particles[j].pos[0]);
	      
	      vtheta_InBin[nParticlesBin] = Radii[i-indexmin] * (cos(theta)*cos(theta)*((-particles[j].vel[0]*particles[j].pos[1]/(particles[j].pos[0]*particles[j].pos[0])) + (particles[j].vel[1]/particles[j].pos[0])));
	      
	      sumvtheta = sumvtheta + vtheta_InBin[nParticlesBin];*/
	      
	      //=========================
	      //mean Z velocity component
	      //=========================
	      
	      vZ_InBin[nParticlesBin] = particles[j].vel[2];
	       
	      sumvZ = sumvZ + vZ_InBin[nParticlesBin];
	      
	      nParticlesBin++;
		}
	      
	      else
		{
		  massInBin = massInBin + galaxy[j].mass;
		  massesInBin[nParticlesBin] =  galaxy[j].mass;
		  pot = pot + galaxy[j].pot;
		  potentialsInBin[nParticlesBin] = galaxy[j].pot;
		  //=========================
		  //mean R velocity component
		  //=========================
		  
		  vR_InBin[nParticlesBin] =   (1.0/Radii[j-indexmin]) * (galaxy[j].pos[0]*galaxy[j].vel[0] + galaxy[j].pos[1]*galaxy[j].vel[1] );
		  
		  sumvR = sumvR + vR_InBin[nParticlesBin];
		  
		  /*
		  //=============================
		  //mean theta velocity component
		  //=============================
		  
		  theta = arctan2(galaxy[j].pos[1], galaxy[j].pos[0]);
		  
		  vtheta_InBin[nParticlesBin] = Radii[i-indexmin] * (cos(theta)*cos(theta)*((-galaxy[j].vel[0]*galaxy[j].pos[1]/(galaxy[j].pos[0]*galaxy[j].pos[0])) + (galaxy[j].vel[1]/galaxy[j].pos[0])));
		  
		  sumvtheta = sumvtheta + vtheta_InBin[nParticlesBin];*/
		  
		  //=========================
		  //mean Z velocity component
		  //=========================
		  
		  vZ_InBin[nParticlesBin] = galaxy[j].vel[2];
		  
		  sumvZ = sumvZ + vZ_InBin[nParticlesBin];
		  
		  nParticlesBin++;
		}
	    }
	}
      
      area= M_PI * (bins[i].xfin*bins[i].xfin - bins[i].xini*bins[i].xini);
      
      meanMass = massInBin/nParticlesBin;
      
      meanPotentialInBin = pot/nParticlesBin;
      
      //==============================================================
      //storing mean potentials in bins to estimate epicycle frecuency
      //==============================================================
      
      potential[i] = meanPotentialInBin;
      cylindrical_radii[i] = bins[i].x;
      
      //===================
      //velocity dispersion
      //===================
      
      double meanvR = 0.0, meanvZ = 0.0;
      
      double a=0.0, c=0.0;
      
      meanvR = sumvR/nParticlesBin;

      v_radial[i] = meanvR;

      //printf("%lf\n",v_radial[i]);
      
      for(j=0;j<nParticlesBin;j++)
	{
	 
	  //a += sqrt((vR_InBin[j]-meanvR)*(vR_InBin[j]-meanvR));
	  a += (vR_InBin[j]-meanvR)*(vR_InBin[j]-meanvR);
	  //b += (vtheta_InBin[j]-meanvtheta)*(vtheta_InBin[j]-meanvtheta);
	  c += (vZ_InBin[j]-meanvZ)*(vZ_InBin[j]-meanvZ);
	}

      //==================================
      //R-velocity dispersion for each bin
      //==================================
      
      v_disp_R[i] = sqrt((1.0/(nParticlesBin-1))*a);

      //v_disp_R[i] = (1.0/(nParticlesBin-1))*a;

      /*
      //======================================
      //theta-velocity dispersion for each bin
      //======================================
      
      v_disp_theta[i] = sqrt((1.0/(nParticlesBin-1))*b);*/

      //==================================
      //R-velocity dispersion for each bin
      //==================================
      
      v_disp_Z[i] = sqrt((1.0/(nParticlesBin-1))*c);
      
      //=======================
      //surface density profile
      //=======================
      
      profile[i].R = bins[i].x; 
      profile[i].Sigma = massInBin/area;
      profile[i].variance = meanMass;

      surf_dens[i] = massInBin/area;
      
      if(nParticlesBin>minParticles)
	{
	  profile[i].accepted = 1;
	  fprintf(pfProfile,"%lf %lf %lf %lf %d\n",profile[i].R,profile[i].Sigma,profile[i].variance,massInBin,nParticlesBin);

	  usedBins ++;
	}

      //===========================
      //write raw data from binning
      //===========================
      
      fprintf(file,"%16.8lf %16.8lf %16.8lf %16.8lf 16.8lf %d\n",profile[i].R,v_disp_R[i],v_radial[i],v_disp_Z[i],profile[i].Sigma,nParticlesBin);
      
    }
  
  nfinalBins = usedBins;
  
  fprintf(pfLog,"Used bins = %d\n",usedBins);
 
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
  free(massesInBin);
  free(potentialsInBin);
  free(index);

  
  return 0;
}

//==================================================
//estimation of volumetric density in spheric shells
//==================================================


int volumetric_density(int limlower, int limuper, int nDisk, char *binning)
{

  int i, j;
  size_t *index;
  double *Radii; 
  double BinVolume, density, meanMass;
  double masaEnBin;
  int nBins, nParticlesBin; // Bins numbers, number of particles per bin
  
  FILE *laufile;

  char out[500];

  index = (size_t *)malloc((size_t)nDisk*sizeof(size_t));
  if(index == NULL){
    printf("Allocation of Index failed\n");
    exit(0);
  }
  
  Radii = (double *)malloc((size_t)nDisk*sizeof(double));
  if(Radii == NULL){
    printf("Allocation of Radii failed\n");
    exit(0);
  }
 
  //*********************
  
  for( i=limlower; i<limuper; i++ )
    if(flag==0)
      {
	Radii[i-limlower] = sqrt( particles[i].pos[0]*particles[i].pos[0] +
				  particles[i].pos[1]*particles[i].pos[1] +  particles[i].pos[2]*particles[i].pos[2]);
      }
    else
      {
	Radii[i-limlower] = sqrt(galaxy[i].pos[0]*galaxy[i].pos[0] +
				 galaxy[i].pos[1]*galaxy[i].pos[1] + galaxy[i].pos[2]*galaxy[i].pos[2]);
      }
  
  gsl_sort_index(index,Radii,1,nDisk);
  
  nBins = ceil( sqrt(nDisk)/factor );
  
  if( binType == 3 ) 
    {
      nBins = ceil( nDisk/particlesPerBin );
      if( (nDisk/particlesPerBin) > 1.0*nBins )
	nBins = nBins +1;
    }
  
  printf("*******************Rmax = %lf, Rmax_arreglo = %lf\n",Rmax, Radii[index[nDisk-1]]);

  if( Radii[index[nDisk-1]] < Rmax )
    Rmax = Radii[index[nDisk-1]];
  
  printf("Rmax = %lf\n",Rmax);
  binningType(binType, nBins, Radii[index[0]], Rmax, nDisk, Radii, particlesPerBin);

  sprintf(out,"volumetric_density_%s_type%d.dat",binning,component);

  laufile = fopen(out,"w");
  
  for( i=0; i<nBins; i++)
    {

      masaEnBin = 0.0;
      nParticlesBin = 0;

      for( j=limlower; j<limuper; j++ )
	{
	  if( (Radii[j-limlower]>=bins[i].xini) && (Radii[j-limlower]<bins[i].xfin) )
	    {
	      if(flag==0)
		{
		  masaEnBin = masaEnBin + particles[j].mass;
		}
	      else
		{
		  masaEnBin = masaEnBin + galaxy[j].mass;
		}
	      
	      nParticlesBin++;
	    }
	}

        
      BinVolume= (1.3333333333)* M_PI * ((bins[i].xfin*bins[i].xfin*bins[i].xfin) - (bins[i].xini*bins[i].xini*bins[i].xini));
      
      meanMass = masaEnBin/nParticlesBin;
      
      density = masaEnBin/BinVolume; 
         
      if(nParticlesBin>minParticles)
	{
	  
	  fprintf(laufile,"%lf %lf %lf %d %lf\n",masaEnBin,density,bins[i].x,nParticlesBin, meanMass);
	}
      //else 
      //printf("nparticlesBin=%d  menor que minParticles = %d\n",nParticlesBin, minParticles);
    }  
  
  free(Radii);
  free(index);
  fclose(laufile);
  return 0;
}
 
