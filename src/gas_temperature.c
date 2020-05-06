int temperature(int indexmin, int indexmax, int nDisk, char *binning)

{
  size_t *indice; // Indice for particles
  int i, j;  // disk particles number 
  int nBins, nPartBin; // Bins numbers, number of particles per bin
  double *Radi; // Radius for particles
  double *TempInBin, *vsoundInBin; // temperature of particles within a bin
  double meanTemp, sumtemp, sumvs, meanvs;
  //char  outfiles[500];
                                                                                
  //====================
  //ALLOCATION OF MEMORY
  //====================
  
  Radi = (double *)malloc((size_t)nDisk*sizeof(double));
  if(Radi == NULL)
    {
      printf("Allocation of Radii failed\n");
      exit(0);
    }
  
  TempInBin = (double *)malloc((size_t)nDisk*sizeof(double));
  if(TempInBin == NULL)
    {
      printf("Allocation of massesInBin failed\n");
      exit(0);
    }
  
   vsoundInBin = (double *)malloc((size_t)nDisk*sizeof(double));
  if(vsoundInBin == NULL)
    {
      printf("Allocation of massesInBin failed\n");
      exit(0);
    }

  indice = (size_t *)malloc((size_t)nDisk*sizeof(size_t));
  if(indice == NULL)
    {
      printf("Allocation of indice failed\n");
      exit(0);
    }

  //============================
  //estimating radii for binning
  //============================

  printf("\t indexmin =%d \t indexmax =%d\n",indexmin,indexmax);
  
  for( i=indexmin; i<indexmax; i++ )
    {
      if(flag == 0)
	{
	  Radi[i-indexmin] = sqrt(particles[i].pos[X]*particles[i].pos[X] +
				  particles[i].pos[Y]*particles[i].pos[Y] );
	}
      else
	{
	  Radi[i-indexmin] = sqrt(galaxy[i].pos[X]*galaxy[i].pos[X] +
				  galaxy[i].pos[Y]*galaxy[i].pos[Y] );
	}
      
    }
  
  gsl_sort_index(indice,Radi,1,nDisk);
  
  nBins = ceil(sqrt(1.0*nDisk)/factor);
  
  printf("\t nBins = %d \t nDisk = %d \t factor = %lf\n",nBins,nDisk,factor);
  
  if( binType == 3 ) 
    {
      nBins = ceil( nDisk/particlesPerBin );
      if( (nDisk/particlesPerBin) > 1.0*nBins )
	nBins = nBins +1;
    }
  
  temperatures = (double *)malloc((size_t)nBins*sizeof(double));

  if(temperatures == NULL)
    {
      printf("Allocation of potential failed\n");
      exit(0);
    }

  vsound = (double *)malloc((size_t)nBins*sizeof(double));

  if(vsound == NULL)
    {
      printf("Allocation of potential failed\n");
      exit(0);
    }

  if( Radi[indice[nDisk-1]] < Rmax )
    {
      Rmax = Radi[indice[nDisk-1]];
      printf("----------Rmax = %lf\n",Rmax);
	
    }
  
  printf("\t Radi[indice[0]] = %lf \t Rmax = %lf\n",Radi[indice[0]],Rmax);
  
  binningType(binType, nBins, Radi[indice[0]], Rmax, nDisk, Radi, particlesPerBin);
  
  //sprintf(outfiles,"temperature_%s_type%d.output",binning,component);

  //pfProfile = fopen(outfiles,"w");

  
  for( i=0; i<nBins; i++) //loop on bins
    {
      nPartBin = 0;
      meanTemp = 0.0;
      meanvs = 0.0;
      sumvs = 0.0;
      sumtemp = 0.0;
      
      for( j=indexmin; j<indexmax; j++ )
	{
	  
	  if( (Radi[j-indexmin]>=bins[i].xini) && (Radi[j-indexmin]<bins[i].xfin) )
	    {
	      
	      //=========================
	      //temperature assignation
	      //=========================
	      if(flag==0)
		{
		  TempInBin[nPartBin] = ((GAMMA_MINUS1)*((gaspro[j].U/KB)*(unitFactor)*mu));
		}
	      else
		{
		  TempInBin[nPartBin] = ((GAMMA_MINUS1)*((galaxy_gaspro[j].U/KB)*(unitFactor)*mu));
		}
	      vsoundInBin[nPartBin] = sqrt(KB*TempInBin[nPartBin]/mp);

	      sumtemp = sumtemp + TempInBin[nPartBin];
	      sumvs = sumvs + vsoundInBin[nPartBin];

	      nPartBin++;	      
	      
	    }
	}

      //=======================
      //mean temperature in bin
      //=======================
    
      meanTemp = sumtemp/nPartBin;
      meanvs = sumvs/nPartBin; // in cm/s
      
      temperatures[i] = meanTemp;
      vsound[i] = meanvs * (1.0/100000);  //from cm/s to km/s
     
    }
  
  
  free(Radi);
  free(TempInBin);
  free(vsoundInBin);
  free(indice);

  
  return 0;
}
