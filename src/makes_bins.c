
int radialBins(double xMin, double xMax, int nBins )
{

  int i;
  double deltaX;

  bins =  (Bin *)malloc((size_t)nBins*sizeof(Bin));
  if(bins == NULL){
    printf("Allocation of bins failed\n");
    exit(0);
  }

  deltaX = (xMax - xMin)/nBins;  

  bins[0].xini = xMin;
  bins[0].xfin = bins[0].xini + deltaX;
  bins[0].x = (bins[0].xini + bins[0].xfin)/2.0;

  for( i=1; i<nBins; i++)
    {
     
      bins[i].xini = bins[i-1].xfin;
      bins[i].xfin = bins[i].xini + deltaX;
      bins[i].x = (bins[i].xini + bins[i].xfin)/2.0;

    }

  return 0;
}

int logRadialBins(double xMin, double xMax, int nBins )
{

  int i;
  double deltaX;
  double xLogIni,xLogFin;
  
  bins =  (Bin *)malloc((size_t)nBins*sizeof(Bin));
  if(bins == NULL){
    printf("Allocation of bins failed\n");
    exit(0);
  }

  xLogIni = log10(xMin);
  xLogFin = log10(xMax);

  deltaX = ( xLogFin - xLogIni )/nBins;

  bins[0].xini = xLogIni;
  bins[0].xfin = bins[0].xini + deltaX;
  bins[0].x = (bins[0].xini + bins[0].xfin)/2.0;

 for( i=1; i<nBins; i++)
   {

     bins[i].xini = bins[i-1].xfin;
     bins[i].xfin = bins[i].xini + deltaX;
     bins[i].x = (bins[i].xini + bins[i].xfin)/2.0;

   }

for( i=0; i<nBins; i++)
   {

     bins[i].xini = pow( 10.0, bins[i].xini );
     bins[i].xfin = pow( 10.0, bins[i].xfin );
     bins[i].x = pow( 10.0, bins[i].x );

   }

 return 0;
}

int nFixedBins(int  nTotalPoints, int deltaN, double R[nTotalPoints])
{

  size_t *index;
  int i, nBins, nIni, nFin;

  nBins = ceil(nTotalPoints/deltaN);
  if( (1.0*nTotalPoints/deltaN) > 1.0*nBins )
    nBins = nBins +1;

  bins =  (Bin *)malloc((size_t)nBins*sizeof(Bin));
  if(bins == NULL){
    printf("Allocation of bins failed\n");
    exit(0);
  }

  index = (size_t *)malloc((size_t)nTotalPoints*sizeof(size_t));
  if(index == NULL){
    printf("Allocation of index failed\n");
    exit(0);
  }

  gsl_sort_index(index,R,1,nTotalPoints);

  printf("deltaN = %d - nBins = %d\n",deltaN, nBins);

  nIni = 0;

  for( i=0; i<nBins; i++ )
    {
      
      nFin = nIni + deltaN;
      if( nFin > nTotalPoints )
	nFin = nTotalPoints-1;
      
      bins[i].xini = R[index[nIni]];  
      bins[i].xfin = R[index[nFin]];

      bins[i].indexIni = nIni;
      bins[i].indexFin = nFin;

      bins[i].x = (bins[i].xini + bins[i].xfin)/2.0;

      printf("%d  %d  %d  %d %.8lf %.8lf\n",i+1,nIni,nFin,nFin-nIni,bins[i].xini,bins[i].xfin);

      nIni = nFin;
    
    }

  return 0;
}

int nLogBins(int  nTotalPoints, int deltaN, double R[nTotalPoints])
{

  size_t *index;
  int i, nBins, nIni, nFin;

  nBins = ceil(nTotalPoints/deltaN);
  if( (1.0*nTotalPoints/deltaN) > 1.0*nBins )
    nBins = nBins +1;

  bins =  (Bin *)malloc((size_t)nBins*sizeof(Bin));
  if(bins == NULL){
    printf("Allocation of bins failed\n");
    exit(0);
  }

  index = (size_t *)malloc((size_t)nTotalPoints*sizeof(size_t));
  if(index == NULL){
    printf("Allocation of index failed\n");
    exit(0);
  }

  gsl_sort_index(index,R,1,nTotalPoints);

  //printf("deltaN = %d\n",deltaN);

  nIni = 0;

  for( i=0; i<nBins; i++ )
    {
      
      nFin = nIni + deltaN;
      if( nFin > nTotalPoints )
	nFin = nTotalPoints-1;


      bins[i].xini = R[index[nIni]];  
      bins[i].xfin = R[index[nFin]];

      bins[i].x = (bins[i].xini + bins[i].xfin)/2.0;

      // printf("%d  %d  %d  %d\n",i+1,nIni,nFin,nFin-nIni);

      nIni = nFin;
    
    }

  return 0;
}
