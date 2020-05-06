
int counterLines(char *infile) 
{

  int NDAT,c;
  static FILE *pf;

  if((pf=fopen(infile,"r"))==NULL)
    printf("no puedo abrir archivo %s\n",infile);

  NDAT=0;

  while((c=fgetc(pf))!= EOF)
    {
      if(c == '\n')
        {
	  ++NDAT;
        }
    }

  fclose(pf);

  return NDAT;

}

int firstDerivative(int nPoints, double x[nPoints], double fx[nPoints], double dfx[nPoints])
{

  int i;
  double meanDeltax;

  meanDeltax = ( fabs(x[1] - x[0]) + fabs(x[2] - x[1]) ) / 2.0;
  
  dfx[0] = ( -3.0*fx[0] + 4.0*fx[1] - fx[2] ) / (2.0*meanDeltax);

  meanDeltax = ( fabs(x[2] - x[1]) + fabs(x[3] - x[2]) ) / 2.0;

  dfx[1] = ( -3.0*fx[1] + 4.0*fx[2] - fx[3] ) / (2.0*meanDeltax);

  for( i=2; i<=nPoints-3; i++)
    {

      meanDeltax = ( fabs(x[i-1]-x[i-2]) + fabs(x[i]-x[i-1]) + fabs(x[i+1]-x[i]) + fabs(x[i+2]-x[i+1]) ) / 4.0;

      dfx[i] = ( -fx[i+2] + 8.0*fx[i+1] - 8.0*fx[i-1] + fx[i-2] ) / (12.0*meanDeltax);
         
    }

  meanDeltax = ( fabs(x[nPoints-3] - x[nPoints-2]) + fabs(x[nPoints-4] - x[nPoints-3]) ) / 2.0;

  dfx[nPoints-2] = ( 3.0*fx[nPoints-2] - 4.0*fx[nPoints-3] + fx[nPoints-4] ) / (2.0*meanDeltax);

  meanDeltax = ( fabs(x[nPoints-2] - x[nPoints-1]) + fabs(x[nPoints-3] - x[nPoints-2]) ) / 2.0;
  
  dfx[nPoints-1] = ( 3.0*fx[nPoints-1] - 4.0*fx[nPoints-2] + fx[nPoints-3] ) / (2.0*meanDeltax);
  
  return 0;

}

int secondDerivative(int nPoints, double x[nPoints], double fx[nPoints], double d2fx[nPoints])
{

  int i;
  double meanDeltax;

  meanDeltax = ( fabs(x[1] - x[0]) + fabs(x[2] - x[1]) + fabs(x[3] - x[2]) ) / 3.0;

  d2fx[0] = ( 2.0*fx[0] - 5.0*fx[1] + 4.0*fx[2] -fx[3]) / (meanDeltax*meanDeltax*meanDeltax);

  meanDeltax = ( fabs(x[2] - x[1]) + fabs(x[3] - x[2]) + fabs(x[4] - x[3]) ) / 3.0;

  d2fx[1] = ( 2.0*fx[1] - 5.0*fx[2] + 4.0*fx[3] -fx[4]) / (meanDeltax*meanDeltax*meanDeltax);

  for( i=2; i<=nPoints-3; i++)
    {

      meanDeltax = ( fabs(x[i-1]-x[i-2]) + fabs(x[i]-x[i-1]) + fabs(x[i+1]-x[i]) + fabs(x[i+2]-x[i+1]) ) / 4.0;

      d2fx[i] = ( -fx[i+2] + 16.0*fx[i+1] - 30.0*fx[i] + 16.0*fx[i-1] - fx[i-2] ) / (12.0*meanDeltax*meanDeltax);

    }

  meanDeltax = ( fabs(x[nPoints-3] - x[nPoints-2]) + fabs(x[nPoints-4] - x[nPoints-3]) + fabs(x[nPoints-5] - x[nPoints-4]) ) / 3.0;
  
  d2fx[nPoints-2] = ( 2.0*fx[nPoints-2] - 5.0*fx[nPoints-3] + 4.0*fx[nPoints-4] -fx[nPoints-5]) / (meanDeltax*meanDeltax*meanDeltax);

  meanDeltax = ( fabs(x[nPoints-2] - x[nPoints-1]) + fabs(x[nPoints-3] - x[nPoints-2]) + fabs(x[nPoints-4] - x[nPoints-3]) ) / 3.0;
  
  d2fx[nPoints-1] = ( 2.0*fx[nPoints-1] - 5.0*fx[nPoints-2] + 4.0*fx[nPoints-3] -fx[nPoints-4]) / (meanDeltax*meanDeltax*meanDeltax);

  
  return 0;

}

// My arctan2: returns angles in all quadrant 
double arctan2(double y, double x)
{

  if( y < 0 )
    return atan2(y,x)+2.0*M_PI;
  else
    return atan2(y,x);
}
