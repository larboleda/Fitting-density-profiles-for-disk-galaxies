
int centerMass_mergercomp(int *Index, int nindex, int limlower, int limuper,CM *centerm)
{
  int i, count;
  double mtot=0.0, cm[3], vcm[3];

  for(i=0; i<3; i++)
    {
      cm[i]=0.0;
      vcm[i]=0.0;
      (*centerm).cm[i] = 0.0;
      (*centerm).vcm[i] = 0.0;
    }
  
  if(Index==NULL)
    {
      
      printf("Computing CM with interval [%d,%d)\n",limlower,limuper);
      
      mtot = 0.0;
      for(i=limlower; i<limuper; i++)
	{
	 
	  cm[0] = cm[0] + galaxy[i].pos[0]*galaxy[i].mass;
	  cm[1] = cm[1] + galaxy[i].pos[1]*galaxy[i].mass;
	  cm[2] = cm[2] + galaxy[i].pos[2]*galaxy[i].mass;
	  
	  vcm[0] = vcm[0] + galaxy[i].vel[0]*galaxy[i].mass;
	  vcm[1] = vcm[1] + galaxy[i].vel[1]*galaxy[i].mass;
	  vcm[2] = vcm[2] + galaxy[i].vel[2]*galaxy[i].mass;
	  
	  mtot = mtot + galaxy[i].mass;
	}
      
      (*centerm).cm[0] = cm[0]/mtot;
      (*centerm).cm[1] = cm[1]/mtot;
      (*centerm).cm[2] = cm[2]/mtot;
      
      (*centerm).vcm[0] = vcm[0]/mtot;
      (*centerm).vcm[1] = vcm[1]/mtot;
      (*centerm).vcm[2] = vcm[2]/mtot;

      // printf("\n (****) cmx = %16.8f \t cmy = %16.8f \t cmz = %16.8f\n",cm[0]/mtot, cm[1]/mtot,cm[2]/mtot);
      
    }
  else
    {
      
      printf("Computing CM with index\n");
      mtot = 0.0;
      count = 0;
      for(i=0; i<nindex; i++)
	{
	 
	  cm[0] = cm[0] + galaxy[Index[i]].pos[0]*galaxy[Index[i]].mass;
	  cm[1] = cm[1] + galaxy[Index[i]].pos[1]*galaxy[Index[i]].mass;
	  cm[2] = cm[2] + galaxy[Index[i]].pos[2]*galaxy[Index[i]].mass;
	  
	  vcm[0] = vcm[0] + galaxy[Index[i]].vel[0]*galaxy[Index[i]].mass;
	  vcm[1] = vcm[1] + galaxy[Index[i]].vel[1]*galaxy[Index[i]].mass;
	  vcm[2] = vcm[2] + galaxy[Index[i]].vel[2]*galaxy[Index[i]].mass;
	  
	  mtot = mtot + galaxy[Index[i]].mass;
	  count++;
	}
      
      (*centerm).cm[0] = cm[0]/mtot;
      (*centerm).cm[1] = cm[1]/mtot;
      (*centerm).cm[2] = cm[2]/mtot;
      
      (*centerm).vcm[0] = vcm[0]/mtot;
      (*centerm).vcm[1] = vcm[1]/mtot;
      (*centerm).vcm[2] = vcm[2]/mtot;
      
      printf("with %d particles\n",count);
      
    }

  (*centerm).M = mtot;  
  
  
  return 0;
}


int compute_angmom_mergercomp(int *Index, int nindex, int limlower, int limuper,CM *centerm)
{

  int i,count;
  double lcm[3];
  
  for(i=0; i<3; i++)
    lcm[i]=0.0;
  
  if(Index==NULL)
    {
      
      printf("Computing LCM with interval (%d,%d)\n",limlower,limuper);
      
      
      for(i=limlower; i<limuper; i++)
	{
	  
	  lcm[0] = lcm[0] + (galaxy[i].pos[1]*galaxy[i].vel[2] 
			     - galaxy[i].pos[2]*galaxy[i].vel[1])*galaxy[i].mass;
	  lcm[1] = lcm[1] - (galaxy[i].pos[0]*galaxy[i].vel[2] 
			     - galaxy[i].pos[2]*galaxy[i].vel[0])*galaxy[i].mass;
	  lcm[2] = lcm[2] + (galaxy[i].pos[0]*galaxy[i].vel[1] 
			     - galaxy[i].pos[1]*galaxy[i].vel[0])*galaxy[i].mass;
	  
	}
      
      (*centerm).lcm[0] = lcm[0];
      (*centerm).lcm[1] = lcm[1];
      (*centerm).lcm[2] = lcm[2];
      
    }
  else
    {
      
      printf("Computing LCM with index\n");
      count = 0;
      for(i=0; i<nindex; i++)
	{
	  lcm[0] = lcm[0] + (galaxy[Index[i]].pos[1]*galaxy[Index[i]].vel[2] 
			     - galaxy[Index[i]].pos[2]*galaxy[Index[i]].vel[1])*galaxy[Index[i]].mass;
	  lcm[1] = lcm[1] - (galaxy[Index[i]].pos[0]*galaxy[Index[i]].vel[2] 
			     - galaxy[Index[i]].pos[2]*galaxy[Index[i]].vel[0])*galaxy[Index[i]].mass;
	  lcm[2] = lcm[2] + (galaxy[Index[i]].pos[0]*galaxy[Index[i]].vel[1] 
			     - galaxy[Index[i]].pos[1]*galaxy[Index[i]].vel[0])*galaxy[Index[i]].mass;
	  
	  count++;
	}
      
      (*centerm).lcm[0] = lcm[0];
      (*centerm).lcm[1] = lcm[1];
      (*centerm).lcm[2] = lcm[2];
      
      
      printf("with %d particles\n",count);
            
    }

  
  return 0;
}

int traslate_ind_mergercomp(CM *zero, int *Index, int nindex)
{
  int i;

  for(i=0; i<nindex; i++)
    {
      
      //printf("%d %d\n",i,Index[i]);
      
      galaxy[Index[i]].pos[0] = galaxy[Index[i]].pos[0] - (*zero).cm[0];
      galaxy[Index[i]].pos[1] = galaxy[Index[i]].pos[1] - (*zero).cm[1]; 
      galaxy[Index[i]].pos[2] = galaxy[Index[i]].pos[2] - (*zero).cm[2]; 
      
      galaxy[Index[i]].vel[0] = galaxy[Index[i]].vel[0] - (*zero).vcm[0];  
      galaxy[Index[i]].vel[1] = galaxy[Index[i]].vel[1] - (*zero).vcm[1];  
      galaxy[Index[i]].vel[2] = galaxy[Index[i]].vel[2] - (*zero).vcm[2];

      //printf("%f %f %f\n",galaxy[Index[i]].pos[0],galaxy[Index[i]].pos[1],galaxy[Index[i]].pos[2]);
      
    }

  return 0;

}


int traslate_mergercomp(CM *zero, int n)
{
  int i;

  for(i=0; i<n; i++)
    {
      
      //printf("%d %d\n",i,i);
      
      galaxy[i].pos[0] = galaxy[i].pos[0] - (*zero).cm[0];
      galaxy[i].pos[1] = galaxy[i].pos[1] - (*zero).cm[1]; 
      galaxy[i].pos[2] = galaxy[i].pos[2] - (*zero).cm[2]; 
      
      galaxy[i].vel[0] = galaxy[i].vel[0] - (*zero).vcm[0];  
      galaxy[i].vel[1] = galaxy[i].vel[1] - (*zero).vcm[1];  
      galaxy[i].vel[2] = galaxy[i].vel[2] - (*zero).vcm[2];

      //printf("%f %f %f\n",galaxy[i].pos[0],galaxy[i].pos[1],galaxy[i].pos[2]);
      
    }

  return 0;

}


//int rotate(CM *vector, int *Index, int nindex)
int rotate_id_mergercomp(float *lcm, int *Index, int nindex)
{
  
  
  int m;
  double thetar=0.0, betar=0.0, r, theta, beta;
  double M1[3][3], M2[3][3];
  float vec1[3], vec2[3], lcmr[3];
  
  
  /* en el plano y,z  --> rotacion al rededor de x*/
  
  r = sqrt( lcm[1]*lcm[1] + lcm[2]*lcm[2] );
  theta = acos(fabs(lcm[1])/r);
  theta = theta*180.0/M_PI;
    
  if((lcm[1] >= 0) && (lcm[2] >= 0))
    thetar = theta + 270.0;
  
  if((lcm[1] >= 0) && (lcm[2] < 0))
    thetar = 180.0 + (90.0 - theta);
  
  if((lcm[1] < 0) && (lcm[2] < 0))
    thetar = theta + 90.0;
  
  if((lcm[1] < 0) && (lcm[2] >= 0))
    thetar = 90.0 - theta;
  
  printf("theta=%lf thetar=%lf\n",theta,thetar);
  thetar  = thetar*M_PI/180.0;
  
  //rotate about x 
  M1[0][0] = 1.0;         M1[0][1] = 0.0;               M1[0][2] = 0.0;
  M1[1][0] = 0.0;         M1[1][1] = cos(thetar);       M1[1][2] = sin(thetar);
  M1[2][0] = 0.0;         M1[2][1] = -sin(thetar);      M1[2][2] = cos(thetar);
  
  
  for(m=0; m<nindex; m++)
    {
      
      vec1[0] = galaxy[Index[m]].pos[0];
      vec1[1] = galaxy[Index[m]].pos[1];
      vec1[2] = galaxy[Index[m]].pos[2];      
      
      vec2[0] = galaxy[Index[m]].vel[0];
      vec2[1] = galaxy[Index[m]].vel[1];
      vec2[2] = galaxy[Index[m]].vel[2];
      
      galaxy[Index[m]].pos[0] = M1[0][0]*vec1[0] + M1[0][1]*vec1[1] + M1[0][2]*vec1[2];
      galaxy[Index[m]].pos[1] = M1[1][0]*vec1[0] + M1[1][1]*vec1[1] + M1[1][2]*vec1[2];
      galaxy[Index[m]].pos[2] = M1[2][0]*vec1[0] + M1[2][1]*vec1[1] + M1[2][2]*vec1[2];

      galaxy[Index[m]].vel[0] = M1[0][0]*vec2[0] + M1[0][1]*vec2[1] + M1[0][2]*vec2[2];
      galaxy[Index[m]].vel[1] = M1[1][0]*vec2[0] + M1[1][1]*vec2[1] + M1[1][2]*vec2[2];
      galaxy[Index[m]].vel[2] = M1[2][0]*vec2[0] + M1[2][1]*vec2[1] + M1[2][2]*vec2[2]; 
      
    } 
  


  //////////////////////////////////////
  /*          Segunda rotacion        */
  //////////////////////////////////////   
  
  
  lcmr[0] = M1[0][0]*lcm[0] + M1[0][1]*lcm[1] + M1[0][2]*lcm[2];
  lcmr[1] = M1[1][0]*lcm[0] + M1[1][1]*lcm[1] + M1[1][2]*lcm[2];
  lcmr[2] = M1[2][0]*lcm[0] + M1[2][1]*lcm[1] + M1[2][2]*lcm[2];

  /* en el plano x,z --> rota al rededor de y */ 
  r = sqrt( lcmr[0]*lcmr[0] + lcmr[2]*lcmr[2] );
  beta = acos(fabs(lcmr[0])/r); 
  beta = beta*180.0/M_PI;
  
  if((lcmr[0] >= 0.0) && (lcmr[2] >= 0.0))
    betar = 90 - beta;
  
  if((lcmr[0] < 0.0) && (lcmr[2] >= 0.0))
    betar = beta + 270.0;
  
  if((lcmr[0] < 0.0) && (lcmr[2] < 0.0))
    betar = (90.0 - beta) + 180.0;
  
  if((lcmr[0] >= 0.0) && (lcmr[2] < 0.0))
    betar = beta + 90.0;
    
  printf("beta=%lf betar=%lf\n",beta,betar);
  betar = betar*M_PI/180.0;
  
  //rotate about Y
  M2[0][0] = cos(betar);      M2[0][1] = 0.0;            M2[0][2] = -sin(betar);
  M2[1][0] = 0.0;             M2[1][1] = 1.0;            M2[1][2] = 0.0;
  M2[2][0] = sin(betar);      M2[2][1] = 0.0;            M2[2][2] = cos(betar);
  

  for(m=0; m<nindex; m++)
    {
      
      vec1[0] = galaxy[Index[m]].pos[0];
      vec1[1] = galaxy[Index[m]].pos[1];
      vec1[2] = galaxy[Index[m]].pos[2];
      
      vec2[0] = galaxy[Index[m]].vel[0];
      vec2[1] = galaxy[Index[m]].vel[1];
      vec2[2] = galaxy[Index[m]].vel[2]; 
      
      galaxy[Index[m]].pos[0] = M2[0][0]*vec1[0] + M2[0][1]*vec1[1] + M2[0][2]*vec1[2];
      galaxy[Index[m]].pos[1] = M2[1][0]*vec1[0] + M2[1][1]*vec1[1] + M2[1][2]*vec1[2];
      galaxy[Index[m]].pos[2] = M2[2][0]*vec1[0] + M2[2][1]*vec1[1] + M2[2][2]*vec1[2];

      galaxy[Index[m]].vel[0] = M2[0][0]*vec2[0] + M2[0][1]*vec2[1] + M2[0][2]*vec2[2];
      galaxy[Index[m]].vel[1] = M2[1][0]*vec2[0] + M2[1][1]*vec2[1] + M2[1][2]*vec2[2];
      galaxy[Index[m]].vel[2] = M2[2][0]*vec2[0] + M2[2][1]*vec2[1] + M2[2][2]*vec2[2];
      
    } 
  

  return 0;
  
}

int rotate_mergercomp(float *lcm, int nindex)
{
  
  
  int m;
  double thetar=0.0, betar=0.0, r, theta, beta;
  double M1[3][3], M2[3][3];
  float vec1[3], vec2[3], lcmr[3];
  
  
  /* en el plano y,z  --> rotacion al rededor de x*/
  
  r = sqrt( lcm[1]*lcm[1] + lcm[2]*lcm[2] );
  theta = acos(fabs(lcm[1])/r);
  theta = theta*180.0/M_PI;
    
  if((lcm[1] >= 0) && (lcm[2] >= 0))
    thetar = theta + 270.0;
  
  if((lcm[1] >= 0) && (lcm[2] < 0))
    thetar = 180.0 + (90.0 - theta);
  
  if((lcm[1] < 0) && (lcm[2] < 0))
    thetar = theta + 90.0;
  
  if((lcm[1] < 0) && (lcm[2] >= 0))
    thetar = 90.0 - theta;
  
  printf("theta=%lf thetar=%lf\n",theta,thetar);
  thetar  = thetar*M_PI/180.0;
  
  //rotate about x 
  M1[0][0] = 1.0;         M1[0][1] = 0.0;               M1[0][2] = 0.0;
  M1[1][0] = 0.0;         M1[1][1] = cos(thetar);       M1[1][2] = sin(thetar);
  M1[2][0] = 0.0;         M1[2][1] = -sin(thetar);      M1[2][2] = cos(thetar);
  
  
  for(m=0; m<nindex; m++)
    {
      
      vec1[0] = galaxy[m].pos[0];
      vec1[1] = galaxy[m].pos[1];
      vec1[2] = galaxy[m].pos[2];      
      
      vec2[0] = galaxy[m].vel[0];
      vec2[1] = galaxy[m].vel[1];
      vec2[2] = galaxy[m].vel[2];
      
      galaxy[m].pos[0] = M1[0][0]*vec1[0] + M1[0][1]*vec1[1] + M1[0][2]*vec1[2];
      galaxy[m].pos[1] = M1[1][0]*vec1[0] + M1[1][1]*vec1[1] + M1[1][2]*vec1[2];
      galaxy[m].pos[2] = M1[2][0]*vec1[0] + M1[2][1]*vec1[1] + M1[2][2]*vec1[2];

      galaxy[m].vel[0] = M1[0][0]*vec2[0] + M1[0][1]*vec2[1] + M1[0][2]*vec2[2];
      galaxy[m].vel[1] = M1[1][0]*vec2[0] + M1[1][1]*vec2[1] + M1[1][2]*vec2[2];
      galaxy[m].vel[2] = M1[2][0]*vec2[0] + M1[2][1]*vec2[1] + M1[2][2]*vec2[2]; 
      
    } 
  


  //////////////////////////////////////
  /*          Segunda rotacion        */
  //////////////////////////////////////   
  
  
  lcmr[0] = M1[0][0]*lcm[0] + M1[0][1]*lcm[1] + M1[0][2]*lcm[2];
  lcmr[1] = M1[1][0]*lcm[0] + M1[1][1]*lcm[1] + M1[1][2]*lcm[2];
  lcmr[2] = M1[2][0]*lcm[0] + M1[2][1]*lcm[1] + M1[2][2]*lcm[2];

  /* en el plano x,z --> rota al rededor de y */ 
  r = sqrt( lcmr[0]*lcmr[0] + lcmr[2]*lcmr[2] );
  beta = acos(fabs(lcmr[0])/r); 
  beta = beta*180.0/M_PI;
  
  if((lcmr[0] >= 0.0) && (lcmr[2] >= 0.0))
    betar = 90 - beta;
  
  if((lcmr[0] < 0.0) && (lcmr[2] >= 0.0))
    betar = beta + 270.0;
  
  if((lcmr[0] < 0.0) && (lcmr[2] < 0.0))
    betar = (90.0 - beta) + 180.0;
  
  if((lcmr[0] >= 0.0) && (lcmr[2] < 0.0))
    betar = beta + 90.0;
    
  printf("beta=%lf betar=%lf\n",beta,betar);
  betar = betar*M_PI/180.0;
  
  //rotate about Y
  M2[0][0] = cos(betar);      M2[0][1] = 0.0;            M2[0][2] = -sin(betar);
  M2[1][0] = 0.0;             M2[1][1] = 1.0;            M2[1][2] = 0.0;
  M2[2][0] = sin(betar);      M2[2][1] = 0.0;            M2[2][2] = cos(betar);
  

  for(m=0; m<nindex; m++)
    {
      
      vec1[0] = galaxy[m].pos[0];
      vec1[1] = galaxy[m].pos[1];
      vec1[2] = galaxy[m].pos[2];
      
      vec2[0] = galaxy[m].vel[0];
      vec2[1] = galaxy[m].vel[1];
      vec2[2] = galaxy[m].vel[2]; 
      
      galaxy[m].pos[0] = M2[0][0]*vec1[0] + M2[0][1]*vec1[1] + M2[0][2]*vec1[2];
      galaxy[m].pos[1] = M2[1][0]*vec1[0] + M2[1][1]*vec1[1] + M2[1][2]*vec1[2];
      galaxy[m].pos[2] = M2[2][0]*vec1[0] + M2[2][1]*vec1[1] + M2[2][2]*vec1[2];

      galaxy[m].vel[0] = M2[0][0]*vec2[0] + M2[0][1]*vec2[1] + M2[0][2]*vec2[2];
      galaxy[m].vel[1] = M2[1][0]*vec2[0] + M2[1][1]*vec2[1] + M2[1][2]*vec2[2];
      galaxy[m].vel[2] = M2[2][0]*vec2[0] + M2[2][1]*vec2[1] + M2[2][2]*vec2[2];
      
    } 
  

  return 0;
  
}

int rotationInclinationPosition_mergercomp(double inclination, double position, int nParticles)
{

  int i;

  double auxVector[3], auxVector2[3];

  for( i=0; i<nParticles; i++ )
    {
      
      // Clockwise rotation of host galaxy its angle inclination around X  
      auxVector[X] = galaxy[i].pos[X];
      auxVector[Y] = galaxy[i].pos[Y]*cos(inclination) + galaxy[i].pos[Z]*sin(inclination);
      auxVector[Z] = -galaxy[i].pos[Y]*sin(inclination) + galaxy[i].pos[Z]*cos(inclination);
      
      auxVector2[X] = galaxy[i].vel[X];
      auxVector2[Y] = galaxy[i].vel[Y]*cos(inclination) + galaxy[i].vel[Z]*sin(inclination);
      auxVector2[Z] = -galaxy[i].vel[Y]*sin(inclination) + galaxy[i].vel[Z]*cos(inclination);

      // Counterclockwise rotation of host galaxy its angle position around Z  
      galaxy[i].pos[X] = auxVector[X]*cos(position) - auxVector[Y]*sin(position); 
      galaxy[i].pos[Y] = auxVector[X]*sin(position) + auxVector[Y]*cos(position); 
      galaxy[i].pos[Z] = auxVector[Z];

      galaxy[i].vel[X] = auxVector2[X]*cos(position) - auxVector2[Y]*sin(position); 
      galaxy[i].vel[Y] = auxVector2[X]*sin(position) + auxVector2[Y]*cos(position); 
      galaxy[i].vel[Z] = auxVector2[Z];
      
    }

  return 0;
}
