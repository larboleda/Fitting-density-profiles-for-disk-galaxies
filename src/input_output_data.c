
//========
//ROUTINES
//========

FILE *fileOpen(char filename[],char mode[]);
int gsl_int_int_sort(int dimension,int *fvector,int *ivector);

int read_header_gadget(char filename[])
{
  FILE *fdata;
  int dummy;
  
  
  //====================================
  //READ DATA FILE
  //====================================

  fdata = fileOpen(filename,"r");
  
  
  //====================================
  //READ HEADER
  //====================================
  
#ifdef GADGET2_FORMAT
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
  returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
#endif

  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  returnRead = fread(&Header,sizeof(io_header),1,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);

  fclose(fdata);

  return 0;
}


  //======================================================================
  //Open a file
  //======================================================================

FILE *fileOpen(char filename[],char mode[])
{
  FILE *f;

  if( !(f=fopen(filename,mode)) )
    {
      fprintf(stderr,"Error openning '%s' for %s\n",filename,mode);
      exit(1);
    }

  return f;
}

int read_gadget(char filename[])
{

  int i, j, type;
  int indexmin, indexmax;
  int nPartWithMass;
  FILE *fdata;
  int dummy;
  
  
  //====================================
  //READ DATA FILE
  //====================================
  
  fdata = fileOpen(filename,"r");
 
  FILE *filel=fopen("bfield.dat","w");
     
  //====================================
  //READ HEADER
  //====================================
  
 #ifdef GADGET2_FORMAT
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
  returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
  printf("------we've read \t %s\n",arr_name);
#endif

  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  returnRead = fread(&Header,sizeof(io_header),1,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);

  //=========================
  //TOTAL NUMBER OF PARTICLES
  //=========================
  
  N_part_total = 0;
  nPartWithMass = 0;
  printf("\nReading snapshot %s with:\n\n",filename);

  for(i=0; i<6; i++)
    {
      N_part[i] = Header.npartTotal[i];
      N_part_total += N_part[i];
      printf("%.8d particles of type %d\n",N_part[i],i);
      if( Header.mass[i]>0.0 )
	nPartWithMass = nPartWithMass +0;
      else
	nPartWithMass = nPartWithMass + Header.npartTotal[i];
    }
  printf("%.8d particles in the snapshot\n",N_part_total);
  printf("Particles with explicit mass = %.9d\n",nPartWithMass);
  
  
  //=================
  //ALLOCATE AND READ
  //=================
  
  particles = (particulas *)malloc((size_t)N_part_total*sizeof(particulas));
  if(particles == NULL){
    printf("Allocation of particles failed\n");
    exit(0);
  }
  
  if( N_part[0] > 0 )
    {
      gaspro = (gas_properties *)malloc((size_t) N_part[0]*sizeof(gas_properties));
      if(gaspro == NULL){
	printf("Allocation of gaspro failed\n");
	exit(0);
      }
    }
  
  printf("\nRead memory blocks:\n");
   
  //============
  //0. POSITIONS
  //===========
  
#ifdef GADGET2_FORMAT
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  returnRead = fread(&arr_name,sizeof(char),4,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  printf("------we've read \t %s\n",arr_name);
#endif

  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
   printf("\ndummy position                      = %.9d",dummy);
  printf(" - expected = %.9lu\n",3*sizeof(float)*N_part_total);
  
  for(i=0; i<N_part_total; i++)
    returnRead = fread(&particles[i].pos[0],sizeof(float),3,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  
  //=============
  //1. VELOCITIES
  //=============
  
#ifdef GADGET2_FORMAT
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  returnRead = fread(&arr_name,sizeof(char),4,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  printf("------we've read \t %s\n",arr_name);
#endif
  

  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  printf("dummy velocity                      = %.9d",dummy);
  printf(" - expected = %.9lu\n",3*sizeof(float)*N_part_total);  
  for(i=0; i<N_part_total; i++)
    returnRead = fread(&particles[i].vel[0],sizeof(float),3,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  
  //======
  //2. IDS
  //======
  
  int idDummy, nOldParticles;
  unsigned int UidDummy;
  nOldParticles = N_part[0]+N_part[1]+N_part[2]+N_part[3]; 
 
  #ifdef GADGET2_FORMAT
    returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
    returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
    returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
    returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
    printf("------we've read \t %s\n",arr_name);
  #endif

  returnRead = fread(&dummy,sizeof(dummy),1,fdata); 
  printf("dummy IDs                           = %.9d",dummy);
  printf(" - expected = %.9lu\n",sizeof(unsigned int)*N_part_total);

  #ifdef LONGIDS
    for(i=0; i<N_part_total; i++)
      returnRead = fread(&particles[i].Id,sizeof(unsigned long long),1,fdata);
  #else
    for(i=0; i<N_part_total; i++)
      {
	if( i<nOldParticles )
	  returnRead = fread(&particles[i].id,sizeof(unsigned int),1,fdata);
	else
	  {
	    returnRead = fread(&idDummy,sizeof(unsigned int),1,fdata);
	    particles[i].id = idDummy;
	    //	  printf("i=%.9d - id int = %.9d\n",i,idDummy);
	  }
	
      }
  #endif //LONGIDS

  returnRead = fread(&dummy,sizeof(dummy),1,fdata);

  #ifdef CODE_GIZMO

  //============
  //3. CHILD IDS
  //============

    #ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      printf("------we've read \t %s\n",arr_name);
    #endif

    returnRead = fread(&dummy,sizeof(dummy),1,fdata); 
    printf("dummy IDs                           = %.9d",dummy);
    printf(" - expected = %.9lu\n",sizeof(unsigned int)*N_part_total);
    
    #ifdef LONGIDS
      for(i=0; i<N_part_total; i++)
	returnRead = fread(&particles[i].ID_child_number,sizeof(unsigned long long),1,fdata);
    #else
      for(i=0; i<N_part_total; i++)
	{
	  if( i<nOldParticles )
	    returnRead = fread(&particles[i].ID_child_number,sizeof(unsigned int),1,fdata);
	  else
	    {
	      returnRead = fread(&idDummy,sizeof(int),1,fdata);
	      particles[i].ID_child_number = UidDummy;
	      //	  printf("i=%.9d - id int = %.9d\n",i,idDummy);
	    }
	}
    #endif //LONGIDS
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  
  //=================
  //4. GENERATION IDS
  //=================
 
    #ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      printf("------we've read \t %s\n",arr_name);
    #endif

    returnRead = fread(&dummy,sizeof(dummy),1,fdata); 
    printf("dummy IDs                           = %.9d",dummy);
    printf(" - expected = %.9lu\n",sizeof(unsigned int)*N_part_total);

    #ifdef LONGIDS
      for(i=0; i<N_part_total; i++)
	returnRead = fread(&particles[i].ID_generation,sizeof(unsigned long long),1,fdata);
    #else
      for(i=0; i<N_part_total; i++)
	{
	  if( i<nOldParticles )
	    returnRead = fread(&particles[i].ID_generation,sizeof(int),1,fdata);
	  else
	    {
	      returnRead = fread(&idDummy,sizeof(int),1,fdata);
	      particles[i].ID_generation = idDummy;
	      //	  printf("i=%.9d - id int = %.9d\n",i,idDummy);
	    }
	}
    #endif //LONGIDS
    returnRead = fread(&dummy,sizeof(dummy),1,fdata);

#endif //CODE GIZMO
   
//=========
//5. MASSES
//=========
    
  if( nPartWithMass>0  ) 
    { 
#ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      printf("------we've read \t %s\n",arr_name);
#endif
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy masses                        = %.9d",dummy);
      printf(" - expected = %.9lu\n",sizeof(float)*nPartWithMass);
    }
  
  N_min=N_max=0;
  
  for(j=0;j<=5;j++)
    {
      N_max=N_max+Header.npartTotal[j];
      printf("*******Header = %d\t N_max = %d \n",Header.npartTotal[j],N_max);

      if( Header.npartTotal[j]>0 )
	{
	  if( Header.mass[j]>0.0 )
	    {
	      for(i=N_min;i<N_max;i++)
		particles[i].mass = Header.mass[j];
	        printf("mass  = %d\n",particles[i].mass);
	    }
	  else
	    {
	      for(i=N_min;i<N_max;i++)
		{
		  returnRead = fread(&particles[i].mass,sizeof(float),1,fdata);
		} 
	    }
	  N_min=N_max;
	}

    }
  
  if( nPartWithMass>0  )  
    returnRead = fread(&dummy,sizeof(dummy),1,fdata);
 
  if(Header.npartTotal[0]!=0)
    {
      //============================================
      //7. Read internal energy for particles of gas
      //============================================
     
#ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      printf("------we've read \t %s\n",arr_name);
#endif
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy internal energy               = %.9d",dummy);
      printf(" - expected = %.9lu\n",sizeof(float)*Header.npartTotal[0]);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].U,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);

      //====================================
      //8. Read density for particles of gas
      //====================================
    
#ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      printf("------we've read \t %s\n",arr_name);
#endif
      
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy density                       = %.9d",dummy);
      printf(" - expected = %.9lu\n",sizeof(float)*Header.npartTotal[0]);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].rho,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);

#ifdef COOLING
      
      //=============================================
      //9. Read electron density for particles of gas
      //=============================================
      
#ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      printf("------we've read \t %s\n",arr_name);
#endif

      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy electron abundance            = %.9d",dummy);
      printf(" - expected = %.9lu\n",sizeof(float)*Header.npartTotal[0]);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].Ne,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);

      //======================================================
      //10. Read neutral hydrigen density for particles of gas
      //======================================================
      
#ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      printf("------we've read \t %s\n",arr_name);
#endif
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy neutral hydrogen abundance    = %.9d",dummy);
      printf(" - expected = %.9lu\n",sizeof(float)*Header.npartTotal[0]);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].Nh,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif //COOLING

      //==================================================
      //11. Read SPH smoothing length for particles of gas
      //==================================================
      
 #ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      printf("------we've read \t %s\n",arr_name);
#endif
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy hsml                          = %.9d",dummy);
      printf(" - expected = %.9lu\n",sizeof(float)*Header.npartTotal[0]);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].h,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);

      //=================================================
      //12. Read star formation rate for particles of gas
      //=================================================
      
#ifdef SFR
      
#ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      printf("------we've read \t %s\n",arr_name);
#endif
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy sfr                           = %.9d",dummy);
      printf(" - expected = %.9lu\n",sizeof(float)*Header.npartTotal[0]);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].sfr,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif
    }
  
  //=============================================
  //13. Read formation time of star for new stars
  //=============================================

  if(Header.npartTotal[4]>0)
    {
      N_min = N_part[0] + N_part[1] + N_part[2] + N_part[3];
      N_max = N_min + N_part[4];

      #ifdef STELLARAGE

      #ifdef GADGET2_FORMAT
        returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
	returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
	returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
	returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
	printf("------we've read \t %s\n",arr_name);
      #endif
    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy stellar age                   = %.9d",dummy);
      printf(" - expected = %.9lu\n",sizeof(float)*Header.npartTotal[4]);
      for(i=N_min;i<N_max;i++)
	returnRead = fread(&particles[i].stellar_age,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      #endif //STELLARAGE
    }

 //==================================
  //Read metallicity for gas and stars
  //==================================
  
  if( (Header.npartTotal[0]>0) || (Header.npartTotal[4]>0) )
    {

      #ifdef METALS 

      #ifdef GADGET2_FORMAT
        returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
	returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
	returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
	returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
	printf("------we've read \t %s\n",arr_name); 
      #endif
     
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy metallicity                   = %.9d",dummy);
      printf(" - expected = %.9lu\n",sizeof(float)*( Header.npartTotal[0] +Header.npartTotal[4] ));
       for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&particles[i].metallicity,sizeof(float),1,fdata);
      
      N_min = N_part[0] + N_part[1] + N_part[2] + N_part[3];
      N_max = N_min + N_part[4];
      
      for(i=N_min;i<N_max;i++)
	returnRead = fread(&particles[i].metallicity,sizeof(float),1,fdata);
      
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      #endif //METALS
    }

  //==================================================
  //23. Read gravitational potential for all particles
  //==================================================
  
#ifdef OUTPUTPOTENTIAL

#ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      printf("------we've read \t %s\n",arr_name);
#endif 
  
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  printf("dummy potential                     = %.9d",dummy);
  printf(" - expected = %.9lu\n",sizeof(float)*N_part_total);
  for(i=0;i<N_part_total;i++)
    returnRead = fread(&particles[i].pot,sizeof(float),1,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif

  //=======================================
  //24. Read acceleration for all particles
  //=======================================
  
#ifdef OUTPUTACCELERATION

#ifdef GADGET2_FORMAT
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
  returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
  printf("------we've read \t %s\n",arr_name);
#endif 
  
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  printf("dummy acceleration                  = %.9d",dummy);
  printf(" - expected = %.9lu\n",3*sizeof(float)*N_part_total);
  for(i=0;i<N_part_total;i++)
    returnRead = fread(&particles[i].acce[0],sizeof(float),3,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif

   //==========================================================
  //36. Read rate of change of internal energy/entropy for gas 
  //==========================================================
  
  if( Header.npartTotal[0] >0 )
    {
#ifdef OUTPUTCHANGEOFENTROPY

#ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      printf("------we've read \t %s\n",arr_name);
#endif 

      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy rate of change of entropy     = %.9d",dummy);
      printf(" - expected = %.9lu\n",sizeof(float)*Header.npartTotal[0]);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].ecr,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif //OUTPUTCHANGEOFENTROPY
    }
  
  //===================================
  //41. Read timestep for all particles
  //===================================
  
#ifdef OUTPUTTIMESTEP

#ifdef GADGET2_FORMAT
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
  returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
  printf("------we've read \t %s\n",arr_name);
#endif 
  
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  printf("dummy time step                     = %.9d",dummy);
  printf(" - expected = %.9lu\n",sizeof(float)*N_part_total);
  for(i=0;i<N_part_total;i++)
    returnRead = fread(&particles[i].timestep,sizeof(float),1,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif //OUTPUTTIMESTEP

  
  
  
  if( Header.npartTotal[0] >0 )
    {
      
#ifdef MAGNETIC_FIELD
      
      //==================
      //42. Magnetic field
      //==================
#ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      printf("------we've read \t %s\n",arr_name);
#endif 

      float Baux[3];
     
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy magnetic field                 = %.9d",dummy);
      printf(" - expected = %.9lu\n",3*sizeof(float)*Header.npartTotal[0]);
      for(i=0;i<Header.npartTotal[0];i++)
	{
	  //returnRead = fread(&gaspro[i].B[0],sizeof(float),3,fdata);
	  returnRead = fread(&Baux[0],sizeof(float),3,fdata);
	  gaspro[i].B[0] = 1.0*Baux[0];
	  gaspro[i].B[1] = 1.0*Baux[1];
	  gaspro[i].B[2] = 1.0*Baux[2];
	  
	}
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
    

      //================================
      //47. Divergence of Magnetic field
      //================================
#ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      printf("------we've read \t %s\n",arr_name);
#endif 
      
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy divergence                    = %.9d",dummy);
      printf(" - expected = %.9lu\n",sizeof(float)*Header.npartTotal[0]);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].divB,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      
      //===========================
      //50. DivBcleaningFunctionPhi
      //===========================
#ifdef GADGET2_FORMAT
        returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
	returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
	returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
	returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
	printf("------we've read \t %s\n",arr_name);
#endif 
	
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy DivBcleaningFunctionPhi                   = %.9d",dummy);
      printf(" - expected = %.9lu\n",sizeof(float)*Header.npartTotal[0]);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].Phi,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);

      //===============================
      //51. DivBcleaningFunctionGradPhi
      //===============================
      
#ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      printf("------we've read \t %s\n",arr_name);
#endif 
      
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy DivBcleaningFunctionPhi                   = %.9d",dummy);
      printf(" - expected = %.9lu\n",3*sizeof(float)*Header.npartTotal[0]);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].GradPhi[0],sizeof(float),3,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      
#endif //MAGNETC_FIELD
      
      /*
	#ifdef OUTPUT_DIV_CURL
	
	#ifdef GADGET2_FORMAT
	returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
	returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
	returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
	returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
	printf("------we've read \t %s\n",arr_name);
	#endif 
	
	returnRead = fread(&dummy,sizeof(dummy),1,fdata);
	printf("dummy divergence                    = %.9d",dummy);
	printf(" - expected = %.9lu\n",sizeof(float)*Header.npartTotal[0]);
	for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].div,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      
      #ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      printf("------we've read \t %s\n",arr_name);
      #endif
      
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy curl                          = %.9d",dummy);
      printf(" - expected = %.9lu\n",sizeof(float)*Header.npartTotal[0]);
      for(i=0;i<Header.npartTotal[0];i++)
      returnRead = fread(&gaspro[i].curl,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);

      #endif
      */

        //============================
      //Read rate of cooling for gas 
      //============================
#ifdef OUTPUTCOOLRATE
         
#ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&arr_name,sizeof(char),4,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);	    
      printf("------we've read \t %s\n",arr_name);
#endif 
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy cooling rate                  = %.9d",dummy);
      printf(" - expected = %.9lu\n",sizeof(float)*Header.npartTotal[0]);
      for(i=0;i<Header.npartTotal[0];i++)
	{
	  returnRead = fread(&gaspro[i].coolr,sizeof(float),1,fdata);
#ifdef MAGNETIC_FIELD
	  fprintf(filel,"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
		  particles[i].pos[0], particles[i].pos[1], particles[i].pos[2],
		  gaspro[i].B[0], gaspro[i].B[1], gaspro[i].B[2], gaspro[i].sfr, gaspro[i].coolr, gaspro[i].rho, gaspro[i].U);
#else
	  fprintf(filel,"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
		  particles[i].pos[0], particles[i].pos[1], particles[i].pos[2],
		  gaspro[i].sfr, gaspro[i].coolr, gaspro[i].rho, gaspro[i].U);	  
#endif
	}
	  
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
     
#endif //OUTPUTCOOLRATE
      
#ifdef OUTPUT_VORTICITY
      
#ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      returnRead = fread(&arr_name,sizeof(char),4,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("------we've read \t %s\n",arr_name);
#endif
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy vorticity                     = %.9d",dummy);
      printf(" - expected = %.9lu\n",3*sizeof(float)*Header.npartTotal[0]);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].vorticity[0],sizeof(float),3,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif //OUTPUT_VORTICITY

#ifdef WINDS

#ifdef GADGET2_FORMAT
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      returnRead = fread(&arr_name,sizeof(char),4,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("------we've read \t %s\n",arr_name);
#endif
      
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      printf("dummy delay time                    = %.9d",dummy);
      printf(" - expected = %.9lu\n",sizeof(float)*Header.npartTotal[0]);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].delaytime,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif  //WINDS

    }
  printf("\n");

  fclose(fdata);
  fclose(filel);
  
  
  //######################################
  //Sorting particles by id
  //######################################

  particulas  *aux;
  gas_properties *auxgas;
  size_t  *p;
  
#ifdef LONGIDS
  unsigned long long *ID;
#else
  unsigned int *ID;
#endif
  
  for(type=0; type<6; type++)
    {
      
      indexmin = indexmax = 0;
      
      for(i=0;i<type;i++)
	{
	  indexmin = indexmin + N_part[i];
	}
      
      indexmax = indexmin + N_part[type];
      
      aux = (particulas *)malloc((size_t)N_part[type]*sizeof(particulas));
      if(particles == NULL)
	{
	  printf("Allocation of aux failed\n");
	  exit(0);
	}
      
      p = (size_t *)malloc((size_t)N_part[type]*sizeof(size_t));
      if(p == NULL)
      {
	printf("Allocation of p failed\n");
	exit(0);
      }
      
#ifdef LONGIDS
      ID = (unsigned long long *)malloc((size_t)N_part[type]*sizeof(unsigned long long));
      if(ID == NULL)
	{
	  printf("Allocation of ID failed\n");
	  exit(0);
	}
#else
      ID = (unsigned int *)malloc((size_t)N_part[type]*sizeof(unsigned int));
      if(ID == NULL)
	{
	  printf("Allocation of ID failed\n");
	  exit(0);
	}
#endif
      for(i=indexmin; i<indexmax; i++)
	{  
	  ID[i-indexmin] = particles[i].id;  
	  aux[i-indexmin] = particles[i];
	}
      gsl_sort_uint_index(p,ID,1,(size_t) N_part[type]); // OJO, con el tipo de ID, no hay funcion para unsigned long long
      
      printf("type %d sorted\n",type); 
      //gsl_sort_index(p,ID,1,(size_t)N_part[type]);
      
      for(i=indexmin; i<indexmax; i++)
	{
	  particles[i] = aux[p[i-indexmin]];   
	}
      
      if( (type == 0) && (N_part[0] > 0) )
	{
	  auxgas = (gas_properties *)malloc((size_t) N_part[0]*sizeof(gas_properties));
	  if(auxgas == NULL){
	    printf("Allocation of auxgas failed\n");
	    exit(0);
	  }
	  
	  for( i=indexmin; i<indexmax; i++)
	    {
	      auxgas[i-indexmin] = gaspro[i];
	    }
	  
	  for( i=indexmin; i<indexmax; i++)
	    {
	      gaspro[i] = auxgas[p[i-indexmin]];
	    }
	  
	  free(auxgas);
	}
      
      free(ID);
      free(aux);
      free(p);
      
    }
  
  n0 = Header.npartTotal[0];
  n1 = Header.npartTotal[1];
  n2 = Header.npartTotal[2];
  n3 = Header.npartTotal[3];
  n4 = Header.npartTotal[4];
  n5 = Header.npartTotal[5];
  
  for( i=0; i<N_part_total; i++)
    {
      if(i < n0) 
	particles[i].type = 0;
      
      if( (i >= n0) && ( i < (n0+n1)) ) 
	particles[i].type = 1;
      
      if( (i >= (n0+n1)) && (i < (n0+n1+n2)) ) 
	particles[i].type = 2;
     
      if( (i >= (n0+n1+n2) ) && (i < (n0+n1+n2+n3)) ) 
	particles[i].type = 3;
      
      if( (i >= (n0+n1+n2+n3)) && (i < (n0+n1+n2+n3+n4)) ) 
	particles[i].type = 4;
      
      if((i >= (n0+n1+n2+n3+n4)) && (i < (n0+n1+n2+n3+n4+n5)) ) 
	particles[i].type = 5;
    }
  
  //==================================
  // Computing total mass by component
  //==================================
 
double totalMasses[7];

 N_min = 0;
 totalMasses[6] = 0.0;
 
for( type=0; type<6; type++)
  {

    N_max = N_min + N_part[type];
    totalMasses[type] = 0.0;

    for(i=N_min;i<N_max;i++)
      {	
	totalMasses[type] = totalMasses[type] + particles[i].mass;  
      }

    N_min = N_max;

    printf("Total mass type %d = %.8e\n",type,totalMasses[type]);
    totalMasses[6] = totalMasses[6] + totalMasses[type];

  }

 printf("Galaxy total mass = %.8e\n",totalMasses[6]);  

 // TAKE CARE!!
 // To free gaspro and particles in the code where this routine is called.  
 
 //free(gaspro);
 //free(partices);
 
 return 0;
}

int print_header(char filename[])
{

////////////////////////////////////////////////////////////////////////
 //SAVE HEADER
////////////////////////////////////////////////////////////////////////

  char tmp[500];
  FILE *fHeader;
  sprintf(tmp,"%s_Header.dat",filename);
  fHeader=fopen(tmp,"w");
  
  fprintf(fHeader,"\n");
  fprintf(fHeader,"Read header from %s\n",filename);
  fprintf(fHeader,"\n");
  fprintf(fHeader,"int Npart[6] : \n");
  fprintf(fHeader,"Npart[0] = %d \n",Header.Npart[0]);
  fprintf(fHeader,"Npart[1] = %d \n",Header.Npart[1]);
  fprintf(fHeader,"Npart[2] = %d \n",Header.Npart[2]);
  fprintf(fHeader,"Npart[3] = %d \n",Header.Npart[3]);
  fprintf(fHeader,"Npart[4] = %d \n",Header.Npart[4]);
  fprintf(fHeader,"Npart[5] = %d \n",Header.Npart[5]);
  fprintf(fHeader,"double mass[6] : \n");
  fprintf(fHeader,"mass[0] = %.10lf \n",Header.mass[0]);
  fprintf(fHeader,"mass[1] = %.10lf \n",Header.mass[1]);
  fprintf(fHeader,"mass[2] = %.10lf \n",Header.mass[2]);
  fprintf(fHeader,"mass[3] = %.10lf \n",Header.mass[3]);
  fprintf(fHeader,"mass[4] = %.10lf \n",Header.mass[4]);
  fprintf(fHeader,"mass[5] = %.10lf \n",Header.mass[5]);
  fprintf(fHeader,"double time = %lf \n",Header.time);
  fprintf(fHeader,"double redshift = %lf \n",Header.redshift);
  fprintf(fHeader,"int flag_sfr = %d \n",Header.flag_sfr);
  fprintf(fHeader,"int flag_feedback = %d \n",Header.flag_feedback);
  fprintf(fHeader,"int npartTotal[6] : \n");
  fprintf(fHeader,"npartTotal[0] = %u \n",Header.npartTotal[0]);
  fprintf(fHeader,"npartTotal[1] = %u \n",Header.npartTotal[1]);
  fprintf(fHeader,"npartTotal[2] = %u \n",Header.npartTotal[2]);
  fprintf(fHeader,"npartTotal[3] = %u \n",Header.npartTotal[3]);
  fprintf(fHeader,"npartTotal[4] = %u \n",Header.npartTotal[4]);
  fprintf(fHeader,"npartTotal[5] = %u \n",Header.npartTotal[5]);
  fprintf(fHeader,"int flag_cooling = %d \n",Header.flag_cooling);
  fprintf(fHeader,"int num_files = %d \n",Header.num_files); 
  fprintf(fHeader,"double BoxSize = %lf \n",Header.BoxSize);
  fprintf(fHeader,"double Omega0 = %lf \n",Header.Omega0);
  fprintf(fHeader,"double OmegaLambda = %lf \n",Header.OmegaLambda);
  fprintf(fHeader,"double HubbleParam = %lf \n",Header.HubbleParam);
  fprintf(fHeader,"int flag_stellarage = %d \n",Header.flag_stellarage);
  fprintf(fHeader,"int flag_metals = %d \n",Header.flag_metals);
  fprintf(fHeader,"int npartTotalHighWord[6] : \n");
  fprintf(fHeader,"npartTotalHighWord[0] = %u \n",Header.npartTotalHighWord[0]);
  fprintf(fHeader,"npartTotalHighWord[1] = %u \n",Header.npartTotalHighWord[1]);
  fprintf(fHeader,"npartTotalHighWord[2] = %u \n",Header.npartTotalHighWord[2]);
  fprintf(fHeader,"npartTotalHighWord[3] = %u \n",Header.npartTotalHighWord[3]);
  fprintf(fHeader,"npartTotalHighWord[4] = %u \n",Header.npartTotalHighWord[4]);
  fprintf(fHeader,"npartTotalHighWord[5] = %u \n",Header.npartTotalHighWord[5]);
  fprintf(fHeader,"int flag_entropy_instead_u = %d \n",Header.flag_entropy_instead_u);
  fprintf(fHeader,"Rest size to 256 Bytes = %lu",sizeof(Header.fill));
  
  fclose(fHeader);

  return 0;
}

int print_data_ascii(char filename[],int component)
{

  int i, indexmin, indexmax;
  char outfile[500];
  /*
  for(i=0; i<component; i++)
    {
      indexmin = indexmax = indexmin + N_part[component];
    }

    indexmax = indexmax + N_part[component];*/
  
  //======================================================
  //para ir al bloque de las particulas de esa componente
  //y leer todo el bloque (si es aislada la lee toda,
  //si es un merger, lee host+sat).
  //======================================================

  if(flag==0)
    {
      indexmin = indexmax= 0;
      for( i=0; i<component; i++)
	{
	  indexmin = indexmin + N_part[i];
	  printf("cuando component==0 no entramos aquí\n");
	}
      indexmax = indexmin + N_part[component];

      printf("ALL BLOCK: \t INDEXMIN = %d \t INDEXMAX = %d \t\n",indexmin,indexmax);
    }
  
  //======================================================
  //para ir al bloque de las particulas de esa componente
  //y leer solo la primera parte del bloque,
  //es decir la primera galaxia -host
  //(solo aplica para cuando hay merger).
  //======================================================

  if(flag==1)
    {
      indexmin = indexmax= 0;
      for( i=0; i<component; i++)
	{
	  //indexmin = indexmax = indexmin + N_part[i];
	  indexmin = indexmin + N_part[i];
	}
      indexmax = indexmin + Nhost[component];
      printf("--GALAXY 1 (host) component %d \t INDEXMIN = %d \t INDEXMAX = %d \t\n",component,indexmin,indexmax);
    }

  //======================================================
  //para ir al bloque de las particulas de esa componente
  //y leer solo la segunda parte del bloque,
  //es decir la segunda galaxia -satellite
  //(solo aplica para cuando hay merger).
  //======================================================
  

  if(flag==2)
    {
      indexmin = indexmax= 0;
      for( i=0; i<component; i++)
	{
	  //indexmin = indexmax = indexmin + N_part[i];
	  indexmin = indexmin + N_part[i];
	}
      indexmin = indexmin + Nhost[component];
      indexmax = indexmin + Nsat[component];
      
      printf("GALAXY 2 (satellite) component %d \t INDEXMIN = %d \t INDEXMAX = %d \t\n",component,indexmin,indexmax);
    }


  if( N_part[component] > 0 )
    {
      sprintf(outfile,"%s.%d",filename,component);
      FILE *outfiles;
      outfiles = fopen(outfile,"w");
      if(outfiles==NULL) printf("Could not open file %s\n",outfile);  
      
      for(i=indexmin;i<indexmax;i++)
	{
	      
#ifdef LONGIDS
	  fprintf(outfiles,"%lu",particles[i].id);
#else
	      fprintf(outfiles,"%u",particles[i].id);
#endif
	      
	      fprintf(outfiles," %f %f %f %f %f %f %.8f",
		      particles[i].pos[0],particles[i].pos[1],particles[i].pos[2],
		      particles[i].vel[0],particles[i].vel[1],particles[i].vel[2],
		      particles[i].mass);
	      
#ifdef PRINT_ALL_ASCII
	      
	      if( component == 0 &&  N_part[0] >0 )
		{
		  fprintf(outfiles," %f %f",gaspro[i].U,gaspro[i].rho);
#ifdef COOLING
		  fprintf(outfiles," %f %f",gaspro[i].Ne,gaspro[i].Nh);
#endif
		  fprintf(outfiles," %f",gaspro[i].h);
#ifdef SFR
		  fprintf(outfiles," %f",gaspro[i].sfr);
#endif
		}
	      
#ifdef STELLARAGE	      
	      if( component == 4 &&  N_part[4] >0 )
		fprintf(outfiles," %f",particles[i].stellar_age);
#endif
	      
#ifdef METALS
	      if( ((component == 0) &&  (N_part[0]>0)) || ((component == 4) &&  (N_part[4]>0)) )
		fprintf(outfiles," %f",particles[i].metallicity);
#endif
	      
#endif // close PRINT_ALL_ASCII
	      
#ifdef OUTPUTPOTENTIAL
	      fprintf(outfiles," %f",particles[i].pot);
#endif
	      
#ifdef PRINT_ALL_ASCII
	      
#ifdef OUTPUTACCELERATION
	      fprintf(outfiles," %f %f %f",particles[i].acce[X],particles[i].acce[Y],particles[i].acce[Z]);
#endif
	      
#ifdef OUTPUTCHANGEOFENTROP
	      if( (component == 0) && (N_part[0] >0) )
		fprintf(outfiles," %f",gaspro[i].ecr);
#endif 
	      
#ifdef OUTPUTTIMESTEP
	      fprintf(outfiles," %f",particles[i].timestep);
#endif
	   
#ifdef OUTPUTCOOLRATE
	      if( (component == 0) && (N_part[0] >0) )
		fprintf(outfiles," %f",gaspro[i].coolr);
#endif
	      
#ifdef OUTPUT_DIV_CURL
	      if( (component == 0) && (N_part[0] >0) )
		fprintf(outfiles," %f %f",gaspro[i].div,gaspro[i].curl);
#endif

#ifdef OUTPUT_VORTICITY
	      fprintf(outfiles," %f %f %f",gaspro[i].vorticity[X],gaspro[i].vorticity[Y],gaspro[i].vorticity[Z]);
#endif
	      
#ifdef WINDS
	      fprintf(outfiles," %f",gaspro[i].delaytime);
#endif
	      
#endif //close PRINT_ALL_ASCII
	      
	      fprintf(outfiles,"\n");	
	      
	}
      
      fclose(outfiles);
    }
  

printf("\nDone.\n\n");

return 0;
}


int print_data_ascii_all_types(char filename[])
{

  int i, type;
  char outfile[500];
  
  N_min = 0;
  
  for( type=0; type<6; type++)
    {
      
      N_max = N_min + N_part[type];  
      
      if( N_part[type] > 0 )
	{
	  sprintf(outfile,"data_merger.%d",type);
	  FILE *outfiles;
	  outfiles = fopen(outfile,"w");
	  if(outfiles==NULL) printf("Could not open file %s\n",outfile);  
	  
	  for(i=N_min;i<N_max;i++)
	    {
	      
#ifdef LONGIDS
	      fprintf(outfiles,"%lu",particles[i].id);
#else
	      fprintf(outfiles,"%u",particles[i].id);
#endif
	      
	      fprintf(outfiles," %f %f %f %f %f %f %.8f",
		      particles[i].pos[0],particles[i].pos[1],particles[i].pos[2],
		      particles[i].vel[0],particles[i].vel[1],particles[i].vel[2],
		      particles[i].mass);

#ifdef PRINT_ALL_ASCII
	      
	      if( type == 0 &&  N_part[0] >0 )
		{
		  fprintf(outfiles," %f %f",gaspro[i].U,gaspro[i].rho);
#ifdef COOLING
		  fprintf(outfiles," %f %f",gaspro[i].Ne,gaspro[i].Nh);
#endif
		  fprintf(outfiles," %f",gaspro[i].h);
#ifdef SFR
		  fprintf(outfiles," %f",gaspro[i].sfr);
#endif
		}

#ifdef STELLARAGE	      
	      if( type == 4 &&  N_part[4] >0 )
		fprintf(outfiles," %f",particles[i].stellar_age);
#endif

#ifdef METALS
	      if( ((type == 0) &&  (N_part[0]>0)) || ((type == 4) &&  (N_part[4]>0)) )
		fprintf(outfiles," %f",particles[i].metallicity);
#endif

#endif // close PRINT_ALL_ASCII
	      
#ifdef OUTPUTPOTENTIAL
	      fprintf(outfiles," %f",particles[i].pot);
#endif

#ifdef PRINT_ALL_ASCII
	      
#ifdef OUTPUTACCELERATION
	      fprintf(outfiles," %f %f %f",particles[i].acce[X],particles[i].acce[Y],particles[i].acce[Z]);
#endif
	      
#ifdef OUTPUTCHANGEOFENTROP
	      if( (type == 0) && (N_part[0] >0) )
		fprintf(outfiles," %f",gaspro[i].ecr);
#endif 
	      
#ifdef OUTPUTTIMESTEP
	      fprintf(outfiles," %f",particles[i].timestep);
#endif
	   
#ifdef OUTPUTCOOLRATE
	      if( (type == 0) && (N_part[0] >0) )
		fprintf(outfiles," %f",gaspro[i].coolr);
#endif
   
#ifdef OUTPUT_DIV_CURL
	      if( (type == 0) && (N_part[0] >0) )
		fprintf(outfiles," %f %f",gaspro[i].div,gaspro[i].curl);
#endif

#ifdef OUTPUT_VORTICITY
	      fprintf(outfiles," %f %f %f",gaspro[i].vorticity[X],gaspro[i].vorticity[Y],gaspro[i].vorticity[Z]);
#endif

#ifdef WINDS
	      fprintf(outfiles," %f",gaspro[i].delaytime);
#endif

#endif //close PRINT_ALL_ASCII

	      fprintf(outfiles,"\n");	
	      
	    }
	  fclose(outfiles);
	}
      
      N_min = N_max;
      
    }
  
  printf("\nDone.\n\n");

  return 0;
}

int write_gadget1(char *filename, char outfile2[])
{
  
  int i,type;
  FILE *fParam,*fHeader,*fGadget;
  char outfile[500];
  int dummy;
    
  
  ////////////////////////////////////////////////////////////////////////
  // READ DATA FILE
  ////////////////////////////////////////////////////////////////////////
  //filename = argv[1];
  fParam = fileOpen(filename,"r");
  
  ////////////////////////////////////////////////////////////////////////
  // NAME OUTPUT FILE
  ////////////////////////////////////////////////////////////////////////
  //returnRead = fscanf(fParam,"%s",outfile);
  //fData = fileOpen(outfile,"r");
  //printf("\nBuilding Gadget binary with format 1 in:\n%s\n\n",outfile);

  //sprintf(outfile2,"%s.gad",outfile);
  fGadget = fopen(outfile2,"w");

  
  ////////////////////////////////////////////////////////////////////////
  // READ HEADER
  ////////////////////////////////////////////////////////////////////////
  returnRead = fscanf(fParam,"%d %d %d %d %d %d",
	 &Header.Npart[0],&Header.Npart[1],&Header.Npart[2],
	 &Header.Npart[3],&Header.Npart[4],&Header.Npart[5]);
  returnRead = fscanf(fParam,"%lf %lf %lf %lf %lf %lf",
	 &Header.mass[0],&Header.mass[1],&Header.mass[2],
	 &Header.mass[3],&Header.mass[4],&Header.mass[5]);
  returnRead = fscanf(fParam,"%lf",&Header.time);
  returnRead = fscanf(fParam,"%lf",&Header.redshift);
  returnRead = fscanf(fParam,"%d",&Header.flag_sfr);
  returnRead = fscanf(fParam,"%d",&Header.flag_feedback);
  returnRead = fscanf(fParam,"%u %u %u %u %u %u",
	 &Header.npartTotal[0],&Header.npartTotal[1],&Header.npartTotal[2],
	 &Header.npartTotal[3],&Header.npartTotal[4],&Header.npartTotal[5]);
  returnRead = fscanf(fParam,"%d",&Header.flag_cooling);
  returnRead = fscanf(fParam,"%d",&Header.num_files);
  returnRead = fscanf(fParam,"%lf",&Header.BoxSize);
  returnRead = fscanf(fParam,"%lf",&Header.Omega0);
  returnRead = fscanf(fParam,"%lf",&Header.OmegaLambda);
  returnRead = fscanf(fParam,"%lf",&Header.HubbleParam);
  returnRead = fscanf(fParam,"%d",&Header.flag_stellarage);
  returnRead = fscanf(fParam,"%d",&Header.flag_metals);
  returnRead = fscanf(fParam,"%u %u %u %u %u %u",
	 &Header.npartTotalHighWord[0],&Header.npartTotalHighWord[1],
	 &Header.npartTotalHighWord[2],&Header.npartTotalHighWord[3],
	 &Header.npartTotalHighWord[4],&Header.npartTotalHighWord[5]);
 returnRead = fscanf(fParam,"%d",&Header.flag_entropy_instead_u);

 fclose(fParam);
  
 ////////////////////////////////////////////////////////////////////////
 // SAVE HEADER
 ////////////////////////////////////////////////////////////////////////
 sprintf(outfile,"%s_Header_wrote.output",outfile2);
 fHeader=fopen(outfile,"w");
 
 fprintf(fHeader,"\n");
 fprintf(fHeader,"Wrote Header for %s\n",outfile2);
 fprintf(fHeader,"\n");
 fprintf(fHeader,"int Npart[6] : \n");
 fprintf(fHeader,"Npart[0] = %d \n",Header.Npart[0]);
 fprintf(fHeader,"Npart[1] = %d \n",Header.Npart[1]);
 fprintf(fHeader,"Npart[2] = %d \n",Header.Npart[2]);
 fprintf(fHeader,"Npart[3] = %d \n",Header.Npart[3]);
 fprintf(fHeader,"Npart[4] = %d \n",Header.Npart[4]);
 fprintf(fHeader,"Npart[5] = %d \n",Header.Npart[5]);
 fprintf(fHeader,"double mass[6] : \n");
 fprintf(fHeader,"mass[0] = %.10e \n",Header.mass[0]);
 fprintf(fHeader,"mass[1] = %.10e \n",Header.mass[1]);
 fprintf(fHeader,"mass[2] = %.10e \n",Header.mass[2]);
 fprintf(fHeader,"mass[3] = %.10e \n",Header.mass[3]);
 fprintf(fHeader,"mass[4] = %.10e \n",Header.mass[4]);
 fprintf(fHeader,"mass[5] = %.10e \n",Header.mass[5]);
 fprintf(fHeader,"double time = %lf \n",Header.time);
 fprintf(fHeader,"double redshift = %lf \n",Header.redshift);
 fprintf(fHeader,"int flag_sfr = %d \n",Header.flag_sfr);
 fprintf(fHeader,"int flag_feedback = %d \n",Header.flag_feedback);
 fprintf(fHeader,"int npartTotal[6] : \n");
 fprintf(fHeader,"npartTotal[0] = %u \n",Header.npartTotal[0]);
 fprintf(fHeader,"npartTotal[1] = %u \n",Header.npartTotal[1]);
 fprintf(fHeader,"npartTotal[2] = %u \n",Header.npartTotal[2]);
 fprintf(fHeader,"npartTotal[3] = %u \n",Header.npartTotal[3]);
 fprintf(fHeader,"npartTotal[4] = %u \n",Header.npartTotal[4]);
 fprintf(fHeader,"npartTotal[5] = %u \n",Header.npartTotal[5]);
 fprintf(fHeader,"int flag_cooling = %d \n",Header.flag_cooling);
 fprintf(fHeader,"int num_files = %d \n",Header.num_files); 
 fprintf(fHeader,"double BoxSize = %lf \n",Header.BoxSize);
 fprintf(fHeader,"double Omega0 = %lf \n",Header.Omega0);
 fprintf(fHeader,"double OmegaLambda = %lf \n",Header.OmegaLambda);
 fprintf(fHeader,"double HubbleParam = %lf \n",Header.HubbleParam);
 fprintf(fHeader,"int flag_stellarage = %d \n",Header.flag_stellarage);
 fprintf(fHeader,"int flag_metals = %d \n",Header.flag_metals);
 fprintf(fHeader,"int npartTotalHighWord[6] : \n");
 fprintf(fHeader,"npartTotalHighWord[0] = %u \n",Header.npartTotalHighWord[0]);
 fprintf(fHeader,"npartTotalHighWord[1] = %u \n",Header.npartTotalHighWord[1]);
 fprintf(fHeader,"npartTotalHighWord[2] = %u \n",Header.npartTotalHighWord[2]);
 fprintf(fHeader,"npartTotalHighWord[3] = %u \n",Header.npartTotalHighWord[3]);
 fprintf(fHeader,"npartTotalHighWord[4] = %u \n",Header.npartTotalHighWord[4]);
 fprintf(fHeader,"npartTotalHighWord[5] = %u \n",Header.npartTotalHighWord[5]);
 fprintf(fHeader,"int flag_entropy_instead_u = %d \n",Header.flag_entropy_instead_u);
 fprintf(fHeader,"Rest size to 256 Bytes = %lu",sizeof(Header.fill));

 fclose(fHeader);   	       
 
 ////////////////////////////////////////////////////////////////////////
 // TOTAL NUMBER OF PARTICLES
 ////////////////////////////////////////////////////////////////////////
 N_part_total = 0;
 for(i=0; i<6; i++)
   {
     N_part[i] = Header.npartTotal[i];
     N_part_total += N_part[i];
   }
 
 ////////////////////////////////////////////////////////////////////////
 // ALLOCATE AND READ
 ////////////////////////////////////////////////////////////////////////
 /*
 particles = (particulas *)malloc((size_t)N_part_total*sizeof(particulas));
 if(particles == NULL){
   printf("Allocation of particles failed\n");
   exit(0);
 }
 
 U = (float *)malloc((size_t)Header.npartTotal[0]*sizeof(float));
 if(particles == NULL){
   printf("Allocation of U failed\n");
   exit(0);
 }
 */

 /*
 ////////////////////////////////////////////////////////////
 // READING PARTICLES DATA
 ///////////////////////////////////////////////////////////
 
 //////////////////////////////////
 // READING GAS
 //////////////////////////////////
 if( Header.npartTotal[0] != 0)
   {
     for( i=0; i<Header.npartTotal[0]; i++)
       {
#ifdef LONGIDS
	 returnRead = fscanf(fData,"%llu %f %f %f %f %f %f %f %f",
		&particles[i].id,
		&particles[i].pos[X],&particles[i].pos[Y],&particles[i].pos[Z],
		&particles[i].vel[X],&particles[i].vel[Y],&particles[i].vel[Z],
		&particles[i].mass,
		&U[i]);
#else
	 returnRead = fscanf(fData,"%u %f %f %f %f %f %f %f %f",
		&particles[i].id,
		&particles[i].pos[X],&particles[i].pos[Y],&particles[i].pos[Z],
		&particles[i].vel[X],&particles[i].vel[Y],&particles[i].vel[Z],
		&particles[i].mass,
		&U[i]);
#endif
       }
     
     imin = Header.npartTotal[0];
   }
 
////////////////////////////////////////////////////////
// READING OTHER PARTICLES TYPE
////////////////////////////////////////////////////////
 for( i=imin; i<N_part_total; i++)
   {
#ifdef LONGIDS
     returnRead = fscanf(fData,"%llu %f %f %f %f %f %f %f",
	    &particles[i].id,
	    &particles[i].pos[X],&particles[i].pos[Y],&particles[i].pos[Z],
	    &particles[i].vel[X],&particles[i].vel[Y],&particles[i].vel[Z],
	    &particles[i].mass);
#else
     returnRead = fscanf(fData,"%u %f %f %f %f %f %f %f",
	    &particles[i].id,
	    &particles[i].pos[X],&particles[i].pos[Y],&particles[i].pos[Z],
	    &particles[i].vel[X],&particles[i].vel[Y],&particles[i].vel[Z],
	    &particles[i].mass);
#endif
   }
 
 fclose(fData);
*/
 //////////////////////////////////
 // WRITING HEADER
 /////////////////////////////////
 dummy = sizeof(Header);
 fwrite(&dummy,sizeof(dummy),1,fGadget);
 fwrite(&Header,sizeof(Header),1,fGadget);
 fwrite(&dummy,sizeof(dummy),1,fGadget);

 ////////////////////////////////////////////////////////////////////////
 // WRITING POSITIONS
 ////////////////////////////////////////////////////////////////////////
 dummy = 3*N_part_total*sizeof(float);
 fwrite(&dummy,sizeof(dummy),1,fGadget); 
 for( i=0; i<N_part_total; i++ ) 
   fwrite(&particles[i].pos,sizeof(float),3,fGadget);
 fwrite(&dummy,sizeof(dummy),1,fGadget);

 ////////////////////////////////////////////////////////////////////////
 // WRITING VELOCITIES
 ////////////////////////////////////////////////////////////////////////
 dummy = 3*N_part_total*sizeof(float);
 fwrite(&dummy,sizeof(dummy),1,fGadget); 
 for( i=0; i<N_part_total; i++ ) 
   fwrite(&particles[i].vel,sizeof(float),3,fGadget);
 fwrite(&dummy,sizeof(dummy),1,fGadget);

 ////////////////////////////////////////////////////////////////////////
 // WRITING IDs
 ////////////////////////////////////////////////////////////////////////
#ifdef LONGIDS
 dummy = N_part_total*sizeof(unsigned long long);	
#else
 dummy = N_part_total*sizeof(unsigned int);		
#endif

 fwrite(&dummy,sizeof(dummy),1,fGadget);
#ifdef LONGIDS
 for(i=0; i<N_part_total; i++)
   fwrite(&particles[i].Id,sizeof(unsigned long long),1,fGadget);
#else
 for(i=0; i<N_part_total; i++)
   fwrite(&particles[i].id,sizeof(unsigned int),1,fGadget);
#endif
 fwrite(&dummy,sizeof(dummy),1,fGadget);
 
 ////////////////////////////////////////////////
 // WRITING MASSES
 ///////////////////////////////////////////////
 dummy = 0; 
 for( type=0; type<6; type++)
   if( Header.npartTotal[type]>0 )
     {
       if( Header.mass[type]>0.0 )
	 continue;
       else
	 dummy = dummy + Header.npartTotal[type]*sizeof(float);
     }
 if(dummy>0)
   {
     fwrite(&dummy,sizeof(dummy),1,fGadget);
     
     N_min = N_max = 0;
     for( type=0; type<6; type++)
       {
	 N_max = N_max + Header.npartTotal[type];
	 if( (Header.npartTotal[type]>0) && (Header.mass[type]>0.0) )
	   continue;
	 else
	   {
	     for( i=N_min; i<N_max; i++)
	       {
		 fwrite(&particles[i].mass,sizeof(float),1,fGadget);
	       }
	     N_min = N_max;
	   }
       } 
     fwrite(&dummy,sizeof(dummy),1,fGadget);
   }

 ////////////////////////////////////////////////////////////
 // WRITING INTERNAL ENERGY FOR GAS
 ///////////////////////////////////////////////////////////
 if( Header.npartTotal[0]>0 )
   {
     dummy = Header.npartTotal[0]*sizeof(float);
     fwrite(&dummy,sizeof(dummy),1,fGadget); 
     for( i=0; i<Header.npartTotal[0]; i++ )
      	 fwrite(&U[i],sizeof(float),1,fGadget);
     fwrite(&dummy,sizeof(dummy),1,fGadget);
   }

 printf("Initial conditions in %s\n\n",outfile2);

 fclose(fGadget);   	      
  
 //free(particles);
 //free(U);
 
 return 0;
}

int write_gadget_ready_data_header(char *outfile)
{
  
  int i,type;
  FILE *fGadget;
  int dummy;
    
  fGadget = fopen(outfile,"w");

 //////////////////////////////////
 // WRITING HEADER
 /////////////////////////////////
 dummy = sizeof(Header);
 fwrite(&dummy,sizeof(dummy),1,fGadget);
 fwrite(&Header,sizeof(Header),1,fGadget);
 fwrite(&dummy,sizeof(dummy),1,fGadget);

 ////////////////////////////////////////////////////////////////////////
 // WRITING POSITIONS
 ////////////////////////////////////////////////////////////////////////
 dummy = 3*N_part_total*sizeof(float);
 fwrite(&dummy,sizeof(dummy),1,fGadget); 
 for( i=0; i<N_part_total; i++ ) 
   fwrite(&particles[i].pos,sizeof(float),3,fGadget);
 fwrite(&dummy,sizeof(dummy),1,fGadget);

 ////////////////////////////////////////////////////////////////////////
 // WRITING VELOCITIES
 ////////////////////////////////////////////////////////////////////////
 dummy = 3*N_part_total*sizeof(float);
 fwrite(&dummy,sizeof(dummy),1,fGadget); 
 for( i=0; i<N_part_total; i++ ) 
   fwrite(&particles[i].vel,sizeof(float),3,fGadget);
 fwrite(&dummy,sizeof(dummy),1,fGadget);

 ////////////////////////////////////////////////////////////////////////
 // WRITING IDs
 ////////////////////////////////////////////////////////////////////////
#ifdef LONGIDS
 dummy = N_part_total*sizeof(unsigned long long);	
#else
 dummy = N_part_total*sizeof(unsigned int);		
#endif

 fwrite(&dummy,sizeof(dummy),1,fGadget);
#ifdef LONGIDS
 for(i=0; i<N_part_total; i++)
   fwrite(&particles[i].Id,sizeof(unsigned long long),1,fGadget);
#else
 for(i=0; i<N_part_total; i++)
   fwrite(&particles[i].id,sizeof(unsigned int),1,fGadget);
#endif
 fwrite(&dummy,sizeof(dummy),1,fGadget);
 
 ////////////////////////////////////////////////
 // WRITING MASSES
 ///////////////////////////////////////////////
 dummy = 0; 
 for( type=0; type<6; type++)
   if( Header.npartTotal[type]>0 )
     {
       if( Header.mass[type]>0.0 )
	 continue;
       else
	 dummy = dummy + Header.npartTotal[type]*sizeof(float);
     }
 if(dummy>0)
   {
     fwrite(&dummy,sizeof(dummy),1,fGadget);
     
     N_min = N_max = 0;
     for( type=0; type<6; type++)
       {
	 N_max = N_max + Header.npartTotal[type];
	 if( (Header.npartTotal[type]>0) && (Header.mass[type]>0.0) )
	   continue;
	 else
	   {
	     for( i=N_min; i<N_max; i++)
	       {
		 fwrite(&particles[i].mass,sizeof(float),1,fGadget);
	       }
	     N_min = N_max;
	   }
       } 
     fwrite(&dummy,sizeof(dummy),1,fGadget);
   }

 ////////////////////////////////////////////////////////////
 // WRITING INTERNAL ENERGY FOR GAS
 ///////////////////////////////////////////////////////////
 if( Header.npartTotal[0]>0 )
   {
     dummy = Header.npartTotal[0]*sizeof(float);
     fwrite(&dummy,sizeof(dummy),1,fGadget); 
     for( i=0; i<Header.npartTotal[0]; i++ )
      	 fwrite(&U[i],sizeof(float),1,fGadget);
     fwrite(&dummy,sizeof(dummy),1,fGadget);
   }

 printf("Initial conditions in %s\n\n",outfile);

 fclose(fGadget);   	      
  
 // TAKE CARE!!
 // To free U and particles in the code where this routine is called.  

 //free(particles);
 //free(U);
 
 return 0;
}
