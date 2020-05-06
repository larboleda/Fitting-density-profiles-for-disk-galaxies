int split_merger()
{
  
  int i, j, k, count, countf, p, type;
  int count_star, count_gas, N_galaxy0;
  int indexmin, indexmax;
  int Nguess;
  double Rsfr;
  double f;

  f = r_scalength[2]; //factor to define aperture for new star selection
  
  char buf[200];
  //FILE *pos=NULL;

  typedef struct
  {
    unsigned int id;
    float pos[3];
    float vel[3];
    float mass;
    int type;
    float  rho;
    float Ne;
    float Nh;
    float U;
    float sfr;
    float coolr;
    float B[3];
    float metallicity;
    float pot;
    
  } aux;

  aux *galaxy_aux;

  //===============
  // HOST
  //===============
  
  if(flag==1) //host
    {
      Nguess = N_part[0] + N_part[4];
      Rsfr = f*r_scalength[0];
      N_galaxy0 = Nhost[0] + Nhost[1] + Nhost[2] + Nhost[3] + Nhost[4] + Nhost[5];
      count_gas = 0;
      count_star = 0;
      count = countf = 0;
      

      printf("\n (***) Nhost[0]=%d \t Nhost[1]=%d \t Nhost[2]=%d \t Nhost[3]=%d \t Nhost[4]=%d\n",Nhost[0],Nhost[1],Nhost[2],Nhost[3],Nhost[4]);

      galaxy = (galaxies *) malloc((size_t)N_galaxy0*sizeof(galaxies));
      
      if(galaxy == NULL)
	{
	  printf("Allocation of galaxy failed\n");
	  exit(0);
	}
      
      galaxy_aux = (aux *)malloc((size_t)N_galaxy0*sizeof(aux));
      
      if(galaxy_aux == NULL)
	{
	  printf("Allocation of galaxy failed\n");
	  exit(0);
	}
      
      //selecting gas particles 
      if(N_part[0]>0)
	{
	  galaxy_gaspro = (galaxy_gas_properties *)malloc((size_t) Nhost[0]*sizeof(galaxy_gas_properties));
	  
	  if(galaxy_gaspro == NULL)
	    {
	      printf("Allocation of galaxy_gaspro failed\n");
	      exit(0);
	    }
	  
	  ////sprintf(buf,"0type_host.dat");
	  
	  //pos = fopen(buf,"w");
	  
	  for(i=0;i<N_part[0];i++)
	    {
	      
	      if(particles[i].id < Nhost[0])
		{
		  galaxy_aux[count_gas].id     = particles[i].id;
		  galaxy_aux[count_gas].pos[0] = particles[i].pos[0];
		  galaxy_aux[count_gas].pos[1] = particles[i].pos[1];
		  galaxy_aux[count_gas].pos[2] = particles[i].pos[2];
		  galaxy_aux[count_gas].vel[0] = particles[i].vel[0];
		  galaxy_aux[count_gas].vel[1] = particles[i].vel[1];
		  galaxy_aux[count_gas].vel[2] = particles[i].vel[2];
		  galaxy_aux[count_gas].mass   = particles[i].mass;
		  galaxy_aux[count_gas].pot    = particles[i].pot;
		  
		  galaxy_aux[count_gas].rho = gaspro[i].rho;
		  galaxy_aux[count_gas].Ne  = gaspro[i].Ne;
		  galaxy_aux[count_gas].Nh  = gaspro[i].Nh;
		  galaxy_aux[count_gas].sfr = gaspro[i].sfr;
		  galaxy_aux[count_gas].coolr = gaspro[i].coolr;
		  galaxy_aux[count_gas].U = gaspro[i].U;
		  galaxy_aux[count_gas].rho = gaspro[i].rho;
#ifdef MAGNETIC_FIELD
		  galaxy_aux[count_gas].B[0] = gaspro[i].B[0];
		  galaxy_aux[count_gas].B[1] = gaspro[i].B[1];
		  galaxy_aux[count_gas].B[2] = gaspro[i].B[2];
#endif	  
		  count_gas += 1;
		}
	    }
	
	  Nhost[0] = count_gas; //actualizando número real de part tipo 0 en host.
	  p=0;
	  
	  double tot_sfr=0.0;
	  FILE *sfr_file = NULL;
	  sfr_file = fopen("sf_rate_host.dat","a");
	  for(k=0; k<count_gas; k++)
	    {
	      galaxy[k].pos[0] = galaxy_aux[k].pos[0];
	      galaxy[k].pos[1] = galaxy_aux[k].pos[1];
	      galaxy[k].pos[2] = galaxy_aux[k].pos[2];
	      galaxy[k].mass   = galaxy_aux[k].mass;
	      galaxy[k].pot    = galaxy_aux[k].pot;
	      
	      galaxy_gaspro[p].rho = galaxy_aux[k].rho;
	      galaxy_gaspro[p].Ne  = galaxy_aux[k].Ne;
	      galaxy_gaspro[p].Nh  = galaxy_aux[k].Nh;
	      galaxy_gaspro[p].sfr = galaxy_aux[k].sfr;
	      galaxy_gaspro[p].coolr = galaxy_aux[k].coolr;
	      galaxy_gaspro[p].U   = galaxy_aux[k].U;
	      galaxy_gaspro[p].rho = galaxy_aux[k].rho;

#ifdef MAGNETIC_FIELD
	      galaxy_gaspro[p].B[0]   = galaxy_aux[k].B[0];
	      galaxy_gaspro[p].B[1]   = galaxy_aux[k].B[1];
	      galaxy_gaspro[p].B[2]   = galaxy_aux[k].B[2];
#endif
             
	      tot_sfr += galaxy_aux[k].sfr;
	      
	      p += 1;
	
	    }
	  
	  printf("\ntype 0 particle's pos&vels for host...DONE\n");
	  // fclose(pos);

	  fprintf(sfr_file,"%16.8lf %16.8lf\n",tot_sfr,Header.time);
	  
	  //realloc;
	  //galaxy_aux = (aux*)realloc(galaxy_aux, count_gas*sizeof(aux));
	  free(galaxy_aux);
	}
      
      
      count = count_gas;
      countf = count + Nhost[1];
      
      for(type=1; type<4; type++) //halo, disk, bulge
	{	  
	 	  
	  //sprintf(buf,"%dtype_host.dat",type);
	  
	  //pos = fopen(buf,"w");
	  
	  p = 0;

	  indexmin = indexmax = 0;
	  
	  for( i=0; i<type; i++)
	    {
	      
	      indexmin = indexmin + N_part[i];
	    }
	  
	  //indexmax = indexmin + Nhost[type]; //ya no lo uso	  
	  
	  for( k=count; k<countf; k++ )
	    {
	      j = indexmin + p;
	      galaxy[k].pos[0] = particles[j].pos[0];
	      galaxy[k].pos[1] = particles[j].pos[1];
	      galaxy[k].pos[2] = particles[j].pos[2];
	      galaxy[k].vel[0] = particles[j].vel[0];
	      galaxy[k].vel[1] = particles[j].vel[1];
	      galaxy[k].vel[2] = particles[j].vel[2];
	      galaxy[k].mass = particles[j].mass;
	      galaxy[k].pot = particles[j].pot;
	      
	      p += 1;
	      
	      //fprintf(pos,"%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",galaxy[k].mass,galaxy[k].pos[0],galaxy[k].pos[1],galaxy[k].pos[2],galaxy[k].vel[0],galaxy[k].vel[1],galaxy[k].vel[2]);
	      
	    }
	  
	  count += Nhost[type];
	  countf += Nhost[type+1];
	  
	  printf("type %d particle's pos&vels for host...DONE\n",type);
	  //fclose(pos);
	}
      
      
      /* cetrar a CM 0*/
      
      if(N_part[4]>0) //new stars
	{
	  //sprintf(buf,"4type_host.dat");
	  
	  //pos = fopen(buf,"w");

	  double r;

	  galaxy_aux = (aux *)malloc((size_t)Nguess*sizeof(aux)); //structure with size  =Nstars in host
	  
	  if(galaxy_aux == NULL)
	    {
	      printf("Allocation of galaxy_aux failed\n");
	      exit(0);
	    }
	 
	  indexmin = indexmax = 0;
	  
	  for( i=0; i<4; i++)
	    {
	      indexmin = indexmin + N_part[i];
	    }
	  
	  indexmax = indexmin + N_part[4];

	  printf("indexmax = %d \t indexmin = %d \n",indexmin,indexmax);
	    
	  printf("\ncounting new stars on host ...\n");
	 
	  for(j=indexmin; j<indexmax; j++)
	    {

	      r = sqrt(particles[j].pos[0]*particles[j].pos[0] + particles[j].pos[1]*particles[j].pos[1] + particles[j].pos[2]*particles[j].pos[2]);
	      
	      if(r<=Rsfr)
		{
		  galaxy_aux[count_star].id = particles[j].id;
		  galaxy_aux[count_star].pos[0] = particles[j].pos[0];
		  galaxy_aux[count_star].pos[1] = particles[j].pos[1];
		  galaxy_aux[count_star].pos[2] = particles[j].pos[2];
		  galaxy_aux[count_star].vel[0] = particles[j].vel[0];
		  galaxy_aux[count_star].vel[1] = particles[j].vel[1];
		  galaxy_aux[count_star].vel[2] = particles[j].vel[2];
		  galaxy_aux[count_star].mass = particles[j].mass;
		  galaxy_aux[count_star].pot = particles[j].pot;
		  
		  count_star += 1;
		}
	      
	    }
	  
	  Nhost[4] = count_star;
	  
	  p = 0;
	  
	  //LAURA
	  count = Nhost[0] + Nhost[1] + Nhost[2] + Nhost[3];
	  countf = count + Nhost[4];
	  
	  printf("count = %d \t countf = %d \n",count,countf);
      
	  for(k=count; k<countf; k++)
	    {
	      galaxy[k].pos[0] = galaxy_aux[p].pos[0];
	      galaxy[k].pos[1] = galaxy_aux[p].pos[1];
	      galaxy[k].pos[2] = galaxy_aux[p].pos[2];
	      galaxy[k].vel[0] = galaxy_aux[p].vel[0];
	      galaxy[k].vel[1] = galaxy_aux[p].vel[1];
	      galaxy[k].vel[2] = galaxy_aux[p].vel[2];
	      galaxy[k].mass = galaxy_aux[p].mass;
	      galaxy[k].pot = galaxy_aux[p].pot;

	      p += 1;
	   
	      //fprintf(pos,"%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",galaxy[k].mass,galaxy[k].pos[0], galaxy[k].pos[1],galaxy[k].pos[2],galaxy[k].vel[0], galaxy[k].vel[1],galaxy[k].vel[2]);
	    }
	  
	  printf("type %d particle's pos&vels for host...DONE\n",type);
	  //fclose(pos);
	  
	 
	  free(galaxy_aux);	  
	}

     
      N_galaxy_tot = Nhost[0] + Nhost[1] + Nhost[2] + Nhost[3] + Nhost[4];
      
      printf("\nNumber of stars in host = %d\n",count_star);
      printf("Total number of particles in host = %d\n",N_galaxy_tot);
      printf("\n (***) Nhost[0]=%d \t Nhost[1]=%d \t Nhost[2]=%d \t Nhost[3]=%d \t Nhost[4]=%d\n",Nhost[0],Nhost[1],Nhost[2],Nhost[3],Nhost[4]);
      
    }
  
  //===============
  // SATELLITE
  //===============

  if(flag==2) //satellite
    {

      Nguess = Nsat[0] + Nsat[4];
      Rsfr = f*r_scalength[1];
      N_galaxy0 = Nsat[0] + Nsat[1] + Nsat[2] + Nsat[3] + Nsat[4] + Nsat[5];
      count_gas = 0;
      count_star = 0;
      count = countf = 0;

      printf("\n (***) Nsat[0]=%d \t Nsat[1]=%d \t Nsat[2]=%d \t Nsat[3]=%d \t Nsat[4]=%d\n",Nsat[0],Nsat[1],Nsat[2],Nsat[3],Nsat[4]);
      
      galaxy = (galaxies *)malloc((size_t) N_galaxy0*sizeof(galaxies));
      
      if(galaxy == NULL)
	{
	  printf("Allocation of galaxy failed\n");
	  exit(0);
	}

       galaxy_aux = (aux *)malloc((size_t)N_galaxy0*sizeof(aux));

      if(galaxy_aux == NULL)
	{
	  printf("Allocation of galaxy failed\n");
	  exit(0);
	}
      
      //selecting gas particles 
      if(N_part[0]>0)    
	{
	  galaxy_gaspro = (galaxy_gas_properties *)malloc((size_t) Nsat[0]*sizeof(galaxy_gas_properties));
	  
	  if(galaxy_gaspro == NULL)
	    {
	      printf("Allocation of galaxy_gaspro failed\n");
	      exit(0);
	    }
	  
	  //sprintf(buf,"0type_sat.dat");
	  
	  //pos = fopen(buf,"w");
	  
	  for(i=0;i<N_part[0];i++)
	    {
	      
	      if(particles[i].id > Nhost[0])
		{
		  galaxy_aux[count_gas].id = particles[i].id;
		  galaxy_aux[count_gas].pos[0] = particles[i].pos[0];
		  galaxy_aux[count_gas].pos[1] = particles[i].pos[1];
		  galaxy_aux[count_gas].pos[2] = particles[i].pos[2];
		  galaxy_aux[count_gas].vel[0] = particles[i].vel[0];
		  galaxy_aux[count_gas].vel[1] = particles[i].vel[1];
		  galaxy_aux[count_gas].vel[2] = particles[i].vel[2];
		  galaxy_aux[count_gas].mass = particles[i].mass;

		  galaxy_aux[count_gas].rho = gaspro[i].rho;
		  galaxy_aux[count_gas].Ne = gaspro[i].Ne;
		  galaxy_aux[count_gas].Nh = gaspro[i].Nh;
		  galaxy_aux[count_gas].U = gaspro[i].U;
		  galaxy_aux[count_gas].sfr = gaspro[i].sfr;
		  galaxy_aux[count_gas].coolr = gaspro[i].coolr;
		  galaxy_aux[count_gas].U = gaspro[i].U;
		  galaxy_aux[count_gas].rho = gaspro[i].rho;
#ifdef MAGNETIC_FIELD
		  galaxy_aux[count_gas].B[0] = gaspro[i].B[0];
		  galaxy_aux[count_gas].B[1] = gaspro[i].B[1];
		  galaxy_aux[count_gas].B[2] = gaspro[i].B[2];
#endif	  
		  
		  count_gas += 1;
		}
	    }
	  
	  Nsat[0] = count_gas;
	  p = 0;
	  double tot_sfr=0.0;
	  FILE *sfr_file = NULL;
	  sfr_file = fopen("sf_rate_sat.dat","a");
	  for(k=0; k<count_gas; k++)
	    {
	      galaxy[k].pos[0] = galaxy_aux[k].pos[0];
	      galaxy[k].pos[1] = galaxy_aux[k].pos[1];
	      galaxy[k].pos[2] = galaxy_aux[k].pos[2];
	      galaxy[k].vel[0] = galaxy_aux[k].vel[0];
	      galaxy[k].vel[1] = galaxy_aux[k].vel[1];
	      galaxy[k].vel[2] = galaxy_aux[k].vel[2];
	      galaxy[k].mass   = galaxy_aux[k].mass;
	      galaxy[k].pot    = galaxy_aux[k].pot;
	      

	      galaxy_gaspro[p].rho = galaxy_aux[k].rho;
	      galaxy_gaspro[p].Ne  = galaxy_aux[k].Ne;
	      galaxy_gaspro[p].Nh  = galaxy_aux[k].Nh;
	      galaxy_gaspro[p].U   = galaxy_aux[k].U;
	      galaxy_gaspro[p].sfr = galaxy_aux[k].sfr;
	      galaxy_gaspro[p].coolr = galaxy_aux[k].coolr;
	      galaxy_gaspro[p].U   = galaxy_aux[k].U;
	      galaxy_gaspro[p].rho = galaxy_aux[k].rho;
	      
#ifdef MAGNETIC_FIELD
	      galaxy_gaspro[p].B[0]   = galaxy_aux[k].B[0];
	      galaxy_gaspro[p].B[1]   = galaxy_aux[k].B[1];
	      galaxy_gaspro[p].B[2]   = galaxy_aux[k].B[2];
#endif
	      tot_sfr += galaxy_aux[k].sfr; 
	      
	      p += 1;
     
	      //fprintf(pos,"%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",galaxy[k].mass, galaxy[k].pos[0], galaxy[k].pos[1], galaxy[k].pos[2],galaxy[k].vel[0], galaxy[k].vel[1], galaxy[k].vel[2]);
	    }

	  
	  printf("\ntype 0 particle's pos&vels for sat...DONE\n");	  
	  // fclose(pos);

	  fprintf(sfr_file,"%16.8lf %16.8lf\n",tot_sfr,Header.time);

	  //realloc;
	  //galaxy_aux = (aux*)realloc(galaxy_aux, count_gas*sizeof(aux));
	  free(galaxy_aux);
	}

      count = count_gas;
      countf = count + Nsat[1];

      for(type=1; type<4; type++) //halo, disk, bulge
	{
	  
	  indexmin = indexmax= 0;
	  
	  //sprintf(buf,"%dtype_sat.dat",type);
	  
	  //pos = fopen(buf,"w");
	  
	  p = 0;
	  
	  for( i=0; i<type; i++)
	    {
	      indexmin = indexmin + N_part[i];
	    }
	  
	  indexmin = indexmin + Nhost[type];
	  indexmax = indexmin + Nsat[type];
	  
	  for( k=count; k<countf; k++ )
	    {
	      j = indexmin + p;
	      galaxy[k].pos[0] = particles[j].pos[0];
	      galaxy[k].pos[1] = particles[j].pos[1];
	      galaxy[k].pos[2] = particles[j].pos[2];
	      galaxy[k].vel[0] = particles[j].vel[0];
	      galaxy[k].vel[1] = particles[j].vel[1];
	      galaxy[k].vel[2] = particles[j].vel[2];
	      galaxy[k].mass = particles[j].mass;
	      galaxy[k].pot = particles[j].pot;
	      p += 1;
	      
	      //fprintf(pos,"%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",galaxy[k].mass,galaxy[k].pos[0],galaxy[k].pos[1],galaxy[k].pos[2],galaxy[k].vel[0],galaxy[k].vel[1],galaxy[k].vel[2]);
	      
	    }
	  
	  printf("type %d particle's pos&vels for sat...DONE\n",type);
	  //fclose(pos);
	  
	  count += Nsat[type];
	  countf += Nsat[type+1];
	  
	}

      //=======================
      //new stars in sat
      //=======================
      
      
      type=4;
      if(N_part[4]>0) //new stars
	{
	  
	  //sprintf(buf,"%dtype_sat.dat",type);
	  
	  //pos = fopen(buf,"w");
	  
	  double r;

	  galaxy_aux = (aux *)malloc((size_t)Nguess*sizeof(aux));
	  
	  if(galaxy_aux == NULL)
	    {
	      printf("Allocation of galaxy_aux failed\n");
	      exit(0);
	    }

	  //============================================
	  //translate all new stars in snap (whole block)
	  //to the sat center of mass
	  //============================================

	  CM cmSat;
	  int n, imin, imax;
	  int *Ind;

	  n = Nsat[0] + Nsat[1] + Nsat[2];
	  imin = 0;
	  imax= Nsat[0] + Nsat[1] + Nsat[2];
	  
	  //==================================
	  //Computing center of mass of sat
	  //==================================
	  
	  //printf("Computing center of mass of satellite\n");
	  centerMass_mergercomp(NULL, n, imin, imax, &cmSat);
	  
	  //==========================================
	  //Traslating new stars to SAT center of mass
	  //==========================================
	  //===========================
	  //counting block of new stars
	  //===========================
	  
	  indexmin = indexmax = 0;
	  
	  for(i=0; i<type; i++)
	    {
	      indexmin = indexmin + N_part[i];
	    }
	  
	  indexmax = indexmin + N_part[type];
	  Ind = (int *)malloc((size_t)N_part[4]*sizeof(int));
	  if(Ind == NULL)
	    {
	      printf("Allocation of Ind failed\n");
	      exit(0);
	    }
	  
	  for(i=indexmin; i<indexmax; i++)
	    {
	      Ind[i-indexmin] =  i;

	      /*
		if((particles[i].id != i) || (particles[i].id - particles[i-1].id != 1))
		{
		prinft("sera esto?\n");
		exit(0);
		}
	      */
	    }
	  
	  //printf("Traslating galaxy to center of mass\n");
	  traslate(&cmSat, Ind, N_part[4]);
	  
	  free(Ind);
	 
	  
	  //printf("counting new stars on sat ...\n");

	  //  FILE *pf_juank=fopen("basura.juank","w");
	  
	  for(j=indexmin; j<indexmax; j++)
	    {
	      r = sqrt(particles[j].pos[0]*particles[j].pos[0] + particles[j].pos[1]*particles[j].pos[1] + particles[j].pos[2]*particles[j].pos[2]);

	      if(r<=Rsfr)
		{
		  galaxy_aux[count_star].id     = particles[j].id;
		  galaxy_aux[count_star].pos[0] = particles[j].pos[0];
		  galaxy_aux[count_star].pos[1] = particles[j].pos[1];
		  galaxy_aux[count_star].pos[2] = particles[j].pos[2];
		  galaxy_aux[count_star].vel[0] = particles[j].vel[0];
		  galaxy_aux[count_star].vel[1] = particles[j].vel[1];
		  galaxy_aux[count_star].vel[2] = particles[j].vel[2];
		  galaxy_aux[count_star].mass   = particles[j].mass;
		  galaxy_aux[count_star].pot    = particles[j].pot;
		  
		  //  fprintf(pf_juank,"%12d %12d %16.8f %16.8f %16.8f\n", j, galaxy_aux[count_star].id,
		  //galaxy_aux[count_star].pos[0]+cmSat.cm[0], galaxy_aux[count_star].pos[1]+cmSat.cm[1],
		  //	  galaxy_aux[count_star].pos[2]+cmSat.cm[2]);
		  
		  
		  count_star += 1;
		}
	      
	    }
	  
	  //fclose(pf_juank);

	  Nsat[4] = count_star;

	  //printf("\n (1.****SPLIT: cmSat) cmx = %16.8f \t cmy = %16.8f \t cmz = %16.8f\n",cmSat.cm[0], cmSat.cm[1], cmSat.cm[3]);
	  
	  p = 0;

	  //LAURA
	  count = Nsat[0] + Nsat[1] + Nsat[2] + Nsat[3];
	  
	  countf = count + Nsat[4];
	  
	  //printf("count = %d \t countf = %d \n",count,countf);
	  
	  for(k=count; k<countf; k++)
	    {
	      galaxy[k].pos[0] = galaxy_aux[p].pos[0] + cmSat.cm[0]; //returning partcl posit.to the original ref. system (not the sat-cm system)
	      galaxy[k].pos[1] = galaxy_aux[p].pos[1] + cmSat.cm[1];
	      galaxy[k].pos[2] = galaxy_aux[p].pos[2] + cmSat.cm[2];
	      galaxy[k].vel[0] = galaxy_aux[p].vel[0];
	      galaxy[k].vel[1] = galaxy_aux[p].vel[1];
	      galaxy[k].vel[2] = galaxy_aux[p].vel[2];
	      galaxy[k].mass   = galaxy_aux[p].mass;
	      galaxy[k].pot    = galaxy_aux[p].pot;
	      p += 1;
	     
	      //fprintf(pos,"%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n",galaxy[k].mass,
	      //      galaxy[k].pos[0]+cmSat.cm[0], galaxy[k].pos[1]+cmSat.cm[1], galaxy[k].pos[2]+cmSat.cm[2],
	      //      galaxy[k].vel[0], galaxy[k].vel[1],galaxy[k].vel[2]);
	      //printf("%16.8f %16.8f %16.8f\n",galaxy[k].pos[0],galaxy[k].pos[1],galaxy[k].pos[2]);
	    }
	  
	  printf("type %d particle's pos&vels for sat...DONE\n",type);
	  // fclose(pos);
      
       
	  free(galaxy_aux);
	  
	}
    
      N_galaxy_tot = Nsat[0] + Nsat[1] + Nsat[2] + Nsat[3] + Nsat[4];
      
      printf("\n (***) Nsat[0]=%d \t Nsat[1]=%d \t Nsat[2]=%d \t Nsat[3]=%d \t Nsat[4]=%d\n",Nsat[0],Nsat[1],Nsat[2],Nsat[3],Nsat[4]);
      
      printf("\nNumber of stars in satellite = %d\n",count_star);
      printf("Total number of particles in satellite = %d\n",N_galaxy_tot); 
    }  
  
  return 0;
}
