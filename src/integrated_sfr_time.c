int integrated_sfr(int limlower, int limuper)
{

  FILE *laufile=NULL;

  laufile = fopen("integrated_sfr_time.dat","a"); //un snap en cada linea                                                               
  int i;
  double  tot_mass;
  double r; //r sphere
  int nstars;

  nstars = 0;

  printf("limlower = %d \t limuper=%d\n",limlower, limuper);
  for(i=limlower; i<limuper; i++)
    {
      if(flag==0)
	{
  
	  tot_mass += particles[i].mass; //summ mass of new stars
	  
	}
      else
	{
	  nstars += 1;
	  tot_mass += galaxy[i].mass;
	  
	}
    }

  double a;
  int snap;
  a = round(Header.time);
  snap = (int)a;
  
  printf("snap = %d \t nstars = %d\n",snap,nstars);
  //tot_mass = tot_mass/Deltat;
 
  fprintf(laufile,"%16.8e %16.8e\n",tot_mass,Header.time);
  
  
  fclose(laufile);

  return 0;
}


int integrated_sfr_gas_selection(int limlower, int limuper)
{
  
  FILE *laufile=NULL;
  FILE *filel = NULL;
  filel = fopen("bfield.dat","w");
  laufile = fopen("gas_integrated_sfr_time.dat","a"); //un snap en cada linea                                                               
  int i;
  double  tot_sfr=0.0,tot_mass=0.0,tot_coolrate=0.0;
  double tot_U=0.0, tot_rho=0.0,tot_bx = 0.0, tot_by = 0.0, tot_bz=0.0;
  double tot_B=0.0;
  double r; //r sphere
  int ngas, nselected;

  double rcomp;
  
  ngas = nselected = 0;

  printf("limlower = %d \t limuper=%d\n",limlower, limuper);
  for(i=limlower; i<limuper; i++)
    {
      if(flag==0)
	{
	 
	  r = sqrt(particles[i].pos[0]*particles[i].pos[0] + particles[i].pos[1]*particles[i].pos[1] + particles[i].pos[2]*particles[i].pos[2]);
	  
	       if(r<=rcomp) //SPHERE
		 {
		   tot_mass += particles[i].mass; //summ mass of new stars
		   tot_sfr += gaspro[i].sfr;
		   tot_coolrate += gaspro[i].coolr;
		   tot_U += gaspro[i].U;
		   tot_rho += gaspro[i].rho;
#ifdef MAGNETIC_FIELD
		   tot_bx += gaspro[i].B[0];
		   tot_by += gaspro[i].B[1];
		   tot_bz += gaspro[i].B[2];
#endif
		 }
	}
      else
	{
	 
	  r = sqrt(galaxy[i].pos[0]*galaxy[i].pos[0] + galaxy[i].pos[1]*galaxy[i].pos[1] + galaxy[i].pos[2]*galaxy[i].pos[2]);
	  
	  ngas += 1;
	  
	  if(flag==1)
	    {
	      rcomp = r_scalength[2]*r_scalength[0];
	    }
	  if(flag==2)
	    {
	      rcomp = r_scalength[2]*r_scalength[1];
	    }
	  
	  if(r<=rcomp)
	    {
	      tot_mass += galaxy[i].mass;
	      tot_sfr += galaxy_gaspro[i].sfr;
	      tot_coolrate += galaxy_gaspro[i].coolr;
	      tot_U += galaxy_gaspro[i].U;
	      tot_rho += galaxy_gaspro[i].rho;
#ifdef MAGNETIC_FIELD
	      tot_bx += galaxy_gaspro[i].B[0];
	      tot_by += galaxy_gaspro[i].B[1];
	      tot_bz += galaxy_gaspro[i].B[2];
	      
	      fprintf(filel,"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
		      galaxy[i].pos[0], galaxy[i].pos[1], galaxy[i].pos[2],
		      galaxy_gaspro[i].B[0], galaxy_gaspro[i].B[1], galaxy_gaspro[i].B[2], galaxy_gaspro[i].sfr, galaxy_gaspro[i].coolr, galaxy_gaspro[i].rho, galaxy_gaspro[i].U);
#else
	      fprintf(filel,"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
		      galaxy[i].pos[0], galaxy[i].pos[1], galaxy[i].pos[2],
		      galaxy_gaspro[i].sfr, galaxy_gaspro[i].coolr, galaxy_gaspro[i].rho, galaxy_gaspro[i].U);	  
#endif
	      nselected += 1;
	    }
	}
      
    }

  
#ifdef MAGNETIC_FIELD
  tot_B = sqrt(tot_bx*tot_bx + tot_by*tot_by + tot_bz*tot_bz);
  
  fprintf(laufile,"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",Header.time,tot_sfr,tot_mass,
	  tot_coolrate,tot_U,tot_rho,tot_bx,tot_by,tot_bz,tot_B);
#else

  fprintf(laufile,"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",Header.time,tot_sfr,tot_mass,tot_coolrate,tot_U,tot_rho);

#endif
  
  printf("\t ngas = %d \t nselected = %d\n",ngas,nselected);
  //tot_mass = tot_mass/Deltat;
  fclose(laufile);
  fclose(filel);
  
  return 0;
}

int integrated_sfr_gas(int limlower, int limuper)
{

  FILE *laufile=NULL;

  laufile = fopen("gas_integrated_sfr_time.dat","a"); //un snap en cada linea                                                               
  int i;
  double  tot_sfr=0.0,tot_mass=0.0,tot_coolrate=0.0;
  double tot_U=0.0, tot_rho=0.0,tot_bx = 0.0, tot_by = 0.0, tot_bz=0.0;
  double tot_B=0.0;
  double r; //r sphere
  int ngas;

  double rcomp;
  
  ngas = 0;

  printf("limlower = %d \t limuper=%d\n",limlower, limuper);
  for(i=limlower; i<limuper; i++)
    {
      if(flag==0)
	{
	
	  tot_mass += particles[i].mass; //summ mass of new stars
	  tot_sfr += gaspro[i].sfr;
	  tot_coolrate += gaspro[i].coolr;
	  tot_U += gaspro[i].U;
	  tot_rho += gaspro[i].rho;

#ifdef MAGNETIC_FIELD
	  tot_bx += gaspro[i].B[0];
	  tot_by += gaspro[i].B[1];
	  tot_bz += gaspro[i].B[2];
#endif
	  ngas += 1;
	  
	}
      else
	{ 
	  tot_mass += galaxy[i].mass;
	  tot_sfr += galaxy_gaspro[i].sfr;
	  tot_coolrate += galaxy_gaspro[i].coolr;
	  tot_U += galaxy_gaspro[i].U;
	  tot_rho += galaxy_gaspro[i].rho;
#ifdef MAGNETIC_FIELD
	  tot_bx += galaxy_gaspro[i].B[0];
	  tot_by += galaxy_gaspro[i].B[1];
	  tot_bz += galaxy_gaspro[i].B[2];
#endif
	  ngas += 1;
	}
    }

#ifdef MAGNETIC_FIELD
  tot_B = sqrt(tot_bx*tot_bx + tot_by*tot_by + tot_bz*tot_bz);
  
  fprintf(laufile,"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",Header.time,tot_sfr,tot_mass,
	  tot_coolrate,tot_U,tot_rho,tot_bx,tot_by,tot_bz,tot_B);
#else
  fprintf(laufile,"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",Header.time,tot_sfr,tot_mass,tot_coolrate,tot_U,tot_rho);

#endif
  
  printf("\t ngas = %d\n",ngas);
  //tot_mass = tot_mass/Deltat;
  fclose(laufile);

  return 0;
}
