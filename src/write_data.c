
int print_data_toascii() //flag 1

{
  int i, type;
  char outfile[500];

  N_min = 0;
   
  for( type=0; type<6; type++)
    {

      if(flag==1)
	{
	
	  N_max = N_min + Nhost[type];	
	  
	  if(Nhost[type] > 0)
	    {
	      sprintf(outfile,"host_%d.dat",type);
	      FILE *outfiles;
	      outfiles = fopen(outfile,"w");
	      if(outfiles==NULL) printf("Could not open file %s\n",outfile);  
	      for(i=N_min;i<N_max;i++)
		{
		  
		  fprintf(outfiles,"%f %f %f %f %f %f %.8f\n",
			  galaxy[i].pos[0],galaxy[i].pos[1],galaxy[i].pos[2],
			  galaxy[i].vel[0],galaxy[i].vel[1],galaxy[i].vel[2],
			  galaxy[i].mass);
		
		  /*
		    #ifdef PRINT_ALL_ASCII
	      
	      if( type == 0 &&  Nsat[0] >0 )
		{
		  fprintf(outfiles," %f %f",galaxy_gaspro[i].U,galaxy_gaspro[i].rho);
#ifdef COOLING
		  fprintf(outfiles," %f %f",galaxy_gaspro[i].Ne,galaxy_gaspro[i].Nh);
#endif
		  fprintf(outfiles," %f",galaxy_gaspro[i].h);
#ifdef SFR
		  fprintf(outfiles," %f",galaxy_gaspro[i].sfr);
#endif
		}

#ifdef STELLARAGE	      
	      if( type == 4 &&  Nsat[4] >0 )
		fprintf(outfiles," %f",galaxy[i].stellar_age);
#endif

#ifdef METALS
	      if( ((type == 0) &&  (Nsat[0]>0)) || ((type == 4) &&  (Nsat[4]>0)) )
		fprintf(outfiles," %f",galaxy[i].metallicity);
#endif

#endif // close PRINT_ALL_ASCII
	      
#ifdef OUTPUTPOTENTIAL
	      fprintf(outfiles," %f",galaxy[i].pot);
#endif

#ifdef PRINT_ALL_ASCII
	      
#ifdef OUTPUTACCELERATION
	      fprintf(outfiles," %f %f %f",galaxy[i].acce[X],galaxy[i].acce[Y],galaxy[i].acce[Z]);
#endif
	      
#ifdef OUTPUTCHANGEOFENTROP
	      if( (type == 0) && (Nsat[0] >0) )
		fprintf(outfiles," %f",galaxy_gaspro[i].ecr);
#endif 
	      
#ifdef OUTPUTTIMESTEP
	      fprintf(outfiles," %f",galaxy[i].timestep);
#endif
	   
#ifdef OUTPUTCOOLRATE
	      if( (type == 0) && (Nsat[0] >0) )
		fprintf(outfiles," %f",galaxy_gaspro[i].coolr);
#endif
   
#ifdef OUTPUT_DIV_CURL
	      if( (type == 0) && (Nsat[0] >0) )
		fprintf(outfiles," %f %f",galaxy_gaspro[i].div,galaxy_gaspro[i].curl);
#endif

#ifdef OUTPUT_VORTICITY
	      fprintf(outfiles," %f %f %f",galaxy_gaspro[i].vorticity[X],galaxy_gaspro[i].vorticity[Y],galaxy_gaspro[i].vorticity[Z]);
#endif

#ifdef WINDS
	      fprintf(outfiles," %f",galaxy_gaspro[i].delaytime);
#endif

#endif //close PRINT_ALL_ASCII*/

	      //fprintf(outfiles,"\n");	
		  
	    }
	      fclose(outfiles);
	    }
	  
	  
	  N_min = N_max;
	 
	  
	}
     

       if(flag==2)
	{
	
	  N_max = N_min + Nsat[type];	
	  
	  if(Nsat[type] > 0)
	    {
	      sprintf(outfile,"satellite_%d.dat",type);
	      FILE *outfiles;
	      outfiles = fopen(outfile,"w");
	      if(outfiles==NULL) printf("Could not open file %s\n",outfile);  
	      for(i=N_min;i<N_max;i++)
		{
		  
		  fprintf(outfiles,"%f %f %f %f %f %f %.8f\n",
			  galaxy[i].pos[0],galaxy[i].pos[1],galaxy[i].pos[2],
			  galaxy[i].vel[0],galaxy[i].vel[1],galaxy[i].vel[2],
			  galaxy[i].mass);
		
		  /*
		    #ifdef PRINT_ALL_ASCII
	      
	      if( type == 0 &&  Nsat[0] >0 )
		{
		  fprintf(outfiles," %f %f",galaxy_gaspro[i].U,galaxy_gaspro[i].rho);
#ifdef COOLING
		  fprintf(outfiles," %f %f",galaxy_gaspro[i].Ne,galaxy_gaspro[i].Nh);
#endif
		  fprintf(outfiles," %f",galaxy_gaspro[i].h);
#ifdef SFR
		  fprintf(outfiles," %f",galaxy_gaspro[i].sfr);
#endif
		}

#ifdef STELLARAGE	      
	      if( type == 4 &&  Nsat[4] >0 )
		fprintf(outfiles," %f",galaxy[i].stellar_age);
#endif

#ifdef METALS
	      if( ((type == 0) &&  (Nsat[0]>0)) || ((type == 4) &&  (Nsat[4]>0)) )
		fprintf(outfiles," %f",galaxy[i].metallicity);
#endif

#endif // close PRINT_ALL_ASCII
	      
#ifdef OUTPUTPOTENTIAL
	      fprintf(outfiles," %f",galaxy[i].pot);
#endif

#ifdef PRINT_ALL_ASCII
	      
#ifdef OUTPUTACCELERATION
	      fprintf(outfiles," %f %f %f",galaxy[i].acce[X],galaxy[i].acce[Y],galaxy[i].acce[Z]);
#endif
	      
#ifdef OUTPUTCHANGEOFENTROP
	      if( (type == 0) && (Nsat[0] >0) )
		fprintf(outfiles," %f",galaxy_gaspro[i].ecr);
#endif 
	      
#ifdef OUTPUTTIMESTEP
	      fprintf(outfiles," %f",galaxy[i].timestep);
#endif
	   
#ifdef OUTPUTCOOLRATE
	      if( (type == 0) && (Nsat[0] >0) )
		fprintf(outfiles," %f",galaxy_gaspro[i].coolr);
#endif
   
#ifdef OUTPUT_DIV_CURL
	      if( (type == 0) && (Nsat[0] >0) )
		fprintf(outfiles," %f %f",galaxy_gaspro[i].div,galaxy_gaspro[i].curl);
#endif

#ifdef OUTPUT_VORTICITY
	      fprintf(outfiles," %f %f %f",galaxy_gaspro[i].vorticity[X],galaxy_gaspro[i].vorticity[Y],galaxy_gaspro[i].vorticity[Z]);
#endif

#ifdef WINDS
	      fprintf(outfiles," %f",galaxy_gaspro[i].delaytime);
#endif

#endif //close PRINT_ALL_ASCII*/

	      //fprintf(outfiles,"\n");	
		  
	    }
	      fclose(outfiles);
	    }
	  
	  
	  N_min = N_max;
	 
	  
	}
  
    }
  if(flag==1)
    {
      printf("\nDone printing data for host.\n\n");
    }
  if(flag==2)
    {
      printf("\nDone printing data for satellite.\n\n");
    }
  
  return 0;
}
