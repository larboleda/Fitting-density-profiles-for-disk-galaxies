int massprofile(double r_compare, int limlower, int limuper, int j)
{

  FILE *laufile=NULL;
  char buff[200];

  sprintf(buff,"%d_total_mass.dat",j);
  laufile = fopen(buff,"a"); //un snap en cada linea                                                                                                                                
  int i;
  double  tot_massr, tot_massR;
  double R; // R disk                                                                                                                           
  double r; //r sphere                                               

  tot_massr = tot_massR = 0.0;

  for(i=limlower; i<limuper; i++)
    {
      if(flag==0)
	{
	  R = sqrt(particles[i].pos[0]*particles[i].pos[0] + particles[i].pos[1]*particles[i].pos[1]);
	  r = sqrt(particles[i].pos[0]*particles[i].pos[0] + particles[i].pos[1]*particles[i].pos[1] + particles[i].pos[2]*particles[i].pos[2]);

	  if(R<=r_compare)
	    tot_massR += particles[i].mass;

	  if(r<=r_compare)
	    tot_massr += particles[i].mass;
	}
      else
	{
	  R = sqrt(galaxy[i].pos[0]*galaxy[i].pos[0] + galaxy[i].pos[1]*galaxy[i].pos[1]);
	  r = sqrt(galaxy[i].pos[0]*galaxy[i].pos[0] + galaxy[i].pos[1]*galaxy[i].pos[1] + galaxy[i].pos[2]*galaxy[i].pos[2]);

	  if(R<=r_compare)
	    tot_massR += galaxy[i].mass;

	  if(r<=r_compare)
	    tot_massr += galaxy[i].mass;
	}
    }

  printf("%lf %lf\n",tot_massR, tot_massr);
  fprintf(laufile,"%16.8e %16.8e %16.8f\n", tot_massR ,Header.time, tot_massr);

  fclose(laufile);

  return 0;
}
