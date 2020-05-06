int indices(int comp)
{
  //====================================================
  //index for the specified galaxy (isolated, host, sat)
  //====================================================
  int i;
  
  //*********
  //Isolated
  //*********

  if(flag==0)
    {
      
      nDisk = N_part[comp]; 
      
      
      for( i=0; i<comp; i++)
	{
	  //indexmin = indexmax = indexmin + N_part[i];
	indexmin = indexmin + N_part[i];
	}
      indexmax = indexmin + N_part[comp];

      printf("Indexmin and Indexmax have been updated for isolated galaxy\n");
      printf("nDisk = %d \t indexmin = %d \t indexmax = %d \n",nDisk, indexmin,indexmax);
    
    }
  
  
  //*********
  //host
  //*********
  
  
  if(flag==1)
    {
      nDisk = Nhost[comp];
      
      for( i=0; i<comp; i++)
	{
	  //indexmin = indexmax = indexmin + N_part[i];
	  indexmin = indexmin + Nhost[i];
	}
      indexmax = indexmin + Nhost[comp];

      printf("Indexmin and Indexmax have been updated for host\n");
      printf("nDisk = %d \t indexmin = %d \t indexmax = %d \n",nDisk, indexmin,indexmax);
      
    }
  
  
  //*********
  //satellite
  //*********
  
  if(flag==2)
    {
      nDisk = Nsat[comp];
      
      for( i=0; i<comp; i++)
	{
	  //indexmin = indexmax = indexmin + N_part[i];
	  indexmin = indexmin + Nsat[i];
	}
      indexmax = indexmin + Nsat[comp];

      printf("Indexmin and Indexmax have been updated for satellite\n");
      printf("nDisk = %d \t indexmin = %d \t indexmax = %d \n",nDisk, indexmin,indexmax);
      
    }

 
  
  return 0;
}
