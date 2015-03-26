// ##################################################################
//
// ham2d.c
//
// Driver code for the 2d unstrctured Hamiltonian path based code
//
// Written by Dr. Jayanarayanan Sitaraman
// #################################################################
// Version 6.0
// Modified by Yong Su Jung(12.20.2014)
//  1. Periodic bc
//  2. WENO 5 reconstruction
//  3. Dual time stepping
//  4. Restart solution 
//
//###################################################################
// Bug fix report
// computeRHS.c (321 , 389 line),2014.11.20
// if(rightCell>-1 || g->test==1) -> if(rightCell>-1)
// computeRHSk.c (289,348) : same
//===============================================================
//
// ##################################################################
#include <stdio.h>
#include <stdlib.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#include <time.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char **argv)
{
  GRID    *g;     // pointer to GRID data structure
  SOLN    *s;     // pointer to SOLN data structure
  int     ngrids; // total number of grids
  int     nsteps,ntotal; // total number of time steps
  int     nwrite; // total number of time steps
  int     i,n,nn;    //
  int     irest,nrest; //restart
  int     idual,ndual;
  double  dt;     // time step
  double  CFL;
  double  l2rho;
  int     msweep;
  char    c;
  char    fname[80];
  char    scheme[20];
  char    timeInteg[20];
  int     order,timeacc,visc,test;
  FILE    *fp;
  clock_t start, end;
  double  cpu_time_used;
    
  //mpi
  int nproc,myid,ierr;

  ierr = MPI_Init(&argc,&argc);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  printf("Hello world! I'm process %i out of %i process\n",myid,nproc);
  ierr = MPI_Barrier(MPI_COMM_WORLD);

  // ==================================================================  
// Read in file inputs from input.ham2d
// ==================================================================
  fp=fopen("input.ham2d","r");
  while((c=fgetc(fp))!='\n'); // skip the 
  while((c=fgetc(fp))!='\n'); // first three
  while((c=fgetc(fp))!='\n'); // lines of input.ham2d
  fscanf(fp,"scheme=%s\n",scheme);
  fscanf(fp,"time integration=%s\n",timeInteg);
  fscanf(fp,"order=%d\n",&order);
  fscanf(fp,"timeacc=%d\n",&timeacc);
  fscanf(fp,"nsteps=%d\n",&nsteps);
  fscanf(fp,"nwrite=%d\n",&nwrite);
  fscanf(fp,"dt=%lf\n",&dt);
  fscanf(fp,"CFL=%lf\n",&CFL);
  fscanf(fp,"msweep=%d\n",&msweep);
  fscanf(fp,"visc=%d\n",&visc);
  fscanf(fp,"testcase=%d\n",&test);
  fscanf(fp,"irest=%d\n",&irest);
  fscanf(fp,"nrest=%d\n",&nrest);
  fscanf(fp,"idual=%d\n",&idual);
  fscanf(fp,"ndual=%d\n",&ndual);
  fclose(fp);

  if(myid==0) 
  {
  trace(nsteps);
  tracef(dt);
  }

  //
  ngrids=1;
  g=(GRID *) malloc(sizeof(GRID)*ngrids);
  s=(SOLN *) malloc(sizeof(SOLN)*ngrids);
  
// ==================================================================
// preprocess grids
// code is written with an overset
// framework in mind (but not implemented yet)
// (ngrid==1) for now
// ==================================================================
  for(i=0;i<ngrids;i++) 
  {
      
      // Read the grid data (i.e., data points and connectivity information
      // from the Matlab mesh-generation)
      readGrid(&g[i],myid,nproc);
      g[i].visc = visc;
      
      // Run the preprocessing steps (computes cell vol., and
      // obtains cell to face and cell to chain connectivity)
      g[i].test = test;
      preprocess(&g[i],myid);

      // Initialize the flow
      initflow(&g[i],&s[i],irest,myid);
      
      // Initialize other grid and solution variables
      g[i].order   = order;   // order of the scheme
      g[i].timeacc = timeacc; // time accuracy
      g[i].CFL     = CFL;     // CFL number for grid
      s[i].cflnum  = CFL;     // CFL number for soln
      g[i].msweep  = msweep;  // number of sweeps 
      s[i].idual   = idual;  // Dual time stepping
      s[i].ndual   = ndual;  // Nsub iteration

      if (strcmp(timeInteg,"bdf1")==0) g[i].timeInteg=BDF1;
      if (strcmp(timeInteg,"bdf2")==0) g[i].timeInteg=BDF2;
  }

  //nn = 1;
  //outputSolution(&g[0],&s[0],nn,myid,nproc); 
  //ierr = MPI_Barrier(MPI_COMM_WORLD);
  //exit(1);
  if(myid==0) tracef(CFL)
// ==================================================================
// Main bulk of the code. Step through the solution for 
// as many steps are required.
// ==================================================================
  if(myid<10)
  sprintf(fname,"./output/sol_his00%i.dat",myid);
  else if(myid<100)
  sprintf(fname,"./output/sol_his0%i.dat",myid);
  else if(myid<1000)
  sprintf(fname,"./output/sol_his%i.dat",myid);


  fp = fopen(fname,"w"); 

  if(myid==0) 
  {   
    printf("#ham2d : using %s scheme for inversion\n",scheme);
    printf("#ham2d : using %d order of spatial accurate for invicid\n",order);
  }
  cpu_time_used=0;

  ierr = MPI_Barrier(MPI_COMM_WORLD);
  ntotal = 0;
  if(irest==1) ntotal = s->nt;
  for (n=ntotal;n<ntotal+nsteps;n++) // loop through the time steps
  { 
    for(i=0;i<ngrids;i++) 
 	  { 
       // start the clock
  	    start = clock();

       // Step through the solution for one time step
  	    stepSolution(scheme,&g[i],&s[i],dt,&l2rho,myid);

       // compute the force on the airfoil
  	    computeForce(&g[i],&s[i]);

       // end the clock
  	    end = clock();

       // compute CPU time used
  	    cpu_time_used += (((double) (end - start)) / CLOCKS_PER_SEC);
     
       if(myid==0) printf("%d %e %2.4f %2.4f %2.4f\n",n,l2rho,s[i].cl,s[i].cd,cpu_time_used);
       fprintf(fp,"%d %e %2.4f %2.4f %2.4f\n",n,l2rho,s[i].cl,s[i].cd,cpu_time_used);
       
	  } 
     ierr = MPI_Barrier(MPI_COMM_WORLD);

     // write output files
     if((n+1)%nwrite==0 || n==ntotal) 
     {
       nn = (n+1)/nwrite;
       outputSolution(&g[0],&s[0],nn,myid,nproc); 
       ierr = MPI_Barrier(MPI_COMM_WORLD);
     }
     
     //// write restart files
     //if((n+1)%nrest==0)
     //{
     //  nn = (n+1)/nrest;
     //  wrest(&g[0],&s[0],n,nn);
     //}
 
     //printf("myid:%d, iter:%d\n",myid,n);

   }
   fclose(fp);
  // Output the solution to screen and files
//  outputSolution(&g[0],&s[0]); 
   printf("finish program,myid:%d\n",myid);
   ierr = MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   return(0);

}
// ##################################################################
// END OF FILE
// ##################################################################
