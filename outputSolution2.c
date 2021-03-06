// ##################################################################
//
// outputSolution.c
//
// Subroutines for output plots
//
// Written by Dr. Jayanarayanan Sitaraman
// ##################################################################
#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#include <mpi.h>
// ==================================================================
// 
// ==================================================================
void outputSolution2(GRID *g,SOLN *s,int nn,int myid, int nproc)
{
  int    i,n;
  int    iface,node1,node2;
  int    icell;
  double x1,y1,x2,y2,rho,rhou,rhov,e,pp,cp;
  FILE   *fp,*fp1;
  char   fname[80]={"output/output/dat"};
  int    iproc,ierr,ncell_prev;
  
  //==============\
  // Surface plot  >
  //==============/


  for(iproc=0;iproc<nproc;iproc++)
  {
    if(iproc==0 && myid==0)
    {
      fp=fopen(fname,"w");
      fprintf(fp,"VARIABLES = \"X\",\"Y\",\"RHO\",\"RHOU\",\"RHOV\",\"E\"\n");
      
      fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N=%d E=%d DATAPACKING=BLOCK\n",g->ntnodes,g->ntcells); 

      fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=CELLCENTERED, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED)\n");
      fclose(fp);

      printf("myid:%d",myid);
    }//myid==0
    ierr = MPI_Barrier(MPI_COMM_WORLD);  
   }


  for(iproc=0;iproc<nproc;iproc++)
  {
    if(iproc==myid)
    {
       //append file
      fp=fopen(fname,"a+");
      // list of x-coordinate positions (after smoothing)
      for(i=0;i<g->nnodes;i++)
        fprintf(fp,"%f\n",g->x[2*i]);
        if(myid==nproc-1) fprintf(fp,"\n");
        fclose(fp);
   }
  } 
   ierr = MPI_Barrier(MPI_COMM_WORLD);       
   
   
   for(iproc=0;iproc<nproc;iproc++)
  {
    if(iproc==myid)
    {
      // append file
      fp=fopen(fname,"a+");
      // list of y-coordinate positions (after smoothing)
      for(i=0;i<g->nnodes;i++)
        fprintf(fp,"%f\n",g->x[2*i+1]);
        
        if(myid==nproc-1) fprintf(fp,"\n");
        fclose(fp);
    }
   }
    


   ierr = MPI_Barrier(MPI_COMM_WORLD);       
   for(iproc=0;iproc<nproc;iproc++)
  {
    if(iproc==myid)
    {
      // append file
      fp=fopen(fname,"a+");
      // list of y-coordinate positions (after smoothing)
        for (i=0;i<g->ncells;i++)
        {
           fprintf(fp,"%f\n",s->q[4*i]);
        }
        
        if(myid==nproc-1) fprintf(fp,"\n");
        fclose(fp);

    }
   }
      
    ierr = MPI_Barrier(MPI_COMM_WORLD);       
   for(iproc=0;iproc<nproc;iproc++)
  {
    if(iproc==myid)
    {
      // append file
      fp=fopen(fname,"a+");
      // list of y-coordinate positions (after smoothing)
        for (i=0;i<g->ncells;i++)
        {
           fprintf(fp,"%f\n",s->q[4*i+1]);
        }
        
        if(myid==nproc-1) fprintf(fp,"\n");
        fclose(fp);
    }   
   }
    
     ierr = MPI_Barrier(MPI_COMM_WORLD);       
   for(iproc=0;iproc<nproc;iproc++)
  {
    if(iproc==myid)
    {
      // append file
      fp=fopen(fname,"a+");
      // list of y-coordinate positions (after smoothing)
        for (i=0;i<g->ncells;i++)
        {
           fprintf(fp,"%f\n",s->q[4*i+2]);
        }
          
        if(myid==nproc-1) fprintf(fp,"\n");
        fclose(fp);
    }
   }
      

     ierr = MPI_Barrier(MPI_COMM_WORLD);       
   for(iproc=0;iproc<nproc;iproc++)
  {
    if(iproc==myid)
    {
      // append file
      fp=fopen(fname,"a+");
      // list of y-coordinate positions (after smoothing)
        for (i=0;i<g->ncells;i++)
        {
           fprintf(fp,"%f\n",s->q[4*i+3]);
        }
        
        if(myid==nproc-1) fprintf(fp,"\n");
        fclose(fp);
    }
   }
   ierr = MPI_Barrier(MPI_COMM_WORLD);  
   
   ncell_prev = 0;
   for(iproc=0;iproc<nproc;iproc++)
  {
    if(iproc==myid)
    {
      // append file
      fp=fopen(fname,"a+");
       // connectivity information for the cells            
        for(i=0;i<g->ncells;i++)
        fprintf(fp,"%d %d %d %d\n",g->conn[4*i]+1+ncell_prev,g->conn[4*i+1]+1+ncell_prev,g->conn[4*i+2]+1+ncell_prev,g->conn[4*i+3]+1+ncell_prev);

        ncell_prev = ncell_prev + g->ncells;
        fclose(fp);
    }

   }


//    }
//  } //iproc

//
//  //=========\
//  // cp curve >
//  //=========/
//  fp=fopen("./output/cp.dat","w");
//  for(i=0;i<g->nbfaces;i++)
//   {
//    iface = g->bfaces[i];
//    node1 = g->faces[6*iface];
//    node2 = g->faces[6*iface+1];
//    x1    = g->x[2*node1];
//    y1    = g->x[2*node1+1];
//    x2    = g->x[2*node2];
//    y2    = g->x[2*node2+1];
//    icell = g->faces[6*iface+2];
//    rho   = s->q[NVAR*icell];
//    rhou  = s->q[NVAR*icell+1];
//    rhov  = s->q[NVAR*icell+2];
//    e     = s->q[NVAR*icell+3];
//    pp    = (gamm-1)*(e-0.5*(rhou*rhou+rhov*rhov)/rho);
//    cp    = (pp-pinf)/(0.5*s->mach*s->mach);	
//
//    fprintf(fp,"%.16e %.16e %.16e\n",(x1+x2)*0.5,(y1+y2)*0.5,cp);
//  }
//  fclose(fp);
//  //fclose(fp1);
//  printf("Writing Output files.......\n");
}

// ==================================================================
// 
// ==================================================================
// ##################################################################
// END OF FILE
// ##################################################################s
