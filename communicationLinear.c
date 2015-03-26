  #include <stdio.h>
  #include "ham2dtypes.h"
  #include "ham2dFunctionDefs.h"
  #include <stdlib.h>
  #define NQ 4
  #include <mpi.h>

  void communicationLinear(GRID *g, SOLN *s,int myid)
  {

  int i,j,ic;
  int iadjp,ist,ilengt,iend;
  int tag,error,ierr;
  int c1,c2,c3;
  int order;
  double *dq1r,*dq2r,*dq3r,*dq4r;
  double *dq1s,*dq2s,*dq3s,*dq4s;
  double *dq1ra,*dq2ra,*dq3ra,*dq4ra;
  double *dq1sa,*dq2sa,*dq3sa,*dq4sa;
  double *dq1rb,*dq2rb,*dq3rb,*dq4rb;
  double *dq1sb,*dq2sb,*dq3sb,*dq4sb;

  order = g->order;
  MPI_Request request[12];
  MPI_Status status[12];

  dq1r = (double *) malloc(sizeof(double)*g->ncommu);
  dq2r = (double *) malloc(sizeof(double)*g->ncommu);
  dq3r = (double *) malloc(sizeof(double)*g->ncommu);
  dq4r = (double *) malloc(sizeof(double)*g->ncommu);

  dq1s = (double *) malloc(sizeof(double)*g->ncommu);
  dq2s = (double *) malloc(sizeof(double)*g->ncommu);
  dq3s = (double *) malloc(sizeof(double)*g->ncommu);
  dq4s = (double *) malloc(sizeof(double)*g->ncommu);
  
  dq1ra = (double *) malloc(sizeof(double)*g->ncommu);
  dq2ra = (double *) malloc(sizeof(double)*g->ncommu);
  dq3ra = (double *) malloc(sizeof(double)*g->ncommu);
  dq4ra = (double *) malloc(sizeof(double)*g->ncommu);

  dq1sa = (double *) malloc(sizeof(double)*g->ncommu);
  dq2sa = (double *) malloc(sizeof(double)*g->ncommu);
  dq3sa = (double *) malloc(sizeof(double)*g->ncommu);
  dq4sa = (double *) malloc(sizeof(double)*g->ncommu);
 
  if(order==5)
  {
  dq1rb = (double *) malloc(sizeof(double)*g->ncommu);
  dq2rb = (double *) malloc(sizeof(double)*g->ncommu);
  dq3rb = (double *) malloc(sizeof(double)*g->ncommu);
  dq4rb = (double *) malloc(sizeof(double)*g->ncommu);

  dq1sb = (double *) malloc(sizeof(double)*g->ncommu);
  dq2sb = (double *) malloc(sizeof(double)*g->ncommu);
  dq3sb = (double *) malloc(sizeof(double)*g->ncommu);
  dq4sb = (double *) malloc(sizeof(double)*g->ncommu);
  }


  //copy to psil
  for(i=0;i<g->nadjp;i++)
  {
    iadjp  = g->iadjp[i];
    ist    = g->istp[i];
    ilengt = g->ilengp[i];
    iend   = ist + ilengt;
    
    for(j=ist;j<iend;j++)
    { 
       c2 = g->irecvconn2[j*6+1];
       c1 = g->irecvconn2[j*6+2];

       if(order==5)
       {
       c3 = g->irecvconn2[j*6];
       
       // third cell
       dq1sb[j] = s->dq[NVAR*c3];
       dq2sb[j] = s->dq[NVAR*c3+1];
       dq3sb[j] = s->dq[NVAR*c3+2];
       dq4sb[j] = s->dq[NVAR*c3+3];
       }

       // second cell 
       dq1s[j] = s->dq[NVAR*c2];
       dq2s[j] = s->dq[NVAR*c2+1];
       dq3s[j] = s->dq[NVAR*c2+2];
       dq4s[j] = s->dq[NVAR*c2+3];
       
       // first cell
       dq1sa[j] = s->dq[NVAR*c1];
       dq2sa[j] = s->dq[NVAR*c1+1];
       dq3sa[j] = s->dq[NVAR*c1+2];
       dq4sa[j] = s->dq[NVAR*c1+3];

    } 
  }


  for(i=0;i<g->nadjp;i++)
  {
  iadjp  = g->iadjp[i];
  ist    = g->istp[i];
  ilengt = g->ilengp[i];

  tag = 10*(2*(myid+1)+2*(iadjp+1));
  //second cell
  error = MPI_Irecv(&dq1r[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[0]);
  tag = tag+1;
  error = MPI_Irecv(&dq2r[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[1]);
tag = tag+1;
  error = MPI_Irecv(&dq3r[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[2]);
tag = tag+1;
  error = MPI_Irecv(&dq4r[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[3]);
 tag = tag+1; 
 //first cell 
    error = MPI_Irecv(&dq1ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[4]);
  tag = tag +1;
  error = MPI_Irecv(&dq2ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[5]);
tag = tag +1;

  error = MPI_Irecv(&dq3ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[6]);
tag = tag +1;

  error = MPI_Irecv(&dq4ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[7]);

if(order==5){
tag = tag +1;
    error = MPI_Irecv(&dq1rb[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[8]);
  tag = tag +1;
  error = MPI_Irecv(&dq2rb[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[9]);
tag = tag +1;

  error = MPI_Irecv(&dq3rb[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[10]);
tag = tag +1;

  error = MPI_Irecv(&dq4rb[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[11]);
  
  }
  }

  for(i=0;i<g->nadjp;i++)
  {
  iadjp  = g->iadjp[i];
  ist    = g->istp[i];
  ilengt = g->ilengp[i];
  tag = 10*(2*(myid+1)+2*(iadjp+1));
 
  //second cell
  error = MPI_Send(&dq1s[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD);
  tag = tag+1;
  error = MPI_Send(&dq2s[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD);
tag = tag+1;
  error = MPI_Send(&dq3s[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD);
tag = tag+1;
  error = MPI_Send(&dq4s[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD);
  //first cell
  tag = tag+1;
  error = MPI_Send(&dq1sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&dq2sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&dq3sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&dq4sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
  
  if(order==5){
  tag = tag +1;
  //third cell
  error = MPI_Send(&dq1sb[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&dq2sb[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&dq3sb[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&dq4sb[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
  }
  }
  
  for(i=0;i<g->nadjp;i++)
  {
    error=MPI_Wait(&request[0],&status[0]);
    error=MPI_Wait(&request[1],&status[1]);
    error=MPI_Wait(&request[2],&status[2]);
    error=MPI_Wait(&request[3],&status[3]);
    error=MPI_Wait(&request[4],&status[4]);
    error=MPI_Wait(&request[5],&status[5]);
    error=MPI_Wait(&request[6],&status[6]);
    error=MPI_Wait(&request[7],&status[7]);
    if(order==5){
    error=MPI_Wait(&request[8],&status[8]);
    error=MPI_Wait(&request[9],&status[9]);
    error=MPI_Wait(&request[10],&status[10]);
    error=MPI_Wait(&request[11],&status[11]);
    }
  }

 ierr = MPI_Barrier(MPI_COMM_WORLD);
 
  // copy from receive data to designated cell 
  for(i=0;i<g->nadjp;i++)
  {
    iadjp  = g->iadjp[i];
    ist    = g->istp[i];
    ilengt = g->ilengp[i];
    iend   = ist + ilengt;
    for(j=ist;j<iend;j++)
    {
      c2 = g->irecv[j*2]; // second cell index            
      c1 = g->irecvconn[j];
      
      // for duplication case
      if(g->idup[c2]==2)
      {
        //second cell
        g->dpsil_deri[c1][0]  = dq1r[j]; 
        g->dpsil_deri[c1][1]  = dq2r[j];   
        g->dpsil_deri[c1][2]  = dq3r[j];   
        g->dpsil_deri[c1][3]  = dq4r[j];
        //fist cell
        g->dpsila_deri[c1][0]  = dq1ra[j]; 
        g->dpsila_deri[c1][1]  = dq2ra[j];   
        g->dpsila_deri[c1][2]  = dq3ra[j];   
        g->dpsila_deri[c1][3]  = dq4ra[j];
        
        if(order==5){
        //third cell
        g->dpsilb_deri[c1][0]  = dq1rb[j]; 
        g->dpsilb_deri[c1][1]  = dq2rb[j];   
        g->dpsilb_deri[c1][2]  = dq3rb[j];   
        g->dpsilb_deri[c1][3]  = dq4rb[j];
        }
      }
      else
      {
        //second cell
        g->psil_deri[c2][0]  = dq1r[j]; 
        g->psil_deri[c2][1]  = dq2r[j];   
        g->psil_deri[c2][2]  = dq3r[j];   
        g->psil_deri[c2][3]  = dq4r[j];
        //first cell
        g->psila_deri[c2][0]  = dq1ra[j]; 
        g->psila_deri[c2][1]  = dq2ra[j];   
        g->psila_deri[c2][2]  = dq3ra[j];   
        g->psila_deri[c2][3]  = dq4ra[j];
        
        if(order==5){
        //third cell
        g->psilb_deri[c2][0]  = dq1rb[j]; 
        g->psilb_deri[c2][1]  = dq2rb[j];   
        g->psilb_deri[c2][2]  = dq3rb[j];   
        g->psilb_deri[c2][3]  = dq4rb[j];
        }

      }

    } 
  }
 
  free(dq1r);
  free(dq2r);
  free(dq3r);
  free(dq4r);

  free(dq1s);
  free(dq2s);
  free(dq3s);
  free(dq4s);

  free(dq1ra);
  free(dq2ra);
  free(dq3ra);
  free(dq4ra);

  free(dq1sa);
  free(dq2sa);
  free(dq3sa);
  free(dq4sa);
 
  if(order==5)
  {
  free(dq1rb);
  free(dq2rb);
  free(dq3rb);
  free(dq4rb);

  free(dq1sb);
  free(dq2sb);
  free(dq3sb);
  free(dq4sb);
  }

 
  // subroutine COMMUNICATION end
  }
