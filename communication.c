  #include <stdio.h>
  #include "ham2dtypes.h"
  #include "ham2dFunctionDefs.h"
  #include <stdlib.h>
  #define NQ 4
  #include <mpi.h>

  void communication(GRID *g, SOLN *s,int myid)
  {

  int i,j,ic;
  int iadjp,ist,ilengt,iend;
  int tag,error,ierr;
  int c1,c2,c3;
  
  int order;
 
  double *q1r,*q2r,*q3r,*q4r;
  double *q1s,*q2s,*q3s,*q4s;
  double *q1ra,*q2ra,*q3ra,*q4ra;
  double *q1sa,*q2sa,*q3sa,*q4sa;
  double *q1rb,*q2rb,*q3rb,*q4rb;
  double *q1sb,*q2sb,*q3sb,*q4sb;

  order = g->order;
  MPI_Request request[12];
  MPI_Status status[12];

  q1r = (double *) malloc(sizeof(double)*g->ncommu);
  q2r = (double *) malloc(sizeof(double)*g->ncommu);
  q3r = (double *) malloc(sizeof(double)*g->ncommu);
  q4r = (double *) malloc(sizeof(double)*g->ncommu);

  q1s = (double *) malloc(sizeof(double)*g->ncommu);
  q2s = (double *) malloc(sizeof(double)*g->ncommu);
  q3s = (double *) malloc(sizeof(double)*g->ncommu);
  q4s = (double *) malloc(sizeof(double)*g->ncommu);
 
  q1ra = (double *) malloc(sizeof(double)*g->ncommu);
  q2ra = (double *) malloc(sizeof(double)*g->ncommu);
  q3ra = (double *) malloc(sizeof(double)*g->ncommu);
  q4ra = (double *) malloc(sizeof(double)*g->ncommu);

  q1sa = (double *) malloc(sizeof(double)*g->ncommu);
  q2sa = (double *) malloc(sizeof(double)*g->ncommu);
  q3sa = (double *) malloc(sizeof(double)*g->ncommu);
  q4sa = (double *) malloc(sizeof(double)*g->ncommu);
 
  if(order==5)
  {
  q1rb = (double *) malloc(sizeof(double)*g->ncommu);
  q2rb = (double *) malloc(sizeof(double)*g->ncommu);
  q3rb = (double *) malloc(sizeof(double)*g->ncommu);
  q4rb = (double *) malloc(sizeof(double)*g->ncommu);

  q1sb = (double *) malloc(sizeof(double)*g->ncommu);
  q2sb = (double *) malloc(sizeof(double)*g->ncommu);
  q3sb = (double *) malloc(sizeof(double)*g->ncommu);
  q4sb = (double *) malloc(sizeof(double)*g->ncommu);
  }
  

  //=============================================
  // copy flow variables of current cell to send
  //=============================================
  for(i=0;i<g->nadjp;i++)
  {
    iadjp  = g->iadjp[i];
    ist    = g->istp[i];
    ilengt = g->ilengp[i];
    iend   = ist + ilengt;
    
    for(j=ist;j<iend;j++)
    { 
       c2  = g->irecvconn2[j*6+1];
       c1  = g->irecvconn2[j*6+2]; 

       if(order==5)
       {
       c3 = g->irecvconn2[j*6]; //third cell
       
       //third cell
       q1sb[j] = s->q[NVAR*c3];
       q2sb[j] = s->q[NVAR*c3+1];
       q3sb[j] = s->q[NVAR*c3+2];
       q4sb[j] = s->q[NVAR*c3+3];
       }
  
       // second cell 
       q1s[j] = s->q[NVAR*c2];
       q2s[j] = s->q[NVAR*c2+1];
       q3s[j] = s->q[NVAR*c2+2];
       q4s[j] = s->q[NVAR*c2+3];
       
       //first cell
       q1sa[j] = s->q[NVAR*c1];
       q2sa[j] = s->q[NVAR*c1+1];
       q3sa[j] = s->q[NVAR*c1+2];
       q4sa[j] = s->q[NVAR*c1+3];
    } 
  }


  for(i=0;i<g->nadjp;i++)
  {
  iadjp  = g->iadjp[i];
  ist    = g->istp[i];
  ilengt = g->ilengp[i];

  tag = 10*(2*(myid+1)+2*(iadjp+1));
  // second cell
  error = MPI_Irecv(&q1r[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[0]);
  tag = tag+1;
  error = MPI_Irecv(&q2r[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[1]);
tag = tag+1;
  error = MPI_Irecv(&q3r[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[2]);
tag = tag+1;
  error = MPI_Irecv(&q4r[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[3]);
  tag = tag+1;
  //first cell
  error = MPI_Irecv(&q1ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[4]);
  tag = tag +1;
  error = MPI_Irecv(&q2ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[5]);
tag = tag +1;
  error = MPI_Irecv(&q3ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[6]);
tag = tag +1;
  error = MPI_Irecv(&q4ra[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[7]);
  
  if(order==5){
  //third cell
  tag = tag+1;
  error = MPI_Irecv(&q1rb[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[8]);
  tag = tag +1;
  error = MPI_Irecv(&q2rb[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[9]);
tag = tag +1;
  error = MPI_Irecv(&q3rb[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[10]);
tag = tag +1;
  error = MPI_Irecv(&q4rb[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD,  &request[11]);
  }
  }

  for(i=0;i<g->nadjp;i++)
  {
  iadjp  = g->iadjp[i];
  ist    = g->istp[i];
  ilengt = g->ilengp[i];
  tag = 10*(2*(myid+1)+2*(iadjp+1));
 
  //second cell
  error = MPI_Send(&q1s[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD);
  tag = tag+1;
  error = MPI_Send(&q2s[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD);
tag = tag+1;
  error = MPI_Send(&q3s[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD);
tag = tag+1;
  error = MPI_Send(&q4s[ist],  ilengt, MPI_DOUBLE_PRECISION,  iadjp,  tag, MPI_COMM_WORLD);
 tag = tag+1; 
   //first cell
  error = MPI_Send(&q1sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&q2sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&q3sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&q4sa[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
  
  if(order==5){
  tag = tag +1;
  //thirsd cell
  error = MPI_Send(&q1sb[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&q2sb[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&q3sb[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
tag = tag +1;

  error = MPI_Send(&q4sb[ist], ilengt, MPI_DOUBLE_PRECISION, iadjp,  tag, MPI_COMM_WORLD);
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
        g->dpsil[c1][0]  = q1r[j]; 
        g->dpsil[c1][1]  = q2r[j];   
        g->dpsil[c1][2]  = q3r[j];   
        g->dpsil[c1][3]  = q4r[j];
        //first cell
        g->dpsila[c1][0]  = q1ra[j]; 
        g->dpsila[c1][1]  = q2ra[j];   
        g->dpsila[c1][2]  = q3ra[j];   
        g->dpsila[c1][3]  = q4ra[j];
        
        if(order==5){
        //third cell
        g->dpsilb[c1][0]  = q1rb[j]; 
        g->dpsilb[c1][1]  = q2rb[j];   
        g->dpsilb[c1][2]  = q3rb[j];   
        g->dpsilb[c1][3]  = q4rb[j];
        }

      }
      else
      {
        //second cell
        g->psil[c2][0]  = q1r[j]; 
        g->psil[c2][1]  = q2r[j];   
        g->psil[c2][2]  = q3r[j];   
        g->psil[c2][3]  = q4r[j];
        //first cell
        g->psila[c2][0]  = q1ra[j]; 
        g->psila[c2][1]  = q2ra[j];   
        g->psila[c2][2]  = q3ra[j];   
        g->psila[c2][3]  = q4ra[j];
        
        if(order==5){
        //third cell
        g->psilb[c2][0]  = q1rb[j]; 
        g->psilb[c2][1]  = q2rb[j];   
        g->psilb[c2][2]  = q3rb[j];   
        g->psilb[c2][3]  = q4rb[j];
        }
 
      }
    } 
  }
  ierr = MPI_Barrier(MPI_COMM_WORLD);

  free(q1r);
  free(q2r);
  free(q3r);
  free(q4r);

  free(q1s);
  free(q2s);
  free(q3s);
  free(q4s);

  free(q1ra);
  free(q2ra);
  free(q3ra);
  free(q4ra);

  free(q1sa);
  free(q2sa);
  free(q3sa);
  free(q4sa);

  if(order==5)
  {
  free(q1rb);
  free(q2rb);
  free(q3rb);
  free(q4rb);

  free(q1sb);
  free(q2sb);
  free(q3sb);
  free(q4sb);
  }

  
  // subroutine COMMUNICATION end
  }
