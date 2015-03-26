// ##################################################################
//
// readGrid.c
//
// Function that reads the grid based on the outputs of the Matlab
// mesh generation code
//
// Data files read from:
//  - coord.dat   (data coordinate points)
//  - conn.dat    (node connectivity information for the triangles)
//  - qedges.dat  (writes the QEdge matrix, 6 cols, refer documentation)
//  - ncolors.dat (number of loops of each colour)
//  - iqloops.dat (index of inner and outer loops of each node)
//  - qloops.dat  (Inner and outer loops  for each node)
//
// Written by Dr. Jayanarayanan Sitaraman
// ##################################################################
#include <stdio.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#include <stdlib.h>
#define NQ 4
#include <mpi.h>

void readGrid(GRID *g,int myid, int nproc)
{
   FILE   *fp;
   char   coord[40],conn[40],qedge[40],qloop[40],iqloop[40],ncolor[40];
   char   c, commu[50], domain[50];
   int    i,l1,l2,l3,j;
   int    workSize;
   double fdummy,tmp;
   int    ierr,iproc;
   int    dummy,c2,c1;
   int    k,kk,ndup; 
   int    *dup;

// total domain
   //if(myid==0) 
   //{
   //  g->ntnodes = 0; 
   //  for(iproc=0;iproc<nproc;iproc++)
   //  {
   //     sprintf(coord,"./QuadData/domain00%i/coord.dat",iproc);
   //     fp=fopen(coord,"r"); 
   //     fscanf(fp,"%d",&(g->tmp));
   //     fclose(fp);
   //     g->ntnodes = g->ntnodes + g->tmp;
   //     
   //     sprintf(conn,"./QuadData/domain00%i/conn.dat",iproc);
   //     fp=fopen(conn,"r");
   //     fscanf(fp,"%d",&(g->tmp));
   //     fclose(fp);
   //     g->ntcells = g->ntcells + g->tmp;
   //  }
   //}
  
// =================================================================
// Initialize dimension
// ================================================================

    if(myid<10)
    {
    sprintf(coord, "./QuadData/domain00%i/coord.dat",myid);
    sprintf(conn,  "./QuadData/domain00%i/conn.dat",myid);
    sprintf(qedge, "./QuadData/domain00%i/qedges.dat",myid);
    sprintf(ncolor,"./QuadData/domain00%i/ncolors.dat",myid);
    sprintf(iqloop,"./QuadData/domain00%i/iqloops.dat",myid);
    sprintf(qloop, "./QuadData/domain00%i/qloops.dat",myid);
    sprintf(commu, "./QuadData/domain00%i/quadCommu00%i.dat",myid,myid);
    sprintf(domain, "./QuadData/domain00%i/subdomain00%i.dat",myid,myid);

    }
    else if(10<=myid<100)
    {
    sprintf(coord, "./QuadData/domain0%i/coord.dat",myid);
    sprintf(conn,  "./QuadData/domain0%i/conn.dat",myid);
    sprintf(qedge, "./QuadData/domain0%i/qedges.dat",myid);
    sprintf(ncolor,"./QuadData/domain0%i/ncolors.dat",myid);
    sprintf(iqloop,"./QuadData/domain0%i/iqloops.dat",myid);
    sprintf(qloop, "./QuadData/domain0%i/qloops.dat",myid);
    sprintf(commu, "./QuadData/domain0%i/quadCommu0%i.dat",myid,myid);
    sprintf(domain, "./QuadData/domain0%i/subdomain0%i.dat",myid,myid);
    }
    else if(100<=myid<1000)
    {
    sprintf(coord, "./QuadData/domain%i/coord.dat",myid);
    sprintf(conn,  "./QuadData/domain%i/conn.dat",myid);
    sprintf(qedge, "./QuadData/domain%i/qedges.dat",myid);
    sprintf(ncolor,"./QuadData/domain%i/ncolors.dat",myid);
    sprintf(iqloop,"./QuadData/domain%i/iqloops.dat",myid);
    sprintf(qloop, "./QuadData/domain%i/qloops.dat",myid);
    sprintf(commu, "./QuadData/domain%i/quadCommu%i.dat",myid,myid);
    sprintf(domain, "./QuadData/domain%i/subdomain%i.dat",myid,myid);
    }

// ==================================================================
// read coordinates
// ==================================================================
   fp=fopen(coord,"r"); 
   fscanf(fp,"%d",&(g->nnodes));
  
   // allocate the required space for the coordinate points
   // for the x and y cartesian position
   g->x=(double *) malloc(sizeof(double)*2*g->nnodes);
   for(i=0;i<g->nnodes;i++)
      fscanf(fp,"%lf %lf %lf\n",&(g->x[2*i]),&(g->x[2*i+1]),&(fdummy));
   fclose(fp);

// ==================================================================
// read connectivity data
// set of three node IDs that make up a quadrilateral
// ==================================================================
   fp=fopen(conn,"r");
   fscanf(fp,"%d",&(g->ncells));
  
   g->conn=(int *) malloc(sizeof(int)*4*g->ncells);
   for(i=0;i<g->ncells;i++)    
      fscanf(fp,"%d %d %d %d\n",&(g->conn[4*i]),&(g->conn[4*i+1]),
     &(g->conn[4*i+2]),&(g->conn[4*i+3]));
    fclose(fp);

// ==================================================================
// read edges (faces)
// This is the Qedges matrix in the Matlab mesh generation
// NOTE: The ordering of the g->faces is different from the 
//       initial Matlab generation
//
// g->faces: 
//  0 - node ID 1 
//  1 - node ID 2
//  2 - left cell ID 
//  3 - left face index
//  4 - right celr ID 
//  5 - right face index
// ==================================================================
   fp=fopen(qedge,"r");
   fscanf(fp,"%d",&(g->nfaces));
   g->faces=(int *) malloc(sizeof(int)*6*g->nfaces);
   for(i=0;i<g->nfaces;i++)
   {
      fscanf(fp,"%d %d %d %d %d %d\n",
       &(g->faces[6*i]),    // node ID of one edge
       &(g->faces[6*i+1]),  // node ID of other edge
       &(g->faces[6*i+2]),  // left  cell
       &(g->faces[6*i+4]),  // right cell
       &(g->faces[6*i+3]),  // left  face index (0-3)
       &(g->faces[6*i+5])); // right face index (0-3)

      //
      // swap if left cell = -1 
      // Therefore as a result all -1s are not in [6*i+4]
      //if (g->faces[6*i+2] == -1 || g->faces[6*i+2] == -5)
      if (g->faces[6*i+2] == -1)


     {
         swap(g->faces[6*i],  g->faces[6*i+1]);
         swap(g->faces[6*i+2],g->faces[6*i+4]);
         swap(g->faces[6*i+3],g->faces[6*i+5]);  
      }
   }
   fclose(fp);

// ==================================================================
// read colours (i.e., number of loops per colour)
// ==================================================================
   fp=fopen(ncolor,"r");
   fscanf(fp,"%d",&(g->ncolors));
   g->chainsPerColor=(int *) malloc(sizeof(int)*g->ncolors);
   for(i=0;i<g->ncolors;i++)
      fscanf(fp,"%d\n",&(g->chainsPerColor[i]));
   fclose(fp);

// ==================================================================
// read chain information (chain = loops)
// ==================================================================
   fp = fopen(iqloop,"r");
   fscanf(fp,"%d\n",&(g->nchains));
   g->nchains--; // not sure why it is reduced by 1
   g->faceStartPerChain=(int *) malloc(sizeof(int)*(g->nchains+1));
   g->nmaxchain=0;

   for(i=0;i<g->nchains+1;i++)
   {
      fscanf(fp,"%d\n",&(g->faceStartPerChain[i]));
      if (i > 0) 
      {
         g->nmaxchain = max(g->nmaxchain,
         g->faceStartPerChain[i]-g->faceStartPerChain[i-1]);
      }
   }
   fclose(fp);
   //
   if(myid==0) trace(g->nmaxchain);   
   //
// ==================================================================
// read quad communication data 
// ==================================================================

  
   fp=fopen(domain,"r"); //read subdomain.dat
   fscanf(fp,"%d %d\n",&(g->ncommu),&(dummy)); // same number of send and recv
  // printf("ncommu:%d\n",g->ncommu); // same number of send and recv
   
   //skip character line   
   do { c=fgetc(fp);} while (c!='\n');
   
   g->icommu = (int *) malloc(sizeof(int)*8*g->ncommu);
   
   g->isend  = (int *) malloc(sizeof(int)*2*g->ncommu);
   g->irecv  = (int *) malloc(sizeof(int)*2*g->ncommu);
   g->irecvconn  = (int *) malloc(sizeof(int)*(g->ncommu));
   g->irecvconn2  = (int *) malloc(sizeof(int)*6*(g->ncommu));
   
   g->idup  = (int *) malloc(sizeof(int)*g->ncells);
   g->psil     = (double **) malloc(sizeof(double *)*(g->ncells));
   g->psila    = (double **) malloc(sizeof(double *)*(g->ncells));
   g->dpsil    = (double **) malloc(sizeof(double *)*(g->ncells));
   g->dpsila   = (double **) malloc(sizeof(double *)*(g->ncells));
   
   g->psilb    = (double **) malloc(sizeof(double *)*(g->ncells));
   g->dpsilb   = (double **) malloc(sizeof(double *)*(g->ncells));
   
   g->psil_deri    = (double **) malloc(sizeof(double *)*(g->ncells));
   g->psila_deri   = (double **) malloc(sizeof(double *)*(g->ncells));   
   g->psilb_deri   = (double **) malloc(sizeof(double *)*(g->ncells));
   g->dpsil_deri   = (double **) malloc(sizeof(double *)*(g->ncells));
   g->dpsila_deri  = (double **) malloc(sizeof(double *)*(g->ncells));
   g->dpsilb_deri  = (double **) malloc(sizeof(double *)*(g->ncells));
   
   for (i=0;i<g->ncells;i++)
   {
     g->psil[i]   = (double *) malloc(sizeof(double)*NVAR);
     g->psila[i]  = (double *) malloc(sizeof(double)*NVAR);
     g->psilb[i]  = (double *) malloc(sizeof(double)*NVAR);
     g->dpsil[i]  = (double *) malloc(sizeof(double)*NVAR);
     g->dpsila[i] = (double *) malloc(sizeof(double)*NVAR);
     g->dpsilb[i] = (double *) malloc(sizeof(double)*NVAR);

     g->psil_deri[i]   = (double *) malloc(sizeof(double)*NVAR);
     g->psila_deri[i]  = (double *) malloc(sizeof(double)*NVAR);
     g->psilb_deri[i]  = (double *) malloc(sizeof(double)*NVAR);
     g->dpsil_deri[i]  = (double *) malloc(sizeof(double)*NVAR);
     g->dpsila_deri[i] = (double *) malloc(sizeof(double)*NVAR);
     g->dpsilb_deri[i] = (double *) malloc(sizeof(double)*NVAR);
   }
   
   // initialize 
   for (i=0;i<g->ncells;i++)
   for (j=0;j<NVAR;j++)
   {
     g->psil[i][j]   = 0.;
     g->psila[i][j]  = 0.;
     g->psilb[i][j]  = 0.;
     
     g->dpsil[i][j]  = 0.;
     g->dpsila[i][j] = 0.;
     g->dpsilb[i][j] = 0.;


     g->psil_deri[i][j]   = 0.;
     g->psila_deri[i][j]  = 0.;
     g->dpsil_deri[i][j]  = 0.;
     g->dpsila_deri[i][j] = 0.;

   }  


   fscanf(fp,"%d\n",&(g->nadjp));

   g->iadjp  = (int *) malloc(sizeof(int)*g->nadjp);
   g->istp   = (int *) malloc(sizeof(int)*g->nadjp);
   g->ilengp = (int *) malloc(sizeof(int)*g->nadjp);
   
      
   for(i=0;i<g->nadjp;i++)
   fscanf(fp,"%d",&(g->iadjp[i]));
   for(i=0;i<g->nadjp;i++)
   fscanf(fp,"%d",&(g->istp[i]));
   for(i=0;i<g->nadjp;i++)
   fscanf(fp,"%d",&(g->ilengp[i]));

   for(i=0;i<g->ncommu;i++)
     fscanf(fp,"%d %d\n",&(g->isend[2*i]),&(g->isend[2*i+1]));
  
   //skip character line   
   do { c=fgetc(fp);} while (c!='\n');

   fscanf(fp,"%d\n",&(g->nadjp));

   for(i=0;i<g->nadjp;i++)
   fscanf(fp,"%d",&(g->iadjp[i]));
   for(i=0;i<g->nadjp;i++)
   fscanf(fp,"%d",&(g->istp[i]));
   for(i=0;i<g->nadjp;i++)
   fscanf(fp,"%d",&(g->ilengp[i]));


   for(i=0;i<g->ncommu;i++)
     fscanf(fp,"%d %d\n",&(g->irecv[2*i]),&(g->irecv[2*i+1]));

   fclose(fp);

   //===========================
   //find duplicate cells in  irecv 
   //==========================
   ndup = 0;
   for(i=0;i<g->ncells;i++) 
   g->idup[i] = 0;
   
   for(i=0;i<g->ncommu;i++)
   {
     c1 = g->irecv[i*2];
     if(g->idup[c1]==1) 
     { 
       ndup = ndup+1;
     }  
     g->idup[c1] = g->idup[c1]+1;
   }
   
   fp=fopen(commu,"r"); 
   
   for(i=0;i<g->ncommu;i++)
     fscanf(fp,"%d %d %d %d %d %d %d %d\n",
     &(g->icommu[8*i]),
     &(g->icommu[8*i+1]),
     &(g->icommu[8*i+2]),
     &(g->icommu[8*i+3]),
     &(g->icommu[8*i+4]),
     &(g->icommu[8*i+5]),
     &(g->icommu[8*i+6]),
     &(g->icommu[8*i+7]));
  
     fclose(fp);

   dup = (int *) malloc(sizeof(int)*g->ncommu);
   for(i=0;i<g->ncommu;i++) dup[i] = 0;
   for(i=0;i<g->ncommu;i++)
   {
     c1  = g->irecv[i*2];
     iproc = g->irecv[i*2+1];
     for(j=0;j<g->ncommu;j++)
     {
        //if(c1 == g->icommu[8*j+1] && iproc == g->icommu[8*j+6] && dup[j]==1)
        //printf("duplicate in irecv: proc:%d,ith:%d\n",myid,i);

        if(c1==g->icommu[8*j+1] && iproc == g->icommu[8*j+6] && dup[j] == 0)
        { 
           g->irecvconn[i] = g->icommu[8*j+2];
           dup[j] = dup[j]+1;
           break;
           //if(myid==16) printf("irecvconn[%d]: %d\n",i,g->irecvconn[i]);
        }
     }
   }

    // ======================
    // make isend global index
    //=======================
    for(i=0;i<g->ncommu;i++)
    {
      c1 = g->isend[2*i];
      iproc = g->isend[2*i+1];
      for(j=0;j<g->ncommu;j++)
      {
        
     //   printf("c1:%d,g->icommu:%d\n",c1,g->icommu[8*j+4]);
        if(c1==g->icommu[8*j+4] && iproc==g->icommu[8*j+6]) 
        {
          g->isend[2*i] = g->icommu[8*j+7];
        }
      }
    }


   //// connectivity inverse
   //for(i=0;i<g->ntcells;i++)
   //{     
   //  g->qconni[3*i]=-1;
   //  g->qconni[3*i+1]=-1;
   //  g->qconni[3*i+2]=-1;
   //}

//   for(i=0;i<g->ncommu;i++)
//   {
//     j = g->icommu[8*i+7];
//     if(g->qconni[3*j]!=-1) 
//     {
//     printf("===> duplicate: myid:%d,ith:%d\n<======",myid,i);
//     //exit(1);
//     }
//     g->qconni[3*j]   = g->icommu[8*i+1]; //second cell 
//     g->qconni[3*j+1] = g->icommu[8*i+2]; //first cell
//     g->qconni[3*j+2] = g->icommu[8*i];   //third cell
//   }
//    

    for(i=0;i<g->ncommu;i++) dup[i] = 0;
    for(i=0;i<g->ncommu;i++)
    {
      c1  = g->isend[2*i]; //global icell cell
      iproc = g->isend[2*i+1]; //processor number
      
      for(j=0;j<g->ncommu;j++)
      {
        
        //if(c1 == g->icommu[8*j+7] && iproc == g->icommu[8*j+6] && dup[j]==1)
        //printf("duplicate in isend: proc:%d,ith:%d\n",myid,i);
       
        if(c1 == g->icommu[8*j+7] && iproc == g->icommu[8*j+6] && dup[j]==0)
        { 
           g->irecvconn2[6*i]   = g->icommu[8*j];
           g->irecvconn2[6*i+1] = g->icommu[8*j+1];
           g->irecvconn2[6*i+2] = g->icommu[8*j+2];
           g->irecvconn2[6*i+3] = g->icommu[8*j+3];
           g->irecvconn2[6*i+4] = g->icommu[8*j+4];
           g->irecvconn2[6*i+5] = g->icommu[8*j+5];
           dup[j] = dup[j]+1;

           //if(myid==17) printf("%d %d %d\n", g->irecvconn2[6*i],g->irecvconn2[6*i+1],g->irecvconn2[6*i+2]);
           break;
        }


      }
    }
    free(dup);


    //for(i=0;i<g->ncommu;i++)
    //printf("%d, %d, %d\n",myid,g->irecvconn[2*i],g->irecvconn[2*i+1]);
   
    if(myid==0) printf("#ham2d: Finished reading communication data.........\n");
// ==================================================================
// Finish reading quad communication data 
// ==================================================================


//   workSize=g->nmaxchain+5;
   workSize=g->nmaxchain+6;

   g->ql    = (double **)  malloc(sizeof(double *)*(workSize));
   g->qr    = (double **)  malloc(sizeof(double *)*(workSize));
   g->dql   = (double **)  malloc(sizeof(double *)*(workSize));
   g->dqr   = (double **)  malloc(sizeof(double *)*(workSize));
   g->f     = (double **)  malloc(sizeof(double *)*(workSize));
   g->fv    = (double **)  malloc(sizeof(double *)*(workSize));
   g->df    = (double **)  malloc(sizeof(double *)*(workSize));
   g->f2    = (double **)  malloc(sizeof(double *)*(workSize));
   g->cindx = (int *)      malloc(sizeof(int)*(workSize));
   g->ctype = (int *)      malloc(sizeof(int)*(workSize));
   g->A     = (double ***) malloc(sizeof(double **)*(workSize));
   g->B     = (double ***) malloc(sizeof(double **)*(workSize));
   g->C     = (double ***) malloc(sizeof(double **)*(workSize));  
   g->F     = (double **)  malloc(sizeof(double *)*(workSize));
   g->Q     = (double **)  malloc(sizeof(double *)*(workSize));
   //
   for (i=0;i<workSize;i++)
   { 
      g->ql[i]  = (double *)  malloc(sizeof(double)*NVAR);
      g->qr[i]  = (double *)  malloc(sizeof(double)*NVAR);
      g->dql[i] = (double *)  malloc(sizeof(double)*NVAR);
      g->dqr[i] = (double *)  malloc(sizeof(double)*NVAR);
      g->f[i]   = (double *)  malloc(sizeof(double)*NVAR);
      g->fv[i]  = (double *)  malloc(sizeof(double)*NVAR);
      g->df[i]  = (double *)  malloc(sizeof(double)*NVAR);
      g->f2[i]  = (double *)  malloc(sizeof(double)*NVAR);
      g->A[i]   = (double **) malloc(sizeof(double *)*NQ);
      g->B[i]   = (double **) malloc(sizeof(double *)*NQ);
      g->C[i]   = (double **) malloc(sizeof(double *)*NQ);
      g->F[i]   = (double *)  malloc(sizeof(double)*NQ);
      g->Q[i]   = (double *)  malloc(sizeof(double)*NQ);
      for(j=0;j<NQ;j++)
     {
        g->A[i][j] = (double *) malloc(sizeof(double)*NQ);
        g->B[i][j] = (double *) malloc(sizeof(double)*NQ);
        g->C[i][j] = (double *) malloc(sizeof(double)*NQ);
     }
   }



// ==================================================================
// read chain information (chain connectivity)
// nchainFaces = nloops (from Matlab)
// ================================================================
  fp=fopen(qloop,"r");
  fscanf(fp,"%d\n",&(g->nchainFaces));
  g->chainConn=(int *)malloc(sizeof(int) * g->nchainFaces);

  for (i=0;i<g->nchainFaces;i++)
     fscanf(fp,"%d\n",&(g->chainConn[i]));
  fclose(fp);
   
  if(myid==0) 
  {
  printf("#ham2d: Finished reading files......................\n");
  trace(g->nnodes);
  trace(g->ncells);
  trace(g->nfaces);
  trace(g->ncolors);
  trace(g->nchains);
  trace(g->nchainFaces);
  }

}
// ##################################################################
// END OF FILE
// ##################################################################
