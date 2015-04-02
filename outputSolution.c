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
void outputSolution(GRID *g,SOLN *s,int nn,int myid, int nproc)
{
  int    i,n;
  int    iface,node1,node2;
  int    icell;
  double x1,y1,x2,y2,rho,rhou,rhov,e,pp,cp;
  FILE   *fp,*fp1;
  char   fname[80];
  int    iproc,ierr;
  //==============\
  // Surface plot  >
  //==============/
  if      (nn < 10  ) {sprintf(fname,"./output/output00%d_%d.dat",nn,myid);}
  else if (nn < 100 ) {sprintf(fname,"./output/output0%d_%d.dat",nn,myid);}
  else if (nn < 1000) {sprintf(fname,"./output/output%d_%d.dat",nn,myid);}

  fp=fopen(fname,"w");

  fprintf(fp,"VARIABLES = \"X\",\"Y\",\"RHO\",\"RHOU\",\"RHOV\",\"E\"\n");
  
  fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N=%d E=%d DATAPACKING=BLOCK\n",g->nnodes,g->ncells); 

  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=CELLCENTERED, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED)\n");

  // list of x-coordinate positions (after smoothing)
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i]);
  fprintf(fp,"\n");

  // list of y-coordinate positions (after smoothing)
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i+1]);
  fprintf(fp,"\n");

  // rho, rhoU, rhoV, E for all the cells
  for(n=0;n<NVAR;n++)
    for (i=0;i<g->ncells;i++)
    {
       fprintf(fp,"%f\n",s->q[4*i+n]);
    }
  fprintf(fp,"\n");

  // connectivity information for the cells
  for(i=0;i<g->ncells;i++)
    fprintf(fp,"%d %d %d %d\n",g->conn[4*i]+1,g->conn[4*i+1]+1,g->conn[4*i+2]+1,g->conn[4*i+3]+1);
  fclose(fp);

// ==================================================================
// Multiple procs (stitch subdomains for tecplot output)
// ==================================================================
	int     ii,j,k,totNumNode,totNumCell,stride,temp;
	int    *allNumNode    = (int *) malloc(sizeof(int)*nproc);
	int    *allNumNodeCum = (int *) malloc(sizeof(int)*nproc);
	int    *allNumCell    = (int *) malloc(sizeof(int)*nproc);
	int    *allNumCellCum = (int *) malloc(sizeof(int)*nproc);
	int    *myConn, *allConn,*subdomainconn;;
	double *allPos, *allSoln;

	MPI_Allgather(&g->nnodes, 1, MPI_INT, allNumNode, 1, MPI_INT, MPI_COMM_WORLD);
	MPI_Allgather(&g->ncells, 1, MPI_INT, allNumCell, 1, MPI_INT, MPI_COMM_WORLD);

	// compute the total number of nodes and cells (by a single proc)
	if (myid == 0)
	{
		totNumNode = totNumCell = 0;
		for (i = 0; i < nproc; i++)
		{
			totNumNode += allNumNode[i];
			totNumCell += allNumCell[i];
		}
	}
	
	// broadcast totNumNode and totNumCell to all procs (needed?)
	MPI_Bcast(&totNumNode,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&totNumCell,1,MPI_INT,0,MPI_COMM_WORLD);

	allNumNodeCum[0] = allNumCellCum[0] = 0;
	for (i = 1; i < nproc; i++)
	{
		allNumNodeCum[i] = allNumNodeCum[i-1] + allNumNode[i-1];
		allNumCellCum[i] = allNumCellCum[i-1] + allNumCell[i-1];
	}

	// Allgatherv the pos, soln, and conn values
	allPos = (double *) malloc(sizeof(double)*2*totNumNode);

	int *recvcounts   = (int *) malloc(sizeof(int)*nproc);
	int *displacement = (int *) malloc(sizeof(int)*nproc);

	//
	// allPos gather
	//
	for (i = 0; i < nproc; i++)
		recvcounts[i] = 2*allNumNode[i];

	displacement[0] = 0;
	for (i = 1; i < nproc; i++)
		displacement[i] = displacement[i-1] + recvcounts[i-1];

	MPI_Allgatherv(g->x, 2*g->nnodes, MPI_DOUBLE, allPos,
		recvcounts, displacement, MPI_DOUBLE, MPI_COMM_WORLD);

	//
	// allSoln gather
	//
	allSoln = (double *) malloc(sizeof(double)*4*totNumCell);
	for (i = 0; i < nproc; i++)
		recvcounts[i] = 4*allNumCell[i];

	displacement[0] = 0;
	for (i = 1; i < nproc; i++)
		displacement[i] = displacement[i-1] + recvcounts[i-1];

	MPI_Allgatherv(s->q, 4*g->ncells, MPI_DOUBLE, allSoln,
		recvcounts, displacement, MPI_DOUBLE, MPI_COMM_WORLD);

	//
	// allConn gather
	//
	myConn  = (int *) malloc(sizeof(int)*4*g->ncells);
	allConn = (int *) malloc(sizeof(int)*4*totNumCell);

	temp = 4*g->ncells;
	for (i = 0; i < temp; i++)
	{
		myConn[i] = g->conn[i] + allNumNodeCum[myid];
	}

	MPI_Allgatherv(myConn, 4*g->ncells, MPI_INT, allConn,
		recvcounts, displacement, MPI_INT, MPI_COMM_WORLD);

// ==================================================================
// Increase the conn number for subsequent domains and update
// the conn values based on subdomainconn.dat (from meshgen)
// ==================================================================
	if (myid == 0)
	{
		int nodeID;
		// rewrite connvalues based on subdomainconn.dat
		subdomainconn = (int *) malloc(sizeof(int)*totNumNode);

		fp = fopen("./QuadData/domain000/subdomainconn.dat","r");
		
		if (fp == NULL)
		{
			printf("WARNING: subdomainconn.dat file missing. Recheck. Stopping.\n");
			MPI_Abort(MPI_COMM_WORLD,33);
		}
		fscanf(fp,"%*d"); // ignore the number

		for (i = 0; i < totNumNode; i++)
			fscanf(fp,"%d",&subdomainconn[i]);

		fclose(fp);

		for (i = 0; i < totNumCell; i++)
		{
			for (j = 0; j < NVAR; j++)
			{
				temp   = 4*i+j;
				nodeID = allConn[temp];

				if (subdomainconn[nodeID] != -1)
					allConn[temp] = subdomainconn[nodeID];

			}
		}


	   //
	   // Write outputs into file
	   //
      fp=fopen("./output/outputall.dat","w");
      
      fprintf(fp,"VARIABLES = \"X\",\"Y\",\"RHO\",\"RHOU\",\"RHOV\",\"E\"\n");
      fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N=%d E=%d DATAPACKING=BLOCK\n",totNumNode,totNumCell); 
      fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=CELLCENTERED, "
      	"4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED)\n");

      // list of x-coordinate positions 
      for(i = 0; i < totNumNode; i++)
        fprintf(fp,"%f\n",allPos[2*i]);
      fprintf(fp,"\n");

      // list of y-coordinate positions
      for(i = 0; i < totNumNode; i++)
        fprintf(fp,"%f\n",allPos[2*i+1]);
      fprintf(fp,"\n");

      // rho, rhoU, rhoV, E for all the cells
      for(n = 0; n < NVAR; n++)
        for (i = 0; i < totNumCell; i++)
        {
           fprintf(fp,"%f\n",allSoln[4*i+n]);
        }
      fprintf(fp,"\n");

      // connectivity information for the cells
      for(i = 0; i < totNumCell; i++)
      {
        fprintf(fp,"%d %d %d %d\n",
        	allConn[4*i]+1,allConn[4*i+1]+1,allConn[4*i+2]+1,allConn[4*i+3]+1);
      }

      fclose(fp);



	}

}

// ==================================================================
// 
// ==================================================================
void outputdq(GRID *g,SOLN *s)
{
  int i,n;
  int iface,node1,node2;
  int icell;
  double x1,y1,x2,y2,rho,rhou,rhov,e,pp,cp;
  FILE *fp;
  char fname[80];
  static int istep0=0;
  //
  sprintf(fname,"dq%d.plt",istep0);
  fp=fopen(fname,"w");
  fprintf(fp,"VARIABLES = \"X\",\"Y\",\"RHO\",\"RHOU\",\"RHOV\",\"E\"\n");
  fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N=%d E=%d DATAPACKING=BLOCK\n",g->nnodes,g->ncells); 
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=CELLCENTERED, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED)\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i+1]);
  fprintf(fp,"\n");
  for(n=0;n<NVAR;n++)
    for(i=0;i<g->ncells;i++)
      fprintf(fp,"%f\n",s->ddq[4*i+n]);
  fprintf(fp,"\n");
  for(i=0;i<g->ncells;i++)
    fprintf(fp,"%d %d %d %d\n",g->conn[4*i]+1,g->conn[4*i+1]+1,g->conn[4*i+2]+1,g->conn[4*i+3]+1);
  fclose(fp);
  istep0++;
}

// ==================================================================
// 
// ==================================================================
void outputr(GRID *g,SOLN *s)
{
  int i,n;
  int iface,node1,node2;
  int icell;
  double x1,y1,x2,y2,rho,rhou,rhov,e,pp,cp;
  FILE *fp;
  char fname[80];
  static int istep0=0;
  //
  sprintf(fname,"r%d.plt",istep0);
  fp=fopen(fname,"w");
  fprintf(fp,"VARIABLES = \"X\",\"Y\",\"RHO\",\"RHOU\",\"RHOV\",\"E\"\n");
  fprintf(fp,"ZONE ZONETYPE=FEQUADRILATERAL N=%d E=%d DATAPACKING=BLOCK\n",g->nnodes,g->ncells); 
  fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=CELLCENTERED, 4=CELLCENTERED, 5=CELLCENTERED, 6=CELLCENTERED)\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[2*i+1]);
  fprintf(fp,"\n");
  for(n=0;n<NVAR;n++)
    for(i=0;i<g->ncells;i++)
      fprintf(fp,"%f\n",s->r[4*i+n]);
  fprintf(fp,"\n");
  for(i=0;i<g->ncells;i++)
    fprintf(fp,"%d %d %d %d\n",g->conn[4*i]+1,g->conn[4*i+1]+1,g->conn[4*i+2]+1,g->conn[4*i+3]+1);
  fclose(fp);
  istep0++;
}
// ##################################################################
// END OF FILE
// ##################################################################s
