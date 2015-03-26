// ##################################################################
//
// computeRHS.c
//
// Compute the RHS, i.e., -R(q^n)
//
// Written by Dr. Jayanarayanan Sitaraman
// ##################################################################
#include <stdio.h>
#include <math.h>
#include "ham2dtypes.h"
#include "ham2dFunctionDefs.h"
#include <mpi.h>


void computeRHS(GRID *g,SOLN *s,double *l2rho, int myid)
{
  //
  int    i,j,k,m,f,n,kp;
  int    f1,f2;
  int    is,ie;
  int    iface;
  int    idv;
  int    chainSize;
  double ds[2];
  double leftState[NVAR];
  double rightState[NVAR];
  double leftState0[NVAR];
  double rightState0[NVAR];
  double consVar[NVAR];
  double flux[NVAR];
  double gm1=gamm-1.0;
  double gamma1=gamm;
  double specRadius;
  double faceVel=0.;
  double dsnorm,nynx,nx2ny;
  double rhoi;
  int    node1,node2,leftCell,rightCell,icell;
  double x1,y1,x2,y2;  
  double pp;
  double th,qt,eps;
  double dscheck[2];
  FILE   *fp;
  int    iface1,iface2,rightCell1,rightCell2,c1,c2,c3;
 
  int nghost,order,iflag,error;
  //  fp = fopen("validate.dat","w");
  // set number of ghost cells
  order = g->order;  
               nghost = 2; // for MUSCL or WENO3
  if(order==5) nghost = 3; // for WENO 

  //
  // zero out residual and spectral radii
  //
  dscheck[0]=dscheck[1]=0;
  for(i=0;i<NVAR*g->ncells;i++) s->r[i]=0.0;
  for(i=0;i<g->ncells;i++) s->sigma[i]=0.0;
  //
  // one loop per chain to evaluate fluxes
  // on all the faces in the chain
  //
  
  for(i=0;i<g->nchains;i++)
  {
    iflag=0;
    // starting and ending indices of loop
    f1 = g->faceStartPerChain[i];
    f2 = g->faceStartPerChain[i+1];

    m = nghost;
 
    // loop through the faces of a loop
    for (f=f1;f<f2;f++)
    {
      iface       = g->chainConn[f];     // face index
      g->cindx[m] = g->faces[6*iface+2]; // left cell index
      m++;
    }

    //
    // add buffer cells to the chain
    //
    if (g->chainConn[f1]==g->chainConn[f2-1])
    {
      // this is a closed chain
      // make it periodic
      //
      iflag = 0;
      f           = f1+1;
      iface       = g->chainConn[f];     // face index
      g->cindx[m] = g->faces[6*iface+2]; // left cell index
      m++;      
      chainSize   = m;

      m = 0;

      for (f=f2-nghost-1;f<f2-1;f++)
      {
        iface       = g->chainConn[f];
        g->cindx[m] = g->faces[6*iface+2];
        m++;
      }

    }

    else 
    {
      //
      // this is a open chain
      // -ve index indicates necessity to create
      // ghost cells
      //

      iflag=1;

      
      //solid bc
      if(g->test == 0){
        if(order==5) //WENO 5
        {     
        m--;    
        g->cindx[m] = -g->cindx[m];
        m++;
        g->cindx[m] = -g->cindx[m-3];
        m++;
        g->cindx[m] = -g->cindx[m-5];

        chainSize = m+1;
        m = 0;
        g->cindx[m] = -g->cindx[m+5];
        m = 1;
        g->cindx[m] = -g->cindx[m+3];
        m = 2;
        g->cindx[m] = -g->cindx[m+1];
        }
        else
        {
        m--;
        g->cindx[m] = -g->cindx[m];
        m++;
        g->cindx[m] = -g->cindx[m-3];
        chainSize = m+1;
        m = 0;
        g->cindx[m] = -g->cindx[m+3];
        m = 1;
        g->cindx[m] = -g->cindx[m+1];
        }
      }
         
          //periodic bc
      if(g->test==1) {
         apply_periodic(&g[0],f1,f2,m);
         chainSize = m+1;          
         if(order==5) chainSize = m+2; 
     }

      //===================================
      // add for parallization
      //======================================
      iface1       = g->chainConn[f1];
      iface2       = g->chainConn[f2-1];
      rightCell1   = g->faces[6*iface1+4];
      rightCell2   = g->faces[6*iface2+4];

      if(rightCell1==-5)
      {  

        if(order==5)
        {
          c3 = g->cindx[5];
          c2 = g->cindx[4];
          c1 = g->cindx[3];
          
          g->cindx[0] = -1*(c2);
          g->cindx[1] = -1*(c2);
          g->cindx[2] = -1*(c2);
        }
        else
        {
        // 2nd cell from the last physical cell        
        c2  = g->cindx[3];
        c1  = g->cindx[2];

        g->cindx[0] = -1*(c2); //
        g->cindx[1] = -1*(c2); //
        }
      }

      if(rightCell2==-5)
      {  

        if(order==5)
        {

          c3 = g->cindx[chainSize-6];
          c2 = g->cindx[chainSize-5];
          c1 = g->cindx[chainSize-4];
          g->cindx[chainSize-1] = -1*(c2); // 
          g->cindx[chainSize-2] = -1*(c2); //
          g->cindx[chainSize-3] = -1*(c2); //

        }
        else
        {
          c2 = g->cindx[chainSize-4];
          c1 = g->cindx[chainSize-3];
          g->cindx[chainSize-1] = -1*(c2); //
          g->cindx[chainSize-2] = -1*(c2); //
        }
      }
      // end add for parallization
      //====================================
//=================================================== 
//============================================
    } // end open/closed chain

    //if(myid==0){
    //if(iflag==1) //open chain 
    //{ 
    //printf("Stopping code\n");
    //printf("chain size=%d\n",chainSize);
    //for(k=0;k<=chainSize-1;k++) 
    //{
    //printf("cell index:%d\n",g->cindx[k]);  
    //}
    //}
    //}


//==================================================

    // go through each of the cells in a chain
    for (j=0;j<chainSize;j++)
    {
      icell = g->cindx[j]; // extract cell index (left)
      if (icell >=0 ) // not negative 
      //if (icell >=0 ||(icell==0&&j==nghost)||(icell==0&&j==chainSize-nghost-1)) // not negative 
      {
        m = NVAR*icell;

        // create temporary array of conservative variables
        // extracted from its location in s->q array
        for(k=0;k<NVAR;k++)
        {
          consVar[k] = s->q[m];
          m++;
        }

        rhoi=1./consVar[0]; // inverse of rho
                
        g->f[j][0] = consVar[0];
        g->f[j][1] = consVar[1];
        g->f[j][2] = consVar[2];
        g->f[j][3] = consVar[3];
      }
      else // icell < 0
      //if(icell<0||(icell==0&&j==nghost-1)||(icell==0&&j==chainSize-nghost))// icell < 0
      {
        //
        // do ghost cells
        // based on whether they are on the solid boundary on that
        if (j < nghost) //0,1,2
        {         
          iface = g->chainConn[f1];
        }
        else
        {
          iface = g->chainConn[f2-1];
        }
        
        rightCell = g->faces[6*iface+4];
       
               
        if (rightCell == -2)  /* this is a face on solid wall */
        {
          node1 = g->faces[6*iface];
          node2 = g->faces[6*iface+1];
      
          x1 = g->x[2*node1];
          y1 = g->x[2*node1+1];
          x2 = g->x[2*node2];
          y2 = g->x[2*node2+1];
          
          ds[0] =  (y2-y1);
          ds[1] = -(x2-x1);

          icell = -icell; // make icell positive
          m     = NVAR*icell;
          
          // temporate variable for conservative quantities
          for(k=0;k<NVAR;k++)
          {
            consVar[k] = s->q[m];
            m++;
          }

          dsnorm = ds[0]*ds[0]+ds[1]*ds[1];
          nynx   = ds[0]*ds[1]/dsnorm;
          nx2ny  = (ds[0]*ds[0]-ds[1]*ds[1])/dsnorm;
          rhoi   = 1./consVar[0];
          
          g->f[j][0] = consVar[0];
          g->f[j][1] = (-consVar[1]*nx2ny-2*consVar[2]*nynx);
          g->f[j][2] = (consVar[2]*nx2ny-2*consVar[1]*nynx);
          g->f[j][3] = consVar[3];            
        }
        //============================================
        // add for parallization
        //=============================================
        else if(rightCell == -5)
        {
          
          icell = -icell; //make positive(second cell index)
          if(j==nghost-2 || j==chainSize-(nghost-1))
          {  
            if(g->idup[icell]==2)
            {
              if(j==nghost-2) c1 = g->cindx[nghost];
              if(j==chainSize-(nghost-1)) c1 = g->cindx[chainSize-(nghost+1)];

              g->f[j][0] = g->dpsil[c1][0];
              g->f[j][1] = g->dpsil[c1][1]; 
              g->f[j][2] = g->dpsil[c1][2]; 
              g->f[j][3] = g->dpsil[c1][3];

            }
            else
            { 
              g->f[j][0] = g->psil[icell][0];
              g->f[j][1] = g->psil[icell][1]; 
              g->f[j][2] = g->psil[icell][2]; 
              g->f[j][3] = g->psil[icell][3];
               
            }

          }
          else if(j==nghost-1 || j==chainSize-(nghost)) //fist cell 
          {
            if(g->idup[icell]==2)
            {
              if(j==nghost-1) c1 = g->cindx[nghost];
              if(j==chainSize-nghost) c1 = g->cindx[chainSize-(nghost+1)];

              g->f[j][0] = g->dpsila[c1][0];
              g->f[j][1] = g->dpsila[c1][1]; 
              g->f[j][2] = g->dpsila[c1][2]; 
              g->f[j][3] = g->dpsila[c1][3];

            }
            else
            {

              g->f[j][0] = g->psila[icell][0];
              g->f[j][1] = g->psila[icell][1]; 
              g->f[j][2] = g->psila[icell][2]; 
              g->f[j][3] = g->psila[icell][3];

            }
          }

          if(order==5)
          {
            if(g->idup[icell]==2)
            {
              if(j==0) c1 = g->cindx[3];
              if(j==chainSize-1) c1 = g->cindx[chainSize-4];

              g->f[j][0] = g->dpsilb[c1][0];
              g->f[j][1] = g->dpsilb[c1][1]; 
              g->f[j][2] = g->dpsilb[c1][2]; 
              g->f[j][3] = g->dpsilb[c1][3];
            }
            else
            { 
              g->f[j][0] = g->psilb[icell][0];
              g->f[j][1] = g->psilb[icell][1]; 
              g->f[j][2] = g->psilb[icell][2]; 
              g->f[j][3] = g->psilb[icell][3];
            }
          }
        }
        //end add for parallization
        //=============================================
        else // if icell<0 and rightCell is not -2, than far field?
        {
          g->f[j][0] = rinf;
          g->f[j][1] = rinf*s->uinf;
          g->f[j][2] = rinf*s->vinf;
          g->f[j][3] = pinf/gm1+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf);
        }        
      } //icell is positive or not
      if(g->f[j][0]==0) printf("ERROR::::\n");
    }

//====================================================================
//==================END of storing flow variables ====================
//===================================================================
    is  = nghost-1;           // start index of chain
    ie  = chainSize-1; // end index of chain
    th  = 1./3;        // third (can be defined global/static)
    qt  = 0.25;        // quarter (same as above)
    if (g->order==1) qt = 0.0; 
    eps = 1e-10;      // epsilon for Koren's limiter


    // Perform MUSCL interpolation with Koren's limiter
    // output is g->ql and g->qr are not computed for the faces 
    // of a chain

    if(order==1 || order==3) 
    {
     muscld(g->f,g->ql,g->qr,g->f2,is,ie,th,qt,eps,chainSize,NVAR);
    }
    
     if(order==5) 
     {
      weno(g->f,g->ql,g->qr,is,ie,eps,chainSize,NVAR); 
     }
    n = is; // start index of chain

    // if idv = 1, it is a closed chain, 0 if open
    idv = (g->chainConn[f1]==g->chainConn[f2-1]);
    // Loop through all the unique faces of a chain
    // Because in closed loops the "starting" face is repeated
    // at the end


    for(f=f1;f<f2-idv;f++)
    {
      iface     = g->chainConn[f];
      node1     = g->faces[6*iface];
      node2     = g->faces[6*iface+1];
      leftCell  = g->faces[6*iface+2];
      rightCell = g->faces[6*iface+4];
      x1        = g->x[2*node1];
      y1        = g->x[2*node1+1];
      x2        = g->x[2*node2];
      y2        = g->x[2*node2+1];
      ds[0]     = (y2-y1);
      ds[1]     = -(x2-x1);

      // loop through the number of variables
      // n is a running index of the chain faces
      for(m=0;m<NVAR;m++)
      {
        // if open chain and is the last element of the chain
        if (f==f2-idv-1 && idv==0) 
        { 
          leftState[m]  = g->ql[n][m];
          rightState[m] = g->qr[n+1][m];
        }
        // any other face
        else
        {
          leftState[m]  = g->qr[n+1][m];
          rightState[m] = g->ql[n][m];
        }

        leftState0[m] = s->q[NVAR*leftCell+m]; //g->ql[j][m];             
        
        if (rightCell > -1) // periodic condition always do!! 
        {
          rightState0[m] = s->q[NVAR*rightCell+m]; //g->qr[j+1][m];
        }  
      }
    
      if (rightCell==-1 && g->test==0) //periodic condition always do not!! 
      { 
        rightState0[0] = rinf;
        rightState0[1] = rinf*s->uinf;
        rightState0[2] = rinf*s->vinf;
        rightState0[3] = pinf/gm1+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf);
      }
      else if (rightCell==-2 && g->test==0) //
      {
        dsnorm         = ds[0]*ds[0]+ds[1]*ds[1];
        nynx           = ds[0]*ds[1]/dsnorm;
        nx2ny          = (ds[0]*ds[0]-ds[1]*ds[1])/dsnorm;
        rightState0[0] = leftState0[0];
        rightState0[1] = -leftState0[1]*nx2ny-2*leftState0[2]*nynx;
        rightState0[2] = leftState0[2]*nx2ny-2*leftState0[1]*nynx;
        rightState0[3] = leftState0[3];
      }
    
      // the '&& 0' condition makes this a false loop always
      if (rightState[1]!=rightState0[1] && 0) 
      {
        trace(leftCell);
        trace(rightCell);
        for(k=0;k<chainSize;k++)
          printf("%d %d %.16e\n",k,g->cindx[k],g->f[k][1]);
        trace(n);
        tracef(leftState[1]);
        tracef(leftState0[1]);
        tracef(rightState[1]);
        tracef(rightState0[1]);
        tracef(s->q[leftCell*NVAR+1])
        trace(i);
        trace(n);
        trace(f);
        exit(0);
      }

      //
      //roeflx(&specRadius,flux,leftState,rightState,faceVel,ds,gm1);

      // Compute Roe flux in 2d one cell at a time
      // Outputs the Roe flux and spectral radius
      flux_roe2d_(ds,leftState,rightState,flux,&specRadius,&gamma1);
      //printf("-----------------------------\n"); 
      //printf("%e %e %e %e\n",flux[0],flux[1],flux[2],flux[3]);
      //printf("%e\n",specRadius);
      //exit(1);   
    
      //
      m = NVAR*leftCell;
      // Compute the residual array for the cells
      for(j=0;j<NVAR;j++)
      {
        s->r[m] -= flux[j];
        m++;

      }

      // accumulate the spectral radius for the left cell
      //if(leftCell==0&&myid==0) printf("lspec:%f,rightcell:%d\n",specRadius,rightCell);
      s->sigma[leftCell] += specRadius;
      

      // if NOT a boundary cell
      
      if (rightCell > -1) //periodic bc always do!! 
      
      {
        m = NVAR*rightCell;
        for(j=0;j<NVAR;j++)
        {
          // Compute the residual array for the cells
          s->r[m] += flux[j];
          m++;          
        }

        // accumulate the spectral radius for the right cell
        s->sigma[rightCell]+=specRadius;
    // if(rightCell==0&&myid==0) printf("rspec:%f,leftCell:%d\n",specRadius,leftCell);
      } 

      n++; // n is a running index of the chain faces (is to ie-1)


    if(i==0) {
  //printf("specRadius:%f\n",specRadius);
  //printf("leftcell:%d,ds[x]:%e, ds[y]:%e\n",leftCell,ds[0],ds[1]);
  //printf("rightcell:%d,ds[x]:%e, ds[y]:%e\n",rightCell,ds[0],ds[1]);
//         m = NVAR*leftCell;
//         printf("leftCell:%d,  res:%e,%e,%e,%e\n",leftCell,s->r[m],s->r[m+1],s->r[m+2],s->r[m+3]);
//         m = NVAR*rightCell;
//         printf("rightCell:%d,    res:%e,%e,%e,%e\n",rightCell,s->r[m],s->r[m+1],s->r[m+2],s->r[m+3]);
          }
    }

  //  exit(1);
  } // end nchains loop


  *l2rho=0.;
  for(i=0;i<g->ncells;i++)
    {
      if ((*l2rho) < fabs(s->r[4*i])) 
    {
      icell  = i; // why have this?
      *l2rho = fabs(s->r[4*i]);
    }
    }

  //if(myid==0)
  //{
  //printf("myid:%d\n",myid);
  //for(i=0;i<g->ncells;i++) printf("sigma[%d]:%lf\n",i,s->sigma[i]);
  //}
  //error = MPI_Barrier(MPI_COMM_WORLD);

  //if(myid==1)
  //{
  // printf("myid:%d\n",myid);
  //for(i=0;i<g->ncells;i++) printf("sigma[%d]:%lf\n",i,s->sigma[i]);
  //}
  //tracef(s->sigma[0]);
  //trace(icell);

//    fclose(fp);

} // end function
// ##################################################################      
// END OF FILE
// ##################################################################  
